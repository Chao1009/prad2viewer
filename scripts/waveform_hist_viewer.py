#!/usr/bin/env python3
"""
Waveform Histogram Viewer
=========================
Self-contained tool: loads an evio file with ``prad2py`` (showing a
cancellable progress dialog), or a previously-saved JSON, and plots the
three per-channel peak histograms shown on the Waveform tab of the
prad2 event monitor:

    - Peak height
    - Peak integral
    - Peak time

Peak finding is a faithful Python port of
``prad2dec/src/WaveAnalyzer.cpp`` (triangular smoothing, iterative
pedestal with outlier rejection, local-maxima search, tail-cut
integration).

Time cut (``-t``)
-----------------
Without ``-t``: every peak above the height threshold fills every
histogram.
With ``-t MIN[,MAX]``: only peaks inside the window fill the time
histogram; the height & integral hists take the best peak (by
integral) in the window per event — the server's logic.

Usage
-----
    python scripts/waveform_hist_viewer.py RUN.evio.00000
    python scripts/waveform_hist_viewer.py RUN.evio.00000 -t 170,190 -n 200000
    python scripts/waveform_hist_viewer.py saved_hist.json
    python scripts/waveform_hist_viewer.py             # File → Open… from menu
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QComboBox, QCheckBox, QCompleter, QFileDialog, QMessageBox,
    QProgressDialog, QSizePolicy, QStatusBar, QToolTip,
)
from PyQt6.QtCore import (
    Qt, QObject, QRectF, QSize, QThread, pyqtSignal,
)
from PyQt6.QtGui import (
    QAction, QPainter, QColor, QPen, QBrush, QFont, QPalette,
)


# ===========================================================================
#  prad2py discovery (mirrors tagger_viewer.py)
# ===========================================================================

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR   = _SCRIPT_DIR.parent

for _cand in (
    _REPO_DIR / "build" / "python",
    _REPO_DIR / "build-release" / "python",
    _REPO_DIR / "build" / "Release" / "python",
):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

try:
    import prad2py                                # type: ignore
    _HAVE_PRAD2PY = True
    _PRAD2PY_ERR  = ""
except Exception as _exc:
    prad2py = None                                # type: ignore
    _HAVE_PRAD2PY = False
    _PRAD2PY_ERR  = f"{type(_exc).__name__}: {_exc}"


# ===========================================================================
#  WaveAnalyzer port — prad2dec/src/WaveAnalyzer.cpp
# ===========================================================================

@dataclass
class WaveConfig:
    resolution:      int   = 2
    threshold:       float = 5.0
    min_threshold:   float = 3.0
    min_peak_ratio:  float = 0.3
    int_tail_ratio:  float = 0.1
    ped_nsamples:    int   = 30
    ped_flatness:    float = 1.0
    ped_max_iter:    int   = 3
    overflow:        int   = 4095
    clk_mhz:         float = 250.0
    max_peaks:       int   = 8


@dataclass
class Peak:
    pos:      int
    left:     int
    right:    int
    height:   float
    integral: float
    time:     float
    overflow: bool


def _smooth(raw: np.ndarray, res: int) -> np.ndarray:
    n = raw.size
    if res <= 1:
        return raw.astype(np.float32)
    out = np.empty(n, dtype=np.float32)
    rf = float(res)
    for i in range(n):
        val  = float(raw[i])
        wsum = 1.0
        for j in range(1, res):
            if j > i or i + j >= n:
                continue
            w = 1.0 - j / (rf + 1.0)
            val  += w * (float(raw[i - j]) + float(raw[i + j]))
            wsum += 2.0 * w
        out[i] = val / wsum
    return out


def _pedestal(buf: np.ndarray, cfg: WaveConfig) -> Tuple[float, float]:
    n = min(cfg.ped_nsamples, buf.size)
    if n <= 0:
        return 0.0, 0.0
    s = buf[:n].astype(np.float64).copy()
    mean = float(s.mean())
    var  = float(s.var())
    rms  = np.sqrt(var) if var > 0 else 0.0
    for _ in range(cfg.ped_max_iter):
        cut = max(rms, cfg.ped_flatness)
        keep = np.abs(s - mean) < cut
        count = int(keep.sum())
        if count == s.size or count < 5:
            break
        s = s[keep]
        mean = float(s.mean())
        var  = float(s.var())
        rms  = np.sqrt(var) if var > 0 else 0.0
    return mean, rms


def _trend_sign(d: float) -> int:
    if abs(d) < 0.1:
        return 0
    return 1 if d > 0 else -1


def _find_peaks(raw: np.ndarray, buf: np.ndarray,
                ped_mean: float, ped_rms: float, thr: float,
                cfg: WaveConfig) -> List[Peak]:
    n = buf.size
    if n < 3:
        return []
    peaks: List[Peak] = []
    pk_range: List[Tuple[int, int]] = []

    i = 1
    while i < n - 1 and len(peaks) < cfg.max_peaks:
        tr1 = _trend_sign(float(buf[i]) - float(buf[i - 1]))
        tr2 = _trend_sign(float(buf[i]) - float(buf[i + 1]))
        if tr1 * tr2 < 0 or (tr1 == 0 and tr2 == 0):
            i += 1
            continue

        flat_end = i
        if tr2 == 0:
            while (flat_end < n - 1 and
                   _trend_sign(float(buf[flat_end]) - float(buf[flat_end + 1])) == 0):
                flat_end += 1
            if (flat_end >= n - 1 or
                _trend_sign(float(buf[flat_end]) - float(buf[flat_end + 1])) <= 0):
                i += 1
                continue
        peak_pos = (i + flat_end) // 2

        left, right = i, flat_end
        while left > 0 and _trend_sign(float(buf[left]) - float(buf[left - 1])) > 0:
            left -= 1
        while (right < n - 1 and
               _trend_sign(float(buf[right]) - float(buf[right + 1])) >= 0):
            right += 1

        span = right - left
        if span <= 0:
            i += 1
            continue

        base = (float(buf[left])  * (right - peak_pos) +
                float(buf[right]) * (peak_pos - left)) / span
        smooth_height = float(buf[peak_pos]) - base
        if smooth_height < thr:
            i = right
            continue

        ped_height = float(buf[peak_pos]) - ped_mean
        if ped_height < thr or ped_height < 3.0 * ped_rms:
            i = right
            continue

        integral = float(buf[peak_pos]) - ped_mean
        tail_cut = ped_height * cfg.int_tail_ratio
        int_left, int_right = peak_pos, peak_pos
        for j in range(peak_pos - 1, left - 1, -1):
            v = float(buf[j]) - ped_mean
            if v < tail_cut or v < ped_rms or v * ped_height < 0:
                int_left = j
                break
            integral += v
            int_left = j
        for j in range(peak_pos + 1, right + 1):
            v = float(buf[j]) - ped_mean
            if v < tail_cut or v < ped_rms or v * ped_height < 0:
                int_right = j
                break
            integral += v
            int_right = j

        raw_pos = peak_pos
        raw_height = float(raw[peak_pos]) - ped_mean
        search = max(1, cfg.resolution) + (flat_end - i) // 2
        for j in range(1, search + 1):
            if peak_pos - j >= 0:
                h = float(raw[peak_pos - j]) - ped_mean
                if h > raw_height:
                    raw_height = h
                    raw_pos = peak_pos - j
            if peak_pos + j < n:
                h = float(raw[peak_pos + j]) - ped_mean
                if h > raw_height:
                    raw_height = h
                    raw_pos = peak_pos + j

        rejected = False
        for k, (lk, rk) in enumerate(pk_range):
            if left <= rk and right >= lk:
                if smooth_height < peaks[k].height * cfg.min_peak_ratio:
                    rejected = True
                    break
        if rejected:
            i = right
            continue

        peaks.append(Peak(
            pos      = raw_pos,
            left     = int_left,
            right    = int_right,
            height   = raw_height,
            integral = integral,
            time     = raw_pos * 1e3 / cfg.clk_mhz,
            overflow = raw[raw_pos] >= cfg.overflow,
        ))
        pk_range.append((left, right))
        i = right

    return peaks


def analyze(samples: np.ndarray, cfg: WaveConfig) -> Tuple[float, float, List[Peak]]:
    n = samples.size
    if n <= 0:
        return 0.0, 0.0, []
    buf = _smooth(samples, cfg.resolution)
    ped_mean, ped_rms = _pedestal(buf, cfg)
    thr = max(cfg.threshold * ped_rms, cfg.min_threshold)
    peaks = _find_peaks(samples, buf, ped_mean, ped_rms, thr, cfg)
    return ped_mean, ped_rms, peaks


# ===========================================================================
#  Accumulator data classes
# ===========================================================================

@dataclass
class Hist1D:
    nbins:   int
    bmin:    float = 0.0
    bstep:   float = 1.0
    bins:    np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    under:   int = 0
    over:    int = 0

    def __post_init__(self):
        if self.bins.size == 0:
            self.bins = np.zeros(self.nbins, dtype=np.int64)

    def fill(self, v: float):
        if v < self.bmin:
            self.under += 1
            return
        b = int((v - self.bmin) / self.bstep)
        if b >= self.nbins:
            self.over += 1
            return
        self.bins[b] += 1

    def to_json(self) -> Dict:
        return {"bins": self.bins.tolist(),
                "underflow": int(self.under), "overflow": int(self.over)}

    @classmethod
    def from_json(cls, d: Dict, bmin: float, bstep: float) -> "Hist1D":
        bins = np.asarray(d.get("bins", []), dtype=np.int64)
        h = cls(nbins=bins.size, bmin=bmin, bstep=bstep)
        h.bins  = bins
        h.under = int(d.get("underflow", 0))
        h.over  = int(d.get("overflow", 0))
        return h


@dataclass
class ChannelHists:
    roc:         int
    slot:        int
    channel:     int
    module:      Optional[str] = None
    events:      int = 0
    peak_events: int = 0
    tcut_events: int = 0
    height:      Optional[Hist1D] = None
    integral:    Optional[Hist1D] = None
    position:    Optional[Hist1D] = None


# ===========================================================================
#  Config / map loaders
# ===========================================================================

def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    with open(path, encoding="utf-8") as f:
        entries = json.load(f)
    return {(int(e["crate"]), int(e["slot"]), int(e["channel"])): e["name"]
            for e in entries}


def load_roc_tag_map(path: Path) -> Dict[int, int]:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    out = {}
    for r in cfg.get("roc_tags", []):
        tag = int(r["tag"], 16) if isinstance(r["tag"], str) else int(r["tag"])
        if r.get("type") == "roc":
            out[tag] = int(r["crate"])
    return out


def load_hist_config(path: Path) -> Dict:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return cfg.get("waveform", {})


def load_trigger_bit_map(path: Path) -> Dict[str, int]:
    if not path.is_file():
        return {}
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return {e["name"]: int(e["bit"]) for e in cfg.get("trigger_bits", [])}


def _mask_from_names(names: List[str], bitmap: Dict[str, int]) -> int:
    m = 0
    for n in names:
        if n in bitmap:
            m |= (1 << bitmap[n])
        else:
            print(f"  warning: trigger bit {n!r} not in trigger_bits.json",
                  file=sys.stderr)
    return m


# ===========================================================================
#  Load worker — runs in a QThread
# ===========================================================================

class LoadWorker(QObject):
    progressed = pyqtSignal(int, int)   # events_seen, channels_seen
    finished   = pyqtSignal(object)     # dict payload
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str,
                 hist_cfg: Dict,
                 daq_map: Dict[Tuple[int, int, int], str],
                 roc_to_crate: Dict[int, int],
                 accept_mask: int,
                 reject_mask: int,
                 max_events: int,
                 wcfg: WaveConfig,
                 time_cut: Optional[Tuple[float, float]],
                 daq_config_path: str = ""):
        super().__init__()
        self._path = evio_path
        self._cfg = hist_cfg
        self._daq_map = daq_map
        self._roc_to_crate = roc_to_crate
        self._accept = accept_mask
        self._reject = reject_mask
        self._max = max_events
        self._wcfg = wcfg
        self._tcut = time_cut
        self._daq_cfg_path = daq_config_path
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            payload = self._run()
            self.finished.emit(payload)
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> Dict:
        if not _HAVE_PRAD2PY:
            raise RuntimeError(
                f"prad2py not importable ({_PRAD2PY_ERR}).  "
                "Build with: cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build")

        h_cfg  = self._cfg.get("height_hist",   {"min": 0, "max": 4000,  "step": 10})
        i_cfg  = self._cfg.get("integral_hist", {"min": 0, "max": 20000, "step": 100})
        p_cfg  = self._cfg.get("time_hist",     {"min": 0, "max": 400,   "step": 4})
        thr_cfg = self._cfg.get("thresholds",   {})
        hist_threshold = float(thr_cfg.get("min_peak_height", 10.0))
        self._wcfg.min_peak_ratio = float(thr_cfg.get("min_secondary_peak_ratio",
                                                      self._wcfg.min_peak_ratio))

        def nbins(c):
            span = c["max"] - c["min"]
            return max(1, int(np.ceil(span / c["step"])))
        h_nbins, i_nbins, p_nbins = nbins(h_cfg), nbins(i_cfg), nbins(p_cfg)

        tcut = self._tcut
        t_min, t_max = (tcut if tcut is not None else (None, None))

        channels: Dict[Tuple[int, int, int], ChannelHists] = {}
        total_events = 0

        dec = prad2py.dec
        cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
               else dec.load_daq_config())
        ch  = dec.EvChannel()
        ch.set_config(cfg)
        st = ch.open(self._path)
        if st != dec.Status.success:
            raise RuntimeError(f"cannot open {self._path}: {st}")

        progress_every = 100
        stop = False
        try:
            while ch.read() == dec.Status.success:
                if self._cancel:
                    stop = True
                    break
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for i in range(ch.get_n_events()):
                    ch.select_event(i)
                    info = ch.info()
                    tb = int(info.trigger_bits)
                    if self._accept and (tb & self._accept) == 0:
                        continue
                    if self._reject and (tb & self._reject):
                        continue
                    fadc_evt = ch.fadc()

                    for r in range(fadc_evt.nrocs):
                        roc = fadc_evt.roc(r)
                        roc_tag = int(roc.tag)
                        if roc_tag not in self._roc_to_crate:
                            continue
                        crate = self._roc_to_crate[roc_tag]
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                samples = slot.channel(c).samples
                                if samples.size < 10:
                                    continue
                                key = (roc_tag, s, c)
                                hits = channels.get(key)
                                if hits is None:
                                    hits = ChannelHists(
                                        roc=roc_tag, slot=s, channel=c,
                                        module=self._daq_map.get((crate, s, c)),
                                        height  =Hist1D(h_nbins, h_cfg["min"], h_cfg["step"]),
                                        integral=Hist1D(i_nbins, i_cfg["min"], i_cfg["step"]),
                                        position=Hist1D(p_nbins, p_cfg["min"], p_cfg["step"]),
                                    )
                                    channels[key] = hits
                                hits.events += 1

                                _, _, peaks = analyze(samples, self._wcfg)

                                if tcut is None:
                                    any_peak = False
                                    for p in peaks:
                                        if p.height < hist_threshold:
                                            continue
                                        any_peak = True
                                        hits.position.fill(p.time)
                                        hits.integral.fill(p.integral)
                                        hits.height.fill(p.height)
                                    if any_peak:
                                        hits.peak_events += 1
                                else:
                                    best_int = -1.0
                                    best_hgt = -1.0
                                    any_peak = False
                                    for p in peaks:
                                        if p.height < hist_threshold:
                                            continue
                                        any_peak = True
                                        if t_min <= p.time <= t_max:
                                            hits.position.fill(p.time)
                                            if p.integral > best_int:
                                                best_int = p.integral
                                                best_hgt = p.height
                                    if any_peak:
                                        hits.peak_events += 1
                                    if best_int >= 0:
                                        hits.tcut_events += 1
                                        hits.integral.fill(best_int)
                                        hits.height.fill(best_hgt)

                    total_events += 1

                    if total_events % progress_every == 0:
                        self.progressed.emit(total_events, len(channels))
                        if self._cancel:
                            stop = True
                            break

                    if self._max and total_events >= self._max:
                        stop = True
                        break
                if stop:
                    break
        finally:
            ch.close()

        # Final tick so the dialog shows the exact numbers
        self.progressed.emit(total_events, len(channels))

        return {
            "source_file":   str(self._path),
            "total_events":  total_events,
            "cancelled":     self._cancel,
            "time_cut":      (None if tcut is None
                              else {"min": float(t_min), "max": float(t_max)}),
            "height_hist":   {"min": h_cfg["min"], "max": h_cfg["max"], "step": h_cfg["step"], "nbins": h_nbins},
            "integral_hist": {"min": i_cfg["min"], "max": i_cfg["max"], "step": i_cfg["step"], "nbins": i_nbins},
            "position_hist": {"min": p_cfg["min"], "max": p_cfg["max"], "step": p_cfg["step"], "nbins": p_nbins},
            "threshold":     hist_threshold,
            "wave_config":   self._wcfg.__dict__.copy(),
            "channels":      channels,
        }


# ===========================================================================
#  Hist1DWidget — QPainter bar chart
# ===========================================================================

class Hist1DWidget(QWidget):
    PAD_L, PAD_R, PAD_T, PAD_B = 58, 14, 26, 34

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(150)
        self._bins: np.ndarray = np.zeros(0, dtype=np.int64)
        self._bmin: float = 0.0
        self._bstep: float = 1.0
        self._under: int = 0
        self._over: int = 0
        self._title: str = ""
        self._xlabel: str = ""
        self._color: QColor = QColor("#58a6ff")
        self._log_y: bool = False
        self._hover_idx: int = -1

        # Per-widget log-Y toggle — floats in the top-right corner.
        self._logy_cb = QCheckBox("log Y", self)
        self._logy_cb.setFont(QFont("Monospace", 9))
        self._logy_cb.setStyleSheet(
            "QCheckBox{color:#8b949e;background:transparent;}"
            "QCheckBox:hover{color:#c9d1d9;}")
        self._logy_cb.toggled.connect(self.set_log_y)
        self._logy_cb.adjustSize()
        self._logy_cb.raise_()

    def set_data(self, bins, bmin: float, bstep: float,
                 under: int = 0, over: int = 0,
                 title: str = "", xlabel: str = "",
                 color: Optional[str] = None):
        self._bins = np.asarray(bins, dtype=np.int64)
        self._bmin = float(bmin)
        self._bstep = float(bstep)
        self._under = int(under)
        self._over  = int(over)
        self._title = title
        self._xlabel = xlabel
        if color:
            self._color = QColor(color)
        self._hover_idx = -1
        self.update()

    def set_log_y(self, on: bool):
        if on != self._log_y:
            self._log_y = on
            self.update()

    def clear(self, title: str = ""):
        self._bins = np.zeros(0, dtype=np.int64)
        self._title = title
        self._under = self._over = 0
        self._hover_idx = -1
        self.update()

    def _plot_rect(self) -> QRectF:
        w, h = self.width(), self.height()
        return QRectF(self.PAD_L, self.PAD_T,
                      max(1.0, w - self.PAD_L - self.PAD_R),
                      max(1.0, h - self.PAD_T - self.PAD_B))

    def resizeEvent(self, ev):
        cb = self._logy_cb
        cb.adjustSize()
        cb.move(self.width() - cb.width() - 6, 4)
        super().resizeEvent(ev)

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        p.fillRect(self.rect(), QColor("#0d1117"))

        r = self._plot_rect()
        p.setPen(QColor("#30363d"))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        if self._title:
            f = QFont("Monospace", 10); f.setBold(True)
            p.setFont(f)
            p.setPen(QColor("#c9d1d9"))
            p.drawText(int(r.left()), int(r.top() - 8), self._title)

        n = self._bins.size
        if n == 0 or self._bins.sum() == 0:
            p.setPen(QColor("#6e7681"))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "(no data)")
            return

        if self._log_y:
            vals = np.where(self._bins > 0, np.log10(self._bins.astype(np.float64)), 0.0)
            ymin, ymax = 0.0, float(vals.max())
        else:
            vals = self._bins.astype(np.float64)
            ymin, ymax = 0.0, float(vals.max())
        if ymax <= 0:
            ymax = 1.0

        bar_w = r.width() / n
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(self._color)
        for i in range(n):
            v = vals[i]
            if v <= 0:
                continue
            h = (v - ymin) / (ymax - ymin) * r.height()
            if h < 1.0:
                continue
            x0 = r.left() + i * bar_w
            y0 = r.bottom() - h
            p.fillRect(QRectF(x0, y0, max(bar_w, 1.0), h), self._color)

        if 0 <= self._hover_idx < n:
            x0 = r.left() + self._hover_idx * bar_w
            p.setPen(QPen(QColor("#ffffff"), 1.2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(QRectF(x0, r.top(), max(bar_w, 1.0), r.height()))

        p.setPen(QColor("#6e7681"))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            if self._log_y:
                val = 10 ** (ymin + frac * (ymax - ymin))
            else:
                val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4), _fmt_count(val))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            x = r.left() + frac * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            val = self._bmin + frac * n * self._bstep
            p.drawText(int(x - 24), int(r.bottom() + 14), f"{val:g}")

        if self._xlabel:
            p.setFont(QFont("Monospace", 9))
            p.drawText(int(r.left() + r.width() / 2 - 50),
                       int(r.bottom() + 28), self._xlabel)

        entries = int(self._bins.sum())
        info = f"N={entries:,}"
        if self._under:
            info += f"  under={self._under:,}"
        if self._over:
            info += f"  over={self._over:,}"
        p.setFont(QFont("Monospace", 9))
        p.setPen(QColor("#8b949e"))
        # Right-align so the text hugs the plot's right edge without
        # ever reaching the log-Y checkbox in the widget corner.
        info_rect = QRectF(r.left(), r.top() - 20,
                           max(1.0, self.width() - self._logy_cb.width() - 20 - r.left()),
                           14)
        p.drawText(info_rect,
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   info)

    def mouseMoveEvent(self, ev):
        r = self._plot_rect()
        n = self._bins.size
        if n == 0 or not r.contains(ev.position()):
            if self._hover_idx != -1:
                self._hover_idx = -1
                QToolTip.hideText()
                self.update()
            return
        idx = int((ev.position().x() - r.left()) / r.width() * n)
        idx = max(0, min(n - 1, idx))
        if idx != self._hover_idx:
            self._hover_idx = idx
            self.update()
        lo = self._bmin + idx * self._bstep
        hi = lo + self._bstep
        QToolTip.showText(ev.globalPosition().toPoint(),
                          f"[{lo:g}, {hi:g})  count={int(self._bins[idx]):,}",
                          self)

    def leaveEvent(self, _ev):
        if self._hover_idx != -1:
            self._hover_idx = -1
            QToolTip.hideText()
            self.update()


def _fmt_count(v: float) -> str:
    av = abs(v)
    if av == 0:
        return "0"
    if av >= 1e6:
        return f"{v/1e6:.2f}M"
    if av >= 1e3:
        return f"{v/1e3:.1f}k"
    if av >= 10:
        return f"{v:.0f}"
    return f"{v:.2g}"


# ===========================================================================
#  Main window
# ===========================================================================

class WaveformHistViewerWindow(QMainWindow):

    def __init__(self,
                 *,
                 hist_config: Dict,
                 daq_map: Dict,
                 roc_to_crate: Dict,
                 accept_mask: int,
                 reject_mask: int,
                 max_events: int,
                 time_cut: Optional[Tuple[float, float]],
                 daq_config_path: str):
        super().__init__()
        self._hist_config    = hist_config
        self._daq_map        = daq_map
        self._roc_to_crate   = roc_to_crate
        self._accept_mask    = accept_mask
        self._reject_mask    = reject_mask
        self._max_events     = max_events
        self._time_cut       = time_cut
        self._daq_config_path = daq_config_path

        self._payload: Optional[Dict] = None        # current in-memory data
        self._items: List[Tuple[str, str]] = []      # (display_label, channel_key)
        self._worker: Optional[LoadWorker] = None
        self._thread: Optional[QThread] = None

        self._apply_dark_palette()
        self._build_ui()
        self._make_menu()

    # -- dark theme --

    def _apply_dark_palette(self):
        pal = self.palette()
        for role, colour in [
            (QPalette.ColorRole.Window,     "#0d1117"),
            (QPalette.ColorRole.WindowText, "#c9d1d9"),
            (QPalette.ColorRole.Base,       "#161b22"),
            (QPalette.ColorRole.Text,       "#c9d1d9"),
            (QPalette.ColorRole.Button,     "#21262d"),
            (QPalette.ColorRole.ButtonText, "#c9d1d9"),
            (QPalette.ColorRole.Highlight,  "#58a6ff"),
        ]:
            pal.setColor(role, QColor(colour))
        self.setPalette(pal)

    # -- UI --

    def _build_ui(self):
        self.setWindowTitle("Waveform Histogram Viewer")
        self.resize(1100, 900)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        self._file_lbl = QLabel("(no file loaded)")
        self._file_lbl.setFont(QFont("Monospace", 10))
        self._file_lbl.setStyleSheet("color:#8b949e;")
        root.addWidget(self._file_lbl)

        picker_row = QHBoxLayout()
        lbl = QLabel("Module:")
        lbl.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
        lbl.setStyleSheet("color:#c9d1d9;")
        picker_row.addWidget(lbl)

        self._combo = QComboBox()
        self._combo.setEditable(True)
        self._combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._combo.setFont(QFont("Monospace", 11))
        self._combo.setMinimumContentsLength(40)
        self._combo.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "selection-background-color:#1f6feb;}")
        comp = self._combo.completer()
        if comp is not None:
            comp.setFilterMode(Qt.MatchFlag.MatchContains)
            comp.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            comp.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        self._combo.currentIndexChanged.connect(self._on_combo_changed)
        picker_row.addWidget(self._combo, stretch=1)
        root.addLayout(picker_row)

        self._info = QLabel("")
        self._info.setFont(QFont("Monospace", 10))
        self._info.setStyleSheet("color:#8b949e;")
        root.addWidget(self._info)

        self._h_height   = Hist1DWidget()
        self._h_integral = Hist1DWidget()
        self._h_position = Hist1DWidget()
        root.addWidget(self._h_height,   stretch=1)
        root.addWidget(self._h_integral, stretch=1)
        root.addWidget(self._h_position, stretch=1)

        self.setStatusBar(QStatusBar())
        self._clear_plots()

    def _make_menu(self):
        mb = self.menuBar()
        mf = mb.addMenu("&File")

        a_evio = QAction("Open &evio…", self)
        a_evio.setShortcut("Ctrl+O")
        a_evio.triggered.connect(self._open_evio_dialog)
        mf.addAction(a_evio)

        a_json = QAction("Open &JSON…", self)
        a_json.setShortcut("Ctrl+J")
        a_json.triggered.connect(self._open_json_dialog)
        mf.addAction(a_json)

        self._a_save = QAction("&Save histograms as JSON…", self)
        self._a_save.setShortcut("Ctrl+S")
        self._a_save.triggered.connect(self._save_json_dialog)
        self._a_save.setEnabled(False)
        mf.addAction(self._a_save)

        mf.addSeparator()
        a_quit = QAction("&Quit", self)
        a_quit.setShortcut("Ctrl+Q")
        a_quit.triggered.connect(self.close)
        mf.addAction(a_quit)

    # -- dispatch file open --

    def open_path(self, path: Path):
        if path.suffix.lower() == ".json":
            self._load_json(path)
        else:
            self._load_evio(path)

    def _open_evio_dialog(self):
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open evio file", str(Path.cwd()),
            "evio files (*.evio *.evio.*);;All files (*)")
        if path_str:
            self._load_evio(Path(path_str))

    def _open_json_dialog(self):
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open histogram JSON", str(Path.cwd()),
            "JSON files (*.json);;All files (*)")
        if path_str:
            self._load_json(Path(path_str))

    def _save_json_dialog(self):
        if not self._payload:
            return
        default = Path(self._payload.get("source_file") or "waveform_hist").name + ".json"
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save histograms as JSON", str(Path.cwd() / default),
            "JSON files (*.json)")
        if not path_str:
            return
        try:
            self._save_json(Path(path_str))
        except Exception as ex:
            QMessageBox.warning(self, "Save failed", f"{path_str}\n\n{ex}")
            return
        self.statusBar().showMessage(f"Saved {path_str}")

    # -- evio load path --

    def _load_evio(self, path: Path):
        if self._thread is not None:
            QMessageBox.information(self, "Busy", "Already loading a file.")
            return
        if not _HAVE_PRAD2PY:
            QMessageBox.critical(
                self, "prad2py missing",
                "prad2py is not importable.\n\n"
                f"{_PRAD2PY_ERR}\n\n"
                "Build with:\n"
                "  cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build")
            return

        pmax = int(self._max_events) if self._max_events > 0 else 1_000_000
        dlg = QProgressDialog(f"Decoding {path.name} …", "Cancel", 0, pmax, self)
        dlg.setWindowTitle("Loading")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        # Show immediately instead of the default 4-second delay — waveform
        # analysis is slow and the first progress signal can be >10 s away.
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()
        self.centralWidget().setEnabled(False)
        self.menuBar().setEnabled(False)

        worker = LoadWorker(
            str(path),
            hist_cfg        = self._hist_config,
            daq_map         = self._daq_map,
            roc_to_crate    = self._roc_to_crate,
            accept_mask     = self._accept_mask,
            reject_mask     = self._reject_mask,
            max_events      = self._max_events,
            wcfg            = WaveConfig(),
            time_cut        = self._time_cut,
            daq_config_path = self._daq_config_path,
        )
        thread = QThread(self)
        worker.moveToThread(thread)

        unlimited = (self._max_events == 0)
        def _on_progress(events: int, channels: int):
            if unlimited:
                m = dlg.maximum()
                if m > 0 and events * 5 >= m * 4:
                    dlg.setMaximum(max(events * 2, m * 2))
            dlg.setValue(events)
            dlg.setLabelText(
                f"Decoding {path.name}\nEvents: {events:,}   Channels: {channels:,}")

        thread.started.connect(worker.run)
        worker.progressed.connect(_on_progress)
        worker.finished.connect(lambda payload: self._on_load_finished(path, payload))
        worker.failed.connect(lambda msg: self._on_load_failed(path, msg))
        dlg.canceled.connect(worker.request_cancel)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(dlg.close)
        thread.finished.connect(self._restore_ui)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)

        self._worker = worker
        self._thread = thread
        thread.start()

    def _restore_ui(self):
        self.centralWidget().setEnabled(True)
        self.menuBar().setEnabled(True)
        self._worker = None
        self._thread = None

    def _on_load_finished(self, path: Path, payload: Dict):
        self._payload = payload
        ch = payload.get("channels", {})
        # Worker returns a dict keyed by (roc, slot, ch) tuples — convert to
        # the JSON-style string key map so _show_channel() can look them up
        # the same way it does for JSON loads.
        string_keyed = {}
        for k, hh in ch.items():
            if isinstance(k, tuple):
                string_keyed[f"{k[0]}_{k[1]}_{k[2]}"] = hh
            else:
                string_keyed[str(k)] = hh
        payload["channels"] = string_keyed
        self._after_data_loaded(path, payload, source_label="(evio)")
        cancelled = bool(payload.get("cancelled"))
        msg = (f"Cancelled after {payload['total_events']:,} events — "
               f"{len(string_keyed):,} channels"
               if cancelled else
               f"Loaded {payload['total_events']:,} events, "
               f"{len(string_keyed):,} channels from {path.name}")
        self.statusBar().showMessage(msg)

    def _on_load_failed(self, path: Path, msg: str):
        QMessageBox.critical(self, "Load failed", f"{path}\n\n{msg}")
        self.statusBar().showMessage(f"Failed to load {path.name}")

    # -- JSON load path --

    def _load_json(self, path: Path):
        try:
            with open(path, encoding="utf-8") as f:
                obj = json.load(f)
        except Exception as ex:
            QMessageBox.warning(self, "Load failed", f"{path}\n\n{ex}")
            return
        if not isinstance(obj, dict) or "channels" not in obj:
            QMessageBox.warning(self, "Load failed",
                                f"{path}\n\nNot a waveform_hist_viewer JSON.")
            return

        hh = obj["height_hist"]; ih = obj["integral_hist"]; ph = obj["position_hist"]
        channels = {}
        for key, c in obj["channels"].items():
            ch = ChannelHists(
                roc=int(c["roc"]), slot=int(c["slot"]), channel=int(c["channel"]),
                module=c.get("module"),
                events=int(c.get("events", 0)),
                peak_events=int(c.get("peak_events", 0)),
                tcut_events=int(c.get("tcut_events", 0)),
                height  =Hist1D.from_json(c.get("height_hist", {}),   hh["min"], hh["step"]),
                integral=Hist1D.from_json(c.get("integral_hist", {}), ih["min"], ih["step"]),
                position=Hist1D.from_json(c.get("position_hist", {}), ph["min"], ph["step"]),
            )
            channels[key] = ch
        obj["channels"] = channels
        self._payload = obj
        self._after_data_loaded(path, obj, source_label="(json)")
        self.statusBar().showMessage(
            f"Loaded {len(channels):,} channels from {path.name}")

    # -- JSON save --

    def _save_json(self, path: Path):
        p = self._payload
        if not p:
            return
        out = {
            "source_file":  p.get("source_file", ""),
            "total_events": int(p.get("total_events", 0)),
            "cancelled":    bool(p.get("cancelled", False)),
            "time_cut":     p.get("time_cut"),
            "height_hist":  p["height_hist"],
            "integral_hist":p["integral_hist"],
            "position_hist":p["position_hist"],
            "threshold":    p.get("threshold"),
            "wave_config":  p.get("wave_config", {}),
            "channels":     {},
        }
        for key, hh in sorted(p["channels"].items()):
            out["channels"][key] = {
                "module":        hh.module,
                "roc":           hh.roc,
                "slot":          hh.slot,
                "channel":       hh.channel,
                "events":        hh.events,
                "peak_events":   hh.peak_events,
                "tcut_events":   hh.tcut_events,
                "height_hist":   hh.height.to_json(),
                "integral_hist": hh.integral.to_json(),
                "position_hist": hh.position.to_json(),
            }
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(out, f)

    # -- post-load common --

    def _after_data_loaded(self, path: Path, payload: Dict, source_label: str):
        n_events = int(payload.get("total_events", 0))
        n_chan   = len(payload["channels"])
        self._file_lbl.setText(
            f"{source_label} {path.name}   events={n_events:,}   channels={n_chan:,}")
        tcut = payload.get("time_cut")
        thr  = payload.get("threshold")
        parts = []
        if thr is not None:
            parts.append(f"threshold={thr:g}")
        if tcut:
            parts.append(f"time_cut=[{tcut['min']:g}, {tcut['max']:g}] ns")
        else:
            parts.append("time_cut=none (every peak counts)")
        self._info.setText("   ".join(parts))

        self._populate_combo()
        self._a_save.setEnabled(True)
        if self._items:
            self._combo.setCurrentIndex(0)     # triggers _on_combo_changed

    def _populate_combo(self):
        self._combo.blockSignals(True)
        self._combo.clear()
        self._items = []
        # Sort: modules with a name first (by module name), then unmapped by
        # roc/slot/ch. Inside each group, stable ordering.
        channels: Dict[str, ChannelHists] = self._payload["channels"]
        named = []
        unnamed = []
        for key, hh in channels.items():
            label_left = (hh.module or "(unmapped)")
            label = f"{label_left:<8}  roc=0x{hh.roc:02X}  s={hh.slot:>2}  ch={hh.channel:>2}"
            item = (label, key)
            (named if hh.module else unnamed).append((label_left, item))
        named.sort(key=lambda x: _natural_sort_key(x[0]))
        unnamed.sort(key=lambda x: x[0])
        all_items = [it for _, it in named] + [it for _, it in unnamed]
        self._items = all_items
        for label, _ in all_items:
            self._combo.addItem(label)
        self._combo.blockSignals(False)

    # -- combo → plot --

    def _on_combo_changed(self, idx: int):
        if not self._payload or idx < 0 or idx >= len(self._items):
            return
        _, key = self._items[idx]
        self._show_channel(key)

    def _show_channel(self, key: str):
        p = self._payload
        hh: Optional[ChannelHists] = p["channels"].get(key)
        if hh is None:
            self._clear_plots()
            return
        module = hh.module or "(unmapped)"
        self._info.setText(
            f"{module}   roc=0x{hh.roc:02X}  slot={hh.slot}  ch={hh.channel}   "
            f"events={hh.events:,}   peaks={hh.peak_events:,}   tcut={hh.tcut_events:,}"
        )
        self._h_height.set_data(hh.height.bins, hh.height.bmin, hh.height.bstep,
                                under=hh.height.under, over=hh.height.over,
                                title=f"{module}  —  Peak Height",
                                xlabel="ADC", color="#e599f7")
        self._h_integral.set_data(hh.integral.bins, hh.integral.bmin, hh.integral.bstep,
                                  under=hh.integral.under, over=hh.integral.over,
                                  title=f"{module}  —  Peak Integral",
                                  xlabel="ADC·sample", color="#00b4d8")
        self._h_position.set_data(hh.position.bins, hh.position.bmin, hh.position.bstep,
                                  under=hh.position.under, over=hh.position.over,
                                  title=f"{module}  —  Peak Time",
                                  xlabel="ns", color="#51cf66")

    def _clear_plots(self):
        self._h_height.clear("Peak Height")
        self._h_integral.clear("Peak Integral")
        self._h_position.clear("Peak Time")

    # -- close --

    def closeEvent(self, ev):
        if self._worker is not None:
            self._worker.request_cancel()
        if self._thread is not None and self._thread.isRunning():
            self._thread.quit()
            self._thread.wait(3000)
        super().closeEvent(ev)


# Natural sort so "W2" < "W10" < "W100".
_NATKEY_RE = re.compile(r"(\d+)")
def _natural_sort_key(s: str):
    return [int(p) if p.isdigit() else p.lower()
            for p in _NATKEY_RE.split(s or "")]


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Waveform-tab histograms — load from evio or JSON and browse.")
    ap.add_argument("path", nargs="?", type=Path,
                    help="Either an evio file (*.evio.*) or a saved JSON "
                         "from a previous run. If omitted, use File → Open…")
    ap.add_argument("-t", "--time-cut", type=str, default=None, metavar="MIN[,MAX]",
                    help="Time-cut window (ns). Without -t every peak fills every "
                         "histogram. With -t MIN,MAX the height & integral hists take "
                         "the best peak (by integral) inside the window, and the time "
                         "hist keeps only peaks inside it. -t MIN is MIN,+infinity.")
    ap.add_argument("-n", "--max-events", type=int, default=1_000_000,
                    help="Stop decoding after N physics events "
                         "(default 1,000,000; pass 0 for unlimited).")
    ap.add_argument("--config", type=Path,
                    default=_REPO_DIR / "database" / "config.json",
                    help="Main config.json (for waveform binning)")
    ap.add_argument("--daq-config", type=Path,
                    default=_REPO_DIR / "database" / "daq_config.json",
                    help="daq_config.json (ROC-tag → crate mapping)")
    ap.add_argument("--daq-map", type=Path,
                    default=_REPO_DIR / "database" / "daq_map.json",
                    help="daq_map.json (module-name lookup)")
    ap.add_argument("--trigger-bits", type=Path,
                    default=_REPO_DIR / "database" / "trigger_bits.json",
                    help="trigger_bits.json (for --accept/--reject-trigger names)")
    ap.add_argument("--accept-trigger", action="append", default=[],
                    metavar="NAME",
                    help="Require at least one of these trigger bits (repeatable).")
    ap.add_argument("--reject-trigger", action="append", default=None,
                    metavar="NAME",
                    help="Drop events with any of these trigger bits (repeatable). "
                         "Default: uses config.json setting.")
    args = ap.parse_args()

    time_cut: Optional[Tuple[float, float]] = None
    if args.time_cut is not None:
        parts = [p.strip() for p in args.time_cut.split(",") if p.strip()]
        try:
            if len(parts) == 1:
                time_cut = (float(parts[0]), float("inf"))
            elif len(parts) == 2:
                time_cut = (float(parts[0]), float(parts[1]))
            else:
                raise ValueError
        except ValueError:
            ap.error(f"--time-cut expects MIN or MIN,MAX (got {args.time_cut!r})")
        if time_cut[0] >= time_cut[1]:
            ap.error(f"--time-cut MIN must be < MAX (got {time_cut})")

    hist_cfg     = load_hist_config(args.config) if args.config.is_file() else {}
    roc_to_crate = load_roc_tag_map(args.daq_config) if args.daq_config.is_file() else {}
    daq_map      = load_daq_map(args.daq_map) if args.daq_map.is_file() else {}
    bit_map      = load_trigger_bit_map(args.trigger_bits)

    accept_names = args.accept_trigger or hist_cfg.get("accept_trigger_bits", []) or []
    if args.reject_trigger is None:
        reject_names = hist_cfg.get("reject_trigger_bits", []) or []
    else:
        reject_names = args.reject_trigger
    accept_mask = _mask_from_names(accept_names, bit_map) if accept_names else 0
    reject_mask = _mask_from_names(reject_names, bit_map) if reject_names else 0

    app = QApplication(sys.argv)
    win = WaveformHistViewerWindow(
        hist_config     = hist_cfg,
        daq_map         = daq_map,
        roc_to_crate    = roc_to_crate,
        accept_mask     = accept_mask,
        reject_mask     = reject_mask,
        max_events      = args.max_events,
        time_cut        = time_cut,
        daq_config_path = str(args.daq_config) if args.daq_config.is_file() else "",
    )
    win.show()
    if args.path is not None:
        win.open_path(args.path)
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
