#!/usr/bin/env python3
"""
Waveform Viewer
===============
Browses an evio file event-by-event.  Opens the file in evio's
random-access mode (no event processing at open time), indexes physics
sub-events, and lets the user step through them.  For each viewed
event, the selected module's waveform is drawn and the four
per-module histograms (peak height, integral, time, n-peaks) accumulate.

The "Process next 10k" button runs a background pass that fills the
current module's histograms without displaying waveforms, for fast
accumulation.

Usage
-----
    python scripts/waveform_viewer.py RUN.evio.00000
    python scripts/waveform_viewer.py             # File → Open…
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QComboBox, QCheckBox, QCompleter, QFileDialog, QMessageBox,
    QProgressDialog, QSizePolicy, QStatusBar, QToolTip, QPushButton, QSpinBox,
    QSplitter,
)
from PyQt6.QtCore import (
    Qt, QObject, QPointF, QRectF, QSize, QThread, pyqtSignal, QTimer,
)
from PyQt6.QtGui import (
    QAction, QKeySequence, QPainter, QColor, QPen, QBrush, QFont, QPolygonF,
    QShortcut,
)

from hycal_geoview import (
    Module as GeoModule, load_modules as load_geo_modules,
    HyCalMapWidget, cmap_qcolor, apply_theme_palette, set_theme,
    available_themes, THEME, themed,
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


def _check_random_access_support() -> Optional[str]:
    """Return None if the loaded prad2py supports random-access, else an error."""
    if not _HAVE_PRAD2PY:
        return _PRAD2PY_ERR
    try:
        ch = prad2py.dec.EvChannel()
        if not hasattr(ch, "open_random_access"):
            return ("prad2py is missing open_random_access — rebuild prad2py "
                    "after the EvChannel.cpp changes:\n"
                    "  cmake -DBUILD_PYTHON=ON -S . -B build && "
                    "cmake --build build --target prad2py")
    except Exception as e:
        return f"{type(e).__name__}: {e}"
    return None


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
#  Histogram accumulator
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

    def reset(self):
        self.bins[:] = 0
        self.under = 0
        self.over = 0

    def to_json(self) -> Dict:
        return {"bins": self.bins.tolist(),
                "underflow": int(self.under), "overflow": int(self.over)}


@dataclass
class ChannelHists:
    roc:         int
    slot:        int
    channel:     int
    module:      Optional[str]
    events:      int = 0
    peak_events: int = 0
    height:      Optional[Hist1D] = None
    integral:    Optional[Hist1D] = None
    position:    Optional[Hist1D] = None
    npeaks:      Optional[Hist1D] = None


# ===========================================================================
#  Config / map loaders
# ===========================================================================

def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    """(crate, slot, channel) -> module_name."""
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
#  Hist bin builder helpers
# ===========================================================================

def _nbins(c: Dict) -> int:
    span = c["max"] - c["min"]
    return max(1, int(np.ceil(span / c["step"])))


def _make_hists(h_cfg: Dict, i_cfg: Dict, p_cfg: Dict, n_cfg: Dict,
                roc: int, slot: int, channel: int,
                module: Optional[str]) -> ChannelHists:
    return ChannelHists(
        roc=roc, slot=slot, channel=channel, module=module,
        height  =Hist1D(_nbins(h_cfg), h_cfg["min"], h_cfg["step"]),
        integral=Hist1D(_nbins(i_cfg), i_cfg["min"], i_cfg["step"]),
        position=Hist1D(_nbins(p_cfg), p_cfg["min"], p_cfg["step"]),
        npeaks  =Hist1D(_nbins(n_cfg), n_cfg["min"], n_cfg["step"]),
    )


# ===========================================================================
#  Indexer — background pass to locate all physics sub-events
# ===========================================================================

class IndexerWorker(QObject):
    """Scans the file once in RA mode to record (evio_idx, sub_idx) per
    physics sub-event.  No waveform decoding — Scan() only."""

    progressed = pyqtSignal(int, int)   # (evio_events_scanned, total_evio_events)
    finished   = pyqtSignal(object)     # {"path": str, "index": list, "total_evio": int}
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str, daq_config_path: str):
        super().__init__()
        self._path = evio_path
        self._daq_cfg_path = daq_config_path
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            self.finished.emit(self._run())
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> Dict:
        dec = prad2py.dec
        cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
               else dec.load_daq_config())
        ch  = dec.EvChannel()
        ch.set_config(cfg)
        st = ch.open_random_access(self._path)
        if st != dec.Status.success:
            raise RuntimeError(f"cannot open {self._path} in RA mode: {st}")

        total_evio = ch.get_random_access_event_count()
        index: List[Tuple[int, int]] = []

        progress_every = max(1, total_evio // 200)
        for ei in range(total_evio):
            if self._cancel:
                break
            if ch.read_event_by_index(ei) != dec.Status.success:
                continue
            if not ch.scan():
                continue
            if ch.get_event_type() != dec.EventType.Physics:
                continue
            n_sub = ch.get_n_events()
            for si in range(n_sub):
                index.append((ei, si))
            if (ei % progress_every) == 0:
                self.progressed.emit(ei + 1, total_evio)

        ch.close()
        self.progressed.emit(total_evio, total_evio)
        return {"path": self._path, "index": index,
                "total_evio": total_evio, "cancelled": self._cancel}


# ===========================================================================
#  Batch processor — fills the current module's hists for the next N events
# ===========================================================================

class BatchWorker(QObject):
    """Reads events start_idx .. start_idx + n - 1 (no display updates) and
    fills histograms.  Runs in its own thread with its own EvChannel handle
    (separate from the UI's).  If ``accum_all`` is True every channel found
    in each event is analysed; otherwise only ``target_key``."""

    progressed = pyqtSignal(int, int, int)  # (done, target, peaks_found)
    finished   = pyqtSignal(int)            # events_processed
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str, daq_config_path: str,
                 index: List[Tuple[int, int]],
                 start_idx: int, count: int,
                 target_key: Tuple[int, int, int],
                 channels: Dict[Tuple[int, int, int], ChannelHists],
                 wcfg: WaveConfig, hist_threshold: float,
                 accept_mask: int, reject_mask: int,
                 accum_all: bool,
                 accumulated: Optional[np.ndarray] = None):
        super().__init__()
        self._path = evio_path
        self._daq_cfg_path = daq_config_path
        self._index = index
        self._start = start_idx
        self._count = count
        self._target_key = target_key
        self._channels = channels
        self._wcfg = wcfg
        self._thr = hist_threshold
        self._accept = accept_mask
        self._reject = reject_mask
        self._accum_all = accum_all
        self._accumulated = accumulated     # shared bool array, mutated in place
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            self.finished.emit(self._run())
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> int:
        dec = prad2py.dec
        cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
               else dec.load_daq_config())
        ch  = dec.EvChannel()
        ch.set_config(cfg)
        st = ch.open_random_access(self._path)
        if st != dec.Status.success:
            raise RuntimeError(f"cannot open {self._path} in RA mode: {st}")

        n_done = 0
        peaks_found = 0
        progress_every = max(1, self._count // 100)
        wcfg = self._wcfg
        thr = self._thr
        channels = self._channels
        accum_all = self._accum_all
        tgt_key = self._target_key

        try:
            for i in range(self._count):
                if self._cancel:
                    break
                phys_idx = self._start + i
                if phys_idx >= len(self._index):
                    break
                # Dedup: skip events already folded in by browsing.
                if (self._accumulated is not None
                        and 0 <= phys_idx < self._accumulated.size
                        and self._accumulated[phys_idx]):
                    n_done += 1
                    continue
                ev_idx, sub_idx = self._index[phys_idx]
                if ch.read_event_by_index(ev_idx) != dec.Status.success:
                    continue
                if not ch.scan():
                    continue
                ch.select_event(sub_idx)
                info = ch.info()
                tb = int(info.trigger_bits)
                if self._accept and (tb & self._accept) == 0:
                    n_done += 1
                    continue
                if self._reject and (tb & self._reject):
                    n_done += 1
                    continue

                fadc_evt = ch.fadc()
                for r in range(fadc_evt.nrocs):
                    roc = fadc_evt.roc(r)
                    roc_tag = int(roc.tag)
                    for s in roc.present_slots():
                        slot = roc.slot(s)
                        for c in slot.present_channels():
                            key = (roc_tag, s, c)
                            if not accum_all and key != tgt_key:
                                continue
                            hits = channels.get(key)
                            if hits is None:
                                continue
                            samples = slot.channel(c).samples
                            if samples.size < 10:
                                continue
                            _, _, peaks = analyze(samples, wcfg)
                            np_kept = 0
                            for p in peaks:
                                if p.height >= thr:
                                    hits.height.fill(p.height)
                                    hits.integral.fill(p.integral)
                                    hits.position.fill(p.time)
                                    np_kept += 1
                                    peaks_found += 1
                            hits.npeaks.fill(np_kept)
                            hits.events += 1
                            if np_kept > 0:
                                hits.peak_events += 1
                if self._accumulated is not None and 0 <= phys_idx < self._accumulated.size:
                    self._accumulated[phys_idx] = True
                n_done += 1

                if (i % progress_every) == 0:
                    self.progressed.emit(n_done, self._count, peaks_found)
        finally:
            ch.close()

        self.progressed.emit(n_done, self._count, peaks_found)
        return n_done


def _find_channel_samples(fadc_evt, roc_tag: int, slot: int, channel: int):
    """Return samples array for (roc, slot, ch), or None if not present."""
    for r in range(fadc_evt.nrocs):
        roc = fadc_evt.roc(r)
        if int(roc.tag) != roc_tag:
            continue
        if slot not in roc.present_slots():
            continue
        slot_data = roc.slot(slot)
        if channel not in slot_data.present_channels():
            continue
        return slot_data.channel(channel).samples
    return None


# ===========================================================================
#  Small themed overlay controls for plot widgets
# ===========================================================================

def _overlay_checkbox_qss() -> str:
    """QSS for a compact checkbox drawn on top of a plot canvas."""
    return (
        f"QCheckBox{{color:{THEME.TEXT_DIM};background:{THEME.PANEL};"
        f"padding:2px 6px;border:1px solid {THEME.BORDER};border-radius:6px;}}"
        f"QCheckBox:hover{{color:{THEME.TEXT};"
        f"border:1px solid {THEME.ACCENT};}}"
        f"QCheckBox:checked{{color:{THEME.TEXT};}}"
        f"QCheckBox::indicator{{width:12px;height:12px;"
        f"border:1px solid {THEME.BORDER};border-radius:3px;"
        f"background:{THEME.BG};}}"
        f"QCheckBox::indicator:hover{{border:1px solid {THEME.ACCENT};}}"
        f"QCheckBox::indicator:checked{{background:{THEME.ACCENT};"
        f"border:1px solid {THEME.ACCENT};}}"
    )


def _overlay_button_qss() -> str:
    """QSS for a compact pushbutton drawn on top of a plot canvas."""
    return (
        f"QPushButton{{color:{THEME.TEXT_DIM};background:{THEME.PANEL};"
        f"padding:2px 8px;border:1px solid {THEME.BORDER};border-radius:6px;"
        f"font:bold 9pt Monospace;}}"
        f"QPushButton:hover{{color:{THEME.TEXT};"
        f"border:1px solid {THEME.ACCENT};}}"
        f"QPushButton:disabled{{color:{THEME.TEXT_MUTED};}}"
    )


# ===========================================================================
#  Hist1DWidget — QPainter bar chart with optional log Y
# ===========================================================================

class Hist1DWidget(QWidget):
    PAD_L, PAD_R, PAD_T, PAD_B = 58, 14, 20, 20

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(140)
        self._bins: np.ndarray = np.zeros(0, dtype=np.int64)
        self._bmin: float = 0.0
        self._bstep: float = 1.0
        self._under: int = 0
        self._over: int = 0
        self._title: str = ""
        self._xlabel: str = ""
        self._color: QColor = QColor(THEME.ACCENT)
        self._log_y: bool = False
        self._hover_idx: int = -1

        self._logy_cb = QCheckBox("log Y", self)
        self._logy_cb.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._logy_cb.setStyleSheet(_overlay_checkbox_qss())
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
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = self._plot_rect()
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        if self._title:
            f = QFont("Monospace", 10); f.setBold(True)
            p.setFont(f)
            p.setPen(QColor(THEME.TEXT))
            p.drawText(int(r.left()), int(r.top() - 8), self._title)

        n = self._bins.size
        if n == 0 or self._bins.sum() == 0:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "(no data)")
            return

        if self._log_y:
            vals = np.where(self._bins > 0,
                            np.log10(self._bins.astype(np.float64)), 0.0)
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
            p.setPen(QPen(QColor(THEME.SELECT_BORDER), 1.2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(QRectF(x0, r.top(), max(bar_w, 1.0), r.height()))

        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            if self._log_y:
                val = 10 ** (ymin + frac * (ymax - ymin))
            else:
                val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4),
                       _fmt_count(val))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            x = r.left() + frac * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            val = self._bmin + frac * n * self._bstep
            p.drawText(int(x - 24), int(r.bottom() + 14), f"{val:g}")

        entries = int(self._bins.sum())
        info = f"N={entries:,}"
        if self._under:
            info += f"  under={self._under:,}"
        if self._over:
            info += f"  over={self._over:,}"
        p.setFont(QFont("Monospace", 9))
        p.setPen(QColor(THEME.TEXT_DIM))
        info_rect = QRectF(r.left(), r.top() - 20,
                           max(1.0, self.width() - self._logy_cb.width()
                               - 20 - r.left()),
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
#  WaveformPlotWidget — draws the current event's raw FADC samples
# ===========================================================================

class WaveformPlotWidget(QWidget):
    PAD_L, PAD_R, PAD_T, PAD_B = 52, 14, 22, 30
    MAX_STACK = 200

    # Same peak colour palette the web frontend uses (resources/viewer.js PC).
    _PEAK_PALETTE = (
        "#00b4d8", "#ff6b6b", "#51cf66", "#ffd43b",
        "#cc5de8", "#ff922b", "#20c997", "#f06595",
    )

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(150)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self._samples: np.ndarray = np.zeros(0, dtype=np.float32)
        self._peaks: List[Peak] = []
        self._ped_mean: float = 0.0
        self._ped_rms: float = 0.0
        self._title: str = ""
        self._clk_mhz: float = 250.0

        # --- stack mode state ---
        self._stack_enabled: bool = False
        self._stack_traces: List[np.ndarray] = []   # bounded by MAX_STACK
        self._stack_key: str = ""                   # (module,roc,slot,ch) key

        # --- overlay controls (top-right) ---
        self._stack_cb = QCheckBox("Stack", self)
        self._stack_cb.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._stack_cb.setStyleSheet(_overlay_checkbox_qss())
        self._stack_cb.setToolTip(
            f"Overlay waveforms across events (up to {self.MAX_STACK}). "
            "Peaks and integral shading hidden in stack mode.")
        self._stack_cb.toggled.connect(self._on_stack_toggled)
        self._stack_cb.adjustSize()

        self._stack_clear_btn = QPushButton("Clear", self)
        self._stack_clear_btn.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._stack_clear_btn.setStyleSheet(_overlay_button_qss())
        self._stack_clear_btn.setToolTip("Drop all stacked waveforms")
        self._stack_clear_btn.clicked.connect(self.clear_stack)
        self._stack_clear_btn.setVisible(False)
        self._stack_clear_btn.adjustSize()

        self._stack_count_lbl = QLabel("", self)
        self._stack_count_lbl.setFont(QFont("Monospace", 9))
        self._stack_count_lbl.setStyleSheet(
            f"color:{THEME.TEXT_DIM};background:transparent;")
        self._stack_count_lbl.setVisible(False)
        self._stack_count_lbl.adjustSize()

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def set_data(self, samples: np.ndarray, peaks: List[Peak],
                 ped_mean: float, ped_rms: float,
                 title: str, clk_mhz: float = 250.0,
                 stack_key: Optional[str] = None):
        samples = np.asarray(samples, dtype=np.float32)
        self._samples = samples
        self._peaks = peaks
        self._ped_mean = ped_mean
        self._ped_rms  = ped_rms
        self._title = title
        self._clk_mhz = clk_mhz

        if self._stack_enabled:
            key = stack_key if stack_key is not None else title
            if key != self._stack_key:
                self._stack_traces = []
                self._stack_key = key
            if samples.size >= 2:
                self._stack_traces.append(samples.copy())
                if len(self._stack_traces) > self.MAX_STACK:
                    self._stack_traces = self._stack_traces[-self.MAX_STACK:]
            self._update_stack_counter()
        self.update()

    def clear(self, title: str = ""):
        self._samples = np.zeros(0, dtype=np.float32)
        self._peaks = []
        self._title = title
        self.update()

    def clear_stack(self):
        """Drop every accumulated trace but keep the current waveform."""
        self._stack_traces = []
        self._stack_key = ""
        self._update_stack_counter()
        self.update()

    def is_stacking(self) -> bool:
        return self._stack_enabled

    # ------------------------------------------------------------------
    #  Internals
    # ------------------------------------------------------------------

    def _on_stack_toggled(self, on: bool):
        self._stack_clear_btn.setVisible(on)
        self._stack_count_lbl.setVisible(on)
        if not on:
            self._stack_traces = []
            self._stack_key = ""
        self._update_stack_counter()
        self._layout_overlays()
        self.update()

    def _update_stack_counter(self):
        self._stack_count_lbl.setText(
            f"{len(self._stack_traces)}/{self.MAX_STACK}")
        self._stack_count_lbl.adjustSize()
        self._layout_overlays()

    def _layout_overlays(self):
        # top-right: [count]  [Clear]  [Stack]
        margin = 6
        x = self.width() - margin
        y = 4
        x -= self._stack_cb.width()
        self._stack_cb.move(x, y)
        if self._stack_clear_btn.isVisible():
            x -= self._stack_clear_btn.width() + 4
            self._stack_clear_btn.move(x, y)
        if self._stack_count_lbl.isVisible():
            x -= self._stack_count_lbl.width() + 6
            self._stack_count_lbl.move(x, y + 2)

    def _plot_rect(self) -> QRectF:
        w, h = self.width(), self.height()
        return QRectF(self.PAD_L, self.PAD_T,
                      max(1.0, w - self.PAD_L - self.PAD_R),
                      max(1.0, h - self.PAD_T - self.PAD_B))

    def resizeEvent(self, ev):
        self._stack_cb.adjustSize()
        self._stack_clear_btn.adjustSize()
        self._stack_count_lbl.adjustSize()
        self._layout_overlays()
        super().resizeEvent(ev)

    # ------------------------------------------------------------------
    #  Painting
    # ------------------------------------------------------------------

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = self._plot_rect()
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        if self._title:
            f = QFont("Monospace", 10); f.setBold(True)
            p.setFont(f)
            p.setPen(QColor(THEME.TEXT))
            suffix = f" — Stacked ({len(self._stack_traces)})" if self._stack_enabled else ""
            p.drawText(int(r.left()), int(r.top() - 6), self._title + suffix)

        if self._stack_enabled:
            self._paint_stacked(p, r)
        else:
            self._paint_single(p, r)

    # --- single-event (default) view ---------------------------------

    def _paint_single(self, p: QPainter, r: QRectF):
        n = self._samples.size
        if n < 2:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(no waveform — click Next to load an event)")
            return

        ymin = float(self._samples.min())
        ymax = float(self._samples.max())
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad_y = (ymax - ymin) * 0.05
        ymin -= pad_y; ymax += pad_y

        def to_sx(i: float) -> float:
            return r.left() + (i / (n - 1)) * r.width()

        def to_sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        # pedestal baseline
        y_ped = to_sy(self._ped_mean) if self._ped_mean != 0 else None
        if y_ped is not None:
            p.setPen(QPen(QColor(THEME.TEXT_DIM), 1, Qt.PenStyle.DashLine))
            p.drawLine(int(r.left()), int(y_ped), int(r.right()), int(y_ped))
            # threshold line (same formula as waveform.js: pm + max(5*pr, 3))
            thr_v = self._ped_mean + max(5.0 * self._ped_rms, 3.0)
            y_thr = to_sy(thr_v)
            p.setPen(QPen(QColor(THEME.TEXT_MUTED), 1, Qt.PenStyle.DotLine))
            p.drawLine(int(r.left()), int(y_thr), int(r.right()), int(y_thr))

        # Fill the integral area (between pedestal and waveform) per peak,
        # colour-coded from _PEAK_PALETTE. Mirrors resources/waveform.js.
        if self._peaks and y_ped is not None:
            for i, pk in enumerate(self._peaks):
                base = QColor(self._PEAK_PALETTE[i % len(self._PEAK_PALETTE)])
                fill = QColor(base); fill.setAlphaF(0.18)
                poly = QPolygonF()
                j = max(0, int(pk.left))
                j_end = min(n - 1, int(pk.right))
                for k in range(j, j_end + 1):
                    poly.append(QPointF(to_sx(k),
                                        to_sy(float(self._samples[k]))))
                # close along the pedestal baseline
                poly.append(QPointF(to_sx(j_end), y_ped))
                poly.append(QPointF(to_sx(j), y_ped))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(fill)
                p.drawPolygon(poly)
                # outline the peak section with the solid palette colour
                p.setPen(QPen(base, 2))
                p.setBrush(Qt.BrushStyle.NoBrush)
                for k in range(j, j_end):
                    p.drawLine(QPointF(to_sx(k),
                                       to_sy(float(self._samples[k]))),
                               QPointF(to_sx(k + 1),
                                       to_sy(float(self._samples[k + 1]))))

        # waveform line (default accent, drawn under peak outlines)
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        for i in range(n - 1):
            p.drawLine(int(to_sx(i)),     int(to_sy(float(self._samples[i]))),
                       int(to_sx(i + 1)), int(to_sy(float(self._samples[i + 1]))))

        # peak markers (diamonds, coloured per peak)
        if self._peaks:
            for i, pk in enumerate(self._peaks):
                if pk.pos < 0 or pk.pos >= n:
                    continue
                col = QColor(self._PEAK_PALETTE[i % len(self._PEAK_PALETTE)])
                p.setPen(QPen(col, 1.2))
                p.setBrush(col)
                cx = to_sx(pk.pos)
                cy = to_sy(float(self._samples[pk.pos]))
                diamond = QPolygonF([
                    QPointF(cx,     cy - 4),
                    QPointF(cx + 4, cy),
                    QPointF(cx,     cy + 4),
                    QPointF(cx - 4, cy),
                ])
                p.drawPolygon(diamond)

        self._paint_axes(p, r, ymin, ymax)

        info = (f"ped={self._ped_mean:.1f}  rms={self._ped_rms:.2f}  "
                f"peaks={len(self._peaks)}")
        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(QRectF(r.left(), r.top() - 20,
                          max(1.0, r.width() - 8), 14),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   info)

    # --- stacked overlay view -----------------------------------------

    def _paint_stacked(self, p: QPainter, r: QRectF):
        traces = self._stack_traces
        if not traces:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(stack is empty — step through events to accumulate)")
            return

        # Compute common y-range across all traces.
        ymin = min(float(w.min()) for w in traces)
        ymax = max(float(w.max()) for w in traces)
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad_y = (ymax - ymin) * 0.05
        ymin -= pad_y; ymax += pad_y

        # Width uses the max length so shorter traces still fit left-aligned.
        n_max = max(w.size for w in traces)

        def to_sx(i: float, n: int) -> float:
            return r.left() + (i / max(1, n - 1)) * r.width()

        def to_sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        # Dimmed stacked traces.
        dim = QColor(THEME.ACCENT); dim.setAlphaF(0.18)
        p.setPen(QPen(dim, 1))
        for w in traces[:-1]:
            n = w.size
            for i in range(n - 1):
                p.drawLine(int(to_sx(i, n)),     int(to_sy(float(w[i]))),
                           int(to_sx(i + 1, n)), int(to_sy(float(w[i + 1]))))

        # Latest trace drawn on top at full colour.
        latest = traces[-1]
        n = latest.size
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        for i in range(n - 1):
            p.drawLine(int(to_sx(i, n)),     int(to_sy(float(latest[i]))),
                       int(to_sx(i + 1, n)), int(to_sy(float(latest[i + 1]))))

        self._paint_axes(p, r, ymin, ymax, n=n_max)

        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(QRectF(r.left(), r.top() - 20,
                          max(1.0, r.width() - 8), 14),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   f"stack={len(traces)}/{self.MAX_STACK}")

    # --- shared axis/tick drawing -------------------------------------

    def _paint_axes(self, p: QPainter, r: QRectF,
                    ymin: float, ymax: float, n: Optional[int] = None):
        if n is None:
            n = self._samples.size
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4), f"{val:.0f}")
        tick_every = max(1, n // 8)
        for i in range(0, n, tick_every):
            x = r.left() + (i / max(1, n - 1)) * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            ns = i * 1e3 / self._clk_mhz
            p.drawText(int(x - 18), int(r.bottom() + 14), f"{ns:g}")

        p.setFont(QFont("Monospace", 9))
        p.drawText(int(r.left() + r.width() / 2 - 10),
                   int(r.bottom() + 26), "ns")


# ===========================================================================
#  WaveformGeoView — small HyCal overview for module selection
# ===========================================================================

# Module types that should get a small name label painted on top of the cell
# (so the tiny LMS / V blocks off to the left of HyCal are identifiable).
_LABEL_TYPES = {"LMS", "SCINT"}

class WaveformGeoView(HyCalMapWidget):
    """Compact HyCal geo view with two colour-coding modes.

    * ``current``  — module colour = max peak integral in the current event.
    * ``overall``  — module colour = occupancy (events-with-peak / accumulated
      events) across all events the user has browsed / batched.

    Modules that have never been seen in an event are drawn in a flat grey
    so "no-data" stays visually distinct from "data, low value".  Clicking
    any module emits moduleClicked with its name.
    """

    MODE_CURRENT = "current"
    MODE_OVERALL = "overall"

    # Resolved at paint time so the active theme wins; see :class:`THEME`.
    @property
    def UNAVAIL_COLOR(self) -> QColor:
        return QColor(THEME.BORDER)

    @property
    def SELECT_COLOR(self) -> QColor:
        return QColor(THEME.SELECT_BORDER)

    def __init__(self, parent=None):
        super().__init__(parent, show_colorbar=True, include_lms=True,
                         margin_top=4, margin_bottom=38,
                         min_size=(220, 260), shrink=0.90)
        self._available: set = set()
        self._selected_name: Optional[str] = None
        self._label_names: set = set()          # filled in set_modules()
        self._mode = self.MODE_CURRENT
        self._current_vals: Dict[str, float] = {}
        self._overall_vals: Dict[str, float] = {}

        # Top-left mode toggle.  Default label matches MODE_CURRENT.
        self._mode_btn = QPushButton("Current", self)
        self._mode_btn.setFixedSize(74, 22)
        _f = QFont("Consolas", 9); _f.setBold(True)
        self._mode_btn.setFont(_f)
        self._mode_btn.setToolTip(
            "Colour coding:\n"
            "  Current — max peak integral in the currently viewed event\n"
            "  Overall — occupancy (events-with-peak / accumulated events)")
        self._mode_btn.setStyleSheet(themed(
            "QPushButton{background:rgba(29,29,31,220);color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:4px;}"
            "QPushButton:hover{background:#28282a;color:#e6edf3;}"))
        self._mode_btn.clicked.connect(self._toggle_mode)

    def set_modules(self, modules):
        super().set_modules(modules)
        self._label_names = {m.name for m in self._modules
                             if m.mod_type in _LABEL_TYPES}

    def set_available(self, names):
        self._available = set(names)
        self.update()

    def set_selected_module(self, name: Optional[str]):
        if name != self._selected_name:
            self._selected_name = name
            self.update()

    def set_current_values(self, vals: Dict[str, float]):
        self._current_vals = vals
        if self._mode == self.MODE_CURRENT:
            self._apply_mode_values()

    def set_overall_values(self, vals: Dict[str, float]):
        self._overall_vals = vals
        if self._mode == self.MODE_OVERALL:
            self._apply_mode_values()

    def _toggle_mode(self):
        self._mode = (self.MODE_OVERALL if self._mode == self.MODE_CURRENT
                      else self.MODE_CURRENT)
        self._mode_btn.setText(
            "Overall" if self._mode == self.MODE_OVERALL else "Current")
        self._apply_mode_values()

    def _apply_mode_values(self):
        if self._mode == self.MODE_OVERALL:
            self.set_values(self._overall_vals)
            self.set_range(0.0, 1.0)          # occupancy fraction
        else:
            self.set_values(self._current_vals)
            self.auto_range()                  # max-integral per-event rescales

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._mode_btn.move(6, 6)

    def _colorbar_center_text(self) -> str:
        return ("occupancy" if self._mode == self.MODE_OVERALL
                else "max peak integral")

    def _paint_modules(self, p):
        # Keep the "not yet seen" grey distinct from the colormap's low end
        # so unused modules don't masquerade as "low value".
        avail = self._available
        u_col = self.UNAVAIL_COLOR
        no_data = self.NO_DATA_COLOR
        stops = self.palette_stops()
        vmin, vmax = self._vmin, self._vmax
        vals = self._values
        for name, rect in self._rects.items():
            if name not in avail:
                p.fillRect(rect, u_col)
                continue
            v = vals.get(name)
            if v is None:
                p.fillRect(rect, no_data)
                continue
            t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
            p.fillRect(rect, cmap_qcolor(t, stops))

    def _paint_overlays(self, p, w, h):
        p.setPen(QColor(THEME.TEXT))
        p.setFont(QFont("Monospace", 7, QFont.Weight.Bold))
        for name in self._label_names:
            r = self._rects.get(name)
            if r is not None:
                p.drawText(r, Qt.AlignmentFlag.AlignCenter, name)

        if self._selected_name and self._selected_name in self._rects:
            p.setPen(QPen(self.SELECT_COLOR, 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._selected_name])
        super()._paint_overlays(p, w, h)   # hover border

    def _tooltip_text(self, name: str) -> str:
        if name not in self._available:
            return f"{name}  (not seen yet)"
        v = self._values.get(name)
        unit = ("occupancy" if self._mode == self.MODE_OVERALL
                else "max integral")
        if v is None:
            return f"{name}  ({unit}: —)"
        return f"{name}  {unit}={v:.3g}"


# ===========================================================================
#  Main window
# ===========================================================================

_NATKEY_RE = re.compile(r"(\d+)")
def _natural_sort_key(s: str):
    return [int(p) if p.isdigit() else p.lower()
            for p in _NATKEY_RE.split(s or "")]


class WaveformViewerWindow(QMainWindow):

    def __init__(self,
                 *,
                 hist_config: Dict,
                 daq_map: Dict,
                 roc_to_crate: Dict,
                 accept_mask: int,
                 reject_mask: int,
                 daq_config_path: str,
                 hycal_modules: Optional[List] = None):
        super().__init__()
        self._hist_config   = hist_config
        self._daq_map       = daq_map
        self._roc_to_crate  = roc_to_crate
        self._accept_mask   = accept_mask
        self._reject_mask   = reject_mask
        self._daq_cfg_path  = daq_config_path
        self._hycal_modules = hycal_modules or []

        # Bin configs — merge user config with defaults, add n-peaks hist
        self._h_cfg = hist_config.get("height_hist",
                                      {"min": 0, "max": 4000,  "step": 10})
        self._i_cfg = hist_config.get("integral_hist",
                                      {"min": 0, "max": 20000, "step": 100})
        self._p_cfg = hist_config.get("time_hist",
                                      {"min": 0, "max": 400,   "step": 4})
        self._n_cfg = {"min": 0, "max": 10, "step": 1}

        thr_cfg = hist_config.get("thresholds", {})
        self._hist_threshold = float(thr_cfg.get("min_peak_height", 10.0))
        self._wcfg = WaveConfig()
        self._wcfg.min_peak_ratio = float(thr_cfg.get(
            "min_secondary_peak_ratio", self._wcfg.min_peak_ratio))

        # File state
        self._evio_path: Optional[Path] = None
        self._index: List[Tuple[int, int]] = []
        self._current_idx: int = -1
        # One bool per physics sub-event — True once its peaks have been
        # folded into self._channels hists.  Lets Prev/Next re-display an
        # event without double-counting.  np.bool = 1 byte/event, so 1 M
        # events ≈ 1 MB, 10 M ≈ 10 MB — negligible.
        self._accumulated: Optional[np.ndarray] = None

        # Per-channel accumulated hists, keyed by (roc, slot, ch)
        self._channels: Dict[Tuple[int, int, int], ChannelHists] = {}
        self._selected_key: Optional[Tuple[int, int, int]] = None

        # Reader state — open EvChannel in RA mode, kept alive across browse
        self._ch: Optional["prad2py.dec.EvChannel"] = None
        self._reader_path: Optional[str] = None

        # Worker threads
        self._idx_worker: Optional[IndexerWorker] = None
        self._idx_thread: Optional[QThread] = None
        self._batch_worker: Optional[BatchWorker] = None
        self._batch_thread: Optional[QThread] = None

        apply_theme_palette(self)
        self._build_ui()
        self._make_menu()

    # -- UI --

    def _build_ui(self):
        self.setWindowTitle("Waveform Viewer")
        self.resize(1500, 1000)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(3)

        self._file_lbl = QLabel("(no file loaded)")
        self._file_lbl.setFont(QFont("Monospace", 10))
        self._file_lbl.setStyleSheet(themed("color:#8b949e;"))
        root.addWidget(self._file_lbl)

        # -- top control bar: navigation on the left, module picker on the right --
        top = QHBoxLayout()

        self._prev_btn = self._small_btn("◀ Prev", self._on_prev)
        self._next_btn = self._small_btn("Next ▶", self._on_next)
        self._prev_btn.setEnabled(False)
        self._next_btn.setEnabled(False)
        top.addWidget(self._prev_btn)
        top.addWidget(self._next_btn)

        top.addSpacing(12)
        top.addWidget(self._mk_label("Event:"))
        self._event_spin = QSpinBox()
        self._event_spin.setFont(QFont("Monospace", 10))
        self._event_spin.setMinimum(0)
        self._event_spin.setMaximum(0)
        self._event_spin.setStyleSheet(themed(
            "QSpinBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:6px;padding:2px 6px;}"))
        self._event_spin.editingFinished.connect(self._on_spin_jump)
        self._event_spin.setEnabled(False)
        top.addWidget(self._event_spin)
        self._total_lbl = QLabel(" / 0")
        self._total_lbl.setFont(QFont("Monospace", 10))
        self._total_lbl.setStyleSheet(themed("color:#8b949e;"))
        top.addWidget(self._total_lbl)

        top.addSpacing(18)
        self._batch_btn = self._small_btn("Process next 10k",
                                          self._on_batch_10k, primary=True)
        self._batch_btn.setEnabled(False)
        top.addWidget(self._batch_btn)
        self._cancel_batch_btn = self._small_btn("Cancel",
                                                 self._on_cancel_batch)
        self._cancel_batch_btn.setVisible(False)
        top.addWidget(self._cancel_batch_btn)

        self._accum_all_cb = QCheckBox("Accumulate all modules")
        self._accum_all_cb.setChecked(True)
        self._accum_all_cb.setToolTip(
            "On: browsing fills histograms for every channel present in the event "
            "(slow ~1-5 s per Next, analysis runs per channel).\n"
            "Off: only the selected channel's hist accumulates (instant browse).")
        self._accum_all_cb.setStyleSheet(themed(
            "QCheckBox{color:#c9d1d9;font:10pt Monospace;}"
            "QCheckBox:hover{color:#e6edf3;}"))
        top.addSpacing(8)
        top.addWidget(self._accum_all_cb)

        self._batch_status = QLabel("")
        self._batch_status.setFont(QFont("Monospace", 10))
        self._batch_status.setStyleSheet(themed("color:#8b949e;"))
        top.addSpacing(8)
        top.addWidget(self._batch_status)

        top.addStretch(1)

        # Right cluster: module dropdown + reset hist.
        mod_lbl = QLabel("Module:")
        mod_lbl.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
        mod_lbl.setStyleSheet(themed("color:#c9d1d9;"))
        top.addWidget(mod_lbl)
        self._combo = QComboBox()
        self._combo.setEditable(True)
        self._combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._combo.setFont(QFont("Monospace", 11))
        self._combo.setMinimumContentsLength(32)
        self._combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:6px;padding:2px 6px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "selection-background-color:#1f6feb;}"))
        comp = self._combo.completer()
        if comp is not None:
            comp.setFilterMode(Qt.MatchFlag.MatchContains)
            comp.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            comp.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        self._combo.currentIndexChanged.connect(self._on_combo_changed)
        top.addWidget(self._combo)
        self._reset_btn = self._small_btn("Reset hist", self._reset_current_hists)
        self._reset_btn.setEnabled(False)
        top.addWidget(self._reset_btn)

        root.addLayout(top)

        self._info = QLabel("")
        self._info.setFont(QFont("Monospace", 10))
        self._info.setStyleSheet(themed("color:#8b949e;"))
        root.addWidget(self._info)

        # -- main split: geo+waveform on the left, 2x2 hists on the right --
        split = QSplitter(Qt.Orientation.Horizontal)

        # Left: geo view (square, top) + waveform plot (bottom)
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        left_lay.setSpacing(4)
        self._geo = WaveformGeoView()
        if self._hycal_modules:
            self._geo.set_modules(self._hycal_modules)
        self._geo.moduleClicked.connect(self._on_geo_clicked)
        self._geo.setSizePolicy(QSizePolicy.Policy.Expanding,
                                QSizePolicy.Policy.Expanding)
        # Geo view is square, so give it more vertical space; the waveform
        # plot is wide-and-short and fills whatever's left below.
        left_lay.addWidget(self._geo, stretch=3)
        self._wave = WaveformPlotWidget()
        left_lay.addWidget(self._wave, stretch=1)
        split.addWidget(left)

        # Right: four histograms stacked vertically (each wide, long in x).
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(0)
        self._h_height   = Hist1DWidget()
        self._h_integral = Hist1DWidget()
        self._h_position = Hist1DWidget()
        self._h_npeaks   = Hist1DWidget()
        for hist in (self._h_height, self._h_integral,
                     self._h_position, self._h_npeaks):
            right_lay.addWidget(hist, stretch=1)
        split.addWidget(right)

        # Even 50/50 split between the geo+waveform column and the hist stack.
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setSizes([750, 750])
        root.addWidget(split, stretch=1)

        self.setStatusBar(QStatusBar())
        self._clear_plots()

        # Keyboard: ← / → to navigate prev / next.
        QShortcut(QKeySequence(Qt.Key.Key_Left),  self, activated=self._on_prev)
        QShortcut(QKeySequence(Qt.Key.Key_Right), self, activated=self._on_next)

    def _small_btn(self, text: str, slot, primary: bool = False) -> QPushButton:
        btn = QPushButton(text)
        bg = "#1f6feb" if primary else "#21262d"
        fg = "#ffffff" if primary else "#c9d1d9"
        btn.setStyleSheet(themed(
            f"QPushButton{{background:{bg};color:{fg};"
            f"border:1px solid #30363d;padding:5px 14px;"
            f"font:bold 10pt Monospace;border-radius:3px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{background:#161b22;color:#484f58;}}"))
        btn.clicked.connect(slot)
        return btn

    def _mk_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 10))
        lbl.setStyleSheet(themed("color:#c9d1d9;"))
        return lbl

    def _make_menu(self):
        mb = self.menuBar()
        mf = mb.addMenu("&File")

        a_open = QAction("Open &evio…", self)
        a_open.setShortcut("Ctrl+O")
        a_open.triggered.connect(self._open_evio_dialog)
        mf.addAction(a_open)

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

    # -- file open --

    def _open_evio_dialog(self):
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open evio file", str(Path.cwd()),
            "evio files (*.evio *.evio.*);;All files (*)")
        if path_str:
            self.open_path(Path(path_str))

    def open_path(self, path: Path):
        err = _check_random_access_support()
        if err:
            QMessageBox.critical(self, "prad2py issue", err)
            return
        if self._idx_thread is not None:
            QMessageBox.information(self, "Busy", "Already indexing.")
            return

        # Tear down any previous reader / hists
        self._close_reader()
        self._channels.clear()
        self._selected_key = None
        self._combo.blockSignals(True); self._combo.clear(); self._combo.blockSignals(False)
        self._index = []
        self._current_idx = -1
        self._accumulated = None
        self._clear_plots()

        # Start indexer
        dlg = QProgressDialog(f"Indexing {path.name} …", "Cancel", 0, 100, self)
        dlg.setWindowTitle("Indexing")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()

        worker = IndexerWorker(str(path), self._daq_cfg_path)
        thread = QThread(self)
        worker.moveToThread(thread)

        def _on_progress(done: int, total: int):
            if total > 0:
                dlg.setMaximum(total)
                dlg.setValue(done)
            dlg.setLabelText(f"Indexing {path.name}\n"
                             f"evio events: {done:,} / {total:,}")

        thread.started.connect(worker.run)
        worker.progressed.connect(_on_progress)
        worker.finished.connect(lambda res: self._on_index_done(path, res))
        worker.failed.connect(lambda msg: self._on_index_failed(path, msg))
        dlg.canceled.connect(worker.request_cancel)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(dlg.close)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(lambda: setattr(self, "_idx_thread", None))

        self._idx_worker = worker
        self._idx_thread = thread
        thread.start()

    def _on_index_done(self, path: Path, res: Dict):
        self._evio_path = path
        self._index = res["index"]
        n_phys = len(self._index)
        cancelled = bool(res.get("cancelled"))
        self._accumulated = np.zeros(n_phys, dtype=bool)

        # Open reader handle for browse use
        ok, err = self._open_reader(str(path))
        if not ok:
            QMessageBox.critical(self, "Open failed",
                                 f"Indexing finished but reader open failed:\n{err}")
            return

        self._file_lbl.setText(
            f"{path.name}   physics events: {n_phys:,}   "
            f"(evio blocks: {res['total_evio']:,})"
            + ("   [indexing cancelled]" if cancelled else ""))
        self._info.setText(
            f"threshold={self._hist_threshold:g}   "
            f"Select a module, then click Next to start browsing.")

        self._event_spin.setMaximum(max(0, n_phys - 1))
        self._event_spin.setValue(0)
        if n_phys > 0:
            self._prev_btn.setEnabled(True)
            self._next_btn.setEnabled(True)
            self._event_spin.setEnabled(True)
            self._batch_btn.setEnabled(True)
            self._a_save.setEnabled(True)
            self._reset_btn.setEnabled(True)
        self._total_lbl.setText(f" / {max(0, n_phys - 1):,}")

        self.statusBar().showMessage(
            f"Indexed {n_phys:,} physics events from {path.name}")

    def _on_index_failed(self, path: Path, msg: str):
        QMessageBox.critical(self, "Indexing failed", f"{path}\n\n{msg}")
        self.statusBar().showMessage(f"Failed to index {path.name}")

    # -- reader (browse handle) --

    def _open_reader(self, path: str) -> Tuple[bool, str]:
        try:
            dec = prad2py.dec
            cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
                   else dec.load_daq_config())
            ch  = dec.EvChannel()
            ch.set_config(cfg)
            st = ch.open_random_access(path)
            if st != dec.Status.success:
                return False, f"status = {st}"
            self._ch = ch
            self._reader_path = path
            return True, ""
        except Exception as e:
            return False, f"{type(e).__name__}: {e}"

    def _close_reader(self):
        if self._ch is not None:
            try:
                self._ch.close()
            except Exception:
                pass
        self._ch = None
        self._reader_path = None

    # -- navigation --

    def _on_prev(self):
        if self._current_idx > 0:
            self._goto(self._current_idx - 1)
        elif self._current_idx == -1 and self._index:
            self._goto(0)

    def _on_next(self):
        if self._current_idx < len(self._index) - 1:
            self._goto(self._current_idx + 1)

    def _on_spin_jump(self):
        v = self._event_spin.value()
        if 0 <= v < len(self._index) and v != self._current_idx:
            self._goto(v)

    def _goto(self, phys_idx: int):
        if self._ch is None or not (0 <= phys_idx < len(self._index)):
            return
        ev_idx, sub_idx = self._index[phys_idx]
        dec = prad2py.dec
        st = self._ch.read_event_by_index(ev_idx)
        if st != dec.Status.success:
            self.statusBar().showMessage(
                f"read_event_by_index({ev_idx}) → {st}")
            return
        if not self._ch.scan():
            self.statusBar().showMessage(f"scan() failed at physics #{phys_idx}")
            return
        self._ch.select_event(sub_idx)
        info = self._ch.info()
        tb = int(info.trigger_bits)

        self._current_idx = phys_idx
        self._event_spin.blockSignals(True)
        self._event_spin.setValue(phys_idx)
        self._event_spin.blockSignals(False)

        # Apply trigger filter: skip updating hists but still show waveform
        trig_ok = True
        if self._accept_mask and (tb & self._accept_mask) == 0: trig_ok = False
        if self._reject_mask and (tb & self._reject_mask):       trig_ok = False

        fadc_evt = self._ch.fadc()
        self._update_channel_list_from_event(fadc_evt)
        self._accumulate_and_display(fadc_evt, info, trig_ok)

    def _update_channel_list_from_event(self, fadc_evt):
        """Add any new (roc, slot, ch) seen in this event to the combo."""
        added = False
        for r in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            if roc_tag not in self._roc_to_crate:
                continue
            crate = self._roc_to_crate[roc_tag]
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = (roc_tag, s, c)
                    if key not in self._channels:
                        module = self._daq_map.get((crate, s, c))
                        self._channels[key] = _make_hists(
                            self._h_cfg, self._i_cfg, self._p_cfg, self._n_cfg,
                            roc_tag, s, c, module)
                        added = True
        if added:
            self._refresh_combo()
            self._geo.set_available({c.module for c in self._channels.values()
                                     if c.module})

    def _refresh_combo(self):
        """Re-populate combo from discovered channels, preserving selection."""
        prev_key = self._selected_key
        items: List[Tuple[str, Tuple[int, int, int], str]] = []  # (sort_key, key, label)
        for key, ch in self._channels.items():
            mod = ch.module or "(unmapped)"
            label = (f"{mod:<8}  roc=0x{ch.roc:02X}  s={ch.slot:>2}  "
                     f"ch={ch.channel:>2}")
            sort_key = (0 if ch.module else 1, _natural_sort_key(mod))
            items.append((sort_key, key, label))
        items.sort(key=lambda x: x[0])

        self._combo.blockSignals(True)
        self._combo.clear()
        self._combo_keys: List[Tuple[int, int, int]] = []
        sel_idx = 0
        for i, (_, key, label) in enumerate(items):
            self._combo.addItem(label)
            self._combo_keys.append(key)
            if key == prev_key:
                sel_idx = i
        if self._combo_keys:
            self._combo.setCurrentIndex(sel_idx)
            self._selected_key = self._combo_keys[sel_idx]
        self._combo.blockSignals(False)

    def _on_combo_changed(self, idx: int):
        if not hasattr(self, "_combo_keys"):
            return
        if 0 <= idx < len(self._combo_keys):
            self._selected_key = self._combo_keys[idx]
            hits = self._channels.get(self._selected_key)
            self._geo.set_selected_module(hits.module if hits else None)
            # Re-render: hists from cache, waveform from current event
            self._display_hists_for_selected()
            if self._current_idx >= 0 and self._ch is not None:
                # Re-pull the waveform for the new module from current event
                fadc_evt = self._ch.fadc()
                self._display_waveform(fadc_evt)

    def _on_geo_clicked(self, name: str):
        """Geo-view click: switch combo to the (first) channel for this module."""
        if not name:
            return
        for i, key in enumerate(getattr(self, "_combo_keys", [])):
            hits = self._channels.get(key)
            if hits and hits.module == name:
                self._combo.setCurrentIndex(i)   # triggers _on_combo_changed
                return
        # Module exists in geometry but hasn't been seen yet — ignore the click
        self.statusBar().showMessage(
            f"{name}: no events seen yet for this module", 2000)

    # -- accumulate + display --

    def _accumulate_and_display(self, fadc_evt, info, trig_ok: bool):
        # Always analyse the selected channel (for waveform display).  When
        # "Accumulate all modules" is on, analyse every channel.  Histogram
        # fills are dedup'd via self._accumulated: an event already folded in
        # still gets re-analysed for display, but isn't counted again.
        sel_peaks: List[Peak] = []
        accum_all = self._accum_all_cb.isChecked()
        threshold = self._hist_threshold
        wcfg = self._wcfg
        sel_key = self._selected_key
        idx = self._current_idx
        already = (self._accumulated is not None and 0 <= idx < self._accumulated.size
                   and bool(self._accumulated[idx]))
        do_fill = trig_ok and not already

        current_vals: Dict[str, float] = {}   # module_name -> max peak integral

        for r in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = (roc_tag, s, c)
                    is_sel = (key == sel_key)
                    if not is_sel and not accum_all:
                        continue
                    hits = self._channels.get(key)
                    if hits is None:
                        continue
                    samples = slot.channel(c).samples
                    if samples.size < 10:
                        continue
                    _, _, peaks = analyze(samples, wcfg)
                    if is_sel:
                        sel_peaks = peaks
                    # Max peak integral (above threshold) for the geo view's
                    # "current" mode.
                    above = [p.integral for p in peaks if p.height >= threshold]
                    if above and hits.module:
                        current_vals[hits.module] = max(above)
                    if not do_fill:
                        continue
                    kept = 0
                    for p in peaks:
                        if p.height >= threshold:
                            hits.height.fill(p.height)
                            hits.integral.fill(p.integral)
                            hits.position.fill(p.time)
                            kept += 1
                    hits.npeaks.fill(kept)
                    hits.events += 1
                    if kept > 0:
                        hits.peak_events += 1

        if do_fill and self._accumulated is not None and 0 <= idx < self._accumulated.size:
            self._accumulated[idx] = True

        self._geo.set_current_values(current_vals)
        # Overall occupancy only needs refreshing when hists actually changed.
        if do_fill:
            self._geo.set_overall_values(self._compute_overall_occupancy())

        if sel_key is None:
            self._set_info_line(info, peaks=None)
            self._wave.clear("(select a module to view its waveform)")
        else:
            self._set_info_line(info, peaks=sel_peaks)
        self._display_hists_for_selected()
        self._display_waveform(fadc_evt)

    def _compute_overall_occupancy(self) -> Dict[str, float]:
        """module_name -> events_with_peak / events_accumulated (skip empty)."""
        out: Dict[str, float] = {}
        for hits in self._channels.values():
            if hits.module and hits.events > 0:
                out[hits.module] = hits.peak_events / hits.events
        return out

    def _set_info_line(self, info, peaks: Optional[List[Peak]]):
        pieces = [
            f"event #{self._current_idx:,}",
            f"tb=0x{int(info.trigger_bits):X}",
            f"evnum={int(info.event_number)}",
        ]
        if peaks is not None:
            pieces.append(f"peaks(this view)={len(peaks)}")
        key = self._selected_key
        if key:
            hits = self._channels.get(key)
            if hits:
                pieces.append(f"accum={hits.events:,}")
                pieces.append(f"w/peak={hits.peak_events:,}")
        self._info.setText("   ".join(pieces))

    def _display_hists_for_selected(self):
        key = self._selected_key
        hits = self._channels.get(key) if key else None
        if hits is None:
            self._clear_hists()
            return
        mod = hits.module or "(unmapped)"
        self._h_height.set_data(hits.height.bins, hits.height.bmin,
                                hits.height.bstep,
                                under=hits.height.under, over=hits.height.over,
                                title=f"{mod}  —  Peak Height [ADC]",
                                color="#e599f7")
        self._h_integral.set_data(hits.integral.bins, hits.integral.bmin,
                                  hits.integral.bstep,
                                  under=hits.integral.under, over=hits.integral.over,
                                  title=f"{mod}  —  Peak Integral [ADC·sample]",
                                  color="#00b4d8")
        self._h_position.set_data(hits.position.bins, hits.position.bmin,
                                  hits.position.bstep,
                                  under=hits.position.under, over=hits.position.over,
                                  title=f"{mod}  —  Peak Time [ns]",
                                  color="#51cf66")
        self._h_npeaks.set_data(hits.npeaks.bins, hits.npeaks.bmin,
                                hits.npeaks.bstep,
                                under=hits.npeaks.under, over=hits.npeaks.over,
                                title=f"{mod}  —  Peaks / Event",
                                color="#ffa657")

    def _display_waveform(self, fadc_evt):
        key = self._selected_key
        if not key:
            self._wave.clear("(select a module)")
            return
        roc_tag, slot, ch = key
        samples = _find_channel_samples(fadc_evt, roc_tag, slot, ch)
        if samples is None or samples.size == 0:
            # Keep the stacker intact — matches resources/waveform.js:105
            # where empty events in stack mode return without touching the
            # plot.
            if self._wave.is_stacking():
                return
            hits = self._channels.get(key)
            mod = hits.module if hits and hits.module else "(unmapped)"
            self._wave.clear(f"{mod} not present in event #{self._current_idx}")
            return
        ped_mean, ped_rms, peaks = analyze(samples, self._wcfg)
        hits = self._channels.get(key)
        mod = hits.module if hits and hits.module else "(unmapped)"
        self._wave.set_data(samples, peaks, ped_mean, ped_rms,
                            title=(f"{mod}   roc=0x{roc_tag:02X}  "
                                   f"slot={slot}  ch={ch}"),
                            clk_mhz=self._wcfg.clk_mhz,
                            stack_key=f"{roc_tag:02X}_{slot}_{ch}")

    def _clear_hists(self):
        self._h_height.clear("Peak Height")
        self._h_integral.clear("Peak Integral")
        self._h_position.clear("Peak Time")
        self._h_npeaks.clear("Peaks / Event")

    def _clear_plots(self):
        self._clear_hists()
        self._wave.clear("(open an evio file and click Next)")

    def _reset_current_hists(self):
        key = self._selected_key
        if not key: return
        hits = self._channels.get(key)
        if not hits: return
        hits.height.reset()
        hits.integral.reset()
        hits.position.reset()
        hits.npeaks.reset()
        hits.events = 0
        hits.peak_events = 0
        self._display_hists_for_selected()
        self._geo.set_overall_values(self._compute_overall_occupancy())
        self.statusBar().showMessage(
            f"Reset histograms for {hits.module or '(unmapped)'}")

    # -- batch 10k --

    def _on_batch_10k(self):
        if self._batch_thread is not None:
            QMessageBox.information(self, "Busy", "Batch already running.")
            return
        if self._selected_key is None:
            QMessageBox.information(self, "No module selected",
                                    "Select a module first (browse to at "
                                    "least one event so the list populates).")
            return
        start_idx = max(0, self._current_idx + 1)
        if start_idx >= len(self._index):
            QMessageBox.information(self, "End of file",
                                    "Already at the last physics event.")
            return
        remaining = len(self._index) - start_idx
        count = min(10_000, remaining)

        hits = self._channels.get(self._selected_key)
        if hits is None:
            return

        worker = BatchWorker(
            evio_path=str(self._evio_path),
            daq_config_path=self._daq_cfg_path,
            index=self._index, start_idx=start_idx, count=count,
            target_key=self._selected_key,
            channels=self._channels, wcfg=self._wcfg,
            hist_threshold=self._hist_threshold,
            accept_mask=self._accept_mask, reject_mask=self._reject_mask,
            accum_all=self._accum_all_cb.isChecked(),
            accumulated=self._accumulated,
        )
        thread = QThread(self)
        worker.moveToThread(thread)

        def _on_progress(done: int, target: int, peaks: int):
            self._batch_status.setText(
                f"batch: {done:,}/{target:,}  peaks={peaks:,}")
            # refresh hists + geo overall map incrementally
            self._display_hists_for_selected()
            self._geo.set_overall_values(self._compute_overall_occupancy())

        def _on_finished(n: int):
            self._current_idx = start_idx + n - 1
            self._event_spin.blockSignals(True)
            self._event_spin.setValue(max(0, self._current_idx))
            self._event_spin.blockSignals(False)
            self._batch_status.setText(
                f"batch done: {n:,} events processed")
            self._display_hists_for_selected()
            self._geo.set_overall_values(self._compute_overall_occupancy())
            # advance one more to show the next waveform
            if self._current_idx + 1 < len(self._index):
                self._goto(self._current_idx + 1)

        def _on_failed(msg: str):
            QMessageBox.critical(self, "Batch failed", msg)
            self._batch_status.setText("batch failed")

        def _cleanup():
            self._batch_btn.setEnabled(True)
            self._cancel_batch_btn.setVisible(False)
            self._batch_worker = None
            self._batch_thread = None

        thread.started.connect(worker.run)
        worker.progressed.connect(_on_progress)
        worker.finished.connect(_on_finished)
        worker.failed.connect(_on_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(_cleanup)

        self._batch_worker = worker
        self._batch_thread = thread
        self._batch_btn.setEnabled(False)
        self._cancel_batch_btn.setVisible(True)
        self._batch_status.setText(f"batch: 0/{count:,}")
        thread.start()

    def _on_cancel_batch(self):
        if self._batch_worker is not None:
            self._batch_worker.request_cancel()
            self._batch_status.setText("batch: cancelling…")

    # -- JSON save --

    def _save_json_dialog(self):
        if not self._channels:
            return
        default = (self._evio_path.name + ".waveform.json"
                   if self._evio_path else "waveform_hist.json")
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save histograms as JSON", str(Path.cwd() / default),
            "JSON files (*.json)")
        if not path_str:
            return
        try:
            out = {
                "source_file": str(self._evio_path) if self._evio_path else "",
                "height_hist":   {"min": self._h_cfg["min"], "max": self._h_cfg["max"],
                                  "step": self._h_cfg["step"]},
                "integral_hist": {"min": self._i_cfg["min"], "max": self._i_cfg["max"],
                                  "step": self._i_cfg["step"]},
                "position_hist": {"min": self._p_cfg["min"], "max": self._p_cfg["max"],
                                  "step": self._p_cfg["step"]},
                "npeaks_hist":   {"min": self._n_cfg["min"], "max": self._n_cfg["max"],
                                  "step": self._n_cfg["step"]},
                "threshold":     self._hist_threshold,
                "wave_config":   self._wcfg.__dict__.copy(),
                "channels":      {},
            }
            for (roc, slot, ch), hits in sorted(self._channels.items()):
                out["channels"][f"{roc}_{slot}_{ch}"] = {
                    "module":      hits.module,
                    "roc":         hits.roc,
                    "slot":        hits.slot,
                    "channel":     hits.channel,
                    "events":      hits.events,
                    "peak_events": hits.peak_events,
                    "height_hist":   hits.height.to_json(),
                    "integral_hist": hits.integral.to_json(),
                    "position_hist": hits.position.to_json(),
                    "npeaks_hist":   hits.npeaks.to_json(),
                }
            Path(path_str).parent.mkdir(parents=True, exist_ok=True)
            with open(path_str, "w", encoding="utf-8") as f:
                json.dump(out, f)
            self.statusBar().showMessage(f"Saved {path_str}")
        except Exception as ex:
            QMessageBox.warning(self, "Save failed", f"{path_str}\n\n{ex}")

    # -- close --

    def closeEvent(self, ev):
        if self._idx_worker is not None:
            self._idx_worker.request_cancel()
        if self._batch_worker is not None:
            self._batch_worker.request_cancel()
        for thr in (self._idx_thread, self._batch_thread):
            if thr is not None and thr.isRunning():
                thr.quit()
                thr.wait(3000)
        self._close_reader()
        super().closeEvent(ev)


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Waveform viewer — browse an evio file event-by-event.")
    ap.add_argument("path", nargs="?", type=Path,
                    help="evio file to open (otherwise use File → Open…).")
    ap.add_argument("--config", type=Path,
                    default=_REPO_DIR / "database" / "config.json",
                    help="Main config.json (for waveform binning).")
    ap.add_argument("--daq-config", type=Path,
                    default=_REPO_DIR / "database" / "daq_config.json",
                    help="daq_config.json (ROC-tag → crate mapping).")
    ap.add_argument("--daq-map", type=Path,
                    default=_REPO_DIR / "database" / "daq_map.json",
                    help="daq_map.json (module-name lookup).")
    ap.add_argument("--hycal-modules", type=Path,
                    default=_REPO_DIR / "database" / "hycal_modules.json",
                    help="hycal_modules.json (module geometry for the geo selector).")
    ap.add_argument("--trigger-bits", type=Path,
                    default=_REPO_DIR / "database" / "trigger_bits.json",
                    help="trigger_bits.json (for --accept/--reject-trigger).")
    ap.add_argument("--accept-trigger", action="append", default=[],
                    metavar="NAME",
                    help="Require at least one of these trigger bits (repeatable).")
    ap.add_argument("--reject-trigger", action="append", default=None,
                    metavar="NAME",
                    help="Drop events with any of these trigger bits (repeatable). "
                         "Default: uses config.json setting.")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

    hist_cfg      = load_hist_config(args.config)      if args.config.is_file()      else {}
    roc_to_crate  = load_roc_tag_map(args.daq_config)  if args.daq_config.is_file()  else {}
    daq_map       = load_daq_map(args.daq_map)         if args.daq_map.is_file()     else {}
    bit_map       = load_trigger_bit_map(args.trigger_bits)
    hycal_modules = (load_geo_modules(args.hycal_modules)
                     if args.hycal_modules.is_file() else [])

    accept_names = args.accept_trigger or hist_cfg.get("accept_trigger_bits", []) or []
    reject_names = (args.reject_trigger if args.reject_trigger is not None
                    else hist_cfg.get("reject_trigger_bits", []) or [])
    accept_mask = _mask_from_names(accept_names, bit_map) if accept_names else 0
    reject_mask = _mask_from_names(reject_names, bit_map) if reject_names else 0

    app = QApplication(sys.argv)
    win = WaveformViewerWindow(
        hist_config     = hist_cfg,
        daq_map         = daq_map,
        roc_to_crate    = roc_to_crate,
        accept_mask     = accept_mask,
        reject_mask     = reject_mask,
        daq_config_path = str(args.daq_config) if args.daq_config.is_file() else "",
        hycal_modules   = hycal_modules,
    )
    win.show()
    if args.path is not None:
        QTimer.singleShot(0, lambda: win.open_path(args.path))
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
