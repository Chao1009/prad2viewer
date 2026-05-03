#!/usr/bin/env python3
"""
HyCal Pedestal Monitor (PyQt6)
==============================
Measures FADC250 pedestals on all 7 HyCal crates via SSH, displays
colour-coded HyCal maps, and reports channels with irregular sigma.

RMS (sigma) is parsed from the faV3peds stdout, not from saved files.
Saved .cnf files contain only pedestal means.

Usage
-----
    python hycal_pedestal_monitor.py            # view existing data
    python hycal_pedestal_monitor.py --sim       # test with simulated data
"""

from __future__ import annotations

import argparse
import json
import random
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QMessageBox,
    QFileDialog, QLineEdit,
)
from PyQt6.QtCore import Qt, QRectF, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QFont

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget, PALETTES, PALETTE_NAMES,
    apply_theme_palette, set_theme, available_themes, THEME,
)


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_map.json"
DAQ_MAP_JSON = DB_DIR / "hycal_map.json"
PEDESTALS_DIR = SCRIPT_DIR / ".." / "pedestals"
ORIGINAL_PED_DIR = Path("/usr/clas12/release/2.0.0/parms/fadc250/peds")

NUM_CRATES = 7
CRATE_NAMES = [f"adchycal{i}" for i in range(1, NUM_CRATES + 1)]
CHANNELS_PER_SLOT = 16

# Display ranges for the maps  (adjust as needed)
DISPLAY_PED_MIN = 50.0
DISPLAY_PED_MAX = 300.0
DISPLAY_DELTA_MIN = -3.0
DISPLAY_DELTA_MAX = 3.0
DISPLAY_RMS_MIN = 0.0
DISPLAY_RMS_MAX = 1.5

# Thresholds for flagging irregular channels  (adjust as needed)
THRESH_PED_MIN = 50.0       # acceptable pedestal mean lower bound
THRESH_PED_MAX = 300.0      # acceptable pedestal mean upper bound
THRESH_DEAD_AVG = 1.0       # avg below this AND rms below THRESH_DEAD_RMS -> DEAD
THRESH_DEAD_RMS = 0.1       # rms below this AND avg below THRESH_DEAD_AVG -> DEAD
THRESH_HIGH_RMS = 1.5       # rms above this -> HIGH RMS
THRESH_DRIFT = 3.0          # |current - configured| above this -> DRIFT

# LMS / V module positions below HyCal  (name -> centre-x)
_BOTTOM_Y = -640.0
_BOTTOM_SZ = 50.0
_LMS_V_XPOS = {
    "LMS1": -170.0, "LMS2": -115.0, "LMS3": -60.0,
    "V1": 5.0, "V2": 60.0, "V3": 115.0, "V4": 170.0,
}


# ===========================================================================
#  Data loading helpers
# ===========================================================================

def prepare_modules(modules: List[Module]) -> List[Module]:
    """Reposition LMS1-3 below HyCal and add V1-V4."""
    result: List[Module] = []
    for m in modules:
        if m.name in _LMS_V_XPOS:
            result.append(Module(m.name, m.mod_type,
                                 _LMS_V_XPOS[m.name], _BOTTOM_Y,
                                 _BOTTOM_SZ, _BOTTOM_SZ))
        else:
            result.append(m)
    for name in ("V1", "V2", "V3", "V4"):
        result.append(Module(name, "Veto",
                             _LMS_V_XPOS[name], _BOTTOM_Y,
                             _BOTTOM_SZ, _BOTTOM_SZ))
    return result


def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    """(crate_index, slot, channel) -> module_name from hycal_map.json.

    Records without a ``daq`` block (boosters, PRad-1 V1-V4) are skipped.
    """
    with open(path) as f:
        data = json.load(f)
    out: Dict[Tuple[int, int, int], str] = {}
    for e in data:
        d = e.get("daq")
        if not d:
            continue
        out[(d["crate"], d["slot"], d["channel"])] = e["n"]
    return out


# ===========================================================================
#  Pedestal .cnf parser  (means only -- RMS is NOT in saved files)
# ===========================================================================

def parse_pedestal_file(filepath: Path) -> Dict[int, List[float]]:
    """Return  slot_number -> [16 pedestal means]."""
    slots: Dict[int, List[float]] = {}
    cur_slot: Optional[int] = None
    vals: List[float] = []
    reading = False

    def _flush():
        nonlocal reading, vals
        if cur_slot is not None and reading and vals:
            slots[cur_slot] = vals[:CHANNELS_PER_SLOT]
        vals = []
        reading = False

    with open(filepath) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("FADC250_CRATE"):
                _flush()
            elif line.startswith("FADC250_SLOT"):
                _flush()
                cur_slot = int(line.split()[1])
            elif line.startswith("FADC250_ALLCH_PED"):
                _flush()
                reading = True
                vals = [float(v) for v in
                        line[len("FADC250_ALLCH_PED"):].split()]
                if len(vals) >= CHANNELS_PER_SLOT:
                    _flush()
            elif reading:
                try:
                    vals.extend(float(v) for v in line.split())
                    if len(vals) >= CHANNELS_PER_SLOT:
                        _flush()
                except ValueError:
                    _flush()
    _flush()
    return slots


def read_all_pedestals(
    ped_dir: Path, suffix: str,
    daq_map: Dict[Tuple[int, int, int], str],
) -> Dict[str, float]:
    """Read all 7 crate files.  Returns  module_name -> pedestal_mean."""
    result: Dict[str, float] = {}
    for ci, cname in enumerate(CRATE_NAMES):
        fp = ped_dir / f"{cname}{suffix}"
        if not fp.exists():
            continue
        for slot, peds in parse_pedestal_file(fp).items():
            for ch, val in enumerate(peds):
                mod = daq_map.get((ci, slot, ch))
                if mod is not None:
                    result[mod] = val
    return result


# ===========================================================================
#  Parse faV3peds stdout for per-channel RMS
# ===========================================================================

_PED_RE = re.compile(
    r"faV3MeasureChannelPedestal:\s*slot\s*(\d+),\s*chan\s*(\d+)\s*=>\s*"
    r"avg\s+([\d.]+),\s*rms\s+([\d.]+),\s*min\s+(\d+),\s*max\s+(\d+)"
)


def parse_measurement_stdout(
    text: str, crate_idx: int,
    daq_map: Dict[Tuple[int, int, int], str],
) -> Dict[str, dict]:
    """Parse stdout lines.  Returns  module_name -> {avg, rms, min, max}."""
    result: Dict[str, dict] = {}
    for m in _PED_RE.finditer(text):
        slot, chan = int(m.group(1)), int(m.group(2))
        mod = daq_map.get((crate_idx, slot, chan))
        if mod is not None:
            result[mod] = {
                "avg": float(m.group(3)), "rms": float(m.group(4)),
                "min": int(m.group(5)),   "max": int(m.group(6)),
            }
    return result


# ===========================================================================
#  Irregular channel detection
# ===========================================================================

def find_irregular_channels(
    measured: Dict[str, dict],
    configured: Dict[str, float],
    daq_map: Dict[Tuple[int, int, int], str],
) -> List[str]:
    """Return formatted lines describing flagged channels."""
    rev: Dict[str, Tuple[str, int, int]] = {}
    for (ci, slot, ch), name in daq_map.items():
        rev[name] = (CRATE_NAMES[ci], slot, ch)

    issues: List[str] = []
    for mod, d in sorted(measured.items(),
                         key=lambda kv: rev.get(kv[0], ("", 0, 0))):
        cname, slot, ch = rev.get(mod, ("???", 0, 0))
        loc = f"{cname} slot {slot:2d} ch {ch:2d}"
        avg, rms = d["avg"], d["rms"]

        if avg < THRESH_DEAD_AVG and rms < THRESH_DEAD_RMS:
            issues.append(f"  DEAD          {mod:<6s}  {loc}  "
                          f"avg={avg:.2f}  rms={rms:.3f}")
        elif avg < THRESH_PED_MIN or avg > THRESH_PED_MAX:
            issues.append(f"  OUT OF RANGE  {mod:<6s}  {loc}  "
                          f"avg={avg:.2f}  rms={rms:.3f}  "
                          f"(valid: {THRESH_PED_MIN:.0f}-{THRESH_PED_MAX:.0f})")
        elif rms > THRESH_HIGH_RMS:
            issues.append(f"  HIGH RMS      {mod:<6s}  {loc}  "
                          f"avg={avg:.2f}  rms={rms:.3f}")

        if mod in configured and not (avg < THRESH_DEAD_AVG and rms < THRESH_DEAD_RMS):
            delta = avg - configured[mod]
            if abs(delta) > THRESH_DRIFT:
                issues.append(f"  DRIFT         {mod:<6s}  {loc}  "
                              f"cur={avg:.2f}  conf={configured[mod]:.2f}  "
                              f"delta={delta:+.2f}")
    return issues


# ===========================================================================
#  Colour helpers
# ===========================================================================

def _time_ago(epoch: float) -> str:
    """Format seconds-since-epoch as a human-readable 'X ago' string."""
    delta = int(time.time() - epoch)
    if delta < 0:
        return "just now"
    if delta < 60:
        return f"{delta}s ago"
    mins = delta // 60
    if mins < 60:
        return f"{mins}min ago"
    hours = mins // 60
    mins_r = mins % 60
    if hours < 24:
        return f"{hours}h {mins_r}min ago" if mins_r else f"{hours}h ago"
    days = hours // 24
    hours_r = hours % 24
    return f"{days}d {hours_r}h ago" if hours_r else f"{days}d ago"


def _ped_mtime(ped_dir: Path, suffix: str) -> Optional[float]:
    """Return the most recent mtime among pedestal files in *ped_dir*, or None."""
    newest = None
    for cname in CRATE_NAMES:
        fp = ped_dir / f"{cname}{suffix}"
        if fp.exists():
            mt = fp.stat().st_mtime
            if newest is None or mt > newest:
                newest = mt
    return newest


# ===========================================================================
#  HyCal map widget
# ===========================================================================

_LABEL_NAMES = {"LMS1", "LMS2", "LMS3", "V1", "V2", "V3", "V4"}


class PedestalMapWidget(HyCalMapWidget):
    """HyCal map with a centred title and labelled LMS/V modules."""

    CB_MAX_WIDTH = 300

    def __init__(self, parent=None):
        super().__init__(parent, include_lms=True, margin_top=30)
        self._title = ""

    def set_data(self, modules: List[Module], values: Dict[str, float],
                 title: str, vmin: float, vmax: float):
        self._title = title
        self.set_modules(modules)
        self._values = values
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def _fmt_value(self, v: float) -> str:
        return f"{v:.1f}"

    def _colorbar_center_text(self) -> str:
        return PALETTE_NAMES[self._palette_idx]   # no [log] flag

    def _paint_before_modules(self, p, w: int, h: int):
        if self._title:
            p.setPen(QColor(THEME.TEXT))
            p.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
            p.drawText(QRectF(0, 4, w, 24),
                       Qt.AlignmentFlag.AlignCenter, self._title)

    def _paint_empty(self, p, w: int, h: int):
        # title-only state (pre-data load)
        if self._title:
            p.setPen(QColor(THEME.TEXT))
            p.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
            p.drawText(QRectF(0, 4, w, 24),
                       Qt.AlignmentFlag.AlignCenter, self._title)
        if not self._values:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Monospace", 12))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter, "No data")

    def _paint_overlays(self, p, w: int, h: int):
        # LMS / V labels
        p.setPen(QColor(THEME.TEXT))
        p.setFont(QFont("Monospace", 7, QFont.Weight.Bold))
        for name in _LABEL_NAMES:
            r = self._rects.get(name)
            if r is not None:
                p.drawText(r, Qt.AlignmentFlag.AlignCenter, name)
        # hover highlight (from base)
        super()._paint_overlays(p, w, h)


# ===========================================================================
#  Measurement thread
# ===========================================================================

class MeasureThread(QThread):
    progress = pyqtSignal(int, str)
    crate_done = pyqtSignal(int, str)
    crate_error = pyqtSignal(int, str)

    def run(self):
        for i, cname in enumerate(CRATE_NAMES):
            self.progress.emit(
                i, f"Measuring {cname} ({i + 1}/{NUM_CRATES})...")
            cmd = (f'ssh {cname} '
                   f'"cd ~/prad2_daq/prad2evviewer/pedestals; '
                   f'faV3peds {cname}_latest.cnf"')
            try:
                r = subprocess.run(cmd, shell=True, capture_output=True,
                                   text=True, timeout=180)
                self.crate_done.emit(i, r.stdout + "\n" + r.stderr)
            except subprocess.TimeoutExpired:
                self.crate_error.emit(i, f"{cname}: TIMEOUT (180 s)")
            except Exception as e:
                self.crate_error.emit(i, f"{cname}: {e}")


# ===========================================================================
#  Main window
# ===========================================================================

class PedestalMonitorWindow(QMainWindow):

    def __init__(self, modules: List[Module],
                 daq_map: Dict[Tuple[int, int, int], str],
                 sim: bool = False):
        super().__init__()
        self._modules = modules
        self._daq_map = daq_map
        self._sim = sim
        self._palette_idx_left = 0
        self._palette_idx_right = 0
        self._right_mode = "delta"   # "delta" or "rms"

        self._configured: Dict[str, float] = {}
        self._latest: Dict[str, float] = {}
        self._measured: Dict[str, dict] = {}

        self._build_ui()
        self._load_data()

    # ---- UI ----

    def _build_ui(self):
        self.setWindowTitle("HyCal Pedestal Monitor")
        self.resize(1600, 1000)
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL PEDESTAL MONITOR")
        lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        lbl.setStyleSheet(f"color:{THEME.ACCENT};")
        top.addWidget(lbl)
        top.addStretch()
        self._measure_btn = self._make_btn(
            "Measure Pedestals", THEME.SUCCESS, self._on_measure)
        top.addWidget(self._measure_btn)
        self._reload_btn = self._make_btn(
            "Reload Files", THEME.TEXT, self._load_data)
        top.addWidget(self._reload_btn)
        self._save_btn = self._make_btn(
            "Save Report", THEME.WARN, self._on_save_report)
        top.addWidget(self._save_btn)
        root.addLayout(top)

        # -- maps --
        maps = QWidget()
        ml = QHBoxLayout(maps)
        ml.setContentsMargins(0, 0, 0, 0)
        ml.setSpacing(8)
        self._map_left = PedestalMapWidget()
        self._map_right = PedestalMapWidget()
        self._map_left.moduleHovered.connect(self._on_hover)
        self._map_right.moduleHovered.connect(self._on_hover)
        self._map_left.paletteClicked.connect(self._cycle_palette_left)
        self._map_right.paletteClicked.connect(self._cycle_palette_right)
        ml.addWidget(self._map_left)
        ml.addWidget(self._map_right)
        root.addWidget(maps, stretch=1)

        # -- range controls --
        rng = QHBoxLayout()
        rng.addWidget(self._slabel("Mean range:"))
        self._left_min = self._sedit(f"{DISPLAY_PED_MIN:.0f}")
        self._left_max = self._sedit(f"{DISPLAY_PED_MAX:.0f}")
        rng.addWidget(self._left_min)
        rng.addWidget(self._slabel("-"))
        rng.addWidget(self._left_max)
        rng.addSpacing(20)

        # Right panel mode toggle
        self._right_mode_btn = QPushButton("Delta")
        self._right_mode_btn.setFixedWidth(75)
        self._right_mode_btn.setStyleSheet(
            f"QPushButton{{background:{THEME.BUTTON};color:{THEME.ACCENT};"
            f"border:1px solid {THEME.BORDER};padding:4px 10px;"
            f"font:bold 13px Monospace;border-radius:8px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")
        self._right_mode_btn.clicked.connect(self._toggle_right_mode)
        rng.addWidget(self._right_mode_btn)
        rng.addWidget(self._slabel("range:"))
        self._right_min = self._sedit(f"{DISPLAY_DELTA_MIN:.1f}")
        self._right_max = self._sedit(f"{DISPLAY_DELTA_MAX:.1f}")
        rng.addWidget(self._right_min)
        rng.addWidget(self._slabel("-"))
        rng.addWidget(self._right_max)
        rng.addSpacing(10)
        rng.addWidget(self._make_btn("Apply", THEME.TEXT,
                                     self._apply_ranges))
        rng.addStretch()
        root.addLayout(rng)

        # -- info bar --
        self._info = QLabel("Hover over a module for details")
        self._info.setFont(QFont("Monospace", 11))
        self._info.setStyleSheet(
            f"QLabel{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"padding:4px 8px;border:1px solid {THEME.BORDER};"
            f"border-radius:8px;}}")
        self._info.setFixedHeight(28)
        root.addWidget(self._info)

        # -- status bar (prominent, for measurement progress) --
        self._status_lbl = QLabel("Ready")
        self._status_lbl.setFont(QFont("Monospace", 12, QFont.Weight.Bold))
        self._status_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._set_status_style("idle")
        self._status_lbl.setFixedHeight(34)
        root.addWidget(self._status_lbl)

        # -- report area --
        self._report = QTextEdit()
        self._report.setReadOnly(True)
        self._report.setFont(QFont("Monospace", 10))
        self._report.setStyleSheet(
            f"QTextEdit{{background:{THEME.PANEL};color:{THEME.TEXT_DIM};"
            f"border:1px solid {THEME.BORDER};border-radius:8px;}}")
        self._report.setMaximumHeight(180)
        root.addWidget(self._report)

    def _set_status_style(self, mode: str):
        # Status rows use a translucent variant of the semantic colour as
        # background to stay within the active theme. Qt's QSS supports
        # rgba(...) directly, but keeping the fill solid via a near-black
        # tint is simpler and theme-neutral: reuse PANEL as the base.
        if mode == "measuring":
            bg, fg = THEME.PANEL, THEME.WARN
            border = f"2px solid {THEME.WARN}"
        elif mode == "done":
            bg, fg = THEME.PANEL, THEME.SUCCESS
            border = f"2px solid {THEME.SUCCESS}"
        elif mode == "error":
            bg, fg = THEME.PANEL, THEME.DANGER
            border = f"2px solid {THEME.DANGER}"
        else:
            bg, fg = THEME.PANEL, THEME.TEXT_DIM
            border = f"1px solid {THEME.BORDER}"
        self._status_lbl.setStyleSheet(
            f"QLabel{{background:{bg};color:{fg};padding:4px;"
            f"border:{border};border-radius:8px;}}")

    # -- helpers --

    def _make_btn(self, text, fg, slot):
        btn = QPushButton(text)
        btn.setStyleSheet(
            f"QPushButton{{background:{THEME.BUTTON};color:{fg};"
            f"border:1px solid {THEME.BORDER};padding:6px 16px;"
            f"font:bold 12px Monospace;border-radius:8px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}"
            f"QPushButton:disabled{{color:{THEME.TEXT_MUTED};}}")
        btn.clicked.connect(slot)
        return btn

    def _slabel(self, text):
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 10))
        lbl.setStyleSheet(f"color:{THEME.TEXT};")
        return lbl

    def _sedit(self, text):
        e = QLineEdit(text)
        e.setFixedWidth(55)
        e.setFont(QFont("Monospace", 10))
        e.setStyleSheet(
            f"QLineEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};border-radius:8px;"
            f"padding:2px 4px;}}")
        e.returnPressed.connect(self._apply_ranges)
        return e

    # ---- palette cycling (independent per map) ----

    def _cycle_palette_left(self):
        self._palette_idx_left = (self._palette_idx_left + 1) % len(PALETTES)
        self._map_left.set_palette(self._palette_idx_left)

    def _cycle_palette_right(self):
        self._palette_idx_right = (self._palette_idx_right + 1) % len(PALETTES)
        self._map_right.set_palette(self._palette_idx_right)

    # ---- right panel mode toggle ----

    def _toggle_right_mode(self):
        if self._right_mode == "delta":
            self._right_mode = "rms"
            self._right_mode_btn.setText("RMS")
            self._right_min.setText(f"{DISPLAY_RMS_MIN:.1f}")
            self._right_max.setText(f"{DISPLAY_RMS_MAX:.1f}")
        else:
            self._right_mode = "delta"
            self._right_mode_btn.setText("Delta")
            self._right_min.setText(f"{DISPLAY_DELTA_MIN:.1f}")
            self._right_max.setText(f"{DISPLAY_DELTA_MAX:.1f}")
        self._update_right_map()
        self._apply_ranges()

    # ---- range editing ----

    def _apply_ranges(self):
        try:
            lmin = float(self._left_min.text())
            lmax = float(self._left_max.text())
            if lmin < lmax:
                self._map_left.set_range(lmin, lmax)
        except ValueError:
            pass
        try:
            rmin = float(self._right_min.text())
            rmax = float(self._right_max.text())
            if rmin < rmax:
                self._map_right.set_range(rmin, rmax)
        except ValueError:
            pass

    # ---- data loading ----

    def _load_data(self):
        if self._sim:
            self._load_sim_data()
            return
        if ORIGINAL_PED_DIR.exists():
            self._configured = read_all_pedestals(
                ORIGINAL_PED_DIR, "_ped.cnf", self._daq_map)
        else:
            self._configured = {}
        if PEDESTALS_DIR.exists():
            self._latest = read_all_pedestals(
                PEDESTALS_DIR, "_latest.cnf", self._daq_map)
        else:
            self._latest = {}

        # If configured files are newer than latest, latest is stale
        mt_conf = _ped_mtime(ORIGINAL_PED_DIR, "_ped.cnf")
        mt_latest = _ped_mtime(PEDESTALS_DIR, "_latest.cnf")
        if mt_conf and mt_latest and mt_conf > mt_latest:
            self._latest.clear()

        n_o, n_l = len(self._configured), len(self._latest)
        age = ""
        if mt_latest is not None and self._latest:
            age = f"    (measured {_time_ago(mt_latest)})"
        elif mt_latest is not None and not self._latest:
            age = f"    (latest stale -- configured files are newer)"
        self._status_lbl.setText(
            f"Loaded {n_o} configured, {n_l} latest channels{age}")
        self._set_status_style("idle")
        self._update_maps()
        self._update_report()

    def _load_sim_data(self):
        rng = random.Random(42)
        self._configured.clear()
        self._latest.clear()
        self._measured.clear()
        all_names: List[str] = []
        for m in self._modules:
            n = m.name
            all_names.append(n)
            o = rng.gauss(160, 25)
            l = o + rng.gauss(0, 1.0)
            self._configured[n] = o
            self._latest[n] = l
            self._measured[n] = {
                "avg": l, "rms": abs(rng.gauss(0.7, 0.15)),
                "min": int(l) - 3, "max": int(l) + 3,
            }
        for n in rng.sample(all_names, min(15, len(all_names))):
            self._configured[n] = 0.0
            self._latest[n] = 0.0
            self._measured[n].update(avg=0.0, rms=0.0)
        for n in rng.sample(all_names, 3):
            val = rng.choice([rng.uniform(10, 40), rng.uniform(320, 500)])
            self._configured[n] = val
            self._latest[n] = val + rng.gauss(0, 0.5)
            self._measured[n].update(avg=self._latest[n],
                                     rms=abs(rng.gauss(0.7, 0.2)))
        for n in rng.sample(all_names, 5):
            if self._measured[n]["avg"] >= THRESH_DEAD_AVG:
                self._measured[n]["rms"] = rng.uniform(1.8, 5.0)
        for n in rng.sample(all_names, 4):
            if self._measured[n]["avg"] >= THRESH_PED_MIN:
                drift = rng.choice([-1, 1]) * rng.uniform(4.0, 12.0)
                self._latest[n] = self._configured[n] + drift
                self._measured[n]["avg"] = self._latest[n]
        self._update_maps()
        self._update_report()

    # ---- update views ----

    def _update_maps(self):
        has_latest = bool(self._latest)
        cur = self._latest if has_latest else self._configured
        label = "Current" if has_latest else "Configured"

        self._map_left.set_data(
            self._modules, cur, f"{label} Pedestal Mean",
            DISPLAY_PED_MIN, DISPLAY_PED_MAX)
        self._map_left.set_palette(self._palette_idx_left)

        self._update_right_map()

    def _update_right_map(self):
        has_latest = bool(self._latest)
        cur = self._latest if has_latest else self._configured

        if self._right_mode == "rms":
            if self._measured:
                rms = {n: d["rms"] for n, d in self._measured.items()}
                self._map_right.set_data(
                    self._modules, rms, "Pedestal RMS (from measurement)",
                    DISPLAY_RMS_MIN, DISPLAY_RMS_MAX)
            else:
                self._map_right.set_data(
                    self._modules, {},
                    "Pedestal RMS (no measurement data)",
                    DISPLAY_RMS_MIN, DISPLAY_RMS_MAX)
        else:
            if has_latest and self._configured:
                delta = {n: cur[n] - self._configured[n]
                         for n in cur if n in self._configured}
                self._map_right.set_data(
                    self._modules, delta,
                    "Mean Difference (Current \u2212 Configured)",
                    DISPLAY_DELTA_MIN, DISPLAY_DELTA_MAX)
            else:
                self._map_right.set_data(
                    self._modules, {},
                    "Mean Difference (no comparison data)",
                    DISPLAY_DELTA_MIN, DISPLAY_DELTA_MAX)
        self._map_right.set_palette(self._palette_idx_right)

    def _update_report(self):
        lines: List[str] = []

        def _stats(label, peds):
            vals = list(peds.values())
            live = [v for v in vals if v >= THRESH_DEAD_AVG]
            dead = sum(1 for v in vals if v < THRESH_DEAD_AVG)
            if not vals:
                lines.append(f"{label}: no data")
                return
            if live:
                avg = sum(live) / len(live)
                lines.append(
                    f"{label}: {len(vals)} ch, {dead} dead, "
                    f"mean={avg:.1f}  min={min(live):.1f}  max={max(live):.1f}")
            else:
                lines.append(f"{label}: {len(vals)} ch, ALL dead")

        if self._configured:
            _stats("Configured", self._configured)
        if self._latest:
            _stats("Current ", self._latest)

        if self._measured:
            issues = find_irregular_channels(
                self._measured, self._configured, self._daq_map)
            lines.append("")
            if issues:
                lines.append(f"IRREGULAR CHANNELS  ({len(issues)} flagged):")
                lines.extend(issues)
            else:
                lines.append("All channels within normal parameters.")
        self._report.setPlainText("\n".join(lines))

    # ---- hover ----

    def _on_hover(self, name: str):
        parts = [name]
        for m in self._modules:
            if m.name == name:
                parts.append(f"({m.mod_type})")
                break
        if name in self._latest:
            parts.append(f"ped: {self._latest[name]:.2f}")
        elif name in self._configured:
            parts.append(f"ped: {self._configured[name]:.2f}")
        if name in self._configured and name in self._latest:
            delta = self._latest[name] - self._configured[name]
            parts.append(f"conf: {self._configured[name]:.2f}")
            parts.append(f"delta: {delta:+.2f}")
        if name in self._measured:
            parts.append(f"rms: {self._measured[name]['rms']:.3f}")
        self._info.setText("    ".join(parts))

    # ---- save ----

    def _on_save_report(self):
        text = self._report.toPlainText().strip()
        if not text:
            QMessageBox.information(self, "Save Report", "Nothing to save.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Irregular Channel Report",
            str(PEDESTALS_DIR / "pedestal_report.txt"),
            "Text files (*.txt);;All files (*)")
        if not path:
            return
        with open(path, "w") as f:
            f.write(text + "\n")
        self._status_lbl.setText(f"Report saved to {path}")
        self._set_status_style("done")

    # ---- measurement ----

    def _on_measure(self):
        reply = QMessageBox.warning(
            self, "Pedestal Measurement",
            "WARNING: Pedestal measurement will INTERRUPT DAQ running!\n"
            "Only proceed when DAQ is IDLE.\n\n"
            "Proceed with measurement on all 7 crates?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No)
        if reply != QMessageBox.StandardButton.Yes:
            return
        self._measure_btn.setEnabled(False)
        self._reload_btn.setEnabled(False)
        self._measure_btn.setText("Measuring...")
        self._measured.clear()
        self._report.clear()

        self._set_status_style("measuring")
        self._status_lbl.setText("MEASURING: starting...")

        self._thread = MeasureThread()
        self._thread.progress.connect(self._on_progress)
        self._thread.crate_done.connect(self._on_crate_done)
        self._thread.crate_error.connect(self._on_crate_error)
        self._thread.finished.connect(self._on_measure_finished)
        self._thread.start()

    def _on_progress(self, idx: int, msg: str):
        self._status_lbl.setText(f"MEASURING: {msg}")

    def _on_crate_done(self, idx: int, stdout: str):
        parsed = parse_measurement_stdout(stdout, idx, self._daq_map)
        self._measured.update(parsed)
        self._report.append(
            f"  {CRATE_NAMES[idx]}: {len(parsed)} channels measured")

    def _on_crate_error(self, idx: int, msg: str):
        self._report.append(f"  ERROR: {msg}")
        self._set_status_style("error")
        self._status_lbl.setText(f"ERROR: {msg}")

    def _on_measure_finished(self):
        self._measure_btn.setEnabled(True)
        self._reload_btn.setEnabled(True)
        self._measure_btn.setText("Measure Pedestals")

        if PEDESTALS_DIR.exists():
            self._latest = read_all_pedestals(
                PEDESTALS_DIR, "_latest.cnf", self._daq_map)

        n = len(self._measured)
        self._status_lbl.setText(
            f"MEASUREMENT COMPLETE: {n} channels measured")
        self._set_status_style("done")

        self._update_maps()
        self._update_report()


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="HyCal Pedestal Monitor")
    ap.add_argument("--sim", action="store_true",
                    help="Use simulated data for testing")
    ap.add_argument("--modules-db", type=Path, default=MODULES_JSON)
    ap.add_argument("--daq-map", type=Path, default=DAQ_MAP_JSON)
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

    modules = prepare_modules(load_modules(args.modules_db))
    daq_map = load_daq_map(args.daq_map)

    app = QApplication(sys.argv)
    win = PedestalMonitorWindow(modules, daq_map, sim=args.sim)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
