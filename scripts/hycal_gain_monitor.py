#!/usr/bin/env python3
"""
HyCal Gain Monitor (PyQt6)
==========================
Visualises LMS-based gain factors across runs for all HyCal modules.
Reads text-based ``prad_{:06d}_LMS.dat`` files produced by the offline
gain analysis, displays a colour-coded HyCal geo map, LMS reference
channel stability charts, and a table of irregular (module, run) entries.

Usage
-----
    python scripts/hycal_gain_monitor.py
"""

from __future__ import annotations

import glob
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QComboBox, QLineEdit, QDoubleSpinBox, QSpinBox,
    QFileDialog, QSplitter, QSizePolicy, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QToolTip,
)
from PyQt6.QtCore import Qt, QRectF, QPointF, pyqtSignal, QSize, QTimer
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QLinearGradient, QPalette,
)


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_modules.json"

LMS_NAMES = ["LMS1", "LMS2", "LMS3"]
LMS_REF_DEFAULT = 1          # index into LMS_NAMES -> "LMS2"
FILE_PATTERN = re.compile(r"prad_(\d{6})_LMS\.dat$")


# ===========================================================================
#  Module database (self-contained, mirrors calibration/scan_utils.py)
# ===========================================================================

class Module:
    __slots__ = ("name", "mod_type", "x", "y", "sx", "sy")

    def __init__(self, name, mod_type, x, y, sx, sy):
        self.name = name
        self.mod_type = mod_type
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy


def load_modules(path: Path = MODULES_JSON) -> List[Module]:
    import json
    with open(path) as f:
        data = json.load(f)
    return [Module(e["n"], e["t"], e["x"], e["y"], e["sx"], e["sy"])
            for e in data]


# ===========================================================================
#  Colour palettes
# ===========================================================================

PALETTES = {
    "blue-orange": [
        (0.00, (10, 42, 110)), (0.25, (30, 90, 180)),
        (0.50, (80, 80, 80)), (0.75, (220, 120, 30)),
        (1.00, (249, 115, 22)),
    ],
    "viridis": [
        (0.00, (68, 1, 84)), (0.25, (59, 82, 139)),
        (0.50, (33, 145, 140)), (0.75, (94, 201, 98)),
        (1.00, (253, 231, 37)),
    ],
    "inferno": [
        (0.00, (0, 0, 4)), (0.25, (120, 28, 109)),
        (0.50, (229, 89, 52)), (0.75, (253, 198, 39)),
        (1.00, (252, 255, 164)),
    ],
    "coolwarm": [
        (0.00, (59, 76, 192)), (0.25, (141, 176, 254)),
        (0.50, (221, 221, 221)), (0.75, (245, 148, 114)),
        (1.00, (180, 4, 38)),
    ],
    "hot": [
        (0.00, (11, 0, 0)), (0.33, (230, 0, 0)),
        (0.66, (255, 210, 0)), (1.00, (255, 255, 255)),
    ],
    "rainbow": [
        (0.00, (30, 58, 95)), (0.25, (59, 130, 246)),
        (0.50, (45, 212, 160)), (0.75, (234, 179, 8)),
        (1.00, (245, 101, 101)),
    ],
}
PALETTE_NAMES = list(PALETTES.keys())

# Separate palette used only in Run-to-Run Drift mode (not cycled by the user)
DRIFT_PALETTE = [
    (0.00, (0, 210, 230)),   # cyan  — large negative drift
    (0.50, (80, 80, 80)),    # grey  — no drift (always maps to 0 in drift mode)
    (1.00, (249, 115, 22)),  # orange — large positive drift
]


def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)


def _cmap_qcolor(t: float, stops) -> QColor:
    t = max(0.0, min(1.0, t))
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1:
            s = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return QColor(_lerp(c0[0], c1[0], s),
                          _lerp(c0[1], c1[1], s),
                          _lerp(c0[2], c1[2], s))
    _, c = stops[-1]
    return QColor(*c)


# ===========================================================================
#  Data structures
# ===========================================================================

@dataclass
class LMSRecord:
    alpha_peak: float
    alpha_sigma: float
    alpha_chi2ndf: float
    lms_peak: float
    lms_sigma: float
    lms_chi2ndf: float


@dataclass
class ModuleRecord:
    lms_peak: float
    lms_sigma: float
    lms_chi2ndf: float
    gain_factors: Tuple[float, float, float]


@dataclass
class RunData:
    run_number: int
    lms: Dict[str, LMSRecord] = field(default_factory=dict)
    modules: Dict[str, ModuleRecord] = field(default_factory=dict)


@dataclass
class IrregularEntry:
    name: str
    mod_type: str
    run_number: int
    gain: float
    mean_gain: float
    std_dev: float
    deviation_sigma: float


@dataclass
class DriftEntry:
    name: str
    mod_type: str
    run_number: int       # current run
    prev_run_number: int  # previous run
    gain_current: float
    gain_prev: float
    rel_change: float     # (gain_current - gain_prev) / gain_prev


@dataclass
class SummaryEntry:
    name: str
    mod_type: str
    drift_count: int      # number of consecutive run pairs with |Δ| > threshold
    max_rel_change: float # largest |Δ| seen (absolute value)
    max_run: int          # run where max drift occurred
    max_prev_run: int     # previous run for that pair


# ===========================================================================
#  File parsing
# ===========================================================================

def parse_dat_file(filepath: str) -> Optional[RunData]:
    """Parse a single prad_NNNNNN_LMS.dat file."""
    basename = os.path.basename(filepath)
    m = FILE_PATTERN.search(basename)
    if not m:
        return None
    run_number = int(m.group(1))
    rd = RunData(run_number=run_number)

    try:
        with open(filepath) as f:
            lines = f.readlines()
    except OSError:
        return None

    if len(lines) < 3:
        return None

    # First 3 lines: LMS reference channels
    for i in range(3):
        parts = lines[i].strip().replace(',', ' ').split()
        if len(parts) < 7:
            continue
        try:
            name = parts[0]
            rd.lms[name] = LMSRecord(
                alpha_peak=float(parts[1]),
                alpha_sigma=float(parts[2]),
                alpha_chi2ndf=float(parts[3]),
                lms_peak=float(parts[4]),
                lms_sigma=float(parts[5]),
                lms_chi2ndf=float(parts[6]),
            )
        except (ValueError, IndexError):
            continue

    # Remaining lines: module data
    for line in lines[3:]:
        parts = line.strip().replace(',', ' ').split()
        if len(parts) < 7:
            continue
        try:
            name = parts[0]
            rd.modules[name] = ModuleRecord(
                lms_peak=float(parts[1]),
                lms_sigma=float(parts[2]),
                lms_chi2ndf=float(parts[3]),
                gain_factors=(float(parts[4]), float(parts[5]), float(parts[6])),
            )
        except (ValueError, IndexError):
            continue

    return rd


def load_all_runs(folder: str) -> List[RunData]:
    """Scan folder for prad_*_LMS.dat files, parse all, sort by run number."""
    pattern = os.path.join(folder, "prad_*_LMS.dat")
    files = sorted(glob.glob(pattern))
    runs: List[RunData] = []
    for f in files:
        rd = parse_dat_file(f)
        if rd is not None:
            runs.append(rd)
    # files are already sorted alphabetically; zero-padded run numbers preserve numerical order
    return runs


# ===========================================================================
#  Outlier detection
# ===========================================================================

def compute_irregular_entries(
    runs: List[RunData],
    ref_idx: int,
    mod_by_name: Dict[str, Module],
    sigma_threshold: float = 3.0,
    min_runs: int = 5,
) -> List[IrregularEntry]:
    """Find (module, run) pairs with outlier gain factors."""

    # Collect gain values per module across all runs
    # module_name -> [(run_number, gain)]
    all_gains: Dict[str, List[Tuple[int, float]]] = {}
    for rd in runs:
        for mname, mrec in rd.modules.items():
            gains = all_gains.setdefault(mname, [])
            gains.append((rd.run_number, mrec.gain_factors[ref_idx]))

    entries: List[IrregularEntry] = []
    for mname, gains_list in all_gains.items():
        if len(gains_list) < min_runs:
            continue
        values = [g for _, g in gains_list]
        mean = sum(values) / len(values)
        variance = sum((v - mean) ** 2 for v in values) / len(values)
        std = math.sqrt(variance) if variance > 0 else 0.0
        if std == 0:
            continue

        mod = mod_by_name.get(mname)
        mod_type = mod.mod_type if mod else "?"

        for run_num, gain in gains_list:
            dev = abs(gain - mean) / std
            if dev > sigma_threshold:
                entries.append(IrregularEntry(
                    name=mname,
                    mod_type=mod_type,
                    run_number=run_num,
                    gain=gain,
                    mean_gain=mean,
                    std_dev=std,
                    deviation_sigma=dev,
                ))

    entries.sort(key=lambda e: (e.name, e.run_number))
    return entries


def compute_drift_entries(
    rd_curr: "RunData",
    rd_prev: "RunData",
    ref_idx: int,
    mod_by_name: Dict[str, "Module"],
    thresh_g: float = 0.10,
    thresh_w: float = 0.05,
) -> List[DriftEntry]:
    """Find modules where gain changed by more than threshold relative to previous run."""
    entries: List[DriftEntry] = []
    for mname, mrec in rd_curr.modules.items():
        prev_mrec = rd_prev.modules.get(mname)
        if prev_mrec is None:
            continue
        g_curr = mrec.gain_factors[ref_idx]
        g_prev = prev_mrec.gain_factors[ref_idx]
        rel_display = math.inf if g_prev == 0 else (g_curr - g_prev) / g_prev
        denom = min(abs(g_curr), abs(g_prev))
        rel_sym = math.inf if denom == 0 else abs(g_curr - g_prev) / denom
        threshold = thresh_g if mname.startswith("G") else thresh_w
        if denom == 0 or rel_sym > threshold:
            mod = mod_by_name.get(mname)
            entries.append(DriftEntry(
                name=mname,
                mod_type=mod.mod_type if mod else "?",
                run_number=rd_curr.run_number,
                prev_run_number=rd_prev.run_number,
                gain_current=g_curr,
                gain_prev=g_prev,
                rel_change=rel_display,
            ))
    entries.sort(key=lambda e: (0 if e.name.startswith("W") else 1,
                                math.isinf(e.rel_change),
                                -(abs(e.gain_current - e.gain_prev) / min(abs(e.gain_current), abs(e.gain_prev))
                                  if min(abs(e.gain_current), abs(e.gain_prev)) != 0 else 0)))
    return entries


def compute_summary(
    runs: List[RunData],
    ref_idx: int,
    mod_by_name: Dict[str, "Module"],
    threshold: float = 0.05,
) -> List[SummaryEntry]:
    """Count run-to-run drift events per module across all consecutive run pairs."""
    stats: Dict[str, List] = {}  # name -> [count, max_rel, max_run, max_prev_run]

    for i in range(1, len(runs)):
        rd_curr = runs[i]
        rd_prev = runs[i - 1]
        curr_run_num = rd_curr.run_number
        prev_run_num = rd_prev.run_number
        # pre-flatten previous run gains to avoid chained lookups in inner loop
        prev_gains: Dict[str, float] = {}
        for mname, mrec in rd_prev.modules.items():
            prev_gains[mname] = mrec.gain_factors[ref_idx]

        for mname, mrec in rd_curr.modules.items():
            g_prev = prev_gains.get(mname)
            if g_prev is None:
                continue
            rel = math.inf if g_prev == 0 else abs(mrec.gain_factors[ref_idx] - g_prev) / g_prev
            if rel > threshold:
                s = stats.get(mname)
                if s is None:
                    stats[mname] = [1, rel, curr_run_num, prev_run_num]
                else:
                    s[0] += 1
                    if rel > s[1]:
                        s[1] = rel
                        s[2] = curr_run_num
                        s[3] = prev_run_num

    entries: List[SummaryEntry] = []
    for mname, (count, max_rel, max_run, max_prev) in stats.items():
        mod = mod_by_name.get(mname)
        entries.append(SummaryEntry(
            name=mname,
            mod_type=mod.mod_type if mod else "?",
            drift_count=count,
            max_rel_change=max_rel,
            max_run=max_run,
            max_prev_run=max_prev,
        ))
    entries.sort(key=lambda e: (
        0 if e.name.startswith("W") else 1,
        0 if e.name.startswith("W") else int(math.isinf(e.max_rel_change)),
        -e.drift_count,
        0 if math.isinf(e.max_rel_change) else -e.max_rel_change,
    ))
    return entries


# ===========================================================================
#  HyCal Gain Map Widget
# ===========================================================================

class HyCalGainMapWidget(QWidget):
    """Colour-coded HyCal module map with zoom/pan and hover tooltips."""

    moduleHovered = pyqtSignal(str)
    moduleClicked = pyqtSignal(str)   # emits name, or "" to deselect
    paletteClicked = pyqtSignal()

    _SHRINK = 0.90
    _CLICK_THRESHOLD = 4

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(400, 400)

        self._modules: List[Module] = []
        self._values: Dict[str, float] = {}
        self._vmin = 0.0
        self._vmax = 1.0
        self._log_scale = False
        self._palette_idx = 0
        self._palette_override = None
        self._legend_mode: Optional[str] = None  # None | "drift" | "summary"
        self._hovered: Optional[str] = None
        self._selected: Optional[str] = None
        self._rects: Dict[str, QRectF] = {}
        self._rect_names_rev: List[str] = []
        self._geo_bounds: Tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0)
        self._cb_rect: Optional[QRectF] = None
        self._layout_dirty = True

        # zoom / pan
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._drag_last: Optional[QPointF] = None
        self._drag_origin: Optional[QPointF] = None
        self._dragging = False

        # overlay reset button (top-right corner)
        self._reset_btn = QPushButton("Reset", self)
        self._reset_btn.setFixedSize(52, 24)
        self._reset_btn.setStyleSheet(
            "QPushButton{background:rgba(22,27,34,200);color:#8b949e;"
            "border:1px solid #30363d;font:bold 9px Consolas;border-radius:3px;}"
            "QPushButton:hover{background:rgba(33,38,45,220);color:#c9d1d9;}")
        self._reset_btn.clicked.connect(self.reset_view)

    # -- public API --

    def set_modules(self, modules: List[Module]):
        self._modules = [m for m in modules if m.mod_type != "LMS"]
        if self._modules:
            self._geo_bounds = (
                min(m.x - m.sx / 2 for m in self._modules),
                max(m.x + m.sx / 2 for m in self._modules),
                min(m.y - m.sy / 2 for m in self._modules),
                max(m.y + m.sy / 2 for m in self._modules),
            )
        self._layout_dirty = True
        self.update()

    def set_gain_data(self, values: Dict[str, float],
                      vmin: float, vmax: float):
        self._values = values
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def set_log_scale(self, on: bool):
        self._log_scale = on
        self.update()

    def set_palette(self, idx: int):
        self._palette_idx = idx % len(PALETTES)
        self._palette_override = None
        self.update()

    def set_palette_override(self, stops):
        """Use a custom stops list instead of the indexed palette. Pass None to clear."""
        self._palette_override = stops
        self.update()

    def set_legend_mode(self, mode: Optional[str]):
        """Set legend overlay: None, 'drift', or 'summary'."""
        if mode != self._legend_mode:
            self._legend_mode = mode
            self.update()

    def set_selected(self, name: Optional[str]):
        self._selected = name
        self.update()

    def reset_view(self):
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._layout_dirty = True
        self.update()

    # -- layout --

    def _recompute_layout(self):
        self._rects.clear()
        if not self._modules:
            self._rect_names_rev = []
            return
        w, h = self.width(), self.height()
        margin, top, bot = 12, 8, 50
        pw, ph = w - 2 * margin, h - top - bot

        x0, x1, y0, y1 = self._geo_bounds

        base_scale = min(pw / (x1 - x0), ph / (y1 - y0))
        sc = base_scale * self._zoom
        dw, dh = (x1 - x0) * sc, (y1 - y0) * sc
        ox = margin + (pw - dw) / 2 + self._pan_x
        oy = top + (ph - dh) / 2 + self._pan_y
        shrink = self._SHRINK

        self._geo_x0 = x0
        self._geo_y1 = y1
        self._geo_sc = sc
        self._geo_ox = ox
        self._geo_oy = oy

        for m in self._modules:
            mw, mh = m.sx * sc * shrink, m.sy * sc * shrink
            cx = ox + (m.x - x0) * sc
            cy = oy + (y1 - m.y) * sc
            self._rects[m.name] = QRectF(cx - mw / 2, cy - mh / 2, mw, mh)
        self._rect_names_rev = list(self._rects)[::-1]
        self._layout_dirty = False

    def resizeEvent(self, event):
        self._layout_dirty = True
        # keep reset button in top-right corner
        self._reset_btn.move(self.width() - self._reset_btn.width() - 6, 6)
        super().resizeEvent(event)

    # -- painting --

    def paintEvent(self, event):
        if self._layout_dirty:
            self._recompute_layout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor("#0a0e14"))

        if not self._rects:
            if not self._values:
                p.setPen(QColor("#555555"))
                p.setFont(QFont("Consolas", 12))
                p.drawText(QRectF(0, 0, w, h),
                           Qt.AlignmentFlag.AlignCenter, "No data loaded")
            p.end()
            return

        if self._palette_override is not None:
            stops = self._palette_override
        else:
            stops = list(PALETTES.values())[self._palette_idx]
        vmin, vmax = self._vmin, self._vmax
        log_scale = self._log_scale
        no_data = QColor("#1a1a2e")

        # precompute log bounds
        if log_scale:
            log_lo = math.log10(max(vmin, 1e-9))
            log_hi = math.log10(max(vmax, vmin * 10, 1e-8))

        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is not None:
                if log_scale:
                    lv = math.log10(max(v, 1e-9))
                    t = (lv - log_lo) / (log_hi - log_lo) if log_hi > log_lo else 0.5
                else:
                    t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                p.fillRect(rect, _cmap_qcolor(t, stops))
            else:
                p.fillRect(rect, no_data)

        # selected highlight (white border)
        if self._selected and self._selected in self._rects:
            p.setPen(QPen(QColor("#ffffff"), 2.5))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._selected])

        # hover highlight
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(QColor("#58a6ff"), 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

        # colour bar
        cb_w = min(300, w - 80)
        cb_h = 14
        cb_x = (w - cb_w) / 2
        cb_y = h - 40
        self._cb_rect = QRectF(cb_x, cb_y, cb_w, cb_h)

        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in stops:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(self._cb_rect, QBrush(grad))
        p.setPen(QPen(QColor("#58a6ff"), 1.0))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(self._cb_rect)

        # range labels + palette name
        p.setPen(QColor("#8b949e"))
        p.setFont(QFont("Consolas", 9))
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, 80, 14),
                   Qt.AlignmentFlag.AlignLeft, f"{vmin:.4f}")
        p.drawText(QRectF(cb_x + cb_w - 80, cb_y + cb_h + 2, 80, 14),
                   Qt.AlignmentFlag.AlignRight, f"{vmax:.4f}")
        pname = "cyan-grey-orange" if self._palette_override is not None else PALETTE_NAMES[self._palette_idx]
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, cb_w, 14),
                   Qt.AlignmentFlag.AlignCenter, pname)

        # legend (just above the colour bar)
        if self._legend_mode == "drift":
            items = [
                (QColor(0, 210, 230),  "gain decreases"),
                (QColor(80, 80, 80),   "stable"),
                (QColor(249, 115, 22), "gain increases"),
            ]
        elif self._legend_mode == "summary":
            items = [
                (QColor(10, 42, 110),  "low drift count"),
                (QColor(249, 115, 22), "high drift count"),
            ]
        elif self._legend_mode == "gain":
            items = [
                (QColor(10, 42, 110),  "low gain"),
                (QColor(249, 115, 22), "high gain"),
            ]
        elif self._legend_mode == "deviation":
            items = [
                (QColor(10, 42, 110),  "below mean"),
                (QColor(80, 80, 80),   "near mean"),
                (QColor(249, 115, 22), "above mean"),
            ]
        else:
            items = []
        if items:
            p.setFont(QFont("Consolas", 9))
            fm = p.fontMetrics()
            swatch = 12
            gap = 5
            item_w = swatch + gap + max(fm.horizontalAdvance(lbl) for _, lbl in items)
            spacing = 18
            total_w = len(items) * item_w + (len(items) - 1) * spacing
            lh = max(swatch, fm.height())
            pad = 4
            lx = (w - total_w) // 2 - pad
            ly = cb_y - lh - 2 * pad - 4
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(10, 14, 20, 200))
            p.drawRoundedRect(QRectF(lx, ly, total_w + 2 * pad, lh + 2 * pad), 4, 4)
            x = lx + pad
            for color, label in items:
                sy = ly + pad + (lh - swatch) // 2
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(color)
                p.drawRect(QRectF(x, sy, swatch, swatch))
                p.setPen(QColor("#c9d1d9"))
                p.drawText(QRectF(x + swatch + gap, ly + pad, fm.horizontalAdvance(label), lh),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                           label)
                x += item_w + spacing

        p.end()

    # -- mouse events (zoom/pan from scan_geoview pattern) --

    def mousePressEvent(self, e):
        if e.button() == Qt.MouseButton.MiddleButton:
            self.reset_view()
            return
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            self._drag_last = e.position()
            self._drag_origin = e.position()
            self._dragging = False

    def mouseReleaseEvent(self, e):
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            if self._dragging:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            elif e.button() == Qt.MouseButton.LeftButton:
                pos = e.position()
                # click on colour bar -> cycle palette
                if self._cb_rect and self._cb_rect.contains(pos):
                    self.paletteClicked.emit()
                else:
                    # click on a module -> select/deselect
                    hit = None
                    for name in self._rect_names_rev:
                        if self._rects[name].contains(pos):
                            hit = name
                            break
                    if hit is not None:
                        new_sel = None if hit == self._selected else hit
                        self._selected = new_sel
                        self.update()
                        self.moduleClicked.emit(new_sel if new_sel else "")
                    elif self._selected is not None:
                        self._selected = None
                        self.update()
                        self.moduleClicked.emit("")
            self._drag_last = None
            self._drag_origin = None
            self._dragging = False

    def mouseMoveEvent(self, e):
        # drag
        if self._drag_last is not None:
            pos = e.position()
            if not self._dragging:
                dx = pos.x() - self._drag_origin.x()
                dy = pos.y() - self._drag_origin.y()
                if dx * dx + dy * dy > self._CLICK_THRESHOLD ** 2:
                    self._dragging = True
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
            if self._dragging:
                self._pan_x += pos.x() - self._drag_last.x()
                self._pan_y += pos.y() - self._drag_last.y()
                self._drag_last = pos
                self._layout_dirty = True
                self.update()
            return

        # hover tooltip
        pos = e.position()
        # hand cursor over colour bar
        if self._cb_rect and self._cb_rect.contains(pos):
            self.setCursor(Qt.CursorShape.PointingHandCursor)
        else:
            self.setCursor(Qt.CursorShape.ArrowCursor)

        found = None
        for name in self._rect_names_rev:
            if self._rects[name].contains(pos):
                found = name
                break
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                v = self._values.get(found)
                tip = f"{found}: {v:.5f}" if v is not None else found
                QToolTip.showText(e.globalPosition().toPoint(), tip, self)
                self.moduleHovered.emit(found)
            else:
                QToolTip.hideText()

    def wheelEvent(self, e):
        factor = 1.15 if e.angleDelta().y() > 0 else 1.0 / 1.15
        new_zoom = max(0.5, min(self._zoom * factor, 20.0))
        if new_zoom == self._zoom:
            return
        pos = e.position()
        ratio = new_zoom / self._zoom
        self._pan_x = pos.x() + (self._pan_x - pos.x()) * ratio
        self._pan_y = pos.y() + (self._pan_y - pos.y()) * ratio
        self._zoom = new_zoom
        self._layout_dirty = True
        self.update()

    def sizeHint(self):
        return QSize(680, 680)


# ===========================================================================
#  LMS Line Chart Widget
# ===========================================================================

class LMSLineChartWidget(QWidget):
    """Line chart with error bars for LMS peak/alpha ratio vs run number."""

    PAD_L, PAD_R, PAD_T, PAD_B = 60, 16, 24, 32

    runClicked = pyqtSignal(int)   # emits actual run number when a point is clicked

    def __init__(self, parent=None):
        super().__init__(parent)
        self._run_numbers: List[int] = []
        self._actual_run_numbers: List[int] = []
        self._ratios: List[float] = []
        self._errors: List[float] = []
        self._title: str = ""
        self._hover_idx: int = -1
        self._highlighted: bool = False
        self._current_run_number: int = -1
        self.setMinimumHeight(100)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

    def set_highlighted(self, on: bool):
        if on != self._highlighted:
            self._highlighted = on
            self.update()

    def set_current_run(self, run_number: int):
        if run_number != self._current_run_number:
            self._current_run_number = run_number
            self.update()

    def set_data(self, run_numbers: List[int], ratios: List[float],
                 errors: List[float], title: str,
                 actual_run_numbers: List[int] = None):
        self._run_numbers = run_numbers
        self._actual_run_numbers = actual_run_numbers if actual_run_numbers is not None else run_numbers
        self._ratios = ratios
        self._errors = errors
        self._title = title
        self._hover_idx = -1
        self.update()

    def _screen_xs(self, w: int) -> List[float]:
        """Return screen x-coordinates for all data points."""
        runs = self._run_numbers
        if not runs:
            return []
        px = self.PAD_L
        pw = w - self.PAD_L - self.PAD_R
        x_min, x_max = runs[0], runs[-1]
        if x_min == x_max:
            x_min -= 1; x_max += 1
        return [px + (r - x_min) / (x_max - x_min) * pw for r in runs]

    def mouseMoveEvent(self, event):
        runs = self._run_numbers
        if not runs:
            return
        sx_list = self._screen_xs(self.width())
        mx = event.position().x()
        best_i, best_d = -1, float("inf")
        for i, sx in enumerate(sx_list):
            d = abs(mx - sx)
            if d < best_d:
                best_d = d
                best_i = i
        new_idx = best_i if best_d < 20 else -1
        if new_idx != self._hover_idx:
            self._hover_idx = new_idx
            self.update()

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        runs = self._run_numbers
        if not runs:
            return
        sx_list = self._screen_xs(self.width())
        mx = event.position().x()
        best_i, best_d = -1, float("inf")
        for i, sx in enumerate(sx_list):
            d = abs(mx - sx)
            if d < best_d:
                best_d = d
                best_i = i
        if best_d < 20 and best_i < len(self._actual_run_numbers):
            self.runClicked.emit(self._actual_run_numbers[best_i])

    def leaveEvent(self, event):
        if self._hover_idx != -1:
            self._hover_idx = -1
            self.update()

    @staticmethod
    def _nice_ticks(lo: float, hi: float, max_ticks: int = 6):
        """Compute nice tick values for an axis range."""
        if hi <= lo:
            return [lo]
        raw = (hi - lo) / max(max_ticks - 1, 1)
        mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
        candidates = [1, 2, 2.5, 5, 10]
        step = mag
        for c in candidates:
            if c * mag >= raw:
                step = c * mag
                break
        start = math.ceil(lo / step) * step
        ticks = []
        v = start
        while v <= hi + step * 0.01:
            ticks.append(v)
            v += step
        return ticks

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor("#0a0e14"))

        # title
        title_color = QColor("#f97316") if self._highlighted else QColor("#58a6ff")
        title_text = (self._title + "  [reference]") if self._highlighted else self._title
        p.setPen(title_color)
        p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
        p.drawText(QRectF(self.PAD_L, 2, w - self.PAD_L - self.PAD_R, 20),
                   Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                   title_text)

        runs = self._run_numbers
        ratios = self._ratios
        errors = self._errors
        if not runs or not ratios:
            p.setPen(QColor("#555555"))
            p.setFont(QFont("Consolas", 10))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter, "No data")
            p.end()
            return

        # plot area
        px = self.PAD_L
        py = self.PAD_T
        pw = w - self.PAD_L - self.PAD_R
        ph = h - self.PAD_T - self.PAD_B
        if pw < 20 or ph < 20:
            p.end()
            return

        # data ranges
        x_min, x_max = runs[0], runs[-1]
        if x_min == x_max:
            x_min -= 1
            x_max += 1

        y_vals = []
        for i, r in enumerate(ratios):
            y_vals.append(r)
            if i < len(errors):
                y_vals.append(r + errors[i])
                y_vals.append(r - errors[i])
        y_lo = min(y_vals)
        y_hi = max(y_vals)
        margin = (y_hi - y_lo) * 0.1 if y_hi > y_lo else 0.05
        y_lo -= margin
        y_hi += margin

        def to_sx(v):
            return px + (v - x_min) / (x_max - x_min) * pw

        def to_sy(v):
            return py + ph - (v - y_lo) / (y_hi - y_lo) * ph

        # grid + axes
        p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
        y_ticks = self._nice_ticks(y_lo, y_hi, 5)
        for yt in y_ticks:
            sy = to_sy(yt)
            p.drawLine(QPointF(px, sy), QPointF(px + pw, sy))

        # axes border
        p.setPen(QPen(QColor("#30363d"), 1))
        p.drawLine(QPointF(px, py), QPointF(px, py + ph))
        p.drawLine(QPointF(px, py + ph), QPointF(px + pw, py + ph))

        # y-axis labels
        p.setPen(QColor("#8b949e"))
        p.setFont(QFont("Consolas", 8))
        for yt in y_ticks:
            sy = to_sy(yt)
            p.drawText(QRectF(0, sy - 8, self.PAD_L - 4, 16),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       f"{yt:.3f}")

        # x-axis labels (skip some if too dense)
        max_labels = max(pw // 60, 2)
        step = max(len(runs) // max_labels, 1)
        for i in range(0, len(runs), step):
            sx = to_sx(runs[i])
            p.drawText(QRectF(sx - 30, py + ph + 2, 60, 18),
                       Qt.AlignmentFlag.AlignCenter, str(runs[i]))

        # error bars
        p.setPen(QPen(QColor("#8b949e"), 1))
        cap = 3
        for i in range(len(runs)):
            if i >= len(ratios):
                break
            sx = to_sx(runs[i])
            r = ratios[i]
            err = errors[i] if i < len(errors) else 0
            sy_top = to_sy(r + err)
            sy_bot = to_sy(r - err)
            p.drawLine(QPointF(sx, sy_top), QPointF(sx, sy_bot))
            p.drawLine(QPointF(sx - cap, sy_top), QPointF(sx + cap, sy_top))
            p.drawLine(QPointF(sx - cap, sy_bot), QPointF(sx + cap, sy_bot))

        series_color = QColor("#f97316") if self._highlighted else QColor("#58a6ff")

        # connecting line
        p.setPen(QPen(series_color, 1.5))
        for i in range(len(runs) - 1):
            if i + 1 >= len(ratios):
                break
            p.drawLine(QPointF(to_sx(runs[i]), to_sy(ratios[i])),
                       QPointF(to_sx(runs[i + 1]), to_sy(ratios[i + 1])))

        # points
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(series_color)
        for i in range(len(runs)):
            if i >= len(ratios):
                break
            radius = 5 if i == self._hover_idx else 3
            p.drawEllipse(QPointF(to_sx(runs[i]), to_sy(ratios[i])), radius, radius)

        # current run red circle
        actual = self._actual_run_numbers
        if self._current_run_number >= 0 and actual:
            for i, rn in enumerate(actual):
                if rn == self._current_run_number and i < len(runs) and i < len(ratios):
                    cx = to_sx(runs[i])
                    cy = to_sy(ratios[i])
                    p.setPen(QPen(QColor("#ff2222"), 2))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawEllipse(QPointF(cx, cy), 8, 8)
                    break

        # hover tooltip
        hi = self._hover_idx
        if 0 <= hi < len(runs) and hi < len(ratios):
            sx = to_sx(runs[hi])
            sy = to_sy(ratios[hi])
            actual_rn = self._actual_run_numbers[hi] if hi < len(self._actual_run_numbers) else runs[hi]
            tip = f"run {actual_rn}\n{ratios[hi]:.4f}"
            p.setFont(QFont("Consolas", 9))
            fm = p.fontMetrics()
            lines = tip.split("\n")
            tw = max(fm.horizontalAdvance(ln) for ln in lines)
            th = fm.height() * len(lines) + 6
            tx = sx + 10
            ty = sy - th // 2
            if tx + tw + 8 > px + pw:
                tx = sx - tw - 14
            if ty < py:
                ty = py
            if ty + th > py + ph:
                ty = py + ph - th
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(20, 30, 45, 220))
            p.drawRoundedRect(QRectF(tx - 4, ty - 2, tw + 8, th), 4, 4)
            p.setPen(QColor("#e6edf3"))
            for j, ln in enumerate(lines):
                p.drawText(QRectF(tx, ty + j * fm.height(), tw, fm.height()),
                           Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                           ln)

        p.end()


# ===========================================================================
#  Irregular channels table
# ===========================================================================

class IrregularTableWidget(QWidget):
    """Table of outlier entries — supports both irregular (deviation) and drift modes."""

    runClicked = pyqtSignal(int)    # emits run number when a row is clicked
    moduleClicked = pyqtSignal(str) # emits module name when a row is clicked

    def __init__(self, parent=None):
        super().__init__(parent)
        self._entries: List[IrregularEntry] = []
        self._drift_entries: List[DriftEntry] = []
        self._summary_entries: List[SummaryEntry] = []
        self._mode: str = "irregular"
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        # filter bar
        fbar = QHBoxLayout()
        fbar.setSpacing(6)

        lbl = QLabel("Search:")
        lbl.setFont(QFont("Consolas", 10))
        lbl.setStyleSheet("color:#c9d1d9;")
        fbar.addWidget(lbl)

        self._search = QLineEdit()
        self._search.setPlaceholderText("module name...")
        self._search.setFixedWidth(120)
        self._search.setFont(QFont("Consolas", 10))
        self._search.setStyleSheet(
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 4px;}")
        self._search.textChanged.connect(self._apply_filter)
        fbar.addWidget(self._search)

        lbl2 = QLabel("Type:")
        lbl2.setFont(QFont("Consolas", 10))
        lbl2.setStyleSheet("color:#c9d1d9;")
        fbar.addWidget(lbl2)

        self._type_filter = QComboBox()
        self._type_filter.addItems(["All", "PbWO4", "PbGlass"])
        self._type_filter.setFixedWidth(100)
        self._type_filter.setFont(QFont("Consolas", 10))
        self._type_filter.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}")
        self._type_filter.currentIndexChanged.connect(
            lambda _: self._apply_filter())
        fbar.addWidget(self._type_filter)

        fbar.addStretch()

        count_lbl = QLabel("")
        count_lbl.setFont(QFont("Consolas", 10))
        count_lbl.setStyleSheet("color:#8b949e;")
        self._count_lbl = count_lbl
        fbar.addWidget(count_lbl)

        layout.addLayout(fbar)

        # table
        self._table = QTableWidget()
        self._table.setColumnCount(6)
        self._table.setHorizontalHeaderLabels(
            ["Module", "Run", "Gain", "Mean", "Std Dev", "Dev (sigma)"])
        self._table.setFont(QFont("Consolas", 10))
        self._table.setStyleSheet(
            "QTableWidget{background:#0d1117;color:#c9d1d9;"
            "gridline-color:#21262d;border:1px solid #30363d;}"
            "QTableWidget::item{padding:2px 6px;}"
            "QHeaderView::section{background:#161b22;color:#58a6ff;"
            "border:1px solid #30363d;font:bold 10pt Consolas;padding:4px;}")
        self._table.setAlternatingRowColors(True)
        pal = self._table.palette()
        pal.setColor(QPalette.ColorRole.AlternateBase, QColor("#131820"))
        self._table.setPalette(pal)
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSortingEnabled(False)
        self._table.cellClicked.connect(self._on_cell_clicked)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(24)

        hdr = self._table.horizontalHeader()
        hdr.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self._table.setColumnWidth(0, 90)
        self._table.setColumnWidth(1, 75)

        layout.addWidget(self._table)

    def set_data(self, entries: List[IrregularEntry]):
        self._entries = entries
        if self._mode != "irregular":
            self._mode = "irregular"
            self._table.setColumnCount(6)
            self._table.setHorizontalHeaderLabels(
                ["Module", "Run", "Gain", "Mean", "Std Dev", "Dev (σ)"])
        else:
            self._mode = "irregular"
        self._apply_filter()

    def set_drift_data(self, entries: List[DriftEntry]):
        self._drift_entries = entries
        if self._mode != "drift":
            self._mode = "drift"
            self._table.setColumnCount(6)
            self._table.setHorizontalHeaderLabels(
                ["Module", "Curr Run", "Prev Run", "Gain (curr)", "Gain (prev)", "Δ (%)"])
        else:
            self._mode = "drift"
        self._apply_filter()

    def set_summary_data(self, entries: List[SummaryEntry]):
        self._summary_entries = entries
        if self._mode != "summary":
            self._mode = "summary"
            self._table.setColumnCount(5)
            self._table.setHorizontalHeaderLabels(
                ["Module", "Drift Counts", "Max |Δ%|", "Worst Run", "Prev Run"])
        else:
            self._mode = "summary"
        self._apply_filter()

    def _apply_filter(self):
        search = self._search.text().strip().upper()
        type_sel = self._type_filter.currentText()

        if self._mode == "drift":
            source = self._drift_entries
        elif self._mode == "summary":
            source = self._summary_entries
        else:
            source = self._entries

        filtered = []
        for e in source:
            if search and search not in e.name.upper():
                continue
            if type_sel == "PbWO4" and e.mod_type != "PbWO4":
                continue
            if type_sel == "PbGlass" and e.mod_type != "PbGlass":
                continue
            filtered.append(e)

        if self._mode == "drift":
            self._count_lbl.setText(f"drifted channels: {len(filtered)} entries")
            self._populate_drift_table(filtered)
        elif self._mode == "summary":
            self._count_lbl.setText(f"problematic channels: {len(filtered)}")
            self._populate_summary_table(filtered)
        else:
            self._count_lbl.setText(f"irregular gains: {len(filtered)} entries")
            self._populate_table(filtered)

    def _populate_table(self, entries: List[IrregularEntry]):
        self._table.setRowCount(len(entries))
        for row, e in enumerate(entries):
            item = QTableWidgetItem(e.name)
            item.setForeground(QColor("#d29922"))
            self._table.setItem(row, 0, item)

            item_run = QTableWidgetItem()
            item_run.setData(Qt.ItemDataRole.DisplayRole, e.run_number)
            self._table.setItem(row, 1, item_run)

            for col, val in [(2, e.gain), (3, e.mean_gain),
                             (4, e.std_dev), (5, e.deviation_sigma)]:
                item_f = QTableWidgetItem()
                item_f.setData(Qt.ItemDataRole.DisplayRole, round(val, 5))
                self._table.setItem(row, col, item_f)


    def _populate_drift_table(self, entries: List[DriftEntry]):
        self._table.setRowCount(len(entries))
        for row, e in enumerate(entries):
            item = QTableWidgetItem(e.name)
            color = QColor("#f85149") if e.rel_change < 0 else QColor("#3fb950")
            item.setForeground(color)
            self._table.setItem(row, 0, item)

            item_curr = QTableWidgetItem()
            item_curr.setData(Qt.ItemDataRole.DisplayRole, e.run_number)
            self._table.setItem(row, 1, item_curr)

            item_prev = QTableWidgetItem()
            item_prev.setData(Qt.ItemDataRole.DisplayRole, e.prev_run_number)
            self._table.setItem(row, 2, item_prev)

            for col, val in [(3, e.gain_current), (4, e.gain_prev)]:
                item_f = QTableWidgetItem()
                item_f.setData(Qt.ItemDataRole.DisplayRole, round(val, 5))
                self._table.setItem(row, col, item_f)

            item_pct = QTableWidgetItem()
            item_pct.setData(Qt.ItemDataRole.DisplayRole, round(e.rel_change * 100, 3))
            self._table.setItem(row, 5, item_pct)


    def _populate_summary_table(self, entries: List[SummaryEntry]):
        self._table.setRowCount(len(entries))
        for row, e in enumerate(entries):
            item = QTableWidgetItem(e.name)
            item.setForeground(QColor("#f85149"))
            self._table.setItem(row, 0, item)

            item_cnt = QTableWidgetItem()
            item_cnt.setData(Qt.ItemDataRole.DisplayRole, e.drift_count)
            self._table.setItem(row, 1, item_cnt)

            item_max = QTableWidgetItem()
            item_max.setData(Qt.ItemDataRole.DisplayRole, round(e.max_rel_change * 100, 3))
            self._table.setItem(row, 2, item_max)

            item_run = QTableWidgetItem()
            item_run.setData(Qt.ItemDataRole.DisplayRole, e.max_run)
            self._table.setItem(row, 3, item_run)

            item_prev = QTableWidgetItem()
            item_prev.setData(Qt.ItemDataRole.DisplayRole, e.max_prev_run)
            self._table.setItem(row, 4, item_prev)


    def select_module(self, name: str):
        """Highlight and scroll to the first row matching name, or clear selection."""
        for row in range(self._table.rowCount()):
            item = self._table.item(row, 0)
            if item is not None and item.text() == name:
                self._table.selectRow(row)
                self._table.scrollToItem(item)
                return
        self._table.clearSelection()

    def _on_cell_clicked(self, row: int, _col: int):
        name_item = self._table.item(row, 0)
        if name_item is not None:
            self.moduleClicked.emit(name_item.text())
        item = self._table.item(row, 1)  # Run column (Curr Run for drift, Run for irregular)
        if item is not None:
            run_num = item.data(Qt.ItemDataRole.DisplayRole)
            if isinstance(run_num, int):
                self.runClicked.emit(run_num)


# ===========================================================================
#  Main window
# ===========================================================================

class GainMonitorWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self._all_modules: List[Module] = []
        self._mod_by_name: Dict[str, Module] = {}
        self._runs: List[RunData] = []
        self._current_run_idx: int = 0
        self._current_ref_idx: int = LMS_REF_DEFAULT
        self._palette_idx = 0
        self._auto_range = True
        self._manual_vmin = 0.9
        self._manual_vmax = 1.1
        self._log_scale = False
        self._selected_module: Optional[str] = None
        self._start_run_idx: int = 0
        self._end_run_idx: int = 0
        self._thresh_g: float = 0.10
        self._thresh_w: float = 0.05
        self._view_mode: int = 2   # 0 = Gain Factor, 1 = Deviation (σ), 2 = Run-to-Run Drift, 3 = Summary
        # pre-computed pairwise diffs for summary mode: (name, mod_type, rel, pair_idx)
        # pair_idx = index of curr run in self._runs; recomputed on load or ref change
        self._pairwise_diffs: List = []
        self._pairwise_ref_idx: int = -1
        self._current_folder: str = ""
        self._deviation_stats_cache: Optional[Dict] = None
        self._deviation_stats_key: Optional[Tuple] = None
        self._file_snapshot: Dict[int, float] = {}   # run_number -> mtime
        self._auto_refresh_timer = QTimer(self)
        self._auto_refresh_timer.timeout.connect(self._auto_refresh_check)

        self._load_geometry()
        self._build_ui()
        legend_map = {0: "gain", 1: "deviation", 2: "drift", 3: "summary"}
        self._map.set_legend_mode(legend_map.get(self._view_mode))
        if self._view_mode == 2:
            self._map.set_palette_override(DRIFT_PALETTE)

    def _load_geometry(self):
        self._all_modules = load_modules()
        self._mod_by_name = {m.name: m for m in self._all_modules}

    @property
    def _active_runs(self) -> List[RunData]:
        return self._runs[self._start_run_idx:self._end_run_idx + 1]

    # ---- UI ----

    def _build_ui(self):
        self.setWindowTitle("HyCal Gain Monitor")
        self.resize(1600, 1000)
        self._apply_dark_palette()

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL GAIN MONITOR")
        lbl.setFont(QFont("Consolas", 14, QFont.Weight.Bold))
        lbl.setStyleSheet("color:#58a6ff;")
        top.addWidget(lbl)
        top.addStretch()

        self._process_btn = self._make_btn(
            "Process Folder...", "#58a6ff", self._on_process_folder)
        top.addWidget(self._process_btn)

        self._refresh_btn = self._make_btn(
            "Refresh", "#3fb950", self._on_refresh)
        self._refresh_btn.setEnabled(False)
        top.addWidget(self._refresh_btn)

        self._auto_refresh_btn = QPushButton("Auto-Refresh: ON")
        self._auto_refresh_btn.setCheckable(True)
        self._auto_refresh_btn.setChecked(True)
        self._auto_refresh_btn.setEnabled(False)
        self._auto_refresh_btn.setFont(QFont("Consolas", 10))
        self._auto_refresh_btn.setFixedHeight(28)
        self._auto_refresh_btn.setStyleSheet(
            "QPushButton{background:#161b22;color:#8b949e;border:1px solid #30363d;"
            "border-radius:3px;padding:0 8px;}"
            "QPushButton:checked{background:#1a3a1a;color:#3fb950;border-color:#3fb950;}"
            "QPushButton:hover{background:#21262d;}")
        self._auto_refresh_btn.toggled.connect(self._on_auto_refresh_toggled)
        top.addWidget(self._auto_refresh_btn)

        top.addWidget(self._slabel("every"))
        self._auto_refresh_interval = QSpinBox()
        self._auto_refresh_interval.setRange(5, 3600)
        self._auto_refresh_interval.setValue(30)
        self._auto_refresh_interval.setSuffix(" s")
        self._auto_refresh_interval.setFixedWidth(72)
        self._auto_refresh_interval.setFont(QFont("Consolas", 10))
        self._auto_refresh_interval.setStyleSheet(
            "QSpinBox{background:#161b22;color:#c9d1d9;border:1px solid #30363d;"
            "border-radius:3px;padding:2px 4px;}"
            "QSpinBox::up-button,QSpinBox::down-button{width:16px;}")
        self._auto_refresh_interval.valueChanged.connect(self._on_refresh_interval_changed)
        top.addWidget(self._auto_refresh_interval)

        self._status_lbl = QLabel("No data loaded")
        self._status_lbl.setFont(QFont("Consolas", 11))
        self._status_lbl.setStyleSheet("color:#8b949e;")
        top.addWidget(self._status_lbl)
        root.addLayout(top)

        # -- body splitter (horizontal: left=map, right=charts+table) --
        body = QSplitter(Qt.Orientation.Horizontal)

        # ---- left panel ----
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(4)

        # controls
        ctrl = QHBoxLayout()
        ctrl.setSpacing(6)

        ctrl.addWidget(self._slabel("Ref:"))
        self._ref_combo = QComboBox()
        self._ref_combo.addItems(LMS_NAMES)
        self._ref_combo.setCurrentIndex(LMS_REF_DEFAULT)
        self._ref_combo.setFixedWidth(90)
        self._ref_combo.setFont(QFont("Consolas", 10))
        self._ref_combo.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}")
        self._ref_combo.currentIndexChanged.connect(self._on_ref_changed)
        ctrl.addWidget(self._ref_combo)

        ctrl.addSpacing(10)
        ctrl.addWidget(self._slabel("View:"))
        self._view_combo = QComboBox()
        self._view_combo.addItems(["Gain Factor", "Deviation (σ)", "Run-to-Run Drift", "Summary"])
        self._view_combo.setFixedWidth(150)
        self._view_combo.setFont(QFont("Consolas", 10))
        self._view_combo.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}")
        self._view_combo.currentIndexChanged.connect(self._on_view_mode_changed)
        ctrl.addWidget(self._view_combo)

        self._view_combo.blockSignals(True)
        self._view_combo.setCurrentIndex(2)
        self._view_combo.blockSignals(False)

        _edit_ss = ("QLineEdit{background:#161b22;color:#c9d1d9;"
                    "border:1px solid #30363d;border-radius:3px;padding:2px 4px;}")

        self._thresh_lbl = self._slabel("G thresh:")
        ctrl.addSpacing(6)
        ctrl.addWidget(self._thresh_lbl)

        self._thresh_g_input = QLineEdit("10.0")
        self._thresh_g_input.setFixedWidth(46)
        self._thresh_g_input.setFont(QFont("Consolas", 10))
        self._thresh_g_input.setStyleSheet(_edit_ss)
        self._thresh_g_input.editingFinished.connect(self._on_drift_threshold_changed)
        ctrl.addWidget(self._thresh_g_input)
        self._thresh_g_pct = self._slabel("%")
        ctrl.addWidget(self._thresh_g_pct)

        self._thresh_w_lbl = self._slabel("  W thresh:")
        ctrl.addWidget(self._thresh_w_lbl)

        self._thresh_w_input = QLineEdit("5.0")
        self._thresh_w_input.setFixedWidth(46)
        self._thresh_w_input.setFont(QFont("Consolas", 10))
        self._thresh_w_input.setStyleSheet(_edit_ss)
        self._thresh_w_input.editingFinished.connect(self._on_drift_threshold_changed)
        ctrl.addWidget(self._thresh_w_input)
        self._thresh_w_pct = self._slabel("%")
        ctrl.addWidget(self._thresh_w_pct)

        ctrl.addSpacing(10)
        ctrl.addWidget(self._slabel("Start:"))
        self._start_combo = QComboBox()
        self._start_combo.setMinimumWidth(100)
        self._start_combo.setFont(QFont("Consolas", 10))
        self._start_combo.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}")
        self._start_combo.currentIndexChanged.connect(self._on_start_run_changed)
        ctrl.addWidget(self._start_combo)

        ctrl.addSpacing(6)
        ctrl.addWidget(self._slabel("End:"))
        self._end_combo = QComboBox()
        self._end_combo.setMinimumWidth(100)
        self._end_combo.setFont(QFont("Consolas", 10))
        self._end_combo.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}")
        self._end_combo.currentIndexChanged.connect(self._on_end_run_changed)
        ctrl.addWidget(self._end_combo)

        ctrl.addSpacing(10)
        ctrl.addWidget(self._slabel("Run:"))

        self._run_combo = QComboBox()
        self._run_combo.setMinimumWidth(100)
        self._run_combo.setFont(QFont("Consolas", 10))
        self._run_combo.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}")
        self._run_combo.currentIndexChanged.connect(self._on_run_changed)
        ctrl.addWidget(self._run_combo)

        self._prev_btn = self._make_btn("<", "#c9d1d9", self._on_prev_run)
        self._prev_btn.setFixedWidth(30)
        ctrl.addWidget(self._prev_btn)

        self._next_btn = self._make_btn(">", "#c9d1d9", self._on_next_run)
        self._next_btn.setFixedWidth(30)
        ctrl.addWidget(self._next_btn)

        ctrl.addStretch()
        left_layout.addLayout(ctrl)

        # range controls
        rng = QHBoxLayout()
        rng.setSpacing(6)

        _EDIT_SS = ("QLineEdit{background:#161b22;color:#c9d1d9;"
                    "border:1px solid #30363d;border-radius:3px;padding:2px 4px;}")

        rng.addWidget(self._slabel("Min:"))
        self._range_min = QLineEdit("0.9")
        self._range_min.setFixedWidth(70)
        self._range_min.setFont(QFont("Consolas", 10))
        self._range_min.setStyleSheet(_EDIT_SS)
        self._range_min.returnPressed.connect(self._on_apply_range)
        rng.addWidget(self._range_min)

        rng.addWidget(self._slabel("Max:"))
        self._range_max = QLineEdit("1.1")
        self._range_max.setFixedWidth(70)
        self._range_max.setFont(QFont("Consolas", 10))
        self._range_max.setStyleSheet(_EDIT_SS)
        self._range_max.returnPressed.connect(self._on_apply_range)
        rng.addWidget(self._range_max)

        self._apply_btn = QPushButton("Apply")
        self._apply_btn.setFixedWidth(55)
        self._apply_btn.clicked.connect(self._on_apply_range)
        rng.addWidget(self._apply_btn)

        self._log_btn = QPushButton("Log")
        self._log_btn.setFixedWidth(45)
        self._log_btn.setCheckable(True)
        self._log_btn.clicked.connect(self._on_log_toggled)
        rng.addWidget(self._log_btn)

        self._auto_btn = QPushButton("Auto")
        self._auto_btn.setFixedWidth(55)
        self._auto_btn.setCheckable(True)
        self._auto_btn.setChecked(True)
        self._auto_btn.clicked.connect(self._on_auto_range)
        rng.addWidget(self._auto_btn)

        # common toggle-button style
        _TOGGLE_SS = (
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 8px;"
            "font:bold 11px Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#30363d;}"
            "QPushButton:checked{background:#1f6feb;color:white;"
            "border-color:#388bfd;}")
        self._apply_btn.setStyleSheet(
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 8px;"
            "font:bold 11px Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#30363d;}")
        self._log_btn.setStyleSheet(_TOGGLE_SS)
        self._auto_btn.setStyleSheet(_TOGGLE_SS)

        rng.addStretch()
        left_layout.addLayout(rng)

        # geo map
        self._map = HyCalGainMapWidget()
        self._map.set_modules(self._all_modules)
        self._map.moduleHovered.connect(self._on_module_hovered)
        self._map.moduleClicked.connect(self._on_module_clicked)
        self._map.paletteClicked.connect(self._on_cycle_palette)
        left_layout.addWidget(self._map, stretch=1)

        # info label
        self._info = QLabel("Hover over a module for details")
        self._info.setFont(QFont("Consolas", 10))
        self._info.setStyleSheet(
            "QLabel{background:#161b22;color:#c9d1d9;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}")
        self._info.setFixedHeight(26)
        left_layout.addWidget(self._info)

        body.addWidget(left)

        # ---- right panel (vertical splitter: charts + table) ----
        right = QSplitter(Qt.Orientation.Vertical)

        charts = QWidget()
        charts_layout = QVBoxLayout(charts)
        charts_layout.setContentsMargins(0, 0, 0, 0)
        charts_layout.setSpacing(2)

        self._charts: List[LMSLineChartWidget] = []
        for i, name in enumerate(LMS_NAMES):
            chart = LMSLineChartWidget()
            chart.set_data([], [], [], f"{name} (lms/alpha)")
            chart.runClicked.connect(self._on_chart_run_clicked)
            self._charts.append(chart)
            charts_layout.addWidget(chart)

        right.addWidget(charts)

        self._irregular_table = IrregularTableWidget()
        self._irregular_table.runClicked.connect(self._on_jump_to_run)
        self._irregular_table.moduleClicked.connect(self._on_module_clicked)
        right.addWidget(self._irregular_table)

        right.setStretchFactor(0, 3)
        right.setStretchFactor(1, 2)

        body.addWidget(right)
        body.setStretchFactor(0, 1)
        body.setStretchFactor(1, 1)

        root.addWidget(body, stretch=1)

    # ---- helpers ----

    def _make_btn(self, text, fg, slot):
        btn = QPushButton(text)
        btn.setStyleSheet(
            f"QPushButton{{background:#21262d;color:{fg};"
            f"border:1px solid #30363d;padding:4px 12px;"
            f"font:bold 11px Consolas;border-radius:3px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{color:#555;}}")
        btn.clicked.connect(slot)
        return btn

    def _slabel(self, text):
        lbl = QLabel(text)
        lbl.setFont(QFont("Consolas", 10))
        lbl.setStyleSheet("color:#c9d1d9;")
        return lbl

    def _apply_dark_palette(self):
        pal = self.palette()
        for role, colour in [
            (QPalette.ColorRole.Window, "#0d1117"),
            (QPalette.ColorRole.WindowText, "#c9d1d9"),
            (QPalette.ColorRole.Base, "#161b22"),
            (QPalette.ColorRole.Text, "#c9d1d9"),
            (QPalette.ColorRole.Button, "#21262d"),
            (QPalette.ColorRole.ButtonText, "#c9d1d9"),
            (QPalette.ColorRole.Highlight, "#58a6ff"),
        ]:
            pal.setColor(role, QColor(colour))
        self.setPalette(pal)

    # ---- keyboard navigation ----

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Left:
            self._on_prev_run()
        elif event.key() == Qt.Key.Key_Right:
            self._on_next_run()
        else:
            super().keyPressEvent(event)

    # ---- slots ----

    def _on_refresh(self):
        if self._current_folder:
            self._load_folder(self._current_folder)

    def _on_auto_refresh_toggled(self, checked: bool):
        self._auto_refresh_btn.setText("Auto-Refresh: ON" if checked else "Auto-Refresh: OFF")
        if checked:
            interval_ms = self._auto_refresh_interval.value() * 1000
            self._auto_refresh_timer.start(interval_ms)
        else:
            self._auto_refresh_timer.stop()

    def _on_refresh_interval_changed(self, value: int):
        if self._auto_refresh_timer.isActive():
            self._auto_refresh_timer.start(value * 1000)

    @staticmethod
    def _take_file_snapshot(folder: str) -> Dict[int, float]:
        """Return {run_number: mtime} for all LMS dat files in folder."""
        snapshot = {}
        for path in Path(folder).glob("prad_*_LMS.dat"):
            m = FILE_PATTERN.match(path.name)
            if m:
                snapshot[int(m.group(1))] = path.stat().st_mtime
        return snapshot

    def _auto_refresh_check(self):
        if not self._current_folder:
            return
        new_snapshot = self._take_file_snapshot(self._current_folder)
        if new_snapshot != self._file_snapshot:
            self._smart_refresh(new_snapshot)

    def _smart_refresh(self, new_snapshot: Dict[int, float]):
        """Reload data while preserving the current run selection and range."""
        current_run_number = (self._runs[self._current_run_idx].run_number
                              if self._runs else None)
        start_run_number = (self._runs[self._start_run_idx].run_number
                            if self._runs else None)
        end_run_number = (self._runs[self._end_run_idx].run_number
                          if self._runs else None)
        old_last_run_number = (self._runs[-1].run_number if self._runs else None)

        new_runs = load_all_runs(self._current_folder)
        if not new_runs:
            return
        self._runs = new_runs
        self._file_snapshot = new_snapshot
        self._deviation_stats_key = None

        # rebuild start/end combos
        self._start_combo.blockSignals(True)
        self._end_combo.blockSignals(True)
        self._start_combo.clear()
        self._end_combo.clear()
        for rd in self._runs:
            self._start_combo.addItem(str(rd.run_number))
            self._end_combo.addItem(str(rd.run_number))

        run_numbers = [rd.run_number for rd in self._runs]
        last = len(self._runs) - 1

        # restore start index
        if start_run_number in run_numbers:
            self._start_run_idx = run_numbers.index(start_run_number)
        else:
            self._start_run_idx = 0
        self._start_combo.setCurrentIndex(self._start_run_idx)

        # restore end index — extend to newest run if it was already at the end before refresh
        if end_run_number == old_last_run_number or end_run_number not in run_numbers:
            self._end_run_idx = last
        else:
            self._end_run_idx = run_numbers.index(end_run_number)
        self._end_combo.setCurrentIndex(self._end_run_idx)

        self._start_combo.blockSignals(False)
        self._end_combo.blockSignals(False)

        new_runs_added = run_numbers[-1] != old_last_run_number
        if new_runs_added:
            # jump to the newest run
            self._current_run_idx = last
        elif current_run_number in run_numbers:
            self._current_run_idx = run_numbers.index(current_run_number)
        else:
            self._current_run_idx = self._end_run_idx

        self._populate_run_combo()
        self._recompute_pairwise_diffs()
        self._status_lbl.setText(f"{len(self._runs)} runs loaded  [auto-refreshed]")
        self._status_lbl.setStyleSheet("color:#3fb950;")
        self._update_all_views()

    def _on_process_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select Data Folder")
        if not folder:
            return
        self._load_folder(folder)

    def _load_folder(self, folder: str):
        self._status_lbl.setText("Loading...")
        self._status_lbl.setStyleSheet("color:#d29922;")
        QApplication.processEvents()

        self._runs = load_all_runs(folder)
        self._current_folder = folder
        if not self._runs:
            self._status_lbl.setText("No data files found")
            self._status_lbl.setStyleSheet("color:#f85149;")
            return

        # populate start/end combos (all runs)
        self._start_combo.blockSignals(True)
        self._start_combo.clear()
        for rd in self._runs:
            self._start_combo.addItem(str(rd.run_number))
        self._start_combo.setCurrentIndex(0)
        self._start_combo.blockSignals(False)

        self._end_combo.blockSignals(True)
        self._end_combo.clear()
        for rd in self._runs:
            self._end_combo.addItem(str(rd.run_number))
        last = len(self._runs) - 1
        self._end_combo.setCurrentIndex(last)
        self._end_combo.blockSignals(False)

        self._start_run_idx = 0
        self._end_run_idx = last
        self._current_run_idx = last
        self._selected_module = None
        self._map.set_selected(None)
        self._deviation_stats_key = None

        self._populate_run_combo()
        self._recompute_pairwise_diffs()
        self._status_lbl.setText(f"{len(self._runs)} runs loaded")
        self._status_lbl.setStyleSheet("color:#3fb950;")
        self._refresh_btn.setEnabled(True)
        self._auto_refresh_btn.setEnabled(True)
        self._file_snapshot = self._take_file_snapshot(folder)
        if self._auto_refresh_btn.isChecked():
            self._auto_refresh_timer.start(self._auto_refresh_interval.value() * 1000)

        self._update_all_views()

    def _populate_run_combo(self):
        self._run_combo.blockSignals(True)
        self._run_combo.clear()
        for rd in self._active_runs:
            self._run_combo.addItem(str(rd.run_number))
        combo_idx = max(0, self._current_run_idx - self._start_run_idx)
        combo_idx = min(combo_idx, self._run_combo.count() - 1)
        self._run_combo.setCurrentIndex(combo_idx)
        self._run_combo.blockSignals(False)

    def _on_start_run_changed(self, index: int):
        if index < 0 or index >= len(self._runs):
            return
        self._start_run_idx = index
        if self._end_run_idx < self._start_run_idx:
            self._end_combo.blockSignals(True)
            self._end_combo.setCurrentIndex(index)
            self._end_combo.blockSignals(False)
            self._end_run_idx = index
        if self._current_run_idx < self._start_run_idx:
            self._current_run_idx = self._start_run_idx
        self._populate_run_combo()
        self._update_all_views()

    def _on_end_run_changed(self, index: int):
        if index < 0 or index >= len(self._runs):
            return
        self._end_run_idx = index
        if self._start_run_idx > self._end_run_idx:
            self._start_combo.blockSignals(True)
            self._start_combo.setCurrentIndex(index)
            self._start_combo.blockSignals(False)
            self._start_run_idx = index
        if self._current_run_idx > self._end_run_idx:
            self._current_run_idx = self._end_run_idx
        self._populate_run_combo()
        self._update_all_views()

    def _on_run_changed(self, index: int):
        if index < 0 or index >= len(self._active_runs):
            return
        if self._view_mode == 3:
            return
        self._current_run_idx = self._start_run_idx + index
        self._update_geo_view()
        if self._view_mode == 2:
            self._update_irregular_table()
        curr_run_number = self._runs[self._current_run_idx].run_number
        for chart in self._charts:
            chart.set_current_run(curr_run_number)

    def _on_ref_changed(self, index: int):
        self._current_ref_idx = index
        if self._view_mode == 3:
            self._recompute_pairwise_diffs()
            self._update_summary_views()
        else:
            self._update_geo_view()
            self._update_irregular_table()
        self._update_line_charts()

    def _on_prev_run(self):
        combo_idx = self._current_run_idx - self._start_run_idx
        if combo_idx > 0:
            self._run_combo.setCurrentIndex(combo_idx - 1)

    def _on_next_run(self):
        combo_idx = self._current_run_idx - self._start_run_idx
        if combo_idx < len(self._active_runs) - 1:
            self._run_combo.setCurrentIndex(combo_idx + 1)

    def _on_apply_range(self):
        try:
            vmin = float(self._range_min.text())
            vmax = float(self._range_max.text())
        except ValueError:
            return
        if vmin >= vmax:
            return
        self._auto_range = False
        self._auto_btn.setChecked(False)
        self._manual_vmin = vmin
        self._manual_vmax = vmax
        self._update_geo_view()

    def _on_log_toggled(self):
        self._log_scale = self._log_btn.isChecked()
        self._map.set_log_scale(self._log_scale)

    def _on_auto_range(self):
        self._auto_range = self._auto_btn.isChecked()
        self._update_geo_view()

    def _on_jump_to_run(self, run_number: int):
        for i, rd in enumerate(self._active_runs):
            if rd.run_number == run_number:
                self._run_combo.setCurrentIndex(i)
                return

    def _on_chart_run_clicked(self, run_number: int):
        if self._view_mode != 3:
            self._on_jump_to_run(run_number)

    def _on_cycle_palette(self):
        if self._view_mode == 2:
            return  # drift mode uses a fixed palette
        self._palette_idx = (self._palette_idx + 1) % len(PALETTES)
        self._map.set_palette(self._palette_idx)

    def _on_view_mode_changed(self, index: int):
        prev_mode = self._view_mode
        self._view_mode = index
        uses_threshold = (index in (2, 3))
        for w in (self._thresh_lbl, self._thresh_g_input, self._thresh_g_pct,
                  self._thresh_w_lbl, self._thresh_w_input, self._thresh_w_pct):
            w.setVisible(uses_threshold)
        self._run_combo.setEnabled(index != 3)
        if prev_mode == 2 and index != 2:
            self._map.set_palette_override(None)
            self._map.set_palette(self._palette_idx)
        if index == 0:
            self._map.set_legend_mode("gain")
        elif index == 1:
            self._map.set_legend_mode("deviation")
        elif index == 2:
            self._map.set_palette_override(DRIFT_PALETTE)
            self._map.set_legend_mode("drift")
        elif index == 3:
            self._map.set_legend_mode("summary")
        if index == 3:
            self._update_summary_views()
        else:
            self._update_geo_view()
            if not (prev_mode in (0, 1) and index in (0, 1)):
                self._update_irregular_table()

    def _on_drift_threshold_changed(self):
        for inp, attr, default in (
            (self._thresh_g_input, "_thresh_g", 0.10),
            (self._thresh_w_input, "_thresh_w", 0.05),
        ):
            inp.blockSignals(True)
            try:
                val = float(inp.text())
                if val <= 0:
                    raise ValueError
                setattr(self, attr, val / 100.0)
            except ValueError:
                inp.setText(f"{getattr(self, attr) * 100:.1f}")
            finally:
                inp.blockSignals(False)
        if self._view_mode == 3:
            self._update_summary_views()
        else:
            self._update_geo_view()
            self._update_irregular_table()

    def _on_module_hovered(self, name: str):
        mod = self._mod_by_name.get(name)
        if not mod:
            return
        active = self._active_runs
        if not active:
            self._info.setText(f"{name} ({mod.mod_type})")
            return
        rd = self._runs[self._current_run_idx]
        mrec = rd.modules.get(name)
        ref = LMS_NAMES[self._current_ref_idx]
        if mrec:
            gain = mrec.gain_factors[self._current_ref_idx]
            base = (f"{name} ({mod.mod_type})  gain[{ref}]: {gain:.5f}"
                    f"  lms_peak: {mrec.lms_peak:.2f}"
                    f"  lms_sigma: {mrec.lms_sigma:.2f}")
            if self._view_mode == 1:
                glist = [r.modules[name].gain_factors[self._current_ref_idx]
                         for r in active if name in r.modules]
                if len(glist) > 1:
                    mean = sum(glist) / len(glist)
                    std = math.sqrt(sum((g - mean) ** 2 for g in glist) / len(glist))
                    dev = (gain - mean) / std if std > 0 else 0.0
                    base += f"  dev: {dev:+.2f}σ"
            self._info.setText(base)
        else:
            self._info.setText(f"{name} ({mod.mod_type})  no data")

    def _on_module_clicked(self, name: str):
        self._selected_module = name if name else None
        self._map.set_selected(self._selected_module)
        if self._selected_module:
            self._irregular_table.select_module(self._selected_module)
        else:
            self._irregular_table._table.clearSelection()
        self._update_line_charts()

    # ---- update views ----

    def _recompute_pairwise_diffs(self):
        """Pre-compute all consecutive-run relative gain changes across ALL runs.
        Stored as flat list of (name, mod_type, rel, pair_idx) where pair_idx is
        the index of the current run in self._runs. Called once on load and on
        ref index change so threshold/range updates just filter this list."""
        ref_idx = self._current_ref_idx
        runs = self._runs
        mod_by_name = self._mod_by_name
        diffs = []
        for i in range(1, len(runs)):
            rd_curr = runs[i]
            rd_prev = runs[i - 1]
            prev_gains = {mname: mrec.gain_factors[ref_idx]
                          for mname, mrec in rd_prev.modules.items()}
            curr_run = rd_curr.run_number
            prev_run = rd_prev.run_number
            for mname, mrec in rd_curr.modules.items():
                g_prev = prev_gains.get(mname)
                if g_prev is None:
                    continue
                g_curr = mrec.gain_factors[ref_idx]
                denom = min(abs(g_curr), abs(g_prev))
                rel = math.inf if denom == 0 else abs(g_curr - g_prev) / denom
                mod = mod_by_name.get(mname)
                diffs.append((mname, mod.mod_type if mod else "?", rel, i,
                               curr_run, prev_run))
        self._pairwise_diffs = diffs
        self._pairwise_ref_idx = ref_idx

    def _update_all_views(self):
        if self._view_mode == 3:
            self._update_summary_views()
        else:
            self._update_geo_view()
            self._update_irregular_table()
        self._update_line_charts()

    def _update_summary_views(self):
        """Filter pre-computed pairwise diffs by active range + threshold."""
        if not self._runs:
            return
        start_idx = self._start_run_idx
        end_idx = self._end_run_idx
        thresh_g = self._thresh_g
        thresh_w = self._thresh_w
        stats: Dict[str, List] = {}
        for name, mod_type, rel, pair_idx, curr_run, prev_run in self._pairwise_diffs:
            if pair_idx <= start_idx or pair_idx > end_idx:
                continue
            threshold = thresh_g if name.startswith("G") else thresh_w
            if rel <= threshold:
                continue
            s = stats.get(name)
            if s is None:
                stats[name] = [1, rel, curr_run, prev_run, mod_type]
            else:
                s[0] += 1
                if rel > s[1]:
                    s[1] = rel
                    s[2] = curr_run
                    s[3] = prev_run

        entries = [SummaryEntry(name=n, mod_type=s[4], drift_count=s[0],
                                max_rel_change=s[1], max_run=s[2], max_prev_run=s[3])
                   for n, s in stats.items()]
        entries.sort(key=lambda e: (
        0 if e.name.startswith("W") else 1,
        0 if e.name.startswith("W") else int(math.isinf(e.max_rel_change)),
        -e.drift_count,
        0 if math.isinf(e.max_rel_change) else -e.max_rel_change,
    ))

        values = {e.name: float(e.drift_count) for e in entries}
        if self._auto_range:
            vmax = max(values.values()) if values else 1.0
            self._range_min.setText("0.0000")
            self._range_max.setText(f"{vmax:.4f}")
            self._map.set_gain_data(values, 0.0, vmax)
        else:
            self._map.set_gain_data(values, self._manual_vmin, self._manual_vmax)
        self._irregular_table.set_summary_data(entries)

    def _update_geo_view(self):
        active = self._active_runs
        if not active:
            return
        rd = self._runs[self._current_run_idx]
        ref_idx = self._current_ref_idx

        if self._view_mode == 1:
            # Deviation (σ): signed (gain − mean) / std across active runs
            cache_key = (self._start_run_idx, self._end_run_idx, ref_idx)
            if self._deviation_stats_key != cache_key:
                all_gains: Dict[str, List[float]] = {}
                for r in active:
                    for mname, mrec in r.modules.items():
                        lst = all_gains.get(mname)
                        if lst is None:
                            all_gains[mname] = [mrec.gain_factors[ref_idx]]
                        else:
                            lst.append(mrec.gain_factors[ref_idx])
                dev_stats: Dict[str, Tuple[float, float]] = {}
                for mname, glist in all_gains.items():
                    mean = sum(glist) / len(glist)
                    var = sum((g - mean) ** 2 for g in glist) / len(glist)
                    dev_stats[mname] = (mean, math.sqrt(var))
                self._deviation_stats_cache = dev_stats
                self._deviation_stats_key = cache_key
            stats = self._deviation_stats_cache

            values: Dict[str, float] = {}
            for mname, mrec in rd.modules.items():
                mean, std = stats.get(mname, (0.0, 0.0))
                if std > 0:
                    values[mname] = (mrec.gain_factors[ref_idx] - mean) / std
                else:
                    values[mname] = 0.0

            if self._auto_range:
                if values:
                    abs_max = max(abs(v) for v in values.values())
                    abs_max = abs_max if abs_max > 0 else 1.0
                    vmin, vmax = -abs_max, abs_max
                else:
                    vmin, vmax = -3.0, 3.0
                self._range_min.setText(f"{vmin:.4f}")
                self._range_max.setText(f"{vmax:.4f}")
            else:
                vmin = self._manual_vmin
                vmax = self._manual_vmax
        elif self._view_mode == 2:
            # Run-to-Run Drift: (gain_current - gain_prev) / gain_prev
            curr_pos = next((i for i, r in enumerate(active)
                             if r.run_number == rd.run_number), None)
            values: Dict[str, float] = {}
            if curr_pos is not None and curr_pos > 0:
                rd_prev = active[curr_pos - 1]
                for mname, mrec in rd.modules.items():
                    prev_mrec = rd_prev.modules.get(mname)
                    if prev_mrec and prev_mrec.gain_factors[ref_idx] != 0:
                        g_prev = prev_mrec.gain_factors[ref_idx]
                        values[mname] = (mrec.gain_factors[ref_idx] - g_prev) / g_prev

            threshold = max(self._thresh_g, self._thresh_w)
            if self._auto_range:
                vmin, vmax = -threshold, threshold
                self._range_min.setText(f"{vmin:.4f}")
                self._range_max.setText(f"{vmax:.4f}")
            else:
                vmin = self._manual_vmin
                vmax = self._manual_vmax

        else:
            # Gain Factor mode (original)
            values: Dict[str, float] = {}
            for mname, mrec in rd.modules.items():
                values[mname] = mrec.gain_factors[ref_idx]

            if self._auto_range:
                if values:
                    sorted_v = sorted(values.values())
                    n = len(sorted_v)
                    lo_idx = max(0, int(n * 0.02))
                    hi_idx = min(n - 1, int(n * 0.98))
                    vmin = sorted_v[lo_idx]
                    vmax = sorted_v[hi_idx]
                    if vmin == vmax:
                        vmin -= 0.01
                        vmax += 0.01
                else:
                    vmin, vmax = 0.9, 1.1
                self._range_min.setText(f"{vmin:.4f}")
                self._range_max.setText(f"{vmax:.4f}")
            else:
                vmin = self._manual_vmin
                vmax = self._manual_vmax

        self._map.set_gain_data(values, vmin, vmax)

    def _update_line_charts(self):
        active = self._active_runs
        if not active:
            return

        first_run = active[0].run_number

        ref_idx = self._current_ref_idx

        if self._selected_module:
            mname = self._selected_module
            # single pass: all 3 channels share the same valid-run check
            chart_data: List[Tuple[List, List, List]] = [([], [], []) for _ in range(3)]
            for idx, rd in enumerate(active):
                mrec = rd.modules.get(mname)
                if mrec is None:
                    continue
                run_num = rd.run_number
                gf = mrec.gain_factors
                for i in range(3):
                    d = chart_data[i]
                    d[0].append(idx + 1)
                    d[1].append(run_num)
                    d[2].append(gf[i])
            curr_run_number = self._runs[self._current_run_idx].run_number
            for i, lms_name in enumerate(LMS_NAMES):
                indices, actual_runs, gains = chart_data[i]
                self._charts[i].set_data(
                    indices, gains, [],
                    f"{mname}  gain[{lms_name}]  (run1={first_run})",
                    actual_runs)
                self._charts[i].set_highlighted(i == ref_idx)
                self._charts[i].set_current_run(curr_run_number)
        else:
            # single pass: collect LMS ratio data for all 3 channels at once
            chart_data2: List[Tuple[List, List, List, List]] = [([], [], [], []) for _ in range(3)]
            for idx, rd in enumerate(active):
                run_num = rd.run_number
                for i, lms_name in enumerate(LMS_NAMES):
                    rec = rd.lms.get(lms_name)
                    if rec is None or rec.alpha_peak == 0 or rec.lms_peak == 0:
                        continue
                    ratio = rec.lms_peak / rec.alpha_peak
                    rel_lms = rec.lms_sigma / rec.lms_peak
                    rel_alpha = rec.alpha_sigma / rec.alpha_peak
                    err = ratio * math.sqrt(rel_lms ** 2 + rel_alpha ** 2)
                    d = chart_data2[i]
                    d[0].append(idx + 1)
                    d[1].append(run_num)
                    d[2].append(ratio)
                    d[3].append(err)
            curr_run_number = self._runs[self._current_run_idx].run_number
            for i, lms_name in enumerate(LMS_NAMES):
                indices, actual_runs, ratios, errors = chart_data2[i]
                self._charts[i].set_data(
                    indices, ratios, errors,
                    f"{lms_name} (lms/alpha)  (run1={first_run})",
                    actual_runs)
                self._charts[i].set_highlighted(i == ref_idx)
                self._charts[i].set_current_run(curr_run_number)

    def _update_irregular_table(self):
        active = self._active_runs
        if not active:
            return
        if self._view_mode == 2:
            rd = self._runs[self._current_run_idx]
            curr_pos = next((i for i, r in enumerate(active)
                             if r.run_number == rd.run_number), None)
            if curr_pos is not None and curr_pos > 0:
                entries = compute_drift_entries(
                    rd, active[curr_pos - 1],
                    self._current_ref_idx, self._mod_by_name,
                    self._thresh_g, self._thresh_w)
            else:
                entries = []
            self._irregular_table.set_drift_data(entries)
        else:
            entries = compute_irregular_entries(
                active, self._current_ref_idx, self._mod_by_name)
            self._irregular_table.set_data(entries)


# ===========================================================================
#  Entry point
# ===========================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", dest="folder", default=None,
                        help="Data folder to load on startup")
    args, qt_args = parser.parse_known_args()

    app = QApplication([sys.argv[0]] + qt_args)
    win = GainMonitorWindow()
    if args.folder:
        win._process_btn.setEnabled(False)
        win._load_folder(args.folder)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
