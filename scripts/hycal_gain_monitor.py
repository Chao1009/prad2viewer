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
import shutil
import subprocess
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import uproot
    _UPROOT_OK = True
except ImportError:
    _UPROOT_OK = False

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QComboBox, QLineEdit, QSpinBox,
    QFileDialog, QSplitter, QSizePolicy, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QMenu,
    QDialog, QFormLayout, QTextEdit, QMessageBox,
)
from PyQt6.QtCore import Qt, QRectF, QPointF, pyqtSignal, QTimer, QProcess
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QFont, QPalette,
)

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget, PALETTES, PALETTE_NAMES,
    apply_theme_palette, set_theme, available_themes, THEME, themed,
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

# Default palette for non-drift modes (matches historical gain-monitor look).
_DEFAULT_PALETTE = "blue-orange"

# Separate palette used only in Run-to-Run Drift mode (not cycled by the user)
DRIFT_PALETTE = [
    (0.00, (0, 210, 230)),   # cyan  — large negative drift
    (0.50, (80, 80, 80)),    # grey  — no drift (always maps to 0 in drift mode)
    (1.00, (249, 115, 22)),  # orange — large positive drift
]


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


# ===========================================================================
#  HyCal Gain Map Widget
# ===========================================================================

_LEGEND_ITEMS = {
    "drift": [
        (QColor(0, 210, 230),  "gain decreases"),
        (QColor(80, 80, 80),   "stable"),
        (QColor(249, 115, 22), "gain increases"),
    ],
    "summary": [
        (QColor(10, 42, 110),  "low drift count"),
        (QColor(249, 115, 22), "high drift count"),
    ],
    "gain": [
        (QColor(10, 42, 110),  "low gain"),
        (QColor(249, 115, 22), "high gain"),
    ],
    "deviation": [
        (QColor(10, 42, 110),  "below mean"),
        (QColor(80, 80, 80),   "near mean"),
        (QColor(249, 115, 22), "above mean"),
    ],
}


class HyCalGainMapWidget(HyCalMapWidget):
    """Gain-monitor specialisation of the shared HyCal map widget.

    Adds a custom palette override (used by Run-to-Run Drift mode), a
    persistent module selection highlight, and a legend overlay above the
    colour bar that explains the active view mode.
    """

    CB_MAX_WIDTH = 300

    def __init__(self, parent=None):
        super().__init__(parent, shrink=0.90, margin_top=8,
                         enable_zoom_pan=True)
        self._palette_override = None
        self._legend_mode: Optional[str] = None
        self._selected: Optional[str] = None

    # -- public API additions --

    def set_gain_data(self, values: Dict[str, float],
                      vmin: float, vmax: float):
        self._values = values
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def set_palette(self, idx_or_name):
        self._palette_override = None
        super().set_palette(idx_or_name)

    def set_palette_override(self, stops):
        """Use a custom stops list instead of the indexed palette. Pass None to clear."""
        self._palette_override = stops
        self.update()

    def set_legend_mode(self, mode: Optional[str]):
        if mode != self._legend_mode:
            self._legend_mode = mode
            self.update()

    def set_selected(self, name: Optional[str]):
        self._selected = name
        self.update()

    # -- base hooks --

    def palette_stops(self):
        if self._palette_override is not None:
            return self._palette_override
        return super().palette_stops()

    def _fmt_value(self, v: float) -> str:
        return f"{v:.4f}"

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        if v is None:
            return name
        return f"{name}: {v:.5f}"

    def _colorbar_center_text(self) -> str:
        if self._palette_override is not None:
            return "cyan-grey-orange"
        return super()._colorbar_center_text()

    def _paint_empty(self, p, w, h):
        if not self._values:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 12))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter, "No data loaded")

    def _paint_overlays(self, p, w, h):
        if self._selected and self._selected in self._rects:
            p.setPen(QPen(QColor(THEME.SELECT_BORDER), 2.5))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._selected])
        super()._paint_overlays(p, w, h)

    def _paint_after_colorbar(self, p, w, h):
        items = _LEGEND_ITEMS.get(self._legend_mode)
        if not items or self._cb_rect is None:
            return
        cb_y = self._cb_rect.y()
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
            p.setPen(QColor(THEME.TEXT))
            p.drawText(QRectF(x + swatch + gap, ly + pad,
                              fm.horizontalAdvance(label), lh),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       label)
            x += item_w + spacing

    def _handle_click(self, pos):
        if self._cb_rect and self._cb_rect.contains(pos):
            self.paletteClicked.emit()
            return
        hit = self._hit(pos)
        if hit is not None:
            new_sel = None if hit == self._selected else hit
            self._selected = new_sel
            self.update()
            self.moduleClicked.emit(new_sel if new_sel else "")
        elif self._selected is not None:
            self._selected = None
            self.update()
            self.moduleClicked.emit("")


def _chart_y_range(values: List[float], errors: List[float]) -> Tuple[float, float]:
    """Return y-axis (lo, hi): mean±20%, expanded if any data point falls outside."""
    mean = sum(values) / len(values)
    y_lo = mean * 0.8
    y_hi = mean * 1.2
    for j, v in enumerate(values):
        err = errors[j] if j < len(errors) else 0.0
        y_lo = min(y_lo, v - err)
        y_hi = max(y_hi, v + err)
    return y_lo, y_hi


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
        self._y_range: Optional[Tuple[float, float]] = None
        self.setMinimumHeight(100)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

    def set_y_range(self, lo: float, hi: float):
        self._y_range = (lo, hi)
        self.update()

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
        p.fillRect(0, 0, w, h, QColor(THEME.CANVAS))

        # title
        title_color = QColor(THEME.HIGHLIGHT) if self._highlighted else QColor(THEME.ACCENT)
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
            p.setPen(QColor(THEME.TEXT_MUTED))
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

        if self._y_range is not None:
            y_lo, y_hi = self._y_range
        else:
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

        if y_hi == y_lo:
            y_lo -= 0.05
            y_hi += 0.05

        def to_sx(v):
            return px + (v - x_min) / (x_max - x_min) * pw

        def to_sy(v):
            return py + ph - (v - y_lo) / (y_hi - y_lo) * ph

        # grid + axes
        p.setPen(QPen(QColor(THEME.BUTTON), 1, Qt.PenStyle.DotLine))
        y_ticks = self._nice_ticks(y_lo, y_hi, 5)
        for yt in y_ticks:
            sy = to_sy(yt)
            p.drawLine(QPointF(px, sy), QPointF(px + pw, sy))

        # axes border
        p.setPen(QPen(QColor(THEME.BORDER), 1))
        p.drawLine(QPointF(px, py), QPointF(px, py + ph))
        p.drawLine(QPointF(px, py + ph), QPointF(px + pw, py + ph))

        # y-axis labels
        p.setPen(QColor(THEME.TEXT_DIM))
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
        p.setPen(QPen(QColor(THEME.TEXT_DIM), 1))
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

        series_color = QColor(THEME.HIGHLIGHT) if self._highlighted else QColor(THEME.ACCENT)

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
                    p.setPen(QPen(QColor(THEME.DANGER), 2))
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
            p.setPen(QColor(THEME.TEXT_STRONG))
            for j, ln in enumerate(lines):
                p.drawText(QRectF(tx, ty + j * fm.height(), tw, fm.height()),
                           Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                           ln)

        p.end()


# ===========================================================================
#  ROOT histogram display widget
# ===========================================================================

_HIST_CACHE_MAX = 20


class RootHistWidget(QWidget):
    """Displays a TH1 histogram read from a fitted LMS ROOT file.

    Left-drag to zoom into an x range; right-click → Unzoom resets to 0–2000.
    """

    PAD_L, PAD_R, PAD_T, PAD_B = 55, 16, 24, 40

    _X_DEFAULT_LO = 0.0
    _X_DEFAULT_HI = 1000.0

    def __init__(self, parent=None):
        super().__init__(parent)
        self._values: List[float] = []
        self._edges: List[float] = []
        self._title: str = ""
        self._gauss: Optional[Tuple[float, float, float, float, float]] = None  # amp,mean,sigma,xmin,xmax
        self._x_lo = self._X_DEFAULT_LO
        self._x_hi = self._X_DEFAULT_HI
        self._drag_start: Optional[float] = None  # data-x where drag began
        self._drag_cur:   Optional[float] = None  # data-x of current cursor
        self.setMinimumHeight(80)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

    # ------------------------------------------------------------------
    def set_histogram(self, values, edges, title: str = "",
                      gauss: Optional[Tuple[float, float, float, float, float]] = None):
        self._values = list(values)
        self._edges = list(edges)
        self._title = title
        self._gauss = gauss
        self._x_lo = self._X_DEFAULT_LO
        self._x_hi = self._X_DEFAULT_HI
        self._drag_start = self._drag_cur = None
        self.update()

    def clear(self):
        self._values = []
        self._edges = []
        self._title = ""
        self._gauss = None
        self._x_lo = self._X_DEFAULT_LO
        self._x_hi = self._X_DEFAULT_HI
        self._drag_start = self._drag_cur = None
        self.update()

    # ------------------------------------------------------------------
    def _plot_rect(self):
        w, h = self.width(), self.height()
        return self.PAD_L, self.PAD_T, w - self.PAD_L - self.PAD_R, h - self.PAD_T - self.PAD_B

    def _sx_to_data(self, sx: float) -> float:
        px, _py, pw, _ph = self._plot_rect()
        if pw <= 0 or self._x_hi == self._x_lo:
            return self._x_lo
        return self._x_lo + (sx - px) / pw * (self._x_hi - self._x_lo)

    # ------------------------------------------------------------------
    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            px, py, pw, ph = self._plot_rect()
            mx, my = event.position().x(), event.position().y()
            if px <= mx <= px + pw and py <= my <= py + ph + self.PAD_B:
                self._drag_start = self._sx_to_data(mx)
                self._drag_cur   = self._drag_start
                self.update()

    def mouseMoveEvent(self, event):
        if self._drag_start is not None:
            self._drag_cur = self._sx_to_data(event.position().x())
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton and self._drag_start is not None:
            d_end = self._sx_to_data(event.position().x())
            d_start = self._drag_start
            self._drag_start = self._drag_cur = None
            span = self._x_hi - self._x_lo
            if abs(d_end - d_start) > span * 0.01:
                self._x_lo = min(d_start, d_end)
                self._x_hi = max(d_start, d_end)
            self.update()

    def contextMenuEvent(self, event):
        menu = QMenu(self)
        menu.setStyleSheet(themed(
            "QMenu{background:#161b22;color:#c9d1d9;border:1px solid #30363d;}"
            "QMenu::item:selected{background:#1f6feb;}"))
        menu.addAction("Unzoom").triggered.connect(self._unzoom)
        menu.exec(event.globalPos())

    def _unzoom(self):
        if self._edges:
            self._x_lo = self._edges[0]
            self._x_hi = self._edges[-1]
        else:
            self._x_lo = self._X_DEFAULT_LO
            self._x_hi = self._X_DEFAULT_HI
        self.update()

    # ------------------------------------------------------------------
    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(THEME.CANVAS))

        if self._title:
            p.setPen(QColor(THEME.ACCENT))
            p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
            p.drawText(QRectF(self.PAD_L, 2, w - self.PAD_L - self.PAD_R, 20),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                       self._title)

        if not self._values or not self._edges or len(self._edges) < 2:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 10))
            p.drawText(QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter, "No histogram")
            p.end()
            return

        px, py, pw, ph = self._plot_rect()
        if pw < 20 or ph < 20:
            p.end()
            return

        x_lo = self._x_lo
        x_hi = self._x_hi if self._x_hi > self._x_lo else self._x_lo + 1

        # y scale from visible bins only
        vis_max = max(
            (v for i, v in enumerate(self._values)
             if i + 1 < len(self._edges)
             and self._edges[i + 1] > x_lo
             and self._edges[i] < x_hi),
            default=0.0)
        y_hi = vis_max * 1.1 if vis_max > 0 else 1.0

        def to_sx(v):
            return px + (v - x_lo) / (x_hi - x_lo) * pw

        def to_sy(v):
            return py + ph * (1.0 - v / y_hi)

        # dotted grid
        p.setPen(QPen(QColor(THEME.BUTTON), 1, Qt.PenStyle.DotLine))
        n_yticks = 5
        for i in range(n_yticks + 1):
            sy = py + ph * i / n_yticks
            p.drawLine(QPointF(px, sy), QPointF(px + pw, sy))

        # bars
        p.setPen(Qt.PenStyle.NoPen)
        for i, v in enumerate(self._values):
            if i + 1 >= len(self._edges):
                break
            b_lo, b_hi = self._edges[i], self._edges[i + 1]
            if b_hi <= x_lo or b_lo >= x_hi:
                continue
            sx1 = max(to_sx(b_lo), px)
            sx2 = min(to_sx(b_hi), px + pw)
            bar_top = to_sy(v)
            bar_h = (py + ph) - bar_top
            if bar_h > 0 and sx2 > sx1:
                p.fillRect(QRectF(sx1, bar_top, sx2 - sx1 - 1, bar_h), QColor(THEME.HIGHLIGHT))

        # Gaussian fit curve
        if self._gauss is not None:
            amp, mean, sigma, g_xmin, g_xmax = self._gauss
            if sigma > 0:
                draw_lo = max(g_xmin, x_lo)
                draw_hi = min(g_xmax, x_hi)
                if draw_hi > draw_lo:
                    n_pts = max(int((draw_hi - draw_lo) / (x_hi - x_lo) * pw), 2)
                    pts = []
                    for k in range(n_pts + 1):
                        gx = draw_lo + k / n_pts * (draw_hi - draw_lo)
                        gy = amp * math.exp(-0.5 * ((gx - mean) / sigma) ** 2)
                        sx = to_sx(gx)
                        sy = to_sy(gy)
                        pts.append(QPointF(sx, sy))
                    p.setPen(QPen(QColor("#00bcd4"), 2))
                    for k in range(len(pts) - 1):
                        p.drawLine(pts[k], pts[k + 1])

        # drag selection overlay
        if self._drag_start is not None and self._drag_cur is not None:
            d_lo = min(self._drag_start, self._drag_cur)
            d_hi = max(self._drag_start, self._drag_cur)
            sx1 = max(to_sx(d_lo), px)
            sx2 = min(to_sx(d_hi), px + pw)
            if sx2 > sx1:
                p.fillRect(QRectF(sx1, py, sx2 - sx1, ph), QColor(255, 255, 100, 50))
                p.setPen(QPen(QColor(255, 255, 100, 180), 1))
                p.drawRect(QRectF(sx1, py, sx2 - sx1, ph))

        # axes
        p.setPen(QPen(QColor(THEME.BORDER), 1))
        p.drawLine(QPointF(px, py), QPointF(px, py + ph))
        p.drawLine(QPointF(px, py + ph), QPointF(px + pw, py + ph))

        # y labels
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Consolas", 8))
        for i in range(n_yticks + 1):
            val = y_hi * (n_yticks - i) / n_yticks
            sy = py + ph * i / n_yticks
            p.drawText(QRectF(0, sy - 8, self.PAD_L - 4, 16),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       f"{val:.0f}")

        # x labels
        x_ticks = LMSLineChartWidget._nice_ticks(x_lo, x_hi, max(pw // 60, 2))
        for xt in x_ticks:
            sx = to_sx(xt)
            p.drawText(QRectF(sx - 25, py + ph + 2, 50, 16),
                       Qt.AlignmentFlag.AlignCenter, f"{xt:.0f}")

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
        lbl.setStyleSheet(themed("color:#c9d1d9;"))
        fbar.addWidget(lbl)

        self._search = QLineEdit()
        self._search.setPlaceholderText("module name...")
        self._search.setFixedWidth(120)
        self._search.setFont(QFont("Consolas", 10))
        self._search.setStyleSheet(themed(
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 4px;}"))
        self._search.textChanged.connect(self._apply_filter)
        fbar.addWidget(self._search)

        lbl2 = QLabel("Type:")
        lbl2.setFont(QFont("Consolas", 10))
        lbl2.setStyleSheet(themed("color:#c9d1d9;"))
        fbar.addWidget(lbl2)

        self._type_filter = QComboBox()
        self._type_filter.addItems(["All", "PbWO4", "PbGlass"])
        self._type_filter.setFixedWidth(100)
        self._type_filter.setFont(QFont("Consolas", 10))
        self._type_filter.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}"))
        self._type_filter.currentIndexChanged.connect(
            lambda _: self._apply_filter())
        fbar.addWidget(self._type_filter)

        fbar.addStretch()

        count_lbl = QLabel("")
        count_lbl.setFont(QFont("Consolas", 10))
        count_lbl.setStyleSheet(themed("color:#8b949e;"))
        self._count_lbl = count_lbl
        fbar.addWidget(count_lbl)

        layout.addLayout(fbar)

        # table
        self._table = QTableWidget()
        self._table.setColumnCount(6)
        self._table.setHorizontalHeaderLabels(
            ["Module", "Run", "Gain", "Mean", "Std Dev", "Dev (sigma)"])
        self._table.setFont(QFont("Consolas", 10))
        self._table.setStyleSheet(themed(
            "QTableWidget{background:#0d1117;color:#c9d1d9;"
            "gridline-color:#21262d;border:1px solid #30363d;}"
            "QTableWidget::item{padding:2px 6px;}"
            "QHeaderView::section{background:#161b22;color:#58a6ff;"
            "border:1px solid #30363d;font:bold 10pt Consolas;padding:4px;}"))
        self._table.setAlternatingRowColors(True)
        pal = self._table.palette()
        pal.setColor(QPalette.ColorRole.AlternateBase, QColor(THEME.ALT_BASE))
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
            item.setForeground(QColor(THEME.WARN))
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
            color = QColor(THEME.DANGER) if e.rel_change < 0 else QColor(THEME.SUCCESS)
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
            item.setForeground(QColor(THEME.DANGER))
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
#  Analyze Data dialog
# ===========================================================================

_SCRIPT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "run_gain_monitor.sh")


class AnalyzeDialog(QDialog):
    """Popup that runs run_gain_monitor.sh and streams its output."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Analyze Data")
        self.resize(700, 500)
        self.setStyleSheet(themed(
            "QDialog{background:#0d1117;color:#c9d1d9;}"
            "QLabel{color:#c9d1d9;font-family:Consolas;font-size:10pt;}"
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;"
            "font-family:Consolas;font-size:10pt;}"
            "QTextEdit{background:#0a0e14;color:#c9d1d9;"
            "border:1px solid #30363d;font-family:Consolas;font-size:9pt;}"
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 12px;"
            "font:bold 10pt Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#30363d;}"
            "QPushButton:disabled{color:#555;}"))

        self._process = QProcess(self)
        self._process.readyReadStandardOutput.connect(self._on_stdout)
        self._process.readyReadStandardError.connect(self._on_stderr)
        self._process.finished.connect(self._on_finished)

        root = QVBoxLayout(self)
        root.setSpacing(8)

        # inputs
        form = QFormLayout()
        form.setSpacing(6)

        self._run_edit = QLineEdit()
        self._run_edit.setPlaceholderText("e.g. 023735")
        self._cpu_edit = QLineEdit("25")

        def _dir_row(placeholder, browse_slot):
            row = QHBoxLayout()
            edit = QLineEdit()
            edit.setPlaceholderText(placeholder)
            btn = QPushButton("Browse…")
            btn.setFixedWidth(80)
            btn.clicked.connect(browse_slot)
            row.addWidget(edit)
            row.addWidget(btn)
            return row, edit

        in_row, self._indir_edit = _dir_row("/data/evio/data", self._browse_indir)
        out_row, self._outdir_edit = _dir_row(
            "/home/clasrun/prad2_daq/gain_monitoring/gain_monitor_output",
            self._browse_outdir)

        form.addRow("Run number:", self._run_edit)
        form.addRow("Number of CPUs:", self._cpu_edit)
        form.addRow("Input directory:", in_row)
        form.addRow("Output directory:", out_row)
        root.addLayout(form)

        # buttons
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Run")
        self._run_btn.setStyleSheet(themed(
            "QPushButton{background:#1f6feb;color:white;border:1px solid #388bfd;"
            "padding:4px 16px;font:bold 10pt Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#388bfd;}"
            "QPushButton:disabled{background:#21262d;color:#555;border-color:#30363d;}"))
        self._run_btn.clicked.connect(self._on_run)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(self._run_btn)
        btn_row.addWidget(self._stop_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

        # output console
        self._console = QTextEdit()
        self._console.setReadOnly(True)
        self._console.document().setMaximumBlockCount(5000)
        root.addWidget(self._console, stretch=1)

    def _browse_indir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if d:
            self._indir_edit.setText(d)

    def _browse_outdir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if d:
            self._outdir_edit.setText(d)

    def _on_run(self):
        run_num = self._run_edit.text().strip()
        n_cpu = self._cpu_edit.text().strip()
        if not run_num:
            self._append("<span style='color:#f85149'>Please enter a run number.</span>")
            return
        if not n_cpu.isdigit() or int(n_cpu) < 1:
            self._append("<span style='color:#f85149'>Number of CPUs must be a positive integer.</span>")
            return
        if not os.path.exists(_SCRIPT_PATH):
            self._append(f"<span style='color:#f85149'>Script not found: {_SCRIPT_PATH}</span>")
            return

        # build environment with optional directory overrides
        env = self._process.processEnvironment()
        from PyQt6.QtCore import QProcessEnvironment
        env = QProcessEnvironment.systemEnvironment()
        indir = self._indir_edit.text().strip()
        outdir = self._outdir_edit.text().strip()
        if indir:
            env.insert("INPUTDIR", indir)
        if outdir:
            env.insert("OUTPUTDIR", outdir)
        self._process.setProcessEnvironment(env)
        self._process.setWorkingDirectory(os.path.dirname(_SCRIPT_PATH))

        self._console.clear()
        extra = ""
        if indir:
            extra += f" INPUTDIR={indir}"
        if outdir:
            extra += f" OUTPUTDIR={outdir}"
        self._append(f"<span style='color:#8b949e'>${extra} {_SCRIPT_PATH} {run_num} {n_cpu}</span><br>")
        self._run_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._process.start("bash", [_SCRIPT_PATH, run_num, n_cpu])

    def _on_stop(self):
        self._process.kill()

    def _on_stdout(self):
        data = self._process.readAllStandardOutput().data().decode(errors="replace")
        self._append(data.replace("\n", "<br>"))

    def _on_stderr(self):
        data = self._process.readAllStandardError().data().decode(errors="replace")
        self._append(f"<span style='color:#f85149'>{data.replace(chr(10), '<br>')}</span>")

    def _on_finished(self, exit_code, exit_status):
        self._run_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)
        color = "#3fb950" if exit_code == 0 else "#f85149"
        self._append(f"<span style='color:{color}'>[Process finished with exit code {exit_code}]</span>")

    def _append(self, html: str):
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)
        self._console.insertHtml(themed(html))
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)

    def closeEvent(self, event):
        if self._process.state() != QProcess.ProcessState.NotRunning:
            self._process.kill()
            self._process.waitForFinished(2000)
        super().closeEvent(event)


# ===========================================================================
#  Get Data dialog
# ===========================================================================

_LOCAL_DATA_BASE  = "/data/evio/data"
_REMOTE_HOST      = "clondaq6"
_REMOTE_DATA_BASE = "/data/stage6"


def _fmt_bytes(b: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if b < 1024:
            return f"{b:.1f} {unit}"
        b //= 1024
    return f"{b:.1f} PB"


def _check_disk_space(remote_host: str, remote_run_dir: str,
                      local_base: str, f_start: int, f_end: int):
    """Return (needed_bytes, free_bytes) for evio files [f_start, f_end].

    SSHes to remote_host and sums the sizes of files whose .evio.NNN suffix
    falls within [f_start, f_end].  Checks free space on the filesystem that
    contains local_base (or its nearest existing ancestor).
    Raises RuntimeError if the SSH call fails.
    """
    result = subprocess.run(
        ["ssh", "-o", "ConnectTimeout=10",
         remote_host, f"ls -l {remote_run_dir}/ 2>/dev/null"],
        capture_output=True, text=True, timeout=30,
    )
    # exit 255 means SSH itself failed to connect; other non-zero codes (e.g. 2
    # when the remote directory doesn't exist yet) are fine — we just get no output
    # and needed stays 0, so the space check passes trivially.
    if result.returncode == 255:
        raise RuntimeError(result.stderr.strip() or "SSH connection failed")

    needed = 0
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) < 9:
            continue
        m = re.search(r'\.evio\.(\d+)$', parts[-1])
        if not m:
            continue
        n = int(m.group(1))
        if f_start <= n <= f_end:
            try:
                needed += int(parts[4])
            except (ValueError, IndexError):
                pass

    # walk up to the nearest existing directory so disk_usage doesn't fail
    check_path = local_base
    while check_path and not os.path.exists(check_path):
        check_path = os.path.dirname(check_path)
    free = shutil.disk_usage(check_path or "/").free
    return needed, free


class GetDataDialog(QDialog):
    """Popup that scps evio files from the DAQ machine for a given run."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Get Data")
        self.resize(700, 520)
        self.setStyleSheet(themed(
            "QDialog{background:#0d1117;color:#c9d1d9;}"
            "QLabel{color:#c9d1d9;font-family:Consolas;font-size:10pt;}"
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;"
            "font-family:Consolas;font-size:10pt;}"
            "QTextEdit{background:#0a0e14;color:#c9d1d9;"
            "border:1px solid #30363d;font-family:Consolas;font-size:9pt;}"
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 12px;"
            "font:bold 10pt Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#30363d;}"
            "QPushButton:disabled{color:#555;}"))

        self._process = QProcess(self)
        self._process.readyReadStandardOutput.connect(self._on_stdout)
        self._process.readyReadStandardError.connect(self._on_stderr)
        self._process.finished.connect(self._on_finished)

        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ---- input fields ----
        form = QFormLayout()
        form.setSpacing(6)

        self._run_edit = QLineEdit()
        self._run_edit.setPlaceholderText("e.g. 023739")
        form.addRow("Run number:", self._run_edit)

        file_range_row = QHBoxLayout()
        self._start_edit = QLineEdit("0")
        self._start_edit.setFixedWidth(80)
        self._end_edit = QLineEdit("99")
        self._end_edit.setFixedWidth(80)
        file_range_row.addWidget(self._start_edit)
        file_range_row.addWidget(QLabel(" to "))
        file_range_row.addWidget(self._end_edit)
        file_range_row.addStretch()
        form.addRow("File number range:", file_range_row)

        def _dir_row(default, slot):
            row = QHBoxLayout()
            edit = QLineEdit(default)
            btn = QPushButton("Browse…")
            btn.setFixedWidth(80)
            btn.clicked.connect(slot)
            row.addWidget(edit)
            row.addWidget(btn)
            return row, edit

        in_row,  self._localbase_edit  = _dir_row(_LOCAL_DATA_BASE,  self._browse_local)
        form.addRow("Local data directory:", in_row)

        host_edit_row = QHBoxLayout()
        self._host_edit = QLineEdit(_REMOTE_HOST)
        self._rembase_edit = QLineEdit(_REMOTE_DATA_BASE)
        host_edit_row.addWidget(self._host_edit)
        host_edit_row.addWidget(QLabel("  base:"))
        host_edit_row.addWidget(self._rembase_edit)
        form.addRow("Remote host:", host_edit_row)

        root.addLayout(form)

        # ---- buttons ----
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Get Data")
        self._run_btn.setStyleSheet(themed(
            "QPushButton{background:#1f6feb;color:white;border:1px solid #388bfd;"
            "padding:4px 16px;font:bold 10pt Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#388bfd;}"
            "QPushButton:disabled{background:#21262d;color:#555;border-color:#30363d;}"))
        self._run_btn.clicked.connect(self._on_get)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(self._run_btn)
        btn_row.addWidget(self._stop_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

        # ---- terminal output ----
        self._console = QTextEdit()
        self._console.setReadOnly(True)
        self._console.document().setMaximumBlockCount(5000)
        root.addWidget(self._console, stretch=1)

    # ------------------------------------------------------------------
    def _browse_local(self):
        d = QFileDialog.getExistingDirectory(self, "Select Local Data Directory")
        if d:
            self._localbase_edit.setText(d)

    def _on_get(self):
        run_num = self._run_edit.text().strip()
        if not run_num:
            self._append("<span style='color:#f85149'>Please enter a run number.</span>")
            return

        start_text = self._start_edit.text().strip() or "0"
        end_text   = self._end_edit.text().strip()   or "9999"
        if not start_text.isdigit() or not end_text.isdigit():
            self._append("<span style='color:#f85149'>File number range must be integers.</span>")
            return
        f_start = int(start_text)
        f_end   = int(end_text)
        if f_end < f_start:
            self._append("<span style='color:#f85149'>End file number must be ≥ start.</span>")
            return

        local_base  = self._localbase_edit.text().strip() or _LOCAL_DATA_BASE
        remote_host = self._host_edit.text().strip() or _REMOTE_HOST
        remote_base = self._rembase_edit.text().strip() or _REMOTE_DATA_BASE

        local_run_dir  = f"{local_base}/prad_{run_num}"
        remote_run_dir = f"{remote_base}/prad_{run_num}"

        # -- pre-check: find files in range that already exist locally --
        existing = []
        if os.path.isdir(local_run_dir):
            import glob as _glob
            for path in sorted(_glob.glob(
                    f"{local_run_dir}/prad_{run_num}.evio.*")):
                m = re.search(r'\.evio\.(\d+)$', os.path.basename(path))
                if m:
                    n = int(m.group(1))
                    if f_start <= n <= f_end:
                        existing.append(os.path.basename(path))

        if existing:
            box = QMessageBox(self)
            box.setWindowTitle("Files Already Present")
            box.setText(
                f"{len(existing)} file(s) in the requested range already exist "
                f"in {local_run_dir}.\nThey will be skipped; only missing files "
                f"will be copied.")
            box.setDetailedText("\n".join(existing))
            box.setStandardButtons(
                QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)
            box.setDefaultButton(QMessageBox.StandardButton.Ok)
            if box.exec() == QMessageBox.StandardButton.Cancel:
                return

        # -- disk space check --
        try:
            needed, free = _check_disk_space(
                remote_host, remote_run_dir, local_base, f_start, f_end)
        except Exception as exc:
            box = QMessageBox(self)
            box.setWindowTitle("Disk Space Check Failed")
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText(f"Could not check remote file sizes:\n{exc}\n\nProceed anyway?")
            box.setStandardButtons(
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            box.setDefaultButton(QMessageBox.StandardButton.No)
            if box.exec() != QMessageBox.StandardButton.Yes:
                return
        else:
            if needed > free:
                box = QMessageBox(self)
                box.setWindowTitle("Insufficient Disk Space")
                box.setIcon(QMessageBox.Icon.Critical)
                box.setText(
                    "Not enough disk space to copy the requested files.\n\n"
                    f"  Required : {_fmt_bytes(needed)}\n"
                    f"  Available: {_fmt_bytes(free)}\n\n"
                    "Free up space and try again.")
                box.setStandardButtons(QMessageBox.StandardButton.Ok)
                box.exec()
                return

        # bash: list remote files in range, skip those already local, scp the rest
        bash_cmd = (
            f"mkdir -p {local_run_dir}\n"
            f"echo 'Local directory: {local_run_dir}'\n"
            f"echo 'Listing remote files...'\n"
            f"ALL_FILES=$(ssh {remote_host} 'ls {remote_run_dir}/' 2>/dev/null | sort)\n"
            f"COPIED=0\n"
            f"ALREADY=0\n"
            f"while IFS= read -r f; do\n"
            f"    NUM=$(echo \"$f\" | grep -oP '\\.evio\\.\\K[0-9]+')\n"
            f"    [ -z \"$NUM\" ] && continue\n"
            f"    N=$((10#$NUM))\n"
            f"    if [ \"$N\" -lt {f_start} ] || [ \"$N\" -gt {f_end} ]; then continue; fi\n"
            f"    if [ -f \"{local_run_dir}/$f\" ]; then\n"
            f"        echo \"  Already exists: $f (skipping)\"\n"
            f"        ALREADY=$((ALREADY+1))\n"
            f"    else\n"
            f"        echo \"  Copying $f\"\n"
            f"        scp {remote_host}:{remote_run_dir}/$f {local_run_dir}/\n"
            f"        COPIED=$((COPIED+1))\n"
            f"    fi\n"
            f"done <<< \"$ALL_FILES\"\n"
            f"echo \"Done. Copied $COPIED file(s), $ALREADY already present.\"\n"
        )

        self._console.clear()
        self._append(
            f"<span style='color:#8b949e'>Run {run_num}, files {f_start}–{f_end}"
            f" → {local_run_dir}</span><br>")
        if existing:
            self._append(
                f"<span style='color:#d29922'>{len(existing)} file(s) skipped "
                f"(already present).</span><br>")
        self._run_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._process.start("bash", ["-c", bash_cmd])

    def _on_stop(self):
        self._process.kill()

    def _on_stdout(self):
        data = self._process.readAllStandardOutput().data().decode(errors="replace")
        self._append(data.replace("\n", "<br>"))

    def _on_stderr(self):
        data = self._process.readAllStandardError().data().decode(errors="replace")
        self._append(f"<span style='color:#f85149'>{data.replace(chr(10), '<br>')}</span>")

    def _on_finished(self, exit_code, exit_status):
        self._run_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)
        color = "#3fb950" if exit_code == 0 else "#f85149"
        self._append(f"<span style='color:{color}'>[Process finished with exit code {exit_code}]</span>")

    def _append(self, html: str):
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)
        self._console.insertHtml(themed(html))
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)

    def closeEvent(self, event):
        if self._process.state() != QProcess.ProcessState.NotRunning:
            self._process.kill()
            self._process.waitForFinished(2000)
        super().closeEvent(event)


# ===========================================================================
#  Do It All dialog  (Get Data + Analyze Data combined)
# ===========================================================================

class DoItAllDialog(QDialog):
    """Runs scp then gain monitor analysis in a single sequential workflow."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Do It All")
        self.resize(750, 600)
        self.setStyleSheet(themed(
            "QDialog{background:#0d1117;color:#c9d1d9;}"
            "QLabel{color:#c9d1d9;font-family:Consolas;font-size:10pt;}"
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;"
            "font-family:Consolas;font-size:10pt;}"
            "QTextEdit{background:#0a0e14;color:#c9d1d9;"
            "border:1px solid #30363d;font-family:Consolas;font-size:9pt;}"
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 12px;"
            "font:bold 10pt Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#30363d;}"
            "QPushButton:disabled{color:#555;}"))

        self._process = QProcess(self)
        self._process.readyReadStandardOutput.connect(self._on_stdout)
        self._process.readyReadStandardError.connect(self._on_stderr)
        self._process.finished.connect(self._on_finished)

        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ---- input fields ----
        form = QFormLayout()
        form.setSpacing(6)

        self._run_edit = QLineEdit()
        self._run_edit.setPlaceholderText("e.g. 023739")
        form.addRow("Run number:", self._run_edit)

        # file range
        file_range_row = QHBoxLayout()
        self._start_edit = QLineEdit("0")
        self._start_edit.setFixedWidth(80)
        self._end_edit = QLineEdit("99")
        self._end_edit.setFixedWidth(80)
        file_range_row.addWidget(self._start_edit)
        file_range_row.addWidget(QLabel(" to "))
        file_range_row.addWidget(self._end_edit)
        file_range_row.addStretch()
        form.addRow("File number range:", file_range_row)

        self._cpu_edit = QLineEdit("25")
        form.addRow("Number of CPUs:", self._cpu_edit)

        def _dir_row(default, slot):
            row = QHBoxLayout()
            edit = QLineEdit(default)
            btn = QPushButton("Browse…")
            btn.setFixedWidth(80)
            btn.clicked.connect(slot)
            row.addWidget(edit)
            row.addWidget(btn)
            return row, edit

        in_row,  self._localbase_edit = _dir_row(_LOCAL_DATA_BASE,  self._browse_local)
        out_row, self._outdir_edit    = _dir_row(
            "/home/clasrun/prad2_daq/gain_monitoring/gain_monitor_output",
            self._browse_out)
        form.addRow("Local data directory:", in_row)
        form.addRow("Output directory:", out_row)

        host_row = QHBoxLayout()
        self._host_edit    = QLineEdit(_REMOTE_HOST)
        self._rembase_edit = QLineEdit(_REMOTE_DATA_BASE)
        host_row.addWidget(self._host_edit)
        host_row.addWidget(QLabel("  base:"))
        host_row.addWidget(self._rembase_edit)
        form.addRow("Remote host:", host_row)

        root.addLayout(form)

        # ---- buttons ----
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Run")
        self._run_btn.setStyleSheet(themed(
            "QPushButton{background:#1f6feb;color:white;border:1px solid #388bfd;"
            "padding:4px 16px;font:bold 10pt Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#388bfd;}"
            "QPushButton:disabled{background:#21262d;color:#555;border-color:#30363d;}"))
        self._run_btn.clicked.connect(self._on_run)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(self._run_btn)
        btn_row.addWidget(self._stop_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

        self._console = QTextEdit()
        self._console.setReadOnly(True)
        self._console.document().setMaximumBlockCount(5000)
        root.addWidget(self._console, stretch=1)

    # ------------------------------------------------------------------
    def _browse_local(self):
        d = QFileDialog.getExistingDirectory(self, "Select Local Data Directory")
        if d:
            self._localbase_edit.setText(d)

    def _browse_out(self):
        d = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if d:
            self._outdir_edit.setText(d)

    def _on_run(self):
        run_num    = self._run_edit.text().strip()
        n_cpu      = self._cpu_edit.text().strip()
        start_text = self._start_edit.text().strip() or "0"
        end_text   = self._end_edit.text().strip()   or "9999"

        if not run_num:
            self._append("<span style='color:#f85149'>Please enter a run number.</span>")
            return
        if not n_cpu.isdigit() or int(n_cpu) < 1:
            self._append("<span style='color:#f85149'>Number of CPUs must be a positive integer.</span>")
            return
        if not start_text.isdigit() or not end_text.isdigit():
            self._append("<span style='color:#f85149'>File number range must be integers.</span>")
            return
        f_start = int(start_text)
        f_end   = int(end_text)
        if f_end < f_start:
            self._append("<span style='color:#f85149'>End file number must be ≥ start.</span>")
            return
        if not os.path.exists(_SCRIPT_PATH):
            self._append(f"<span style='color:#f85149'>Script not found: {_SCRIPT_PATH}</span>")
            return

        local_base  = self._localbase_edit.text().strip() or _LOCAL_DATA_BASE
        remote_host = self._host_edit.text().strip()      or _REMOTE_HOST
        remote_base = self._rembase_edit.text().strip()   or _REMOTE_DATA_BASE
        outdir      = self._outdir_edit.text().strip()    or \
            "/home/clasrun/prad2_daq/gain_monitoring/gain_monitor_output"

        local_run_dir  = f"{local_base}/prad_{run_num}"
        remote_run_dir = f"{remote_base}/prad_{run_num}"
        script_dir     = os.path.dirname(_SCRIPT_PATH)

        # pre-check for existing local files
        existing = []
        if os.path.isdir(local_run_dir):
            import glob as _glob
            for path in sorted(_glob.glob(f"{local_run_dir}/prad_{run_num}.evio.*")):
                m = re.search(r'\.evio\.(\d+)$', os.path.basename(path))
                if m:
                    n = int(m.group(1))
                    if f_start <= n <= f_end:
                        existing.append(os.path.basename(path))

        if existing:
            box = QMessageBox(self)
            box.setWindowTitle("Files Already Present")
            box.setText(
                f"{len(existing)} file(s) in the requested range already exist "
                f"in {local_run_dir}.\nThey will be skipped.")
            box.setDetailedText("\n".join(existing))
            box.setStandardButtons(
                QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)
            box.setDefaultButton(QMessageBox.StandardButton.Ok)
            if box.exec() == QMessageBox.StandardButton.Cancel:
                return

        # -- disk space check --
        try:
            needed, free = _check_disk_space(
                remote_host, remote_run_dir, local_base, f_start, f_end)
        except Exception as exc:
            box = QMessageBox(self)
            box.setWindowTitle("Disk Space Check Failed")
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText(f"Could not check remote file sizes:\n{exc}\n\nProceed anyway?")
            box.setStandardButtons(
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            box.setDefaultButton(QMessageBox.StandardButton.No)
            if box.exec() != QMessageBox.StandardButton.Yes:
                return
        else:
            if needed > free:
                box = QMessageBox(self)
                box.setWindowTitle("Insufficient Disk Space")
                box.setIcon(QMessageBox.Icon.Critical)
                box.setText(
                    "Not enough disk space to copy the requested files.\n\n"
                    f"  Required : {_fmt_bytes(needed)}\n"
                    f"  Available: {_fmt_bytes(free)}\n\n"
                    "Free up space and try again.")
                box.setStandardButtons(QMessageBox.StandardButton.Ok)
                box.exec()
                return

        # Combined bash: step 1 = scp, step 2 = analyze
        bash_cmd = (
            f"set -e\n"
            f"\n"
            f"# ── Step 1: Get Data ──────────────────────────────────────\n"
            f"mkdir -p {local_run_dir}\n"
            f"echo '=== Step 1: Copying evio files ==='\n"
            f"echo 'Local directory: {local_run_dir}'\n"
            f"echo 'Listing remote files...'\n"
            f"ALL_FILES=$(ssh {remote_host} 'ls {remote_run_dir}/' 2>/dev/null | sort)\n"
            f"COPIED=0; ALREADY=0\n"
            f"while IFS= read -r f; do\n"
            f"    NUM=$(echo \"$f\" | grep -oP '\\.evio\\.\\K[0-9]+')\n"
            f"    [ -z \"$NUM\" ] && continue\n"
            f"    N=$((10#$NUM))\n"
            f"    if [ \"$N\" -lt {f_start} ] || [ \"$N\" -gt {f_end} ]; then continue; fi\n"
            f"    if [ -f \"{local_run_dir}/$f\" ]; then\n"
            f"        echo \"  Already exists: $f (skipping)\"\n"
            f"        ALREADY=$((ALREADY+1))\n"
            f"    else\n"
            f"        echo \"  Copying $f\"\n"
            f"        scp {remote_host}:{remote_run_dir}/$f {local_run_dir}/\n"
            f"        COPIED=$((COPIED+1))\n"
            f"    fi\n"
            f"done <<< \"$ALL_FILES\"\n"
            f"echo \"Step 1 done. Copied $COPIED file(s), $ALREADY already present.\"\n"
            f"\n"
            f"# ── Step 2: Analyze Data ──────────────────────────────────\n"
            f"echo ''\n"
            f"echo '=== Step 2: Running gain monitor analysis ==='\n"
            f"cd {script_dir}\n"
            f"INPUTDIR={local_base} OUTPUTDIR={outdir} bash {_SCRIPT_PATH} {run_num} {n_cpu}\n"
        )

        self._console.clear()
        self._append(
            f"<span style='color:#8b949e'>Run {run_num} | files {f_start}–{f_end}"
            f" | {n_cpu} CPUs</span><br>")
        self._run_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._process.start("bash", ["-c", bash_cmd])

    def _on_stop(self):
        self._process.kill()

    def _on_stdout(self):
        data = self._process.readAllStandardOutput().data().decode(errors="replace")
        self._append(data.replace("\n", "<br>"))

    def _on_stderr(self):
        data = self._process.readAllStandardError().data().decode(errors="replace")
        self._append(f"<span style='color:#f85149'>{data.replace(chr(10), '<br>')}</span>")

    def _on_finished(self, exit_code, exit_status):
        self._run_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)
        color = "#3fb950" if exit_code == 0 else "#f85149"
        self._append(f"<span style='color:{color}'>[Process finished with exit code {exit_code}]</span>")

    def _append(self, html: str):
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)
        self._console.insertHtml(themed(html))
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)

    def closeEvent(self, event):
        if self._process.state() != QProcess.ProcessState.NotRunning:
            self._process.kill()
            self._process.waitForFinished(2000)
        super().closeEvent(event)


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
        self._palette_idx = PALETTE_NAMES.index(_DEFAULT_PALETTE)
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
        self._hist_cache: OrderedDict = OrderedDict()  # (run_number, module) -> (values, edges)
        self._auto_refresh_timer = QTimer(self)
        self._auto_refresh_timer.timeout.connect(self._auto_refresh_check)

        self._load_geometry()
        self._build_ui()
        legend_map = {0: "gain", 1: "deviation", 2: "drift", 3: "summary"}
        self._map.set_legend_mode(legend_map.get(self._view_mode))
        if self._view_mode == 2:
            self._map.set_palette_override(DRIFT_PALETTE)

    def _load_geometry(self):
        self._all_modules = load_modules(MODULES_JSON)
        self._mod_by_name = {m.name: m for m in self._all_modules}

    @property
    def _active_runs(self) -> List[RunData]:
        return self._runs[self._start_run_idx:self._end_run_idx + 1]

    # ---- UI ----

    def _build_ui(self):
        self.setWindowTitle("HyCal Gain Monitor")
        self.resize(1800, 1000)   # 18:10
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL GAIN MONITOR")
        lbl.setFont(QFont("Consolas", 14, QFont.Weight.Bold))
        lbl.setStyleSheet(themed("color:#58a6ff;"))
        top.addWidget(lbl)
        top.addStretch()

        self._process_btn = self._make_btn(
            "Process Folder...", THEME.ACCENT, self._on_process_folder)
        top.addWidget(self._process_btn)

        self._analyze_btn = self._make_btn(
            "Analyze Data", THEME.SUCCESS, self._on_analyze_data)
        top.addWidget(self._analyze_btn)

        self._getdata_btn = self._make_btn(
            "Get Data", THEME.WARN, self._on_get_data)
        top.addWidget(self._getdata_btn)

        self._doitall_btn = self._make_btn(
            "Do It All", "#bc8cff", self._on_do_it_all)
        top.addWidget(self._doitall_btn)

        self._refresh_btn = self._make_btn(
            "Refresh", THEME.SUCCESS, self._on_refresh)
        self._refresh_btn.setEnabled(False)
        top.addWidget(self._refresh_btn)

        self._auto_refresh_btn = QPushButton("Auto-Refresh: ON")
        self._auto_refresh_btn.setCheckable(True)
        self._auto_refresh_btn.setChecked(True)
        self._auto_refresh_btn.setEnabled(False)
        self._auto_refresh_btn.setFont(QFont("Consolas", 10))
        self._auto_refresh_btn.setFixedHeight(28)
        # Themed via f-string so the :checked state's green tint comes from
        # THEME.SUCCESS / THEME.BUTTON_HOVER and tracks dark/light correctly.
        self._auto_refresh_btn.setStyleSheet(
            f"QPushButton{{background:{THEME.PANEL};color:{THEME.TEXT_DIM};"
            f"border:1px solid {THEME.BORDER};border-radius:3px;padding:0 8px;}}"
            f"QPushButton:checked{{background:{THEME.BUTTON_HOVER};"
            f"color:{THEME.SUCCESS};border-color:{THEME.SUCCESS};}}"
            f"QPushButton:hover{{background:{THEME.BUTTON};}}")
        self._auto_refresh_btn.toggled.connect(self._on_auto_refresh_toggled)
        top.addWidget(self._auto_refresh_btn)

        top.addWidget(self._slabel("every"))
        self._auto_refresh_interval = QSpinBox()
        self._auto_refresh_interval.setRange(5, 3600)
        self._auto_refresh_interval.setValue(30)
        self._auto_refresh_interval.setSuffix(" s")
        self._auto_refresh_interval.setFixedWidth(72)
        self._auto_refresh_interval.setFont(QFont("Consolas", 10))
        self._auto_refresh_interval.setStyleSheet(themed(
            "QSpinBox{background:#161b22;color:#c9d1d9;border:1px solid #30363d;"
            "border-radius:3px;padding:2px 4px;}"
            "QSpinBox::up-button,QSpinBox::down-button{width:16px;}"))
        self._auto_refresh_interval.valueChanged.connect(self._on_refresh_interval_changed)
        top.addWidget(self._auto_refresh_interval)

        self._status_lbl = QLabel("No data loaded")
        self._status_lbl.setFont(QFont("Consolas", 11))
        self._status_lbl.setStyleSheet(themed("color:#8b949e;"))
        top.addWidget(self._status_lbl)
        root.addLayout(top)

        # -- body splitter (horizontal: table | map | charts+reserved) --
        body = QSplitter(Qt.Orientation.Horizontal)
        self._body = body

        # ---- left panel: report table ----
        self._irregular_table = IrregularTableWidget()
        self._irregular_table.runClicked.connect(self._on_jump_to_run)
        self._irregular_table.moduleClicked.connect(self._on_module_clicked)
        body.addWidget(self._irregular_table)

        # ---- middle panel: HyCal map ----
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
        self._ref_combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}"))
        self._ref_combo.currentIndexChanged.connect(self._on_ref_changed)
        ctrl.addWidget(self._ref_combo)

        ctrl.addSpacing(10)
        ctrl.addWidget(self._slabel("View:"))
        self._view_combo = QComboBox()
        self._view_combo.addItems(["Gain Factor", "Deviation (σ)", "Run-to-Run Drift", "Summary"])
        self._view_combo.setFixedWidth(150)
        self._view_combo.setFont(QFont("Consolas", 10))
        self._view_combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}"))
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
        self._thresh_g_input.setStyleSheet(themed(_edit_ss))
        self._thresh_g_input.editingFinished.connect(self._on_drift_threshold_changed)
        ctrl.addWidget(self._thresh_g_input)
        self._thresh_g_pct = self._slabel("%")
        ctrl.addWidget(self._thresh_g_pct)

        self._thresh_w_lbl = self._slabel("  W thresh:")
        ctrl.addWidget(self._thresh_w_lbl)

        self._thresh_w_input = QLineEdit("5.0")
        self._thresh_w_input.setFixedWidth(46)
        self._thresh_w_input.setFont(QFont("Consolas", 10))
        self._thresh_w_input.setStyleSheet(themed(_edit_ss))
        self._thresh_w_input.editingFinished.connect(self._on_drift_threshold_changed)
        ctrl.addWidget(self._thresh_w_input)
        self._thresh_w_pct = self._slabel("%")
        ctrl.addWidget(self._thresh_w_pct)

        ctrl.addSpacing(10)
        ctrl.addWidget(self._slabel("Start:"))
        self._start_combo = QComboBox()
        self._start_combo.setMinimumWidth(100)
        self._start_combo.setFont(QFont("Consolas", 10))
        self._start_combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}"))
        self._start_combo.currentIndexChanged.connect(self._on_start_run_changed)
        ctrl.addWidget(self._start_combo)

        ctrl.addSpacing(6)
        ctrl.addWidget(self._slabel("End:"))
        self._end_combo = QComboBox()
        self._end_combo.setMinimumWidth(100)
        self._end_combo.setFont(QFont("Consolas", 10))
        self._end_combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}"))
        self._end_combo.currentIndexChanged.connect(self._on_end_run_changed)
        ctrl.addWidget(self._end_combo)

        ctrl.addSpacing(10)
        ctrl.addWidget(self._slabel("Run:"))

        self._run_combo = QComboBox()
        self._run_combo.setMinimumWidth(100)
        self._run_combo.setFont(QFont("Consolas", 10))
        self._run_combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox::drop-down{border:none;width:18px;}"
            "QComboBox::down-arrow{border-left:4px solid transparent;"
            "border-right:4px solid transparent;border-top:5px solid #8b949e;"
            "margin-right:4px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;selection-background-color:#1f6feb;}"))
        self._run_combo.currentIndexChanged.connect(self._on_run_changed)
        ctrl.addWidget(self._run_combo)

        self._prev_btn = self._make_btn("<", THEME.TEXT, self._on_prev_run)
        self._prev_btn.setFixedWidth(30)
        ctrl.addWidget(self._prev_btn)

        self._next_btn = self._make_btn(">", THEME.TEXT, self._on_next_run)
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
        self._range_min.setStyleSheet(themed(_EDIT_SS))
        self._range_min.returnPressed.connect(self._on_apply_range)
        rng.addWidget(self._range_min)

        rng.addWidget(self._slabel("Max:"))
        self._range_max = QLineEdit("1.1")
        self._range_max.setFixedWidth(70)
        self._range_max.setFont(QFont("Consolas", 10))
        self._range_max.setStyleSheet(themed(_EDIT_SS))
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
        self._apply_btn.setStyleSheet(themed(
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 8px;"
            "font:bold 11px Consolas;border-radius:3px;}"
            "QPushButton:hover{background:#30363d;}"))
        self._log_btn.setStyleSheet(themed(_TOGGLE_SS))
        self._auto_btn.setStyleSheet(themed(_TOGGLE_SS))

        rng.addStretch()
        left_layout.addLayout(rng)

        # geo map
        self._map = HyCalGainMapWidget()
        self._map.set_modules(self._all_modules)
        self._map.set_palette(self._palette_idx)
        self._map.moduleHovered.connect(self._on_module_hovered)
        self._map.moduleClicked.connect(self._on_module_clicked)
        self._map.paletteClicked.connect(self._on_cycle_palette)
        left_layout.addWidget(self._map, stretch=1)

        # info label
        self._info = QLabel("Hover over a module for details")
        self._info.setFont(QFont("Consolas", 10))
        self._info.setStyleSheet(themed(
            "QLabel{background:#161b22;color:#c9d1d9;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}"))
        self._info.setFixedHeight(26)
        left_layout.addWidget(self._info)

        body.addWidget(left)

        # ---- right panel (vertical splitter: charts top, reserved bottom) ----
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
        self._hist_widget = RootHistWidget()
        right.addWidget(self._hist_widget)
        right.setStretchFactor(0, 3)
        right.setStretchFactor(1, 2)

        body.addWidget(right)
        body.setStretchFactor(0, 1)  # table
        body.setStretchFactor(1, 1)  # HyCal map
        body.setStretchFactor(2, 5)  # charts + reserved
        QTimer.singleShot(0, lambda: self._body.setSizes([200, 350, 1050]))

        root.addWidget(body, stretch=1)

    # ---- helpers ----

    def _make_btn(self, text, fg, slot):
        btn = QPushButton(text)
        btn.setStyleSheet(themed(
            f"QPushButton{{background:#21262d;color:{fg};"
            f"border:1px solid #30363d;padding:4px 12px;"
            f"font:bold 11px Consolas;border-radius:3px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{color:#555;}}"))
        btn.clicked.connect(slot)
        return btn

    def _slabel(self, text):
        lbl = QLabel(text)
        lbl.setFont(QFont("Consolas", 10))
        lbl.setStyleSheet(themed("color:#c9d1d9;"))
        return lbl

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

        # Only parse files that are new or whose mtime changed — avoids re-reading
        # the entire folder on every auto-refresh tick when the dataset is large.
        changed_runs = {rnum for rnum, mtime in new_snapshot.items()
                        if self._file_snapshot.get(rnum) != mtime}
        if not changed_runs:
            return
        runs_by_number = {rd.run_number: rd for rd in self._runs}
        for rnum in changed_runs:
            path = os.path.join(self._current_folder,
                                f"prad_{rnum:06d}_LMS.dat")
            rd = parse_dat_file(path)
            if rd is not None:
                runs_by_number[rnum] = rd
        # Remove runs whose files have disappeared from the snapshot
        for rnum in list(runs_by_number):
            if rnum not in new_snapshot:
                del runs_by_number[rnum]
        new_runs = sorted(runs_by_number.values(), key=lambda r: r.run_number)
        if not new_runs:
            return
        self._runs = new_runs
        self._file_snapshot = new_snapshot
        self._deviation_stats_key = None
        self._hist_cache.clear()

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
        self._status_lbl.setStyleSheet(themed("color:#3fb950;"))
        self._update_all_views()

    def _on_process_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select Data Folder")
        if not folder:
            return
        self._load_folder(folder)

    def _on_analyze_data(self):
        dlg = AnalyzeDialog(self)
        dlg.exec()

    def _on_get_data(self):
        dlg = GetDataDialog(self)
        dlg.exec()

    def _on_do_it_all(self):
        dlg = DoItAllDialog(self)
        dlg.exec()

    def _load_folder(self, folder: str):
        self._status_lbl.setText("Loading...")
        self._status_lbl.setStyleSheet(themed("color:#d29922;"))
        QApplication.processEvents()

        self._runs = load_all_runs(folder)
        self._current_folder = folder
        if not self._runs:
            self._status_lbl.setText("No data files found")
            self._status_lbl.setStyleSheet(themed("color:#f85149;"))
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
        self._status_lbl.setStyleSheet(themed("color:#3fb950;"))
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
        self._load_lms_hist(self._selected_module)

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
        else:
            self._load_lms_hist(self._selected_module, run_number)

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
            self._hist_widget.clear()
        else:
            self._update_geo_view()
            if not (prev_mode in (0, 1) and index in (0, 1)):
                self._update_irregular_table()
            self._load_lms_hist(self._selected_module)

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
        self._load_lms_hist(self._selected_module)

    def _load_lms_hist(self, module_name: Optional[str],
                       run_number: Optional[int] = None):
        # In summary mode a run_number must be supplied explicitly (from chart click).
        if self._view_mode == 3 and run_number is None:
            self._hist_widget.clear()
            return
        if not module_name or not self._current_folder or not self._runs:
            self._hist_widget.clear()
            return
        if not _UPROOT_OK:
            self._hist_widget.clear()
            return
        if run_number is None:
            run_number = self._runs[self._current_run_idx].run_number
        key = (run_number, module_name)
        if key in self._hist_cache:
            self._hist_cache.move_to_end(key)
            values, edges, gauss = self._hist_cache[key]
        else:
            root_path = os.path.join(self._current_folder,
                                     f"prad_{run_number:06d}_LMS_fitted.root")
            if not os.path.exists(root_path):
                self._hist_widget.clear()
                return
            try:
                with uproot.open(root_path) as rf:
                    hist_key = f"{module_name}_LMS"
                    if hist_key not in rf:
                        self._hist_widget.clear()
                        return
                    h = rf[hist_key]
                    values = h.values()
                    edges = h.axis().edges()
                    gauss = None
                    fns = h.member("fFunctions")
                    for fn in fns:
                        if fn.classname == "TF1":
                            params = fn.member("fFormula").member("fClingParameters")
                            if len(params) >= 3:
                                gauss = (float(params[0]), float(params[1]),
                                         float(params[2]),
                                         float(fn.member("fXmin")),
                                         float(fn.member("fXmax")))
                            break
            except Exception:
                self._hist_widget.clear()
                return
            self._hist_cache[key] = (values, edges, gauss)
            if len(self._hist_cache) > _HIST_CACHE_MAX:
                self._hist_cache.popitem(last=False)
        self._hist_widget.set_histogram(values, edges,
                                        f"{module_name}  run {run_number}",
                                        gauss=gauss)

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
                if gains:
                    self._charts[i].set_y_range(*_chart_y_range(gains, []))
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
                if ratios:
                    self._charts[i].set_y_range(*_chart_y_range(ratios, errors))
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
    parser.add_argument("--theme", choices=available_themes(), default="dark",
                        help="Colour theme (default: dark)")
    args, qt_args = parser.parse_known_args()

    set_theme(args.theme)

    app = QApplication([sys.argv[0]] + qt_args)
    win = GainMonitorWindow()
    if args.folder:
        win._process_btn.setEnabled(False)
        win._load_folder(args.folder)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
