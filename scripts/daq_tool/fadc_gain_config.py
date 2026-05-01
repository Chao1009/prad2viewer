#!/usr/bin/env python3
"""
FADC Gain Config Editor (PyQt6)
================================
Interactive HyCal geo-view editor for the FADC250 gain config
(``adchycal_gain.cnf``).  Always opens a GUI: load a calibration JSON,
edit per-channel gains by clicking or dragging on the HyCal map, switch
to *Closing* mode to null gains (effectively masking channels off), and
save the resulting trigger config.

Workflow inside the GUI
-----------------------
* **Load Calibration…** – read a JSON list of ``{name, factor, ...}``
  entries; gains for matched modules become the displayed colormap.
* **Paint mode**:
    * *Set gain* – left-click / drag paints the spinbox value.
    * *Closing (null = 0)* – left-click / drag zeroes channels (mask off).
  Drag the mouse with the left button held to bulk-paint a region.
* **Safe Cap** – clamp every channel's gain to a configurable
  ``[Min, Max]`` (defaults ``[0.0, 0.15]``); gains stay non-negative.
* **Reset** – revert all manual edits to the loaded base.
* **Load .cnf… / Save .cnf…** – open / save an ``adchycal_gain.cnf`` file.

Optional CLI shortcuts (still always open the GUI):
    python fadc_gain_config.py
    python fadc_gain_config.py -c database/calibration/adc_to_mev_factors_cosmic.json
    python fadc_gain_config.py -o /path/to/adchycal_gain.cnf
    python fadc_gain_config.py -i existing.cnf
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent

NUM_CRATES = 7
CRATE_NAMES = [f"adchycal{i}" for i in range(1, NUM_CRATES + 1)]
CHANNELS_PER_SLOT = 16

DEFAULT_UNMAPPED_GAIN = 1.0
DEFAULT_LMS_GAIN = 1.0
DEFAULT_SCINT_GAIN = 1.0


# ---------------------------------------------------------------------------
#  Database auto-discovery
# ---------------------------------------------------------------------------

def find_database_dir(explicit: Optional[str] = None) -> Path:
    if explicit:
        p = Path(explicit).resolve()
        if not p.is_dir():
            sys.exit(f"error: --database path does not exist: {p}")
        return p

    candidates = [
        SCRIPT_DIR / ".." / ".." / "database",
        Path.cwd() / "database",
        Path.cwd(),
    ]
    for c in candidates:
        if (c / "hycal_daq_map.json").is_file() and (c / "hycal_modules.json").is_file():
            return c.resolve()
    sys.exit("error: could not locate database directory "
             "(looked for hycal_daq_map.json + hycal_modules.json)")


def load_modules(db_dir: Path) -> Dict[str, str]:
    """Return {module_name: module_type} for all HyCal modules."""
    with open(db_dir / "hycal_modules.json") as f:
        mods = json.load(f)
    return {m["n"]: m["t"] for m in mods}


def load_daq_map(db_dir: Path) -> List[Tuple[str, int, int, int]]:
    """Return list of (name, crate, slot, channel)."""
    with open(db_dir / "hycal_daq_map.json") as f:
        entries = json.load(f)
    return [(e["name"], e["crate"], e["slot"], e["channel"]) for e in entries]


# ---------------------------------------------------------------------------
#  Gain source
# ---------------------------------------------------------------------------

def load_calibration(path: Path) -> Dict[str, float]:
    """Return {module_name: gain_factor} from a calibration JSON file."""
    with open(path) as f:
        data = json.load(f)
    out: Dict[str, float] = {}
    for entry in data:
        name = entry.get("name")
        factor = entry.get("factor")
        if name is None or factor is None:
            continue
        out[name] = float(factor)
    return out


def resolve_gain(name: str,
                 mod_type: Optional[str],
                 cal: Dict[str, float],
                 pbwo4_gain: float,
                 pbglass_gain: float) -> float:
    if name in cal:
        return cal[name]
    if mod_type == "PbWO4":
        return pbwo4_gain
    if mod_type == "PbGlass":
        return pbglass_gain
    if mod_type == "LMS":
        return DEFAULT_LMS_GAIN
    # V1-V4 scintillators and anything else
    return DEFAULT_SCINT_GAIN


def safe_cap_gains(gains: Dict[str, float],
                   min_cap: float = 0.0,
                   max_cap: float = 0.15) -> Dict[str, float]:
    """Clamp each gain in ``gains`` to ``[max(0, min_cap), max_cap]``.

    The lower bound is floored at 0 — gain factors are physical and must
    never be negative.  Raises ``ValueError`` if ``max_cap`` is below the
    (non-negative) lower bound.
    """
    lo = max(0.0, float(min_cap))
    hi = float(max_cap)
    if hi < lo:
        raise ValueError(f"max_cap ({hi}) must be >= min_cap ({lo})")
    return {k: lo if v < lo else (hi if v > hi else v)
            for k, v in gains.items()}


# ---------------------------------------------------------------------------
#  Config text rendering / parsing
# ---------------------------------------------------------------------------

def format_gain(g: float) -> str:
    return f"{g:.6f}"


def render_cnf(daq: List[Tuple[str, int, int, int]],
               gains_by_name: Dict[str, float],
               header_comments: Optional[List[str]] = None,
               ) -> Tuple[str, int]:
    """Build the ``.cnf`` text from a per-module gain dict.

    Returns ``(text, num_unmapped_channels)``.  Channels for which no
    module is present in the DAQ map at the given (crate, slot, channel)
    get the fallback :data:`DEFAULT_UNMAPPED_GAIN`.
    """
    slots: Dict[Tuple[int, int], Dict[int, Tuple[str, float]]] = {}
    for name, crate, slot, ch in daq:
        if crate < 0 or slot < 0 or ch < 0:
            continue
        gain = gains_by_name.get(name, DEFAULT_UNMAPPED_GAIN)
        slots.setdefault((crate, slot), {})[ch] = (name, gain)

    lines: List[str] = ["# adchycal_gain.cnf",
                        "# Generated by fadc_gain_config.py"]
    if header_comments:
        lines.extend(header_comments)
    lines.append("")

    unmapped = 0
    for ci in range(NUM_CRATES):
        crate_slots = sorted(s for (c, s) in slots if c == ci)
        if not crate_slots:
            continue
        lines.append(f"FAV3_CRATE {CRATE_NAMES[ci]}")
        for slot in crate_slots:
            ch_map = slots[(ci, slot)]
            gains: List[str] = []
            names: List[str] = []
            for ch in range(CHANNELS_PER_SLOT):
                entry = ch_map.get(ch)
                if entry is None:
                    gains.append(format_gain(DEFAULT_UNMAPPED_GAIN))
                    names.append(f"ch{ch}:unmapped")
                    unmapped += 1
                else:
                    name, g = entry
                    gains.append(format_gain(g))
                    names.append(name)
            lines.append(f"# slot {slot}: {', '.join(names)}")
            lines.append(f"FAV3_SLOT {slot}")
            lines.append(f"FAV3_ALLCH_GAIN {' '.join(gains)}")
        lines.append("FAV3_CRATE end")
        lines.append("")

    return "\n".join(lines), unmapped


def parse_cnf_text(text: str,
                   daq: List[Tuple[str, int, int, int]]) -> Dict[str, float]:
    """Parse ``.cnf`` text, return ``{module_name: gain}`` for mapped channels."""
    name_to_ci = {n: i for i, n in enumerate(CRATE_NAMES)}
    daq_lookup: Dict[Tuple[int, int, int], str] = {}
    for name, crate, slot, ch in daq:
        if crate >= 0 and slot >= 0 and ch >= 0:
            daq_lookup[(crate, slot, ch)] = name

    gains: Dict[str, float] = {}
    current_crate_idx = -1
    current_slot = -1

    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        kw = parts[0]
        if kw == "FAV3_CRATE":
            cname = parts[1] if len(parts) > 1 else ""
            current_crate_idx = (
                -1 if cname == "end" else name_to_ci.get(cname, -1))
        elif kw == "FAV3_SLOT" and len(parts) > 1:
            try:
                current_slot = int(parts[1])
            except ValueError:
                current_slot = -1
        elif (kw == "FAV3_ALLCH_GAIN"
              and current_crate_idx >= 0 and current_slot >= 0):
            for ch, val in enumerate(parts[1:1 + CHANNELS_PER_SLOT]):
                name = daq_lookup.get((current_crate_idx, current_slot, ch))
                if name:
                    try:
                        gains[name] = float(val)
                    except ValueError:
                        pass
    return gains


# ---------------------------------------------------------------------------
#  GUI: HyCal geo-view editor
# ---------------------------------------------------------------------------

# LMS / V module display positions, slotted in a row below HyCal so that
# they remain visible on the geo-view (matches trigger_mask_editor).
_BOTTOM_Y = -640.0
_BOTTOM_SZ = 50.0
_LMS_V_XPOS = {
    "LMS1": -200.0, "LMS2": -145.0, "LMS3": -90.0, "LMSP": -35.0,
    "V1":     35.0, "V2":     90.0, "V3":   145.0, "V4":  200.0,
}
_LABEL_NAMES = set(_LMS_V_XPOS.keys())


class _ModuleInfo:
    __slots__ = ("name", "mod_type", "x", "y", "sx", "sy",
                 "crate", "slot", "channel")

    def __init__(self, name, mod_type, x, y, sx, sy,
                 crate=-1, slot=-1, channel=-1):
        self.name = name
        self.mod_type = mod_type
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy
        self.crate = crate
        self.slot = slot
        self.channel = channel


def _load_module_info(db_dir: Path) -> List[_ModuleInfo]:
    """Load HyCal modules joined with the DAQ map for the GUI."""
    with open(db_dir / "hycal_modules.json") as f:
        mods_json = json.load(f)
    with open(db_dir / "hycal_daq_map.json") as f:
        daq_json = json.load(f)

    daq_by_name: Dict[str, Tuple[int, int, int]] = {
        d["name"]: (d["crate"], d["slot"], d["channel"]) for d in daq_json
    }

    modules: List[_ModuleInfo] = []
    for m in mods_json:
        name = m["n"]
        if name in _LMS_V_XPOS:
            x, y, sx, sy = _LMS_V_XPOS[name], _BOTTOM_Y, _BOTTOM_SZ, _BOTTOM_SZ
        else:
            x, y, sx, sy = m["x"], m["y"], m["sx"], m["sy"]
        crate, slot, ch = daq_by_name.get(name, (-1, -1, -1))
        modules.append(_ModuleInfo(name, m["t"], x, y, sx, sy, crate, slot, ch))

    have = {m.name for m in modules}
    for name in _LMS_V_XPOS:
        if name in have:
            continue
        crate, slot, ch = daq_by_name.get(name, (-1, -1, -1))
        if crate >= 0:
            t = "Scintillator" if name.startswith("V") else "LMS"
            modules.append(_ModuleInfo(name, t, _LMS_V_XPOS[name],
                                       _BOTTOM_Y, _BOTTOM_SZ, _BOTTOM_SZ,
                                       crate, slot, ch))
    return modules


# PyQt + hycal_geoview imports (script is always GUI-only).
sys.path.insert(0, str(SCRIPT_DIR.parent))
from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QColor, QPen, QFont, QDoubleValidator
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QTextEdit, QSplitter, QFileDialog,
    QDoubleSpinBox, QGroupBox, QFormLayout, QToolTip, QInputDialog,
    QMessageBox,
)
from hycal_geoview import (
    HyCalMapWidget as _HyCalMapBase,
    Module as _Module,
    ColorRangeControl,
    THEME, apply_theme_palette, set_theme, available_themes, themed,
    cmap_qcolor,
)


def _btn_style(checked_color: Optional[str] = None) -> str:
    base = themed(
        f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT};"
        f"border:1px solid {THEME.BORDER};padding:6px 14px;"
        f"font:10pt;border-radius:8px;}}"
        f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")
    if checked_color:
        base += themed(
            f"QPushButton:checked{{background:{checked_color};"
            f"color:{THEME.TEXT};border:1px solid {checked_color};}}")
    return base


# ---------------------------------------------------------------------------
#  HyCal geo-view widget
# ---------------------------------------------------------------------------
# Three interaction modes (selected via the editor's right-panel buttons):
#   * Edit (default): click on a module emits ``moduleEditRequested`` with
#     the current gain.  The editor opens a popup dialog to set a new value.
#   * Mask: drag-paint zeros (channel closed).
#   * Set:  drag-paint the value carried by ``_paint_value``.
# At drag end ``paintCommitted`` fires with the batch of
# ``(name, prior_override_or_None)`` tuples so the editor can record the
# action on its undo stack.

PAINT_MODE_EDIT = "edit"
PAINT_MODE_MASK = "mask"
PAINT_MODE_SET = "set"


class _HyCalGainMap(_HyCalMapBase):
    moduleEditRequested = pyqtSignal(str, float)   # name, current gain
    paintCommitted = pyqtSignal(list)              # [(name, prior_value_or_None), ...]

    def __init__(self, modules: List[_ModuleInfo], parent=None):
        super().__init__(parent, shrink=0.92, margin_top=10,
                         margin_bottom=40, include_lms=True,
                         show_colorbar=True, min_size=(500, 500))
        self._mod_map: Dict[str, _ModuleInfo] = {m.name: m for m in modules}
        self._gains: Dict[str, float] = {}
        self._overrides: Dict[str, float] = {}

        self._paint_mode: str = PAINT_MODE_EDIT
        self._paint_value: float = 0.0
        self._paint_dragging = False
        self._drag_visited: Set[str] = set()
        self._drag_batch: List[Tuple[str, Optional[float]]] = []

        base_modules = [_Module(m.name, m.mod_type, m.x, m.y, m.sx, m.sy)
                        for m in modules]
        self.set_modules(base_modules)
        self.set_range(0.0, 1.0)

    # ---- public API ----

    def set_paint_mode(self, mode: str) -> None:
        if mode not in (PAINT_MODE_EDIT, PAINT_MODE_MASK, PAINT_MODE_SET):
            return
        self._paint_mode = mode
        self._paint_dragging = False
        self._drag_visited.clear()
        self._drag_batch = []
        self.setCursor(
            Qt.CursorShape.CrossCursor if mode != PAINT_MODE_EDIT
            else Qt.CursorShape.ArrowCursor)

    def set_paint_value(self, v: float) -> None:
        self._paint_value = float(v)

    @property
    def paint_mode(self) -> str:
        return self._paint_mode

    def set_gains(self, gains: Dict[str, float],
                  overrides: Optional[Dict[str, float]] = None) -> None:
        self._gains = dict(gains)
        self._overrides = dict(overrides) if overrides is not None else {}
        # Share the dict with the base widget's _values so the colormap
        # picks up our edits without separate set_values calls.
        self.set_values(self._gains)
        self.update()

    @property
    def gains(self) -> Dict[str, float]:
        return self._gains

    @property
    def overrides(self) -> Dict[str, float]:
        return self._overrides

    # ---- painting ----

    def _paint_modules(self, p):
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        null_color = QColor(THEME.DANGER)
        vmin, vmax = self._vmin, self._vmax
        for name, rect in self._rects.items():
            m = self._mod_map.get(name)
            if m is None or m.crate < 0:
                p.fillRect(rect, no_data)
                continue
            v = self._gains.get(name, 0.0)
            if v == 0.0:
                p.fillRect(rect, null_color)
            else:
                t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                t = max(0.0, min(1.0, t))
                p.fillRect(rect, cmap_qcolor(t, stops))

    def _paint_overlays(self, p, w, h):
        # White border around modules whose gain was set via the GUI.
        sel_pen = QPen(QColor(THEME.SELECT_BORDER), 1.5)
        sel_pen.setCosmetic(True)
        p.setPen(sel_pen)
        p.setBrush(Qt.BrushStyle.NoBrush)
        for name in self._overrides:
            rect = self._rects.get(name)
            if rect is not None:
                p.drawRect(rect)
        # LMS / V text labels
        p.setPen(QColor(THEME.TEXT))
        p.setFont(QFont("Monospace", 7, QFont.Weight.Bold))
        for name in _LABEL_NAMES:
            r = self._rects.get(name)
            if r is not None:
                p.drawText(r, Qt.AlignmentFlag.AlignCenter, name)
        # Hover highlight
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(QColor(THEME.ACCENT), 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

    def _paint_after_colorbar(self, p, w, h):
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 9))
        n_masked = sum(1 for v in self._gains.values() if v == 0.0)
        info = f"Edits: {len(self._overrides)}    Masked: {n_masked}"
        if self._paint_mode == PAINT_MODE_MASK:
            info += "    [MASK]"
        elif self._paint_mode == PAINT_MODE_SET:
            info += f"    [SET={self._paint_value:.6g}]"
        p.drawText(QRectF(8, h - 18, w - 16, 16),
                   Qt.AlignmentFlag.AlignLeft, info)

    # ---- hit / mouse ----

    def _hit(self, pos) -> Optional[str]:
        for name in self._rect_names_rev:
            if self._rects[name].contains(pos):
                m = self._mod_map.get(name)
                if m and m.crate >= 0:
                    return name
        return None

    def _apply_paint(self, name: str) -> None:
        """Apply the current paint mode's value to ``name``.  Records the
        prior override (or ``None``) onto the active drag batch so the
        editor can undo the action."""
        m = self._mod_map.get(name)
        if not m or m.crate < 0:
            return
        if self._paint_mode == PAINT_MODE_MASK:
            v = 0.0
        elif self._paint_mode == PAINT_MODE_SET:
            v = self._paint_value
        else:
            return
        # No-op if the cell already holds this exact override.
        if (name in self._overrides
                and self._overrides[name] == v
                and self._gains.get(name) == v):
            return
        prior = self._overrides.get(name)
        self._drag_batch.append((name, prior))
        self._gains[name] = v
        self._overrides[name] = v

    def _tooltip_text(self, name: str) -> str:
        m = self._mod_map.get(name)
        v = self._gains.get(name, 0.0)
        tip = f"{name}: {v:.6g}"
        if v == 0.0:
            tip += "  [masked]"
        if name in self._overrides:
            tip += "  (edit)"
        if m and m.crate >= 0:
            tip += f"\ncrate={m.crate} slot={m.slot} ch={m.channel}"
        if self._paint_mode == PAINT_MODE_EDIT:
            tip += "\n(click to edit gain)"
        elif self._paint_mode == PAINT_MODE_SET:
            tip += f"\n(click/drag to set {self._paint_value:.6g})"
        else:
            tip += "\n(click/drag to mask)"
        return tip

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        pos = event.position()
        # Defer to the base inline-range edit feature first — this would
        # otherwise be swallowed by our paint-mode dispatch below.
        if self._check_inline_range_edit_click(pos):
            return
        # Click on the colour bar cycles palettes (base widget feature)
        if self._cb_rect and self._cb_rect.contains(pos):
            self.cycle_palette()
            return
        found = self._hit(pos)
        if not found:
            return
        if self._paint_mode == PAINT_MODE_EDIT:
            # Editor will show a popup dialog
            self.moduleEditRequested.emit(found, self._gains.get(found, 0.0))
        else:
            self._paint_dragging = True
            self._drag_visited = {found}
            self._drag_batch = []
            self._apply_paint(found)
            self.update()

    def mouseMoveEvent(self, event):
        pos = event.position()
        found = self._hit(pos)
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                QToolTip.showText(event.globalPosition().toPoint(),
                                  self._tooltip_text(found), self)
                self.moduleHovered.emit(found)
            else:
                QToolTip.hideText()
        if (self._paint_mode != PAINT_MODE_EDIT and self._paint_dragging
                and found and found not in self._drag_visited):
            self._drag_visited.add(found)
            self._apply_paint(found)
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        if self._paint_dragging and self._drag_batch:
            self.paintCommitted.emit(list(self._drag_batch))
        self._paint_dragging = False
        self._drag_visited.clear()
        self._drag_batch = []

    def wheelEvent(self, event):
        event.ignore()


class _GainEditor(QMainWindow):
    def __init__(self,
                 modules: List[_ModuleInfo],
                 daq: List[Tuple[str, int, int, int]],
                 mod_types: Dict[str, str],
                 db_dir: Path,
                 cal: Optional[Dict[str, float]] = None,
                 cal_path: Optional[Path] = None,
                 pbwo4_gain: float = 1.0,
                 pbglass_gain: float = 1.0,
                 initial_overrides: Optional[Dict[str, float]] = None,
                 output_path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("FADC Gain Editor")
        self.resize(1100, 800)

        self._modules = modules
        self._daq = daq
        self._mod_types = mod_types
        self._db_dir = db_dir
        self._cal: Dict[str, float] = dict(cal or {})
        self._cal_path: Optional[Path] = cal_path
        self._pbwo4 = pbwo4_gain
        self._pbglass = pbglass_gain
        self._output_path = output_path
        self._mod_map: Dict[str, _ModuleInfo] = {m.name: m for m in modules}

        apply_theme_palette(self)
        # Window-scoped stylesheet so QLabel / QRadioButton / QDoubleSpinBox
        # (which don't reliably pick up the QPalette on Windows native style)
        # render text against the dark surfaces correctly.
        self.setStyleSheet(themed(
            f"QLabel{{color:{THEME.TEXT};background:transparent;}}"
            f"QDoubleSpinBox,QSpinBox,QLineEdit{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT};border:1px solid {THEME.BORDER};"
            f"border-radius:4px;padding:2px 6px;}}"
            f"QGroupBox{{color:{THEME.TEXT};background:transparent;"
            f"border:1px solid {THEME.BORDER};border-radius:6px;"
            f"margin-top:10px;padding-top:10px;}}"
            f"QGroupBox::title{{subcontrol-origin:margin;left:10px;"
            f"padding:0 6px;color:{THEME.TEXT_DIM};}}"
            f"QSplitter::handle{{background:{THEME.BORDER};}}"
        ))

        self._base_gains = self._compute_base_gains()
        overrides = self._diff_overrides(initial_overrides or {})

        # Undo stack: each entry is a batch (list) of (name, prior_value_or_None)
        # tuples.  ``prior_value=None`` means "module had no override before"
        # — undoing reverts it to the base gain.
        self._history: List[List[Tuple[str, Optional[float]]]] = []

        # Colormap range control built later in _build_right_panel; until
        # then any _notify_range_values() call is a no-op.
        self._range_ctrl: Optional[ColorRangeControl] = None

        self._map = _HyCalGainMap(modules)
        merged = dict(self._base_gains)
        merged.update(overrides)
        self._map.set_gains(merged, overrides)
        # Pre-fit so the colorbar isn't 0..1 at first paint.
        self._fit_initial_range()

        self._build_right_panel()

        self._map.moduleEditRequested.connect(self._on_module_edit_requested)
        self._map.paintCommitted.connect(self._on_paint_committed)
        self._map.moduleHovered.connect(self._on_hover)

        self._status = QLabel("Click or drag modules to apply gain")
        self._status.setStyleSheet(themed(
            f"color:{THEME.TEXT_DIM};font:10pt Monospace;padding:4px;"))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._map)
        splitter.addWidget(self._right)
        splitter.setSizes([700, 400])

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter, 1)
        layout.addWidget(self._status)
        self.setCentralWidget(central)

        self._update_cal_label()
        self._refresh_text()

    # ---- helpers ----

    def _compute_base_gains(self) -> Dict[str, float]:
        base: Dict[str, float] = {}
        for name, crate, slot, ch in self._daq:
            if crate < 0 or slot < 0 or ch < 0:
                continue
            mt = self._mod_types.get(name)
            base[name] = resolve_gain(name, mt, self._cal,
                                      self._pbwo4, self._pbglass)
        return base

    def _diff_overrides(self,
                        candidate: Dict[str, float]) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for name, val in candidate.items():
            base = self._base_gains.get(name)
            if base is None or abs(base - val) > 1e-12:
                out[name] = val
        return out

    def _fit_initial_range(self) -> None:
        """One-time min/max-of-nonzero fit before the range control exists."""
        vals = [v for v in self._map.gains.values() if v > 0]
        if vals:
            vmin = min(vals)
            vmax = max(vals)
            if vmin == vmax:
                vmax = vmin + max(abs(vmin) * 0.1, 1e-3)
            self._map.set_range(vmin, vmax)
        else:
            self._map.set_range(0.0, 1.0)

    def _notify_range_values(self) -> None:
        """Tell the range control the gain dict changed.  Re-fits if the
        Auto button is in persistent (pinned) mode; otherwise no-op."""
        if self._range_ctrl is not None:
            self._range_ctrl.notify_values_changed(self._map.gains)

    def _rebuild_from_base(self,
                           keep_overrides: bool = True) -> None:
        """Recompute base gains and re-merge with overrides on top."""
        self._base_gains = self._compute_base_gains()
        overrides = self._map.overrides if keep_overrides else {}
        merged = dict(self._base_gains)
        merged.update(overrides)
        self._map.set_gains(merged, overrides)
        self._notify_range_values()
        self._refresh_text()

    # ---- right panel ----

    def _build_right_panel(self) -> None:
        self._right = QWidget()
        v = QVBoxLayout(self._right)
        v.setContentsMargins(8, 8, 8, 8)
        v.setSpacing(8)

        # ---- Source group ----
        src_grp = QGroupBox("Source")
        sform = QFormLayout(src_grp)

        cal_row = QHBoxLayout()
        self._cal_label = QLabel("(no calibration loaded)")
        self._cal_label.setStyleSheet(themed(
            f"color:{THEME.TEXT_DIM};font:9pt Monospace;"))
        cal_row.addWidget(self._cal_label, 1)
        btn_load_cal = QPushButton("Load Calibration…")
        btn_load_cal.setStyleSheet(_btn_style())
        btn_load_cal.clicked.connect(self._load_calibration_file)
        cal_row.addWidget(btn_load_cal)
        sform.addRow(cal_row)

        self._sb_pbwo4 = QDoubleSpinBox()
        self._sb_pbwo4.setRange(0.0, 100.0)
        self._sb_pbwo4.setDecimals(6)
        self._sb_pbwo4.setValue(self._pbwo4)
        self._sb_pbwo4.valueChanged.connect(self._on_default_changed)
        sform.addRow(QLabel("PbWO4 default:"), self._sb_pbwo4)

        self._sb_pbglass = QDoubleSpinBox()
        self._sb_pbglass.setRange(0.0, 100.0)
        self._sb_pbglass.setDecimals(6)
        self._sb_pbglass.setValue(self._pbglass)
        self._sb_pbglass.valueChanged.connect(self._on_default_changed)
        sform.addRow(QLabel("PbGlass default:"), self._sb_pbglass)

        v.addWidget(src_grp)

        # ---- Color Range group ----
        # Reusable widget from hycal_geoview: min/max edits + Auto button.
        # auto_fit="minmax_nonzero" ignores zero-valued (masked) channels.
        # Click the Auto button for a one-shot fit; double-click to keep
        # auto-fitting after every edit.
        range_grp = QGroupBox("Color Range")
        rlayout = QHBoxLayout(range_grp)
        rlayout.setContentsMargins(8, 4, 8, 4)
        self._range_ctrl = ColorRangeControl(
            self._map,
            auto_fit="minmax_nonzero",
            orientation="horizontal",
        )
        rlayout.addWidget(self._range_ctrl)
        v.addWidget(range_grp)

        # ---- Mask group ----
        # Default interaction is "edit": clicking a module opens a popup
        # dialog to set its gain.  Toggling "Mask" or "Set" switches to a
        # drag-paint mode — Mask zeroes channels, Set assigns the value
        # in the line edit.  Reset clears all GUI edits; Undo reverts the
        # most recent action.
        mask_grp = QGroupBox("Mask")
        mlayout = QVBoxLayout(mask_grp)
        mlayout.setSpacing(6)

        # Build buttons + value edit; layout is split into two rows below.
        self._set_value_edit = QLineEdit("0.150000")
        self._set_value_edit.setMaximumWidth(110)
        self._set_value_edit.setValidator(
            QDoubleValidator(0.0, 100.0, 6, self._set_value_edit))
        self._set_value_edit.editingFinished.connect(self._on_set_value_changed)

        self._btn_set = QPushButton("Set")
        self._btn_set.setStyleSheet(
            _btn_style(checked_color=THEME.ACCENT_STRONG))
        self._btn_set.setCheckable(True)
        self._btn_set.setToolTip(
            "Toggle set mode — click or drag modules to apply the value")
        self._btn_set.toggled.connect(self._on_set_toggled)

        self._btn_set_all = QPushButton("Set All")
        self._btn_set_all.setStyleSheet(_btn_style())
        self._btn_set_all.setToolTip(
            "Apply the Set value to every DAQ-mapped channel")
        self._btn_set_all.clicked.connect(self._on_set_all)

        self._btn_mask = QPushButton("Mask")
        self._btn_mask.setStyleSheet(_btn_style(checked_color=THEME.DANGER))
        self._btn_mask.setCheckable(True)
        self._btn_mask.setToolTip(
            "Toggle mask mode — click or drag modules to close (gain = 0)")
        self._btn_mask.toggled.connect(self._on_mask_toggled)

        self._btn_mask_all = QPushButton("Mask All")
        self._btn_mask_all.setStyleSheet(_btn_style())
        self._btn_mask_all.setToolTip(
            "Set every DAQ-mapped channel's gain to 0")
        self._btn_mask_all.clicked.connect(self._on_mask_all)

        self._btn_undo = QPushButton("Undo")
        self._btn_undo.setStyleSheet(_btn_style())
        self._btn_undo.setToolTip("Revert the most recent edit")
        self._btn_undo.clicked.connect(self._undo)

        self._btn_reset = QPushButton("Reset")
        self._btn_reset.setStyleSheet(_btn_style())
        self._btn_reset.setToolTip(
            "Discard all manual edits, revert to loaded base")
        self._btn_reset.clicked.connect(self._reset_overrides)

        # Row 1: gain value edit + Set paint mode + Set All bulk apply.
        set_row = QHBoxLayout()
        set_row.addWidget(self._set_value_edit)
        set_row.addSpacing(6)
        set_row.addWidget(self._btn_set)
        set_row.addWidget(self._btn_set_all)
        set_row.addStretch()
        mlayout.addLayout(set_row)

        # Row 2: Mask paint + Mask All + Undo + Reset.
        mask_row = QHBoxLayout()
        mask_row.addWidget(self._btn_mask)
        mask_row.addWidget(self._btn_mask_all)
        mask_row.addStretch()
        mask_row.addWidget(self._btn_undo)
        mask_row.addWidget(self._btn_reset)
        mlayout.addLayout(mask_row)

        v.addWidget(mask_grp)

        # ---- Safe Cap group ----
        # One-shot clamp of every DAQ-mapped channel to [min, max].  Min is
        # floored at 0 (gains are non-negative), max defaults to 0.15.
        cap_grp = QGroupBox("Safe Cap")
        clayout = QHBoxLayout(cap_grp)
        clayout.setContentsMargins(8, 4, 8, 4)
        clayout.setSpacing(6)

        self._sb_cap_min = QDoubleSpinBox()
        self._sb_cap_min.setRange(0.0, 100.0)
        self._sb_cap_min.setDecimals(6)
        self._sb_cap_min.setValue(0.0)
        self._sb_cap_min.setToolTip(
            "Lower bound (clamped to >= 0; gains are non-negative)")

        self._sb_cap_max = QDoubleSpinBox()
        self._sb_cap_max.setRange(0.0, 100.0)
        self._sb_cap_max.setDecimals(6)
        self._sb_cap_max.setValue(0.15)
        self._sb_cap_max.setToolTip("Upper bound")

        self._btn_cap = QPushButton("Apply Safe Cap")
        self._btn_cap.setStyleSheet(_btn_style())
        self._btn_cap.setToolTip(
            "Clamp every DAQ-mapped channel's gain to [Min, Max]")
        self._btn_cap.clicked.connect(self._on_apply_safe_cap)

        clayout.addWidget(QLabel("Min:"))
        clayout.addWidget(self._sb_cap_min)
        clayout.addSpacing(6)
        clayout.addWidget(QLabel("Max:"))
        clayout.addWidget(self._sb_cap_max)
        clayout.addStretch()
        clayout.addWidget(self._btn_cap)

        v.addWidget(cap_grp)

        # ---- File row ----
        btn_load_cnf = QPushButton("Load .cnf…")
        btn_load_cnf.setStyleSheet(_btn_style())
        btn_load_cnf.clicked.connect(self._load_cnf_file)

        btn_save = QPushButton("Save .cnf…")
        btn_save.setStyleSheet(_btn_style())
        btn_save.clicked.connect(self._save_as)

        row = QHBoxLayout()
        row.addStretch()
        row.addWidget(btn_load_cnf)
        row.addWidget(btn_save)
        v.addLayout(row)

        # ---- Live .cnf preview ----
        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setStyleSheet(themed(
            f"QTextEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"font:9pt Monospace;border:1px solid {THEME.BORDER};}}"))
        v.addWidget(self._text, 1)

    # ---- handlers ----

    def _update_cal_label(self) -> None:
        if self._cal_path is not None:
            self._cal_label.setText(
                f"{self._cal_path.name}  ({len(self._cal)} entries)")
        elif self._cal:
            self._cal_label.setText(f"(in-memory, {len(self._cal)} entries)")
        else:
            self._cal_label.setText("(no calibration loaded)")

    def _on_mask_toggled(self, on: bool) -> None:
        if on and self._btn_set.isChecked():
            # Mask and Set are mutually exclusive — silently un-check Set.
            self._btn_set.blockSignals(True)
            self._btn_set.setChecked(False)
            self._btn_set.blockSignals(False)
        if on:
            self._map.set_paint_mode(PAINT_MODE_MASK)
            self._status.setText(
                "Mask mode — click or drag modules to close (gain = 0)")
        elif not self._btn_set.isChecked():
            self._map.set_paint_mode(PAINT_MODE_EDIT)
            self._status.setText(
                "Edit mode — click a module to set its gain")

    def _on_set_toggled(self, on: bool) -> None:
        if on:
            v = self._read_set_value()
            if v is None:
                self._status.setText(
                    "Set: enter a numeric gain value first")
                self._btn_set.blockSignals(True)
                self._btn_set.setChecked(False)
                self._btn_set.blockSignals(False)
                return
            if self._btn_mask.isChecked():
                self._btn_mask.blockSignals(True)
                self._btn_mask.setChecked(False)
                self._btn_mask.blockSignals(False)
            self._map.set_paint_value(v)
            self._map.set_paint_mode(PAINT_MODE_SET)
            self._status.setText(
                f"Set mode — click or drag modules to apply gain = {v:.6g}")
        elif not self._btn_mask.isChecked():
            self._map.set_paint_mode(PAINT_MODE_EDIT)
            self._status.setText(
                "Edit mode — click a module to set its gain")

    def _on_set_value_changed(self) -> None:
        v = self._read_set_value()
        if v is None:
            return
        self._map.set_paint_value(v)
        if self._btn_set.isChecked():
            self._status.setText(f"Set value = {v:.6g}")

    def _read_set_value(self) -> Optional[float]:
        try:
            return float(self._set_value_edit.text())
        except ValueError:
            return None

    def _on_mask_all(self) -> None:
        if QMessageBox.question(
                self, "Mask all channels?",
                "Set every DAQ-mapped channel's gain to 0?\n"
                "Use Undo to revert.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
        ) != QMessageBox.StandardButton.Yes:
            return
        self._bulk_apply(0.0, "Masked all")

    def _on_set_all(self) -> None:
        v = self._read_set_value()
        if v is None:
            self._status.setText("Set All: enter a numeric gain value first")
            return
        if QMessageBox.question(
                self, "Set all channels?",
                f"Set every DAQ-mapped channel's gain to {v:.6g}?\n"
                "Use Undo to revert.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
        ) != QMessageBox.StandardButton.Yes:
            return
        self._bulk_apply(v, f"Set all to {v:.6g}")

    def _bulk_apply(self, value: float, status_msg: str) -> None:
        """Apply ``value`` to every DAQ-mapped channel; record one undo batch."""
        gains = dict(self._map.gains)
        overrides = dict(self._map.overrides)
        batch: List[Tuple[str, Optional[float]]] = []
        for name, m in self._mod_map.items():
            if m.crate < 0:
                continue
            if (gains.get(name) == value and overrides.get(name) == value):
                continue
            batch.append((name, overrides.get(name)))
            gains[name] = value
            overrides[name] = value
        if not batch:
            self._status.setText("Nothing to change")
            return
        self._history.append(batch)
        self._map.set_gains(gains, overrides)
        self._notify_range_values()
        self._refresh_text()
        self._status.setText(f"{status_msg} — {len(batch)} channel(s)")

    def _on_apply_safe_cap(self) -> None:
        try:
            capped = safe_cap_gains(
                self._map.gains,
                self._sb_cap_min.value(),
                self._sb_cap_max.value())
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid range", str(exc))
            return

        min_cap = max(0.0, self._sb_cap_min.value())
        max_cap = self._sb_cap_max.value()
        if QMessageBox.question(
                self, "Apply safe cap?",
                f"Clamp every DAQ-mapped channel's gain to "
                f"[{min_cap:.6g}, {max_cap:.6g}]?\nUse Undo to revert.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
        ) != QMessageBox.StandardButton.Yes:
            return

        gains = dict(self._map.gains)
        overrides = dict(self._map.overrides)
        batch: List[Tuple[str, Optional[float]]] = []
        for name, m in self._mod_map.items():
            if m.crate < 0:
                continue
            new_v = capped.get(name)
            if new_v is None or new_v == gains.get(name):
                continue
            batch.append((name, overrides.get(name)))
            gains[name] = new_v
            overrides[name] = new_v
        if not batch:
            self._status.setText("Safe Cap: all channels already in range")
            return
        self._history.append(batch)
        self._map.set_gains(gains, overrides)
        self._notify_range_values()
        self._refresh_text()
        self._status.setText(
            f"Safe-capped to [{min_cap:.6g}, {max_cap:.6g}] — "
            f"{len(batch)} channel(s)")

    def _on_default_changed(self, _) -> None:
        self._pbwo4 = self._sb_pbwo4.value()
        self._pbglass = self._sb_pbglass.value()
        self._rebuild_from_base(keep_overrides=True)

    def _on_module_edit_requested(self, name: str, current: float) -> None:
        """Open a popup dialog to set ``name``'s gain.  Apply / Enter
        commits and closes; Cancel / Esc / X discards."""
        m = self._mod_map.get(name)
        crate_str = (f"  ({CRATE_NAMES[m.crate]} slot {m.slot} ch {m.channel})"
                     if m and m.crate >= 0 else "")
        new_val, ok = QInputDialog.getDouble(
            self, f"Edit gain — {name}",
            f"{name}{crate_str}\nGain (current: {current:.6g}):",
            current, 0.0, 100.0, 6)
        if not ok:
            return
        prior = self._map.overrides.get(name)
        self._history.append([(name, prior)])

        gains = dict(self._map.gains)
        overrides = dict(self._map.overrides)
        gains[name] = new_val
        overrides[name] = new_val
        self._map.set_gains(gains, overrides)

        self._notify_range_values()
        self._refresh_text()
        self._status.setText(f"Set {name} = {new_val:.6g}")

    def _on_paint_committed(self,
                            batch: List[Tuple[str, Optional[float]]]) -> None:
        """Drag-paint (mask) finished — record batch onto the undo stack."""
        if not batch:
            return
        self._history.append(list(batch))
        self._notify_range_values()
        self._refresh_text()
        self._status.setText(f"Masked {len(batch)} module(s)")

    def _undo(self) -> None:
        if not self._history:
            self._status.setText("Nothing to undo")
            return
        batch = self._history.pop()
        gains = dict(self._map.gains)
        overrides = dict(self._map.overrides)
        for name, prior in batch:
            if prior is None:
                overrides.pop(name, None)
                gains[name] = self._base_gains.get(name, 0.0)
            else:
                overrides[name] = prior
                gains[name] = prior
        self._map.set_gains(gains, overrides)
        self._notify_range_values()
        self._refresh_text()
        self._status.setText(f"Undone {len(batch)} edit(s)")

    def _on_hover(self, name: str) -> None:
        m = self._mod_map.get(name)
        if not m or m.crate < 0:
            self._status.setText(f"{name}  (no DAQ mapping)")
            return
        v = self._map.gains.get(name, 0.0)
        flags = []
        if name in self._map.overrides:
            flags.append("edit")
        if v == 0.0:
            flags.append("closed")
        tag = ("  [" + ", ".join(flags) + "]") if flags else ""
        self._status.setText(
            f"{name}  ({CRATE_NAMES[m.crate]} slot {m.slot} ch {m.channel}) "
            f" gain={v:.6g}{tag}")

    def _refresh_text(self) -> None:
        cal_line = (f"# calibration : {self._cal_path}"
                    if self._cal_path
                    else "# calibration : (none)")
        text, _ = render_cnf(
            self._daq, self._map.gains,
            [cal_line,
             f"# edits applied : {len(self._map.overrides)}",
             f"# defaults : PbWO4={self._pbwo4}, PbGlass={self._pbglass}"])
        self._text.setPlainText(text)

    def _reset_overrides(self) -> None:
        self._map.set_gains(self._base_gains, {})
        self._history.clear()
        self._notify_range_values()
        self._refresh_text()
        self._status.setText("All edits cleared")

    def _load_calibration_file(self) -> None:
        cal_dir = self._db_dir / "calibration"
        start = str(cal_dir if cal_dir.is_dir() else self._db_dir)
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Calibration JSON", start,
            "JSON (*.json);;All Files (*)")
        if not path:
            return
        try:
            cal = load_calibration(Path(path))
        except Exception as exc:
            self._status.setText(f"Calibration load failed: {exc}")
            return
        self._cal = cal
        self._cal_path = Path(path)
        self._update_cal_label()
        # Drop manual edits — the user just changed the base.
        self._history.clear()
        self._rebuild_from_base(keep_overrides=False)
        self._status.setText(
            f"Loaded calibration {self._cal_path.name}: {len(cal)} entries")

    def _save_as(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save FADC Gain Config",
            self._output_path or "adchycal_gain.cnf",
            "Config Files (*.cnf);;All Files (*)")
        if not path:
            return
        self._output_path = path
        with open(path, "w") as f:
            f.write(self._text.toPlainText())
        self._status.setText(f"Saved to {path}")

    def _load_cnf_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load FADC Gain Config", "",
            "Config Files (*.cnf);;All Files (*)")
        if not path:
            return
        with open(path) as f:
            text = f.read()
        loaded = parse_cnf_text(text, self._daq)
        overrides = self._diff_overrides(loaded)
        merged = dict(self._base_gains)
        merged.update(loaded)
        self._history.clear()
        self._map.set_gains(merged, overrides)
        self._notify_range_values()
        self._refresh_text()
        self._status.setText(
            f"Loaded {Path(path).name}: {len(loaded)} channels, "
            f"{len(overrides)} differ from base")


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interactive FADC gain config editor (always opens GUI)")
    parser.add_argument("-c", "--calibration",
                        help="Pre-load this calibration JSON in the editor "
                             "(equivalent to clicking Load Calibration…)")
    parser.add_argument("--pbwo4-gain", type=float, default=1.0,
                        help="Initial PbWO4 default gain (default: 1.0)")
    parser.add_argument("--pbglass-gain", type=float, default=1.0,
                        help="Initial PbGlass default gain (default: 1.0)")
    parser.add_argument("-d", "--database",
                        help="Database directory (default: auto-search)")
    parser.add_argument("-o", "--output", default="adchycal_gain.cnf",
                        help="Default save path (default: adchycal_gain.cnf)")
    parser.add_argument("-i", "--input",
                        help="Pre-load this .cnf as starting edits")
    parser.add_argument("--theme", default="dark",
                        choices=available_themes(),
                        help="GUI colour theme (default: dark)")
    args = parser.parse_args()

    db_dir = find_database_dir(args.database)
    print(f"database : {db_dir}")

    mod_types = load_modules(db_dir)
    daq = load_daq_map(db_dir)
    print(f"modules  : {len(mod_types)}   daq entries: {len(daq)}")

    cal: Dict[str, float] = {}
    cal_path: Optional[Path] = None
    if args.calibration:
        p = Path(args.calibration)
        if not p.is_absolute() and not p.is_file():
            alt = db_dir / "calibration" / p.name
            if alt.is_file():
                p = alt
        if not p.is_file():
            sys.exit(f"error: calibration file not found: {args.calibration}")
        cal = load_calibration(p)
        cal_path = p
        print(f"cal file : {cal_path}  ({len(cal)} entries)")

    initial_overrides: Dict[str, float] = {}
    if args.input:
        input_path = Path(args.input)
        if not input_path.is_file():
            sys.exit(f"error: --input file not found: {input_path}")
        with open(input_path) as f:
            initial_overrides = parse_cnf_text(f.read(), daq)
        print(f"input cnf: {input_path}  ({len(initial_overrides)} channels)")

    modules_info = _load_module_info(db_dir)

    set_theme(args.theme)
    app = QApplication.instance() or QApplication(sys.argv)
    win = _GainEditor(modules_info, daq, mod_types, db_dir,
                      cal=cal, cal_path=cal_path,
                      pbwo4_gain=args.pbwo4_gain,
                      pbglass_gain=args.pbglass_gain,
                      initial_overrides=initial_overrides,
                      output_path=args.output)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
