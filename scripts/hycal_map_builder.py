#!/usr/bin/env python3
"""
HyCal Map Builder (PyQt6)
=========================
Simple HyCal geometry viewer that colour-maps user data loaded from
JSON or plain-text files.

Data formats
------------
* JSON  : {"<module_name>": {"<field>": <value>, ...}, ...}

          Values may also be a list of history entries; the last entry
          of each list is used (so gain_equalization_results.json-style
          per-module history files work directly). Nested dicts are
          flattened with dot notation, e.g. fit.slope / edge.percentage.
          Non-numeric fields (strings like timestamps) are ignored.

* Text  : whitespace / comma / tab delimited rows

            <module_name> <val1> <val2> ...

          Lines starting with '#' are ignored. If the first non-comment
          row has a non-numeric second column it is treated as a header
          naming the columns; otherwise columns get default names
          (col1, col2, ...).

Usage
-----
    python scripts/hycal_map_builder.py                        # empty map
    python scripts/hycal_map_builder.py mydata.json            # auto-load
    python scripts/hycal_map_builder.py mydata.txt --field rms
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QComboBox, QSlider, QCheckBox,
    QFileDialog, QMessageBox,
)
from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QColor, QFont, QPen

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget, cmap_qcolor,
    PALETTE_NAMES, apply_theme_palette, set_theme,
    available_themes, THEME, ColorRangeControl,
)


# ===========================================================================
#  Paths
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent

def _find_modules_json() -> Path:
    candidates = [
        SCRIPT_DIR / ".." / "database" / "hycal_map.json",
        Path.cwd() / "database" / "hycal_map.json",
        Path.cwd() / "hycal_map.json",
    ]
    for c in candidates:
        if c.is_file():
            return c.resolve()
    return (SCRIPT_DIR / ".." / "database" / "hycal_map.json").resolve()

MODULES_JSON = _find_modules_json()


# ===========================================================================
#  Data loading
# ===========================================================================

def load_data_file(path: Path) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    """
    Returns (data, fields) where
      data[field][module_name] -> float
      fields is the ordered list of field names
    """
    text = path.read_text()
    try:
        obj = json.loads(text)
        if isinstance(obj, dict):
            return _data_from_json_dict(obj)
    except json.JSONDecodeError:
        pass
    return _data_from_text(text)


def _flatten_dict(d: Dict, prefix: str = "") -> Dict[str, object]:
    out: Dict[str, object] = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten_dict(v, key + "."))
        else:
            out[key] = v
    return out


def _data_from_json_dict(obj: Dict) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    # Normalise each module's value to a flat dict of scalars.
    #   {name: {key: val, ...}}              -> as-is
    #   {name: [{entry}, {entry}, ...]}      -> last entry (history)
    # Nested dicts are flattened with dot-joined keys.
    per_module: Dict[str, Dict[str, object]] = {}
    for name, entry in obj.items():
        if isinstance(entry, list):
            if not entry:
                continue
            entry = entry[-1]
        if not isinstance(entry, dict):
            continue
        flat = _flatten_dict(entry)
        if flat:
            per_module[str(name)] = flat

    fields: List[str] = []
    seen = set()
    for flat in per_module.values():
        for k in flat:
            if k not in seen:
                seen.add(k)
                fields.append(k)

    data: Dict[str, Dict[str, float]] = {f: {} for f in fields}
    for name, flat in per_module.items():
        for k, v in flat.items():
            try:
                data[k][name] = float(v)
            except (TypeError, ValueError):
                pass

    # Drop fields that ended up with no numeric values (e.g. timestamps).
    data = {k: v for k, v in data.items() if v}
    fields = [f for f in fields if f in data]
    return data, fields


_SPLIT = re.compile(r"[,\s\t]+")

def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False

def _data_from_text(text: str) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    rows: List[List[str]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p for p in _SPLIT.split(line) if p]
        if parts:
            rows.append(parts)
    if not rows:
        return {}, []

    first = rows[0]
    ncols = max(len(r) for r in rows) - 1
    if ncols <= 0:
        return {}, []

    header = not all(_is_number(c) for c in first[1:])
    if header:
        fields = [f or f"col{i+1}" for i, f in enumerate(first[1:])]
        data_rows = rows[1:]
    else:
        fields = [f"col{i+1}" for i in range(ncols)]
        data_rows = rows

    # pad/trim field list to ncols
    while len(fields) < ncols:
        fields.append(f"col{len(fields)+1}")
    fields = fields[:ncols]

    data: Dict[str, Dict[str, float]] = {f: {} for f in fields}
    for row in data_rows:
        if len(row) < 2:
            continue
        name = row[0]
        for i, field in enumerate(fields):
            idx = i + 1
            if idx >= len(row):
                break
            try:
                data[field][name] = float(row[idx])
            except ValueError:
                pass
    return data, fields


# ===========================================================================
#  HyCal map widget  (PbGlass alpha overlay + zoom/pan)
# ===========================================================================

def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "\u2014"
    if v == 0:
        return "0"
    return f"{v:.6g}"


class MapBuilderWidget(HyCalMapWidget):
    """HyCal map with adjustable PbGlass transparency and field-name colourbar.

    LMS and Veto modules are partitioned from the loaded module list and can
    be toggled independently via :meth:`set_lms_visible` /
    :meth:`set_scint_visible`. The canvas is always centred on geometric
    origin ``(0, 0)``, so toggling those modules does not shift the view.
    """

    _BOUNDS_TYPES = ("PbWO4", "PbGlass")   # modules that define the HyCal footprint

    def __init__(self, parent=None):
        super().__init__(parent, margin_top=8, enable_zoom_pan=True,
                         min_size=(500, 500))
        self._hycal: List[Module] = []
        self._lms: List[Module] = []
        self._scint: List[Module] = []
        self._show_lms = False
        self._show_scint = False
        self._pbglass_names: set = set()
        self._pbglass_alpha: float = 1.0
        self._field_label = ""

    def set_modules(self, modules: List[Module]):
        self._hycal = [m for m in modules if m.mod_type in self._BOUNDS_TYPES]
        self._lms   = [m for m in modules if m.mod_type == "LMS"]
        self._scint = [m for m in modules if m.mod_type == "Veto"]
        self._pbglass_names = {m.name for m in self._hycal
                               if m.mod_type == "PbGlass"}
        self._rebuild_layout()

    def set_lms_visible(self, on: bool):
        if bool(on) == self._show_lms:
            return
        self._show_lms = bool(on)
        self._rebuild_layout()

    def set_scint_visible(self, on: bool):
        if bool(on) == self._show_scint:
            return
        self._show_scint = bool(on)
        self._rebuild_layout()

    def _rebuild_layout(self):
        display = list(self._hycal)
        if self._show_lms:
            display.extend(self._lms)
        if self._show_scint:
            display.extend(self._scint)
        self._modules = display
        # Symmetric bounds around origin so the canvas stays centred on
        # (0, 0) even when LMS/Veto sit off to the side.
        if display:
            rx = max(max(abs(m.x - m.sx / 2), abs(m.x + m.sx / 2))
                     for m in display)
            ry = max(max(abs(m.y - m.sy / 2), abs(m.y + m.sy / 2))
                     for m in display)
            self._geo_bounds = (-rx, rx, -ry, ry)
        self._layout_dirty = True
        self.update()

    def set_values(self, values: Dict[str, float], label: str = ""):
        self._field_label = label
        super().set_values(values)

    def set_pbglass_alpha(self, a: float):
        self._pbglass_alpha = max(0.0, min(1.0, a))
        self.update()

    def _fmt_value(self, v: float) -> str:
        return _fmt(v)

    def _colorbar_center_text(self) -> str:
        mid = PALETTE_NAMES[self._palette_idx]
        if self._log_scale:
            mid += "  [log]"
        if self._field_label:
            mid = self._field_label + "  \u2014  " + mid
        return mid

    def _paint_empty(self, p, w, h):
        p.setPen(self.EMPTY_TEXT)
        p.setFont(QFont("Consolas", 12))
        p.drawText(QRectF(0, 0, w, h),
                   Qt.AlignmentFlag.AlignCenter, "No modules loaded")

    def _paint_modules(self, p):
        stops = self.palette_stops()
        vmin, vmax = self._vmin, self._vmax
        no_data = self.NO_DATA_COLOR
        log_scale = self._log_scale
        if log_scale:
            log_lo = math.log10(max(vmin, 1e-9))
            log_hi = math.log10(max(vmax, vmin * 10, 1e-8))

        glass_alpha = self._pbglass_alpha
        # Module frame — TEXT_DIM is defined to be readable on the active
        # theme's canvas, so it's visible in both dark and light modes.
        frame_col_base = QColor(THEME.TEXT_DIM)
        pbglass_names = self._pbglass_names

        for name, rect in self._rects.items():
            is_glass = name in pbglass_names
            a = glass_alpha if is_glass else 1.0

            v = self._values.get(name)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                fill = QColor(no_data)
            else:
                if log_scale:
                    lv = math.log10(max(v, 1e-9))
                    t = (lv - log_lo) / (log_hi - log_lo) if log_hi > log_lo else 0.5
                else:
                    t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                fill = cmap_qcolor(t, stops)

            if a < 1.0:
                fill = QColor(fill)
                fill.setAlphaF(a)
            if a > 0.0:
                p.fillRect(rect, fill)
                frame = QColor(frame_col_base)
                frame.setAlphaF(a * 0.8)
                p.setPen(QPen(frame, 1.0))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(rect)


# ===========================================================================
#  Main window
# ===========================================================================

class MapBuilderWindow(QMainWindow):

    def __init__(self, modules: List[Module],
                 data_file: Optional[Path] = None,
                 initial_field: Optional[str] = None):
        super().__init__()
        self._modules = modules
        self._data: Dict[str, Dict[str, float]] = {}
        self._fields: List[str] = []
        self._current_field: Optional[str] = None
        self._data_path: Optional[Path] = None

        self._build_ui()
        self._map.set_modules(modules)

        if data_file is not None:
            self._load_file(data_file, preferred_field=initial_field)

    # -- ui --

    def _build_ui(self):
        self.setWindowTitle("HyCal Map Builder")
        self.resize(900, 960)
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar: title + file --
        top = QHBoxLayout()
        title = QLabel("HYCAL MAP BUILDER")
        title.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        title.setStyleSheet(f"color:{THEME.ACCENT};")
        top.addWidget(title)

        top.addStretch()

        self._file_lbl = QLabel("(no file loaded)")
        self._file_lbl.setFont(QFont("Monospace", 10))
        self._file_lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};")
        top.addWidget(self._file_lbl)

        top.addWidget(self._make_btn("Open File...", THEME.ACCENT,
                                     self._open_file))
        root.addLayout(top)

        # -- map --
        self._map = MapBuilderWidget()
        self._map.paletteClicked.connect(self._cycle_palette)
        self._map.moduleHovered.connect(self._on_hover)
        root.addWidget(self._map, stretch=1)

        # -- controls: field + palette + auto + range + log --
        ctrl = QHBoxLayout()
        ctrl.addWidget(self._styled_label("Field:"))
        self._field_box = QComboBox()
        self._field_box.setMinimumWidth(160)
        self._field_box.setFont(QFont("Monospace", 11))
        self._field_box.setStyleSheet(
            f"QComboBox{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};border-radius:8px;"
            f"padding:2px 6px;}}"
            f"QComboBox QAbstractItemView{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT};}}")
        self._field_box.currentTextChanged.connect(self._on_field_changed)
        ctrl.addWidget(self._field_box)

        ctrl.addSpacing(12)
        ctrl.addWidget(self._make_btn("Palette \u25B6", THEME.TEXT,
                                      self._cycle_palette))

        ctrl.addSpacing(12)
        # Reusable control: min/max edits, Auto button (single-click=oneshot,
        # double-click=pin), Log toggle.  auto_fit="percentile" with (2, 98)
        # matches the long-tail distributions in calibration histograms.
        # Starts pinned so switching fields auto-rescales until the user opts out.
        self._range_ctrl = ColorRangeControl(
            self._map,
            auto_fit="percentile",
            auto_fit_percentile=(2.0, 98.0),
            include_log=True,
            start_pinned=True,
        )
        ctrl.addWidget(self._range_ctrl)

        ctrl.addSpacing(12)
        ctrl.addWidget(self._styled_label("PbGlass \u03B1:"))
        self._alpha_slider = QSlider(Qt.Orientation.Horizontal)
        self._alpha_slider.setRange(0, 100)
        self._alpha_slider.setValue(100)
        self._alpha_slider.setFixedWidth(120)
        self._alpha_slider.setStyleSheet(
            f"QSlider::groove:horizontal{{background:{THEME.BORDER};height:4px;border-radius:2px;}}"
            f"QSlider::sub-page:horizontal{{background:{THEME.ACCENT};height:4px;border-radius:2px;}}"
            f"QSlider::handle:horizontal{{background:{THEME.ACCENT};width:12px;"
            f"margin:-5px 0;border-radius:6px;}}"
            f"QSlider::handle:horizontal:hover{{background:{THEME.ACCENT_STRONG};}}")
        self._alpha_slider.valueChanged.connect(self._on_alpha_changed)
        ctrl.addWidget(self._alpha_slider)
        self._alpha_lbl = QLabel("100%")
        self._alpha_lbl.setFont(QFont("Monospace", 10))
        self._alpha_lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};")
        self._alpha_lbl.setFixedWidth(40)
        ctrl.addWidget(self._alpha_lbl)

        ctrl.addSpacing(12)
        self._lms_chk = self._styled_checkbox("LMS", self._on_lms_toggled)
        self._scint_chk = self._styled_checkbox("Veto", self._on_scint_toggled)
        ctrl.addWidget(self._lms_chk)
        ctrl.addWidget(self._scint_chk)

        ctrl.addStretch()
        root.addLayout(ctrl)

        # -- info / stats --
        info_row = QHBoxLayout()
        info_row.setSpacing(6)
        self._info = QLabel("Hover over a module")
        self._info.setFont(QFont("Monospace", 11))
        self._info.setStyleSheet(
            f"QLabel{{background:{THEME.PANEL};color:{THEME.TEXT};padding:4px 8px;"
            f"border:1px solid {THEME.BORDER};border-radius:8px;}}")
        self._info.setFixedHeight(28)
        info_row.addWidget(self._info, stretch=1)

        self._stats_lbl = QLabel("")
        self._stats_lbl.setFont(QFont("Monospace", 11))
        self._stats_lbl.setStyleSheet(
            f"QLabel{{background:{THEME.PANEL};color:{THEME.TEXT_DIM};padding:4px 8px;"
            f"border:1px solid {THEME.BORDER};border-radius:8px;}}")
        self._stats_lbl.setFixedHeight(28)
        info_row.addWidget(self._stats_lbl)
        root.addLayout(info_row)

    def _make_btn(self, text: str, fg: str, slot) -> QPushButton:
        btn = QPushButton(text)
        f = QFont("Monospace", 11); f.setBold(True)
        btn.setFont(f)
        btn.setStyleSheet(
            f"QPushButton{{background:{THEME.BUTTON};color:{fg};"
            f"border:1px solid {THEME.BORDER};padding:5px 14px;"
            f"border-radius:8px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")
        btn.clicked.connect(slot)
        return btn

    def _styled_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 11))
        lbl.setStyleSheet(f"color:{THEME.TEXT};")
        return lbl

    def _styled_checkbox(self, text: str, slot) -> QCheckBox:
        cb = QCheckBox(text)
        cb.setChecked(False)
        cb.setFont(QFont("Monospace", 11))
        cb.setStyleSheet(
            f"QCheckBox{{color:{THEME.TEXT};spacing:6px;}}"
            f"QCheckBox::indicator{{width:13px;height:13px;border-radius:3px;"
            f"border:1px solid {THEME.BORDER};background:{THEME.PANEL};}}"
            f"QCheckBox::indicator:hover{{border:1px solid {THEME.ACCENT};}}"
            f"QCheckBox::indicator:checked{{background:{THEME.ACCENT};"
            f"border:1px solid {THEME.ACCENT};}}")
        cb.toggled.connect(slot)
        return cb

    # -- actions --

    def _open_file(self):
        start_dir = str(self._data_path.parent) if self._data_path else str(Path.cwd())
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open data file", start_dir,
            "Data files (*.json *.txt *.dat *.csv *.tsv);;All files (*)")
        if path_str:
            self._load_file(Path(path_str))

    def _load_file(self, path: Path, preferred_field: Optional[str] = None):
        try:
            data, fields = load_data_file(path)
        except Exception as ex:
            QMessageBox.warning(self, "Load failed", f"{path}\n\n{ex}")
            return
        if not fields:
            QMessageBox.warning(self, "No data", f"No usable fields in {path}")
            return

        self._data = data
        self._fields = fields
        self._data_path = path
        self._file_lbl.setText(path.name)
        self._file_lbl.setToolTip(str(path))

        self._field_box.blockSignals(True)
        self._field_box.clear()
        self._field_box.addItems(fields)
        if preferred_field and preferred_field in fields:
            self._field_box.setCurrentText(preferred_field)
        self._field_box.blockSignals(False)
        self._current_field = self._field_box.currentText()
        self._refresh_values()

    def _on_field_changed(self, field: str):
        if not field:
            return
        self._current_field = field
        self._refresh_values()

    def _refresh_values(self):
        if not self._current_field:
            return
        values = self._data.get(self._current_field, {})
        self._map.set_values(values, label=self._current_field)
        # Pin handles re-fit when on; otherwise the user's manual range stays.
        self._range_ctrl.notify_values_changed(values)
        self._update_stats()

    def _update_stats(self):
        if not self._current_field:
            self._stats_lbl.setText("")
            return
        vals = [v for v in self._data.get(self._current_field, {}).values()
                if v is not None and not (isinstance(v, float) and math.isnan(v))]
        if not vals:
            self._stats_lbl.setText("no data")
            return
        arr = np.asarray(vals, dtype=float)
        self._stats_lbl.setText(
            f"N={arr.size}  mean={_fmt(float(arr.mean()))}"
            f"  rms={_fmt(float(arr.std()))}")

    def _cycle_palette(self):
        self._map.set_palette(self._map.palette_idx() + 1)

    def _on_alpha_changed(self, v: int):
        self._map.set_pbglass_alpha(v / 100.0)
        self._alpha_lbl.setText(f"{v}%")

    def _on_lms_toggled(self, on: bool):
        self._map.set_lms_visible(on)

    def _on_scint_toggled(self, on: bool):
        self._map.set_scint_visible(on)

    def _on_hover(self, name: str):
        parts = [name]
        for m in self._modules:
            if m.name == name:
                parts.append(f"({m.mod_type})")
                break
        if self._current_field:
            v = self._data.get(self._current_field, {}).get(name)
            if v is not None:
                parts.append(f"{self._current_field} = {_fmt(v)}")
        self._info.setText("    ".join(parts))


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="HyCal geo-view map builder")
    ap.add_argument("data_file", nargs="?", type=Path,
                    help="JSON or plain-text data file to load at startup")
    ap.add_argument("--field", type=str, default=None,
                    help="Initial field to display (default: first field)")
    ap.add_argument("--modules", type=Path, default=MODULES_JSON,
                    help=f"Path to hycal_map.json (default: {MODULES_JSON})")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

    if not args.modules.is_file():
        print(f"ERROR: hycal modules file not found: {args.modules}",
              file=sys.stderr)
        sys.exit(1)

    modules = load_modules(args.modules)
    print(f"Loaded {len(modules)} HyCal modules from {args.modules}")

    app = QApplication(sys.argv)
    win = MapBuilderWindow(modules,
                           data_file=args.data_file,
                           initial_field=args.field)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
