#!/usr/bin/env python3
"""
HyCal Trigger Mask Editor (PyQt6)
==================================
Visual editor for FAV3 trigger masks.  Displays a HyCal geo view with
LMS / V modules below.  Click or drag to toggle channels off/on.
Generates trigger mask config text (only disabled channels are written).

Usage
-----
    python trigger_mask_editor.py
    python trigger_mask_editor.py -o output.cnf          # auto-save path
    python trigger_mask_editor.py -i existing.cnf         # load existing mask
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# This tool lives in scripts/daq_tool/; import hycal_geoview from the
# parent scripts/ directory at runtime.
_SCRIPTS_DIR = Path(__file__).resolve().parent.parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QSplitter, QFileDialog,
    QGroupBox, QFormLayout,
)
from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QColor, QPen, QFont

from hycal_geoview import (
    HyCalMapWidget as _HyCalMapBase,
    THEME, apply_theme_palette, set_theme, available_themes, themed,
)


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / ".." / "database"
HYCAL_MAP_JSON = DB_DIR / "hycal_map.json"

NUM_CRATES = 7
CRATE_NAMES = [f"adchycal{i}" for i in range(1, NUM_CRATES + 1)]
CHANNELS_PER_SLOT = 16

# LMS / V module display positions below HyCal
_BOTTOM_Y = -640.0
_BOTTOM_SZ = 50.0
_LMS_V_XPOS = {
    "LMS1": -200.0, "LMS2": -145.0, "LMS3": -90.0,
    "V1": 35.0, "V2": 90.0, "V3": 145.0, "V4": 200.0,
}
_LABEL_NAMES = set(_LMS_V_XPOS.keys())

# Colours — resolved against the active :class:`THEME` at paint time.
def _col_on()       -> QColor: return QColor(THEME.SUCCESS)        # enabled
def _col_off()      -> QColor: return QColor(THEME.DANGER)         # disabled
def _col_on_glass() -> QColor: return QColor(THEME.SUCCESS).darker(140)
def _col_hover()    -> QColor: return QColor(THEME.ACCENT)
def _col_no_daq()   -> QColor: return QColor(THEME.NO_DATA)
def _col_text()     -> QColor: return QColor(THEME.TEXT)


def _btn_style() -> str:
    return themed(
        f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT};"
        f"border:1px solid {THEME.BORDER};padding:6px 14px;"
        f"font:10pt;border-radius:8px;}}"
        f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")


# ===========================================================================
#  Data loading
# ===========================================================================

class ModuleInfo:
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


def load_data() -> List[ModuleInfo]:
    with open(HYCAL_MAP_JSON) as f:
        entries = json.load(f)

    modules: List[ModuleInfo] = []
    for m in entries:
        name = m["n"]
        d = m.get("daq") or {}
        g = m.get("geo") or {}
        # Reposition LMS below HyCal for the editor display
        if name in _LMS_V_XPOS:
            x, y, sx, sy = _LMS_V_XPOS[name], _BOTTOM_Y, _BOTTOM_SZ, _BOTTOM_SZ
        else:
            x, y, sx, sy = g.get("x", 0.0), g.get("y", 0.0), g.get("sx", 0.0), g.get("sy", 0.0)
        crate = d.get("crate", -1)
        slot  = d.get("slot",  -1)
        ch    = d.get("channel", -1)
        modules.append(ModuleInfo(name, m["t"], x, y, sx, sy, crate, slot, ch))

    return modules


# ===========================================================================
#  Trigger mask I/O
# ===========================================================================

def generate_trigger_mask(modules: List[ModuleInfo],
                          disabled: Set[str]) -> str:
    """Generate trigger mask text.  Only crate/slot combos with at least
    one disabled channel are written.  Unmapped channels (no module in the
    DAQ map for that slot position) are always masked off."""
    # Group by (crate, slot) -> {channel: module_name}
    crate_slots: Dict[Tuple[int, int], Dict[int, str]] = {}
    for m in modules:
        if m.crate < 0:
            continue
        crate_slots.setdefault((m.crate, m.slot), {})[m.channel] = m.name

    # Build output per crate
    lines: List[str] = []
    for ci in range(NUM_CRATES):
        slot_lines: List[str] = []
        for slot in sorted(set(s for (c, s) in crate_slots if c == ci)):
            ch_map = crate_slots.get((ci, slot), {})
            mask = []
            off_names: List[str] = []
            has_disabled = False
            for ch in range(CHANNELS_PER_SLOT):
                mod_name = ch_map.get(ch)
                if mod_name is None:
                    mask.append("0")
                    off_names.append(f"ch{ch}:unmapped")
                    has_disabled = True
                elif mod_name in disabled:
                    mask.append("0")
                    off_names.append(mod_name)
                    has_disabled = True
                else:
                    mask.append("1")
            if has_disabled:
                slot_lines.append(f"# off: {', '.join(off_names)}")
                slot_lines.append(f"FAV3_SLOT {slot}")
                slot_lines.append(f"FAV3_TRG_MASK {' '.join(mask)}")

        if slot_lines:
            lines.append(f"FAV3_CRATE {CRATE_NAMES[ci]}")
            lines.extend(slot_lines)
            lines.append("FAV3_CRATE end")
            lines.append("")

    return "\n".join(lines)


def parse_trigger_mask(text: str, modules: List[ModuleInfo]) -> Set[str]:
    """Parse existing trigger mask text, return set of disabled module names."""
    # Build reverse lookup: (crate_name, slot, channel) -> module_name
    name_to_ci = {n: i for i, n in enumerate(CRATE_NAMES)}
    daq_lookup: Dict[Tuple[int, int, int], str] = {}
    for m in modules:
        if m.crate >= 0:
            daq_lookup[(m.crate, m.slot, m.channel)] = m.name

    disabled: Set[str] = set()
    current_crate_idx = -1

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("FAV3_CRATE"):
            cname = line.split()[1] if len(line.split()) > 1 else ""
            if cname == "end":
                current_crate_idx = -1
            else:
                current_crate_idx = name_to_ci.get(cname, -1)
        elif line.startswith("FAV3_SLOT"):
            current_slot = int(line.split()[1])
        elif line.startswith("FAV3_TRG_MASK") and current_crate_idx >= 0:
            parts = line.split()[1:]
            for ch, val in enumerate(parts):
                if val == "0":
                    mod = daq_lookup.get((current_crate_idx, current_slot, ch))
                    if mod:
                        disabled.add(mod)
    return disabled


# ===========================================================================
#  HyCal Map Widget
# ===========================================================================

class HyCalMapWidget(_HyCalMapBase):
    """Interactive HyCal geo view for trigger mask editing."""

    module_hovered = pyqtSignal(str)   # kept snake_case for script-level callers
    mask_changed = pyqtSignal()

    def __init__(self, modules: List[ModuleInfo], parent=None):
        super().__init__(parent, shrink=0.92, margin_top=10,
                         margin_bottom=32, include_lms=True,
                         show_colorbar=False, min_size=(500, 500))
        self._mod_map: Dict[str, ModuleInfo] = {m.name: m for m in modules}
        self._disabled: Set[str] = set()

        # drag-paint state (distinct from the base's pan-drag state)
        self._paint_dragging = False
        self._paint_mode: Optional[bool] = None  # True = disabling, False = enabling
        self._drag_visited: Set[str] = set()

        self.set_modules(modules)

    @property
    def disabled(self) -> Set[str]:
        return self._disabled

    @disabled.setter
    def disabled(self, s: Set[str]):
        self._disabled = set(s)
        self.update()

    # -- painting: fill-by-state + LMS/V labels + stats line --

    def _paint_modules(self, p):
        disabled = self._disabled
        mod_map = self._mod_map
        col_on       = _col_on()
        col_on_glass = _col_on_glass()
        col_off      = _col_off()
        col_no_daq   = _col_no_daq()
        for name, rect in self._rects.items():
            m = mod_map.get(name)
            if m and m.crate < 0:
                p.fillRect(rect, col_no_daq)
            elif name in disabled:
                p.fillRect(rect, col_off)
            elif m and m.mod_type == "PbGlass":
                p.fillRect(rect, col_on_glass)
            else:
                p.fillRect(rect, col_on)

    def _paint_overlays(self, p, w, h):
        # LMS / V labels
        p.setPen(_col_text())
        p.setFont(QFont("Monospace", 7, QFont.Weight.Bold))
        for name in _LABEL_NAMES:
            r = self._rects.get(name)
            if r is not None:
                p.drawText(r, Qt.AlignmentFlag.AlignCenter, name)
        # hover highlight (use mask-editor's hover colour)
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(_col_hover(), 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

    def _paint_after_colorbar(self, p, w, h):
        p.setFont(QFont("Monospace", 9))
        sw = 10
        y = h - 22
        swatches = (
            (_col_on(),       "enabled"),
            (_col_on_glass(), "glass"),
            (_col_off(),      "disabled"),
            (_col_no_daq(),   "no DAQ"),
        )

        x = 10
        fm = p.fontMetrics()
        for col, label in swatches:
            p.fillRect(QRectF(x, y, sw, sw), col)
            p.setPen(QColor(THEME.BORDER))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(QRectF(x, y, sw, sw))
            p.setPen(QColor(THEME.TEXT_DIM))
            text_w = fm.horizontalAdvance(label)
            p.drawText(QRectF(x + sw + 4, y - 2, text_w + 4, sw + 6),
                       Qt.AlignmentFlag.AlignLeft
                       | Qt.AlignmentFlag.AlignVCenter,
                       label)
            x += sw + 4 + text_w + 14

        total = sum(1 for m in self._mod_map.values() if m.crate >= 0)
        n_off = len(self._disabled)
        n_on = total - n_off
        p.setPen(QColor(THEME.TEXT))
        p.drawText(QRectF(0, y - 2, w - 10, sw + 6),
                   Qt.AlignmentFlag.AlignRight
                   | Qt.AlignmentFlag.AlignVCenter,
                   f"{n_on} on  /  {n_off} off    ({total} total)")

    # -- hit test: only DAQ-mapped modules are hittable --

    def _hit(self, pos) -> Optional[str]:
        for name, rect in self._rects.items():
            if rect.contains(pos):
                m = self._mod_map.get(name)
                if m and m.crate >= 0:
                    return name
        return None

    # -- mouse: drag-paint replaces base zoom/pan --

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        found = self._hit(event.position())
        if found:
            self._paint_dragging = True
            self._paint_mode = found not in self._disabled  # True = disable
            self._drag_visited = {found}
            if self._paint_mode:
                self._disabled.add(found)
            else:
                self._disabled.discard(found)
            self.update()
            self.mask_changed.emit()

    def mouseMoveEvent(self, event):
        pos = event.position()
        found = self._hit(pos)
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                self.module_hovered.emit(found)
        if self._paint_dragging and found and found not in self._drag_visited:
            self._drag_visited.add(found)
            if self._paint_mode:
                self._disabled.add(found)
            else:
                self._disabled.discard(found)
            self.update()
            self.mask_changed.emit()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._paint_dragging = False
            self._paint_mode = None
            self._drag_visited.clear()

    def wheelEvent(self, event):
        # wheel events should not do anything here (no zoom)
        event.ignore()


# ===========================================================================
#  Main Window
# ===========================================================================

class TriggerMaskEditor(QMainWindow):
    def __init__(self, modules: List[ModuleInfo],
                 initial_disabled: Set[str] = None,
                 output_path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("HyCal Trigger Mask Editor")
        self.resize(1600, 900)
        self._modules = modules
        self._output_path = output_path
        self._mod_map: Dict[str, ModuleInfo] = {m.name: m for m in modules}
        self._total_mapped = sum(1 for m in modules if m.crate >= 0)

        apply_theme_palette(self)
        # Window-scoped stylesheet so QLabel / QGroupBox text renders against
        # the dark surfaces correctly on Windows native style (where palette
        # propagation through QGroupBox / QSplitter is unreliable).
        self.setStyleSheet(themed(
            f"QLabel{{color:{THEME.TEXT};background:transparent;}}"
            f"QGroupBox{{color:{THEME.TEXT};background:transparent;"
            f"border:1px solid {THEME.BORDER};border-radius:6px;"
            f"margin-top:10px;padding-top:10px;}}"
            f"QGroupBox::title{{subcontrol-origin:margin;left:10px;"
            f"padding:0 6px;color:{THEME.TEXT_DIM};}}"
            f"QSplitter::handle{{background:{THEME.BORDER};}}"
        ))

        self._map = HyCalMapWidget(modules)
        if initial_disabled:
            self._map.disabled = initial_disabled

        self._status = QLabel("Click or drag modules to toggle trigger mask")
        self._status.setStyleSheet(themed(
            f"color:{THEME.TEXT_DIM};font:10pt Monospace;padding:6px;"))

        self._right = self._build_right_panel()

        # Signals
        self._map.module_hovered.connect(self._on_hover)
        self._map.mask_changed.connect(self._on_mask_changed)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._map)
        splitter.addWidget(self._right)
        splitter.setSizes([900, 700])

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter, 1)
        layout.addWidget(self._status)
        self.setCentralWidget(central)

        self._refresh_text()
        self._refresh_stats()

    def _build_right_panel(self) -> QWidget:
        right = QWidget()
        v = QVBoxLayout(right)
        v.setContentsMargins(8, 8, 8, 8)
        v.setSpacing(8)

        # ---- Status group ----
        stat_grp = QGroupBox("Status")
        sf = QFormLayout(stat_grp)
        sf.setHorizontalSpacing(12)
        sf.setVerticalSpacing(4)
        self._stat_total = QLabel("—")
        self._stat_enabled = QLabel("—")
        self._stat_disabled = QLabel("—")
        for lbl in (self._stat_total, self._stat_enabled, self._stat_disabled):
            lbl.setStyleSheet(themed(
                f"color:{THEME.TEXT};font:10pt Monospace;"))
            lbl.setTextFormat(Qt.TextFormat.RichText)
        sf.addRow("Total:",    self._stat_total)
        sf.addRow("Enabled:",  self._stat_enabled)
        sf.addRow("Disabled:", self._stat_disabled)
        v.addWidget(stat_grp)

        # ---- Bulk actions ----
        bulk_grp = QGroupBox("Bulk")
        bv = QHBoxLayout(bulk_grp)
        btn_enable_all = QPushButton("Enable All")
        btn_enable_all.setStyleSheet(_btn_style())
        btn_enable_all.clicked.connect(self._enable_all)
        btn_disable_all = QPushButton("Disable All")
        btn_disable_all.setStyleSheet(_btn_style())
        btn_disable_all.clicked.connect(self._disable_all)
        bv.addWidget(btn_enable_all)
        bv.addWidget(btn_disable_all)
        v.addWidget(bulk_grp)

        # ---- File ops ----
        file_grp = QGroupBox("File")
        fv = QHBoxLayout(file_grp)
        btn_load = QPushButton("Load…")
        btn_load.setStyleSheet(_btn_style())
        btn_load.clicked.connect(self._load_file)
        btn_save = QPushButton("Save As…")
        btn_save.setStyleSheet(_btn_style())
        btn_save.clicked.connect(self._save_as)
        fv.addWidget(btn_load)
        fv.addWidget(btn_save)
        v.addWidget(file_grp)

        # ---- Preview ----
        preview_lbl = QLabel("Preview")
        preview_lbl.setStyleSheet(themed(
            f"color:{THEME.TEXT_DIM};font:9pt;padding-top:4px;"))
        v.addWidget(preview_lbl)

        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setStyleSheet(themed(
            f"QTextEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"font:9pt Monospace;border:1px solid {THEME.BORDER};"
            f"border-radius:6px;}}"))
        v.addWidget(self._text, 1)

        return right

    def _on_mask_changed(self) -> None:
        self._refresh_text()
        self._refresh_stats()

    def _refresh_stats(self) -> None:
        n_off = len(self._map.disabled)
        n_on = self._total_mapped - n_off
        self._stat_total.setText(f"{self._total_mapped}")
        self._stat_enabled.setText(
            f'<span style="color:{THEME.SUCCESS}">{n_on}</span>')
        self._stat_disabled.setText(
            f'<span style="color:{THEME.DANGER}">{n_off}</span>')

    def _on_hover(self, name: str):
        m = self._mod_map.get(name)
        if m and m.crate >= 0:
            state = "OFF" if name in self._map.disabled else "ON"
            self._status.setText(
                f"{name}  ({CRATE_NAMES[m.crate]} slot {m.slot} ch {m.channel})  [{state}]")

    def _refresh_text(self):
        text = generate_trigger_mask(self._modules, self._map.disabled)
        if text:
            self._text.setPlainText(text)
        else:
            self._text.setPlainText("# All channels enabled — no mask needed")

    def _enable_all(self):
        self._map.disabled = set()
        self._map.update()
        self._on_mask_changed()

    def _disable_all(self):
        self._map.disabled = {m.name for m in self._modules if m.crate >= 0}
        self._map.update()
        self._on_mask_changed()

    def _save_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Trigger Mask", self._output_path or "trigger_mask.cnf",
            "Config Files (*.cnf);;All Files (*)")
        if path:
            self._output_path = path
            with open(path, "w") as f:
                f.write(self._text.toPlainText())
            self._status.setText(f"Saved to {path}")

    def _load_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Trigger Mask", "",
            "Config Files (*.cnf);;All Files (*)")
        if path:
            with open(path) as f:
                text = f.read()
            disabled = parse_trigger_mask(text, self._modules)
            self._map.disabled = disabled
            self._map.update()
            self._on_mask_changed()
            self._status.setText(f"Loaded {path} — {len(disabled)} channels disabled")


# ===========================================================================
#  Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(description="HyCal Trigger Mask Editor")
    parser.add_argument("-i", "--input", help="Load existing trigger mask file")
    parser.add_argument("-o", "--output", help="Default save path")
    parser.add_argument("--theme", choices=available_themes(), default="dark",
                        help="Colour theme (default: dark)")
    args = parser.parse_args()

    set_theme(args.theme)

    modules = load_data()

    initial_disabled: Set[str] = set()
    if args.input:
        with open(args.input) as f:
            initial_disabled = parse_trigger_mask(f.read(), modules)

    app = QApplication(sys.argv)
    win = TriggerMaskEditor(modules, initial_disabled, args.output)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
