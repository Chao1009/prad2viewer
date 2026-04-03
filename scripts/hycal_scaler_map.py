#!/usr/bin/env python3
"""
HyCal FADC Scaler Map (PyQt6)
=============================
Polls EPICS scaler channels (B_DET_HYCAL_FADC_<name>) for every HyCal
module and displays a live colour-coded geo map.

Usage
-----
    python scripts/hycal_scaler_map.py              # real EPICS
    python scripts/hycal_scaler_map.py --sim         # simulation (random)
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QSizePolicy,
)
from PyQt6.QtCore import Qt, QRectF, QTimer, pyqtSignal
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QLinearGradient, QPalette,
)


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_modules.json"

SCALER_PREFIX = "B_DET_HYCAL_FADC_"
POLL_INTERVAL_MS = 10_000   # 10 seconds


# ===========================================================================
#  Module database
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


def load_modules(path: Path) -> List[Module]:
    with open(path) as f:
        data = json.load(f)
    return [Module(e["n"], e["t"], e["x"], e["y"], e["sx"], e["sy"])
            for e in data]


# ===========================================================================
#  EPICS interfaces
# ===========================================================================

class RealScalerEPICS:
    """Read scaler PVs via pyepics."""

    def __init__(self, modules: List[Module]):
        import epics as _epics
        self._pvs: Dict[str, object] = {}
        for m in modules:
            if m.mod_type in ("PbWO4", "PbGlass"):
                pv = _epics.PV(SCALER_PREFIX + m.name, connection_timeout=2.0)
                self._pvs[m.name] = pv

    def get(self, name: str) -> Optional[float]:
        pv = self._pvs.get(name)
        if pv and pv.connected:
            return pv.get()
        return None

    def connection_count(self) -> Tuple[int, int]:
        n = sum(1 for pv in self._pvs.values() if pv.connected)
        return n, len(self._pvs)


class SimulatedScalerEPICS:
    """Return random values for testing."""

    def __init__(self, modules: List[Module]):
        self._rng = random.Random(0)
        self._names = [m.name for m in modules
                       if m.mod_type in ("PbWO4", "PbGlass")]

    def get(self, name: str) -> Optional[float]:
        return self._rng.uniform(0, 1000)

    def connection_count(self) -> Tuple[int, int]:
        return len(self._names), len(self._names)


# ===========================================================================
#  Colour helpers
# ===========================================================================

def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)


_VIRIDIS = [
    (0.00, (68,   1,  84)),
    (0.25, (59,  82, 139)),
    (0.50, (33, 145, 140)),
    (0.75, (94, 201,  98)),
    (1.00, (253, 231, 37)),
]


def _cmap_qcolor(t: float) -> QColor:
    t = max(0.0, min(1.0, t))
    for i in range(len(_VIRIDIS) - 1):
        t0, c0 = _VIRIDIS[i]
        t1, c1 = _VIRIDIS[i + 1]
        if t <= t1:
            s = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return QColor(_lerp(c0[0], c1[0], s),
                          _lerp(c0[1], c1[1], s),
                          _lerp(c0[2], c1[2], s))
    _, c = _VIRIDIS[-1]
    return QColor(*c)


# ===========================================================================
#  HyCal map widget
# ===========================================================================

class HyCalMapWidget(QWidget):
    module_hovered = pyqtSignal(str)

    _SHRINK = 0.92

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(500, 500)

        self._modules: List[Module] = []
        self._values: Dict[str, float] = {}
        self._vmin = 0.0
        self._vmax = 1000.0
        self._hovered: Optional[str] = None
        self._rects: Dict[str, QRectF] = {}
        self._layout_dirty = True

    def set_modules(self, modules: List[Module]):
        self._modules = modules
        self._layout_dirty = True
        self.update()

    def set_values(self, values: Dict[str, float]):
        self._values = values
        self.update()

    def set_range(self, vmin: float, vmax: float):
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def auto_range(self):
        vals = list(self._values.values())
        if vals:
            self._vmin = min(vals)
            self._vmax = max(vals)
            if self._vmin == self._vmax:
                self._vmax = self._vmin + 1.0
        self.update()
        return self._vmin, self._vmax

    # -- layout --

    def _recompute_layout(self):
        self._rects.clear()
        det = [m for m in self._modules if m.mod_type != "LMS"]
        if not det:
            return
        w, h = self.width(), self.height()
        margin, top, bot = 12, 10, 50
        pw, ph = w - 2 * margin, h - top - bot
        x0 = min(m.x - m.sx / 2 for m in det)
        x1 = max(m.x + m.sx / 2 for m in det)
        y0 = min(m.y - m.sy / 2 for m in det)
        y1 = max(m.y + m.sy / 2 for m in det)
        sc = min(pw / (x1 - x0), ph / (y1 - y0))
        dw, dh = (x1 - x0) * sc, (y1 - y0) * sc
        ox = margin + (pw - dw) / 2
        oy = top + (ph - dh) / 2
        shrink = self._SHRINK
        for m in det:
            mw, mh = m.sx * sc * shrink, m.sy * sc * shrink
            cx = ox + (m.x - x0) * sc
            cy = oy + (y1 - m.y) * sc
            self._rects[m.name] = QRectF(cx - mw / 2, cy - mh / 2, mw, mh)
        self._layout_dirty = False

    def resizeEvent(self, event):
        self._layout_dirty = True
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
            p.end()
            return

        vmin, vmax = self._vmin, self._vmax
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is not None:
                t = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
                p.fillRect(rect, _cmap_qcolor(t))
            else:
                p.fillRect(rect, QColor("#15181d"))

        # Hover highlight
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(QColor("#58a6ff"), 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

        # Colour bar
        cb_w = min(400, w - 80)
        cb_h = 14
        cb_x = (w - cb_w) / 2
        cb_y = h - 40
        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in _VIRIDIS:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(QRectF(cb_x, cb_y, cb_w, cb_h), QBrush(grad))
        p.setPen(QPen(QColor("#555"), 0.5))
        p.drawRect(QRectF(cb_x, cb_y, cb_w, cb_h))
        p.setPen(QColor("#8b949e"))
        p.setFont(QFont("Monospace", 9))
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, 60, 14),
                   Qt.AlignmentFlag.AlignLeft, f"{vmin:.0f}")
        p.drawText(QRectF(cb_x + cb_w - 60, cb_y + cb_h + 2, 60, 14),
                   Qt.AlignmentFlag.AlignRight, f"{vmax:.0f}")
        p.end()

    # -- mouse --

    def mouseMoveEvent(self, event):
        pos = event.position()
        found = None
        for name, rect in self._rects.items():
            if rect.contains(pos):
                found = name
                break
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                self.module_hovered.emit(found)


# ===========================================================================
#  Main window
# ===========================================================================

class ScalerMapWindow(QMainWindow):

    def __init__(self, modules: List[Module], epics_source, simulation: bool):
        super().__init__()
        self._modules = modules
        self._ep = epics_source
        self._simulation = simulation
        self._scalable = [m for m in modules
                          if m.mod_type in ("PbWO4", "PbGlass")]
        self._values: Dict[str, float] = {}
        self._polling = True

        self._build_ui()

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._refresh)
        self._timer.start(POLL_INTERVAL_MS)
        self._refresh()

    def _build_ui(self):
        self.setWindowTitle("HyCal Scaler Map" +
                            ("  [SIMULATION]" if self._simulation
                             else "  [REALTIME]"))
        self.resize(800, 860)
        self._apply_dark_palette()

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL SCALER MAP")
        lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        lbl.setStyleSheet("color:#58a6ff;")
        top.addWidget(lbl)

        mode = "SIMULATION" if self._simulation else "REALTIME"
        mode_clr = "#d29922" if self._simulation else "#3fb950"
        mode_lbl = QLabel(mode)
        mode_lbl.setFont(QFont("Monospace", 10, QFont.Weight.Bold))
        mode_lbl.setStyleSheet(f"color:{mode_clr};")
        top.addWidget(mode_lbl)
        top.addStretch()

        self._poll_btn = self._make_btn("Polling: ON", "#3fb950",
                                        self._toggle_polling)
        top.addWidget(self._poll_btn)
        top.addWidget(self._make_btn("Refresh Now", "#c9d1d9",
                                     self._refresh))
        root.addLayout(top)

        # -- map --
        self._map = HyCalMapWidget()
        self._map.set_modules(self._modules)
        self._map.module_hovered.connect(self._on_hover)
        root.addWidget(self._map, stretch=1)

        # -- range controls --
        ctrl = QHBoxLayout()
        ctrl.addWidget(self._styled_label("Range:"))

        self._min_edit = self._styled_edit("0")
        self._max_edit = self._styled_edit("1000")
        ctrl.addWidget(self._min_edit)
        ctrl.addWidget(self._styled_label("-"))
        ctrl.addWidget(self._max_edit)
        ctrl.addWidget(self._make_btn("Apply", "#c9d1d9",
                                      self._apply_range))
        ctrl.addWidget(self._make_btn("Auto", "#d29922",
                                      self._auto_range))
        ctrl.addStretch()

        self._conn_lbl = QLabel("EPICS: --")
        self._conn_lbl.setFont(QFont("Monospace", 10))
        self._conn_lbl.setStyleSheet("color:#8b949e;")
        ctrl.addWidget(self._conn_lbl)
        root.addLayout(ctrl)

        # -- info bar --
        self._info = QLabel("Hover over a module")
        self._info.setFont(QFont("Monospace", 11))
        self._info.setStyleSheet(
            "QLabel{background:#161b22;color:#c9d1d9;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}")
        self._info.setFixedHeight(28)
        root.addWidget(self._info)

    # -- helpers --

    def _make_btn(self, text: str, fg: str, slot) -> QPushButton:
        btn = QPushButton(text)
        btn.setStyleSheet(
            f"QPushButton{{background:#21262d;color:{fg};"
            f"border:1px solid #30363d;padding:5px 14px;"
            f"font:bold 11px Monospace;border-radius:4px;}}"
            f"QPushButton:hover{{background:#30363d;}}")
        btn.clicked.connect(slot)
        return btn

    def _styled_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 11))
        lbl.setStyleSheet("color:#c9d1d9;")
        return lbl

    def _styled_edit(self, text: str) -> QLineEdit:
        e = QLineEdit(text)
        e.setFixedWidth(70)
        e.setFont(QFont("Monospace", 11))
        e.setStyleSheet(
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}")
        e.returnPressed.connect(self._apply_range)
        return e

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

    # -- actions --

    def _refresh(self):
        for m in self._scalable:
            v = self._ep.get(m.name)
            if v is not None:
                self._values[m.name] = float(v)
        self._map.set_values(self._values)

        n_ok, n_total = self._ep.connection_count()
        fg = "#3fb950" if n_ok == n_total else (
             "#d29922" if n_ok > 0 else "#8b949e")
        self._conn_lbl.setText(f"EPICS: {n_ok}/{n_total}")
        self._conn_lbl.setStyleSheet(f"color:{fg};font:10px Monospace;")

        if self._values:
            lo = min(self._values.values())
            hi = max(self._values.values())
            self.statusBar().showMessage(
                f"Data: {lo:.0f} .. {hi:.0f}    "
                f"Channels: {len(self._values)}")

    def _toggle_polling(self):
        self._polling = not self._polling
        if self._polling:
            self._timer.start(POLL_INTERVAL_MS)
            self._poll_btn.setText("Polling: ON")
            self._poll_btn.setStyleSheet(
                self._poll_btn.styleSheet().replace("#f85149", "#3fb950"))
        else:
            self._timer.stop()
            self._poll_btn.setText("Polling: OFF")
            self._poll_btn.setStyleSheet(
                self._poll_btn.styleSheet().replace("#3fb950", "#f85149"))

    def _apply_range(self):
        try:
            vmin = float(self._min_edit.text())
            vmax = float(self._max_edit.text())
            if vmin < vmax:
                self._map.set_range(vmin, vmax)
        except ValueError:
            pass

    def _auto_range(self):
        vmin, vmax = self._map.auto_range()
        self._min_edit.setText(f"{vmin:.0f}")
        self._max_edit.setText(f"{vmax:.0f}")

    def _on_hover(self, name: str):
        parts = [name]
        for m in self._modules:
            if m.name == name:
                parts.append(f"({m.mod_type})")
                break
        v = self._values.get(name)
        if v is not None:
            parts.append(f"{v:.1f}")
        self._info.setText("    ".join(parts))


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="HyCal FADC Scaler Map")
    ap.add_argument("--sim", action="store_true",
                    help="Simulation mode (random values, no EPICS)")
    ap.add_argument("--database", type=Path, default=MODULES_JSON,
                    help="Path to hycal_modules.json")
    args = ap.parse_args()

    modules = load_modules(args.database)
    print(f"Loaded {len(modules)} modules")

    if args.sim:
        ep = SimulatedScalerEPICS(modules)
    else:
        try:
            ep = RealScalerEPICS(modules)
        except ImportError:
            print("ERROR: pyepics not available. Use --sim or install pyepics.")
            sys.exit(1)

    app = QApplication(sys.argv)
    win = ScalerMapWindow(modules, ep, simulation=args.sim)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
