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
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

#local path for testing on farm
#sys.path.append('/home/wrightso/.local/bin/*')

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit,
)
from PyQt6.QtCore import QTimer
from PyQt6.QtGui import QFont

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget, PALETTES,
    apply_theme_palette, set_theme, available_themes, THEME,
)


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_modules.json"

SCALER_PV = "B_DET_HYCAL_FADC_{label}:c"
POLL_INTERVAL_MS = 2_500   # 1 seconds


# ===========================================================================
#  EPICS interfaces
# ===========================================================================

class RealScalerEPICS:
    """Read scaler PVs via pyepics."""

    def __init__(self, modules: List[Module]):
        import epics as _epics
        self._pvs: Dict[str, object] = {}
        for m in modules:
            if m.mod_type in ("PbWO4", "PbGlass", "LMS"):
                pv = _epics.PV(SCALER_PV.format(label=m.name), connection_timeout=2.0)
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
                       if m.mod_type in ("PbWO4", "PbGlass", "LMS")]

    def get(self, name: str) -> Optional[float]:
        return self._rng.uniform(0, 1000)

    def connection_count(self) -> Tuple[int, int]:
        return len(self._names), len(self._names)


# ===========================================================================
#  HyCal map widget  (subclass customises min size + vmax default)
# ===========================================================================

class ScalerMapWidget(HyCalMapWidget):
    """Simple value → colour map with palette cycle and log scale."""

    def __init__(self, parent=None):
        super().__init__(parent, min_size=(500, 500), include_lms=True)
        self._vmax = 1000.0   # sensible default for kHz rates

    def _fmt_value(self, v: float) -> str:
        return f"{v:.0f}"


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
                          if m.mod_type in ("PbWO4", "PbGlass", "LMS")]
        self._values: Dict[str, float] = {}
        self._polling = True
        self._palette_idx = 0
        self._auto_range_on = True

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
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL SCALER MAP")
        lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        lbl.setStyleSheet(f"color:{THEME.ACCENT};")
        top.addWidget(lbl)

        mode = "SIMULATION" if self._simulation else "REALTIME"
        mode_clr = THEME.WARN if self._simulation else THEME.SUCCESS
        mode_lbl = QLabel(mode)
        mode_lbl.setFont(QFont("Monospace", 10, QFont.Weight.Bold))
        mode_lbl.setStyleSheet(f"color:{mode_clr};")
        top.addWidget(mode_lbl)
        top.addStretch()

        self._poll_btn = self._make_btn("Polling: ON", THEME.SUCCESS,
                                        self._toggle_polling)
        top.addWidget(self._poll_btn)
        top.addWidget(self._make_btn("Refresh Now", THEME.TEXT,
                                     self._refresh))
        root.addLayout(top)

        # -- map --
        # Only show modules that actually have scaler PVs; SCINT (V1–V4)
        # has no scaler rate, so rendering them here would misleadingly
        # grey them out. LMS is still filtered by include_lms=False.
        self._map = ScalerMapWidget()
        self._map.set_modules(self._scalable)
        self._map.moduleHovered.connect(self._on_hover)
        self._map.paletteClicked.connect(self._cycle_palette)
        root.addWidget(self._map, stretch=1)

        # -- range controls --
        ctrl = QHBoxLayout()
        ctrl.addWidget(self._styled_label("Range:"))

        self._min_edit = self._styled_edit("0")
        self._max_edit = self._styled_edit("1000")
        ctrl.addWidget(self._min_edit)
        ctrl.addWidget(self._styled_label("-"))
        ctrl.addWidget(self._max_edit)
        ctrl.addWidget(self._make_btn("Apply", THEME.TEXT,
                                      self._apply_range))
        self._auto_btn = self._make_btn("Auto Scale", THEME.WARN,
                                        self._toggle_auto_range)
        ctrl.addWidget(self._auto_btn)
        self._update_auto_btn()
        self._log_btn = self._make_btn("Log: OFF", THEME.TEXT_DIM,
                                       self._toggle_log)
        ctrl.addWidget(self._log_btn)
        ctrl.addStretch()

        self._conn_lbl = QLabel("EPICS: --")
        self._conn_lbl.setFont(QFont("Monospace", 10))
        self._conn_lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};")
        ctrl.addWidget(self._conn_lbl)
        root.addLayout(ctrl)

        # -- info bar --
        self._info = QLabel("Hover over a module")
        self._info.setFont(QFont("Monospace", 11))
        self._info.setStyleSheet(
            f"QLabel{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"padding:4px 8px;border:1px solid {THEME.BORDER};"
            f"border-radius:8px;}}")
        self._info.setFixedHeight(28)
        root.addWidget(self._info)

    # -- helpers --

    def _make_btn(self, text: str, fg: str, slot) -> QPushButton:
        btn = QPushButton(text)
        btn.setStyleSheet(
            f"QPushButton{{background:{THEME.BUTTON};color:{fg};"
            f"border:1px solid {THEME.BORDER};padding:5px 14px;"
            f"font:bold 11px Monospace;border-radius:8px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")
        btn.clicked.connect(slot)
        return btn

    def _styled_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 11))
        lbl.setStyleSheet(f"color:{THEME.TEXT};")
        return lbl

    def _styled_edit(self, text: str) -> QLineEdit:
        e = QLineEdit(text)
        e.setFixedWidth(70)
        e.setFont(QFont("Monospace", 11))
        e.setStyleSheet(
            f"QLineEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};border-radius:8px;"
            f"padding:2px 6px;}}")
        e.returnPressed.connect(self._apply_range)
        return e

    # -- actions --

    def _refresh(self):
        W_totalSum = 0
        topSum = 0
        botSum = 0
        leftSum = 0
        rightSum = 0
        for m in self._scalable:
            v = self._ep.get(m.name)
            if v is not None:
                self._values[m.name] = float(v)
                if "W" in m.name and int(m.name.strip("W")!=1125):
                    W_totalSum += v
                    if(int(m.name.strip("W"))<578):
                        topSum += v
                    else:
                        botSum += v
                    if(int(((float(m.name.strip("W"))/34.0)%1)*100)<=50):
                        leftSum += v
                    else:
                        rightSum += v


        self._map.set_values(self._values)

        #Convert Sums of rates to kHz
        W_totalSum = W_totalSum/1000.0
        y_asym = (topSum-botSum)/1000.0
        x_asym = (rightSum-leftSum)/1000.0
        #Get Center of the Rate Relative to the center of the beam hole
        #x_COM = x_asym/(20.5*17)
        #y_COM = y_asym/(20.5*17)

        if self._auto_range_on and self._values:
            self._do_auto_range()

        n_ok, n_total = self._ep.connection_count()
        fg = THEME.SUCCESS if n_ok == n_total else (
             THEME.WARN if n_ok > 0 else THEME.TEXT_DIM)
        self._conn_lbl.setText(f"EPICS: {n_ok}/{n_total}")
        self._conn_lbl.setStyleSheet(f"color:{fg};font:10px Monospace;")

        if self._values:
            lo = min(self._values.values())/1000.0
            hi = max(self._values.values())/1000.0
            self.statusBar().showMessage(
                f"Data: {lo:.0f}kHz .. {hi:.0f}kHz  "
                f"Channels: {len(self._values)}  "
                f"PbWO4 Total: {W_totalSum:.2f}kHz  "
                f"Ave: {W_totalSum/1152:3f}kHz  "
                f"Asym (kHz): [{x_asym:.3f}, {y_asym:.3f}]")
                #f"CoR (mm): [{x_COM:.3f},{y_COM:.3f}]")

    def _toggle_polling(self):
        self._polling = not self._polling
        if self._polling:
            self._timer.start(POLL_INTERVAL_MS)
            self._poll_btn.setText("Polling: ON")
            self._poll_btn.setStyleSheet(
                self._poll_btn.styleSheet().replace(THEME.DANGER, THEME.SUCCESS))
        else:
            self._timer.stop()
            self._poll_btn.setText("Polling: OFF")
            self._poll_btn.setStyleSheet(
                self._poll_btn.styleSheet().replace(THEME.SUCCESS, THEME.DANGER))

    def _apply_range(self):
        try:
            vmin = float(self._min_edit.text())
            vmax = float(self._max_edit.text())
            if vmin < vmax:
                self._map.set_range(vmin, vmax)
                self._auto_range_on = False
                self._update_auto_btn()
        except ValueError:
            pass

    def _toggle_auto_range(self):
        self._auto_range_on = not self._auto_range_on
        self._update_auto_btn()
        if self._auto_range_on:
            self._do_auto_range()

    def _do_auto_range(self):
        vmin, vmax = self._map.auto_range()
        self._min_edit.setText(f"{vmin:.0f}")
        self._max_edit.setText(f"{vmax:.0f}")

    def _update_auto_btn(self):
        if self._auto_range_on:
            self._auto_btn.setStyleSheet(
                f"QPushButton{{background:{THEME.WARN};color:{THEME.BG};"
                f"border:1px solid {THEME.WARN};padding:5px 14px;"
                f"font:bold 11px Monospace;border-radius:8px;}}"
                f"QPushButton:hover{{background:{THEME.WARN};}}")
        else:
            self._auto_btn.setStyleSheet(
                f"QPushButton{{background:{THEME.BUTTON};color:{THEME.WARN};"
                f"border:1px solid {THEME.BORDER};padding:5px 14px;"
                f"font:bold 11px Monospace;border-radius:8px;}}"
                f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")

    def _cycle_palette(self):
        self._palette_idx = (self._palette_idx + 1) % len(PALETTES)
        self._map.set_palette(self._palette_idx)

    def _toggle_log(self):
        on = not self._map.is_log_scale()
        self._map.set_log_scale(on)
        if on:
            self._log_btn.setText("Log: ON")
            self._log_btn.setStyleSheet(
                self._log_btn.styleSheet().replace(THEME.TEXT_DIM, THEME.ACCENT))
        else:
            self._log_btn.setText("Log: OFF")
            self._log_btn.setStyleSheet(
                self._log_btn.styleSheet().replace(THEME.ACCENT, THEME.TEXT_DIM))

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
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

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
