#!/usr/bin/env python3
"""
HyCal Snake Scan -- Module Scanner (PyQt6)
==========================================
PyQt6 GUI that drives the HyCal transporter in a snake pattern so the
beam centres on each scanned module, dwells for a configurable time,
then advances to the next module.

Includes a live FADC scaler overlay so the beam spot is visible as a
hot region on the HyCal map.

Usage
-----
    python hycal_snake_scan.py                          # simulation
    python hycal_snake_scan.py --expert                  # expert operator
    python hycal_snake_scan.py --observer                # read-only monitor

Coordinate system
-----------------
    ptrans_x, ptrans_y = (-126.75, 10.11)  -->  beam at HyCal centre (0,0)
    ptrans_x = BEAM_CENTER_X + module_x
    ptrans_y = BEAM_CENTER_Y - module_y

Writable PVs (the ONLY PVs this tool writes to):
    ptrans_x.VAL / ptrans_y.VAL    -- absolute set-point
    ptrans_x.SPMG / ptrans_y.SPMG  -- motor mode  Stop(0) Pause(1) Move(2) Go(3)

Requirements
------------
    Python 3.8+, PyQt6
    pyepics  (only for --expert / --observer mode)
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Dict, List, Optional

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QGroupBox, QPushButton, QLabel, QComboBox, QSpinBox,
    QDoubleSpinBox, QTextEdit, QProgressBar, QMessageBox, QSplitter,
    QSizePolicy, QFrame, QDialog, QLineEdit, QScrollArea, QSlider,
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QFont

from scan_utils import (
    C, Module, load_modules, module_to_ptrans, ptrans_to_module,
    ptrans_in_limits, filter_scan_modules, DARK_QSS,
    BEAM_CENTER_X, BEAM_CENTER_Y, DEFAULT_DB_PATH,
)
from scan_epics import (
    SPMG, SPMG_LABELS, epics_move_to, epics_stop,
)
from scan_engine import (
    ScanState, ScanEngine, build_scan_path, estimate_scan_time,
    DEFAULT_DWELL, DEFAULT_POS_THRESHOLD, DEFAULT_BEAM_THRESHOLD,
    DEFAULT_VELO_X, DEFAULT_VELO_Y, MAX_LG_LAYERS,
)
from scan_geoview import HyCalScanMapWidget, PALETTES, PALETTE_NAMES
from scan_gui_common import (
    PATHS_FILE, SCALER_POLL_MS, POLL_MS,
    open_session_log, format_log_line, append_log_line,
    build_position_check_panel, EncoderDriftChecker,
    load_profiles, setup_motor_epics, setup_scaler_epics,
)


# ============================================================================
#  MODULE INFO DIALOG
# ============================================================================

class ModuleInfoDialog(QDialog):
    """Pop-up showing module details with a Move To button.

    Call :meth:`setModule` to refresh the content for a different module
    without closing and re-opening the dialog.
    """

    _FIELDS = ("Scaler", "Name", "Type", "Sector", "Row/Col", "Size", "HyCal", "Ptrans", "In limits")

    def __init__(self, ep, log_fn, parent=None):
        super().__init__(parent)
        self._mod: Optional[Module] = None
        self._ep = ep
        self._log = log_fn
        self.setStyleSheet(DARK_QSS)
        self.setFixedWidth(360)

        lo = QVBoxLayout(self)

        grid = QGridLayout()
        grid.setSpacing(4)
        self._value_labels: Dict[str, QLabel] = {}
        for r, label in enumerate(self._FIELDS):
            lk = QLabel(f"{label}:")
            lk.setStyleSheet(f"color: {C.DIM}; font: 13pt 'Consolas';")
            lk.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(lk, r, 0)
            lv = QLabel("--")
            grid.addWidget(lv, r, 1)
            self._value_labels[label] = lv
        lo.addLayout(grid)

        lo.addSpacing(8)

        btn_row = QHBoxLayout()
        self._btn_move = QPushButton("Move To")
        self._btn_move.setProperty("cssClass", "accent")
        self._btn_move.clicked.connect(self._doMove)
        btn_row.addWidget(self._btn_move)
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_row.addWidget(btn_close)
        lo.addLayout(btn_row)

    def setModule(self, mod: Module, scaler_value: Optional[float] = None):
        self._mod = mod
        px, py = module_to_ptrans(mod.x, mod.y)
        in_limits = ptrans_in_limits(px, py)
        vals = {
            "Scaler":    f"{scaler_value:.1f}" if scaler_value is not None else "--",
            "Name":      mod.name,
            "Type":      mod.mod_type,
            "Sector":    mod.sector or "--",
            "Row/Col":   f"{mod.row} / {mod.col}" if mod.row else "--",
            "Size":      f"{mod.sx:.2f} x {mod.sy:.2f} mm",
            "HyCal":     f"({mod.x:.2f}, {mod.y:.2f}) mm",
            "Ptrans":    f"({px:.2f}, {py:.2f}) mm",
            "In limits": "Yes" if in_limits else "No",
        }
        for label, lv in self._value_labels.items():
            lv.setText(vals.get(label, "--"))
            if label == "In limits" and not in_limits:
                lv.setStyleSheet(f"color: {C.RED}; font: bold 13pt 'Consolas';")
            elif label == "Scaler" and scaler_value is not None:
                lv.setStyleSheet(f"color: {C.GREEN}; font: bold 13pt 'Consolas';")
            else:
                lv.setStyleSheet("")
        self._btn_move.setText(f"Move To {mod.name}")
        self._btn_move.setEnabled(in_limits)
        self.setWindowTitle(f"Module {mod.name}")

    def _doMove(self):
        mod = self._mod
        if not mod: return
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if epics_move_to(self._ep, px, py):
            win = self.parent()
            if hasattr(win, '_setTarget'):
                win._setTarget(px, py, mod.name)
        else:
            self._log(f"BLOCKED: ptrans({px:.3f}, {py:.3f}) outside limits", level="error")


# ============================================================================
#  MAIN WINDOW
# ============================================================================


class SnakeScanWindow(QMainWindow):
    _logSignal = pyqtSignal(str, str)
    AUTOGEN = "(autogen)"
    NONE = "(none)"

    def __init__(self, motor_ep, scaler_ep, simulation, all_modules,
                 profiles=None, observer=False):
        super().__init__()
        self.ep = motor_ep
        self.scaler_ep = scaler_ep
        self.simulation = simulation
        self.observer = observer
        self.all_modules = all_modules
        self._profiles = profiles or {}
        self._active_profile = self.NONE
        self._lg_layers = 0

        glass = [m for m in all_modules if m.mod_type == "PbGlass"]
        self._lg_sx = glass[0].sx if glass else 38.15
        self._lg_sy = glass[0].sy if glass else 38.15

        self._mod_by_name = {m.name: m for m in all_modules}
        self._log_lines = []
        self._log_file = open_session_log("snake_scan", self.simulation, self.observer)

        self.scan_modules = []
        self.engine = ScanEngine(motor_ep, self.scan_modules, self._log)
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(self.engine.path)}
        self._scan_names = {m.name for m in self.scan_modules}
        self._selected_start_idx = 0
        self._selected_mod_name = None
        self._mod_dlg: Optional[ModuleInfoDialog] = None
        self._status_labels: Dict[str, QLabel] = {}

        self._encoder_checker = EncoderDriftChecker()

        # target position — set once when a move is commanded
        self._target_px: Optional[float] = None
        self._target_py: Optional[float] = None
        self._target_name: str = ""
        self._last_scan_idx: int = -1  # track scan engine moves

        self._logSignal.connect(self._appendLog)
        self._buildUI()

        if self.observer:
            self._disableControls()
        if not self.simulation and not self.observer:
            disc = self.ep.disconnected_pvs()
            if disc:
                self._disableControls()
                QMessageBox.critical(self, "PV Connection Error",
                    "Not connected:\n" + "\n".join(f"  {p}" for p in disc))

        # main poll timer (5 Hz)
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._poll)
        self._timer.start(POLL_MS)

        # scaler poll timer
        self._scaler_timer = QTimer(self)
        self._scaler_timer.timeout.connect(self._pollScalers)
        self._scaler_timer.start(SCALER_POLL_MS)
        self._pollScalers()  # initial read

    # =======================================================================
    #  Layout — 16:9
    #
    #  [TOP BAR: title | mode | ===BEAM=== | state              ]
    #  [LEFT half                 | RIGHT half                   ]
    #  [  HyCal geo view          |  Controls (scrollable)       ]
    #  [  legend row              |  ─────────────────────────── ]
    #  [  scaler controls         |  Event Log                   ]
    # =======================================================================

    def _buildUI(self):
        if self.observer:       suffix = "  [OBSERVER]"
        elif self.simulation:   suffix = "  [SIMULATION]"
        else:                   suffix = "  [EXPERT OPERATOR]"
        self.setWindowTitle("HyCal Snake Scan" + suffix)
        self.setStyleSheet(DARK_QSS)
        self.resize(1600, 900)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── top bar ──────────────────────────────────────────────────────
        top = QWidget()
        top.setFixedHeight(48)
        top.setStyleSheet("background: #0d1520;")
        tl = QHBoxLayout(top)
        tl.setContentsMargins(12, 0, 12, 0)

        lbl = QLabel("HYCAL SNAKE SCAN")
        lbl.setStyleSheet(f"color: {C.GREEN}; font: bold 17pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl)

        if self.observer:       mt, mf = "OBSERVER", C.ORANGE
        elif self.simulation:   mt, mf = "SIMULATION", C.YELLOW
        else:                   mt, mf = "EXPERT", C.GREEN
        lbl_mode = QLabel(mt)
        lbl_mode.setStyleSheet(f"color: {mf}; font: bold 13pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl_mode)
        tl.addSpacing(16)

        # prominent beam current
        beam_frame = QFrame()
        beam_frame.setStyleSheet(
            "QFrame { background: #161b22; border: 1px solid #30363d; border-radius: 4px; }")
        beam_frame.setFixedHeight(36)
        bf_layout = QHBoxLayout(beam_frame)
        bf_layout.setContentsMargins(10, 0, 10, 0)
        bf_layout.setSpacing(6)
        beam_icon = QLabel("BEAM")
        beam_icon.setStyleSheet("color: #8b949e; font: bold 12pt 'Consolas'; background: transparent; border: none;")
        bf_layout.addWidget(beam_icon)
        self._lbl_beam_val = QLabel("-- nA")
        self._lbl_beam_val.setStyleSheet(
            f"color: {C.GREEN}; font: bold 18pt 'Consolas'; background: transparent; border: none;")
        self._lbl_beam_val.setMinimumWidth(140)
        bf_layout.addWidget(self._lbl_beam_val)
        self._lbl_beam_status = QLabel("")
        self._lbl_beam_status.setStyleSheet(
            "color: transparent; font: bold 13pt 'Consolas'; background: transparent; border: none;")
        bf_layout.addWidget(self._lbl_beam_status)
        tl.addWidget(beam_frame)

        tl.addStretch()

        self._lbl_state = QLabel("IDLE")
        self._lbl_state.setStyleSheet(
            f"color: {C.DIM}; font: bold 15pt 'Consolas'; background: transparent;")
        tl.addWidget(self._lbl_state)

        root.addWidget(top)

        # ── main body — horizontal splitter ──────────────────────────────
        body_splitter = QSplitter(Qt.Orientation.Horizontal)
        body_splitter.setContentsMargins(6, 4, 6, 6)

        # --- LEFT half: map + legend + scaler controls ---
        left = QWidget()
        left_lo = QVBoxLayout(left)
        left_lo.setContentsMargins(0, 0, 0, 0)
        left_lo.setSpacing(2)

        self._canvas_label = QLabel()
        self._canvas_label.setStyleSheet(f"color: {C.ACCENT}; font: bold 13pt 'Consolas';")
        left_lo.addWidget(self._canvas_label)

        self._map = HyCalScanMapWidget(self.all_modules)
        self._map.moduleClicked.connect(self._onCanvasClick)
        left_lo.addWidget(self._map, stretch=1)

        # reset button overlaid at bottom-right corner, outside the map drawing area
        self._btn_reset_view = QPushButton("Reset", self._map)
        self._btn_reset_view.setFixedSize(56, 28)
        self._btn_reset_view.setStyleSheet(
            f"QPushButton{{background:rgba(22,27,34,220);color:{C.DIM};"
            f"border:1px solid #30363d;border-radius:2px;padding:0;"
            f"font:12pt Consolas;}}"
            f"QPushButton:hover{{color:{C.TEXT};border-color:{C.ACCENT};}}")
        self._btn_reset_view.clicked.connect(self._map.resetView)
        self._map.installEventFilter(self)

        # legend row
        leg = QHBoxLayout()
        leg.setSpacing(4); leg.setContentsMargins(0, 0, 0, 0)
        for label, colour in [("Todo", C.MOD_TODO), ("Skipped", C.MOD_SKIPPED),
                               ("Moving", C.MOD_CURRENT), ("Dwell", C.MOD_DWELL),
                               ("Done", C.MOD_DONE), ("Error", C.MOD_ERROR),
                               ("Start", C.MOD_SELECTED), ("PbGlass", C.MOD_GLASS)]:
            sw = QLabel(); sw.setFixedSize(10, 10)
            sw.setStyleSheet(f"background: {colour}; border: none;")
            leg.addWidget(sw)
            ll = QLabel(label); ll.setStyleSheet(f"color: {C.DIM}; font: 12pt 'Consolas';")
            leg.addWidget(ll)
        leg.addStretch()
        left_lo.addLayout(leg)

        # scaler controls row
        sc_row = QHBoxLayout()
        sc_row.setSpacing(4); sc_row.setContentsMargins(0, 2, 0, 0)

        self._btn_scaler_toggle = QPushButton("Scalers: ON")
        self._btn_scaler_toggle.setStyleSheet(self._scaler_btn_ss(True))
        self._btn_scaler_toggle.setFixedHeight(28)
        self._btn_scaler_toggle.clicked.connect(self._toggleScaler)
        sc_row.addWidget(self._btn_scaler_toggle)

        self._btn_scaler_auto = QPushButton("Auto")
        self._btn_scaler_auto.setFixedHeight(28)
        self._btn_scaler_auto.clicked.connect(self._toggleScalerAuto)
        self._scaler_auto_on = True
        self._updateScalerAutoBtn()
        sc_row.addWidget(self._btn_scaler_auto)

        self._scaler_min_edit = QLineEdit("0")
        self._scaler_min_edit.setFixedWidth(50); self._scaler_min_edit.setFixedHeight(28)
        self._scaler_min_edit.setFont(QFont("Consolas", 8))
        self._scaler_min_edit.setStyleSheet(
            "QLineEdit{background:#161b22;color:#c9d1d9;border:1px solid #30363d;border-radius:2px;padding:1px 4px;}")
        self._scaler_min_edit.returnPressed.connect(self._applyScalerRange)
        sc_row.addWidget(self._scaler_min_edit)

        sc_row.addWidget(QLabel("-"))

        self._scaler_max_edit = QLineEdit("1000")
        self._scaler_max_edit.setFixedWidth(50); self._scaler_max_edit.setFixedHeight(28)
        self._scaler_max_edit.setFont(QFont("Consolas", 8))
        self._scaler_max_edit.setStyleSheet(self._scaler_min_edit.styleSheet())
        self._scaler_max_edit.returnPressed.connect(self._applyScalerRange)
        sc_row.addWidget(self._scaler_max_edit)

        btn_apply = QPushButton("Apply"); btn_apply.setFixedHeight(28)
        btn_apply.clicked.connect(self._applyScalerRange)
        sc_row.addWidget(btn_apply)

        self._btn_scaler_log = QPushButton("Log: OFF"); self._btn_scaler_log.setFixedHeight(28)
        self._btn_scaler_log.setStyleSheet(self._small_btn_ss(C.DIM))
        self._btn_scaler_log.clicked.connect(self._toggleScalerLog)
        sc_row.addWidget(self._btn_scaler_log)

        self._btn_palette = QPushButton(); self._btn_palette.setFixedSize(90, 28)
        self._btn_palette.setToolTip("Click to cycle colour palette")
        self._btn_palette.clicked.connect(self._cycleScalerPalette)
        self._updatePaletteBtn()
        sc_row.addWidget(self._btn_palette)

        sc_row.addStretch()
        left_lo.addLayout(sc_row)

        body_splitter.addWidget(left)

        # --- RIGHT half: controls (top 2/3) + log (bottom 1/3) ---
        right_splitter = QSplitter(Qt.Orientation.Vertical)

        # top area: two columns — left: scan/direct/position, right: motor status
        ctrl_columns = QWidget()
        ctrl_cols_lo = QHBoxLayout(ctrl_columns)
        ctrl_cols_lo.setContentsMargins(0, 0, 0, 0)
        ctrl_cols_lo.setSpacing(4)

        # left column: scan control + direct control + position check
        left_col_scroll = QScrollArea()
        left_col_scroll.setWidgetResizable(True)
        left_col_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_col_widget = QWidget()
        left_col_lo = QVBoxLayout(left_col_widget)
        left_col_lo.setSpacing(4)
        left_col_lo.setContentsMargins(0, 0, 0, 0)
        self._buildScanControl(left_col_lo)
        self._buildDirectControl(left_col_lo)
        left_col_lo.addStretch()
        left_col_scroll.setWidget(left_col_widget)
        ctrl_cols_lo.addWidget(left_col_scroll, stretch=1)

        # right column: position check + motor status
        right_col_scroll = QScrollArea()
        right_col_scroll.setWidgetResizable(True)
        right_col_scroll.setFrameShape(QFrame.Shape.NoFrame)
        right_col_widget = QWidget()
        right_col_lo = QVBoxLayout(right_col_widget)
        right_col_lo.setSpacing(4)
        right_col_lo.setContentsMargins(0, 0, 0, 0)
        self._buildPositionCheck(right_col_lo)
        self._buildMotorStatus(right_col_lo)
        self._buildScalerControl(right_col_lo)
        right_col_lo.addStretch()
        right_col_scroll.setWidget(right_col_widget)
        ctrl_cols_lo.addWidget(right_col_scroll, stretch=1)

        right_splitter.addWidget(ctrl_columns)

        # event log
        log_group = QGroupBox("Event Log")
        log_layout = QVBoxLayout(log_group)
        log_layout.setContentsMargins(4, 4, 4, 4)
        self._log_text = QTextEdit()
        self._log_text.setReadOnly(True)
        log_layout.addWidget(self._log_text)
        right_splitter.addWidget(log_group)

        right_splitter.setStretchFactor(0, 3)  # controls
        right_splitter.setStretchFactor(1, 2)  # log

        body_splitter.addWidget(right_splitter)
        body_splitter.setStretchFactor(0, 1)  # left half
        body_splitter.setStretchFactor(1, 1)  # right half

        root.addWidget(body_splitter, stretch=1)
        self._updateCanvasLabel()

    # -- scaler control helpers ----------------------------------------------

    @staticmethod
    def _scaler_btn_ss(on):
        fg = C.GREEN if on else C.RED
        return (f"QPushButton{{background:#21262d;color:{fg};"
                f"border:1px solid #30363d;padding:1px 8px;"
                f"font:bold 12pt Consolas;border-radius:2px;}}"
                f"QPushButton:hover{{background:#30363d;}}")

    @staticmethod
    def _small_btn_ss(fg):
        return (f"QPushButton{{background:#21262d;color:{fg};"
                f"border:1px solid #30363d;padding:1px 8px;"
                f"font:bold 12pt Consolas;border-radius:2px;}}"
                f"QPushButton:hover{{background:#30363d;}}")

    def _toggleScaler(self):
        on = not self._map._scaler_enabled
        self._map.setScalerEnabled(on)
        self._btn_scaler_toggle.setText("Scalers: ON" if on else "Scalers: OFF")
        self._btn_scaler_toggle.setStyleSheet(self._scaler_btn_ss(on))

    def _toggleScalerAuto(self):
        self._scaler_auto_on = not self._scaler_auto_on
        self._map.setScalerAutoRange(self._scaler_auto_on)
        self._updateScalerAutoBtn()
        if self._scaler_auto_on:
            vmin, vmax = self._map.scalerRange()
            self._scaler_min_edit.setText(f"{vmin:.0f}")
            self._scaler_max_edit.setText(f"{vmax:.0f}")

    def _updateScalerAutoBtn(self):
        if self._scaler_auto_on:
            self._btn_scaler_auto.setStyleSheet(
                "QPushButton{background:#d29922;color:#0d1117;"
                "border:1px solid #d29922;padding:1px 8px;"
                "font:bold 12pt Consolas;border-radius:2px;}"
                "QPushButton:hover{background:#e0a82b;}")
        else:
            self._btn_scaler_auto.setStyleSheet(self._small_btn_ss(C.YELLOW))

    def _applyScalerRange(self):
        try:
            vmin = float(self._scaler_min_edit.text())
            vmax = float(self._scaler_max_edit.text())
            if vmin < vmax:
                self._map.setScalerRange(vmin, vmax)
                self._scaler_auto_on = False
                self._map.setScalerAutoRange(False)
                self._updateScalerAutoBtn()
        except ValueError:
            pass

    def _toggleScalerLog(self):
        on = not self._map._scaler_log
        self._map.setScalerLogScale(on)
        self._btn_scaler_log.setText("Log: ON" if on else "Log: OFF")
        self._btn_scaler_log.setStyleSheet(
            self._small_btn_ss(C.ACCENT if on else C.DIM))

    def _cycleScalerPalette(self):
        self._map.cyclePalette()
        self._updatePaletteBtn()

    def _updatePaletteBtn(self):
        """Set the palette button background to a CSS linear-gradient of the current palette."""
        idx = self._map._palette_idx
        stops = list(PALETTES.values())[idx]
        css_stops = ", ".join(
            f"rgb({r},{g},{b}) {int(t * 100)}%" for t, (r, g, b) in stops)
        self._btn_palette.setStyleSheet(
            f"QPushButton{{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,{self._palette_grad_stops(stops)});"
            f"border:1px solid #30363d;border-radius:2px;color:#c9d1d9;"
            f"font:bold 11pt Consolas;padding:0 4px;}}"
            f"QPushButton:hover{{border-color:#58a6ff;}}")
        self._btn_palette.setText(PALETTE_NAMES[idx])

    @staticmethod
    def _palette_grad_stops(stops):
        parts = []
        for t, (r, g, b) in stops:
            parts.append(f"stop:{t:.2f} rgb({r},{g},{b})")
        return ",".join(parts)

    def _pollScalers(self):
        vals = self.scaler_ep.get_all()
        if vals:
            self._map.setScalerValues(vals)
            if self._scaler_auto_on:
                vmin, vmax = self._map.scalerRange()
                self._scaler_min_edit.setText(f"{vmin:.0f}")
                self._scaler_max_edit.setText(f"{vmax:.0f}")

    # -- Scan Control --------------------------------------------------------

    def _buildScanControl(self, parent):
        sc = QGroupBox("Scan Control")
        lo = QVBoxLayout(sc)

        r = QHBoxLayout(); r.addWidget(QLabel("Path:"))
        self._profile_combo = QComboBox()
        self._profile_combo.addItems([self.NONE, self.AUTOGEN] + sorted(self._profiles.keys()))
        self._profile_combo.setCurrentText(self.NONE)
        self._profile_combo.activated.connect(self._onPathProfileChanged)
        r.addWidget(self._profile_combo, stretch=1); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("LG layers (0-2):"))
        self._lg_spin = QSpinBox(); self._lg_spin.setRange(0, MAX_LG_LAYERS)
        self._lg_spin.valueChanged.connect(self._onLgLayersChanged)
        r.addWidget(self._lg_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Start:"))
        names = [m.name for m in self.engine.path]
        self._start_combo = QComboBox(); self._start_combo.setEditable(True)
        self._start_combo.addItems(names); self._start_combo.setMinimumWidth(80)
        self._start_combo.activated.connect(self._onStartSelected)
        r.addWidget(self._start_combo)
        r.addWidget(QLabel("Count:"))
        self._count_spin = QSpinBox(); self._count_spin.setRange(0, len(names))
        self._count_spin.setSpecialValueText("all")
        self._count_spin.valueChanged.connect(lambda _: self._drawPathPreview())
        r.addWidget(self._count_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Dwell (s):"))
        self._dwell_spin = QDoubleSpinBox(); self._dwell_spin.setRange(1, 9999)
        self._dwell_spin.setValue(DEFAULT_DWELL); self._dwell_spin.setDecimals(0)
        r.addWidget(self._dwell_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Pos. threshold (mm):"))
        self._thresh_spin = QDoubleSpinBox(); self._thresh_spin.setRange(0.01, 10.0)
        self._thresh_spin.setValue(DEFAULT_POS_THRESHOLD); self._thresh_spin.setSingleStep(0.1)
        self._thresh_spin.setDecimals(2)
        r.addWidget(self._thresh_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Beam threshold (nA):"))
        self._beam_thresh_spin = QDoubleSpinBox(); self._beam_thresh_spin.setRange(0.0, 1000.0)
        self._beam_thresh_spin.setValue(DEFAULT_BEAM_THRESHOLD)
        self._beam_thresh_spin.setSingleStep(0.1); self._beam_thresh_spin.setDecimals(2)
        self._beam_thresh_spin.setSpecialValueText("off")
        r.addWidget(self._beam_thresh_spin); lo.addLayout(r)

        bf = QHBoxLayout()
        self._btn_start = QPushButton("Start Scan"); self._btn_start.setProperty("cssClass", "green")
        self._btn_start.clicked.connect(self._cmdStart); bf.addWidget(self._btn_start)
        self._btn_pause = QPushButton("Pause"); self._btn_pause.setProperty("cssClass", "warn")
        self._btn_pause.clicked.connect(self._cmdPause); bf.addWidget(self._btn_pause)
        self._btn_stop = QPushButton("Stop"); self._btn_stop.setProperty("cssClass", "danger")
        self._btn_stop.clicked.connect(self._cmdStop); bf.addWidget(self._btn_stop)
        lo.addLayout(bf)

        bf2 = QHBoxLayout()
        self._btn_skip = QPushButton("Skip Module"); self._btn_skip.clicked.connect(self._cmdSkip)
        bf2.addWidget(self._btn_skip)
        self._btn_ack = QPushButton("Ack Error"); self._btn_ack.setProperty("cssClass", "warn")
        self._btn_ack.clicked.connect(self._cmdAckError); bf2.addWidget(self._btn_ack)
        lo.addLayout(bf2)

        r = QHBoxLayout()
        self._lbl_progress = QLabel("Progress: --/--"); r.addWidget(self._lbl_progress)
        self._progress_bar = QProgressBar(); self._progress_bar.setMaximumWidth(140)
        r.addWidget(self._progress_bar); lo.addLayout(r)

        self._lbl_current = QLabel("Current:  --"); lo.addWidget(self._lbl_current)
        self._lbl_eta = QLabel("ETA:      --")
        self._lbl_eta.setStyleSheet(f"color: {C.DIM};"); lo.addWidget(self._lbl_eta)
        self._lbl_dwell_cd = QLabel("")
        self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};"); lo.addWidget(self._lbl_dwell_cd)

        parent.addWidget(sc)

    def _buildDirectControl(self, parent):
        dc = QGroupBox("Direct Control"); lo = QVBoxLayout(dc)
        self._btn_move = QPushButton("Move to Starting Point")
        self._btn_move.clicked.connect(self._cmdMoveToModule); lo.addWidget(self._btn_move)
        self._btn_reset = QPushButton("Reset to Beam Center")
        self._btn_reset.setProperty("cssClass", "accent")
        self._btn_reset.clicked.connect(self._cmdResetCenter); lo.addWidget(self._btn_reset)
        parent.addWidget(dc)

    def _buildMotorStatus(self, parent):
        ms = QGroupBox("Motor Status"); lo = QVBoxLayout(ms)
        self._motor_state_labels = {}
        for title, axis, fields in [
            ("X Motor", "x", [
                ("Encoder", "x_encoder"), ("RBV", "x_rbv"), ("VAL", "x_val"),
                ("MOVN", "x_movn"), ("SPMG", "x_spmg"), ("VELO", "x_velo"),
                ("ACCL", "x_accl"), ("TDIR", "x_tdir"), ("MSTA", "x_msta"), ("ATHM", "x_athm")]),
            ("Y Motor", "y", [
                ("Encoder", "y_encoder"), ("RBV", "y_rbv"), ("VAL", "y_val"),
                ("MOVN", "y_movn"), ("SPMG", "y_spmg"), ("VELO", "y_velo"),
                ("ACCL", "y_accl"), ("TDIR", "y_tdir"), ("MSTA", "y_msta"), ("ATHM", "y_athm")])]:
            # title row with inline state badge
            tr = QHBoxLayout()
            tl = QLabel(title)
            tl.setStyleSheet(f"color: {C.ACCENT}; font: bold 13pt 'Consolas';")
            tr.addWidget(tl)
            sl = QLabel("Idle")
            sl.setStyleSheet(f"color: {C.DIM}; font: bold 12pt 'Consolas'; "
                             f"background: #21262d; border: 1px solid #30363d; "
                             f"border-radius: 3px; padding: 1px 6px;")
            tr.addWidget(sl)
            tr.addStretch()
            self._motor_state_labels[axis] = sl
            lo.addLayout(tr)
            g = QGridLayout(); g.setSpacing(2)
            half = (len(fields) + 1) // 2
            for i, (label, key) in enumerate(fields):
                c = 0 if i < half else 2; r = i % half
                ln = QLabel(f"{label}:")
                ln.setStyleSheet(f"color: {C.DIM}; font: 12pt 'Consolas';")
                ln.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                g.addWidget(ln, r, c)
                lv = QLabel("--"); lv.setMinimumWidth(90)
                g.addWidget(lv, r, c + 1)
                self._status_labels[key] = lv
            lo.addLayout(g)
        parent.addWidget(ms)

    def _buildPositionCheck(self, parent):
        labels = build_position_check_panel(parent)
        self._lbl_expected = labels["target"]
        self._lbl_actual = labels["actual"]
        self._lbl_error = labels["diff"]
        self._lbl_drift = labels["drift"]

    def _buildScalerControl(self, parent):
        sc = QGroupBox("Scalers"); lo = QVBoxLayout(sc)

        r = QHBoxLayout()
        btn_refresh = QPushButton("Refresh Now")
        btn_refresh.clicked.connect(self._pollScalers)
        r.addWidget(btn_refresh)
        r.addStretch()
        lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Poll:"))
        self._scaler_interval_slider = QSlider(Qt.Orientation.Horizontal)
        self._scaler_interval_slider.setRange(20, 100)  # 2.0s - 10.0s in 0.1s steps
        self._scaler_interval_slider.setValue(50)        # default 5.0s
        self._scaler_interval_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self._scaler_interval_slider.setTickInterval(10)
        self._scaler_interval_slider.valueChanged.connect(self._onScalerIntervalChanged)
        r.addWidget(self._scaler_interval_slider)
        self._lbl_scaler_interval = QLabel("5.0 s")
        self._lbl_scaler_interval.setMinimumWidth(40)
        r.addWidget(self._lbl_scaler_interval)
        lo.addLayout(r)

        parent.addWidget(sc)

    def _onScalerIntervalChanged(self, val):
        sec = val / 10.0
        self._lbl_scaler_interval.setText(f"{sec:.1f} s")
        self._scaler_timer.setInterval(int(sec * 1000))

    def _disableControls(self):
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_skip, self._btn_ack, self._btn_move, self._btn_reset):
            w.setEnabled(False)
        for w in (self._start_combo, self._profile_combo, self._lg_spin,
                  self._count_spin, self._dwell_spin, self._thresh_spin,
                  self._beam_thresh_spin):
            w.setEnabled(False)

    # -- canvas helpers ------------------------------------------------------

    def _updateCanvasLabel(self):
        n_pwo4 = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        base = f"Scan Path: {n_pwo4} PbWO4 + {n_lg} LG" if n_lg else f"Scan Path: {n_pwo4} PbWO4"
        if not self.scan_modules:
            base = "Scan Path: none"
        elif self.engine.path and 0 <= self._selected_start_idx < len(self.engine.path):
            start_name = self.engine.path[self._selected_start_idx].name
            base += f"  start: {start_name}"
        self._canvas_label.setText(f" {base} ")

    def _updateCanvas(self):
        if self.observer:
            colors = {}
            for m in self.all_modules:
                if m.mod_type == "LMS": continue
                if m.mod_type == "PbGlass": colors[m.name] = C.MOD_GLASS
                elif m.mod_type == "PbWO4": colors[m.name] = C.MOD_PWO4_BG
            for m in self.scan_modules:
                colors[m.name] = C.MOD_TODO
            self._map.setModuleColors(colors)
            rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
            self._map.setMarkerPosition(*ptrans_to_module(rx, ry))
            self._map.update(); return

        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        idle = eng.state in (ScanState.IDLE, ScanState.COMPLETED)

        # Build scan-state colours.  When scaler overlay is on, these are
        # drawn as borders (heat map fill stays visible).  When off, they
        # are drawn as fills.
        colors = {}
        for m in self.all_modules:
            if m.name in self._scan_names or m.mod_type == "LMS": continue
            colors[m.name] = C.MOD_EXCLUDED if running else (
                C.MOD_GLASS if m.mod_type == "PbGlass" else C.MOD_PWO4_BG if m.mod_type == "PbWO4" else C.MOD_LMS)

        count = self._count_spin.value()
        si = self._selected_start_idx
        ei = min(si + count, len(eng.path)) if count > 0 else len(eng.path)
        if not idle:
            si = eng.current_idx; ei = getattr(eng, '_end_idx', len(eng.path))
        for i, mod in enumerate(eng.path):
            if i == eng.current_idx and eng.state == ScanState.DWELLING:
                colors[mod.name] = C.MOD_DWELL
            elif i == eng.current_idx and eng.state in (ScanState.MOVING, ScanState.PAUSED):
                colors[mod.name] = C.MOD_CURRENT
            elif i in eng.error_modules: colors[mod.name] = C.MOD_ERROR
            elif i in eng.completed: colors[mod.name] = C.MOD_DONE
            elif idle and i == self._selected_start_idx: colors[mod.name] = C.MOD_SELECTED
            elif i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
            else: colors[mod.name] = C.MOD_TODO
        self._map.setModuleColors(colors)
        if idle:
            self._drawPathPreview(); self._map.setDashPreview([])
        elif running:
            self._map.setPathPreview([])
            self._map.setDashPreview([self._map.modCenter(eng.path[i]) for i in range(eng.current_idx + 1, ei)])
        rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
        self._map.setMarkerPosition(*ptrans_to_module(rx, ry))
        self._map.update()

    def _drawPathPreview(self):
        path = self.engine.path
        s = self._selected_start_idx
        if s >= len(path): self._map.setPathPreview([]); return
        c = self._count_spin.value()
        e = min(s + c, len(path)) if c > 0 else len(path)
        self._map.setPathPreview([self._map.modCenter(path[i]) for i in range(s, e)])

    def _onCanvasClick(self, name):
        if self._selected_mod_name == name:
            self._selected_mod_name = None; self._map.setHighlight(None)
            self._updateCanvasLabel(); return
        self._selected_mod_name = name
        if name in self._scan_name_to_idx:
            self._selected_start_idx = self._scan_name_to_idx[name]
            idx = self._start_combo.findText(name)
            if idx >= 0: self._start_combo.setCurrentIndex(idx)
            self._drawPathPreview()
        self._map.setHighlight(name); self._updateCanvasLabel()

        idle = self.engine.state in (ScanState.IDLE, ScanState.COMPLETED)
        if idle and not self.observer:
            mod = self._mod_by_name.get(name)
            if mod:
                if self._mod_dlg is None:
                    self._mod_dlg = ModuleInfoDialog(self.ep, self._log, parent=self)
                sv = self._map._scaler_values.get(mod.name)
                self._mod_dlg.setModule(mod, sv)
                self._mod_dlg.show()
                self._mod_dlg.raise_()

    # -- commands ------------------------------------------------------------

    def _onStartSelected(self, _):
        name = self._start_combo.currentText()
        for i, m in enumerate(self.engine.path):
            if m.name == name:
                self._selected_start_idx = i
                self._drawPathPreview()
                self._updateCanvasLabel()
                break

    def _onPathProfileChanged(self, _):
        name = self._profile_combo.currentText()
        if name == self._active_profile: return
        self._active_profile = name
        if name == self.AUTOGEN:
            self._lg_spin.setEnabled(True); self._onLgLayersChanged(force=True); return
        self._lg_spin.setEnabled(False)
        if name == self.NONE:
            self._setPath([])
            self._log("Path: none (direct control only)")
            return
        mod_by_name = {m.name: m for m in self.all_modules}
        path_mods = [mod_by_name[n] for n in self._profiles.get(name, []) if n in mod_by_name]
        if not path_mods:
            self._log(f"Profile '{name}' empty", level="error"); return
        self._setPath(path_mods)
        self._log(f"Path profile: {name} ({len(path_mods)} modules)")

    def _onLgLayersChanged(self, value=0, force=False):
        if self._active_profile != self.AUTOGEN: return
        nl = self._lg_spin.value()
        if nl == self._lg_layers and not force: return
        self._lg_layers = nl
        mods = filter_scan_modules(self.all_modules, nl, self._lg_sx, self._lg_sy)
        # generate the snake path once at autogen — from now on the order is fixed
        path, n_unopt = build_scan_path(mods)
        if n_unopt:
            self._log(f"WARNING: {n_unopt} modules with unoptimized path", level="warn")
        self._setPath(path)
        np_ = sum(1 for m in path if m.mod_type == "PbWO4")
        ng = sum(1 for m in path if m.mod_type == "PbGlass")
        self._log(f"LG layers: {nl} ({np_} PbWO4 + {ng} PbGlass = {len(path)})")

    def _setPath(self, path):
        """Set the scan path. ``path`` is the final ordered list of modules."""
        self.scan_modules = path
        self._scan_names = {m.name for m in path}
        self.engine = ScanEngine(self.ep, path, self._log)
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(path)}
        self._selected_start_idx = 0
        ns = [m.name for m in path]
        self._start_combo.clear(); self._start_combo.addItems(ns)
        self._count_spin.setMaximum(len(ns)); self._count_spin.setValue(0)
        if not path:
            self._map.setPathPreview([]); self._map.setDashPreview([])
            self._map.setHighlight(None); self._selected_mod_name = None
        self._updateCanvasLabel()

    def _cmdStart(self):
        self._onStartSelected(0)
        path = self.engine.path; s = self._selected_start_idx; c = self._count_spin.value()
        e = min(s + c, len(path)) if c > 0 else len(path)
        oob = [path[i].name for i in range(s, e) if not ptrans_in_limits(*module_to_ptrans(path[i].x, path[i].y))]
        if oob:
            ns = ", ".join(oob[:5]) + (f" ... ({len(oob)} total)" if len(oob) > 5 else "")
            self._log(f"BLOCKED: {len(oob)} modules outside limits: {ns}", level="error")
            QMessageBox.critical(self, "Out of Bounds", f"{len(oob)} outside limits:\n{ns}"); return
        self.engine.dwell_time = self._dwell_spin.value()
        self.engine.pos_threshold = self._thresh_spin.value()
        self.engine.beam_threshold = self._beam_thresh_spin.value()
        self.engine.start(self._selected_start_idx, count=c)

    def _cmdPause(self):
        eng = self.engine
        if eng.state == ScanState.PAUSED:
            eng.resume_scan(); self._btn_pause.setText("Pause")
        elif eng.state in (ScanState.MOVING, ScanState.DWELLING):
            eng.pause_scan(); self._btn_pause.setText("Resume")

    def _cmdStop(self):
        if self.engine.state != ScanState.IDLE:
            self.engine.stop_scan(); self._btn_pause.setText("Pause")
        else:
            epics_stop(self.ep); self._log("Motors stopped")

    def _cmdSkip(self):      self.engine.skip_module()
    def _cmdAckError(self):  self.engine.acknowledge_error()

    def _setTarget(self, px, py, name=""):
        self._target_px = px
        self._target_py = py
        self._target_name = name

    def _cmdMoveToModule(self):
        self._onStartSelected(0)
        if not self.engine.path: return
        mod = self.engine.path[self._selected_start_idx]
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if epics_move_to(self.ep, px, py):
            self._setTarget(px, py, mod.name)
        else:
            self._log(f"BLOCKED: outside limits", level="error")

    def _cmdResetCenter(self):
        self._log(f"Resetting to beam centre ptrans({BEAM_CENTER_X}, {BEAM_CENTER_Y})")
        if epics_move_to(self.ep, BEAM_CENTER_X, BEAM_CENTER_Y):
            self._setTarget(BEAM_CENTER_X, BEAM_CENTER_Y, "Beam Center")

    # -- event filter (reposition overlay button on map resize) ---------------

    def eventFilter(self, obj, event):
        if obj is self._map and event.type() == event.Type.Resize:
            btn = self._btn_reset_view
            btn.move(self._map.width() - btn.width() - 2,
                     self._map.height() - btn.height() - 2)
        return super().eventFilter(obj, event)

    # -- logging -------------------------------------------------------------

    def _log(self, msg, level="info"):
        line = format_log_line(msg, level)
        self._log_lines.append(line)
        if self._log_file and not self._log_file.closed:
            self._log_file.write(line + "\n"); self._log_file.flush()
        self._logSignal.emit(line, level)

    def _appendLog(self, line, level):
        append_log_line(self._log_text, line, level)

    # -- polling (5 Hz) ------------------------------------------------------

    def _poll(self):
        # detect scan engine moving to a new module
        eng = self.engine
        if eng.state in (ScanState.MOVING, ScanState.DWELLING) and eng.current_idx != self._last_scan_idx:
            self._last_scan_idx = eng.current_idx
            mod = eng.current_module
            if mod:
                px, py = module_to_ptrans(mod.x, mod.y)
                self._setTarget(px, py, mod.name)
        elif eng.state == ScanState.IDLE:
            self._last_scan_idx = -1
        self._updateStatus()
        self._updateCanvas()
        self._updateScanInfo()
        self._updateButtons()
        self._updateBeamDisplay()
        self._checkEncoder()

    def _updateBeamDisplay(self):
        bc = self.ep.get("beam_cur", None)
        if bc is None:
            self._lbl_beam_val.setText("-- nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.DIM}; font: bold 18pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("")
            return
        thresh = self._beam_thresh_spin.value()
        tripped = self.engine.beam_tripped
        if tripped:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.RED}; font: bold 18pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("TRIP")
            self._lbl_beam_status.setStyleSheet(f"color: {C.RED}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
        elif thresh > 0 and bc < thresh:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.YELLOW}; font: bold 18pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("LOW")
            self._lbl_beam_status.setStyleSheet(f"color: {C.YELLOW}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
        else:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.GREEN}; font: bold 18pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("")

    def _checkEncoder(self):
        self._encoder_checker.update(self.ep, self._log, self._lbl_drift)

    def _updateStatus(self):
        for key, lbl in self._status_labels.items():
            val = self.ep.get(key, "--")
            if val == "--" or val is None:
                lbl.setText("--"); lbl.setStyleSheet(f"color: {C.DIM};"); continue
            if key.endswith("_msta"):        txt = f"0x{int(val):X}"
            elif key.endswith("_spmg"):      txt = f"{SPMG_LABELS.get(int(val), '?')}({int(val)})"
            elif key.endswith(("_movn", "_athm", "_tdir")): txt = str(int(val))
            elif isinstance(val, float):     txt = f"{val:.3f}"
            else:                            txt = str(val)
            fg = C.TEXT
            if key.endswith("_movn") and int(val) == 1: fg = C.YELLOW
            elif key.endswith("_spmg"):
                sv = int(val)
                fg = C.RED if sv == SPMG.STOP else C.ORANGE if sv == SPMG.PAUSE else C.GREEN
            lbl.setText(txt); lbl.setStyleSheet(f"color: {fg};")

        # update motor state badges
        for axis, sl in self._motor_state_labels.items():
            spmg = self.ep.get(f"{axis}_spmg", None)
            movn = self.ep.get(f"{axis}_movn", None)
            if spmg is not None and int(spmg) == SPMG.STOP:
                sl.setText("Stop"); fg, bg = C.RED, "#3d1214"
            elif spmg is not None and int(spmg) == SPMG.PAUSE:
                sl.setText("Pause"); fg, bg = C.ORANGE, "#3d2a0e"
            elif movn is not None and int(movn) == 1:
                sl.setText("Moving"); fg, bg = C.YELLOW, "#3d3010"
            else:
                sl.setText("Idle"); fg, bg = C.DIM, "#21262d"
            sl.setStyleSheet(f"color: {fg}; font: bold 12pt 'Consolas'; "
                             f"background: {bg}; border: 1px solid #30363d; "
                             f"border-radius: 3px; padding: 1px 6px;")
        # position check — target is set when a move is commanded
        rx, ry = self.ep.get("x_rbv", 0.0), self.ep.get("y_rbv", 0.0)
        self._lbl_actual.setText(f"Actual: ({rx:.3f}, {ry:.3f})")
        px, py = self._target_px, self._target_py
        if px is not None and py is not None:
            err = math.sqrt((rx - px)**2 + (ry - py)**2)
            name_html = f' <b style="color:{C.ACCENT}">{self._target_name}</b>' if self._target_name else ""
            self._lbl_expected.setText(f"Target:   ({px:.3f}, {py:.3f}){name_html}")
            # ETA from distance and motor velocities
            vx = self.ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
            vy = self.ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
            dx, dy = abs(rx - px), abs(ry - py)
            eta_sec = max(dx / vx if vx > 0 else 0, dy / vy if vy > 0 else 0)
            if eta_sec >= 60:
                eta_str = f" ({int(eta_sec)//60}m {int(eta_sec)%60}s)"
            elif eta_sec >= 1:
                eta_str = f" ({eta_sec:.0f}s)"
            else:
                eta_str = ""
            scanning = self.engine.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
            if scanning:
                ef = C.RED if err > self.engine.pos_threshold else C.GREEN
                self._lbl_error.setText(f"Diff:   {err:.3f} mm{eta_str}")
                self._lbl_error.setStyleSheet(f"color: {ef}; font: bold 13pt 'Consolas';")
            else:
                self._lbl_error.setText(f"Diff:   {err:.3f} mm{eta_str}")
                self._lbl_error.setStyleSheet(f"color: {C.DIM}; font: bold 13pt 'Consolas';")
        else:
            self._lbl_expected.setText("Target: --")
            self._lbl_error.setText("Diff:   --")
            self._lbl_error.setStyleSheet(f"color: {C.DIM}; font: bold 13pt 'Consolas';")

    def _updateScanInfo(self):
        eng = self.engine
        sc = {ScanState.IDLE: C.DIM, ScanState.MOVING: C.YELLOW, ScanState.DWELLING: C.GREEN,
              ScanState.PAUSED: C.ORANGE, ScanState.ERROR: C.RED, ScanState.COMPLETED: C.ACCENT}
        self._lbl_state.setText(eng.state)
        self._lbl_state.setStyleSheet(f"color: {sc.get(eng.state, C.DIM)}; font: bold 15pt 'Consolas'; background: transparent;")
        done = len(eng.completed)
        s = getattr(eng, '_start_idx', 0); e = getattr(eng, '_end_idx', len(eng.path))
        total = e - s
        self._lbl_progress.setText(f"Progress: {done}/{total}")
        self._progress_bar.setMaximum(max(total, 1)); self._progress_bar.setValue(done)
        mod = eng.current_module
        self._lbl_current.setText(f"Current:  {mod.name}" if mod else "Current:  --")
        if eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED):
            eta = eng.eta_seconds; h, rem = divmod(int(eta), 3600); m, s = divmod(rem, 60)
            self._lbl_eta.setText(f"ETA:      {h}h {m:02d}m {s:02d}s")
        elif eng.state == ScanState.IDLE and eng.path:
            vx = self.ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
            vy = self.ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
            eta = estimate_scan_time(eng.path, self._selected_start_idx,
                                     self._count_spin.value(), self._dwell_spin.value(), vx, vy)
            if eta > 0:
                h, rem = divmod(int(eta), 3600); m, s = divmod(rem, 60)
                self._lbl_eta.setText(f"ETA:      ~{h}h {m:02d}m {s:02d}s")
            else: self._lbl_eta.setText("ETA:      --")
        else: self._lbl_eta.setText("ETA:      --")
        if eng.beam_tripped:
            self._lbl_dwell_cd.setText("Dwell:    BEAM TRIP -- waiting")
            self._lbl_dwell_cd.setStyleSheet(f"color: {C.RED}; font: bold 13pt 'Consolas';")
        elif eng.state == ScanState.DWELLING:
            self._lbl_dwell_cd.setText(f"Dwell:    {eng.dwell_remaining:.1f}s remaining")
            self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};")
        else:
            self._lbl_dwell_cd.setText(""); self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};")

    def _updateButtons(self):
        if self.observer: return
        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        has_path = len(eng.path) > 0
        self._btn_start.setEnabled(not running and has_path)
        self._btn_pause.setEnabled(running)
        self._btn_stop.setEnabled(True)
        self._btn_skip.setEnabled(eng.state == ScanState.DWELLING)
        self._btn_ack.setEnabled(eng.state == ScanState.ERROR)
        self._start_combo.setEnabled(not running and has_path)
        self._count_spin.setEnabled(not running and has_path)
        self._profile_combo.setEnabled(not running)
        self._lg_spin.setEnabled(not running and self._active_profile == self.AUTOGEN)
        # dwell/thresholds are only applied at scan start, so lock them while running
        self._dwell_spin.setEnabled(not running)
        self._thresh_spin.setEnabled(not running)
        self._beam_thresh_spin.setEnabled(not running)
        self._btn_move.setEnabled(not running and has_path)
        self._btn_reset.setEnabled(not running)

    def closeEvent(self, e):
        self._timer.stop()
        self._scaler_timer.stop()
        if self._log_file and not self._log_file.closed:
            self._log_file.close()
        self._log_file = None
        super().closeEvent(e)


# ============================================================================
#  MAIN
# ============================================================================


def main():
    parser = argparse.ArgumentParser(description="HyCal Snake Scan")
    parser.add_argument("--expert", action="store_true")
    parser.add_argument("--observer", action="store_true")
    parser.add_argument("--database", default=DEFAULT_DB_PATH)
    parser.add_argument("--paths", default=PATHS_FILE)
    args = parser.parse_args()

    all_modules = load_modules(args.database)
    by_type: Dict[str, int] = {}
    for m in all_modules:
        by_type[m.mod_type] = by_type.get(m.mod_type, 0) + 1
    print(f"Loaded {len(all_modules)} modules from {args.database}")
    for t, n in sorted(by_type.items()):
        print(f"  {t}: {n}")

    profiles = load_profiles(args.paths)
    if profiles:
        print(f"Loaded {len(profiles)} path profiles")

    observer = args.observer
    simulation = not args.expert and not observer

    motor_ep = setup_motor_epics(observer, simulation)
    scaler_ep = setup_scaler_epics(simulation, all_modules)

    app = QApplication(sys.argv)
    win = SnakeScanWindow(motor_ep, scaler_ep, simulation, all_modules,
                          profiles, observer=observer)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
