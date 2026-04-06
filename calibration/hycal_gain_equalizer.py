#!/usr/bin/env python3
"""
HyCal Gain Equalizer (PyQt6)
=============================
Automatic gain equalization for HyCal crystal modules.  Moves the beam
to each module, collects peak height histograms from prad2_server, finds
the right edge of the Bremsstrahlung spectrum, and adjusts HV via
prad2hvd until the edge converges to a target ADC value.

Shares scan_utils, scan_epics, scan_engine, scan_geoview, and
gain_scanner modules with hycal_snake_scan.

Usage
-----
    python hycal_gain_equalizer.py                     # simulation
    python hycal_gain_equalizer.py --expert             # expert operator
    python hycal_gain_equalizer.py --observer            # read-only monitor
"""

from __future__ import annotations

import argparse
import html as html_mod
import json
import math
import os
import sys
from datetime import datetime
from typing import Dict, List, Optional

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QGroupBox, QPushButton, QLabel, QComboBox, QSpinBox,
    QDoubleSpinBox, QTextEdit, QProgressBar, QMessageBox, QSplitter,
    QSizePolicy, QFrame, QLineEdit, QScrollArea, QSlider,
)
from PyQt6.QtCore import Qt, QRectF, QTimer, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QPainter, QPen

from scan_utils import (
    C, Module, load_modules, module_to_ptrans, ptrans_to_module,
    ptrans_in_limits, filter_scan_modules, DARK_QSS,
    BEAM_CENTER_X, BEAM_CENTER_Y, DEFAULT_DB_PATH,
)
from scan_epics import (
    SPMG, SPMG_LABELS, MotorEPICS, ObserverEPICS, SimulatedMotorEPICS,
    ScalerPVGroup, SimulatedScalerEPICS,
    epics_move_to, epics_stop,
)
from scan_engine import (
    build_scan_path,
    DEFAULT_POS_THRESHOLD, DEFAULT_BEAM_THRESHOLD,
    DEFAULT_VELO_X, DEFAULT_VELO_Y, MAX_LG_LAYERS,
)
from scan_geoview import HyCalScanMapWidget, PALETTES, PALETTE_NAMES
from gain_scanner import (
    GainScanEngine, GainScanState, ServerClient, HVClient,
)


# ============================================================================
#  HISTOGRAM WIDGET
# ============================================================================

class HistogramWidget(QWidget):
    """Lightweight bar chart for peak height histogram display."""

    PAD_L, PAD_R, PAD_T, PAD_B = 50, 12, 28, 24

    def __init__(self, parent=None):
        super().__init__(parent)
        self._bins: List[int] = []
        self._target_bin: Optional[int] = None  # vertical target line
        self._edge_bin: Optional[int] = None     # detected edge marker
        self._title: str = ""
        self._info: str = ""                     # e.g. "V=1525.0  edge=3200"
        self.setMinimumHeight(140)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def setData(self, bins: List[int], target_bin: Optional[int] = None,
                edge_bin: Optional[int] = None):
        self._bins = bins
        self._target_bin = target_bin
        self._edge_bin = edge_bin
        self.update()

    def setTitle(self, text: str):
        self._title = text; self.update()

    def setInfo(self, text: str):
        self._info = text; self.update()

    def clear(self):
        self._bins = []; self._target_bin = None; self._edge_bin = None
        self._title = ""; self._info = ""
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor("#0d1117"))

        L, R, T, B = self.PAD_L, self.PAD_R, self.PAD_T, self.PAD_B
        pw, ph = w - L - R, h - T - B
        if pw < 10 or ph < 10 or not self._bins:
            # title and info even when empty
            p.setPen(QColor(C.ACCENT))
            p.setFont(QFont("Consolas", 11, QFont.Weight.Bold))
            p.drawText(QRectF(L, 2, pw, T - 2), Qt.AlignmentFlag.AlignLeft, self._title)
            p.setPen(QColor(C.DIM))
            p.setFont(QFont("Consolas", 10))
            p.drawText(QRectF(L, 2, pw, T - 2), Qt.AlignmentFlag.AlignRight, self._info)
            p.end(); return

        import math as _math
        bins = self._bins
        n = len(bins)
        vmax = max(bins) if bins else 1
        if vmax == 0: vmax = 1
        log_vmax = _math.log10(max(vmax, 1))

        # title + info
        p.setPen(QColor(C.ACCENT))
        p.setFont(QFont("Consolas", 11, QFont.Weight.Bold))
        p.drawText(QRectF(L, 2, pw, T - 2), Qt.AlignmentFlag.AlignLeft, self._title)
        p.setPen(QColor(C.DIM))
        p.setFont(QFont("Consolas", 10))
        p.drawText(QRectF(L, 2, pw, T - 2), Qt.AlignmentFlag.AlignRight, self._info)

        # axes
        p.setPen(QPen(QColor("#30363d"), 1))
        p.drawLine(L, T, L, T + ph)
        p.drawLine(L, T + ph, L + pw, T + ph)

        # bars (log y scale)
        bar_w = pw / n
        p.setPen(Qt.PenStyle.NoPen)
        for i, v in enumerate(bins):
            if v <= 0: continue
            frac = _math.log10(v) / log_vmax if log_vmax > 0 else 0
            bh = frac * ph
            x = L + i * bar_w
            y = T + ph - bh
            p.fillRect(QRectF(x, y, max(bar_w - 0.5, 0.5), bh), QColor(C.ACCENT))

        # target line (red dashed vertical)
        if self._target_bin is not None and 0 <= self._target_bin < n:
            tx = L + (self._target_bin + 0.5) * bar_w
            p.setPen(QPen(QColor(C.RED), 1.5, Qt.PenStyle.DashLine))
            p.drawLine(int(tx), T, int(tx), T + ph)

        # edge marker (green solid vertical)
        if self._edge_bin is not None and 0 <= self._edge_bin < n:
            ex = L + (self._edge_bin + 0.5) * bar_w
            p.setPen(QPen(QColor(C.GREEN), 2))
            p.drawLine(int(ex), T, int(ex), T + ph)

        # y-axis labels (log scale: 1, 10, 100, ...)
        p.setPen(QColor(C.DIM))
        p.setFont(QFont("Consolas", 9))
        p.drawText(QRectF(0, T - 2, L - 4, 14),
                   Qt.AlignmentFlag.AlignRight, f"{vmax}")
        p.drawText(QRectF(0, T + ph - 7, L - 4, 14),
                   Qt.AlignmentFlag.AlignRight, "1")
        # grid lines at powers of 10
        p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
        decade = 10
        while decade < vmax:
            frac = _math.log10(decade) / log_vmax
            gy = T + ph - frac * ph
            p.drawLine(L + 1, int(gy), L + pw, int(gy))
            p.setPen(QColor(C.DIM))
            p.setFont(QFont("Consolas", 8))
            p.drawText(QRectF(0, gy - 7, L - 4, 14),
                       Qt.AlignmentFlag.AlignRight, f"{decade}")
            p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
            decade *= 10

        p.end()


# ============================================================================
#  MAIN WINDOW
# ============================================================================

SCALER_POLL_MS = 5_000
PATHS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paths.json")


class GainEqualizerWindow(QMainWindow):
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
        self._log_lines: List[str] = []

        log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
        os.makedirs(log_dir, exist_ok=True)
        self._log_file = open(os.path.join(log_dir,
            datetime.now().strftime("gain_eq_%Y%m%d_%H%M%S.log")), "w")

        self.scan_modules: List[Module] = []
        self._scan_names: set = set()
        self._scan_name_to_idx: Dict[str, int] = {}
        self._ordered_path: List[Module] = []
        self._use_profile_order = False
        self._selected_start_idx = 0
        self._selected_mod_name: Optional[str] = None
        self._mod_dlg = None
        self._gain_engine: Optional[GainScanEngine] = None

        self._enc_offset_x: Optional[float] = None
        self._enc_offset_y: Optional[float] = None

        self._target_px: Optional[float] = None
        self._target_py: Optional[float] = None
        self._target_name: str = ""

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

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._poll)
        self._timer.start(200)

        self._scaler_timer = QTimer(self)
        self._scaler_timer.timeout.connect(self._pollScalers)
        self._scaler_timer.start(SCALER_POLL_MS)
        self._pollScalers()

    # =======================================================================
    #  Layout
    # =======================================================================

    def _buildUI(self):
        if self.observer:       suffix = "  [OBSERVER]"
        elif self.simulation:   suffix = "  [SIMULATION]"
        else:                   suffix = "  [EXPERT OPERATOR]"
        self.setWindowTitle("HyCal Gain Equalizer" + suffix)
        self.setStyleSheet(DARK_QSS)
        self.resize(1600, 900)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── top bar ──
        top = QWidget(); top.setFixedHeight(48)
        top.setStyleSheet("background: #0d1520;")
        tl = QHBoxLayout(top); tl.setContentsMargins(12, 0, 12, 0)

        lbl = QLabel("HYCAL GAIN EQUALIZER")
        lbl.setStyleSheet(f"color: {C.GREEN}; font: bold 17pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl)

        if self.observer:       mt, mf = "OBSERVER", C.ORANGE
        elif self.simulation:   mt, mf = "SIMULATION", C.YELLOW
        else:                   mt, mf = "EXPERT", C.GREEN
        lbl_mode = QLabel(mt)
        lbl_mode.setStyleSheet(f"color: {mf}; font: bold 13pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl_mode); tl.addSpacing(16)

        # beam current
        beam_frame = QFrame()
        beam_frame.setStyleSheet("QFrame{background:#161b22;border:1px solid #30363d;border-radius:4px;}")
        beam_frame.setFixedHeight(36)
        bf_lo = QHBoxLayout(beam_frame); bf_lo.setContentsMargins(10, 0, 10, 0); bf_lo.setSpacing(6)
        bf_lo.addWidget(QLabel("BEAM"))
        self._lbl_beam_val = QLabel("-- nA")
        self._lbl_beam_val.setStyleSheet(f"color:{C.GREEN};font:bold 18pt 'Consolas';background:transparent;border:none;")
        self._lbl_beam_val.setMinimumWidth(140)
        bf_lo.addWidget(self._lbl_beam_val)
        self._lbl_beam_status = QLabel("")
        bf_lo.addWidget(self._lbl_beam_status)
        tl.addWidget(beam_frame); tl.addStretch()

        self._lbl_state = QLabel("IDLE")
        self._lbl_state.setStyleSheet(f"color:{C.DIM};font:bold 15pt 'Consolas';background:transparent;")
        tl.addWidget(self._lbl_state)
        root.addWidget(top)

        # ── body ──
        body_splitter = QSplitter(Qt.Orientation.Horizontal)
        body_splitter.setContentsMargins(6, 4, 6, 6)

        # LEFT: map + legend + scaler controls
        left = QWidget(); left_lo = QVBoxLayout(left)
        left_lo.setContentsMargins(0, 0, 0, 0); left_lo.setSpacing(2)

        self._canvas_label = QLabel()
        self._canvas_label.setStyleSheet(f"color:{C.ACCENT};font:bold 13pt 'Consolas';")
        left_lo.addWidget(self._canvas_label)

        self._map = HyCalScanMapWidget(self.all_modules)
        self._map.moduleClicked.connect(self._onCanvasClick)
        left_lo.addWidget(self._map, stretch=1)

        # reset view button
        self._btn_reset_view = QPushButton("Reset", self._map)
        self._btn_reset_view.setFixedSize(56, 28)
        self._btn_reset_view.setStyleSheet(
            f"QPushButton{{background:rgba(22,27,34,220);color:{C.DIM};"
            f"border:1px solid #30363d;border-radius:2px;padding:0;"
            f"font:10pt Consolas;}}"
            f"QPushButton:hover{{color:{C.TEXT};border-color:{C.ACCENT};}}")
        self._btn_reset_view.clicked.connect(self._map.resetView)
        self._map.installEventFilter(self)

        # legend
        leg = QHBoxLayout(); leg.setSpacing(4); leg.setContentsMargins(0, 0, 0, 0)
        for label, colour in [("Converged", C.GREEN), ("Failed", C.RED),
                               ("In progress", C.YELLOW), ("Todo", C.MOD_TODO),
                               ("Skipped", C.MOD_SKIPPED)]:
            sw = QLabel(); sw.setFixedSize(10, 10)
            sw.setStyleSheet(f"background:{colour};border:none;")
            leg.addWidget(sw)
            ll = QLabel(label); ll.setStyleSheet(f"color:{C.DIM};font:12pt 'Consolas';")
            leg.addWidget(ll)
        leg.addStretch()
        left_lo.addLayout(leg)

        # scaler controls
        sc_row = QHBoxLayout(); sc_row.setSpacing(4); sc_row.setContentsMargins(0, 2, 0, 0)

        self._btn_scaler_toggle = QPushButton("Scalers: ON")
        self._btn_scaler_toggle.setStyleSheet(self._toggle_btn_ss(True))
        self._btn_scaler_toggle.setFixedHeight(28)
        self._btn_scaler_toggle.clicked.connect(self._toggleScaler)
        sc_row.addWidget(self._btn_scaler_toggle)

        self._btn_scaler_auto = QPushButton("Auto"); self._btn_scaler_auto.setFixedHeight(28)
        self._btn_scaler_auto.clicked.connect(self._toggleScalerAuto)
        self._scaler_auto_on = True; self._updateScalerAutoBtn()
        sc_row.addWidget(self._btn_scaler_auto)

        self._scaler_min_edit = self._small_edit("0"); sc_row.addWidget(self._scaler_min_edit)
        sc_row.addWidget(QLabel("-"))
        self._scaler_max_edit = self._small_edit("1000"); sc_row.addWidget(self._scaler_max_edit)

        btn_apply = QPushButton("Apply"); btn_apply.setFixedHeight(28)
        btn_apply.clicked.connect(self._applyScalerRange); sc_row.addWidget(btn_apply)

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

        # RIGHT: control panel | histogram | event log (vertical splitter)
        right_splitter = QSplitter(Qt.Orientation.Vertical)

        # row 1: control panel (scrollable)
        ctrl_scroll = QScrollArea()
        ctrl_scroll.setWidgetResizable(True)
        ctrl_scroll.setFrameShape(QFrame.Shape.NoFrame)
        ctrl_w = QWidget(); ctrl_lo = QVBoxLayout(ctrl_w)
        ctrl_lo.setSpacing(4); ctrl_lo.setContentsMargins(0, 0, 0, 0)
        self._buildPathControl(ctrl_lo)
        self._buildGainControl(ctrl_lo)
        self._buildPositionCheck(ctrl_lo)
        ctrl_lo.addStretch()
        ctrl_scroll.setWidget(ctrl_w)
        right_splitter.addWidget(ctrl_scroll)

        # row 2: peak height histogram
        hist_group = QGroupBox("Peak Height Histogram")
        hist_lo = QVBoxLayout(hist_group); hist_lo.setContentsMargins(4, 4, 4, 4)
        self._histogram = HistogramWidget()
        hist_lo.addWidget(self._histogram)
        right_splitter.addWidget(hist_group)

        # row 3: event log
        log_group = QGroupBox("Event Log")
        log_lo = QVBoxLayout(log_group); log_lo.setContentsMargins(4, 4, 4, 4)
        self._log_text = QTextEdit(); self._log_text.setReadOnly(True)
        log_lo.addWidget(self._log_text)
        right_splitter.addWidget(log_group)

        right_splitter.setStretchFactor(0, 2)  # controls
        right_splitter.setStretchFactor(1, 3)  # histogram
        right_splitter.setStretchFactor(2, 2)  # log

        body_splitter.addWidget(right_splitter)
        body_splitter.setStretchFactor(0, 1)
        body_splitter.setStretchFactor(1, 1)
        root.addWidget(body_splitter, stretch=1)

        self._updateCanvasLabel()

    # -- small widget helpers -----------------------------------------------

    @staticmethod
    def _toggle_btn_ss(on):
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

    def _small_edit(self, text):
        e = QLineEdit(text); e.setFixedWidth(50); e.setFixedHeight(28)
        e.setFont(QFont("Consolas", 10))
        e.setStyleSheet("QLineEdit{background:#161b22;color:#c9d1d9;"
                         "border:1px solid #30363d;border-radius:2px;padding:1px 4px;}")
        e.returnPressed.connect(self._applyScalerRange)
        return e

    # -- control panel builders ---------------------------------------------

    def _buildPathControl(self, parent):
        pc = QGroupBox("Scan Path"); lo = QVBoxLayout(pc)
        self._path_group = pc

        r = QHBoxLayout(); r.addWidget(QLabel("Path:"))
        self._profile_combo = QComboBox()
        self._profile_combo.addItems([self.NONE, self.AUTOGEN] + sorted(self._profiles.keys()))
        self._profile_combo.setCurrentText(self.NONE)
        self._profile_combo.activated.connect(self._onPathProfileChanged)
        r.addWidget(self._profile_combo, stretch=1); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("LG layers:"))
        self._lg_spin = QSpinBox(); self._lg_spin.setRange(0, MAX_LG_LAYERS)
        self._lg_spin.valueChanged.connect(self._onLgLayersChanged)
        r.addWidget(self._lg_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Start:"))
        self._start_combo = QComboBox(); self._start_combo.setEditable(True)
        self._start_combo.setMinimumWidth(80)
        self._start_combo.activated.connect(self._onStartSelected)
        r.addWidget(self._start_combo)
        r.addWidget(QLabel("Count:"))
        self._count_spin = QSpinBox(); self._count_spin.setRange(0, 0)
        self._count_spin.setSpecialValueText("all")
        self._count_spin.valueChanged.connect(lambda _: self._drawPathPreview())
        r.addWidget(self._count_spin); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Pos. threshold (mm):"))
        self._thresh_spin = QDoubleSpinBox(); self._thresh_spin.setRange(0.01, 10.0)
        self._thresh_spin.setValue(DEFAULT_POS_THRESHOLD)
        self._thresh_spin.setSingleStep(0.1); self._thresh_spin.setDecimals(2)
        r.addWidget(self._thresh_spin); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Beam threshold (nA):"))
        self._beam_thresh_spin = QDoubleSpinBox(); self._beam_thresh_spin.setRange(0.0, 1000.0)
        self._beam_thresh_spin.setValue(DEFAULT_BEAM_THRESHOLD)
        self._beam_thresh_spin.setSingleStep(0.1); self._beam_thresh_spin.setDecimals(2)
        self._beam_thresh_spin.setSpecialValueText("off")
        r.addWidget(self._beam_thresh_spin); lo.addLayout(r)

        parent.addWidget(pc)

    def _buildGainControl(self, parent):
        ge = QGroupBox("Gain Equalizer"); lo = QVBoxLayout(ge)

        r = QHBoxLayout(); r.addWidget(QLabel("Server:"))
        self._ge_server_edit = QLineEdit("http://clondaq6:5051")
        r.addWidget(self._ge_server_edit); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("HV:"))
        self._ge_hv_edit = QLineEdit("ws://clonpc19:8765")
        r.addWidget(self._ge_hv_edit); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("HV Password:"))
        self._ge_hv_pw = QLineEdit(); self._ge_hv_pw.setEchoMode(QLineEdit.EchoMode.Password)
        r.addWidget(self._ge_hv_pw); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Target ADC:"))
        self._ge_target = QSpinBox(); self._ge_target.setRange(500, 4000)
        self._ge_target.setValue(3200); r.addWidget(self._ge_target)
        r.addWidget(QLabel("Min counts:"))
        self._ge_counts = QSpinBox(); self._ge_counts.setRange(100, 1000000)
        self._ge_counts.setValue(10000); self._ge_counts.setSingleStep(1000)
        r.addWidget(self._ge_counts); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Max iter:"))
        self._ge_maxiter = QSpinBox(); self._ge_maxiter.setRange(1, 50)
        self._ge_maxiter.setValue(8); r.addWidget(self._ge_maxiter)
        r.addWidget(QLabel("Tolerance:"))
        self._ge_tol = QSpinBox(); self._ge_tol.setRange(10, 500)
        self._ge_tol.setValue(50); r.addWidget(self._ge_tol); lo.addLayout(r)

        bf = QHBoxLayout()
        self._btn_start = QPushButton("Start")
        self._btn_start.setProperty("cssClass", "green")
        self._btn_start.clicked.connect(self._cmdStart); bf.addWidget(self._btn_start)
        self._btn_pause = QPushButton("Pause")
        self._btn_pause.setProperty("cssClass", "warn")
        self._btn_pause.clicked.connect(self._cmdPause); bf.addWidget(self._btn_pause)
        self._btn_stop = QPushButton("Stop")
        self._btn_stop.setProperty("cssClass", "danger")
        self._btn_stop.clicked.connect(self._cmdStop); bf.addWidget(self._btn_stop)
        lo.addLayout(bf)

        bf2 = QHBoxLayout()
        self._btn_skip = QPushButton("Skip Module")
        self._btn_skip.clicked.connect(self._cmdSkip); bf2.addWidget(self._btn_skip)
        bf2.addStretch()
        lo.addLayout(bf2)

        self._lbl_progress = QLabel("Progress: --/--"); lo.addWidget(self._lbl_progress)
        self._progress_bar = QProgressBar(); lo.addWidget(self._progress_bar)
        self._lbl_ge_status = QLabel("Idle")
        self._lbl_ge_status.setStyleSheet(f"color:{C.DIM};"); lo.addWidget(self._lbl_ge_status)
        self._lbl_ge_detail = QLabel("")
        self._lbl_ge_detail.setStyleSheet(f"color:{C.DIM};"); lo.addWidget(self._lbl_ge_detail)

        parent.addWidget(ge)

    def _buildPositionCheck(self, parent):
        pe = QGroupBox("Position Check"); lo = QVBoxLayout(pe)
        self._lbl_expected = QLabel("Target: --"); lo.addWidget(self._lbl_expected)
        self._lbl_actual = QLabel("Actual: --"); lo.addWidget(self._lbl_actual)
        self._lbl_error = QLabel("Diff:   --")
        self._lbl_error.setStyleSheet(f"font:bold 13pt 'Consolas';"); lo.addWidget(self._lbl_error)
        self._lbl_drift = QLabel("Drift:   --"); lo.addWidget(self._lbl_drift)
        parent.addWidget(pe)

    def _disableControls(self):
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_skip):
            w.setEnabled(False)

    # -- commands -----------------------------------------------------------

    def _cmdStart(self):
        if self._gain_engine and self._gain_engine.state not in (
                GainScanState.IDLE, GainScanState.COMPLETED):
            return
        if not self.scan_modules:
            self._log("Select a path first", level="error"); return
        # sync start index from combo selection
        self._onStartSelected(0)

        ro = self.simulation
        if not ro and not self._ge_hv_pw.text().strip():
            self._log("Expert mode requires HV password", level="error")
            QMessageBox.warning(self, "HV Password Required",
                                "Enter the prad2hvd password before starting in expert mode.")
            return

        try:
            server = ServerClient(self._ge_server_edit.text().strip(),
                                  log_fn=self._log, read_only=ro)
            key_map = server.build_key_map()
            mode = "read-only" if ro else "read-write"
            self._log(f"Server connected ({mode}), {len(key_map)} DAQ channels")
        except Exception as e:
            self._log(f"Server error: {e}", level="error"); return

        try:
            hv = HVClient(self._ge_hv_edit.text().strip(),
                          log_fn=self._log, read_only=ro)
            hv.connect(password=self._ge_hv_pw.text())
        except Exception as e:
            self._log(f"HV error: {e}", level="error")
            if not ro:
                QMessageBox.critical(self, "HV Connection Failed", str(e))
            return

        eng = GainScanEngine(
            motor_ep=self.ep, server=server, hv=hv,
            modules=self.scan_modules, log_fn=self._log, key_map=key_map)
        # use profile order if a named path is selected
        if self._use_profile_order:
            eng.path = list(self._ordered_path)
        eng.target_adc = self._ge_target.value()
        eng.min_counts = self._ge_counts.value()
        eng.max_iterations = self._ge_maxiter.value()
        eng.convergence_tol = self._ge_tol.value()
        eng.beam_threshold = self._beam_thresh_spin.value()
        eng.pos_threshold = self._thresh_spin.value()
        self._gain_engine = eng
        eng.start(self._selected_start_idx, count=self._count_spin.value())

    def _cmdPause(self):
        eng = self._gain_engine
        if not eng: return
        if eng._paused:
            eng.resume(); self._btn_pause.setText("Pause")
        else:
            eng.pause(); self._btn_pause.setText("Resume")

    def _cmdStop(self):
        if self._gain_engine:
            self._gain_engine.stop()
            self._btn_pause.setText("Pause")
        else:
            epics_stop(self.ep); self._log("Motors stopped")

    def _cmdSkip(self):
        if self._gain_engine:
            self._gain_engine.skip_module()

    def _setTarget(self, px, py, name=""):
        self._target_px = px; self._target_py = py; self._target_name = name

    def _cmdMoveToModule(self):
        self._onStartSelected(0)
        names = [m.name for m in (self._gain_engine.path if self._gain_engine else [])]
        if not names: return
        path = self._ordered_path
        if self._selected_start_idx >= len(path): return
        mod = path[self._selected_start_idx]
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if epics_move_to(self.ep, px, py):
            self._setTarget(px, py, mod.name)
        else:
            self._log("BLOCKED: outside limits", level="error")

    def _cmdResetCenter(self):
        self._log(f"Resetting to beam centre ptrans({BEAM_CENTER_X}, {BEAM_CENTER_Y})")
        if epics_move_to(self.ep, BEAM_CENTER_X, BEAM_CENTER_Y):
            self._setTarget(BEAM_CENTER_X, BEAM_CENTER_Y, "Beam Center")

    # -- path management (shared logic) -------------------------------------

    def _onStartSelected(self, _):
        name = self._start_combo.currentText()
        path = self._ordered_path
        for i, m in enumerate(path):
            if m.name == name:
                self._selected_start_idx = i; self._drawPathPreview(); break

    def _onPathProfileChanged(self, _):
        name = self._profile_combo.currentText()
        if name == self._active_profile: return
        self._active_profile = name
        if name == self.AUTOGEN:
            self._lg_spin.setEnabled(True); self._onLgLayersChanged(force=True); return
        self._lg_spin.setEnabled(False)
        if name == self.NONE:
            self.scan_modules = []; self._scan_names = set()
            self._scan_name_to_idx = {}; self._selected_start_idx = 0
            self._start_combo.clear(); self._count_spin.setMaximum(0); self._count_spin.setValue(0)
            self._map.setPathPreview([]); self._map.setDashPreview([])
            self._map.setHighlight(None); self._selected_mod_name = None
            self._updateCanvasLabel(); self._log("Path: none"); return
        mod_by_name = {m.name: m for m in self.all_modules}
        path_mods = [mod_by_name[n] for n in self._profiles.get(name, []) if n in mod_by_name]
        if not path_mods:
            self._log(f"Profile '{name}' empty", level="error"); return
        self._setPath(path_mods, use_profile_order=True)
        self._log(f"Path profile: {name} ({len(path_mods)} modules)")

    def _onLgLayersChanged(self, value=0, force=False):
        if self._active_profile != self.AUTOGEN: return
        nl = self._lg_spin.value()
        if nl == self._lg_layers and not force: return
        self._lg_layers = nl
        mods = filter_scan_modules(self.all_modules, nl, self._lg_sx, self._lg_sy)
        self._setPath(mods)
        np_ = sum(1 for m in mods if m.mod_type == "PbWO4")
        ng = sum(1 for m in mods if m.mod_type == "PbGlass")
        self._log(f"LG layers: {nl} ({np_} PbWO4 + {ng} PbGlass = {len(mods)})")

    def _setPath(self, mods, use_profile_order=False):
        self.scan_modules = mods
        self._scan_names = {m.name for m in mods}
        self._use_profile_order = use_profile_order
        if use_profile_order:
            path = mods  # preserve order from paths.json
        else:
            path, _ = build_scan_path(mods)
        self._ordered_path = path
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(path)}
        self._selected_start_idx = 0
        ns = [m.name for m in path]
        self._start_combo.clear(); self._start_combo.addItems(ns)
        self._count_spin.setMaximum(len(ns)); self._count_spin.setValue(0)
        self._updateCanvasLabel()

    # -- canvas -------------------------------------------------------------

    def _updateCanvasLabel(self):
        n_pwo4 = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        base = f"Scan Path: {n_pwo4} PbWO4 + {n_lg} LG" if n_lg else f"Scan Path: {n_pwo4} PbWO4"
        if not self.scan_modules: base = "Scan Path: none"
        self._canvas_label.setText(f" {base} ")

    def _updateCanvas(self):
        eng = self._gain_engine
        colors: Dict[str, str] = {}

        # non-scan modules: show type colour only when scalers are off
        if not self._map._scaler_enabled:
            for m in self.all_modules:
                if m.name in self._scan_names or m.mod_type == "LMS": continue
                colors[m.name] = (C.MOD_GLASS if m.mod_type == "PbGlass"
                                  else C.MOD_PWO4_BG if m.mod_type == "PbWO4" else C.MOD_LMS)

        # scan path modules
        path = self._ordered_path
        if eng and eng.state not in (GainScanState.IDLE, GainScanState.COMPLETED):
            for i, mod in enumerate(path):
                if i == eng.current_idx and eng.state in (GainScanState.MOVING,):
                    colors[mod.name] = C.YELLOW
                elif i == eng.current_idx:
                    colors[mod.name] = C.ACCENT
                elif i in eng.converged:
                    colors[mod.name] = C.GREEN
                elif i in eng.failed:
                    colors[mod.name] = C.RED
                else:
                    colors[mod.name] = C.MOD_TODO
        else:
            si = self._selected_start_idx
            count = self._count_spin.value()
            ei = min(si + count, len(path)) if count > 0 else len(path)
            if eng:
                for i, mod in enumerate(path):
                    if i in eng.converged: colors[mod.name] = C.GREEN
                    elif i in eng.failed: colors[mod.name] = C.RED
                    elif i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
                    elif i == si: colors[mod.name] = C.MOD_SELECTED
                    else: colors[mod.name] = C.MOD_TODO
            else:
                for i, mod in enumerate(path):
                    if i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
                    elif i == si: colors[mod.name] = C.MOD_SELECTED
                    else: colors[mod.name] = C.MOD_TODO

        self._map.setModuleColors(colors)
        self._drawPathPreview()
        rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
        self._map.setMarkerPosition(*ptrans_to_module(rx, ry))
        self._map.update()

    def _drawPathPreview(self):
        path = self._ordered_path
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
        self._map.setHighlight(name); self._updateCanvasLabel()

    # -- scaler controls (shared logic) -------------------------------------

    def _toggleScaler(self):
        on = not self._map._scaler_enabled
        self._map.setScalerEnabled(on)
        self._btn_scaler_toggle.setText("Scalers: ON" if on else "Scalers: OFF")
        self._btn_scaler_toggle.setStyleSheet(self._toggle_btn_ss(on))

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
        self._btn_scaler_log.setStyleSheet(self._small_btn_ss(C.ACCENT if on else C.DIM))

    def _cycleScalerPalette(self):
        self._map.cyclePalette(); self._updatePaletteBtn()

    def _updatePaletteBtn(self):
        idx = self._map._palette_idx
        stops = list(PALETTES.values())[idx]
        parts = [f"stop:{t:.2f} rgb({r},{g},{b})" for t, (r, g, b) in stops]
        self._btn_palette.setStyleSheet(
            f"QPushButton{{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,{','.join(parts)});"
            f"border:1px solid #30363d;border-radius:2px;color:#c9d1d9;"
            f"font:bold 11pt Consolas;padding:0 4px;}}"
            f"QPushButton:hover{{border-color:#58a6ff;}}")
        self._btn_palette.setText(PALETTE_NAMES[idx])

    def _pollScalers(self):
        vals = self.scaler_ep.get_all()
        if vals:
            self._map.setScalerValues(vals)
            if self._scaler_auto_on:
                vmin, vmax = self._map.scalerRange()
                self._scaler_min_edit.setText(f"{vmin:.0f}")
                self._scaler_max_edit.setText(f"{vmax:.0f}")

    # -- logging ------------------------------------------------------------

    def _log(self, msg, level="info"):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {level.upper().ljust(5)} {msg}"
        self._log_lines.append(line)
        if self._log_file and not self._log_file.closed:
            self._log_file.write(line + "\n"); self._log_file.flush()
        self._logSignal.emit(line, level)

    def _appendLog(self, line, level):
        colors = {"info": C.TEXT, "warn": C.YELLOW, "error": C.RED}
        c = colors.get(level, C.DIM)
        self._log_text.append(
            f'<span style="color:{c};font-family:Consolas;font-size:13pt;">'
            f'{html_mod.escape(line)}</span>')

    # -- polling (5 Hz) -----------------------------------------------------

    def _poll(self):
        self._updateGainStatus()
        self._updatePositionCheck()
        self._updateCanvas()
        self._updateBeamDisplay()
        self._checkEncoder()

    def _updateGainStatus(self):
        eng = self._gain_engine
        if eng is None: return
        running = eng.state not in (GainScanState.IDLE, GainScanState.COMPLETED)
        self._path_group.setVisible(not running)
        self._btn_start.setEnabled(not running)
        for w in (self._ge_server_edit, self._ge_hv_edit, self._ge_hv_pw,
                  self._ge_target, self._ge_counts, self._ge_maxiter, self._ge_tol):
            w.setEnabled(not running)
        self._btn_pause.setEnabled(running)
        self._btn_stop.setEnabled(running)
        self._btn_skip.setEnabled(running)
        sc = {GainScanState.IDLE: C.DIM, GainScanState.MOVING: C.YELLOW,
              GainScanState.COLLECTING: C.ACCENT, GainScanState.ANALYZING: C.ACCENT,
              GainScanState.ADJUSTING: C.ORANGE, GainScanState.CONVERGED: C.GREEN,
              GainScanState.FAILED: C.RED, GainScanState.COMPLETED: C.GREEN}
        self._lbl_state.setText(eng.state)
        self._lbl_state.setStyleSheet(f"color:{sc.get(eng.state, C.DIM)};font:bold 15pt 'Consolas';background:transparent;")
        self._lbl_ge_status.setText(eng.state)
        self._lbl_ge_status.setStyleSheet(f"color:{sc.get(eng.state, C.DIM)};")

        done = len(eng.converged) + len(eng.failed)
        total = getattr(eng, '_end_idx', len(eng.path)) - getattr(eng, '_start_idx', 0)
        self._lbl_progress.setText(f"Progress: {done}/{total}")
        self._progress_bar.setMaximum(max(total, 1)); self._progress_bar.setValue(done)

        mod = eng.current_module
        parts = []
        if mod: parts.append(mod.name)
        if eng.current_iteration > 0:
            parts.append(f"iter {eng.current_iteration}/{eng.max_iterations}")
        if eng.last_edge_adc is not None:
            parts.append(f"edge={eng.last_edge_adc:.0f}")
        if eng.last_dv is not None:
            parts.append(f"ΔV={eng.last_dv:+.0f}")
        if eng.state == GainScanState.COLLECTING:
            parts.append(f"counts={eng.module_counts}")
            if eng.collect_rate > 0:
                parts.append(f"{eng.collect_rate:.0f} Hz")
        parts.append(f"[{len(eng.converged)}ok {len(eng.failed)}fail]")
        self._lbl_ge_detail.setText("  ".join(parts))

        if eng.state == GainScanState.MOVING and mod:
            px, py = module_to_ptrans(mod.x, mod.y)
            self._setTarget(px, py, mod.name)

        # update histogram display
        mod_name = mod.name if mod else ""
        target_bin = int((eng.target_adc - eng.analyzer.bin_min) / eng.analyzer.bin_step) \
            if eng.analyzer.bin_step > 0 else None
        self._histogram.setTitle(mod_name)
        info_parts = []
        if eng.last_vset is not None:
            info_parts.append(f"V={eng.last_vset:.1f}")
        if eng.last_edge_adc is not None:
            info_parts.append(f"edge={eng.last_edge_adc:.0f}")
        if eng.last_dv is not None:
            info_parts.append(f"ΔV={eng.last_dv:+.0f}")
        if eng.state == GainScanState.COLLECTING and eng.collect_rate > 0:
            info_parts.append(f"{eng.collect_rate:.0f} Hz")
        self._histogram.setInfo("  ".join(info_parts))

        # fetch live histogram during collection for preview (~every 2s)
        if eng.state == GainScanState.COLLECTING and mod:
            import time as _time
            now = _time.time()
            if now - getattr(self, '_last_hist_fetch', 0) > 2.0:
                self._last_hist_fetch = now
                key = eng.key_map.get(mod.name)
                if key and eng.module_counts > 0:
                    try:
                        hist = eng.server.get_height_histogram(key, quiet=True)
                        live_bins = hist.get("bins", [])
                        if live_bins:
                            self._histogram.setData(live_bins, target_bin, None)
                    except Exception:
                        pass
        elif eng.last_bins:
            self._histogram.setData(eng.last_bins, target_bin, eng.last_edge_bin)

    def _updatePositionCheck(self):
        # position check
        rx, ry = self.ep.get("x_rbv", 0.0), self.ep.get("y_rbv", 0.0)
        self._lbl_actual.setText(f"Actual: ({rx:.3f}, {ry:.3f})")
        px, py = self._target_px, self._target_py
        if px is not None and py is not None:
            err = math.sqrt((rx - px)**2 + (ry - py)**2)
            name_html = f' <b style="color:{C.ACCENT}">{self._target_name}</b>' if self._target_name else ""
            self._lbl_expected.setText(f"Target: ({px:.3f}, {py:.3f}){name_html}")
            vx = self.ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
            vy = self.ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
            dx, dy = abs(rx - px), abs(ry - py)
            eta = max(dx / vx if vx > 0 else 0, dy / vy if vy > 0 else 0)
            eta_str = f" ({int(eta)//60}m {int(eta)%60}s)" if eta >= 60 else (f" ({eta:.0f}s)" if eta >= 1 else "")
            self._lbl_error.setText(f"Diff:   {err:.3f} mm{eta_str}")
        else:
            self._lbl_expected.setText("Target: --")
            self._lbl_error.setText("Diff:   --")

    def _updateBeamDisplay(self):
        bc = self.ep.get("beam_cur", None)
        if bc is None:
            self._lbl_beam_val.setText("-- nA"); return
        thresh = self._beam_thresh_spin.value()
        if thresh > 0 and bc < thresh:
            fg = C.RED if (self._gain_engine and self._gain_engine._paused) else C.YELLOW
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color:{fg};font:bold 18pt 'Consolas';background:transparent;border:none;")
        else:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color:{C.GREEN};font:bold 18pt 'Consolas';background:transparent;border:none;")

    ENCODER_DRIFT_WARN = 0.5
    ENCODER_DRIFT_ERR  = 1.5

    def _checkEncoder(self):
        enc_x = self.ep.get("x_encoder", None)
        enc_y = self.ep.get("y_encoder", None)
        rbv_x = self.ep.get("x_rbv", None)
        rbv_y = self.ep.get("y_rbv", None)
        if enc_x is None or enc_y is None or rbv_x is None or rbv_y is None: return
        if self._enc_offset_x is None:
            self._enc_offset_x = enc_x - rbv_x
            self._enc_offset_y = enc_y - rbv_y
            self._log(f"Encoder calibrated: offset X={self._enc_offset_x:.4f} Y={self._enc_offset_y:.4f}")
            return
        dx = abs((enc_x - self._enc_offset_x) - rbv_x)
        dy = abs((enc_y - self._enc_offset_y) - rbv_y)
        fx = C.RED if dx > self.ENCODER_DRIFT_ERR else (C.YELLOW if dx > self.ENCODER_DRIFT_WARN else C.GREEN)
        fy = C.RED if dy > self.ENCODER_DRIFT_ERR else (C.YELLOW if dy > self.ENCODER_DRIFT_WARN else C.GREEN)
        self._lbl_drift.setText(
            f'Drift:   X <span style="color:{fx}">{dx:.4f}</span>  '
            f'Y <span style="color:{fy}">{dy:.4f}</span>')

    def eventFilter(self, obj, event):
        if obj is self._map and event.type() == event.Type.Resize:
            btn = self._btn_reset_view
            btn.move(self._map.width() - btn.width() - 2,
                     self._map.height() - btn.height() - 2)
        return super().eventFilter(obj, event)

    def closeEvent(self, e):
        self._timer.stop()
        self._scaler_timer.stop()
        if self._gain_engine:
            self._gain_engine.stop()
            t = getattr(self._gain_engine, '_thread', None)
            if t and t.is_alive():
                t.join(timeout=2.0)
            # close HV WebSocket to unblock reader thread
            if hasattr(self._gain_engine, 'hv'):
                self._gain_engine.hv.close()
        if self._log_file and not self._log_file.closed:
            self._log_file.close()
        self._log_file = None
        super().closeEvent(e)
        # force exit in case daemon threads are still blocked
        import os
        os._exit(0)


# ============================================================================
#  MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="HyCal Gain Equalizer")
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

    profiles: Dict[str, List[str]] = {}
    if os.path.exists(args.paths):
        with open(args.paths) as f:
            profiles = json.load(f)
        print(f"Loaded {len(profiles)} path profiles")

    observer = args.observer
    simulation = not args.expert and not observer

    if observer:        motor_ep = ObserverEPICS()
    elif simulation:    motor_ep = SimulatedMotorEPICS()
    else:               motor_ep = MotorEPICS(writable=True)

    n_ok, n_total = motor_ep.connect()
    if not simulation:
        print(f"EPICS: {n_ok}/{n_total} PVs connected")

    if simulation:  scaler_ep = SimulatedScalerEPICS(all_modules)
    else:           scaler_ep = ScalerPVGroup(all_modules)
    scaler_ep.connect()

    app = QApplication(sys.argv)
    win = GainEqualizerWindow(motor_ep, scaler_ep, simulation, all_modules,
                              profiles, observer=observer)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
