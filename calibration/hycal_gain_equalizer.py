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
import os
import sys
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
    SPMG, SPMG_LABELS, epics_move_to, epics_stop,
)
from scan_engine import (
    build_scan_path,
    DEFAULT_POS_THRESHOLD, DEFAULT_BEAM_THRESHOLD, MAX_LG_LAYERS,
)
from scan_geoview import HyCalScanMapWidget, PALETTES, PALETTE_NAMES
from gain_scanner import (
    GainScanEngine, GainScanState, ServerClient, HVClient,
)
from scan_gui_common import (
    PATHS_FILE, SCALER_POLL_MS, POLL_MS,
    open_session_log, format_log_line, append_log_line,
    build_position_check_panel, update_position_check, EncoderDriftChecker,
    load_profiles, setup_motor_epics, setup_scaler_epics,
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
        self._target_bin: Optional[int] = None
        self._edge_bin: Optional[int] = None
        self._bin_min: float = 0.0      # ADC offset of bin 0
        self._bin_step: float = 1.0     # ADC width per bin
        self._title: str = ""
        self._info: str = ""
        self._log_y: bool = True  # log or linear y scale
        self.setMinimumHeight(140)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def setLogY(self, on: bool):
        self._log_y = on; self.update()

    def setBinning(self, bin_min: float, bin_step: float):
        """Set the ADC mapping (bin index → ADC value).

        Used by the x-axis tick labels.  Safe to call once when the
        engine is created — the analyzer's binning is fixed for the
        lifetime of a scan.
        """
        self._bin_min = float(bin_min)
        self._bin_step = float(bin_step)
        self.update()

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
        use_log = self._log_y
        log_vmax = _math.log10(max(vmax, 1)) if use_log else 0

        def _bar_frac(v):
            if use_log:
                return _math.log10(v) / log_vmax if log_vmax > 0 else 0
            return v / vmax

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

        # bars
        bar_w = pw / n
        p.setPen(Qt.PenStyle.NoPen)
        for i, v in enumerate(bins):
            if v <= 0: continue
            bh = _bar_frac(v) * ph
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

        # y-axis labels + grid
        p.setPen(QColor(C.DIM))
        p.setFont(QFont("Consolas", 9))
        p.drawText(QRectF(0, T - 2, L - 4, 14),
                   Qt.AlignmentFlag.AlignRight, f"{vmax}")
        p.drawText(QRectF(0, T + ph - 7, L - 4, 14),
                   Qt.AlignmentFlag.AlignRight, "1" if use_log else "0")
        if use_log:
            decade = 10
            while decade < vmax:
                frac = _math.log10(decade) / log_vmax
                gy = T + ph - frac * ph
                p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
                p.drawLine(L + 1, int(gy), L + pw, int(gy))
                p.setPen(QColor(C.DIM))
                p.setFont(QFont("Consolas", 8))
                p.drawText(QRectF(0, gy - 7, L - 4, 14),
                           Qt.AlignmentFlag.AlignRight, f"{decade}")
                decade *= 10
        else:
            n_grid = 4
            for gi in range(1, n_grid + 1):
                frac = gi / (n_grid + 1)
                gy = T + ph - frac * ph
                gval = int(vmax * frac)
                p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
                p.drawLine(L + 1, int(gy), L + pw, int(gy))
                p.setPen(QColor(C.DIM))
                p.setFont(QFont("Consolas", 8))
                p.drawText(QRectF(0, gy - 7, L - 4, 14),
                           Qt.AlignmentFlag.AlignRight, f"{gval}")

        # x-axis tick labels in ADC units (from bin_min / bin_step)
        if self._bin_step > 0:
            adc_min = self._bin_min
            adc_max = self._bin_min + n * self._bin_step
            span = adc_max - adc_min
            # target tick count adapts to plot width so labels don't collide.
            # Budget is ~80 px per label so even when "nice" rounding lands
            # on a smaller step (and we end up with more ticks than the
            # target), adjacent 50-px label boxes still don't overlap.
            target_ticks = max(2, int(pw / 80))
            raw = span / target_ticks
            if raw > 0:
                mag = 10 ** _math.floor(_math.log10(raw))
                norm = raw / mag
                if   norm < 1.5: nice = 1
                elif norm < 3.0: nice = 2
                elif norm < 7.0: nice = 5
                else:            nice = 10
                tick_step = nice * mag
                # round adc_min UP to a multiple of tick_step (with a tiny
                # epsilon so values that land just barely above an integer
                # multiple aren't pushed to the next tick)
                first_tick = _math.ceil(adc_min / tick_step - 1e-9) * tick_step
                p.setFont(QFont("Consolas", 8))
                adc = first_tick
                while adc <= adc_max + 1e-6:
                    bin_idx = (adc - adc_min) / self._bin_step
                    if 0 <= bin_idx <= n:
                        x = L + bin_idx * bar_w
                        p.setPen(QPen(QColor("#30363d"), 1))
                        p.drawLine(int(x), T + ph, int(x), T + ph + 3)
                        p.setPen(QColor(C.DIM))
                        label = f"{int(adc)}" if tick_step >= 1 else f"{adc:g}"
                        # keep edge labels inside the widget bounds:
                        # left ticks → left-align, right ticks → right-align,
                        # interior ticks → centered
                        if x < L + 25:
                            rect = QRectF(x - 2, T + ph + 4, 50, 14)
                            align = Qt.AlignmentFlag.AlignLeft
                        elif x > L + pw - 25:
                            rect = QRectF(x - 48, T + ph + 4, 50, 14)
                            align = Qt.AlignmentFlag.AlignRight
                        else:
                            rect = QRectF(x - 25, T + ph + 4, 50, 14)
                            align = Qt.AlignmentFlag.AlignCenter
                        p.drawText(rect, align, label)
                    adc += tick_step

        p.end()


# ============================================================================
#  MAIN WINDOW
# ============================================================================

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
        self._log_file = open_session_log("gain_eq", self.simulation, self.observer)

        self.scan_modules: List[Module] = []
        self._scan_names: set = set()
        self._scan_name_to_idx: Dict[str, int] = {}
        self._selected_start_idx = 0
        self._selected_mod_name: Optional[str] = None
        self._mod_dlg = None
        self._gain_engine: Optional[GainScanEngine] = None

        self._encoder_checker = EncoderDriftChecker()

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
        self._timer.start(POLL_MS)

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
        self.resize(1600, 1000)

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

        # row 1: control panel (scrollable, vertical only)
        ctrl_scroll = QScrollArea()
        ctrl_scroll.setWidgetResizable(True)
        ctrl_scroll.setFrameShape(QFrame.Shape.NoFrame)
        ctrl_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        ctrl_w = QWidget(); ctrl_lo = QVBoxLayout(ctrl_w)
        ctrl_lo.setSpacing(4); ctrl_lo.setContentsMargins(0, 0, 0, 0)
        self._buildPathControl(ctrl_lo)
        self._buildGainControl(ctrl_lo)
        self._buildControlPanel(ctrl_lo)
        self._buildPositionCheck(ctrl_lo)
        ctrl_lo.addStretch()
        ctrl_scroll.setWidget(ctrl_w)
        right_splitter.addWidget(ctrl_scroll)

        # row 2: peak height histogram
        hist_group = QGroupBox("Peak Height Histogram")
        hist_lo = QVBoxLayout(hist_group); hist_lo.setContentsMargins(4, 4, 4, 4)
        self._histogram = HistogramWidget()
        hist_lo.addWidget(self._histogram)
        self._hist_group = hist_group
        hist_group.setVisible(False)
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

        # Row 1: Path + LG layers
        r = QHBoxLayout(); r.addWidget(QLabel("Path:"))
        self._profile_combo = QComboBox()
        self._profile_combo.addItems([self.NONE, self.AUTOGEN] + sorted(self._profiles.keys()))
        self._profile_combo.setCurrentText(self.NONE)
        self._profile_combo.activated.connect(self._onPathProfileChanged)
        r.addWidget(self._profile_combo, stretch=1)
        r.addWidget(QLabel("LG:"))
        self._lg_spin = QSpinBox(); self._lg_spin.setRange(0, MAX_LG_LAYERS)
        self._lg_spin.setMaximumWidth(60)
        self._lg_spin.valueChanged.connect(self._onLgLayersChanged)
        r.addWidget(self._lg_spin); lo.addLayout(r)

        # Row 2: Start + Count
        r = QHBoxLayout(); r.addWidget(QLabel("Start:"))
        self._start_combo = QComboBox(); self._start_combo.setEditable(True)
        self._start_combo.setMinimumWidth(70)
        self._start_combo.setMaximumWidth(100)
        self._start_combo.activated.connect(self._onStartSelected)
        r.addWidget(self._start_combo)
        r.addStretch()
        r.addWidget(QLabel("Count:"))
        self._count_spin = QSpinBox(); self._count_spin.setRange(0, 0)
        self._count_spin.setSpecialValueText("all")
        self._count_spin.setMaximumWidth(110)
        self._count_spin.setMinimumWidth(100)
        self._count_spin.valueChanged.connect(lambda _: self._drawPathPreview())
        r.addWidget(self._count_spin); lo.addLayout(r)

        # Row 3: Thresholds: Pos (mm) + Curr (nA) — Curr right-aligned
        r = QHBoxLayout()
        r.addWidget(QLabel("Thres.  Pos. (mm)"))
        self._thresh_spin = QDoubleSpinBox(); self._thresh_spin.setRange(0.01, 10.0)
        self._thresh_spin.setValue(DEFAULT_POS_THRESHOLD)
        self._thresh_spin.setSingleStep(0.1); self._thresh_spin.setDecimals(2)
        self._thresh_spin.setMaximumWidth(100)
        r.addWidget(self._thresh_spin)
        r.addStretch()
        r.addWidget(QLabel("Curr. (nA)"))
        self._beam_thresh_spin = QDoubleSpinBox(); self._beam_thresh_spin.setRange(0.0, 1000.0)
        self._beam_thresh_spin.setValue(DEFAULT_BEAM_THRESHOLD)
        self._beam_thresh_spin.setSingleStep(0.1); self._beam_thresh_spin.setDecimals(2)
        self._beam_thresh_spin.setSpecialValueText("off")
        self._beam_thresh_spin.setMaximumWidth(100)
        r.addWidget(self._beam_thresh_spin)
        lo.addLayout(r)

        parent.addWidget(pc)

    def _buildGainControl(self, parent):
        ge = QGroupBox("Gain Equalization"); lo = QVBoxLayout(ge)
        self._gain_group = ge

        r = QHBoxLayout(); r.addWidget(QLabel("Server:"))
        self._ge_server_edit = QLineEdit("http://clondaq6:5051")
        r.addWidget(self._ge_server_edit); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("HV:"))
        self._ge_hv_edit = QLineEdit("http://clonpc19:8765")
        r.addWidget(self._ge_hv_edit); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("HV Password:"))
        self._ge_hv_pw = QLineEdit("prad2_admin"); self._ge_hv_pw.setEchoMode(QLineEdit.EchoMode.Password)
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

        r = QHBoxLayout()
        r.addWidget(QLabel("Edge frac %:"))
        self._ge_edge_frac = QDoubleSpinBox(); self._ge_edge_frac.setRange(0.1, 20.0)
        self._ge_edge_frac.setValue(5.0); self._ge_edge_frac.setSingleStep(0.5)
        self._ge_edge_frac.setDecimals(1); r.addWidget(self._ge_edge_frac)
        self._ge_log_y = QPushButton("LogY: ON")
        self._ge_log_y.setCheckable(True); self._ge_log_y.setChecked(True)
        self._ge_log_y.clicked.connect(self._toggleLogY)
        self._updateLogYBtnStyle(True)
        r.addWidget(self._ge_log_y); lo.addLayout(r)

        parent.addWidget(ge)

    def _buildControlPanel(self, parent):
        cp = QGroupBox("Control"); lo = QVBoxLayout(cp)

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
        self._btn_redo = QPushButton("Redo Current")
        self._btn_redo.setProperty("cssClass", "accent")
        self._btn_redo.clicked.connect(self._cmdRedo); bf2.addWidget(self._btn_redo)
        self._btn_skip = QPushButton("Skip Current")
        self._btn_skip.clicked.connect(self._cmdSkip); bf2.addWidget(self._btn_skip)
        lo.addLayout(bf2)

        self._lbl_progress = QLabel("Progress: --/--"); lo.addWidget(self._lbl_progress)
        self._progress_bar = QProgressBar(); lo.addWidget(self._progress_bar)
        self._lbl_ge_status = QLabel("Idle")
        self._lbl_ge_status.setStyleSheet(f"color:{C.DIM};"); lo.addWidget(self._lbl_ge_status)
        self._lbl_ge_detail = QLabel("")
        self._lbl_ge_detail.setStyleSheet(f"color:{C.DIM};"); lo.addWidget(self._lbl_ge_detail)

        parent.addWidget(cp)

    def _updateLogYBtnStyle(self, on: bool):
        if on:
            self._ge_log_y.setStyleSheet(
                "QPushButton{background:#1f6feb;color:white;"
                "border:1px solid #388bfd;border-radius:3px;padding:5px 12px;"
                "font:bold 13pt 'Consolas';}"
                "QPushButton:hover{background:#388bfd;}")
        else:
            self._ge_log_y.setStyleSheet("")

    def _buildPositionCheck(self, parent):
        labels = build_position_check_panel(parent)
        self._lbl_expected = labels["target"]
        self._lbl_actual = labels["actual"]
        self._lbl_error = labels["diff"]
        self._lbl_drift = labels["drift"]
        self._pos_labels = labels

    def _disableControls(self):
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_redo, self._btn_skip):
            w.setEnabled(False)

    # -- commands -----------------------------------------------------------

    def _toggleLogY(self):
        on = self._ge_log_y.isChecked()
        self._ge_log_y.setText("LogY: ON" if on else "LogY: OFF")
        self._updateLogYBtnStyle(on)
        self._histogram.setLogY(on)

    def _cmdStart(self):
        if self._gain_engine and self._gain_engine.state not in (
                GainScanState.IDLE, GainScanState.COMPLETED, GainScanState.FAILED):
            return

        # resume from failure or stop — reuse existing engine, retry from current module
        eng = self._gain_engine
        if eng and eng._has_run and eng.state in (GainScanState.FAILED, GainScanState.IDLE):
            self._clearHistogramData()
            eng.start()
            return

        if not self.scan_modules:
            # no path loaded — fall back to single-module mode if the user
            # has clicked a scannable module on the canvas
            if not self._selected_mod_name:
                self._log("Select a path, or click a module on the canvas "
                          "to equalize a single module", level="error")
                return
            mod = self._mod_by_name.get(self._selected_mod_name)
            if mod is None or mod.mod_type not in ("PbWO4", "PbGlass"):
                self._log(f"Module {self._selected_mod_name!r} is not scannable "
                          f"(needs to be PbWO4 or PbGlass)", level="error")
                return
            self._setPath([mod])
            self._log(f"Single-module start: {mod.name}")
        # sync start index from combo selection
        self._onStartSelected(0)

        ro = self.simulation
        server_url = self._ge_server_edit.text().strip()
        hv_url = self._ge_hv_edit.text().strip()
        hv_pw = self._ge_hv_pw.text()

        if not ro and not hv_pw.strip():
            self._log("Expert mode requires HV password", level="error")
            QMessageBox.warning(self, "HV Password Required",
                                "Enter the prad2hvd password before starting in expert mode.")
            return

        # pre-flight: verify server is reachable and build key map
        try:
            test_server = ServerClient(server_url, log_fn=self._log, read_only=ro)
            key_map = test_server.build_key_map()
            mode = "read-only" if ro else "read-write"
            self._log(f"Server OK ({mode}), {len(key_map)} DAQ channels")
        except Exception as e:
            self._log(f"Server error: {e}", level="error"); return

        # pre-flight: verify HV is reachable (connections are per-module during scan)
        if not ro:
            try:
                test_hv = HVClient(hv_url, log_fn=self._log, read_only=False)
                test_hv.connect(password=hv_pw)
                test_hv.close()
                self._log("HV pre-flight OK")
            except Exception as e:
                self._log(f"HV error: {e}", level="error")
                QMessageBox.critical(self, "HV Connection Failed", str(e))
                return

        eng = GainScanEngine(
            motor_ep=self.ep,
            server_url=server_url, hv_url=hv_url,
            hv_password=hv_pw, read_only=ro,
            modules=self.scan_modules, log_fn=self._log, key_map=key_map,
            report_prefix="SIM_" if self.simulation else "")
        eng.target_adc = self._ge_target.value()
        eng.min_counts = self._ge_counts.value()
        eng.max_iterations = self._ge_maxiter.value()
        eng.convergence_tol = self._ge_tol.value()
        eng.beam_threshold = self._beam_thresh_spin.value()
        eng.pos_threshold = self._thresh_spin.value()
        eng.analyzer.edge_fraction = self._ge_edge_frac.value() / 100.0
        eng.analyzer.use_log_cumul = self._ge_log_y.isChecked()
        eng.use_log_y = self._ge_log_y.isChecked()
        self._gain_engine = eng
        # tell the histogram widget how to map bin index → ADC value for x-ticks
        self._histogram.setBinning(eng.analyzer.bin_min, eng.analyzer.bin_step)
        eng.start(self._selected_start_idx, count=self._count_spin.value())

    def _cmdPause(self):
        eng = self._gain_engine
        if not eng: return
        if eng._paused:
            eng.resume(); self._btn_pause.setText("Pause")
        else:
            eng.pause(); self._btn_pause.setText("Resume")

    def _cmdStop(self):
        eng = self._gain_engine
        if not eng:
            epics_stop(self.ep); self._log("Motors stopped")
            return
        # if currently running → stop the scan (will be in resumable state)
        running = eng.state not in (
            GainScanState.IDLE, GainScanState.COMPLETED, GainScanState.FAILED)
        if running:
            eng.stop()
            self._btn_pause.setText("Pause")
            return
        # already stopped (FAILED / IDLE-with-has_run / COMPLETED) →
        # discard engine, reset to clean state for a fresh start
        self._gain_engine = None
        self._log("Reset to fresh start")

    def _cmdRedo(self):
        if self._gain_engine:
            self._gain_engine.redo_module()
            self._clearHistogramData()

    def _cmdSkip(self):
        if self._gain_engine:
            self._gain_engine.skip_module()
            self._clearHistogramData()

    def _clearHistogramData(self):
        """Drop any visible histogram bars (target just changed).

        Keeps the title/info banners; only the bars are wiped so the
        operator doesn't see prior-module data while the next collection
        spins up.  Resets the live-fetch timer so the next poll fetches
        immediately once new counts start arriving.
        """
        self._histogram.setData([], None, None)
        self._last_hist_fetch = 0

    def _setTarget(self, px, py, name=""):
        self._target_px = px; self._target_py = py; self._target_name = name

    def _cmdMoveToModule(self):
        self._onStartSelected(0)
        names = [m.name for m in (self._gain_engine.path if self._gain_engine else [])]
        if not names: return
        path = self.scan_modules
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
        path = self.scan_modules
        for i, m in enumerate(path):
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
            self._log("Path: none")
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
        path, _ = build_scan_path(mods)
        self._setPath(path)
        np_ = sum(1 for m in path if m.mod_type == "PbWO4")
        ng = sum(1 for m in path if m.mod_type == "PbGlass")
        self._log(f"LG layers: {nl} ({np_} PbWO4 + {ng} PbGlass = {len(path)})")

    def _setPath(self, path):
        """Set the scan path. ``path`` is the final ordered list of modules."""
        self.scan_modules = path
        self._scan_names = {m.name for m in path}
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(path)}
        self._selected_start_idx = 0
        ns = [m.name for m in path]
        self._start_combo.clear(); self._start_combo.addItems(ns)
        self._count_spin.setMaximum(len(ns)); self._count_spin.setValue(0)
        if not path:
            self._map.setPathPreview([]); self._map.setDashPreview([])
            self._map.setModuleColors({}); self._map.update()
            self._map.setHighlight(None); self._selected_mod_name = None
        self._updateCanvasLabel()

    # -- canvas -------------------------------------------------------------

    def _updateCanvasLabel(self):
        n_pwo4 = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        base = f"Scan Path: {n_pwo4} PbWO4 + {n_lg} LG" if n_lg else f"Scan Path: {n_pwo4} PbWO4"
        if not self.scan_modules:
            base = "Scan Path: none"
        elif self.scan_modules and 0 <= self._selected_start_idx < len(self.scan_modules):
            start_name = self.scan_modules[self._selected_start_idx].name
            base += f"  start: {start_name}"
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
        path = self.scan_modules
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
        # while running, show a dashed preview of modules still ahead;
        # idle / completed / failed → solid preview of the planned path.
        if eng and eng.state not in (GainScanState.IDLE,
                                     GainScanState.COMPLETED,
                                     GainScanState.FAILED):
            ei = getattr(eng, '_end_idx', len(eng.path))
            self._map.setPathPreview([])
            self._map.setDashPreview([
                self._map.modCenter(eng.path[i])
                for i in range(eng.current_idx + 1, ei)
            ])
        else:
            self._drawPathPreview()
            self._map.setDashPreview([])
        rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
        self._map.setMarkerPosition(*ptrans_to_module(rx, ry))
        self._map.update()

    def _drawPathPreview(self):
        path = self.scan_modules
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
        line = format_log_line(msg, level)
        self._log_lines.append(line)
        if self._log_file and not self._log_file.closed:
            self._log_file.write(line + "\n"); self._log_file.flush()
        self._logSignal.emit(line, level)

    def _appendLog(self, line, level):
        append_log_line(self._log_text, line, level)

    # -- polling (5 Hz) -----------------------------------------------------

    def _poll(self):
        self._updateGainStatus()
        self._updatePositionCheck()
        self._updateCanvas()
        self._updateBeamDisplay()
        self._checkEncoder()

    def _updateGainStatus(self):
        eng = self._gain_engine
        if eng is None:
            # clean state: no engine
            self._path_group.setVisible(True)
            self._gain_group.setVisible(True)
            self._hist_group.setVisible(False)
            self._btn_start.setText("Start")
            self._btn_start.setEnabled(True)
            self._btn_pause.setEnabled(False)
            self._btn_stop.setText("Stop")
            self._btn_stop.setEnabled(False)
            self._btn_redo.setEnabled(False)
            self._btn_redo.setText("Redo Current")
            self._btn_skip.setEnabled(False)
            self._btn_skip.setText("Skip Current")
            return

        running = eng.state not in (
            GainScanState.IDLE, GainScanState.COMPLETED, GainScanState.FAILED)
        # "resumable": scan was started, then stopped/failed mid-way
        resumable = (not running) and eng._has_run and \
                    eng.state != GainScanState.COMPLETED

        # detect transition to a new module → drop stale histogram bars
        # (covers natural progression and the post-skip transition once
        # the engine processes the skip flag)
        prev_idx = getattr(self, '_hist_last_idx', None)
        if running and eng.current_idx != prev_idx:
            self._clearHistogramData()
        self._hist_last_idx = eng.current_idx if running else None

        # panel visibility
        self._path_group.setVisible(not running and not resumable)
        self._gain_group.setVisible(not running and not resumable)
        has_data = bool(eng.last_bins) or bool(eng.iteration_history)
        self._hist_group.setVisible(running or has_data)

        # Start button: "Resume" while resumable, "Start" otherwise
        self._btn_start.setText("Resume" if resumable else "Start")
        self._btn_start.setEnabled(not running)
        # Stop button: "Reset" while resumable, "Stop" while running
        self._btn_stop.setText("Reset" if resumable else "Stop")
        self._btn_stop.setEnabled(running or resumable)
        self._btn_pause.setEnabled(running)
        self._btn_redo.setEnabled(running)
        self._btn_skip.setEnabled(running)

        # show current module name on the per-module action buttons
        cur = eng.current_module
        cur_name = cur.name if (running and cur is not None) else ""
        self._btn_redo.setText(f"Redo Current ({cur_name})" if cur_name else "Redo Current")
        self._btn_skip.setText(f"Skip Current ({cur_name})" if cur_name else "Skip Current")

        # sync combo to current_idx so operator sees resume point
        if resumable and 0 <= eng.current_idx < len(eng.path):
            name = eng.path[eng.current_idx].name
            if self._start_combo.currentText() != name:
                idx = self._start_combo.findText(name)
                if idx >= 0:
                    self._start_combo.setCurrentIndex(idx)
                self._selected_start_idx = eng.current_idx
                self._updateCanvasLabel()
        for w in (self._ge_server_edit, self._ge_hv_edit, self._ge_hv_pw,
                  self._ge_target, self._ge_counts, self._ge_maxiter, self._ge_tol,
                  self._ge_edge_frac, self._ge_log_y):
            w.setEnabled(not running and not resumable)
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

        if mod:
            px, py = module_to_ptrans(mod.x, mod.y)
            if self._target_name != mod.name:
                self._setTarget(px, py, mod.name)

        # update histogram display
        mod_name = mod.name if mod else ""
        target_bin = int((eng.target_adc - eng.analyzer.bin_min) / eng.analyzer.bin_step) \
            if eng.analyzer.bin_step > 0 else None
        # title: module name + HV readings
        title_parts = [mod_name]
        if eng.last_vmon is not None:
            title_parts.append(f"VMon={eng.last_vmon:.1f}")
        if eng.last_vset is not None:
            title_parts.append(f"VSet={eng.last_vset:.1f}")
        self._histogram.setTitle("  ".join(title_parts))
        info_parts = []
        if eng.last_edge_adc is not None:
            info_parts.append(f"edge={eng.last_edge_adc:.0f}")
        if eng.last_dv is not None:
            info_parts.append(f"ΔV={eng.last_dv:+.0f}")
        if eng.state == GainScanState.COLLECTING and eng.collect_rate > 0:
            info_parts.append(f"{eng.collect_rate:.0f} Hz")
        self._histogram.setInfo("  ".join(info_parts))

        # fetch live histogram during collection for preview (~every 2s)
        if eng.state == GainScanState.COLLECTING and mod and eng.server:
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
        update_position_check(self._pos_labels, self.ep,
                              self._target_px, self._target_py,
                              self._target_name)

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

    def _checkEncoder(self):
        self._encoder_checker.update(self.ep, self._log, self._lbl_drift)

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
        if self._log_file and not self._log_file.closed:
            self._log_file.close()
        self._log_file = None
        super().closeEvent(e)


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

    profiles = load_profiles(args.paths)
    if profiles:
        print(f"Loaded {len(profiles)} path profiles")

    observer = args.observer
    simulation = not args.expert and not observer

    motor_ep = setup_motor_epics(observer, simulation)
    scaler_ep = setup_scaler_epics(simulation, all_modules)

    app = QApplication(sys.argv)
    win = GainEqualizerWindow(motor_ep, scaler_ep, simulation, all_modules,
                              profiles, observer=observer)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
