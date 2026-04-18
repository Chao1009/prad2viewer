#!/usr/bin/env python3
"""
GEM Event Viewer (PyQt6) — step through an EVIO file event-by-event,
running full GEM reconstruction (pedestal + CM + ZS + clustering) in
process via ``prad2py``.

Features:
  * Pre-scans the file on open to build an event index (progress dialog).
  * Prev / Next / Goto event# + slider for navigation.
  * Primary threshold sliders (ZS, CM, min-cluster-hits) above the canvas.
  * Collapsible "Advanced tuning" dock for every other knob on
    GemSystem / GemCluster — values are live; each change re-runs
    reconstruction on cached SSP data (no EVIO I/O).

Usage:
    python gem/gem_event_viewer.py [file.evio.00000]

If an EVIO path is given on the command line the viewer starts scanning
it immediately; otherwise use File → Open EVIO.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# prad2py auto-discovery — walk up from this script to find build/python/.
# ---------------------------------------------------------------------------

_SCRIPT_DIR = Path(__file__).resolve().parent
_probe = _SCRIPT_DIR
for _ in range(5):
    _probe = _probe.parent
    for _sub in ("build/python", "build-release/python", "build/Release/python"):
        _cand = _probe / _sub
        if _cand.is_dir() and str(_cand) not in sys.path:
            sys.path.insert(0, str(_cand))

try:
    import prad2py
    from prad2py import dec, det
    HAVE_PRAD2PY = True
    PRAD2PY_ERROR = ""
except Exception as _exc:  # noqa: BLE001
    prad2py = None  # type: ignore
    dec = None  # type: ignore
    det = None  # type: ignore
    HAVE_PRAD2PY = False
    PRAD2PY_ERROR = f"{type(_exc).__name__}: {_exc}"


from PyQt6.QtCore import (  # noqa: E402
    QObject, Qt, QThread, QTimer, pyqtSignal,
)
from PyQt6.QtGui import QAction, QKeySequence  # noqa: E402
from PyQt6.QtWidgets import (  # noqa: E402
    QApplication,
    QComboBox,
    QDockWidget,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QSlider,
    QSpinBox,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

# Sibling imports — this file lives in gem/ alongside the helpers.
# These modules in turn import gem_strip_map which requires prad2py.det, so
# wrap the whole block so the GUI can still *start* and show an error
# dialog when the module is missing.
try:
    from gem_layout import build_strip_layout, load_gem_map  # noqa: E402
    from gem_view import (  # noqa: E402
        build_apv_map,
        build_det_list_from_gemsys,
        build_zs_apvs_from_gemsys,
        draw_event,
        process_zs_hits,
    )
except Exception as _sib_exc:  # noqa: BLE001
    build_strip_layout = load_gem_map = None  # type: ignore
    build_apv_map = build_det_list_from_gemsys = build_zs_apvs_from_gemsys = None  # type: ignore
    draw_event = process_zs_hits = None  # type: ignore
    if HAVE_PRAD2PY:
        PRAD2PY_ERROR = (PRAD2PY_ERROR + "\n" if PRAD2PY_ERROR else "") + \
                        f"sibling import: {type(_sib_exc).__name__}: {_sib_exc}"
    HAVE_PRAD2PY = False


# ---------------------------------------------------------------------------
# Event index
# ---------------------------------------------------------------------------


@dataclass
class EventMeta:
    """Metadata for one physics sub-event in the source EVIO file.

    ``record_idx`` / ``subevt_idx`` locate the event in the file's
    record structure — needed by Stepper to re-read it.  Identity
    fields (event_number / trigger_*) come from the EventInfo bank and
    are shown in the UI.
    """
    record_idx: int
    subevt_idx: int
    event_number: int
    trigger_number: int
    trigger_bits: int


class ScanWorker(QObject):
    """Builds an EventMeta list for every Physics sub-event in a file.

    Runs on a QThread so the UI stays responsive.  Emits ``progress``
    every ``PROGRESS_EVERY`` physics events seen and ``finished`` once
    done.  Cancel via ``request_cancel()``.
    """

    PROGRESS_EVERY = 2000

    progress = pyqtSignal(int, int)           # (physics_seen, records_seen)
    finished = pyqtSignal(object, float)      # (List[EventMeta], elapsed_seconds)
    failed   = pyqtSignal(str)

    def __init__(self, path: str, daq_config_path: str):
        super().__init__()
        self._path = path
        self._daq = daq_config_path
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            cfg = dec.load_daq_config(self._daq)
            ch = dec.EvChannel()
            ch.set_config(cfg)
            st = ch.open(self._path)
            if st != dec.Status.success:
                raise RuntimeError(f"cannot open {self._path}: {st}")
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return

        events: List[EventMeta] = []
        record_idx = 0
        start = time.monotonic()
        try:
            while ch.read() == dec.Status.success:
                if self._cancel:
                    break
                if not ch.scan():
                    record_idx += 1
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    record_idx += 1
                    continue
                n_sub = ch.get_n_events()
                for i in range(n_sub):
                    ch.select_event(i)
                    info = ch.info()
                    events.append(EventMeta(
                        record_idx=record_idx,
                        subevt_idx=i,
                        event_number=int(info.event_number),
                        trigger_number=int(info.trigger_number),
                        trigger_bits=int(info.trigger_bits),
                    ))
                    if len(events) % self.PROGRESS_EVERY == 0:
                        self.progress.emit(len(events), record_idx + 1)
                record_idx += 1
        except Exception as exc:  # noqa: BLE001
            ch.close()
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        ch.close()
        elapsed = time.monotonic() - start
        self.progress.emit(len(events), record_idx)
        self.finished.emit(events, elapsed)


# ---------------------------------------------------------------------------
# Stepper — reads one event by index from the EVIO file
# ---------------------------------------------------------------------------


class Stepper:
    """Positioned EVIO reader.  Fetches SspEventData for any event in the
    file, using a cache of recently-visited events for fast Prev/Next.

    EVIO is sequential, so random access to event N sometimes requires a
    re-open + scan-up-to-N.  We only re-open when the requested event is
    behind our current read position; forward jumps are a linear walk.
    """

    CACHE_SIZE = 16

    def __init__(self, path: str, daq_config_path: str):
        self._path = path
        self._daq = daq_config_path
        self._ch: Optional[object] = None
        # Index of the record the reader is currently positioned on
        # (after a successful scan).  -1 before the first Read().
        self._position = -1
        self._cache: Dict[int, object] = {}  # event_idx -> SspEventData
        self._cache_order: List[int] = []

    # --- lifecycle -------------------------------------------------------

    def open(self):
        cfg = dec.load_daq_config(self._daq)
        self._ch = dec.EvChannel()
        self._ch.set_config(cfg)
        st = self._ch.open(self._path)
        if st != dec.Status.success:
            raise RuntimeError(f"cannot open {self._path}: {st}")
        self._position = -1

    def close(self):
        if self._ch is not None:
            self._ch.close()
            self._ch = None
        self._position = -1

    # --- fetch -----------------------------------------------------------

    def get_ssp(self, evmeta: EventMeta, event_idx: int):
        """Return the SspEventData for ``evmeta`` (event_idx is the ordinal
        into the EventMeta list, used as cache key)."""
        if event_idx in self._cache:
            self._cache_order.remove(event_idx)
            self._cache_order.append(event_idx)
            return self._cache[event_idx]

        if self._ch is None:
            self.open()

        target = evmeta.record_idx
        if self._position > target:
            # Target is behind us — re-open and walk forward.
            self.close()
            self.open()

        # Walk forward until the reader is positioned on the target record.
        # (If _position == target we're already there — a second sub-event
        # on the same record.)
        while self._position < target:
            st = self._ch.read()
            if st != dec.Status.success:
                raise RuntimeError(
                    f"EOF before reaching record {target}: {st}")
            self._position += 1
            if not self._ch.scan():
                # Non-physics / unscannable record — keep walking.
                continue

        self._ch.select_event(evmeta.subevt_idx)
        ssp = self._ch.gem()

        # LRU cache insert (keep CACHE_SIZE most recent SSP payloads).
        self._cache[event_idx] = ssp
        self._cache_order.append(event_idx)
        while len(self._cache_order) > self.CACHE_SIZE:
            old = self._cache_order.pop(0)
            self._cache.pop(old, None)

        return ssp


# ---------------------------------------------------------------------------
# Matplotlib canvas
# ---------------------------------------------------------------------------


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self):
        self.fig = Figure(figsize=(10, 4), constrained_layout=True)
        super().__init__(self.fig)


# ---------------------------------------------------------------------------
# Advanced tuning dock
# ---------------------------------------------------------------------------


class AdvancedDock(QDockWidget):
    """Collapsible dock with every clustering / XY-match knob exposed as
    spinboxes.  Emits ``changed`` whenever any value changes."""

    changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__("Advanced tuning", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea
                             | Qt.DockWidgetArea.RightDockWidgetArea)
        self.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )

        root = QWidget()
        self.setWidget(root)
        layout = QVBoxLayout(root)
        layout.setContentsMargins(6, 6, 6, 6)

        # --- Clustering ----------------------------------------------
        cg = QGroupBox("Clustering")
        cf = QFormLayout(cg)
        self.max_cluster_hits = self._mkspin(1, 100, 20)
        self.consecutive_thres = self._mkspin(0, 10, 1)
        self.split_thres = self._mkdspin(0.0, 100.0, 6.0, 0.1)
        self.cross_talk_width = self._mkspin(1, 64, 16)
        cf.addRow("max_cluster_hits", self.max_cluster_hits)
        cf.addRow("consecutive_thres", self.consecutive_thres)
        cf.addRow("split_thres (ADC)", self.split_thres)
        cf.addRow("cross_talk_width", self.cross_talk_width)
        layout.addWidget(cg)

        # --- XY matching ---------------------------------------------
        xg = QGroupBox("XY matching")
        xf = QFormLayout(xg)
        self.match_mode = QComboBox()
        self.match_mode.addItems(["0 — sorted pairing", "1 — cartesian + cuts"])
        self.match_mode.setCurrentIndex(1)
        self.match_adc_asym = self._mkdspin(0.0, 1.0, 0.8, 0.05)
        self.match_time_diff = self._mkdspin(0.0, 200.0, 50.0, 1.0)
        self.ts_period = self._mkdspin(1.0, 100.0, 25.0, 0.5)
        xf.addRow("match_mode", self.match_mode)
        xf.addRow("match_adc_asymmetry", self.match_adc_asym)
        xf.addRow("match_time_diff (ns)", self.match_time_diff)
        xf.addRow("ts_period (ns)", self.ts_period)
        layout.addWidget(xg)

        layout.addStretch(1)

        # Wire every editor → changed
        for w in (self.max_cluster_hits, self.consecutive_thres,
                  self.cross_talk_width):
            w.valueChanged.connect(lambda *_: self.changed.emit())
        for w in (self.split_thres, self.match_adc_asym,
                  self.match_time_diff, self.ts_period):
            w.valueChanged.connect(lambda *_: self.changed.emit())
        self.match_mode.currentIndexChanged.connect(lambda *_: self.changed.emit())

    @staticmethod
    def _mkspin(lo: int, hi: int, val: int) -> QSpinBox:
        sb = QSpinBox(); sb.setRange(lo, hi); sb.setValue(val); return sb

    @staticmethod
    def _mkdspin(lo: float, hi: float, val: float, step: float) -> QDoubleSpinBox:
        sb = QDoubleSpinBox(); sb.setRange(lo, hi)
        sb.setSingleStep(step); sb.setDecimals(3); sb.setValue(val)
        return sb

    def apply_to(self, cluster: "det.GemCluster"):
        """Write current dock values into a GemCluster's ClusterConfig."""
        cfg = cluster.get_config()
        cfg.max_cluster_hits = int(self.max_cluster_hits.value())
        cfg.consecutive_thres = int(self.consecutive_thres.value())
        cfg.split_thres = float(self.split_thres.value())
        cfg.cross_talk_width = int(self.cross_talk_width.value())
        cfg.match_mode = int(self.match_mode.currentIndex())
        cfg.match_adc_asymmetry = float(self.match_adc_asym.value())
        cfg.match_time_diff = float(self.match_time_diff.value())
        cfg.ts_period = float(self.ts_period.value())
        cluster.set_config(cfg)


# ---------------------------------------------------------------------------
# Config discovery helpers
# ---------------------------------------------------------------------------


def _find_first(candidates: List[Path]) -> Optional[Path]:
    for p in candidates:
        if p.is_file():
            return p
    return None


def default_daq_config() -> Optional[Path]:
    root = _SCRIPT_DIR.parent.parent
    return _find_first([
        root / "database" / "daq_config.json",
        Path.cwd() / "database" / "daq_config.json",
        Path.cwd() / "daq_config.json",
    ])


def default_gem_map() -> Optional[Path]:
    root = _SCRIPT_DIR.parent.parent
    return _find_first([
        root / "database" / "gem_map.json",
        Path.cwd() / "database" / "gem_map.json",
        Path.cwd() / "gem_map.json",
    ])


# NOTE: no default_gem_ped() auto-discovery.  Pedestals are per-run
# calibration products and a wrong file is worse than none — we require
# the caller to pick the ped file explicitly via --gem-ped / the File menu.


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------


class GemEventViewer(QMainWindow):
    REDRAW_DEBOUNCE_MS = 120

    def __init__(self,
                 initial_evio: Optional[str] = None,
                 daq_config_path: Optional[str] = None,
                 gem_map_path: Optional[str] = None,
                 gem_ped_path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("GEM Event Viewer")
        self.resize(1400, 820)

        self._daq_config_path = str(daq_config_path or default_daq_config() or "")
        self._gem_map_path    = str(gem_map_path or default_gem_map() or "")
        # Pedestals: no auto-discovery — loaded only if the user passes
        # --gem-ped or picks one via File → Choose gem_ped.json.
        self._gem_ped_path    = str(gem_ped_path or "")

        # Latched after a full-readout event is seen without a ped file;
        # used to show a persistent banner in the status bar and to avoid
        # re-showing the warning dialog on every event.
        self._ped_warning_shown = False

        # Geometry (loaded once per gem_map change)
        self._detectors: Dict[int, dict] = {}
        self._hole: Optional[dict] = None
        self._gem_raw: dict = {}
        self._apv_map: dict = {}

        # GEM reconstruction objects (live)
        self._gsys: Optional[det.GemSystem] = None
        self._gcl: Optional[det.GemCluster] = None

        # Gem-map defaults, captured after Init() so "Reset defaults" can
        # restore whatever the JSON specified rather than wholly-untuned.
        self._default_zs = 5.0
        self._default_cm = 20.0
        self._default_cluster_cfg = None

        # EVIO file state
        self._evio_path = ""
        self._events: List[EventMeta] = []
        self._current = -1
        self._stepper: Optional[Stepper] = None
        self._last_ssp = None  # cached SspEventData for threshold re-runs

        # Worker thread for scanning
        self._scan_thread: Optional[QThread] = None
        self._scan_worker: Optional[ScanWorker] = None
        self._progress: Optional[QProgressDialog] = None

        self._build_ui()
        self._load_geometry_and_gemsys()

        # Debounce timer for slider changes — avoids rebuilding on every tick
        # of a drag.
        self._redraw_timer = QTimer(self)
        self._redraw_timer.setSingleShot(True)
        self._redraw_timer.timeout.connect(self._re_reconstruct_current)

        if initial_evio:
            QTimer.singleShot(50, lambda: self._open_evio(initial_evio))

    # -----------------------------------------------------------------
    # UI construction
    # -----------------------------------------------------------------

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)

        # --- Top row: file + nav ---
        top = QHBoxLayout()
        self.file_label = QLabel("(no file loaded)")
        self.file_label.setStyleSheet("color: #444;")
        top.addWidget(self.file_label, 1)

        self.btn_prev = QPushButton("◀ Prev"); self.btn_prev.setShortcut("Left")
        self.btn_next = QPushButton("Next ▶"); self.btn_next.setShortcut("Right")
        self.btn_prev.clicked.connect(lambda: self._step(-1))
        self.btn_next.clicked.connect(lambda: self._step(+1))
        self.btn_prev.setEnabled(False); self.btn_next.setEnabled(False)
        top.addWidget(self.btn_prev)
        top.addWidget(self.btn_next)

        top.addWidget(QLabel("Goto #"))
        self.goto_spin = QSpinBox()
        self.goto_spin.setRange(0, 0)
        self.goto_spin.setEnabled(False)
        self.goto_spin.editingFinished.connect(self._on_goto)
        top.addWidget(self.goto_spin)

        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, 0)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self._on_slider)
        top.addWidget(self.slider, 2)

        self.btn_advanced = QPushButton("Advanced…")
        self.btn_advanced.setCheckable(True)
        self.btn_advanced.toggled.connect(self._toggle_advanced)
        top.addWidget(self.btn_advanced)

        root.addLayout(top)

        # --- Second row: primary threshold sliders ---
        thr = QHBoxLayout()
        thr.addWidget(QLabel("ZS σ:"))
        self.zs_slider, self.zs_spin = self._mkfloat_slider(2.0, 15.0, 5.0, 0.1)
        thr.addWidget(self.zs_slider); thr.addWidget(self.zs_spin)

        thr.addSpacing(16)
        thr.addWidget(QLabel("CM thr:"))
        self.cm_slider, self.cm_spin = self._mkfloat_slider(5.0, 50.0, 20.0, 0.5)
        thr.addWidget(self.cm_slider); thr.addWidget(self.cm_spin)

        thr.addSpacing(16)
        thr.addWidget(QLabel("min cluster hits:"))
        self.mch_spin = QSpinBox(); self.mch_spin.setRange(1, 10); self.mch_spin.setValue(1)
        thr.addWidget(self.mch_spin)

        thr.addSpacing(16)
        self.btn_reset = QPushButton("Reset defaults")
        self.btn_reset.clicked.connect(self._reset_defaults)
        thr.addWidget(self.btn_reset)

        thr.addStretch(1)
        root.addLayout(thr)

        # --- Canvas + matplotlib toolbar ---
        self.canvas = MplCanvas()
        mpl_toolbar = NavigationToolbar2QT(self.canvas, self)
        root.addWidget(mpl_toolbar)
        root.addWidget(self.canvas, 1)

        # --- Status bar ---
        self.setStatusBar(QStatusBar(self))
        self._status = QLabel("")
        self.statusBar().addPermanentWidget(self._status, 1)
        # Red warning badge — only visible when full-readout data is loaded
        # without a pedestal file.  Right-aligned, fixed width.
        self._ped_badge = QLabel("")
        self._ped_badge.setStyleSheet(
            "color: #b00020; font-weight: bold; padding: 0 8px;")
        self._ped_badge.hide()
        self.statusBar().addPermanentWidget(self._ped_badge)
        self._set_status("Open an EVIO file via File → Open EVIO… (Ctrl+O).")

        # --- Advanced dock (hidden by default) ---
        self.adv_dock = AdvancedDock(self)
        self.adv_dock.hide()
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.adv_dock)
        self.adv_dock.visibilityChanged.connect(
            lambda vis: self.btn_advanced.setChecked(bool(vis)))
        self.adv_dock.changed.connect(self._on_threshold_change)

        # --- Wire primary threshold widgets ---
        self.zs_spin.valueChanged.connect(lambda *_: self._on_threshold_change())
        self.cm_spin.valueChanged.connect(lambda *_: self._on_threshold_change())
        self.mch_spin.valueChanged.connect(lambda *_: self._on_threshold_change())

        # --- Menu bar ---
        self._build_menu()

    def _build_menu(self):
        mb = self.menuBar()
        m_file = mb.addMenu("&File")

        act_open = QAction("&Open EVIO…", self)
        act_open.setShortcut(QKeySequence("Ctrl+O"))
        act_open.triggered.connect(self._pick_evio)
        m_file.addAction(act_open)

        act_map = QAction("Choose &gem_map.json…", self)
        act_map.triggered.connect(self._pick_gem_map)
        m_file.addAction(act_map)

        act_ped = QAction("Choose gem_&ped.json…", self)
        act_ped.triggered.connect(self._pick_gem_ped)
        m_file.addAction(act_ped)

        m_file.addSeparator()
        act_quit = QAction("&Quit", self)
        act_quit.setShortcut(QKeySequence("Ctrl+Q"))
        act_quit.triggered.connect(self.close)
        m_file.addAction(act_quit)

        m_view = mb.addMenu("&View")
        self.act_adv = QAction("Show &Advanced tuning", self, checkable=True)
        self.act_adv.triggered.connect(self._toggle_advanced)
        m_view.addAction(self.act_adv)

    @staticmethod
    def _mkfloat_slider(lo: float, hi: float, val: float, step: float
                        ) -> Tuple[QSlider, QDoubleSpinBox]:
        """Build a (slider, spinbox) pair for a float parameter.

        The slider holds integer ticks of ``step`` resolution; the spinbox
        shows the float value.  They stay synchronized — edits on either
        widget reach both.  Return (slider, spinbox) so the caller can add
        them to a layout in whatever order it likes.
        """
        n_ticks = int(round((hi - lo) / step))
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(0, n_ticks)
        slider.setValue(int(round((val - lo) / step)))
        spin = QDoubleSpinBox()
        spin.setRange(lo, hi); spin.setSingleStep(step); spin.setDecimals(2)
        spin.setValue(val)
        spin.setMaximumWidth(90)

        def s2b(tick: int):
            spin.blockSignals(True)
            spin.setValue(lo + tick * step)
            spin.blockSignals(False)

        def b2s(v: float):
            slider.blockSignals(True)
            slider.setValue(int(round((v - lo) / step)))
            slider.blockSignals(False)

        slider.valueChanged.connect(s2b)
        spin.valueChanged.connect(b2s)
        return slider, spin

    # -----------------------------------------------------------------
    # Geometry + GemSystem init
    # -----------------------------------------------------------------

    def _load_geometry_and_gemsys(self):
        if not HAVE_PRAD2PY:
            self._fatal_prad2py_missing()
            return
        if not self._gem_map_path or not os.path.isfile(self._gem_map_path):
            QMessageBox.warning(
                self, "gem_map.json not found",
                "Could not locate gem_map.json — use File → Choose gem_map.json to set it.")
            return

        try:
            layers, apvs, hole, raw = load_gem_map(self._gem_map_path)
            self._detectors = build_strip_layout(layers, apvs, hole, raw)
            self._apv_map = build_apv_map(apvs)
            self._hole = hole
            self._gem_raw = raw
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Bad gem_map.json", str(exc))
            return

        try:
            self._gsys = det.GemSystem()
            self._gsys.init(self._gem_map_path)
            if self._gem_ped_path and os.path.isfile(self._gem_ped_path):
                self._gsys.load_pedestals(self._gem_ped_path)
            self._gcl = det.GemCluster()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "GemSystem init failed", str(exc))
            self._gsys = self._gcl = None
            return

        # Capture gem-map defaults so the "Reset defaults" button can
        # restore the original configuration — GemSystem thresholds come
        # from the JSON map, GemCluster config is the library default.
        self._default_zs = float(self._gsys.zero_sup_threshold)
        self._default_cm = float(self._gsys.common_mode_threshold)
        self._default_cluster_cfg = det.GemCluster().get_config()

        # Pull current system defaults into the widgets.
        self.zs_spin.blockSignals(True)
        self.cm_spin.blockSignals(True)
        self.zs_spin.setValue(self._default_zs)
        self.cm_spin.setValue(self._default_cm)
        self.zs_spin.blockSignals(False)
        self.cm_spin.blockSignals(False)

        cfg = self._gcl.get_config()
        self.adv_dock.max_cluster_hits.setValue(int(cfg.max_cluster_hits))
        self.adv_dock.consecutive_thres.setValue(int(cfg.consecutive_thres))
        self.adv_dock.split_thres.setValue(float(cfg.split_thres))
        self.adv_dock.cross_talk_width.setValue(int(cfg.cross_talk_width))
        self.adv_dock.match_mode.setCurrentIndex(int(cfg.match_mode))
        self.adv_dock.match_adc_asym.setValue(float(cfg.match_adc_asymmetry))
        self.adv_dock.match_time_diff.setValue(float(cfg.match_time_diff))
        self.adv_dock.ts_period.setValue(float(cfg.ts_period))
        self.mch_spin.setValue(int(cfg.min_cluster_hits))

        ped_status = ("ped: " + os.path.basename(self._gem_ped_path)
                      if self._gem_ped_path and os.path.isfile(self._gem_ped_path)
                      else "no pedestals loaded (required for full-readout data)")
        self._set_status(
            f"GEM system ready — {self._gsys.get_n_detectors()} detectors, "
            f"{self._gsys.get_n_apvs()} APVs, {ped_status}.")

        # Ped-loaded state changed — reset the "already warned" latch so the
        # next full-readout event can warn again if still missing peds.
        self._ped_warning_shown = bool(
            self._gem_ped_path and os.path.isfile(self._gem_ped_path))

    # -----------------------------------------------------------------
    # File picker / pre-scan
    # -----------------------------------------------------------------

    def _pick_evio(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open EVIO file", "",
            "EVIO files (*.evio *.evio.*);;All files (*)")
        if path:
            self._open_evio(path)

    def _pick_gem_map(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose gem_map.json", self._gem_map_path or "",
            "JSON (*.json);;All files (*)")
        if path:
            self._gem_map_path = path
            self._load_geometry_and_gemsys()
            if self._events and self._current >= 0:
                self._show_event(self._current)

    def _pick_gem_ped(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose gem_ped.json", self._gem_ped_path or "",
            "JSON (*.json);;All files (*)")
        if path:
            self._gem_ped_path = path
            self._load_geometry_and_gemsys()
            if self._events and self._current >= 0:
                self._show_event(self._current)

    def _open_evio(self, path: str):
        if not HAVE_PRAD2PY:
            self._fatal_prad2py_missing()
            return
        if self._gsys is None:
            QMessageBox.warning(self, "GEM system not ready",
                                "Load a gem_map.json before opening an EVIO file.")
            return
        if not os.path.isfile(path):
            QMessageBox.warning(self, "File not found", path)
            return

        # Tear down previous run if any.
        self._close_current_run()

        self._evio_path = path
        self.file_label.setText(f"Scanning: {os.path.basename(path)}")
        self._set_status("Pre-scanning EVIO file for event index…")

        # Progress dialog
        self._progress = QProgressDialog(
            f"Scanning {os.path.basename(path)}…", "Cancel", 0, 0, self)
        self._progress.setWindowTitle("Building event index")
        self._progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self._progress.setMinimumDuration(0)
        self._progress.setAutoClose(False)
        self._progress.setAutoReset(False)
        self._progress.canceled.connect(self._cancel_scan)
        self._progress.show()

        # Worker + thread
        self._scan_thread = QThread(self)
        self._scan_worker = ScanWorker(path, self._daq_config_path)
        self._scan_worker.moveToThread(self._scan_thread)
        self._scan_thread.started.connect(self._scan_worker.run)
        self._scan_worker.progress.connect(self._on_scan_progress)
        self._scan_worker.finished.connect(self._on_scan_done)
        self._scan_worker.failed.connect(self._on_scan_failed)
        self._scan_thread.start()

    def _cancel_scan(self):
        if self._scan_worker is not None:
            self._scan_worker.request_cancel()

    def _on_scan_progress(self, phys_seen: int, records_seen: int):
        if self._progress is not None:
            self._progress.setLabelText(
                f"Scanning {os.path.basename(self._evio_path)}…\n"
                f"{phys_seen:,} physics events found in {records_seen:,} records")

    def _on_scan_failed(self, msg: str):
        if self._progress is not None:
            self._progress.close(); self._progress = None
        self._tear_down_worker()
        QMessageBox.critical(self, "Scan failed", msg)
        self.file_label.setText("(no file loaded)")
        self._set_status("Scan failed — see dialog.")

    def _on_scan_done(self, events: List[EventMeta], elapsed: float):
        if self._progress is not None:
            self._progress.close(); self._progress = None
        self._tear_down_worker()

        self._events = events
        if not events:
            self.file_label.setText(
                f"{os.path.basename(self._evio_path)} — no physics events")
            self._set_status("No physics events in file.")
            return

        self._stepper = Stepper(self._evio_path, self._daq_config_path)
        try:
            self._stepper.open()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Reopen failed", str(exc))
            return

        # Configure navigation widgets
        last = len(events) - 1
        ev_nums = [e.event_number for e in events]
        lo_ev, hi_ev = min(ev_nums), max(ev_nums)
        self.slider.blockSignals(True); self.slider.setRange(0, last); self.slider.setValue(0); self.slider.setEnabled(True); self.slider.blockSignals(False)
        self.goto_spin.blockSignals(True); self.goto_spin.setRange(lo_ev, hi_ev); self.goto_spin.setValue(events[0].event_number); self.goto_spin.setEnabled(True); self.goto_spin.blockSignals(False)
        self.btn_prev.setEnabled(True); self.btn_next.setEnabled(True)

        self.file_label.setText(
            f"{os.path.basename(self._evio_path)} — {len(events):,} physics events "
            f"(scanned in {elapsed:.1f} s)")
        self._show_event(0)

    def _tear_down_worker(self):
        if self._scan_thread is not None:
            self._scan_thread.quit()
            self._scan_thread.wait(2000)
            self._scan_thread = None
        self._scan_worker = None

    def _close_current_run(self):
        if self._stepper is not None:
            self._stepper.close()
            self._stepper = None
        self._events = []
        self._current = -1
        self._last_ssp = None
        self.slider.setEnabled(False); self.goto_spin.setEnabled(False)
        self.btn_prev.setEnabled(False); self.btn_next.setEnabled(False)

    # -----------------------------------------------------------------
    # Event navigation
    # -----------------------------------------------------------------

    def _step(self, delta: int):
        if not self._events:
            return
        tgt = max(0, min(len(self._events) - 1, self._current + delta))
        if tgt != self._current:
            self._show_event(tgt)

    def _on_goto(self):
        """User typed an event number in the Goto spinbox.  Find the closest
        event (numbers may not be contiguous) and jump to it."""
        if not self._events:
            return
        want = int(self.goto_spin.value())
        # binary-ish search — events are sorted by event_number
        best = 0; best_d = abs(self._events[0].event_number - want)
        for i, e in enumerate(self._events):
            d = abs(e.event_number - want)
            if d < best_d:
                best_d = d; best = i
            if e.event_number >= want:
                break
        if best != self._current:
            self._show_event(best)

    def _on_slider(self, value: int):
        if not self._events:
            return
        if value != self._current:
            self._show_event(value)

    def _show_event(self, event_idx: int):
        if self._stepper is None or not self._events:
            return
        if not (0 <= event_idx < len(self._events)):
            return
        evmeta = self._events[event_idx]
        self._set_status(f"Fetching event #{evmeta.event_number}…")
        QApplication.processEvents()
        try:
            ssp = self._stepper.get_ssp(evmeta, event_idx)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Read failed", str(exc))
            return
        self._current = event_idx
        self._last_ssp = ssp

        # Reflect in nav widgets without looping back.
        self.slider.blockSignals(True); self.slider.setValue(event_idx); self.slider.blockSignals(False)
        self.goto_spin.blockSignals(True); self.goto_spin.setValue(evmeta.event_number); self.goto_spin.blockSignals(False)

        self._check_pedestal_requirement(ssp)
        self._re_reconstruct_current()

    def _check_pedestal_requirement(self, ssp) -> None:
        """If this event is full-readout (firmware did not run online ZS)
        and we have no pedestal file, warn the user — loudly once, and via
        a persistent status-bar indicator afterwards.

        Full-readout data without pedestals cannot produce meaningful hits:
        the default noise is 5000 ADC, nothing clears the ZS threshold, and
        every event looks empty.  Better to say so explicitly than leave
        the user staring at a silent canvas.
        """
        have_ped = bool(self._gem_ped_path and os.path.isfile(self._gem_ped_path))
        if have_ped:
            self._ped_badge.hide()
            return
        # Scan APVs for any in full-readout mode.  The authoritative signal
        # is nstrips == 128 (firmware sent every channel); has_online_cm
        # alone is unreliable — the MPD can emit CM debug headers while
        # still sending all 128 strips raw.
        APV_STRIP_SIZE = 128
        full_readout = False
        for m in range(ssp.nmpds):
            mpd = ssp.mpd(m)
            if not mpd.present:
                continue
            for a in range(16):  # ssp::MAX_APVS_PER_MPD
                apv = mpd.apv(a)
                if apv.present and apv.nstrips == APV_STRIP_SIZE:
                    full_readout = True
                    break
            if full_readout:
                break
        if not full_readout:
            # Data is online-ZS; peds not needed.
            self._ped_badge.hide()
            return

        # Full-readout + no pedestal = every event will look empty.
        self._ped_badge.setText("⚠ NO PEDESTAL FILE — full-readout data will reconstruct empty")
        self._ped_badge.show()
        if not self._ped_warning_shown:
            self._ped_warning_shown = True
            QMessageBox.warning(
                self, "Pedestal file required",
                "This file contains <b>full-readout</b> GEM data "
                "(no online zero-suppression), but no pedestal file is "
                "loaded.<br><br>"
                "Without pedestals, zero-suppression uses a default noise "
                "value and nothing will survive → every event will appear "
                "empty.<br><br>"
                "Either:<br>"
                "• Generate pedestals:  "
                "<code>gem_dump -m ped &lt;file.evio&gt; -o gem_ped.json</code><br>"
                "• Then load them via File → Choose gem_ped.json…")

    # -----------------------------------------------------------------
    # Reconstruction + drawing
    # -----------------------------------------------------------------

    def _on_threshold_change(self):
        # Debounce: coalesce multiple rapid changes (e.g. slider drag) into
        # one reconstruction pass.
        self._redraw_timer.start(self.REDRAW_DEBOUNCE_MS)

    def _reset_defaults(self):
        if self._gsys is None or self._gcl is None:
            return
        cfg = self._default_cluster_cfg
        if cfg is None:
            cfg = det.GemCluster().get_config()
        self._gcl.set_config(cfg)

        self.zs_spin.blockSignals(True); self.zs_spin.setValue(self._default_zs); self.zs_spin.blockSignals(False)
        self.cm_spin.blockSignals(True); self.cm_spin.setValue(self._default_cm); self.cm_spin.blockSignals(False)
        self.mch_spin.blockSignals(True); self.mch_spin.setValue(int(cfg.min_cluster_hits)); self.mch_spin.blockSignals(False)

        self.adv_dock.max_cluster_hits.blockSignals(True); self.adv_dock.max_cluster_hits.setValue(int(cfg.max_cluster_hits)); self.adv_dock.max_cluster_hits.blockSignals(False)
        self.adv_dock.consecutive_thres.blockSignals(True); self.adv_dock.consecutive_thres.setValue(int(cfg.consecutive_thres)); self.adv_dock.consecutive_thres.blockSignals(False)
        self.adv_dock.split_thres.blockSignals(True); self.adv_dock.split_thres.setValue(float(cfg.split_thres)); self.adv_dock.split_thres.blockSignals(False)
        self.adv_dock.cross_talk_width.blockSignals(True); self.adv_dock.cross_talk_width.setValue(int(cfg.cross_talk_width)); self.adv_dock.cross_talk_width.blockSignals(False)
        self.adv_dock.match_mode.blockSignals(True); self.adv_dock.match_mode.setCurrentIndex(int(cfg.match_mode)); self.adv_dock.match_mode.blockSignals(False)
        self.adv_dock.match_adc_asym.blockSignals(True); self.adv_dock.match_adc_asym.setValue(float(cfg.match_adc_asymmetry)); self.adv_dock.match_adc_asym.blockSignals(False)
        self.adv_dock.match_time_diff.blockSignals(True); self.adv_dock.match_time_diff.setValue(float(cfg.match_time_diff)); self.adv_dock.match_time_diff.blockSignals(False)
        self.adv_dock.ts_period.blockSignals(True); self.adv_dock.ts_period.setValue(float(cfg.ts_period)); self.adv_dock.ts_period.blockSignals(False)

        self._re_reconstruct_current()

    def _re_reconstruct_current(self):
        if self._gsys is None or self._gcl is None:
            return
        if self._last_ssp is None:
            return

        # Push widget values into the GemSystem / GemCluster.
        self._gsys.zero_sup_threshold = float(self.zs_spin.value())
        self._gsys.common_mode_threshold = float(self.cm_spin.value())
        cfg = self._gcl.get_config()
        cfg.min_cluster_hits = int(self.mch_spin.value())
        self._gcl.set_config(cfg)
        self.adv_dock.apply_to(self._gcl)

        # Run reconstruction on cached SSP.
        t0 = time.monotonic()
        self._gsys.clear()
        self._gsys.process_event(self._last_ssp)
        self._gsys.reconstruct(self._gcl)
        elapsed_ms = (time.monotonic() - t0) * 1000.0

        # Build the structures gem_view expects.
        det_list = build_det_list_from_gemsys(self._gsys)
        zs_apvs = build_zs_apvs_from_gemsys(self._gsys)
        det_hits = process_zs_hits(zs_apvs, self._apv_map,
                                   self._detectors, self._hole, self._gem_raw)

        evmeta = self._events[self._current]
        title = (f"GEM Event Viewer — "
                 f"ev #{evmeta.event_number}  "
                 f"trig #{evmeta.trigger_number}  "
                 f"bits 0x{evmeta.trigger_bits:X}")
        draw_event(self.canvas.fig, self._detectors, det_list, det_hits,
                   self._hole, title=title)
        self.canvas.draw_idle()

        n_2d = sum(len(d.get("hits_2d", [])) for d in det_list)
        self._set_status(
            f"#{self._current + 1}/{len(self._events)}  "
            f"ev={evmeta.event_number}  trig={evmeta.trigger_number}  "
            f"bits=0x{evmeta.trigger_bits:X}  "
            f"2D hits: {n_2d}   reco: {elapsed_ms:.1f} ms")

    # -----------------------------------------------------------------
    # Misc
    # -----------------------------------------------------------------

    def _toggle_advanced(self, checked: bool):
        self.adv_dock.setVisible(checked)
        self.act_adv.setChecked(checked)
        self.btn_advanced.setChecked(checked)

    def _set_status(self, text: str):
        self._status.setText(text)

    def _fatal_prad2py_missing(self):
        QMessageBox.critical(
            self, "prad2py not available",
            "The prad2py pybind11 module could not be imported.\n\n"
            "Build it with:\n"
            "    cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build\n\n"
            f"Details:\n{PRAD2PY_ERROR}")

    def closeEvent(self, ev):
        # Cancel any in-flight pre-scan before we tear down the run.
        if self._scan_worker is not None:
            self._scan_worker.request_cancel()
            self._tear_down_worker()
        self._close_current_run()
        super().closeEvent(ev)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Interactive PyQt6 GEM event viewer — step through an "
                    "EVIO file with live threshold tuning.")
    parser.add_argument("evio", nargs="?", help="EVIO file to open on start.")
    parser.add_argument("-D", "--daq-config", default=None,
                        help="Override daq_config.json path.")
    parser.add_argument("-G", "--gem-map", default=None,
                        help="Override gem_map.json path.")
    parser.add_argument("-P", "--gem-ped", default=None,
                        help="Pedestal file (required for full-readout data; "
                             "ignored for online-ZS production data).")
    args = parser.parse_args()

    app = QApplication(sys.argv)
    win = GemEventViewer(
        initial_evio=args.evio,
        daq_config_path=args.daq_config,
        gem_map_path=args.gem_map,
        gem_ped_path=args.gem_ped,
    )
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
