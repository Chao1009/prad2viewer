#!/usr/bin/env python3
"""
gem_hycal_match_viewer.py — interactive PyQt6 GEM↔HyCal matching viewer.

Open an EVIO file event-by-event and inspect the GEM↔HyCal coincidence:
  * front view (X-Y) — HyCal geometry with cluster centroids, GEM hits
    projected through the target onto the HyCal plane (one mark per
    detector, color-coded), and matching circles drawn at N·sigma_total;
  * side view (Z-Y) — target / GEM planes / HyCal face with hit markers
    and the HyCal→GEM line for each matched HC cluster;
  * match table — one row per (HC cluster × GEM detector) with residual
    and sigma_total;
  * toolbar — First/Prev/Next/Goto/Last for navigation, plus a
    "Next matched" search controlled by two thresholds (N hits per
    detector, K detectors with ≥N hits).  The search is a foreground scan
    with a cancellable progress dialog.
  * show/hide — checkboxes per detector + HyCal cluster overlay.

Usage:
    python scripts/gem_hycal_match_viewer.py <file.evio.00000> [-r RUN]

Configuration is read from `database/`:
  * monitor_config.json        (waveform binning, trigger filter)
  * daq_config.json            (DAQ + raw decoding)
  * reconstruction_config.json (runinfo pointer + matching constants)
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# prad2py + analysis._common discovery — walk up from this script to find
# build/python/ and add analysis/pyscripts/ to sys.path.
# ---------------------------------------------------------------------------

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR = _SCRIPT_DIR.parent
_probe = _SCRIPT_DIR
for _ in range(5):
    _probe = _probe.parent
    for _sub in ("build/python", "build-release/python", "build/Release/python"):
        _cand = _probe / _sub
        if _cand.is_dir() and str(_cand) not in sys.path:
            sys.path.insert(0, str(_cand))

# analysis/pyscripts/_common.py — parametric matching helpers shared with
# the offline TSV/CSV writer (gem_hycal_matching.py).
_ANA_PY = _REPO_DIR / "analysis" / "pyscripts"
if _ANA_PY.is_dir() and str(_ANA_PY) not in sys.path:
    sys.path.insert(0, str(_ANA_PY))

try:
    from prad2py import dec, det
    HAVE_PRAD2PY = True
    PRAD2PY_ERROR = ""
except Exception as _exc:
    dec = None  # type: ignore
    det = None  # type: ignore
    HAVE_PRAD2PY = False
    PRAD2PY_ERROR = f"{type(_exc).__name__}: {_exc}"

try:
    import _common as C  # type: ignore
    HAVE_COMMON = True
except Exception as _exc:
    C = None  # type: ignore
    HAVE_COMMON = False
    PRAD2PY_ERROR = (PRAD2PY_ERROR + "\n" if PRAD2PY_ERROR else "") + \
                    f"_common import: {type(_exc).__name__}: {_exc}"

from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QAction, QBrush, QColor, QFont, QKeySequence, QPainter, QPen
from PyQt6.QtWidgets import (
    QApplication, QCheckBox, QDoubleSpinBox, QFileDialog,
    QGroupBox, QHBoxLayout, QLabel, QMainWindow, QMessageBox, QProgressDialog,
    QPushButton, QSizePolicy, QSpinBox, QSplitter, QStatusBar, QTableWidget,
    QTableWidgetItem, QToolBar, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

# hycal_geoview is a sibling — provides the HyCal front-view widget and the
# theme system reused across event viewers.
from hycal_geoview import (
    HyCalMapWidget, apply_theme_palette, available_themes,
    load_modules as load_geo_modules, set_theme,
)

# Per-detector colour palette — same as the web viewer (resources/gem.js).
GEM_COLORS = [
    QColor("#ff7f0e"),  # GEM0 — orange
    QColor("#1f77b4"),  # GEM1 — blue
    QColor("#2ca02c"),  # GEM2 — green
    QColor("#d62728"),  # GEM3 — red
]
GEM_NAMES = ["GEM0", "GEM1", "GEM2", "GEM3"]


# ============================================================================
# Data structures
# ============================================================================

class HCCluster:
    __slots__ = ("idx", "x", "y", "z", "energy", "center_id", "nblocks",
                 "lab_x", "lab_y", "lab_z")

    def __init__(self, idx, h, lab):
        self.idx = idx
        self.x = float(h.x)
        self.y = float(h.y)
        self.z = 0.0
        self.energy = float(h.energy)
        self.center_id = int(h.center_id)
        self.nblocks = int(h.nblocks)
        self.lab_x, self.lab_y, self.lab_z = lab


class GEMHit:
    __slots__ = ("det_id", "lx", "ly", "lab_x", "lab_y", "lab_z",
                 "x_charge", "y_charge", "x_size", "y_size")

    def __init__(self, det_id, lx, ly, lab, x_charge, y_charge, x_size, y_size):
        self.det_id = det_id
        self.lx = float(lx)
        self.ly = float(ly)
        self.lab_x, self.lab_y, self.lab_z = lab
        self.x_charge = float(x_charge)
        self.y_charge = float(y_charge)
        self.x_size = int(x_size)
        self.y_size = int(y_size)


class Match:
    __slots__ = ("hc_idx", "det_id", "gem_idx", "proj_x", "proj_y",
                 "residual", "sigma_total")

    def __init__(self, hc_idx, det_id, gem_idx, proj_x, proj_y, residual, sigma_total):
        self.hc_idx = hc_idx
        self.det_id = det_id
        self.gem_idx = gem_idx
        self.proj_x = float(proj_x)
        self.proj_y = float(proj_y)
        self.residual = float(residual)
        self.sigma_total = float(sigma_total)


class EventResult:
    """Everything we need to render one event."""
    def __init__(self):
        self.event_num = 0
        self.trigger_bits = 0
        self.hc: List[HCCluster] = []
        self.gem: List[List[GEMHit]] = [[], [], [], []]
        self.matches: List[Match] = []   # one entry per (hc, det) best match


# ============================================================================
# Matching loop (parametric sigma — mirrors gem_hycal_matching.py / .C)
# ============================================================================

def compute_matches(hc: List[HCCluster], gem: List[List[GEMHit]],
                    pr_A: float, pr_B: float, pr_C: float,
                    gem_pos_res: List[float], match_nsigma: float
                    ) -> List[Match]:
    """For each (HC cluster × GEM detector) pair, find the closest GEM hit
    inside `match_nsigma · σ_total` at the GEM plane.  At most one match per
    (HC, det)."""
    out: List[Match] = []
    for k, h in enumerate(hc):
        if h.lab_z <= 0:
            continue
        sig_face = C.hycal_pos_resolution(pr_A, pr_B, pr_C, h.energy)
        for d in range(4):
            gl = gem[d]
            if not gl:
                continue
            z_gem = gl[0].lab_z
            if z_gem <= 0:
                continue
            scale = z_gem / h.lab_z
            proj_x = h.lab_x * scale
            proj_y = h.lab_y * scale
            sig_hc_at_gem = sig_face * scale
            sig_gem = gem_pos_res[d] if d < len(gem_pos_res) else 0.1
            sig_total = math.sqrt(sig_hc_at_gem**2 + sig_gem**2)
            cut = match_nsigma * sig_total
            best_gi = -1
            best_dr = cut
            for gi, g in enumerate(gl):
                dx = g.lab_x - proj_x
                dy = g.lab_y - proj_y
                dr = math.sqrt(dx*dx + dy*dy)
                if dr <= best_dr:
                    best_dr = dr
                    best_gi = gi
            if best_gi >= 0:
                out.append(Match(k, d, best_gi, proj_x, proj_y, best_dr, sig_total))
    return out


def event_passes(matches: List[Match], min_hits_per_det: int, min_dets: int) -> bool:
    """Event qualifies if at least `min_dets` GEM detectors each have at
    least `min_hits_per_det` matched hits (counted across all HC clusters)."""
    counts = [0, 0, 0, 0]
    for m in matches:
        counts[m.det_id] += 1
    n_pass = sum(1 for c in counts if c >= min_hits_per_det)
    return n_pass >= min_dets


# ============================================================================
# Reconstruction pipeline (uses prad2py, matches gem_hycal_matching.py setup)
# ============================================================================

class Pipeline:
    """Wraps `_common.setup_pipeline` so the viewer reconstructs identically
    to the offline TSV/CSV writer.  Owns the matching constants too."""

    def __init__(self, db_dir: Path, run_num: int, evio_path: Path):
        self.db_dir = db_dir
        self.run_num = run_num
        self.evio_path = evio_path
        self.match_nsigma = 3.0
        # Initialize via the shared helper.  daq_config="" → installed default;
        # we override the env var so the helper looks in the user's db_dir.
        os.environ.setdefault("PRAD2_DATABASE_DIR", str(db_dir))
        self._p = C.setup_pipeline(
            evio_path=str(evio_path),
            run_num=run_num,
        )
        self.daq_cfg      = self._p.cfg
        self.hycal        = self._p.hycal
        self.hc_clusterer = self._p.hc_clusterer
        self.gem_sys      = self._p.gem_sys
        self.gem_clusterer = self._p.gem_clusterer
        self.wave_ana     = self._p.wave_ana
        self.geo          = self._p.geo

        # Matching config (parametric sigma).  Push A,B,C into HyCalSystem
        # so PositionResolution(E) works the same as the monitor does.
        (A, B, Cc), gpr = C.load_matching_config()
        self.match_A, self.match_B, self.match_C = A, B, Cc
        self.gem_pos_res = (list(gpr) + [0.1] * 4)[:4] if gpr else [0.1] * 4
        try:
            self.hycal.set_position_resolution_params(A, B, Cc)
        except AttributeError:
            pass  # binding not yet built — we still use values from Python

    # -----------------------------------------------------------------------
    # Per-event reco — returns an EventResult.
    # -----------------------------------------------------------------------
    def reconstruct(self, fadc_evt, ssp_evt) -> EventResult:
        ev = EventResult()
        ev.event_num = int(fadc_evt.info.event_number)
        ev.trigger_bits = int(fadc_evt.info.trigger_bits)

        # HyCal: waveform → energy → cluster (same logic as gem_hycal_matching.py).
        hc_raw = C.reconstruct_hycal(self._p, fadc_evt)
        for k, h in enumerate(hc_raw):
            z_local = det.shower_depth(h.center_id, h.energy)
            lab = (C.transform_hycal(h.x, h.y, z_local, self.geo)
                   if self.geo else (float(h.x), float(h.y), 0.0))
            ev.hc.append(HCCluster(k, h, lab))

        # GEM: pedestal + CM + ZS → 1D + 2D
        C.reconstruct_gem(self._p, ssp_evt)
        for d in range(min(4, self.gem_sys.get_n_detectors())):
            for g in self.gem_sys.get_hits(d):
                lab = (C.transform_gem(g.x, g.y, 0.0, d, self.geo)
                       if self.geo else (float(g.x), float(g.y), 0.0))
                ev.gem[d].append(GEMHit(d, g.x, g.y, lab,
                                        g.x_charge, g.y_charge,
                                        g.x_size, g.y_size))

        # Matching
        ev.matches = compute_matches(ev.hc, ev.gem,
                                     self.match_A, self.match_B, self.match_C,
                                     self.gem_pos_res, self.match_nsigma)
        return ev


# ============================================================================
# HyCal front-view subclass — overlay GEM-projected hits + matching markers.
# ============================================================================

class FrontView(HyCalMapWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._show_hc = True
        self._show_gem = [True, True, True, True]
        self._show_matches = True
        self._evt: Optional[EventResult] = None
        self._z_hc = 6225.0   # default; updated when geometry loads
        self._z_gem = [5400.0] * 4

    def set_event(self, evt: Optional[EventResult]):
        self._evt = evt
        self.update()

    def set_zs(self, z_hc: float, z_gem: List[float]):
        if z_hc > 0:
            self._z_hc = z_hc
        for i, z in enumerate(z_gem[:4]):
            if z > 0:
                self._z_gem[i] = z

    def set_show_hc(self, on: bool):
        self._show_hc = on; self.update()

    def set_show_gem(self, det_id: int, on: bool):
        if 0 <= det_id < 4:
            self._show_gem[det_id] = on; self.update()

    def set_show_matches(self, on: bool):
        self._show_matches = on; self.update()

    def _paint_overlays(self, p: QPainter, w: int, h: int):
        super()._paint_overlays(p, w, h)
        if not self._evt:
            return
        ev = self._evt

        # HyCal cluster centroids (HyCal-local x,y).
        if self._show_hc:
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            for c in ev.hc:
                pt = self.geo_to_canvas(c.x, c.y)
                p.setPen(QPen(QColor("#ffffff"), 2.0))
                p.setBrush(QBrush(QColor(255, 255, 255, 90)))
                p.drawEllipse(pt, 6.0, 6.0)
                p.setPen(QPen(QColor("#000000"), 1.0))
                font = QFont(p.font()); font.setPointSize(8); p.setFont(font)
                p.drawText(pt + QPointF(8, -2), f"HC{c.idx}: {c.energy:.0f} MeV")

        # GEM hits projected through the target onto HyCal-local x,y.
        # px = lab_x · (z_hc / z_gem),  py = lab_y · (z_hc / z_gem).
        # HyCal-local equals lab x,y for an untilted HyCal centred at (0,0,z_hc),
        # which is the standard PRad-II geometry; if hycal_x/y or tilts are non-
        # zero the runinfo loader has already absorbed them into lab coords.
        for d in range(4):
            if not self._show_gem[d]:
                continue
            color = GEM_COLORS[d]
            for g in ev.gem[d]:
                z_gem = g.lab_z if g.lab_z > 0 else self._z_gem[d]
                if z_gem <= 0:
                    continue
                scale = self._z_hc / z_gem
                px = g.lab_x * scale
                py = g.lab_y * scale
                qp = self.geo_to_canvas(px, py)
                p.setPen(QPen(color, 1.4))
                p.setBrush(QBrush(QColor(color.red(), color.green(), color.blue(), 130)))
                p.drawEllipse(qp, 4.0, 4.0)

        # Matching circle around HC cluster — drawn at HC plane.
        # Radius = match_nsigma · σ_total at HC plane (project σ_GEM up via
        # z_hc/z_gem).  We don't have nsigma here, but the match record
        # carries σ_total directly at the GEM plane; at HC it scales by
        # z_hc/z_gem.
        if self._show_matches and ev.matches:
            for m in ev.matches:
                if not self._show_gem[m.det_id]:
                    continue
                if m.hc_idx >= len(ev.hc):
                    continue
                c = ev.hc[m.hc_idx]
                # Best-match GEM hit projected to HC plane → endpoint of line.
                gl = ev.gem[m.det_id]
                if m.gem_idx >= len(gl):
                    continue
                g = gl[m.gem_idx]
                z_gem = g.lab_z if g.lab_z > 0 else self._z_gem[m.det_id]
                scale = self._z_hc / z_gem if z_gem > 0 else 1.0
                px = g.lab_x * scale
                py = g.lab_y * scale
                a = self.geo_to_canvas(c.x, c.y)
                b = self.geo_to_canvas(px, py)
                pen = QPen(GEM_COLORS[m.det_id], 1.6)
                pen.setStyle(Qt.PenStyle.DashLine)
                p.setPen(pen)
                p.drawLine(a, b)


# ============================================================================
# Side view (Z-Y) — matplotlib canvas
# ============================================================================

class SideView(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self._fig = Figure(figsize=(6, 4), tight_layout=True)
        super().__init__(self._fig)
        if parent is not None:
            self.setParent(parent)
        self._ax = self._fig.add_subplot(111)
        self._ax.set_xlabel("z (mm)")
        self._ax.set_ylabel("y (mm)")
        self._evt: Optional[EventResult] = None
        self._z_hc = 6225.0
        self._z_gem = [5400.0] * 4
        self._x_size_gem = [600.0] * 4
        self._y_size_gem = [600.0] * 4
        self._show_hc = True
        self._show_gem = [True, True, True, True]
        self._show_matches = True

    def set_geom(self, z_hc: float, z_gem: List[float],
                 y_size_gem: List[float]):
        if z_hc > 0:
            self._z_hc = z_hc
        for i, z in enumerate(z_gem[:4]):
            if z > 0:
                self._z_gem[i] = z
        for i, y in enumerate(y_size_gem[:4]):
            if y > 0:
                self._y_size_gem[i] = y

    def set_event(self, evt: Optional[EventResult]):
        self._evt = evt
        self.redraw()

    def set_show_hc(self, on: bool):    self._show_hc = on;       self.redraw()
    def set_show_matches(self, on: bool):
        self._show_matches = on; self.redraw()

    def set_show_gem(self, det_id: int, on: bool):
        if 0 <= det_id < 4:
            self._show_gem[det_id] = on; self.redraw()

    def redraw(self):
        ax = self._ax
        ax.clear()
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(True, color="#888", alpha=0.2, linewidth=0.5)

        # Detector frames (dashed) — GEM and HyCal active areas in y.
        for d in range(4):
            yh = self._y_size_gem[d] / 2
            color = (GEM_COLORS[d].redF(), GEM_COLORS[d].greenF(), GEM_COLORS[d].blueF())
            ax.plot([self._z_gem[d], self._z_gem[d]], [-yh, yh],
                    "--", color=color, alpha=0.6, linewidth=1.0)
            ax.text(self._z_gem[d], yh + 20, GEM_NAMES[d],
                    color=color, fontsize=8, ha="center")
        ax.axvline(self._z_hc, color="#cccccc", linestyle="--", linewidth=1.0)
        ax.text(self._z_hc, 0, "HyCal", color="#cccccc", fontsize=8,
                ha="center", va="bottom", rotation=90)
        ax.axvline(0.0, color="#888", linestyle=":", linewidth=0.8)
        ax.text(0, 0, "T", color="#888", fontsize=8, ha="right")

        evt = self._evt
        if evt:
            # HyCal cluster markers (z_hc, lab_y)
            if self._show_hc:
                for c in evt.hc:
                    ax.plot([c.lab_z], [c.lab_y], "s", color="white",
                            markersize=8, markeredgecolor="black",
                            markeredgewidth=0.8)
            # GEM hits per detector
            for d in range(4):
                if not self._show_gem[d]:
                    continue
                color = (GEM_COLORS[d].redF(), GEM_COLORS[d].greenF(), GEM_COLORS[d].blueF())
                xs = [g.lab_z for g in evt.gem[d]]
                ys = [g.lab_y for g in evt.gem[d]]
                if xs:
                    ax.plot(xs, ys, "o", color=color, markersize=5,
                            markeredgecolor="black", markeredgewidth=0.4,
                            alpha=0.85)
            # Matched HC↔GEM lines
            if self._show_matches:
                for m in evt.matches:
                    if not self._show_gem[m.det_id]:
                        continue
                    if m.hc_idx >= len(evt.hc):
                        continue
                    c = evt.hc[m.hc_idx]
                    gl = evt.gem[m.det_id]
                    if m.gem_idx >= len(gl):
                        continue
                    g = gl[m.gem_idx]
                    color = (GEM_COLORS[m.det_id].redF(),
                             GEM_COLORS[m.det_id].greenF(),
                             GEM_COLORS[m.det_id].blueF())
                    ax.plot([g.lab_z, c.lab_z], [g.lab_y, c.lab_y],
                            "-", color=color, linewidth=1.0, alpha=0.7)

        # Y range — pad around the largest detector size.
        y_max = max(max(self._y_size_gem) / 2, 50.0)
        ax.set_ylim(-y_max * 1.1, y_max * 1.1)
        ax.set_xlim(-200, self._z_hc * 1.05)
        self.draw_idle()


# ============================================================================
# Main window
# ============================================================================

class GemHycalMatchViewer(QMainWindow):
    def __init__(self, evio_path: Optional[Path] = None,
                 db_dir: Optional[Path] = None,
                 run_num: int = -1):
        super().__init__()
        self.setWindowTitle("GEM↔HyCal Matching Viewer")
        self.resize(1500, 900)

        if not HAVE_PRAD2PY or not HAVE_COMMON:
            QMessageBox.critical(self, "prad2py / _common not available",
                                 PRAD2PY_ERROR or "Cannot find prad2py or "
                                 "analysis/pyscripts/_common.py — build the "
                                 "Python bindings and re-run.")
            sys.exit(1)

        self._db_dir = Path(db_dir or os.environ.get(
            "PRAD2_DATABASE_DIR", _REPO_DIR / "database")).resolve()
        self._run_num = run_num
        self._pipeline: Optional[Pipeline] = None
        self._physics_index: List[Tuple[int, int]] = []  # (record_idx, sub_idx)
        self._cur_idx = -1
        self._cur_event: Optional[EventResult] = None

        self._build_ui()

        if evio_path is not None:
            self._open_file(Path(evio_path))

    # ---- UI ---------------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        outer = QVBoxLayout(central); outer.setContentsMargins(4, 4, 4, 4)

        self._build_toolbar()
        self._build_search_row(outer)
        self._build_visibility_row(outer)
        self._build_views(outer)
        self._build_match_table(outer)

        self._status = QStatusBar()
        self.setStatusBar(self._status)
        self._update_status_bar()

        # File menu
        m = self.menuBar().addMenu("&File")
        act_open = QAction("&Open EVIO…", self)
        act_open.setShortcut(QKeySequence.StandardKey.Open)
        act_open.triggered.connect(self._on_open)
        m.addAction(act_open)
        act_quit = QAction("&Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        m.addAction(act_quit)

        apply_theme_palette(self)

    def _build_toolbar(self):
        tb = QToolBar()
        tb.setMovable(False)
        self.addToolBar(tb)

        def addBtn(text, slot, shortcut=None, tip=None):
            act = QAction(text, self)
            if shortcut:
                act.setShortcut(shortcut)
            if tip:
                act.setToolTip(tip)
            act.triggered.connect(slot)
            tb.addAction(act)
            return act

        addBtn("⏮ First", self._first, "Home", "Jump to first physics event")
        addBtn("◀ Prev",  self._prev,  "Left", "Previous physics event")
        tb.addSeparator()
        tb.addWidget(QLabel(" Event #"))
        self._sb_idx = QSpinBox()
        self._sb_idx.setRange(0, 0)
        self._sb_idx.setKeyboardTracking(False)
        self._sb_idx.editingFinished.connect(
            lambda: self._goto(self._sb_idx.value()))
        tb.addWidget(self._sb_idx)
        addBtn("Goto", lambda: self._goto(self._sb_idx.value()))
        tb.addSeparator()
        addBtn("Next ▶", self._next, "Right", "Next physics event")
        addBtn("Last ⏭", self._last, "End", "Jump to last physics event")

    def _build_search_row(self, outer: QVBoxLayout):
        box = QGroupBox("Next-matched search")
        lay = QHBoxLayout(box); lay.setContentsMargins(8, 4, 8, 4)
        lay.addWidget(QLabel("≥"))
        self._sb_N = QSpinBox(); self._sb_N.setRange(1, 99); self._sb_N.setValue(1)
        lay.addWidget(self._sb_N)
        lay.addWidget(QLabel("matched hits per detector,"))
        lay.addWidget(QLabel("≥"))
        self._sb_K = QSpinBox(); self._sb_K.setRange(1, 4); self._sb_K.setValue(2)
        lay.addWidget(self._sb_K)
        lay.addWidget(QLabel("detectors satisfied"))
        btn = QPushButton("Find next ▶▶"); btn.setShortcut("Shift+Right")
        btn.clicked.connect(self._find_next_matched)
        lay.addWidget(btn)
        lay.addStretch(1)

        # nsigma override
        lay.addWidget(QLabel("  nσ:"))
        self._sb_ns = QDoubleSpinBox()
        self._sb_ns.setRange(0.5, 20.0); self._sb_ns.setSingleStep(0.5)
        self._sb_ns.setValue(3.0)
        self._sb_ns.valueChanged.connect(self._on_nsigma_changed)
        lay.addWidget(self._sb_ns)

        outer.addWidget(box)

    def _build_visibility_row(self, outer: QVBoxLayout):
        row = QHBoxLayout()
        row.setContentsMargins(8, 0, 8, 0)
        self._cb_hc = QCheckBox("HyCal"); self._cb_hc.setChecked(True)
        self._cb_hc.toggled.connect(self._on_show_hc)
        row.addWidget(self._cb_hc)
        self._cb_gem: List[QCheckBox] = []
        for d in range(4):
            cb = QCheckBox(GEM_NAMES[d]); cb.setChecked(True)
            cb.setStyleSheet(
                f"color: {GEM_COLORS[d].name()}; font-weight: bold;")
            cb.toggled.connect(lambda on, dd=d: self._on_show_gem(dd, on))
            self._cb_gem.append(cb)
            row.addWidget(cb)
        self._cb_match = QCheckBox("Matches"); self._cb_match.setChecked(True)
        self._cb_match.toggled.connect(self._on_show_matches)
        row.addWidget(self._cb_match)
        row.addStretch(1)
        outer.addLayout(row)

    def _build_views(self, outer: QVBoxLayout):
        split = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(split, stretch=1)

        self._front = FrontView()
        self._front.setSizePolicy(QSizePolicy.Policy.Expanding,
                                  QSizePolicy.Policy.Expanding)
        split.addWidget(self._front)

        self._side = SideView()
        split.addWidget(self._side)
        split.setSizes([700, 700])

    def _build_match_table(self, outer: QVBoxLayout):
        self._tbl = QTableWidget(0, 7)
        self._tbl.setHorizontalHeaderLabels(
            ["HC#", "GEM", "GEM x (mm)", "GEM y (mm)",
             "residual (mm)", "σ_total (mm)", "ratio"])
        self._tbl.horizontalHeader().setStretchLastSection(True)
        self._tbl.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._tbl.setMaximumHeight(160)
        outer.addWidget(self._tbl)

    # ---- File handling ----------------------------------------------------
    def _on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open EVIO file", str(_REPO_DIR),
            "EVIO files (*.evio *.evio.*);;All files (*)")
        if path:
            self._open_file(Path(path))

    def _open_file(self, path: Path):
        if not path.is_file():
            QMessageBox.warning(self, "EVIO not found", str(path))
            return
        try:
            self._pipeline = Pipeline(self._db_dir, self._run_num, path)
        except (Exception, SystemExit) as exc:
            QMessageBox.critical(self, "Pipeline init failed",
                                 f"{type(exc).__name__}: {exc}")
            return

        # Index physics events with a progress dialog.  The channel is
        # released afterwards — every navigation call opens a fresh channel.
        ch = dec.EvChannel()
        ch.set_config(self._pipeline.daq_cfg)
        if ch.open_auto(str(path)) != dec.Status.success:
            QMessageBox.critical(self, "Cannot open EVIO", str(path))
            return
        self._physics_index = self._index_physics_events(ch)
        if not self._physics_index:
            QMessageBox.information(self, "No physics events",
                                    "Scanned the file but found no physics records.")
            return

        # Push detector geometry into the views.
        if self._pipeline.geo:
            geo = self._pipeline.geo
            self._front.set_zs(geo.hycal_z, geo.gem_z)
            y_size = []
            for d in range(min(4, self._pipeline.gem_sys.get_n_detectors())):
                dets = self._pipeline.gem_sys.get_detectors()
                y_size.append(float(dets[d].planes[1].size))
            while len(y_size) < 4:
                y_size.append(600.0)
            self._side.set_geom(geo.hycal_z, geo.gem_z, y_size)

        # Push HyCal modules into the front view.
        modules_rel = self._pipeline.daq_cfg.modules_file or "hycal_modules.json"
        modules_path = Path(modules_rel)
        if not modules_path.is_absolute():
            modules_path = self._db_dir / modules_rel
        if modules_path.is_file():
            try:
                modules = load_geo_modules(modules_path)
                self._front.set_modules(modules)
            except Exception as exc:
                self.statusBar().showMessage(f"hycal_modules.json: {exc}")

        self._sb_idx.setRange(0, len(self._physics_index) - 1)
        self.setWindowTitle(f"GEM↔HyCal Matching Viewer — {path.name}")
        self._goto(0)

    def _index_physics_events(self, ch) -> List[Tuple[int, int]]:
        idx: List[Tuple[int, int]] = []
        dlg = QProgressDialog("Indexing physics events…", "Cancel", 0, 0, self)
        dlg.setWindowModality(Qt.WindowModality.ApplicationModal)
        dlg.setMinimumDuration(250)
        dlg.show()
        rec = 0
        last_ui = time.monotonic()
        while ch.read() == dec.Status.success:
            rec += 1
            if not ch.scan():
                continue
            if ch.get_event_type() != dec.EventType.Physics:
                continue
            for sub in range(ch.get_n_events()):
                idx.append((rec, sub))
            now = time.monotonic()
            if now - last_ui > 0.1:
                dlg.setLabelText(f"Indexing physics events… {len(idx):,} found")
                QApplication.processEvents()
                last_ui = now
                if dlg.wasCanceled():
                    break
        dlg.close()
        # Rewind for sub-event re-decoding.  EvChannel doesn't expose seek;
        # we keep reading sequentially and decode lazily by record id.
        # Instead, _goto() opens a fresh channel for each random access.
        return idx

    # ---- Navigation -------------------------------------------------------
    def _first(self): self._goto(0)
    def _prev(self):  self._goto(self._cur_idx - 1)
    def _next(self):  self._goto(self._cur_idx + 1)
    def _last(self):  self._goto(len(self._physics_index) - 1)

    def _goto(self, idx: int):
        if not self._physics_index:
            return
        idx = max(0, min(idx, len(self._physics_index) - 1))
        if idx == self._cur_idx and self._cur_event is not None:
            return
        self._cur_idx = idx
        self._cur_event = self._decode_at(idx)
        if self._cur_event is not None:
            self._front.set_event(self._cur_event)
            self._side.set_event(self._cur_event)
            self._populate_match_table(self._cur_event)
        self._sb_idx.blockSignals(True)
        self._sb_idx.setValue(idx); self._sb_idx.blockSignals(False)
        self._update_status_bar()

    def _decode_at(self, idx: int) -> Optional[EventResult]:
        """Fresh channel + sequential read up to the target record/sub-event."""
        if idx < 0 or idx >= len(self._physics_index):
            return None
        target_rec, target_sub = self._physics_index[idx]
        ch = dec.EvChannel()
        ch.set_config(self._pipeline.daq_cfg)
        if ch.open_auto(str(self._pipeline.evio_path)) != dec.Status.success:
            return None
        rec = 0
        while ch.read() == dec.Status.success:
            rec += 1
            if rec != target_rec:
                continue
            if not ch.scan():
                return None
            decoded = ch.decode_event(target_sub, with_ssp=True)
            if not decoded.get("ok"):
                return None
            return self._pipeline.reconstruct(decoded["event"], decoded["ssp"])
        return None

    # ---- Search -----------------------------------------------------------
    def _find_next_matched(self):
        if not self._physics_index:
            return
        N = self._sb_N.value(); K = self._sb_K.value()
        start = self._cur_idx + 1
        end   = len(self._physics_index)
        dlg = QProgressDialog(
            f"Searching for an event with ≥{K} detectors at ≥{N} hit(s)…",
            "Cancel", start, end, self)
        dlg.setWindowModality(Qt.WindowModality.ApplicationModal)
        dlg.setMinimumDuration(250)
        # Open a fresh channel and walk records sequentially; this is much
        # faster than calling _decode_at() per event (which reopens the file).
        ch = dec.EvChannel()
        ch.set_config(self._pipeline.daq_cfg)
        if ch.open_auto(str(self._pipeline.evio_path)) != dec.Status.success:
            return
        rec = 0
        target_rec, target_sub = self._physics_index[start] if start < end else (None, None)
        cursor = start
        last_ui = time.monotonic()
        found = -1
        while ch.read() == dec.Status.success:
            rec += 1
            if cursor >= end:
                break
            if target_rec is None:
                break
            if rec != target_rec:
                continue
            if not ch.scan():
                # Advance past every (rec, sub) on this skipped record.
                while cursor < end and self._physics_index[cursor][0] == rec:
                    cursor += 1
                if cursor < end:
                    target_rec, target_sub = self._physics_index[cursor]
                continue
            # Walk every sub-event in this record while we still need them.
            while cursor < end and self._physics_index[cursor][0] == rec:
                _, sub = self._physics_index[cursor]
                decoded = ch.decode_event(sub, with_ssp=True)
                if decoded.get("ok"):
                    evr = self._pipeline.reconstruct(decoded["event"], decoded["ssp"])
                    if event_passes(evr.matches, N, K):
                        found = cursor
                        self._cur_event = evr
                        break
                cursor += 1
                now = time.monotonic()
                if now - last_ui > 0.1:
                    dlg.setValue(cursor); QApplication.processEvents()
                    last_ui = now
                    if dlg.wasCanceled():
                        cursor = end
                        break
            if found >= 0:
                break
            if cursor < end:
                target_rec, target_sub = self._physics_index[cursor]
        dlg.close()
        if found >= 0:
            self._cur_idx = found
            self._sb_idx.blockSignals(True)
            self._sb_idx.setValue(found); self._sb_idx.blockSignals(False)
            self._front.set_event(self._cur_event)
            self._side.set_event(self._cur_event)
            self._populate_match_table(self._cur_event)
            self._update_status_bar()
        else:
            self.statusBar().showMessage(
                "No matching event found before EOF.", 5000)

    # ---- Visibility -------------------------------------------------------
    def _on_show_hc(self, on: bool):
        self._front.set_show_hc(on); self._side.set_show_hc(on)

    def _on_show_gem(self, det_id: int, on: bool):
        self._front.set_show_gem(det_id, on); self._side.set_show_gem(det_id, on)

    def _on_show_matches(self, on: bool):
        self._front.set_show_matches(on); self._side.set_show_matches(on)

    def _on_nsigma_changed(self, v: float):
        if self._pipeline:
            self._pipeline.match_nsigma = float(v)
        # Re-run matching on the current event without re-decoding.
        if self._cur_event:
            self._cur_event.matches = compute_matches(
                self._cur_event.hc, self._cur_event.gem,
                self._pipeline.match_A, self._pipeline.match_B,
                self._pipeline.match_C, self._pipeline.gem_pos_res,
                self._pipeline.match_nsigma)
            self._front.set_event(self._cur_event)
            self._side.set_event(self._cur_event)
            self._populate_match_table(self._cur_event)

    # ---- Match table ------------------------------------------------------
    def _populate_match_table(self, evt: EventResult):
        self._tbl.setRowCount(len(evt.matches))
        for r, m in enumerate(evt.matches):
            ratio = m.residual / m.sigma_total if m.sigma_total > 0 else float("inf")
            cells = [str(m.hc_idx), GEM_NAMES[m.det_id],
                     f"{evt.gem[m.det_id][m.gem_idx].lab_x:.2f}",
                     f"{evt.gem[m.det_id][m.gem_idx].lab_y:.2f}",
                     f"{m.residual:.2f}", f"{m.sigma_total:.2f}",
                     f"{ratio:.2f}σ"]
            for c, txt in enumerate(cells):
                item = QTableWidgetItem(txt)
                if c == 1:
                    item.setForeground(GEM_COLORS[m.det_id])
                self._tbl.setItem(r, c, item)

    # ---- Status -----------------------------------------------------------
    def _update_status_bar(self):
        if not self._physics_index:
            self._status.showMessage("No file loaded.")
            return
        if self._cur_event is None:
            self._status.showMessage(
                f"event {self._cur_idx + 1} / {len(self._physics_index)}  (decoding…)")
            return
        ev = self._cur_event
        # Count matched detectors (≥1 hit).
        per_det = [0, 0, 0, 0]
        for m in ev.matches:
            per_det[m.det_id] += 1
        n_dets = sum(1 for c in per_det if c > 0)
        msg = (f"event {self._cur_idx + 1}/{len(self._physics_index)}  "
               f"(#{ev.event_num})  trig=0x{ev.trigger_bits:08X}  "
               f"HC={len(ev.hc)}  matches={len(ev.matches)} on {n_dets} det "
               f"[{per_det[0]},{per_det[1]},{per_det[2]},{per_det[3]}]")
        self._status.showMessage(msg)


# ============================================================================
# Entry point
# ============================================================================

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("evio", nargs="?", help="EVIO file to open on startup")
    ap.add_argument("--db", default=None,
                    help="database directory (default: $PRAD2_DATABASE_DIR or "
                         "<repo>/database)")
    ap.add_argument("-r", "--run", type=int, default=-1,
                    help="run number for runinfo lookup (default: sniff filename)")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="colour theme")
    args = ap.parse_args(argv)

    set_theme(args.theme)
    app = QApplication(sys.argv)
    w = GemHycalMatchViewer(
        evio_path=Path(args.evio) if args.evio else None,
        db_dir=Path(args.db) if args.db else None,
        run_num=args.run,
    )
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
