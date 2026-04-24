#!/usr/bin/env python3
"""
ep_calib_viewer.py  —  epCalib output visualization GUI
========================================================
Reads Physics_calib/{run}/calib_iter{N}.json  +  CalibResult_iter{N}.root

ROOT file layout assumed:
  module_energy/h_{name}   — per-module energy TH1F
  ratio_all                — ratio distribution TH1F
  h2_energy_theta          — 2D energy-theta TH2F
  h2_Nevents_moduleMap     — 2D (col, row) event map TH2F
  measured_peak            — distribution of measured peak values
  recon_sigma              — distribution of sigma values

Usage:
  python scripts/ep_calib_viewer.py [Physics_calib_dir]
  python scripts/ep_calib_viewer.py build/Physics_calib
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
from PyQt6.QtCore import Qt, QRectF, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QFont
from PyQt6.QtWidgets import (
    QApplication, QButtonGroup, QCheckBox, QComboBox, QFileDialog,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton,
    QSizePolicy, QSplitter, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

try:
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ── find hycal_geoview  ───────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from hycal_geoview import (          # noqa: E402  (after sys.path tweak)
    HyCalMapWidget, Module, load_modules, cmap_qcolor,
    PALETTE_NAMES, apply_theme_palette, set_theme, available_themes, THEME,
)

# ── database path  ────────────────────────────────────────────────────────────
_DB_CANDIDATES = [
    _SCRIPT_DIR.parent / "database" / "hycal_modules.json",
    _SCRIPT_DIR.parent / "build"    / "database" / "hycal_modules.json",
]
MODULES_JSON = next((p for p in _DB_CANDIDATES if p.is_file()), None)


# =============================================================================
#  Data structures
# =============================================================================

class ModuleMetrics:
    """Per-module computed metrics for one (run, iteration) pair."""
    __slots__ = ("name", "stats", "peak", "sigma", "chi2",
                 "factor", "base_energy")

    def __init__(self, name: str):
        self.name        = name
        self.stats       = 0.0    # event count (histogram integral)
        self.peak        = 0.0    # measured peak (MeV)
        self.sigma       = 0.0    # Gaussian sigma (MeV)
        self.chi2        = 0.0    # chi2/ndf of Gaussian fit
        self.factor      = 1.0    # calibration factor from JSON
        self.base_energy = 0.0    # expected peak from JSON (MeV)


class IterData:
    """All data for one (run, iteration) pair."""
    def __init__(self):
        self.metrics:        Dict[str, ModuleMetrics] = {}
        self.root_path:      Optional[Path]           = None
        self.factors:        Dict[str, dict]          = {}   # raw JSON entries
        self.expected_peaks: Dict[str, float]         = {}   # from .dat ExpectedPeak
        self._hist_cache:    Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    @property
    def json_path(self) -> Optional[Path]:
        """Derive the calib_iter{N}.json path from root_path."""
        if self.root_path is None:
            return None
        return self.root_path.with_name(
            self.root_path.name
                .replace("CalibResult_iter", "calib_iter")
                .replace(".root", ".json")
        )

    @property
    def dat_path(self) -> Optional[Path]:
        """Derive the fitting_parameters_iter{N}.dat path from root_path."""
        if self.root_path is None:
            return None
        return self.root_path.with_name(
            self.root_path.name
                .replace("CalibResult_iter", "fitting_parameters_iter")
                .replace(".root", ".dat")
        )


# =============================================================================
#  I/O helpers
# =============================================================================

def scan_calib_dir(base: Path) -> Dict[str, Dict[int, IterData]]:
    """Return {run_str: {iter_num: IterData}} discovered under *base*."""
    result: Dict[str, Dict[int, IterData]] = {}
    if not base.is_dir():
        return result
    for run_dir in sorted(base.iterdir()):
        if not run_dir.is_dir():
            continue
        iters: Dict[int, IterData] = {}
        for jf in sorted(run_dir.glob("calib_iter*.json")):
            m = re.search(r"calib_iter(\d+)\.json", jf.name)
            if not m:
                continue
            it   = int(m.group(1))
            data = IterData()
            with open(jf) as f:
                entries = json.load(f)
            data.factors = {e["name"]: e for e in entries}
            rf = run_dir / f"CalibResult_iter{it}.root"
            if rf.is_file():
                data.root_path = rf
            # parse fitting_parameters_iter{N}.dat for ExpectedPeak
            df = run_dir / f"fitting_parameters_iter{it}.dat"
            if df.is_file():
                data.expected_peaks = _parse_dat_expected_peaks(df)
            iters[it] = data
        if iters:
            result[run_dir.name] = iters
    return result


def _parse_dat_expected_peaks(dat_path: Path) -> Dict[str, float]:
    """Parse fitting_parameters_iter{N}.dat and return {module: ExpectedPeak}."""
    result: Dict[str, float] = {}
    try:
        for line in dat_path.read_text().splitlines():
            tokens = line.split()
            # skip header and blank lines
            if len(tokens) < 2 or tokens[0] == "Module":
                continue
            try:
                result[tokens[0]] = float(tokens[1])
            except ValueError:
                pass
    except Exception:
        pass
    return result


def _gaussian(x, amp, mu, sigma):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _fit_histogram(counts: np.ndarray, edges: np.ndarray
                   ) -> Tuple[float, float, float]:
    """Return (peak_MeV, sigma_MeV, chi2_ndf). Falls back to moment estimates."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw      = edges[1] - edges[0]
    total   = counts.sum()
    if total < 5:
        return 0.0, 0.0, 0.0

    # initial guess from histogram moments
    mu0    = (centers * counts).sum() / total
    sig0   = max(math.sqrt(((centers - mu0) ** 2 * counts).sum() / total), bw)
    amp0   = counts.max()

    # restrict fit window to ±2.5 sigma around peak
    lo, hi = mu0 - 2.5 * sig0, mu0 + 2.5 * sig0
    mask   = (centers >= lo) & (centers <= hi) & (counts > 0)
    if mask.sum() < 5:
        return mu0, sig0, 0.0

    xf, yf = centers[mask], counts[mask].astype(float)

    if HAS_SCIPY:
        try:
            popt, _ = curve_fit(_gaussian, xf, yf,
                                p0=[amp0, mu0, sig0],
                                bounds=([0, lo, bw],
                                        [amp0 * 3, hi, hi - lo]),
                                maxfev=3000)
            amp_f, mu_f, sig_f = popt
            resid = yf - _gaussian(xf, *popt)
            chi2  = (resid ** 2 / np.maximum(yf, 1)).sum()
            ndf   = max(len(xf) - 3, 1)
            return mu_f, abs(sig_f), chi2 / ndf
        except Exception:
            pass

    return mu0, sig0, 0.0


def compute_metrics(data: IterData) -> None:
    """Fill *data.metrics* from ROOT histograms + JSON factors."""
    factors = data.factors

    # --- histogram-based metrics ---
    if HAS_UPROOT and data.root_path is not None:
        with uproot.open(str(data.root_path)) as f:
            # collect available module_energy histogram names
            me_keys = {k.split(";")[0].replace("module_energy/", ""): k
                       for k in f.keys() if k.startswith("module_energy/")}

            for hname, raw_key in me_keys.items():
                # hname like "h_W1" -> module name "W1"
                mod_name = hname[len("h_"):]
                counts, edges = f[raw_key].to_numpy()
                data._hist_cache[mod_name] = (counts, edges)
                mm = data.metrics.setdefault(mod_name, ModuleMetrics(mod_name))
                mm.stats = float(counts.sum())
                mm.peak, mm.sigma, mm.chi2 = _fit_histogram(counts, edges)

    # --- JSON factor data ---
    for name, entry in factors.items():
        mm = data.metrics.setdefault(name, ModuleMetrics(name))
        mm.factor = float(entry.get("factor") or 1.0)
        # Prefer ExpectedPeak from .dat (per-run, per-module value);
        # fall back to base_energy in JSON only if .dat was not loaded.
        if name in data.expected_peaks:
            mm.base_energy = data.expected_peaks[name]
        else:
            mm.base_energy = float(entry.get("base_energy") or 0.0)

    # --- compute delta_E = (peak - base_energy) / base_energy ---
    # handled on-demand in _refresh_map()


# =============================================================================
#  Map widget
# =============================================================================

class CalibMapWidget(HyCalMapWidget):
    """HyCalMapWidget specialised for calibration quality display."""

    MODE_DELTA_E = "ΔE/E"
    MODE_SIGMA   = "σ (MeV)"
    MODE_CHI2    = "χ²/ndf"
    MODE_STATS   = "Statistics"
    MODE_FACTOR  = "Calib Factor"

    def __init__(self, parent=None):
        super().__init__(parent, enable_zoom_pan=True, min_size=(460, 460))
        self._label = ""

    def set_map_label(self, label: str):
        self._label = label

    def _colorbar_center_text(self) -> str:
        return self._label or PALETTE_NAMES[self._palette_idx]

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        base = name if v is None else f"{name}: {self._fmt_value(v)}"
        return base


# =============================================================================
#  Matplotlib canvas
# =============================================================================

class MplCanvas(FigureCanvas):
    def __init__(self, width=4, height=3, dpi=90):
        fc = THEME.CANVAS
        self.fig = Figure(figsize=(width, height), dpi=dpi,
                          facecolor=fc, tight_layout=True)
        self.ax  = self.fig.add_subplot(111)
        self.ax.set_facecolor(fc)
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

    def _style_ax(self, title="", xlabel="", ylabel=""):
        ax = self.ax
        ax.set_title(title, color=THEME.TEXT, fontsize=12)
        ax.set_xlabel(xlabel, color=THEME.TEXT_DIM, fontsize=10)
        ax.set_ylabel(ylabel, color=THEME.TEXT_DIM, fontsize=10)
        ax.tick_params(colors=THEME.TEXT_DIM, labelsize=9)
        for sp in ax.spines.values():
            sp.set_edgecolor(THEME.BORDER)

    def clear_ax(self):
        self.ax.cla()
        self.ax.set_facecolor(THEME.CANVAS)

    def redraw(self):
        self.fig.canvas.draw_idle()


# =============================================================================
#  Module detail panel
# =============================================================================

class ModuleDetailPanel(QWidget):
    refitApplied = pyqtSignal(str)   # emitted with module name after Apply+Save

    def __init__(self, parent=None):
        super().__init__(parent)
        self._root_path:       Optional[Path]          = None
        self._mm:              Optional[ModuleMetrics] = None
        self._iter_data:       Optional[IterData]      = None
        self._cur_module_name: str   = ""
        self._refit_peak:      float = 0.0
        self._refit_sigma:     float = 0.0
        self._refit_chi2:      float = 0.0
        self._refit_new_factor:float = 0.0
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        self._name_lbl = QLabel("← click a module on the map")
        self._name_lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        self._name_lbl.setStyleSheet(f"color: {THEME.ACCENT};")
        lay.addWidget(self._name_lbl)

        self._info_lbl = QLabel("")
        self._info_lbl.setFont(QFont("Monospace", 11))
        self._info_lbl.setStyleSheet(f"color: {THEME.TEXT};")
        self._info_lbl.setWordWrap(True)
        lay.addWidget(self._info_lbl)

        self._canvas = MplCanvas(width=4, height=3)
        lay.addWidget(self._canvas)

        # ── Manual re-fit controls ──────────────────────────────────────
        rf_box = QGroupBox("Manual Gauss Re-fit")
        rf_lay = QVBoxLayout(rf_box)
        rf_lay.setContentsMargins(6, 12, 6, 4)
        rf_lay.setSpacing(4)

        range_row = QHBoxLayout()
        range_row.setSpacing(4)
        range_row.addWidget(QLabel("Fit range:"))
        self._rf_xmin = QLineEdit()
        self._rf_xmin.setFixedWidth(72)
        self._rf_xmin.setPlaceholderText("Xmin")
        range_row.addWidget(self._rf_xmin)
        range_row.addWidget(QLabel("–"))
        self._rf_xmax = QLineEdit()
        self._rf_xmax.setFixedWidth(72)
        self._rf_xmax.setPlaceholderText("Xmax")
        range_row.addWidget(self._rf_xmax)
        range_row.addSpacing(16)
        range_row.addWidget(QLabel("Init μ₀:"))
        self._rf_mu0 = QLineEdit()
        self._rf_mu0.setFixedWidth(72)
        self._rf_mu0.setPlaceholderText("μ₀ (MeV)")
        range_row.addWidget(self._rf_mu0)
        range_row.addStretch()
        rf_lay.addLayout(range_row)

        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        self._refit_btn = QPushButton("Run Refit")
        self._refit_btn.setEnabled(False)
        self._refit_btn.clicked.connect(self._do_refit)
        btn_row.addWidget(self._refit_btn)
        self._mean_btn = QPushButton("Use Mean (χ²=1)")
        self._mean_btn.setEnabled(False)
        self._mean_btn.setToolTip(
            "Compute weighted mean/σ from histogram (in fit range if set)\n"
            "and use as peak value; χ²/ndf is fixed at 1.0")
        self._mean_btn.clicked.connect(self._do_use_mean)
        btn_row.addWidget(self._mean_btn)
        self._apply_btn = QPushButton("Apply && Save JSON")
        self._apply_btn.setEnabled(False)
        self._apply_btn.clicked.connect(self._apply_and_save)
        btn_row.addWidget(self._apply_btn)
        btn_row.addStretch()
        rf_lay.addLayout(btn_row)

        self._refit_status = QLabel("")
        self._refit_status.setFont(QFont("Monospace", 10))
        self._refit_status.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        self._refit_status.setWordWrap(True)
        rf_lay.addWidget(self._refit_status)

        lay.addWidget(rf_box)

    def set_root_path(self, path: Optional[Path]):
        self._root_path = path

    def show_module(self, name: str, mm: Optional[ModuleMetrics]):
        self._mm = mm
        self._cur_module_name  = name
        self._refit_peak       = 0.0
        self._refit_new_factor = 0.0
        self._apply_btn.setEnabled(False)
        self._refit_status.setText("")
        self._refit_btn.setEnabled(bool(name) and HAS_SCIPY and HAS_UPROOT)
        self._mean_btn.setEnabled(bool(name) and HAS_UPROOT)

        # Pre-populate fit range from current fit result
        if mm and mm.peak > 0 and mm.sigma > 0:
            self._rf_xmin.setText(f"{mm.peak - 3 * mm.sigma:.1f}")
            self._rf_xmax.setText(f"{mm.peak + 3 * mm.sigma:.1f}")
            self._rf_mu0.setText(f"{mm.peak:.1f}")
        elif mm and mm.base_energy > 0:
            self._rf_mu0.setText(f"{mm.base_energy:.1f}")

        self._name_lbl.setText(f"Module: {name}")
        if mm:
            delta_e = ""
            if mm.base_energy > 0 and mm.peak > 0:
                de = (mm.peak - mm.base_energy) / mm.base_energy * 100
                delta_e = f"  ΔE/E={de:+.2f}%"
            self._info_lbl.setText(
                f"Factor: {mm.factor:.4f}   Base energy: {mm.base_energy:.1f} MeV\n"
                f"Events: {mm.stats:.0f}   "
                f"Peak: {mm.peak:.1f} MeV   σ: {mm.sigma:.1f} MeV   "
                f"χ²/ndf: {mm.chi2:.3f}{delta_e}"
            )
        else:
            self._info_lbl.setText("")

        self._canvas.clear_ax()
        ax = self._canvas.ax

        hist_data = self._read_module_hist(name)
        if hist_data is not None:
            counts, edges = hist_data
            centers = 0.5 * (edges[:-1] + edges[1:])
            # crop to non-zero range + some margin
            nz = np.where(counts > 0)[0]
            if nz.size:
                lo_i = max(0, nz[0] - 5)
                hi_i = min(len(counts), nz[-1] + 6)
                counts  = counts[lo_i:hi_i]
                centers = centers[lo_i:hi_i]
                edges_c = edges[lo_i:hi_i + 1]
            else:
                edges_c = edges

            ax.bar(centers, counts, width=np.diff(edges_c),
                   color=THEME.ACCENT, alpha=0.75, linewidth=0)

            # overlay Gaussian fit
            if mm and mm.peak > 0 and mm.sigma > 0 and HAS_SCIPY:
                x  = np.linspace(edges_c[0], edges_c[-1], 500)
                bw = edges_c[1] - edges_c[0] if len(edges_c) > 1 else 1.0
                area = counts.sum() * bw
                amp  = area / (mm.sigma * math.sqrt(2 * math.pi))
                gaus = amp * np.exp(-0.5 * ((x - mm.peak) / mm.sigma) ** 2)
                ax.plot(x, gaus, color=THEME.WARN, linewidth=1.5,
                        label=f"Gauss fit  μ={mm.peak:.0f}")

            if mm and mm.base_energy > 0:
                ax.axvline(mm.base_energy, color=THEME.SUCCESS,
                           linewidth=1.2, linestyle="--",
                           label=f"Expected {mm.base_energy:.0f} MeV")

            ax.legend(fontsize=9, labelcolor=THEME.TEXT,
                      facecolor=THEME.PANEL, edgecolor=THEME.BORDER)
        else:
            msg = ("uproot not installed\npip install uproot"
                   if not HAS_UPROOT else "No histogram found in ROOT file")
            ax.text(0.5, 0.5, msg, transform=ax.transAxes,
                    ha="center", va="center",
                    color=THEME.TEXT_DIM, fontsize=9)

        self._canvas._style_ax(f"{name}  energy histogram",
                                "Energy (MeV)", "Counts")
        self._canvas.redraw()

    def _read_module_hist(self, name: str
                          ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        # Fast path: return in-memory cache populated during compute_metrics
        if self._iter_data is not None and name in self._iter_data._hist_cache:
            return self._iter_data._hist_cache[name]
        if not HAS_UPROOT or self._root_path is None:
            return None
        key = f"module_energy/h_{name}"
        try:
            with uproot.open(str(self._root_path)) as f:
                avail = {k.split(";")[0] for k in f.keys()}
                if key not in avail:
                    return None
                result = f[key].to_numpy()
                if self._iter_data is not None:
                    self._iter_data._hist_cache[name] = result
                return result
        except Exception:
            return None

    # ──────────────────────────────────────────────────────────────────────

    def set_iter_data(self, data: Optional["IterData"]) -> None:
        """Store reference to current IterData (needed for JSON write-back)."""
        self._iter_data = data

    def _do_use_mean(self) -> None:
        """Use histogram weighted mean as peak value; chi2/ndf is fixed at 1.0."""
        name = self._cur_module_name
        if not name or self._root_path is None:
            self._refit_status.setText("No module selected")
            return

        hist_data = self._read_module_hist(name)
        if hist_data is None:
            self._refit_status.setText("Histogram not available in ROOT file")
            return

        counts, edges = hist_data
        centers = 0.5 * (edges[:-1] + edges[1:])

        # Apply fit range if provided
        xmin_txt = self._rf_xmin.text().strip()
        xmax_txt = self._rf_xmax.text().strip()
        try:
            xmin = float(xmin_txt) if xmin_txt else None
        except ValueError:
            xmin = None
        try:
            xmax = float(xmax_txt) if xmax_txt else None
        except ValueError:
            xmax = None

        mask = counts > 0
        if xmin is not None:
            mask &= centers >= xmin
        if xmax is not None:
            mask &= centers <= xmax

        c_sel = counts[mask]
        x_sel = centers[mask]
        total = c_sel.sum()
        if total < 1:
            self._refit_status.setText("No data in selected range")
            return

        mu0  = float((x_sel * c_sel).sum() / total)
        var  = float(((x_sel - mu0) ** 2 * c_sel).sum() / total)
        bw   = float(edges[1] - edges[0]) if len(edges) > 1 else 1.0
        sig0 = math.sqrt(max(var, bw ** 2 / 12.0))

        old_peak   = self._mm.peak   if (self._mm and self._mm.peak > 0.0) else 0.0
        old_factor = self._mm.factor if self._mm else 1.0
        new_factor = old_factor * mu0 / old_peak if old_peak > 0.0 else 0.0

        self._refit_peak       = mu0
        self._refit_sigma      = sig0
        self._refit_chi2       = 1.0
        self._refit_new_factor = new_factor

        range_note = ""
        if xmin is not None or xmax is not None:
            lo_s = f"{xmin:.1f}" if xmin is not None else "-inf"
            hi_s = f"{xmax:.1f}" if xmax is not None else "+inf"
            range_note = f"  [range: {lo_s} – {hi_s}]"

        msg = (
            f"Mean:  μ = {mu0:.2f} MeV     "
            f"σ = {sig0:.2f} MeV     χ²/ndf = 1.000  (fixed){range_note}"
        )
        if new_factor > 0.0:
            msg += f"\nNew factor: {new_factor:.5f}   (was {old_factor:.5f})"
        self._refit_status.setText(msg)
        self._apply_btn.setEnabled(new_factor > 0.0)
        self._redraw_refit_overlay(name, counts, edges, mu0, sig0)

    def _do_refit(self) -> None:
        """Gaussian refit using user-specified range and initial peak center."""
        if not HAS_SCIPY or not HAS_UPROOT:
            self._refit_status.setText("scipy and uproot required")
            return
        name = self._cur_module_name
        if not name or self._root_path is None:
            self._refit_status.setText("No module selected")
            return
        try:
            xmin = float(self._rf_xmin.text())
            xmax = float(self._rf_xmax.text())
            mu0  = float(self._rf_mu0.text())
        except ValueError:
            self._refit_status.setText("Invalid input — enter numeric values in all three fields")
            return
        if xmin >= xmax:
            self._refit_status.setText("Fit range error: Xmin must be < Xmax")
            return

        hist_data = self._read_module_hist(name)
        if hist_data is None:
            self._refit_status.setText("Histogram not available in ROOT file")
            return

        counts, edges = hist_data
        centers = 0.5 * (edges[:-1] + edges[1:])
        bw = float(edges[1] - edges[0]) if len(edges) > 1 else 1.0

        mask = (centers >= xmin) & (centers <= xmax) & (counts > 0)
        if mask.sum() < 4:
            self._refit_status.setText(
                f"Only {mask.sum()} non-empty bins in fit range — widen the range")
            return

        xf   = centers[mask]
        yf   = counts[mask].astype(float)
        amp0 = float(yf.max())
        sig0 = max((xmax - xmin) / 5.0, bw)

        try:
            popt, pcov = curve_fit(
                _gaussian, xf, yf,
                p0=[amp0, mu0, sig0],
                bounds=([0.0, xmin, bw], [amp0 * 3.0, xmax, xmax - xmin]),
                maxfev=6000,
            )
            amp_f, mu_f, sig_f = popt
            resid    = yf - _gaussian(xf, *popt)
            chi2     = float((resid ** 2 / np.maximum(yf, 1.0)).sum())
            ndf      = max(len(xf) - 3, 1)
            chi2_ndf = chi2 / ndf

            self._refit_peak  = float(mu_f)
            self._refit_sigma = float(abs(sig_f))
            self._refit_chi2  = chi2_ndf

            old_peak   = self._mm.peak   if (self._mm and self._mm.peak > 0.0) else 0.0
            old_factor = self._mm.factor if self._mm else 1.0
            self._refit_new_factor = (
                old_factor * mu_f / old_peak if old_peak > 0.0 else 0.0
            )

            msg = (
                f"Refit:  μ = {mu_f:.2f} MeV     "
                f"σ = {abs(sig_f):.2f} MeV     χ²/ndf = {chi2_ndf:.3f}"
            )
            if self._refit_new_factor > 0.0:
                msg += (
                    f"\nNew factor: {self._refit_new_factor:.5f}"
                    f"   (was {old_factor:.5f})"
                )
            self._refit_status.setText(msg)
            self._apply_btn.setEnabled(self._refit_new_factor > 0.0)
            self._redraw_refit_overlay(
                name, counts, edges, float(mu_f), float(abs(sig_f)))

        except Exception as exc:
            self._refit_status.setText(f"Fit failed: {exc}")
            self._refit_peak       = 0.0
            self._refit_new_factor = 0.0

    def _redraw_refit_overlay(
            self,
            name: str,
            counts: np.ndarray,
            edges: np.ndarray,
            mu_f: float,
            sig_f: float,
    ) -> None:
        """Redraw histogram with old (dashed) + new (solid) Gaussian overlaid."""
        self._canvas.clear_ax()
        ax = self._canvas.ax

        centers = 0.5 * (edges[:-1] + edges[1:])
        nz = np.where(counts > 0)[0]
        if nz.size:
            lo_i = max(0, nz[0] - 5)
            hi_i = min(len(counts), nz[-1] + 6)
            c, ce, ee = counts[lo_i:hi_i], centers[lo_i:hi_i], edges[lo_i:hi_i + 1]
        else:
            c, ce, ee = counts, centers, edges

        bw = float(ee[1] - ee[0]) if len(ee) > 1 else 1.0
        ax.bar(ce, c, width=bw, color=THEME.ACCENT, alpha=0.6, linewidth=0)

        x = np.linspace(float(ee[0]), float(ee[-1]), 600)
        # original fit — dim dashed line
        mm = self._mm
        if mm and mm.peak > 0.0 and mm.sigma > 0.0:
            area = float(c.sum()) * bw
            amp  = area / (mm.sigma * math.sqrt(2 * math.pi))
            ax.plot(x, amp * np.exp(-0.5 * ((x - mm.peak) / mm.sigma) ** 2),
                    color=THEME.TEXT_DIM, linewidth=1.0, linestyle="--",
                    label=f"Original  μ={mm.peak:.0f}")
        # new refit — highlighted solid line
        area = float(c.sum()) * bw
        amp  = area / (sig_f * math.sqrt(2 * math.pi))
        ax.plot(x, amp * np.exp(-0.5 * ((x - mu_f) / sig_f) ** 2),
                color=THEME.WARN, linewidth=2.0,
                label=f"Refit  μ={mu_f:.0f}")

        if mm and mm.base_energy > 0.0:
            ax.axvline(mm.base_energy, color=THEME.SUCCESS,
                       linewidth=1.2, linestyle="--",
                       label=f"Expected {mm.base_energy:.0f} MeV")

        ax.legend(fontsize=9, labelcolor=THEME.TEXT,
                  facecolor=THEME.PANEL, edgecolor=THEME.BORDER)
        self._canvas._style_ax(f"{name}  energy histogram (refit)",
                                "Energy (MeV)", "Counts")
        self._canvas.redraw()

    def _apply_and_save(self) -> None:
        """Update in-memory IterData metrics and write calibration JSON to disk."""
        name = self._cur_module_name
        if (
            not name
            or self._iter_data is None
            or self._refit_peak <= 0.0
            or self._refit_new_factor <= 0.0
        ):
            return
        mm = self._iter_data.metrics.get(name)
        if mm is None:
            return

        old_factor = mm.factor
        mm.peak   = self._refit_peak
        mm.sigma  = self._refit_sigma
        mm.chi2   = self._refit_chi2
        mm.factor = self._refit_new_factor

        # keep factors dict in sync so JSON reflects the update
        entry = self._iter_data.factors.get(name)
        if entry is not None:
            entry["factor"] = self._refit_new_factor

        # write JSON to disk
        jpath = self._iter_data.json_path
        if jpath is not None:
            try:
                entries = list(self._iter_data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                save_note = f"  →  saved {jpath.name}"
            except Exception as exc:
                self._refit_status.setText(f"Save failed: {exc}")
                return
        else:
            save_note = "  (JSON path not found — not saved to disk)"

        # write .dat to disk
        dpath = self._iter_data.dat_path
        dat_note = ""
        if dpath is not None and dpath.is_file():
            try:
                lines = dpath.read_text().splitlines(keepends=True)
                new_lines = []
                updated = False
                for line in lines:
                    tokens = line.split()
                    if tokens and tokens[0] == name:
                        # reconstruct line preserving column widths
                        ratio = (mm.base_energy / mm.peak
                                 if mm.peak > 0.0 else 0.0)
                        new_line = (
                            f"{name:<16}{mm.base_energy:<16.6g}"
                            f"{mm.peak:<16.6g}{ratio:<16.6g}"
                            f"{mm.sigma:<16.6g}{mm.chi2:<16.6g}\n"
                        )
                        new_lines.append(new_line)
                        updated = True
                    else:
                        new_lines.append(line)
                if updated:
                    dpath.write_text("".join(new_lines))
                    dat_note = f", {dpath.name}"
                else:
                    dat_note = f" (module {name!r} not found in .dat)"
            except Exception as exc:
                dat_note = f" (.dat write failed: {exc})"
        elif dpath is not None:
            dat_note = " (.dat file not found)"

        self._mm = mm
        delta_e = ""
        if mm.base_energy > 0 and mm.peak > 0:
            de = (mm.peak - mm.base_energy) / mm.base_energy * 100
            delta_e = f"  ΔE/E={de:+.2f}%"
        self._info_lbl.setText(
            f"Factor: {mm.factor:.5f}  (was {old_factor:.5f})   "
            f"Base energy: {mm.base_energy:.1f} MeV\n"
            f"Events: {mm.stats:.0f}   "
            f"Peak: {mm.peak:.1f} MeV   σ: {mm.sigma:.1f} MeV   "
            f"χ²/ndf: {mm.chi2:.3f}{delta_e}"
        )
        first_line = self._refit_status.text().split("\n")[0]
        self._refit_status.setText(
            first_line + f"\nApplied.{save_note}{dat_note}")
        self._apply_btn.setEnabled(False)
        self.refitApplied.emit(name)


# =============================================================================
#  Global stats panel
# =============================================================================

class GlobalStatsPanel(QWidget):

    _RIGHT_VIEWS = ["Energy vs θ", "One Cluster Energy", "Resolution (σ)"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self._root_path: Optional[Path] = None
        self._right_view_idx = 0
        self._cached_hists: dict = {}
        self._build_ui()

    def _build_ui(self):
        lay = QHBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.setSpacing(4)
        self._cv_ratio  = MplCanvas(width=3, height=2.4)
        self._cv_etheta = MplCanvas(width=3, height=2.4)
        self._etheta_cbar = None
        lay.addWidget(self._cv_ratio)

        # right panel: toggle button + canvas
        right_box = QWidget()
        right_lay = QVBoxLayout(right_box)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(2)
        self._right_toggle = QPushButton(self._RIGHT_VIEWS[self._right_view_idx])
        self._right_toggle.setToolTip("Click to switch between views")
        self._right_toggle.clicked.connect(self._cycle_right_view)
        right_lay.addWidget(self._right_toggle)
        right_lay.addWidget(self._cv_etheta)
        lay.addWidget(right_box)

    def _cycle_right_view(self):
        self._right_view_idx = (self._right_view_idx + 1) % len(self._RIGHT_VIEWS)
        self._right_toggle.setText(self._RIGHT_VIEWS[self._right_view_idx])
        self._draw_right()

    def set_root_path(self, path: Optional[Path]):
        self._root_path = path
        self._cached_hists.clear()
        if HAS_UPROOT and path is not None:
            try:
                with uproot.open(str(path)) as f:
                    avail = {k.split(";")[0] for k in f.keys()}
                    for key in ("ratio_all", "h2_energy_theta",
                                "one_cluster_energy", "recon_sigma"):
                        if key in avail:
                            try:
                                self._cached_hists[key] = f[key].to_numpy()
                            except Exception:
                                pass
            except Exception:
                pass
        self._draw_ratio()
        self._draw_right()

    def _draw_right(self):
        idx = self._right_view_idx
        if idx == 0:
            self._draw_etheta()
        elif idx == 1:
            self._draw_one_cluster_energy()
        else:
            self._draw_resolution()

    def _open_root(self):
        if not HAS_UPROOT or self._root_path is None:
            return None
        try:
            return uproot.open(str(self._root_path))
        except Exception:
            return None

    def _draw_ratio(self):
        ax = self._cv_ratio.ax
        self._cv_ratio.clear_ax()
        data = self._cached_hists.get("ratio_all")
        if data is not None:
            try:
                counts, edges = data
                centers = 0.5 * (edges[:-1] + edges[1:])
                nz = counts > 0
                if nz.any():
                    ax.bar(centers[nz], counts[nz],
                           width=(edges[1] - edges[0]),
                           color=THEME.SUCCESS, alpha=0.8, linewidth=0)
                ax.axvline(1.0, color=THEME.WARN, linewidth=1.2,
                           linestyle="--", label="ratio = 1")
                lo = max(0, centers[nz][0] - 0.3) if nz.any() else 0
                hi = min(4, centers[nz][-1] + 0.3) if nz.any() else 4
                ax.set_xlim(lo, hi)
                ax.legend(fontsize=9, labelcolor=THEME.TEXT,
                          facecolor=THEME.PANEL, edgecolor=THEME.BORDER)
            except Exception:
                pass
        else:
            ax.text(0.5, 0.5, "ratio_all\nnot found",
                    transform=ax.transAxes, ha="center", va="center",
                    color=THEME.TEXT_DIM, fontsize=9)
        self._cv_ratio._style_ax("Peak ratio distribution",
                                  "Ratio (exp / meas)", "Modules")
        self._cv_ratio.redraw()

    def _draw_etheta(self):
        ax = self._cv_etheta.ax
        # remove previous colorbar before clearing axes
        if self._etheta_cbar is not None:
            self._etheta_cbar.remove()
            self._etheta_cbar = None
        self._cv_etheta.clear_ax()
        data = self._cached_hists.get("h2_energy_theta")
        if data is not None:
            try:
                counts, xedges, yedges = data
                pcm = ax.pcolormesh(xedges, yedges, counts.T,
                                    cmap="viridis", shading="flat")
                self._etheta_cbar = self._cv_etheta.fig.colorbar(pcm, ax=ax, pad=0.01)
            except Exception:
                ax.text(0.5, 0.5, "h2_energy_theta\nnot found",
                        transform=ax.transAxes, ha="center", va="center",
                        color=THEME.TEXT_DIM, fontsize=9)
        self._cv_etheta._style_ax("Energy vs θ", "θ (deg)", "Energy (MeV)")
        self._cv_etheta.redraw()

    def _draw_1d_hist(self, hist_key: str, title: str, xlabel: str,
                      color: str, axvline: Optional[float] = None):
        """Generic helper: draw a 1D histogram from the ROOT file."""
        if self._etheta_cbar is not None:
            self._etheta_cbar.remove()
            self._etheta_cbar = None
        ax = self._cv_etheta.ax
        self._cv_etheta.clear_ax()
        data = self._cached_hists.get(hist_key)
        if data is not None:
            try:
                counts, edges = data
                centers = 0.5 * (edges[:-1] + edges[1:])
                nz = counts > 0
                if nz.any():
                    ax.bar(centers, counts, width=(edges[1] - edges[0]),
                           color=color, alpha=0.85, linewidth=0)
                    lo = centers[nz][0] - (edges[1] - edges[0])
                    hi = centers[nz][-1] + (edges[1] - edges[0])
                    ax.set_xlim(lo, hi)
                if axvline is not None:
                    ax.axvline(axvline, color=THEME.WARN,
                               linewidth=1.2, linestyle="--")
            except Exception:
                ax.text(0.5, 0.5, f"{hist_key}\nnot found",
                        transform=ax.transAxes, ha="center", va="center",
                        color=THEME.TEXT_DIM, fontsize=9)
        self._cv_etheta._style_ax(title, xlabel, "Counts")
        self._cv_etheta.redraw()

    def _draw_one_cluster_energy(self):
        self._draw_1d_hist("one_cluster_energy",
                           "One Cluster Energy", "Energy (MeV)",
                           color=THEME.ACCENT)

    def _draw_resolution(self):
        self._draw_1d_hist("recon_sigma",
                           "Resolution (σ)", "σ (MeV)",
                           color=THEME.SUCCESS)


# =============================================================================
#  Background worker
# =============================================================================

class _MetricsWorker(QThread):
    """Run compute_metrics() in a background thread to keep the UI responsive."""
    finished = pyqtSignal(object)   # emits the IterData when done

    def __init__(self, data: "IterData", parent=None):
        super().__init__(parent)
        self._data = data

    def run(self) -> None:
        compute_metrics(self._data)
        self.finished.emit(self._data)


# =============================================================================
#  Main window
# =============================================================================

class EpCalibViewerWindow(QMainWindow):

    _MODES = [
        CalibMapWidget.MODE_DELTA_E,
        CalibMapWidget.MODE_STATS,
        CalibMapWidget.MODE_SIGMA,
        CalibMapWidget.MODE_CHI2,
        CalibMapWidget.MODE_FACTOR,
    ]

    def __init__(self, initial_dir: Optional[Path] = None):
        super().__init__()
        self._calib_dir: Optional[Path] = None
        self._scan_data: Dict[str, Dict[int, IterData]] = {}
        self._modules:   List[Module] = []
        self._cur_data:  Optional[IterData] = None
        self._map_mode = CalibMapWidget.MODE_DELTA_E
        self._worker:    Optional[_MetricsWorker] = None

        self._build_ui()
        self.setWindowTitle("epCalib Viewer")
        self.resize(1380, 900)
        self._apply_stylesheet()

        if MODULES_JSON:
            self._modules = load_modules(MODULES_JSON)
            self._map.set_modules(self._modules)

        if initial_dir:
            self._load_calib_dir(initial_dir)

    # ── stylesheet ────────────────────────────────────────────────────────────

    def _apply_stylesheet(self):
        t = THEME
        self.setStyleSheet(f"""
            QMainWindow, QWidget {{ background: {t.BG}; color: {t.TEXT}; }}
            QComboBox {{ background: {t.PANEL}; border: 1px solid {t.BORDER};
                         color: {t.TEXT}; padding: 2px 6px; }}
            QComboBox QAbstractItemView {{ background: {t.PANEL}; color: {t.TEXT}; }}
            QPushButton {{ background: {t.BUTTON}; border: 1px solid {t.BORDER};
                           color: {t.TEXT}; padding: 4px 10px; border-radius: 4px; }}
            QPushButton:hover   {{ background: {t.BUTTON_HOVER}; }}
            QPushButton:checked {{ background: {t.ACCENT_STRONG};
                                   border: 1px solid {t.ACCENT_BORDER}; }}
            QGroupBox {{ color: {t.TEXT_DIM}; font-size: 12px;
                         border: 1px solid {t.BORDER}; margin-top: 6px; }}
            QGroupBox::title {{ subcontrol-origin: margin; left: 8px; }}
            QLineEdit {{ background: {t.PANEL}; border: 1px solid {t.BORDER};
                         color: {t.TEXT}; padding: 2px 4px; border-radius: 3px; }}
            QLineEdit:focus {{ border: 1px solid {t.ACCENT_BORDER}; }}
            QCheckBox {{ color: {t.TEXT}; spacing: 5px; }}
            QCheckBox::indicator {{ width: 14px; height: 14px;
                                    background: {t.PANEL}; border: 1px solid {t.BORDER};
                                    border-radius: 3px; }}
            QCheckBox::indicator:checked {{ background: {t.ACCENT_STRONG};
                                            border: 1px solid {t.ACCENT_BORDER}; }}
            QSplitter::handle {{ background: {t.PANEL}; }}
            QStatusBar {{ color: {t.TEXT_DIM}; font-size: 11px; }}
        """)

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(5)

        # ── toolbar ───────────────────────────────────────────────────────────
        top = QHBoxLayout()
        top.setSpacing(8)
        title = QLabel("epCalib Viewer")
        title.setFont(QFont("Monospace", 15, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {THEME.ACCENT};")
        top.addWidget(title)
        self._dir_btn = QPushButton("Open Physics_calib dir…")
        self._dir_btn.clicked.connect(self._browse_dir)
        top.addWidget(self._dir_btn)
        self._refresh_btn = QPushButton("⟳ Refresh")
        self._refresh_btn.setToolTip("Rescan the current directory for new/updated files")
        self._refresh_btn.clicked.connect(self._refresh_dir)
        self._refresh_btn.setEnabled(False)
        top.addWidget(self._refresh_btn)
        self._dir_lbl = QLabel("(no directory)")
        self._dir_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 12px;")
        top.addWidget(self._dir_lbl)
        top.addStretch()
        top.addWidget(QLabel("Run:"))
        self._run_cb = QComboBox()
        self._run_cb.setMinimumWidth(100)
        self._run_cb.currentTextChanged.connect(self._on_run_changed)
        top.addWidget(self._run_cb)
        top.addWidget(QLabel("Iteration:"))
        self._iter_cb = QComboBox()
        self._iter_cb.setMinimumWidth(70)
        self._iter_cb.currentTextChanged.connect(self._on_iter_changed)
        top.addWidget(self._iter_cb)
        root.addLayout(top)

        # ── map mode buttons ──────────────────────────────────────────────────
        mode_bar = QHBoxLayout()
        mode_bar.addWidget(QLabel("Map:"))
        self._mode_grp = QButtonGroup(self)
        self._mode_grp.setExclusive(True)
        for mode in self._MODES:
            btn = QPushButton(mode)
            btn.setCheckable(True)
            if mode == self._map_mode:
                btn.setChecked(True)
            btn.clicked.connect(lambda _, m=mode: self._set_map_mode(m))
            self._mode_grp.addButton(btn)
            mode_bar.addWidget(btn)
        mode_bar.addStretch()
        self._bad_lbl = QLabel("")
        self._bad_lbl.setStyleSheet(f"color: {THEME.DANGER}; font-size: 11px;")
        mode_bar.addWidget(self._bad_lbl)
        root.addLayout(mode_bar)

        # ── Z-axis range / log controls ───────────────────────────────────────
        z_bar = QHBoxLayout()
        z_bar.setSpacing(6)

        self._logz_btn = QCheckBox("Log Z")
        self._logz_btn.setChecked(False)
        self._logz_btn.toggled.connect(self._toggle_logz)
        z_bar.addWidget(self._logz_btn)

        z_bar.addWidget(QLabel("Zmin:"))
        self._zmin_edit = QLineEdit()
        self._zmin_edit.setFixedWidth(80)
        self._zmin_edit.setPlaceholderText("auto")
        self._zmin_edit.returnPressed.connect(self._apply_z_range)
        z_bar.addWidget(self._zmin_edit)

        z_bar.addWidget(QLabel("Zmax:"))
        self._zmax_edit = QLineEdit()
        self._zmax_edit.setFixedWidth(80)
        self._zmax_edit.setPlaceholderText("auto")
        self._zmax_edit.returnPressed.connect(self._apply_z_range)
        z_bar.addWidget(self._zmax_edit)

        self._zauto_btn = QPushButton("Auto range")
        self._zauto_btn.clicked.connect(self._do_auto_range)
        z_bar.addWidget(self._zauto_btn)

        z_bar.addStretch()
        root.addLayout(z_bar)

        # ── main splitter ─────────────────────────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(4)
        root.addWidget(splitter, stretch=1)

        # left: HyCal map
        map_box = QGroupBox("HyCal Map  (scroll to zoom · drag to pan · click module)")
        ml = QVBoxLayout(map_box)
        ml.setContentsMargins(2, 14, 2, 2)
        self._map = CalibMapWidget()
        self._map.moduleClicked.connect(self._on_module_clicked)
        self._map.moduleHovered.connect(self._on_module_hovered)
        ml.addWidget(self._map)
        splitter.addWidget(map_box)

        # right panel
        right = QSplitter(Qt.Orientation.Vertical)
        right.setHandleWidth(4)
        splitter.addWidget(right)

        detail_box = QGroupBox("Module Detail")
        dl = QVBoxLayout(detail_box)
        dl.setContentsMargins(2, 14, 2, 2)
        self._detail = ModuleDetailPanel()
        self._detail.refitApplied.connect(self._on_refit_applied)
        dl.addWidget(self._detail)
        right.addWidget(detail_box)

        stats_box = QGroupBox("Global Statistics")
        sl = QVBoxLayout(stats_box)
        sl.setContentsMargins(2, 14, 2, 2)
        self._stats = GlobalStatsPanel()
        sl.addWidget(self._stats)
        right.addWidget(stats_box)

        splitter.setSizes([580, 670])
        right.setSizes([490, 290])

        # status bar
        self._hover_lbl = QLabel("")
        self.statusBar().addWidget(self._hover_lbl)
        if not HAS_UPROOT:
            self.statusBar().showMessage(
                "uproot not installed — histograms unavailable.  "
                "Run: pip install uproot", 0)
        elif not HAS_SCIPY:
            self.statusBar().showMessage(
                "scipy not installed — Gaussian fitting disabled.  "
                "Run: pip install scipy", 0)

    # ── directory loading ─────────────────────────────────────────────────────

    def _browse_dir(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select Physics_calib directory",
            str(self._calib_dir or Path.cwd()))
        if path:
            self._load_calib_dir(Path(path))

    def _refresh_dir(self):
        if self._calib_dir:
            cur_run = self._run_cb.currentText()
            cur_iter = self._iter_cb.currentText()
            self._load_calib_dir(self._calib_dir)
            # try to restore previous run/iter selection
            if cur_run in self._scan_data:
                self._run_cb.setCurrentText(cur_run)
                self._on_run_changed(cur_run)
                if cur_iter:
                    self._iter_cb.setCurrentText(cur_iter)
                    self._on_iter_changed(cur_iter)

    def _load_calib_dir(self, path: Path):
        self._calib_dir = path
        self._dir_lbl.setText(str(path))
        self._scan_data = scan_calib_dir(path)
        self._refresh_btn.setEnabled(True)

        self._run_cb.blockSignals(True)
        self._run_cb.clear()
        for run in sorted(self._scan_data.keys()):
            self._run_cb.addItem(run)
        self._run_cb.blockSignals(False)

        if self._scan_data:
            self._run_cb.setCurrentIndex(0)
            self._on_run_changed(self._run_cb.currentText())

    # ── combo callbacks ───────────────────────────────────────────────────────

    def _on_run_changed(self, run: str):
        if run not in self._scan_data:
            return
        self._iter_cb.blockSignals(True)
        self._iter_cb.clear()
        for it in sorted(self._scan_data[run].keys()):
            self._iter_cb.addItem(str(it))
        self._iter_cb.blockSignals(False)
        if self._iter_cb.count():
            self._iter_cb.setCurrentIndex(self._iter_cb.count() - 1)
            self._on_iter_changed(self._iter_cb.currentText())

    def _on_iter_changed(self, it_str: str):
        if not it_str:
            return
        run = self._run_cb.currentText()
        it  = int(it_str)
        data = self._scan_data.get(run, {}).get(it)
        if data is None:
            return
        self._cur_data = data
        self._detail.set_root_path(data.root_path)
        self._detail.set_iter_data(data)
        self._stats.set_root_path(data.root_path)
        if not data.metrics:
            # stop any previous worker still running
            if self._worker is not None and self._worker.isRunning():
                self._worker.quit()
                self._worker.wait()
            self.statusBar().showMessage(
                f"Computing metrics for run {run} iter {it} …")
            self._worker = _MetricsWorker(data, self)
            self._worker.finished.connect(self._on_metrics_ready)
            self._worker.start()
        else:
            self._refresh_map()

    def _on_metrics_ready(self, data: IterData) -> None:
        """Called from background worker when compute_metrics finishes."""
        self.statusBar().clearMessage()
        if data is self._cur_data:
            self._refresh_map()
        self._worker = None

    # ── map mode ──────────────────────────────────────────────────────────────

    def _set_map_mode(self, mode: str):
        self._map_mode = mode
        self._do_auto_range()   # reset range when switching modes

    def _toggle_logz(self, on: bool):
        self._map.set_log_scale(on)
        self._map.update()

    def _apply_z_range(self):
        """Read Zmin/Zmax edits and apply to map; empty = auto."""
        try:
            zmin = float(self._zmin_edit.text())
        except ValueError:
            zmin = None
        try:
            zmax = float(self._zmax_edit.text())
        except ValueError:
            zmax = None
        if zmin is None or zmax is None:
            self._do_auto_range()
            return
        if zmin >= zmax:
            return
        self._map.set_range(zmin, zmax)
        self._map.update()

    def _do_auto_range(self):
        """Auto-range from current map values and update the edit boxes."""
        if self._cur_data is None:
            return
        self._refresh_map()
        vmin, vmax = self._map.auto_range()
        self._zmin_edit.setText(f"{vmin:.6g}")
        self._zmax_edit.setText(f"{vmax:.6g}")

    def _refresh_map(self):
        data = self._cur_data
        if data is None:
            return

        values: Dict[str, float] = {}
        mode = self._map_mode

        if mode == CalibMapWidget.MODE_DELTA_E:
            for name, mm in data.metrics.items():
                if mm.base_energy > 0 and mm.peak > 0:
                    values[name] = abs(mm.peak - mm.base_energy) / mm.base_energy
            label = "ΔE/E = |peak − expected| / expected"

        elif mode == CalibMapWidget.MODE_STATS:
            for name, mm in data.metrics.items():
                values[name] = mm.stats
            label = "Event count per module"

        elif mode == CalibMapWidget.MODE_SIGMA:
            for name, mm in data.metrics.items():
                if mm.sigma > 0:
                    values[name] = mm.sigma
            label = "Gaussian σ (MeV)"

        elif mode == CalibMapWidget.MODE_CHI2:
            for name, mm in data.metrics.items():
                if mm.chi2 > 0:
                    values[name] = mm.chi2
            label = "χ²/ndf  (Gaussian fit)"

        elif mode == CalibMapWidget.MODE_FACTOR:
            for name, mm in data.metrics.items():
                values[name] = mm.factor
            label = "Calibration factor"

        self._map.set_map_label(label)
        self._map.set_values(values)
        # only auto-range when both edit boxes are empty (user hasn't set limits)
        zmin_txt = self._zmin_edit.text().strip()
        zmax_txt = self._zmax_edit.text().strip()
        if values and (not zmin_txt or not zmax_txt):
            vmin, vmax = self._map.auto_range()
            self._zmin_edit.setText(f"{vmin:.6g}")
            self._zmax_edit.setText(f"{vmax:.6g}")
        elif values and zmin_txt and zmax_txt:
            # re-apply the user-defined range (values didn't change the limits)
            try:
                self._map.set_range(float(zmin_txt), float(zmax_txt))
            except ValueError:
                pass

        # count "bad" modules: chi2 > 2 or |ΔE/E| > 2%
        n_bad = sum(
            1 for mm in data.metrics.values()
            if mm.chi2 > 2.0 or (
                mm.base_energy > 0 and mm.peak > 0
                and abs(mm.peak - mm.base_energy) / mm.base_energy > 0.02
            )
        )
        total = len(data.metrics)
        self._bad_lbl.setText(
            f"  Flagged (χ²>2 or |ΔE/E|>2%): {n_bad}/{total}" if total else "")

    # ── module interaction ────────────────────────────────────────────────────

    def _on_module_clicked(self, name: str):
        if not name or self._cur_data is None:
            return
        mm = self._cur_data.metrics.get(name)
        self._detail.show_module(name, mm)

    def _on_module_hovered(self, name: str):
        if not name or self._cur_data is None:
            self._hover_lbl.setText("")
            return
        mm = self._cur_data.metrics.get(name)
        if mm:
            de = ""
            if mm.base_energy > 0 and mm.peak > 0:
                de = f"  ΔE/E={(mm.peak - mm.base_energy)/mm.base_energy*100:+.2f}%"
            self._hover_lbl.setText(
                f"{name}  |  factor={mm.factor:.4f}  "
                f"events={mm.stats:.0f}  "
                f"peak={mm.peak:.1f} MeV  σ={mm.sigma:.1f} MeV  "
                f"χ²/ndf={mm.chi2:.3f}{de}")
        else:
            self._hover_lbl.setText(name)

    def _on_refit_applied(self, _name: str) -> None:
        """Called after a manual refit is saved; refresh map colors."""
        self._refresh_map()


# =============================================================================
#  Entry point
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description="epCalib output GUI viewer")
    ap.add_argument("calib_dir", nargs="?", type=Path,
                    help="Path to Physics_calib directory")
    ap.add_argument("--theme", choices=available_themes(), default="dark")
    args = ap.parse_args()

    set_theme(args.theme)

    app = QApplication(sys.argv)
    app.setApplicationName("epCalib Viewer")
    win = EpCalibViewerWindow(initial_dir=args.calib_dir)
    apply_theme_palette(win)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
