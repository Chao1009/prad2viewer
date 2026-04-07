#!/usr/bin/env python3
"""
cosmic_eq_analyze — pass 1 of the offline cosmic-data gain-equalization workflow.

Reads one or more ROOT files produced by ``replay_rawdata -p`` (which contain
peak_height[nch][MAX_PEAKS] and peak_integral[nch][MAX_PEAKS] plus the FP
trigger word per event). All input files are chained together and treated as
one logical event stream. For each HyCal channel:

  1. Filter events: only those with the sum trigger bit set (bit 8).
  2. Build a peak-height histogram and a peak-integral histogram.
  3. Run a peak search + Gaussian fit on each → mean / sigma.
  4. Optionally query prad2hvd for the current VSet / VMon for that module.
  5. Append ONE new entry per channel to the unified history JSON file.

Each call appends exactly one iteration regardless of how many input files
are provided.

Pure-Python dependencies (no ROOT install needed):
    uproot, awkward, numpy, scipy, matplotlib
    pip install uproot awkward numpy scipy matplotlib

Typical use (one EVIO file → one ROOT file):
    replay_rawdata prad_023527.evio -o prad_023527.root -p
    python3 cosmic_eq_analyze.py prad_023527.root --history cosmic_history.json

Multiple EVIO splits (one ROOT file per split, all chained for one iteration):
    for f in /data/stage6/prad_023600/prad_023600.evio.000??; do
        replay_rawdata $f -o /tmp/$(basename $f).root -p
    done
    python3 cosmic_eq_analyze.py /tmp/prad_023600.evio.*.root \\
            --history cosmic_history.json

Glob form (shell expansion or --glob):
    python3 cosmic_eq_analyze.py --glob '/tmp/prad_023600.evio.*.root' \\
            --history cosmic_history.json
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import uproot
import awkward as ak
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from cosmic_eq_common import (
    HVClient, filter_modules_by_type, load_daq_map, load_history,
    load_hycal_modules, n_iterations, natural_module_sort_key, save_history,
)

SUM_TRIGGER_BIT  = 8
SUM_TRIGGER_MASK = 1 << SUM_TRIGGER_BIT


# ============================================================================
#  Argument parsing
# ============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="HyCal per-channel peak-height/integral fits → cosmic history JSON")
    p.add_argument("input_root", nargs="*",
                   help="One or more ROOT files from replay_rawdata -p (chained)")
    p.add_argument("--glob", default=None,
                   help="Shell-style glob for input files (alternative to positional)")
    p.add_argument("--history", required=True,
                   help="Path to cosmic history JSON (created/appended)")
    p.add_argument("--pdf", default=None,
                   help="Optional PDF for per-channel fit canvases (default: <history>_iterN.pdf)")
    p.add_argument("-n", "--max-events", type=int, default=-1,
                   help="Process at most N events")
    p.add_argument("--exclusive", action="store_true",
                   help="Require ONLY the sum bit (reject events with any other bit set)")
    p.add_argument("--types", default="PbGlass,PbWO4",
                   help="Module types to include (default: PbGlass,PbWO4; 'all' for everything)")
    # peak-height histogram
    p.add_argument("--ph-min",  type=float, default=0.0,    help="Peak-height hist min ADC")
    p.add_argument("--ph-max",  type=float, default=4096.0, help="Peak-height hist max ADC")
    p.add_argument("--ph-bins", type=int,   default=256,    help="Peak-height hist bins")
    # peak-integral histogram
    p.add_argument("--pi-min",  type=float, default=0.0,     help="Peak-integral hist min")
    p.add_argument("--pi-max",  type=float, default=80000.0, help="Peak-integral hist max")
    p.add_argument("--pi-bins", type=int,   default=400,     help="Peak-integral hist bins")
    # fit knobs
    p.add_argument("--min-entries", type=int, default=200,
                   help="Skip channels with fewer entries (default: 200)")
    p.add_argument("--peak-prom", type=float, default=0.10,
                   help="find_peaks prominence as fraction of histogram max (default: 0.10)")
    p.add_argument("--fit-window", type=float, default=2.5,
                   help="Fit window = peak ± N * rough_sigma (default: 2.5)")
    # HV
    p.add_argument("--hv-url", default="http://clonpc19:8765",
                   help="prad2hvd URL (default: http://clonpc19:8765)")
    p.add_argument("--no-hv", action="store_true",
                   help="Don't query prad2hvd; record VSet/VMon as null")
    # db overrides
    p.add_argument("--daq-map", default=None)
    p.add_argument("--modules-db", default=None)
    p.add_argument("-q", "--quiet", action="store_true")
    return p.parse_args()


def resolve_inputs(args: argparse.Namespace) -> List[str]:
    """Combine positional file list with --glob; expand any shell-style globs."""
    import glob as _glob
    files: List[str] = []
    for entry in args.input_root:
        if any(c in entry for c in "*?["):
            matches = sorted(_glob.glob(entry))
            if not matches:
                print(f"WARNING: positional glob '{entry}' matched nothing")
            files.extend(matches)
        else:
            files.append(entry)
    if args.glob:
        matches = sorted(_glob.glob(args.glob))
        if not matches:
            print(f"WARNING: --glob '{args.glob}' matched nothing")
        files.extend(matches)
    if not files:
        sys.exit("ERROR: no input ROOT files (give one or more positional args, "
                 "or use --glob '...')")
    seen, unique = set(), []
    for f in files:
        if f not in seen:
            seen.add(f); unique.append(f)
    return unique


# ============================================================================
#  Lightweight Histogram (replaces TH1F)
# ============================================================================

class Histogram:
    __slots__ = ("nbins", "lo", "hi", "width", "counts")

    def __init__(self, nbins: int, lo: float, hi: float):
        self.nbins  = int(nbins)
        self.lo     = float(lo)
        self.hi     = float(hi)
        self.width  = (self.hi - self.lo) / self.nbins
        self.counts = np.zeros(self.nbins, dtype=np.int64)

    def fill_many(self, values: np.ndarray) -> None:
        v = np.asarray(values, dtype=np.float64)
        v = v[np.isfinite(v)]
        v = v[(v >= self.lo) & (v < self.hi)]
        if v.size == 0:
            return
        bins = ((v - self.lo) / self.width).astype(np.int64)
        np.add.at(self.counts, bins, 1)

    @property
    def entries(self) -> int:
        return int(self.counts.sum())

    def bin_centers(self) -> np.ndarray:
        return self.lo + (np.arange(self.nbins) + 0.5) * self.width


# ============================================================================
#  Branch-name discovery (Replay.cpp uses 'hycal.<x>' branch names)
# ============================================================================

BRANCH_CANDIDATES = {
    "trigger":      ["trigger"],
    "nch":          ["hycal.nch", "nch"],
    "crate":        ["hycal.crate", "crate"],
    "slot":         ["hycal.slot", "slot"],
    "channel":      ["hycal.channel", "channel"],
    "npeaks":       ["hycal.npeaks", "npeaks"],
    "peak_height":  ["hycal.peak_height", "peak_height"],
    "peak_integral": ["hycal.peak_integral", "peak_integral"],
}


def resolve_branches(tree) -> Dict[str, str]:
    """Map our logical names to whichever branch name actually exists."""
    keys = set(tree.keys())
    out: Dict[str, str] = {}
    for logical, candidates in BRANCH_CANDIDATES.items():
        for c in candidates:
            if c in keys:
                out[logical] = c
                break
        else:
            sys.exit(f"ERROR: tree has none of {candidates} (looking for '{logical}')")
    return out


# ============================================================================
#  Vectorised histogram building (uproot + awkward + numpy)
# ============================================================================

def build_histograms(input_files: List[str],
                     daq_map: Dict[Tuple[int, int, int], str],
                     module_set: set,
                     args: argparse.Namespace
                    ) -> Tuple[Dict[str, Histogram], Dict[str, Histogram]]:
    hists_ph: Dict[str, Histogram] = {}
    hists_pi: Dict[str, Histogram] = {}

    # Pre-build a numpy lookup table: lookup[crate, slot, channel] → channel ID
    # (-1 = not in module set). Replaces a dict lookup in the inner loop.
    MAX_CRATE, MAX_SLOT, MAX_CHAN = 32, 32, 16
    lookup = np.full((MAX_CRATE, MAX_SLOT, MAX_CHAN), -1, dtype=np.int32)
    id_to_name: List[str] = []
    for (c, s, ch), name in daq_map.items():
        if name not in module_set:
            continue
        if 0 <= c < MAX_CRATE and 0 <= s < MAX_SLOT and 0 <= ch < MAX_CHAN:
            lookup[c, s, ch] = len(id_to_name)
            id_to_name.append(name)

    def get_or_make_ph(name: str) -> Histogram:
        h = hists_ph.get(name)
        if h is None:
            h = Histogram(args.ph_bins, args.ph_min, args.ph_max)
            hists_ph[name] = h
        return h

    def get_or_make_pi(name: str) -> Histogram:
        h = hists_pi.get(name)
        if h is None:
            h = Histogram(args.pi_bins, args.pi_min, args.pi_max)
            hists_pi[name] = h
        return h

    n_passed = 0
    n_skipped = 0
    n_total = 0
    t0 = time.time()
    print(f"Reading {len(input_files)} file(s)")

    for path in input_files:
        try:
            f = uproot.open(path)
        except Exception as e:
            print(f"WARNING: cannot open {path}: {e}")
            continue
        if "events" not in f:
            print(f"WARNING: no 'events' tree in {path}")
            continue
        tree = f["events"]
        bnames = resolve_branches(tree)
        n_in_file = tree.num_entries
        if not args.quiet:
            print(f"  {path}  ({n_in_file} events)")

        for chunk in tree.iterate(list(bnames.values()), step_size="100 MB",
                                   library="ak"):
            n_chunk = len(chunk)

            triggers = np.asarray(chunk[bnames["trigger"]])
            ev_mask = (triggers & SUM_TRIGGER_MASK) != 0
            if args.exclusive:
                ev_mask &= (triggers & ~SUM_TRIGGER_MASK) == 0

            # max-events cap
            if args.max_events > 0:
                room = max(0, args.max_events - n_total)
                if room < n_chunk:
                    cap = np.zeros(n_chunk, dtype=bool)
                    cap[:room] = True
                    ev_mask &= cap
                    n_total += room
                else:
                    n_total += n_chunk
            else:
                n_total += n_chunk
            n_passed  += int(np.sum(ev_mask))
            n_skipped += int(np.sum(~ev_mask))

            if not ev_mask.any():
                if args.max_events > 0 and n_total >= args.max_events:
                    break
                continue

            crate = chunk[bnames["crate"]][ev_mask]
            slot  = chunk[bnames["slot"]][ev_mask]
            chnl  = chunk[bnames["channel"]][ev_mask]
            npks  = chunk[bnames["npeaks"]][ev_mask]
            ph    = chunk[bnames["peak_height"]][ev_mask]
            pi    = chunk[bnames["peak_integral"]][ev_mask]

            # Flatten across events. peak_height/integral are jagged in nch
            # and fixed in the inner dim — we want index [0] (highest peak).
            flat_crate = ak.flatten(crate).to_numpy().astype(np.int64)
            flat_slot  = ak.flatten(slot).to_numpy().astype(np.int64)
            flat_chnl  = ak.flatten(chnl).to_numpy().astype(np.int64)
            flat_npks  = ak.flatten(npks).to_numpy().astype(np.int64)
            flat_ph    = ak.flatten(ph[:, :, 0]).to_numpy().astype(np.float64)
            flat_pi    = ak.flatten(pi[:, :, 0]).to_numpy().astype(np.float64)

            # filter: npeaks > 0 and (crate,slot,chan) within bounds
            ok = ((flat_npks > 0)
                  & (flat_crate >= 0) & (flat_crate < MAX_CRATE)
                  & (flat_slot  >= 0) & (flat_slot  < MAX_SLOT)
                  & (flat_chnl  >= 0) & (flat_chnl  < MAX_CHAN))
            flat_crate = flat_crate[ok]
            flat_slot  = flat_slot[ok]
            flat_chnl  = flat_chnl[ok]
            flat_ph    = flat_ph[ok]
            flat_pi    = flat_pi[ok]

            if flat_crate.size == 0:
                continue

            # One vectorised lookup → channel IDs (or -1)
            ids = lookup[flat_crate, flat_slot, flat_chnl]
            valid = ids >= 0
            ids = ids[valid]
            ph_vals = flat_ph[valid]
            pi_vals = flat_pi[valid]
            if ids.size == 0:
                continue

            # Group by channel and bulk-fill each histogram
            unique_ids, inverse = np.unique(ids, return_inverse=True)
            for i, cid in enumerate(unique_ids):
                name = id_to_name[cid]
                cmask = inverse == i
                get_or_make_ph(name).fill_many(ph_vals[cmask])
                get_or_make_pi(name).fill_many(pi_vals[cmask])

            if args.max_events > 0 and n_total >= args.max_events:
                break

        if args.max_events > 0 and n_total >= args.max_events:
            break

    elapsed = time.time() - t0
    print(f"Processed {n_total} events ({n_passed} passed sum filter, "
          f"{n_skipped} rejected) in {elapsed:.1f}s; filled {len(hists_ph)} channels.")
    return hists_ph, hists_pi


# ============================================================================
#  Per-channel peak search + Gaussian fit
# ============================================================================

class FitOutcome:
    __slots__ = ("mean", "sigma", "chi2_ndf", "status",
                 "fit_lo", "fit_hi", "amp")

    def __init__(self, mean: float, sigma: float, chi2_ndf: float,
                 status: str, fit_lo: float = 0.0, fit_hi: float = 0.0,
                 amp: float = 0.0):
        self.mean     = mean
        self.sigma    = sigma
        self.chi2_ndf = chi2_ndf
        self.status   = status
        self.fit_lo   = fit_lo
        self.fit_hi   = fit_hi
        self.amp      = amp


def _gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_one(hist: Histogram, args: argparse.Namespace) -> FitOutcome:
    entries = hist.entries
    if entries < args.min_entries:
        return FitOutcome(0.0, 0.0, 0.0, "LOW_STATS")

    centers = hist.bin_centers()
    y = hist.counts.astype(np.float64)

    ymax = float(y.max())
    if ymax <= 0:
        return FitOutcome(0.0, 0.0, 0.0, "NO_PEAK")
    peaks, _ = find_peaks(y, prominence=args.peak_prom * ymax)
    if peaks.size == 0:
        peak_idx = int(np.argmax(y))
    else:
        peak_idx = int(peaks[np.argmax(y[peaks])])

    peak_x = centers[peak_idx]
    peak_y = y[peak_idx]
    if peak_y <= 0:
        return FitOutcome(peak_x, 0.0, 0.0, "NO_PEAK")

    # rough sigma from FWHM around peak
    half = peak_y / 2.0
    lo_idx = peak_idx
    while lo_idx > 0 and y[lo_idx] > half:
        lo_idx -= 1
    hi_idx = peak_idx
    while hi_idx < hist.nbins - 1 and y[hi_idx] > half:
        hi_idx += 1
    fwhm = centers[hi_idx] - centers[lo_idx]
    rough_sigma = max(fwhm / 2.355, 2 * hist.width)

    fit_lo = max(peak_x - args.fit_window * rough_sigma, hist.lo)
    fit_hi = min(peak_x + args.fit_window * rough_sigma, hist.hi)
    mask = (centers >= fit_lo) & (centers <= fit_hi)
    if mask.sum() < 4:
        return FitOutcome(peak_x, rough_sigma, 0.0, "FIT_FAIL", fit_lo, fit_hi)

    # Poisson errors (sqrt(N), floored at 1 for empty bins)
    sigma_y = np.sqrt(np.maximum(y[mask], 1.0))
    try:
        popt, _ = curve_fit(
            _gauss, centers[mask], y[mask],
            p0=[peak_y, peak_x, rough_sigma],
            sigma=sigma_y, absolute_sigma=False, maxfev=2000)
    except Exception:
        return FitOutcome(peak_x, rough_sigma, 0.0, "FIT_FAIL", fit_lo, fit_hi)

    A, mu, sigma = popt
    sigma = abs(sigma)
    fit_y = _gauss(centers[mask], A, mu, sigma)
    chi2 = float(np.sum(((y[mask] - fit_y) / sigma_y) ** 2))
    ndf  = max(int(mask.sum()) - 3, 1)
    return FitOutcome(float(mu), float(sigma), chi2 / ndf, "OK",
                      fit_lo, fit_hi, float(A))


# ============================================================================
#  PDF output (matplotlib)
# ============================================================================

def write_pdf(out_pdf: str,
              ordered: List[str],
              hists_ph: Dict[str, Histogram],
              hists_pi: Dict[str, Histogram],
              fits_ph: Dict[str, FitOutcome],
              fits_pi: Dict[str, FitOutcome],
              rows_per_page: int = 4) -> None:
    plt.rcParams.update({"font.size": 8})
    with PdfPages(out_pdf) as pdf:
        i = 0
        fig, axes = plt.subplots(rows_per_page, 2, figsize=(10, 11))
        for name in ordered:
            row = i % rows_per_page
            if row == 0 and i > 0:
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                fig, axes = plt.subplots(rows_per_page, 2, figsize=(10, 11))

            _draw_panel(axes[row][0], name, "peak height",
                        hists_ph.get(name), fits_ph.get(name), "steelblue")
            _draw_panel(axes[row][1], name, "peak integral",
                        hists_pi.get(name), fits_pi.get(name), "seagreen")
            i += 1

        if i > 0:
            tail = i % rows_per_page
            if tail != 0:
                for r in range(tail, rows_per_page):
                    axes[r][0].axis("off")
                    axes[r][1].axis("off")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def _draw_panel(ax, name: str, label: str,
                hist: Optional[Histogram],
                fit:  Optional[FitOutcome],
                color: str) -> None:
    if hist is None or hist.entries == 0:
        ax.axis("off")
        return
    centers = hist.bin_centers()
    ax.bar(centers, hist.counts, width=hist.width, color=color,
           edgecolor="none", linewidth=0)
    if fit is not None and fit.status == "OK":
        xs = np.linspace(fit.fit_lo, fit.fit_hi, 200)
        ys = _gauss(xs, fit.amp, fit.mean, fit.sigma)
        ax.plot(xs, ys, "r-", lw=1.5)
        ax.text(0.97, 0.95,
                f"mu={fit.mean:.2f}\nsig={fit.sigma:.2f}\nchi2/ndf={fit.chi2_ndf:.2f}\nN={hist.entries}",
                transform=ax.transAxes, ha="right", va="top",
                family="monospace",
                bbox=dict(boxstyle="round,pad=0.2",
                          facecolor="white", alpha=0.8, edgecolor="0.5"))
    elif fit is not None:
        ax.text(0.97, 0.95, f"{fit.status}\nN={hist.entries}",
                transform=ax.transAxes, ha="right", va="top",
                color="red", family="monospace",
                bbox=dict(boxstyle="round,pad=0.2",
                          facecolor="white", alpha=0.8, edgecolor="0.5"))
    ax.set_yscale("log")
    ax.set_title(f"{name}  {label}", fontsize=9)
    ax.set_xlabel("ADC", fontsize=8)
    ax.tick_params(labelsize=7)


# ============================================================================
#  Main
# ============================================================================

def main() -> int:
    args = parse_args()
    input_files = resolve_inputs(args)
    print(f"Input files ({len(input_files)}):")
    for f in input_files:
        print(f"  {f}")

    daq_map = load_daq_map(args.daq_map) if args.daq_map else load_daq_map()
    modules = (load_hycal_modules(args.modules_db)
               if args.modules_db else load_hycal_modules())

    all_names = sorted({n for n in daq_map.values()})
    if args.types.lower() == "all":
        wanted_names = all_names
    else:
        wanted_types = [t.strip() for t in args.types.split(",") if t.strip()]
        wanted_names = filter_modules_by_type(all_names, modules, wanted_types)
    print(f"Selected {len(wanted_names)} modules of types: {args.types}")
    module_set = set(wanted_names)

    hists_ph, hists_pi = build_histograms(input_files, daq_map, module_set, args)

    ordered = sorted(set(hists_ph.keys()) | set(hists_pi.keys()),
                     key=natural_module_sort_key)
    fits_ph: Dict[str, FitOutcome] = {}
    fits_pi: Dict[str, FitOutcome] = {}
    for name in ordered:
        if name in hists_ph: fits_ph[name] = fit_one(hists_ph[name], args)
        if name in hists_pi: fits_pi[name] = fit_one(hists_pi[name], args)
        if not args.quiet:
            r_ph = fits_ph.get(name); r_pi = fits_pi.get(name)
            line = f"  {name:<8}"
            if r_ph: line += f"  PH mu={r_ph.mean:>8.2f} sig={r_ph.sigma:>6.2f} {r_ph.status}"
            if r_pi: line += f"  PI mu={r_pi.mean:>10.1f} sig={r_pi.sigma:>8.1f} {r_pi.status}"
            print(line)

    # Optional HV query
    hv_data: Dict[str, Optional[dict]] = {}
    if not args.no_hv:
        hv = HVClient(args.hv_url)
        try:
            hv.authenticate()
        except Exception as e:
            print(f"WARNING: HV authentication failed ({e}); proceeding read-only")
        for name in ordered:
            try:
                hv_data[name] = hv.get_voltage(name)
            except Exception as e:
                print(f"WARNING: HV GET {name} failed: {e}")
                hv_data[name] = None

    # Append to history
    channels = load_history(args.history)
    iter_idx = n_iterations(channels) + 1
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    input_basenames = [os.path.basename(p) for p in input_files]
    input_field: Any = (input_basenames[0] if len(input_basenames) == 1
                        else input_basenames)

    for name in ordered:
        r_ph = fits_ph.get(name)
        r_pi = fits_pi.get(name)
        h_ph = hists_ph.get(name)
        count = h_ph.entries if h_ph is not None else 0

        hv_entry = hv_data.get(name)
        vset = hv_entry["vset"] if hv_entry else None
        vmon = hv_entry["vmon"] if hv_entry else None

        if r_ph is None or r_pi is None:
            status = "MISSING"
        elif r_ph.status == "OK" and r_pi.status == "OK":
            status = "OK"
        else:
            status = f"{r_ph.status if r_ph else '-'}/{r_pi.status if r_pi else '-'}"

        entry = {
            "iter":               iter_idx,
            "timestamp":          timestamp,
            "input":              input_field,
            "count":              count,
            "peak_height_mean":   round(r_ph.mean,  3) if r_ph else None,
            "peak_height_sigma":  round(r_ph.sigma, 3) if r_ph else None,
            "peak_height_chi2":   round(r_ph.chi2_ndf, 3) if r_ph else None,
            "peak_integral_mean": round(r_pi.mean,  2) if r_pi else None,
            "peak_integral_sigma": round(r_pi.sigma, 2) if r_pi else None,
            "peak_integral_chi2": round(r_pi.chi2_ndf, 3) if r_pi else None,
            "VMon":               round(vmon, 2) if vmon is not None else None,
            "VSet":               round(vset, 2) if vset is not None else None,
            "status":             status,
        }
        # HV addressing — required by prad2hvd /api/load_settings
        if hv_entry is not None:
            if "crate"   in hv_entry: entry["hv_crate"]   = hv_entry["crate"]
            if "slot"    in hv_entry: entry["hv_slot"]    = hv_entry["slot"]
            if "channel" in hv_entry: entry["hv_channel"] = hv_entry["channel"]
        channels.setdefault(name, []).append(entry)

    save_history(args.history, channels)
    print(f"Appended iteration {iter_idx} to {args.history}")

    pdf_path = args.pdf or f"{os.path.splitext(args.history)[0]}_iter{iter_idx}.pdf"
    write_pdf(pdf_path, ordered, hists_ph, hists_pi, fits_ph, fits_pi)
    print(f"Wrote {pdf_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
