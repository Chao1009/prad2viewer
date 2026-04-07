#!/usr/bin/env python3
"""
gain_eq_analyze — pass 1 of the offline gain-equalization workflow.

Reads a ROOT file produced by ``replay_rawdata -p`` (which contains
peak_height[nch][MAX_PEAKS] and peak_integral[nch][MAX_PEAKS] plus the FP
trigger word per event). For each HyCal channel:

  1. Filter events: only those with the sum trigger bit set (bit 8).
  2. Build a peak-height histogram and a peak-integral histogram.
  3. Run TSpectrum + Gaussian fit on both → mean / sigma.
  4. Optionally query prad2hvd for the current VSet / VMon for that module.
  5. Append a new entry to the channel's iteration list in the unified
     history JSON file.

Each call appends one new iteration. The same JSON file is fed to
gain_eq_distribution.py and gain_eq_propose.py.

Typical use:
    # iter 1
    replay_rawdata prad_023527.evio -o prad_023527.root -p
    python3 gain_eq_analyze.py prad_023527.root --history gain_history.json

    # iter 2 (after applying new VSet and taking new data)
    replay_rawdata prad_023600.evio -o prad_023600.root -p
    python3 gain_eq_analyze.py prad_023600.root --history gain_history.json
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import ROOT  # noqa: E402

from gain_eq_common import (
    HVClient, filter_modules_by_type, latest_iteration, load_daq_map,
    load_history, load_hycal_modules, n_iterations, natural_module_sort_key,
    save_history,
)

SUM_TRIGGER_BIT  = 8
SUM_TRIGGER_MASK = 1 << SUM_TRIGGER_BIT


# ============================================================================
#  Argument parsing
# ============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="HyCal per-channel peak-height/integral fits → gain history JSON")
    p.add_argument("input_root",
                   help="ROOT file from replay_rawdata -p")
    p.add_argument("--history", required=True,
                   help="Path to gain history JSON (created/appended)")
    p.add_argument("--pdf", default=None,
                   help="Optional PDF for per-channel fit canvases (default: <history>_iterN.pdf)")
    p.add_argument("--save-hists", default=None,
                   help="Optional ROOT file to save all per-channel histograms")
    p.add_argument("-n", "--max-events", type=int, default=-1,
                   help="Process at most N events")
    p.add_argument("--exclusive", action="store_true",
                   help="Require ONLY the sum bit (reject events with any other bit set)")
    p.add_argument("--types", default="PbGlass,PWO",
                   help="Comma-separated module types to include (default: PbGlass,PWO). "
                        "Use 'all' to include everything.")
    # peak height histogram
    p.add_argument("--ph-min",  type=float, default=0.0,    help="Peak-height hist min ADC")
    p.add_argument("--ph-max",  type=float, default=4096.0, help="Peak-height hist max ADC")
    p.add_argument("--ph-bins", type=int,   default=256,    help="Peak-height hist bins")
    # peak integral histogram
    p.add_argument("--pi-min",  type=float, default=0.0,     help="Peak-integral hist min")
    p.add_argument("--pi-max",  type=float, default=80000.0, help="Peak-integral hist max")
    p.add_argument("--pi-bins", type=int,   default=400,     help="Peak-integral hist bins")
    # fit knobs
    p.add_argument("--min-entries", type=int, default=200,
                   help="Skip channels with fewer entries (default: 200)")
    p.add_argument("--peak-sigma", type=float, default=2.0,
                   help="TSpectrum sigma in bins (default: 2)")
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


# ============================================================================
#  Histogram building
# ============================================================================

def build_histograms(input_root: str,
                     daq_map: Dict[Tuple[int, int, int], str],
                     module_set: set,
                     args: argparse.Namespace
                    ) -> Tuple[Dict[str, "ROOT.TH1F"], Dict[str, "ROOT.TH1F"]]:
    """Loop over events; fill peak-height + peak-integral hists per channel."""
    f = ROOT.TFile.Open(input_root)
    if not f or f.IsZombie():
        sys.exit(f"ERROR: cannot open {input_root}")
    tree = f.Get("events")
    if not tree:
        sys.exit(f"ERROR: 'events' tree not found in {input_root}")

    n_total = tree.GetEntries()
    n_to_process = n_total if args.max_events < 0 else min(n_total, args.max_events)
    print(f"Reading {n_to_process}/{n_total} events from {input_root}")

    hists_ph: Dict[str, ROOT.TH1F] = {}
    hists_pi: Dict[str, ROOT.TH1F] = {}

    def get_or_make_ph(name: str) -> ROOT.TH1F:
        h = hists_ph.get(name)
        if h is None:
            h = ROOT.TH1F(f"hph_{name}",
                          f"{name} peak height;peak height [ADC];events",
                          args.ph_bins, args.ph_min, args.ph_max)
            h.SetDirectory(0)
            hists_ph[name] = h
        return h

    def get_or_make_pi(name: str) -> ROOT.TH1F:
        h = hists_pi.get(name)
        if h is None:
            h = ROOT.TH1F(f"hpi_{name}",
                          f"{name} peak integral;peak integral [ADC];events",
                          args.pi_bins, args.pi_min, args.pi_max)
            h.SetDirectory(0)
            hists_pi[name] = h
        return h

    n_passed = 0
    n_skipped_trig = 0
    t0 = time.time()
    for i in range(n_to_process):
        tree.GetEntry(i)
        trig = int(tree.trigger)

        if not (trig & SUM_TRIGGER_MASK):
            n_skipped_trig += 1
            continue
        if args.exclusive and (trig & ~SUM_TRIGGER_MASK):
            n_skipped_trig += 1
            continue

        nch = int(tree.nch)
        for k in range(nch):
            crate = int(tree.crate[k])
            slot  = int(tree.slot[k])
            ch    = int(tree.channel[k])
            name  = daq_map.get((crate, slot, ch))
            if not name or name not in module_set:
                continue
            npeaks = int(tree.npeaks[k])
            if npeaks <= 0:
                continue
            ph = float(tree.peak_height[k][0])
            pi = float(tree.peak_integral[k][0])
            get_or_make_ph(name).Fill(ph)
            get_or_make_pi(name).Fill(pi)

        n_passed += 1
        if not args.quiet and (i + 1) % 50000 == 0:
            rate = (i + 1) / max(time.time() - t0, 1e-3)
            print(f"  ... {i + 1}/{n_to_process}  ({rate:.0f} ev/s)")

    print(f"Done. {n_passed} events passed sum-trigger filter, "
          f"{n_skipped_trig} rejected. Filled {len(hists_ph)} channels.")
    f.Close()
    return hists_ph, hists_pi


# ============================================================================
#  Per-channel TSpectrum + Gaussian fit
# ============================================================================

class FitOutcome:
    __slots__ = ("mean", "sigma", "chi2_ndf", "status", "fit_lo", "fit_hi")

    def __init__(self, mean: float, sigma: float, chi2_ndf: float,
                 status: str, fit_lo: float = 0.0, fit_hi: float = 0.0):
        self.mean     = mean
        self.sigma    = sigma
        self.chi2_ndf = chi2_ndf
        self.status   = status
        self.fit_lo   = fit_lo
        self.fit_hi   = fit_hi


def fit_one(hist: "ROOT.TH1F", args: argparse.Namespace) -> FitOutcome:
    entries = int(hist.GetEntries())
    if entries < args.min_entries:
        return FitOutcome(0.0, 0.0, 0.0, "LOW_STATS")

    spec = ROOT.TSpectrum(10)
    nfound = spec.Search(hist, args.peak_sigma, "nodraw nobackground", 0.10)
    if nfound <= 0:
        peak_x = hist.GetBinCenter(hist.GetMaximumBin())
    else:
        xs = spec.GetPositionX()
        best_x = xs[0]
        best_h = hist.GetBinContent(hist.FindBin(xs[0]))
        for j in range(1, nfound):
            h_j = hist.GetBinContent(hist.FindBin(xs[j]))
            if h_j > best_h:
                best_h = h_j
                best_x = xs[j]
        peak_x = best_x

    # rough sigma from FWHM around peak
    peak_bin = hist.FindBin(peak_x)
    peak_y = hist.GetBinContent(peak_bin)
    if peak_y <= 0:
        return FitOutcome(peak_x, 0.0, 0.0, "NO_PEAK")
    half = peak_y / 2.0
    lo = peak_bin
    while lo > 1 and hist.GetBinContent(lo) > half:
        lo -= 1
    hi = peak_bin
    while hi < hist.GetNbinsX() and hist.GetBinContent(hi) > half:
        hi += 1
    fwhm = hist.GetBinCenter(hi) - hist.GetBinCenter(lo)
    rough_sigma = max(fwhm / 2.355, 2 * hist.GetBinWidth(1))

    fit_lo = max(peak_x - args.fit_window * rough_sigma, hist.GetXaxis().GetXmin())
    fit_hi = min(peak_x + args.fit_window * rough_sigma, hist.GetXaxis().GetXmax())

    fn = ROOT.TF1(f"g_{hist.GetName()}", "gaus", fit_lo, fit_hi)
    fn.SetParameters(peak_y, peak_x, rough_sigma)
    res = hist.Fit(fn, "RQNS")

    if not res.Get() or res.Status() != 0:
        return FitOutcome(peak_x, rough_sigma, 0.0, "FIT_FAIL", fit_lo, fit_hi)

    mean  = fn.GetParameter(1)
    sigma = abs(fn.GetParameter(2))
    chi2  = fn.GetChisquare()
    ndf   = max(fn.GetNDF(), 1)

    # attach a fit function to the histogram for the PDF report
    fn_save = ROOT.TF1(f"fit_{hist.GetName()}", "gaus", fit_lo, fit_hi)
    fn_save.SetParameters(fn.GetParameter(0), mean, sigma)
    fn_save.SetLineColor(ROOT.kRed)
    fn_save.SetLineWidth(2)
    hist.GetListOfFunctions().Add(fn_save)

    return FitOutcome(mean, sigma, chi2 / ndf, "OK", fit_lo, fit_hi)


# ============================================================================
#  PDF output (one row per channel = peak-height + peak-integral side by side)
# ============================================================================

def write_pdf(out_pdf: str,
              ordered: List[str],
              hists_ph: Dict[str, "ROOT.TH1F"],
              hists_pi: Dict[str, "ROOT.TH1F"],
              rows_per_page: int = 4) -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(1110)
    ROOT.gStyle.SetOptFit(111)

    canvas = ROOT.TCanvas("c_gain", "gain fits", 1200, 900)
    canvas.Divide(2, rows_per_page)
    canvas.Print(f"{out_pdf}[")

    pad_idx = 0
    for name in ordered:
        h_ph = hists_ph.get(name)
        h_pi = hists_pi.get(name)
        if h_ph is None and h_pi is None:
            continue
        # peak height
        pad_idx += 1
        canvas.cd(pad_idx)
        ROOT.gPad.SetLogy(True)
        if h_ph: h_ph.Draw()
        # peak integral
        pad_idx += 1
        canvas.cd(pad_idx)
        ROOT.gPad.SetLogy(True)
        if h_pi: h_pi.Draw()

        if pad_idx >= 2 * rows_per_page:
            canvas.Print(out_pdf)
            canvas.Clear()
            canvas.Divide(2, rows_per_page)
            pad_idx = 0

    if pad_idx > 0:
        canvas.Print(out_pdf)
    canvas.Print(f"{out_pdf}]")


# ============================================================================
#  Main
# ============================================================================

def main() -> int:
    args = parse_args()

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

    # Build histograms
    hists_ph, hists_pi = build_histograms(args.input_root, daq_map, module_set, args)

    # Fit
    ordered = sorted(set(hists_ph.keys()) | set(hists_pi.keys()),
                     key=natural_module_sort_key)
    fits_ph: Dict[str, FitOutcome] = {}
    fits_pi: Dict[str, FitOutcome] = {}
    for name in ordered:
        if name in hists_ph: fits_ph[name] = fit_one(hists_ph[name], args)
        if name in hists_pi: fits_pi[name] = fit_one(hists_pi[name], args)
        if not args.quiet:
            r_ph = fits_ph.get(name)
            r_pi = fits_pi.get(name)
            line = f"  {name:<8}"
            if r_ph:
                line += (f"  PH μ={r_ph.mean:>8.2f} σ={r_ph.sigma:>6.2f}"
                         f" {r_ph.status}")
            if r_pi:
                line += (f"  PI μ={r_pi.mean:>10.1f} σ={r_pi.sigma:>8.1f}"
                         f" {r_pi.status}")
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

    # Append to history JSON
    channels = load_history(args.history)
    iter_idx = n_iterations(channels) + 1
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    input_basename = os.path.basename(args.input_root)

    for name in ordered:
        r_ph = fits_ph.get(name)
        r_pi = fits_pi.get(name)
        h_ph = hists_ph.get(name)
        if h_ph is None:
            count = 0
        else:
            count = int(h_ph.GetEntries())

        hv_entry = hv_data.get(name)
        vset = hv_entry["vset"] if hv_entry else None
        vmon = hv_entry["vmon"] if hv_entry else None

        # combined status — OK only if both fits succeeded
        if r_ph is None or r_pi is None:
            status = "MISSING"
        elif r_ph.status == "OK" and r_pi.status == "OK":
            status = "OK"
        else:
            status = f"{r_ph.status if r_ph else '-'}/{r_pi.status if r_pi else '-'}"

        entry = {
            "iter":               iter_idx,
            "timestamp":          timestamp,
            "input":              input_basename,
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
        channels.setdefault(name, []).append(entry)

    save_history(args.history, channels)
    print(f"Appended iteration {iter_idx} to {args.history}")

    # PDF
    pdf_path = args.pdf or f"{os.path.splitext(args.history)[0]}_iter{iter_idx}.pdf"
    write_pdf(pdf_path, ordered, hists_ph, hists_pi)
    print(f"Wrote {pdf_path}")

    # Optional histograms ROOT file
    if args.save_hists:
        rf = ROOT.TFile.Open(args.save_hists, "RECREATE")
        for name in ordered:
            if name in hists_ph: hists_ph[name].Write()
            if name in hists_pi: hists_pi[name].Write()
        rf.Close()
        print(f"Wrote {args.save_hists}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
