#!/usr/bin/env python3
"""
cosmic_eq_distribution — pass 2 of the offline gain-equalization workflow.

Reads the unified gain history JSON, picks one iteration (latest by default),
and produces a histogram of per-channel mean values across all good fits.
This is what shows the bimodal distribution: an "already equalized" group at
high mean and an "underequalized" group at low mean.

Outputs:
    <out>.pdf — overview canvas with:
                  · histogram of peak_height_mean across channels
                  · histogram of peak_integral_mean
                  · scatter peak_height_mean vs VSet (when HV recorded)
                  · per-iteration overlay (if multiple iters present)
    Stdout    — summary stats: n_channels, mean, median, group split estimate.

Typical use:
    python3 cosmic_eq_distribution.py gain_history.json
    python3 cosmic_eq_distribution.py gain_history.json --iter 2
    python3 cosmic_eq_distribution.py gain_history.json --max-sigma 200 --type W
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, List, Optional

import ROOT  # noqa: E402

from cosmic_eq_common import (
    filter_modules_by_type, load_history, load_hycal_modules,
    natural_module_sort_key,
)


# ============================================================================
#  Argument parsing
# ============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot the distribution of fitted means across channels")
    p.add_argument("history", help="Path to gain history JSON")
    p.add_argument("--iter", type=int, default=-1,
                   help="Iteration index (1-based). Default -1 = latest.")
    p.add_argument("--types", default="PbGlass,PWO",
                   help="Module types to include (default: PbGlass,PWO; 'all' for everything)")
    p.add_argument("--max-sigma", type=float, default=0.0,
                   help="Quality cut: drop channels with peak_height_sigma above this. "
                        "0 = no cut.")
    p.add_argument("--max-chi2", type=float, default=0.0,
                   help="Quality cut: drop channels with peak_height_chi2 above this. "
                        "0 = no cut.")
    p.add_argument("-o", "--output", default=None,
                   help="Output PDF path (default: <history>_dist_iterN.pdf)")
    p.add_argument("--ph-bins", type=int, default=80, help="Hist bins for PH means")
    p.add_argument("--ph-range", default="0,4096", help="PH mean range 'min,max'")
    p.add_argument("--pi-bins", type=int, default=80, help="Hist bins for PI means")
    p.add_argument("--pi-range", default="0,80000", help="PI mean range 'min,max'")
    p.add_argument("--modules-db", default=None)
    return p.parse_args()


# ============================================================================
#  Iteration extraction
# ============================================================================

def extract_iteration(channels: Dict[str, List[dict]],
                      iter_idx: int,
                      wanted: set,
                      max_sigma: float,
                      max_chi2: float) -> List[dict]:
    """Pull one iteration entry per channel, applying quality cuts.

    iter_idx <= 0 → latest.
    """
    out: List[dict] = []
    for name, entries in channels.items():
        if name not in wanted:
            continue
        if not entries:
            continue
        if iter_idx <= 0:
            entry = entries[-1]
        else:
            # match by 1-based "iter" field, not by index — entries can be out of order
            match = [e for e in entries if e.get("iter") == iter_idx]
            if not match:
                continue
            entry = match[-1]

        if entry.get("peak_height_mean") is None:
            continue
        sig = entry.get("peak_height_sigma") or 0.0
        chi = entry.get("peak_height_chi2") or 0.0
        if max_sigma > 0 and sig > max_sigma:
            continue
        if max_chi2 > 0 and chi > max_chi2:
            continue

        record = dict(entry)
        record["name"] = name
        out.append(record)
    out.sort(key=lambda r: natural_module_sort_key(r["name"]))
    return out


# ============================================================================
#  Plotting
# ============================================================================

def parse_range(s: str, default: tuple) -> tuple:
    try:
        a, b = s.split(",")
        return float(a), float(b)
    except Exception:
        return default


def make_overview(records: List[dict],
                  iter_idx: int,
                  args: argparse.Namespace,
                  out_pdf: str) -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(1110)

    ph_lo, ph_hi = parse_range(args.ph_range, (0, 4096))
    pi_lo, pi_hi = parse_range(args.pi_range, (0, 80000))

    h_ph = ROOT.TH1F("h_ph_mean",
                      f"Peak-height mean per channel (iter {iter_idx});"
                      f"peak_height_mean [ADC];channels",
                      args.ph_bins, ph_lo, ph_hi)
    h_pi = ROOT.TH1F("h_pi_mean",
                      f"Peak-integral mean per channel (iter {iter_idx});"
                      f"peak_integral_mean [ADC];channels",
                      args.pi_bins, pi_lo, pi_hi)

    g_ph_vs_v = ROOT.TGraph()
    g_ph_vs_v.SetTitle(f"Peak-height mean vs VSet (iter {iter_idx});"
                       f"VSet [V];peak_height_mean [ADC]")
    g_ph_vs_v.SetMarkerStyle(20)
    g_ph_vs_v.SetMarkerSize(0.6)

    n_with_hv = 0
    for r in records:
        h_ph.Fill(r["peak_height_mean"])
        if r.get("peak_integral_mean") is not None:
            h_pi.Fill(r["peak_integral_mean"])
        v = r.get("VSet")
        if v is not None:
            g_ph_vs_v.SetPoint(n_with_hv, v, r["peak_height_mean"])
            n_with_hv += 1

    canvas = ROOT.TCanvas("c_dist", "gain distribution", 1400, 900)
    canvas.Divide(2, 2)

    canvas.cd(1)
    h_ph.SetFillColorAlpha(ROOT.kBlue - 7, 0.6)
    h_ph.Draw()

    canvas.cd(2)
    h_pi.SetFillColorAlpha(ROOT.kGreen - 7, 0.6)
    h_pi.Draw()

    canvas.cd(3)
    if n_with_hv > 0:
        g_ph_vs_v.Draw("AP")
    else:
        t = ROOT.TLatex(0.5, 0.5, "no VSet recorded for this iteration")
        t.SetNDC(True); t.SetTextAlign(22); t.Draw()

    # Summary text in pad 4
    canvas.cd(4)
    summ = compute_stats(records)
    t = ROOT.TLatex()
    t.SetNDC(True)
    t.SetTextSize(0.045)
    t.SetTextFont(82)
    y = 0.88
    for line in summ:
        t.DrawLatex(0.05, y, line)
        y -= 0.06

    canvas.Update()
    canvas.Print(out_pdf)


# ============================================================================
#  Per-iteration overlay (when ≥ 2 iterations exist)
# ============================================================================

def make_iteration_overlay(channels: Dict[str, List[dict]],
                           wanted: set,
                           args: argparse.Namespace,
                           out_pdf: str) -> None:
    n_iters = max((len(v) for v in channels.values()), default=0)
    if n_iters < 2:
        return

    ph_lo, ph_hi = parse_range(args.ph_range, (0, 4096))
    canvas = ROOT.TCanvas("c_overlay", "iteration overlay", 1200, 800)
    canvas.SetLogy(False)

    leg = ROOT.TLegend(0.7, 0.7, 0.95, 0.92)
    hists = []
    colors = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2,
              ROOT.kMagenta + 1, ROOT.kOrange + 1, ROOT.kCyan + 2]

    for it in range(1, n_iters + 1):
        recs = extract_iteration(channels, it, wanted, args.max_sigma, args.max_chi2)
        if not recs:
            continue
        h = ROOT.TH1F(f"h_overlay_{it}",
                      "Peak-height mean per iteration;peak_height_mean [ADC];channels",
                      args.ph_bins, ph_lo, ph_hi)
        h.SetDirectory(0)
        for r in recs:
            h.Fill(r["peak_height_mean"])
        col = colors[(it - 1) % len(colors)]
        h.SetLineColor(col); h.SetLineWidth(2)
        hists.append(h)
        leg.AddEntry(h, f"iter {it}  (n={int(h.GetEntries())})", "l")

    if not hists:
        return
    hists[0].Draw("HIST")
    for h in hists[1:]:
        h.Draw("HIST SAME")
    leg.Draw()
    canvas.Update()
    canvas.Print(out_pdf)


# ============================================================================
#  Stats
# ============================================================================

def compute_stats(records: List[dict]) -> List[str]:
    n = len(records)
    if n == 0:
        return ["no records"]
    means = sorted(r["peak_height_mean"] for r in records)
    median = means[n // 2]
    avg = sum(means) / n
    mn, mx = means[0], means[-1]

    # crude bimodal split: bin 1/3 from min toward max
    split = mn + (mx - mn) / 3.0
    n_low  = sum(1 for m in means if m <  split)
    n_high = sum(1 for m in means if m >= split)
    avg_low  = (sum(m for m in means if m <  split) / n_low ) if n_low  else 0.0
    avg_high = (sum(m for m in means if m >= split) / n_high) if n_high else 0.0

    return [
        f"channels:  {n}",
        f"mean:      {avg:.1f}",
        f"median:    {median:.1f}",
        f"min..max:  {mn:.1f} .. {mx:.1f}",
        "",
        f"split @    {split:.1f}",
        f"low group:  {n_low:>3}  avg {avg_low:.1f}",
        f"high group: {n_high:>3}  avg {avg_high:.1f}",
        f"gap:        {avg_high - avg_low:.1f}",
    ]


# ============================================================================
#  Main
# ============================================================================

def main() -> int:
    args = parse_args()

    channels = load_history(args.history)
    if not channels:
        sys.exit(f"ERROR: no channels in {args.history}")

    modules = (load_hycal_modules(args.modules_db)
               if args.modules_db else load_hycal_modules())
    all_in_history = list(channels.keys())
    if args.types.lower() == "all":
        wanted = set(all_in_history)
    else:
        wanted_types = [t.strip() for t in args.types.split(",") if t.strip()]
        wanted = set(filter_modules_by_type(all_in_history, modules, wanted_types))

    iter_idx = args.iter
    if iter_idx <= 0:
        iter_idx = max((e.get("iter", 0) for v in channels.values() for e in v),
                       default=0)
    if iter_idx <= 0:
        sys.exit("ERROR: no iterations in history")
    print(f"Using iteration {iter_idx}")

    records = extract_iteration(channels, iter_idx, wanted,
                                 args.max_sigma, args.max_chi2)
    print(f"After cuts: {len(records)} channels")

    out_pdf = (args.output or
               f"{os.path.splitext(args.history)[0]}_dist_iter{iter_idx}.pdf")
    if os.path.exists(out_pdf):
        os.remove(out_pdf)

    # write the multi-page output: overview, then per-iteration overlay
    ROOT.gROOT.SetBatch(True)
    canvas_seed = ROOT.TCanvas("seed", "seed", 100, 100)
    canvas_seed.Print(f"{out_pdf}[")  # open
    canvas_seed.Close()

    make_overview(records, iter_idx, args, out_pdf)
    make_iteration_overlay(channels, wanted, args, out_pdf)

    closer = ROOT.TCanvas("close", "close", 100, 100)
    closer.Print(f"{out_pdf}]")  # close
    closer.Close()

    print(f"Wrote {out_pdf}")
    print()
    for line in compute_stats(records):
        print(f"  {line}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
