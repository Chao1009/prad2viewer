#!/usr/bin/env python3
"""
cosmic_eq_distribution — pass 2 of the offline cosmic-data gain-equalization
workflow.

Reads the unified cosmic history JSON, picks one iteration (latest by default),
and produces a histogram of per-channel mean values across all good fits.
This is what shows the bimodal distribution: an "already equalized" group at
high mean and an "underequalized" group at low mean.

Outputs:
    <out>.pdf — multi-page PDF with:
                  · histogram of peak_height_mean across channels
                  · histogram of peak_integral_mean
                  · scatter peak_*_mean vs VSet (when HV recorded)
                  · per-iteration overlay (if multiple iters present)
    Stdout    — summary stats: n_channels, mean, median, group split estimate.

Pure-Python dependencies (no ROOT install needed):
    numpy, matplotlib

Typical use:
    python3 cosmic_eq_distribution.py cosmic_history.json
    python3 cosmic_eq_distribution.py cosmic_history.json --iter 2
    python3 cosmic_eq_distribution.py cosmic_history.json --max-sigma 4 --types W
    python3 cosmic_eq_distribution.py cosmic_history.json --metric peak_integral
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, List, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    p.add_argument("history", help="Path to cosmic history JSON")
    p.add_argument("--iter", type=int, default=-1,
                   help="Iteration index (1-based). Default -1 = latest.")
    p.add_argument("--metric", choices=("peak_height", "peak_integral"),
                   default="peak_height",
                   help="Which fitted mean drives the iteration overlay and stats "
                        "(default: peak_height). The 2x2 overview always shows both.")
    p.add_argument("--types", default="PbGlass,PbWO4",
                   help="Module types to include (default: PbGlass,PbWO4; 'all' for everything)")
    p.add_argument("--max-sigma", type=float, default=0.0,
                   help="Quality cut: drop channels with <metric>_sigma above this. 0 = no cut.")
    p.add_argument("--max-chi2", type=float, default=0.0,
                   help="Quality cut: drop channels with <metric>_chi2 above this. 0 = no cut.")
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

def metric_keys(metric: str):
    if metric == "peak_integral":
        return "peak_integral_mean", "peak_integral_sigma", "peak_integral_chi2"
    return "peak_height_mean", "peak_height_sigma", "peak_height_chi2"


def extract_iteration(channels: Dict[str, List[dict]],
                      iter_idx: int,
                      wanted: set,
                      max_sigma: float,
                      max_chi2: float,
                      metric: str = "peak_height") -> List[dict]:
    """Pull one iteration entry per channel, applying quality cuts on the
    chosen metric. iter_idx <= 0 → latest."""
    mean_key, sigma_key, chi2_key = metric_keys(metric)
    out: List[dict] = []
    for name, entries in channels.items():
        if name not in wanted:
            continue
        if not entries:
            continue
        if iter_idx <= 0:
            entry = entries[-1]
        else:
            match = [e for e in entries if e.get("iter") == iter_idx]
            if not match:
                continue
            entry = match[-1]

        if entry.get(mean_key) is None:
            continue
        sig = entry.get(sigma_key) or 0.0
        chi = entry.get(chi2_key) or 0.0
        if max_sigma > 0 and sig > max_sigma:
            continue
        if max_chi2 > 0 and chi > max_chi2:
            continue

        record = dict(entry)
        record["name"] = name
        record["_metric_mean"] = entry[mean_key]
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


def make_overview(records: List[dict], iter_idx: int,
                  args: argparse.Namespace, pdf: PdfPages) -> None:
    ph_lo, ph_hi = parse_range(args.ph_range, (0, 4096))
    pi_lo, pi_hi = parse_range(args.pi_range, (0, 80000))

    ph_means = np.array([r["peak_height_mean"]
                          for r in records
                          if r.get("peak_height_mean") is not None], dtype=float)
    pi_means = np.array([r["peak_integral_mean"]
                          for r in records
                          if r.get("peak_integral_mean") is not None], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Cosmic-eq overview — iteration {iter_idx}", fontsize=12)

    # Pad 1: peak-height mean histogram
    ax = axes[0][0]
    ax.hist(ph_means, bins=args.ph_bins, range=(ph_lo, ph_hi),
            color="steelblue", alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("peak_height_mean [ADC]")
    ax.set_ylabel("channels")
    ax.set_title("Peak-height mean per channel")
    ax.grid(True, alpha=0.3)

    # Pad 2: peak-integral mean histogram
    ax = axes[0][1]
    ax.hist(pi_means, bins=args.pi_bins, range=(pi_lo, pi_hi),
            color="seagreen", alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("peak_integral_mean [ADC]")
    ax.set_ylabel("channels")
    ax.set_title("Peak-integral mean per channel")
    ax.grid(True, alpha=0.3)

    # Pad 3: scatter mean vs VSet (using the chosen metric)
    ax = axes[1][0]
    metric_label = ("peak_height_mean" if args.metric == "peak_height"
                    else "peak_integral_mean")
    vsets, ms = [], []
    for r in records:
        v = r.get("VSet")
        if v is not None and r.get("_metric_mean") is not None:
            vsets.append(v); ms.append(r["_metric_mean"])
    if vsets:
        ax.scatter(vsets, ms, s=12, color="purple", alpha=0.7)
        ax.set_xlabel("VSet [V]")
        ax.set_ylabel(f"{metric_label} [ADC]")
        ax.set_title(f"{metric_label} vs VSet")
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, "no VSet recorded for this iteration",
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")

    # Pad 4: stats summary
    ax = axes[1][1]
    ax.axis("off")
    summ = compute_stats(records)
    text = "\n".join(summ)
    ax.text(0.05, 0.95, text, transform=ax.transAxes,
            ha="left", va="top", family="monospace", fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def make_iteration_overlay(channels: Dict[str, List[dict]],
                           wanted: set,
                           args: argparse.Namespace,
                           pdf: PdfPages) -> None:
    n_iters = max((len(v) for v in channels.values()), default=0)
    if n_iters < 2:
        return

    ph_lo, ph_hi = parse_range(args.ph_range, (0, 4096))
    metric_label = ("peak_height_mean" if args.metric == "peak_height"
                    else "peak_integral_mean")

    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.tab10.colors

    drew_any = False
    for it in range(1, n_iters + 1):
        recs = extract_iteration(channels, it, wanted, args.max_sigma,
                                  args.max_chi2, args.metric)
        if not recs:
            continue
        means = np.array([r["_metric_mean"] for r in recs], dtype=float)
        if means.size == 0:
            continue
        ax.hist(means, bins=args.ph_bins, range=(ph_lo, ph_hi),
                histtype="step", linewidth=2,
                color=colors[(it - 1) % len(colors)],
                label=f"iter {it}  (n={means.size})")
        drew_any = True

    if not drew_any:
        plt.close(fig)
        return

    ax.set_xlabel(f"{metric_label} [ADC]")
    ax.set_ylabel("channels")
    ax.set_title(f"Iteration overlay — {metric_label}")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


# ============================================================================
#  Stats
# ============================================================================

def compute_stats(records: List[dict]) -> List[str]:
    n = len(records)
    if n == 0:
        return ["no records"]
    means = sorted(r["_metric_mean"] for r in records)
    median = means[n // 2]
    avg = sum(means) / n
    mn, mx = means[0], means[-1]

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
                                 args.max_sigma, args.max_chi2, args.metric)
    print(f"After cuts: {len(records)} channels")

    out_pdf = (args.output or
               f"{os.path.splitext(args.history)[0]}_dist_iter{iter_idx}.pdf")

    with PdfPages(out_pdf) as pdf:
        make_overview(records, iter_idx, args, pdf)
        make_iteration_overlay(channels, wanted, args, pdf)

    print(f"Wrote {out_pdf}")
    print()
    for line in compute_stats(records):
        print(f"  {line}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
