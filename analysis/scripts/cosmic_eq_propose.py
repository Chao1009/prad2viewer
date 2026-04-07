#!/usr/bin/env python3
"""
cosmic_eq_propose — pass 3 of the offline gain-equalization workflow.

Reads the unified gain history JSON (created by cosmic_eq_analyze) and writes
a batch VSet JSON in the prad2hvmon_settings_v1 format that can be POSTed
to prad2hvd /api/load_settings.

Iteration logic per channel:
  - 1 iteration in history: apply a constant ΔV (default +100 V) to bring
    the response into the linear regime where we can measure it.
  - ≥ 2 iterations: compute the per-channel response constant
        rc = (mean₂ - mean₁) / (VSet₂ - VSet₁)         [ADC / V]
    using the latest two iterations, then propose
        ΔV = (target_mean - latest_mean) / rc
    clipped to ±max_dv.

Channels above the target by more than `tolerance` are left alone.
Voltage limits from voltage_limits.json are enforced — proposed VSet is
clamped, and channels that hit the limit are reported.

Outputs:
    <out>.json — prad2hvmon_settings_v1 batch settings (the artifact you
                 actually feed to prad2hvd, manually or via --apply)
    Stdout     — table of {name, current_vset, current_mean, rc, dV, new_vset}
                 with reasons for any skipped channels.

Typical use:
    # iteration 1: blanket +100 V to underequalized channels
    python3 cosmic_eq_propose.py gain_history.json \\
            --target-mean 2500 --output vset_iter1.json

    # iteration 2 and beyond: history now has ≥ 2 points → linear correction
    python3 cosmic_eq_propose.py gain_history.json \\
            --target-mean 2500 --output vset_iter2.json

    # Apply directly via prad2hvd HTTP API (asks for password)
    python3 cosmic_eq_propose.py gain_history.json \\
            --target-mean 2500 --output vset_iter1.json --apply
"""

from __future__ import annotations

import argparse
import getpass
import json
import os
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from cosmic_eq_common import (
    HVClient, build_batch_settings, filter_modules_by_type,
    load_history, load_hycal_modules, load_voltage_limits,
    natural_module_sort_key, voltage_limit_for,
)


# ============================================================================
#  Argument parsing
# ============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Propose new HyCal VSet values from gain history → batch JSON")
    p.add_argument("history", help="Path to gain history JSON")
    p.add_argument("-o", "--output", required=True,
                   help="Output batch settings JSON path")
    p.add_argument("--target-mean", type=float, required=True,
                   help="Target value for the chosen --metric (ADC). "
                        "This is the peak_height_mean OR peak_integral_mean "
                        "you want every channel driven toward.")
    p.add_argument("--metric", choices=("peak_height", "peak_integral"),
                   default="peak_height",
                   help="Which fitted mean to drive: 'peak_height' (default) "
                        "uses peak_height_mean/sigma; 'peak_integral' uses "
                        "peak_integral_mean/sigma instead.")
    p.add_argument("--types", default="PbGlass,PbWO4",
                   help="Module types to include (default: PbGlass,PbWO4)")
    p.add_argument("--initial-step", type=float, default=100.0,
                   help="Constant ΔV applied on the first iteration (default: +100)")
    p.add_argument("--initial-threshold", type=float, default=0.0,
                   help="First-iteration only: apply the initial step ONLY when "
                        "|target_mean - current_mean| > threshold. Bidirectional — "
                        "channels significantly below target step UP, channels "
                        "significantly above target step DOWN, channels within "
                        "± threshold of target are left alone. "
                        "Default 0 = step in either direction (current behavior).")
    p.add_argument("--max-dv", type=float, default=80.0,
                   help="Cap on per-channel ΔV in subsequent iterations (default: 80)")
    p.add_argument("--tolerance", type=float, default=50.0,
                   help="Symmetric leave-alone window around target (default: 50). "
                        "Channels with mean in [target - tolerance, target + tolerance] "
                        "are not touched. Override per-side with --tolerance-low / --tolerance-high.")
    p.add_argument("--tolerance-low", type=float, default=None,
                   help="Asymmetric: how far BELOW target is acceptable (defaults to --tolerance). "
                        "Channels with target - mean in [0, tolerance-low] are left alone.")
    p.add_argument("--tolerance-high", type=float, default=None,
                   help="Asymmetric: how far ABOVE target is acceptable (defaults to --tolerance). "
                        "Channels with mean - target in [0, tolerance-high] are left alone.")
    p.add_argument("--min-rc", type=float, default=0.5,
                   help="Reject response constants below this (ADC/V) — too small to trust")
    p.add_argument("--max-rc", type=float, default=20.0,
                   help="Reject response constants above this (ADC/V) — likely a glitch")
    p.add_argument("--min-good-fit", action="store_true",
                   help="Only propose for channels whose latest entry has status==OK")
    p.add_argument("--max-sigma", type=float, default=0.0,
                   help="Skip channels whose latest <metric>_sigma exceeds this (0 = no cut). "
                        "Applies to whichever metric --metric selects.")
    p.add_argument("--hv-url", default="http://clonpc19:8765",
                   help="prad2hvd URL")
    p.add_argument("--query-hv", action="store_true",
                   help="Query prad2hvd for current VSet/crate addressing instead of using "
                        "the values stored in the history file")
    p.add_argument("--apply", action="store_true",
                   help="POST the batch JSON directly via /api/load_settings (asks for password)")
    p.add_argument("--no-limits", action="store_true",
                   help="Skip voltage_limits.json clamping (default: clamp)")
    p.add_argument("--exclude-list", default=None,
                   help="Path to a text file of module names to leave alone "
                        "(one per line, '#' comments allowed). These channels "
                        "are reported but never get a new VSet.")
    p.add_argument("--modules-db", default=None)
    return p.parse_args()


def load_exclude_list(path: str) -> set:
    """Read a text file of module names; return a set.

    Format: one name per line. Blank lines and lines starting with '#' are
    ignored. Inline '#' comments are stripped. Whitespace-trimmed.
    """
    out: set = set()
    with open(path) as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            for token in line.split():
                out.add(token)
    return out


# ============================================================================
#  Proposal logic
# ============================================================================

class Proposal:
    __slots__ = ("name", "current_vset", "current_mean", "rc", "dv",
                 "new_vset", "reason", "iter_count")

    def __init__(self, name: str):
        self.name         = name
        self.current_vset: Optional[float] = None
        self.current_mean: Optional[float] = None
        self.rc:           Optional[float] = None
        self.dv:           Optional[float] = None
        self.new_vset:     Optional[float] = None
        self.reason:       str             = ""
        self.iter_count:   int             = 0


def metric_keys(metric: str) -> Tuple[str, str, str]:
    """Return the (mean, sigma, chi2) field names for the chosen metric."""
    if metric == "peak_integral":
        return "peak_integral_mean", "peak_integral_sigma", "peak_integral_chi2"
    return "peak_height_mean", "peak_height_sigma", "peak_height_chi2"


def propose_for_channel(name: str,
                        entries: List[dict],
                        args: argparse.Namespace,
                        v_limit: float,
                        exclude: set) -> Proposal:
    """Decide what new VSet to propose for one channel."""
    p = Proposal(name)
    p.iter_count = len(entries)
    if name in exclude:
        p.reason = "excluded by list"
        if entries:
            latest = entries[-1]
            p.current_vset = latest.get("VSet")
            mean_key, _, _ = metric_keys(args.metric)
            p.current_mean = latest.get(mean_key)
        return p
    if not entries:
        p.reason = "no entries"
        return p

    mean_key, sigma_key, _ = metric_keys(args.metric)
    latest = entries[-1]
    mean = latest.get(mean_key)
    vset = latest.get("VSet")

    if mean is None:
        p.reason = "no fit mean"
        return p
    if vset is None:
        p.reason = "no VSet recorded"
        return p

    if args.min_good_fit and latest.get("status") != "OK":
        p.reason = f"status={latest.get('status')}"
        return p

    sig = latest.get(sigma_key) or 0.0
    if args.max_sigma > 0 and sig > args.max_sigma:
        p.reason = f"sigma={sig:.1f} > {args.max_sigma:.1f}"
        return p

    p.current_mean = mean
    p.current_vset = vset

    diff = args.target_mean - mean   # > 0 → mean below target, < 0 → above
    if diff >= 0 and diff <= args.tolerance_low:
        p.reason = f"within low tolerance (diff={diff:+.1f})"
        return p
    if diff < 0 and -diff <= args.tolerance_high:
        p.reason = f"within high tolerance (diff={diff:+.1f})"
        return p

    # Decide ΔV
    if len(entries) < 2:
        # First-iteration constant step (signed by direction toward target).
        # If --initial-threshold is set, only step channels whose |diff| from
        # target exceeds the threshold (bidirectional): channels well below
        # target step UP, channels well above step DOWN, channels within
        # ± threshold of target are left at their current voltage.
        if args.initial_threshold > 0 and abs(diff) <= args.initial_threshold:
            p.reason = (f"|diff|={abs(diff):.1f} <= initial threshold "
                        f"{args.initial_threshold:.1f}")
            return p
        step = abs(args.initial_step)
        p.dv = step if diff > 0 else -step
        p.reason = "initial constant step"
    else:
        prev = entries[-2]
        m1 = prev.get(mean_key)
        v1 = prev.get("VSet")
        if m1 is None or v1 is None:
            # fall back to constant step
            step = abs(args.initial_step)
            p.dv = step if diff > 0 else -step
            p.reason = "prev entry missing data → constant step"
        else:
            dv_prev = vset - v1
            dm_prev = mean - m1
            if abs(dv_prev) < 1e-3:
                p.reason = "no ΔV between iterations"
                return p
            rc = dm_prev / dv_prev
            p.rc = rc
            if not (args.min_rc <= abs(rc) <= args.max_rc):
                p.reason = f"rc={rc:.2f} out of [{args.min_rc},{args.max_rc}]"
                return p
            dv = diff / rc
            # cap
            if dv >  args.max_dv: dv =  args.max_dv
            if dv < -args.max_dv: dv = -args.max_dv
            p.dv = dv
            p.reason = "linear-response correction"

    new_vset = vset + p.dv

    # Clamp to voltage limit
    if v_limit > 0 and new_vset > v_limit:
        new_vset = v_limit
        p.reason += f" (clamped to limit {v_limit:.0f}V)"
    if new_vset < 0:
        new_vset = 0.0
        p.reason += " (clamped at 0V)"

    p.new_vset = round(new_vset, 2)
    if abs(p.new_vset - vset) < 0.05:
        p.reason += " — no change after clamp"
        p.new_vset = None  # don't include in batch

    return p


# ============================================================================
#  Build batch JSON
# ============================================================================

def build_channel_entries(proposals: List[Proposal],
                          hv_lookup: Dict[str, dict],
                          channels: Dict[str, List[Dict[str, Any]]]
                         ) -> Tuple[List[Dict[str, Any]], int]:
    """Build the per-channel dicts that go into prad2hvmon_settings_v1.

    prad2hvd's /api/load_settings matches channels by HV crate/slot/channel
    addressing — NOT by name. Without these fields the daemon reports the
    entry as 'skipped'. We try, in order:
        1. live hv_lookup  (populated when --query-hv or --apply was used)
        2. hv_crate/hv_slot/hv_channel stored in the latest history entry
           (recorded by cosmic_eq_patch_hv or by cosmic_eq_analyze's HV query)

    Returns (entries, n_missing_addr).
    """
    out: List[Dict[str, Any]] = []
    n_missing_addr = 0
    for p in proposals:
        if p.new_vset is None:
            continue
        entry: Dict[str, Any] = {
            "name":   p.name,
            "params": {"V0Set": p.new_vset},
        }

        # source 1: live HV lookup
        hv_info = hv_lookup.get(p.name)
        crate = slot = chan = None
        if hv_info is not None:
            crate = hv_info.get("crate")
            slot  = hv_info.get("slot")
            chan  = hv_info.get("channel")

        # source 2: history-stored addressing (set by patch_hv / analyze)
        if crate is None or slot is None or chan is None:
            history = channels.get(p.name, [])
            if history:
                latest = history[-1]
                if crate is None: crate = latest.get("hv_crate")
                if slot  is None: slot  = latest.get("hv_slot")
                if chan  is None: chan  = latest.get("hv_channel")

        if crate is not None: entry["crate"]   = crate
        if slot  is not None: entry["slot"]    = slot
        if chan  is not None: entry["channel"] = chan

        if crate is None or slot is None or chan is None:
            n_missing_addr += 1

        out.append(entry)
    return out, n_missing_addr


# ============================================================================
#  Main
# ============================================================================

def main() -> int:
    args = parse_args()

    # Asymmetric tolerance defaults: fall back to symmetric --tolerance
    if args.tolerance_low is None:
        args.tolerance_low = args.tolerance
    if args.tolerance_high is None:
        args.tolerance_high = args.tolerance

    channels = load_history(args.history)
    if not channels:
        sys.exit(f"ERROR: no channels in {args.history}")

    modules = (load_hycal_modules(args.modules_db)
               if args.modules_db else load_hycal_modules())
    all_in_history = list(channels.keys())
    if args.types.lower() == "all":
        wanted = list(all_in_history)
    else:
        wanted_types = [t.strip() for t in args.types.split(",") if t.strip()]
        wanted = filter_modules_by_type(all_in_history, modules, wanted_types)
    wanted.sort(key=natural_module_sort_key)
    print(f"Considering {len(wanted)} modules of types: {args.types}")

    # voltage limits
    if args.no_limits:
        limits = []
    else:
        limits = load_voltage_limits()
        if not limits:
            print("WARNING: voltage_limits.json not found — proceeding without limits")

    # exclude list (channels with known-bad fits — touched VSet would be wrong)
    exclude: set = set()
    if args.exclude_list:
        exclude = load_exclude_list(args.exclude_list)
        in_scope = sum(1 for n in wanted if n in exclude)
        print(f"Loaded {len(exclude)} names from {args.exclude_list} "
              f"({in_scope} match the current type filter)")

    # Optional HV refresh
    hv_lookup: Dict[str, dict] = {}
    if args.query_hv or args.apply:
        hv = HVClient(args.hv_url)
        for name in wanted:
            try:
                info = hv.get_voltage(name)
                if info is not None:
                    hv_lookup[name] = info
                    if args.query_hv:
                        # overwrite the latest entry's VSet/VMon with live values
                        if channels.get(name):
                            channels[name][-1]["VSet"] = info.get("vset")
                            channels[name][-1]["VMon"] = info.get("vmon")
            except Exception as e:
                print(f"WARNING: HV GET {name} failed: {e}")

    # Build proposals
    proposals: List[Proposal] = []
    for name in wanted:
        v_lim = voltage_limit_for(name, limits) if limits else 0.0
        proposals.append(propose_for_channel(name, channels.get(name, []),
                                              args, v_lim, exclude))

    # Pretty print
    print()
    print(f"{'name':<8} {'iter':>4} {'cur_VSet':>9} {'cur_mean':>9} "
          f"{'rc':>6} {'dV':>7} {'new_VSet':>9}  reason")
    print("-" * 100)
    n_set = 0
    for p in proposals:
        cv  = f"{p.current_vset:.1f}" if p.current_vset is not None else "—"
        cm  = f"{p.current_mean:.1f}" if p.current_mean is not None else "—"
        rc  = f"{p.rc:+.2f}"           if p.rc           is not None else "—"
        dv  = f"{p.dv:+.1f}"           if p.dv           is not None else "—"
        nv  = f"{p.new_vset:.1f}"      if p.new_vset     is not None else "—"
        if p.new_vset is not None:
            n_set += 1
        print(f"{p.name:<8} {p.iter_count:>4} {cv:>9} {cm:>9} {rc:>6} {dv:>7} {nv:>9}  "
              f"{p.reason}")

    # Tally rejection reasons (helpful when n_set is unexpectedly small/zero)
    from collections import Counter
    tally: Counter = Counter()
    for p in proposals:
        if p.new_vset is not None:
            tally["proposed"] += 1
        else:
            # bucket the reason by its first phrase (everything before " (" or " — ")
            r = p.reason.split(" (")[0].split(" — ")[0].strip() or "unknown"
            tally[r] += 1
    print()
    print("Rejection breakdown:")
    for reason, count in tally.most_common():
        marker = "✓" if reason == "proposed" else " "
        print(f"  {marker} {count:>5}  {reason}")

    # Build batch JSON
    chan_entries, n_missing_addr = build_channel_entries(proposals, hv_lookup, channels)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    settings  = build_batch_settings(timestamp, chan_entries)

    with open(args.output, "w") as f:
        json.dump(settings, f, indent=2)
    print()
    print(f"Wrote {args.output}  —  {n_set} channels with new VSet")
    if n_missing_addr > 0:
        print(f"WARNING: {n_missing_addr} entries lack hv_crate/hv_slot/hv_channel — "
              f"prad2hvd /api/load_settings will SKIP them. "
              f"Run cosmic_eq_patch_hv.py first to populate addressing in the history "
              f"JSON, or pass --query-hv to fetch it live now.")

    # Optional direct apply
    if args.apply:
        if n_set == 0:
            print("Nothing to apply.")
            return 0
        password = os.environ.get("PRAD2HVD_PASSWORD") or getpass.getpass(
            "prad2hvd expert password: ")
        hv = HVClient(args.hv_url, password=password)
        try:
            granted = hv.authenticate()
        except Exception as e:
            sys.exit(f"ERROR: auth failed: {e}")
        if granted < 2:
            sys.exit(f"ERROR: insufficient privileges (granted={granted}, need 2)")
        try:
            resp = hv.load_settings(settings)
        except Exception as e:
            sys.exit(f"ERROR: /api/load_settings failed: {e}")
        print(f"prad2hvd response: {json.dumps(resp)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
