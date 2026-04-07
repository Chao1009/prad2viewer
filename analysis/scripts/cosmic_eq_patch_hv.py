#!/usr/bin/env python3
"""
cosmic_eq_patch_hv — fill in VMon/VSet for the latest iteration entry of each
channel in a gain history JSON.

Use case: a colleague sends you a JSON with fit results (peak_height_mean,
peak_integral_mean, count, ...) but no VMon/VSet because they don't have
network access to prad2hvd. This script queries prad2hvd for each channel
and writes the live VMon/VSet into the latest iteration entry, leaving all
other fields untouched. After patching, the file is ready for
cosmic_eq_propose.py.

By default the patch is idempotent: a channel whose latest entry already has
both VMon and VSet is skipped. Use --force to overwrite.

The script also tolerates and lightly normalises a few common variations:
  - missing 'iter'      → backfilled as the entry's 1-based position
  - missing 'timestamp' → added (now) only when patching that entry
  - missing 'status'    → set to "OK" only when patching, never overwritten

Default behavior writes back in place with a .bak alongside. Use --output to
write to a different file (then no .bak is created).

Typical use:
    python3 cosmic_eq_patch_hv.py colleague_results.json
    python3 cosmic_eq_patch_hv.py colleague_results.json --output gain_history.json
    python3 cosmic_eq_patch_hv.py gain_history.json --force --hv-url http://clonpc19:8765
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional

from cosmic_eq_common import HVClient, HISTORY_FORMAT, natural_module_sort_key


# ============================================================================
#  Argument parsing
# ============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Patch VMon/VSet into the latest iteration of a gain history JSON")
    p.add_argument("input", help="History JSON to patch")
    p.add_argument("-o", "--output", default=None,
                   help="Write to this path instead of in-place "
                        "(when omitted, input is rewritten and a .bak is created)")
    p.add_argument("--hv-url", default="http://clonpc19:8765",
                   help="prad2hvd URL")
    p.add_argument("--force", action="store_true",
                   help="Overwrite VMon/VSet even if already present")
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would change but do not write the file")
    p.add_argument("--no-format-stamp", action="store_true",
                   help="Do not add the 'format: cosmic_eq_history_v1' envelope when missing")
    return p.parse_args()


# ============================================================================
#  Loader / saver tolerant of un-enveloped colleague files
# ============================================================================

def load_tolerant(path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Load a history JSON. Accepts either:
       - the canonical {"format": "cosmic_eq_history_v1", "channels": {...}} envelope
       - a bare {channel_name: [iter, ...]} dict (colleague-style)
    Returns the channels dict in either case.
    """
    with open(path) as f:
        data = json.load(f)
    if not isinstance(data, dict):
        sys.exit(f"ERROR: {path} is not a JSON object")
    if "channels" in data and isinstance(data["channels"], dict):
        return data["channels"]
    # bare form — every value should be a list of entries
    bare_ok = all(isinstance(v, list) for v in data.values())
    if not bare_ok:
        sys.exit(f"ERROR: {path} is not a recognised history JSON layout")
    return data


def save_history(path: str,
                 channels: Dict[str, List[Dict[str, Any]]],
                 add_envelope: bool = True) -> None:
    if add_envelope:
        body = {"format": HISTORY_FORMAT, "channels": channels}
    else:
        body = channels
    with open(path, "w") as f:
        json.dump(body, f, indent=2)


# ============================================================================
#  Main
# ============================================================================

def main() -> int:
    args = parse_args()

    channels = load_tolerant(args.input)
    if not channels:
        sys.exit(f"ERROR: no channels found in {args.input}")

    # Set up output path & backup before any HV calls so we fail fast on bad paths.
    if args.dry_run:
        out_path = None
    elif args.output:
        out_path = args.output
    else:
        out_path = args.input
        bak = args.input + ".bak"
        shutil.copy2(args.input, bak)
        print(f"Backed up {args.input} → {bak}")

    hv = HVClient(args.hv_url)

    n_total    = 0
    n_patched  = 0
    n_skipped  = 0
    n_unchanged = 0
    n_missing  = 0

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    for name in sorted(channels.keys(), key=natural_module_sort_key):
        entries = channels[name]
        if not entries:
            continue
        n_total += 1
        latest = entries[-1]

        has_vset = latest.get("VSet") is not None
        has_vmon = latest.get("VMon") is not None
        has_addr = (latest.get("hv_crate") is not None
                    and latest.get("hv_slot") is not None
                    and latest.get("hv_channel") is not None)
        if has_vset and has_vmon and has_addr and not args.force:
            n_unchanged += 1
            continue

        try:
            info = hv.get_voltage(name)
        except Exception as e:
            print(f"  {name:<8}  HV GET failed: {e}")
            n_skipped += 1
            continue
        if info is None:
            print(f"  {name:<8}  not found in prad2hvd")
            n_missing += 1
            continue

        new_vset = info.get("vset")
        new_vmon = info.get("vmon")
        if new_vset is None:
            print(f"  {name:<8}  prad2hvd returned no vset")
            n_skipped += 1
            continue

        # Apply the patch
        old_vset = latest.get("VSet")
        old_vmon = latest.get("VMon")
        latest["VSet"] = round(float(new_vset), 2)
        if new_vmon is not None:
            latest["VMon"] = round(float(new_vmon), 2)
        # HV crate/slot/channel addressing is what /api/load_settings needs to
        # match channels — store it so cosmic_eq_propose can write it later.
        if "crate" in info:   latest["hv_crate"]   = info["crate"]
        if "slot" in info:    latest["hv_slot"]    = info["slot"]
        if "channel" in info: latest["hv_channel"] = info["channel"]
        latest.setdefault("status", "OK")
        latest.setdefault("timestamp", timestamp)
        n_patched += 1
        if old_vset is None:
            print(f"  {name:<8}  VSet={latest['VSet']:.1f}  VMon={latest.get('VMon', '—')}  "
                  f"addr={latest.get('hv_crate')}/{latest.get('hv_slot')}/{latest.get('hv_channel')}")
        else:
            print(f"  {name:<8}  VSet {old_vset} → {latest['VSet']}  "
                  f"VMon {old_vmon} → {latest.get('VMon', '—')}")

    # Backfill iter field for entries that lack it (1-based, by position).
    for name, entries in channels.items():
        for i, e in enumerate(entries):
            e.setdefault("iter", i + 1)

    print()
    print(f"Channels considered : {n_total}")
    print(f"  patched           : {n_patched}")
    print(f"  unchanged (had HV): {n_unchanged}")
    print(f"  skipped (HV err)  : {n_skipped}")
    print(f"  not in HV system  : {n_missing}")

    if args.dry_run:
        print("--dry-run: no file written")
        return 0

    save_history(out_path, channels,
                 add_envelope=not args.no_format_stamp)
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
