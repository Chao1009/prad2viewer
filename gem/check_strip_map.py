#!/usr/bin/env python3
"""
Cross-check strip mapping between three implementations:
  1. PRadAnalyzer  (PRad-I, SRS electronics, no hybrid board)
  2. mpd_gem_view_ssp (PRad-II, MPD electronics, hybrid board)
  3. Our implementation (GemSystem::buildStripMap, configurable)

Verifies that our config-driven pipeline produces identical plane-wide
strip numbers as the reference code for all 128 channels of every APV.

Usage:
    python check_strip_map.py [path/to/gem_map.json]
"""

import argparse
import json
import sys
import os
from gem_strip_map import map_strip

APV_SIZE = 128


# =========================================================================
# Reference implementations (hardcoded logic from original codebases)
# =========================================================================

def map_strip_pradanalyzer(ch, plane_type, plane_index, orient, plane_orient=1):
    """PRadAnalyzer PRadGEMAPV::MapStrip — SRS, no hybrid board."""
    strip = 32 * (ch % 4) + 8 * (ch // 4) - 31 * (ch // 16)

    if plane_type == 0 and plane_index == 11:
        strip = (48 - (strip + 1) // 2) if (strip & 1) else (48 + strip // 2)
    else:
        strip = (32 - (strip + 1) // 2) if (strip & 1) else (32 + strip // 2)
    strip &= 0x7f
    local = strip

    if orient != plane_orient:
        strip = 127 - strip

    if plane_type == 0 and plane_index == 11:
        strip += -16 + APV_SIZE * (plane_index - 1)
    else:
        strip += APV_SIZE * plane_index
    return local, strip


def map_strip_mpd(ch, plane_type, plane_index, orient):
    """mpd_gem_view_ssp GEMAPV::MapStripPRad — MPD with hybrid board."""
    strip = 32 * (ch % 4) + 8 * (ch // 4) - 31 * (ch // 16)
    strip = strip + 1 + strip % 4 - 5 * ((strip // 4) % 2)

    if plane_type == 0 and plane_index == 11:
        strip = (48 - (strip + 1) // 2) if (strip & 1) else (48 + strip // 2)
    else:
        strip = (32 - (strip + 1) // 2) if (strip & 1) else (32 + strip // 2)
    strip &= 0x7f
    local = strip

    if orient == 1:
        strip = 127 - strip

    if plane_type == 0 and plane_index == 11:
        strip += -16 + APV_SIZE * (plane_index - 1)
    else:
        strip += APV_SIZE * plane_index
    return local, strip


# =========================================================================
# Our implementation (config-driven, from GemSystem::buildStripMap)
# =========================================================================

def map_strip_ours(ch, plane_index, orient, pin_rotate=0, shared_pos=-1,
                   hybrid_board=True, apv_channels=128, readout_center=32):
    """Our implementation — delegates to shared gem_strip_map.map_strip."""
    return map_strip(ch, plane_index, orient,
                     pin_rotate=pin_rotate, shared_pos=shared_pos,
                     hybrid_board=hybrid_board,
                     apv_channels=apv_channels, readout_center=readout_center)


# =========================================================================
# Comparison logic
# =========================================================================

def check_apv(plane_type, plane_index, orient, pin_rotate=0, shared_pos=-1,
              hybrid_board=True, verbose=False):
    """Check all 128 channels. Returns (mpd_fail, prad_fail, details)."""
    mismatches_mpd = 0
    mismatches_prad = 0
    details = []

    for ch in range(APV_SIZE):
        _, mpd_plane = map_strip_mpd(ch, plane_type, plane_index, orient)
        _, our_plane = map_strip_ours(ch, plane_index, orient,
                                       pin_rotate=pin_rotate,
                                       shared_pos=shared_pos,
                                       hybrid_board=hybrid_board)
        if mpd_plane != our_plane:
            mismatches_mpd += 1
            if verbose and mismatches_mpd <= 5:
                details.append(f"    ch={ch}: mpd={mpd_plane} ours={our_plane} (diff={our_plane - mpd_plane})")

    # also check vs PRadAnalyzer (SRS = no hybrid board)
    for ch in range(APV_SIZE):
        _, prad_plane = map_strip_pradanalyzer(ch, plane_type, plane_index, orient)
        _, our_plane = map_strip_ours(ch, plane_index, orient,
                                       pin_rotate=pin_rotate,
                                       shared_pos=shared_pos,
                                       hybrid_board=False)
        if prad_plane != our_plane:
            mismatches_prad += 1

    return mismatches_mpd, mismatches_prad, details


def main():
    parser = argparse.ArgumentParser(
        description="Cross-check GEM strip mapping vs PRadAnalyzer + "
                    "mpd_gem_view_ssp reference implementations.")
    parser.add_argument("-G", "--gem-map", default=None,
                        help="GEM map JSON path (default: auto-search).")
    args = parser.parse_args()

    gem_map_path = args.gem_map
    if not gem_map_path:
        for candidate in ["database/gem_map.json",
                          "../database/gem_map.json",
                          "../../database/gem_map.json",
                          "gem_map.json"]:
            if os.path.exists(candidate):
                gem_map_path = candidate
                break
        else:
            print("Error: cannot find gem_map.json (pass -G <path>)")
            sys.exit(1)

    with open(gem_map_path, encoding="utf-8") as f:
        raw = json.load(f)

    apvs = [e for e in raw["apvs"] if "crate" in e]

    print(f"Checking {len(apvs)} APVs from {gem_map_path}")
    print(f"Comparing against mpd_gem_view_ssp (hybrid board) and PRadAnalyzer (SRS)")
    print()

    # ---- check all APVs from config ----
    total_mpd_fail = 0
    total_prad_fail = 0

    print(f"{'det':>4} {'plane':>6} {'pos':>4} {'orient':>7} {'pinrot':>7} {'shared':>7} {'match':>6}  {'vs MPD':>8} {'vs PRAna':>9}")
    print("-" * 72)

    for apv in apvs:
        det = apv["det"]
        plane_str = apv.get("plane", "X")
        plane_type = 1 if plane_str in ("Y", "1") else 0
        pos = apv["pos"]
        orient = apv["orient"]
        pin_rotate = apv.get("pin_rotate", 0)
        shared_pos = apv.get("shared_pos", -1)
        hybrid_board = apv.get("hybrid_board", True)
        match = apv.get("match", "")

        mpd_fail, prad_fail, details = check_apv(
            plane_type, pos, orient,
            pin_rotate=pin_rotate,
            shared_pos=shared_pos,
            hybrid_board=hybrid_board,
            verbose=True)

        mpd_status = "OK" if mpd_fail == 0 else f"FAIL({mpd_fail})"
        prad_status = "OK" if prad_fail == 0 else f"FAIL({prad_fail})"

        special = ""
        if pin_rotate != 0 or shared_pos >= 0 or match:
            special = " *"

        sp_str = str(shared_pos) if shared_pos >= 0 else "-"
        print(f"{det:4d} {plane_str:>6} {pos:4d} {orient:7d} {pin_rotate:7d} {sp_str:>7} {match:>6}  {mpd_status:>8} {prad_status:>9}{special}")

        for d in details:
            print(d)

        total_mpd_fail += mpd_fail
        total_prad_fail += prad_fail

    # ---- summary ----
    print()
    print("=" * 72)
    print(f"vs mpd_gem_view_ssp (MPD, hybrid board): ", end="")
    if total_mpd_fail == 0:
        print("ALL PASS")
    else:
        print(f"{total_mpd_fail} channel mismatches")

    print(f"vs PRadAnalyzer (SRS, no hybrid board):  ", end="")
    if total_prad_fail == 0:
        print("ALL PASS")
    else:
        print(f"{total_prad_fail} channel mismatches")
    print("=" * 72)

    # ---- also test edge cases not in config ----
    print("\nAdditional edge case checks:")

    edge_cases = [
        ("Y plane pos=0 orient=1",           1,  0, 1,  0, -1),
        ("Y plane pos=23 orient=1",          1, 23, 1,  0, -1),
        ("X pos=0 orient=0",                 0,  0, 0,  0, -1),
        ("X pos=9 orient=0 (last normal)",   0,  9, 0,  0, -1),
        ("X pos=10 orient=0 (hole neighbor)",0, 10, 0,  0, -1),
        ("X pos=11 orient=0 (pin_rot=16)",   0, 11, 0, 16, 10),
        ("X pos=11 orient=1 (pin_rot=16)",   0, 11, 1, 16, 10),
    ]

    all_edge_ok = True
    for label, pt, pi, ori, pr, sp in edge_cases:
        mpd_fail, _, details = check_apv(pt, pi, ori, pin_rotate=pr,
                                          shared_pos=sp, verbose=True)
        status = "OK" if mpd_fail == 0 else f"FAIL({mpd_fail})"
        print(f"  {label:45s} {status}")
        for d in details:
            print(d)
        if mpd_fail > 0:
            all_edge_ok = False

    print()
    if total_mpd_fail == 0 and all_edge_ok:
        print("ALL CHECKS PASSED")
    else:
        print("SOME CHECKS FAILED")
    return 0 if (total_mpd_fail == 0 and all_edge_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
