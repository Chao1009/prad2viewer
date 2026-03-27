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

import json
import sys
import os

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

def map_strip_ours(ch, plane_index, orient, readout_offset=32, strip_offset=0,
                   hybrid_board=True):
    """Our implementation — configurable pipeline."""
    strip = 32 * (ch % 4) + 8 * (ch // 4) - 31 * (ch // 16)

    if hybrid_board:
        strip = strip + 1 + strip % 4 - 5 * ((strip // 4) % 2)

    if readout_offset > 0:
        if strip & 1:
            strip = readout_offset - (strip + 1) // 2
        else:
            strip = readout_offset + strip // 2
    strip &= 0x7f
    local = strip

    if orient == 1:
        strip = 127 - strip

    strip += strip_offset + plane_index * APV_SIZE
    return local, strip


# =========================================================================
# Comparison logic
# =========================================================================

def check_apv(plane_type, plane_index, orient, readout_offset=32,
              strip_offset=0, hybrid_board=True, verbose=False):
    """Check all 128 channels. Returns (ok, details_string)."""
    mismatches_mpd = 0
    mismatches_prad = 0
    details = []

    for ch in range(APV_SIZE):
        _, mpd_plane = map_strip_mpd(ch, plane_type, plane_index, orient)
        _, our_plane = map_strip_ours(ch, plane_index, orient,
                                       readout_offset=readout_offset,
                                       strip_offset=strip_offset,
                                       hybrid_board=hybrid_board)
        if mpd_plane != our_plane:
            mismatches_mpd += 1
            if verbose and mismatches_mpd <= 5:
                details.append(f"    ch={ch}: mpd={mpd_plane} ours={our_plane} (diff={our_plane - mpd_plane})")

    # also check vs PRadAnalyzer (plane-wide only, different local is expected)
    for ch in range(APV_SIZE):
        _, prad_plane = map_strip_pradanalyzer(ch, plane_type, plane_index, orient)
        _, our_plane = map_strip_ours(ch, plane_index, orient,
                                       readout_offset=readout_offset,
                                       strip_offset=strip_offset,
                                       hybrid_board=False)  # SRS = no hybrid
        if prad_plane != our_plane:
            mismatches_prad += 1

    return mismatches_mpd, mismatches_prad, details


def main():
    # find gem_map.json
    if len(sys.argv) > 1:
        gem_map_path = sys.argv[1]
    else:
        for candidate in ["database/gem_map.json", "../database/gem_map.json", "gem_map.json"]:
            if os.path.exists(candidate):
                gem_map_path = candidate
                break
        else:
            print("Usage: python check_strip_map.py [path/to/gem_map.json]")
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

    print(f"{'det':>4} {'plane':>6} {'pos':>4} {'orient':>7} {'ro_off':>7} {'st_off':>7} {'match':>6}  {'vs MPD':>8} {'vs PRAna':>9}")
    print("-" * 72)

    for apv in apvs:
        det = apv["det"]
        plane_str = apv.get("plane", "X")
        plane_type = 1 if plane_str in ("Y", "1") else 0
        pos = apv["pos"]
        orient = apv["orient"]
        readout_offset = apv.get("readout_offset", 32)
        strip_offset = apv.get("strip_offset", 0)
        hybrid_board = apv.get("hybrid_board", True)
        match = apv.get("match", "")

        mpd_fail, prad_fail, details = check_apv(
            plane_type, pos, orient,
            readout_offset=readout_offset,
            strip_offset=strip_offset,
            hybrid_board=hybrid_board,
            verbose=True)

        mpd_status = "OK" if mpd_fail == 0 else f"FAIL({mpd_fail})"
        prad_status = "OK" if prad_fail == 0 else f"FAIL({prad_fail})"

        special = ""
        if strip_offset != 0 or readout_offset != 32 or match:
            special = " *"

        print(f"{det:4d} {plane_str:>6} {pos:4d} {orient:7d} {readout_offset:7d} {strip_offset:7d} {match:>6}  {mpd_status:>8} {prad_status:>9}{special}")

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
        ("Y plane pos=0 orient=1", 1, 0, 1, 32, 0),
        ("Y plane pos=23 orient=1", 1, 23, 1, 32, 0),
        ("X pos=0 orient=0", 0, 0, 0, 32, 0),
        ("X pos=9 orient=0 (last normal)", 0, 9, 0, 32, 0),
        ("X pos=10 orient=0 (hole neighbor)", 0, 10, 0, 32, 0),
        ("X pos=11 orient=0 (correct config)", 0, 11, 0, 48, -144),
        ("X pos=11 orient=1 (correct config)", 0, 11, 1, 48, -144),
    ]

    all_edge_ok = True
    for label, pt, pi, ori, ro, so in edge_cases:
        mpd_fail, _, details = check_apv(pt, pi, ori, readout_offset=ro,
                                          strip_offset=so, verbose=True)
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
