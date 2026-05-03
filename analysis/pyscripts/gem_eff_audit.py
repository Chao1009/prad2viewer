#!/usr/bin/env python3
"""
gem_eff_audit.py — offline audit of the GEM tracking-efficiency monitor.

Purpose
-------
The online monitor's `gem_eff_max_chi2 = 10` was tuned against a buggy
identity-transform code path (DetectorTransform cache wasn't being
invalidated when the runinfo block wrote the GEM poses).  After the fix,
the χ²/dof distribution and per-detector efficiency numbers are no longer
comparable to whatever was on screen before.  This script reads EVIO
files offline, runs the same line-fit algorithm, and reports:

  1. Distribution of χ²/dof for accepted tracks (PNG histogram).
  2. Per-detector efficiency in five modes — to expose the seed-bias
     that makes GEM 0/1 look ~95% and GEM 2/3 ~50%:
       * current:   seeds only from GEM 0 and GEM 1 (mirrors
                    AppState::runGemEfficiency).  Seed detector is
                    *automatically* matched, biasing its efficiency up.
       * all-seeds: seeds from every detector that has hits.
       * unbiased:  no seed at all.  For each HyCal cluster, enumerate
                    every GEM-subset of size ≥3 and every cluster combo
                    within those subsets, fit HyCal + combo, take the
                    lowest-χ²/dof fit that passes the gate.  At most one
                    track per HyCal cluster.  Per-detector efficiency =
                    tracks-with-this-detector / total-good-tracks.  This
                    is the right metric — every detector competes on
                    equal footing.
       * target:    seed line is (target → HyCal cluster); no GEM in the
                    seed.  Project to each GEM, take closest hit within
                    window, fit HyCal + matched GEMs (≥3).  Same per-
                    detector counting as `unbiased`.  Selects only
                    target-pointing tracks — comparing target vs
                    unbiased per detector estimates the upstream-halo
                    fraction.  Target xyz comes from runinfo (the same
                    "target" array the C++ replay uses).
       * loo:       leave-one-out — for each test detector D, fit the
                    track from HyCal + the OTHER three GEMs (≥2 of them
                    matched), project to D, count hit-or-not within
                    match window.  Independent unbiased estimator;
                    should agree with `unbiased` if both algorithms
                    are sound.

  Event filter (applied before any of the above): require ≥1 HyCal
  cluster AND ≥3 GEM detectors with at least one cluster — the minimum
  to fit a 4-parameter line with non-trivial dof.
  3. Sweep of efficiency vs `max_chi2_per_dof` so a sensible cut can be
     picked from the data instead of a guess.

All matching uses the same physics as the C++ monitor:
    σ_HC@gem = σ_HC(E_GeV) · |z_gem / z_hc|
    σ_total  = sqrt(σ_HC@gem² + σ_GEM[d]²)
    cut      = match_nsigma · σ_total

Usage
-----
    python analysis/pyscripts/gem_eff_audit.py <evio_path> <out_dir> \
        [--max-events N] [--match-nsigma 3.0] [--max-chi2 10.0] \
        [--max-hits-per-det 3]

`<evio_path>` accepts a glob, directory, or single split (same as the
other analysis scripts).  `<out_dir>` is created if missing; the script
writes `summary.txt`, `chi2_per_dof.png`, `efficiency_vs_chi2.png`, and
`residuals_loo.png`.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass, field
from itertools import combinations, product
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import _common as C
from _common import dec, det  # noqa: F401  (det is imported for future use)


# ---------------------------------------------------------------------------
# Geometry-aware tracking primitives — Python port of the C++ originals
# ---------------------------------------------------------------------------

def seed_line(x1: float, y1: float, z1: float,
              x2: float, y2: float, z2: float
              ) -> Tuple[float, float, float, float]:
    """Two-point line in lab frame: x(z)=ax+bx·z, y(z)=ay+by·z.
    Caller must ensure z1 != z2."""
    dz = z2 - z1
    if abs(dz) < 1e-6:
        return (x1, 0.0, y1, 0.0)
    bx = (x2 - x1) / dz
    ax = x1 - bx * z1
    by = (y2 - y1) / dz
    ay = y1 - by * z1
    return (ax, bx, ay, by)


def fit_weighted_line(z: Sequence[float], x: Sequence[float],
                      y: Sequence[float], w: Sequence[float]
                      ) -> Optional[Tuple[float, float, float, float, float]]:
    """Independent (z,x) and (z,y) weighted LSQ — 4-parameter line.
    Mirrors fitWeightedLine in app_state.cpp.  Returns
    (ax, bx, ay, by, chi2_per_dof) or None if the normal-equations
    determinant is degenerate.  dof = 2N − 4."""
    N = len(z)
    if N < 2:
        return None
    Sw = Sz = Szz = Sx = Sxz = Sy = Syz = 0.0
    for wi, zi, xi, yi in zip(w, z, x, y):
        Sw  += wi
        Sz  += wi * zi
        Szz += wi * zi * zi
        Sx  += wi * xi
        Sxz += wi * xi * zi
        Sy  += wi * yi
        Syz += wi * yi * zi
    delta = Sw * Szz - Sz * Sz
    if abs(delta) < 1e-9:
        return None
    bx = (Sw * Sxz - Sz * Sx) / delta
    ax = (Sx - bx * Sz) / Sw
    by = (Sw * Syz - Sz * Sy) / delta
    ay = (Sy - by * Sz) / Sw
    dof = 2 * N - 4
    chi2 = 0.0
    if dof > 0:
        for wi, zi, xi, yi in zip(w, z, x, y):
            dxp = (ax + bx * zi) - xi
            dyp = (ay + by * zi) - yi
            chi2 += wi * (dxp * dxp + dyp * dyp)
        chi2_per_dof = chi2 / dof
    else:
        chi2_per_dof = 0.0
    return (ax, bx, ay, by, chi2_per_dof)


def find_closest(hits: Sequence[Tuple[float, float, float]],
                 pred_x: float, pred_y: float, cut: float
                 ) -> Tuple[int, float]:
    """Return (index, dr) of the lab-frame hit closest to (pred_x, pred_y)
    within `cut`, or (-1, +inf) if none qualify."""
    best_idx = -1
    best_d2 = cut * cut
    for i, h in enumerate(hits):
        dx = h[0] - pred_x
        dy = h[1] - pred_y
        d2 = dx * dx + dy * dy
        if d2 < best_d2:
            best_d2 = d2
            best_idx = i
    if best_idx < 0:
        return -1, float("inf")
    return best_idx, math.sqrt(best_d2)


# ---------------------------------------------------------------------------
# Tracking-efficiency runs
# ---------------------------------------------------------------------------

@dataclass
class TrackingParams:
    match_nsigma:      float
    max_chi2:          float
    max_hits_per_det:  int
    gem_pos_res:       List[float]   # per-detector (mm), len 4

@dataclass
class TrackResult:
    chi2_per_dof: float
    matched:      List[bool]                      # len n_dets
    cand:         List[Optional[Tuple[float, float, float]]]
    fit:          Tuple[float, float, float, float]   # ax, bx, ay, by


def _try_seed(seed_d: int, seed_idx: int,
              hcx: float, hcy: float, hcz: float,
              sigma_hc: float,
              hits_by_det: List[List[Tuple[float, float, float]]],
              gem_z: Sequence[float],
              candidate_dets: Sequence[int],
              params: TrackingParams,
              min_match: int) -> Optional[TrackResult]:
    """Build a track using seed_d/seed_idx.  Match candidates only on the
    detectors listed in `candidate_dets` (always include seed_d).  Need at
    least `min_match` matched detectors for a valid fit (counts seed)."""
    if seed_d not in candidate_dets:
        return None
    hits_s = hits_by_det[seed_d]
    if seed_idx >= len(hits_s):
        return None
    g0 = hits_s[seed_idx]

    ax_s, bx_s, ay_s, by_s = seed_line(hcx, hcy, hcz, g0[0], g0[1], g0[2])

    n_dets = len(hits_by_det)
    matched = [False] * n_dets
    cand: List[Optional[Tuple[float, float, float]]] = [None] * n_dets
    matched[seed_d] = True
    cand[seed_d] = g0

    for d in candidate_dets:
        if d == seed_d:
            continue
        zd = gem_z[d]
        pred_x = ax_s + bx_s * zd
        pred_y = ay_s + by_s * zd
        s_hc_at_gem = sigma_hc * abs(zd / hcz) if hcz != 0 else sigma_hc
        s_gem = params.gem_pos_res[d] if d < len(params.gem_pos_res) else 0.1
        s_total = math.sqrt(s_hc_at_gem * s_hc_at_gem + s_gem * s_gem)
        cut = params.match_nsigma * s_total
        idx, _ = find_closest(hits_by_det[d], pred_x, pred_y, cut)
        if idx >= 0:
            matched[d] = True
            cand[d] = hits_by_det[d][idx]

    if sum(matched) < min_match:
        return None

    # Weighted fit of HyCal + every matched GEM.
    z_arr: List[float] = [hcz]
    x_arr: List[float] = [hcx]
    y_arr: List[float] = [hcy]
    w_arr: List[float] = [1.0 / (sigma_hc * sigma_hc)]
    for d in range(n_dets):
        if not matched[d] or cand[d] is None:
            continue
        h = cand[d]
        z_arr.append(h[2])
        x_arr.append(h[0])
        y_arr.append(h[1])
        s_gem = params.gem_pos_res[d] if d < len(params.gem_pos_res) else 0.1
        w_arr.append(1.0 / (s_gem * s_gem))

    fit = fit_weighted_line(z_arr, x_arr, y_arr, w_arr)
    if fit is None:
        return None
    ax, bx, ay, by, chi2_per_dof = fit
    if chi2_per_dof > params.max_chi2:
        return None
    return TrackResult(chi2_per_dof, matched, cand, (ax, bx, ay, by))


def best_track(hcx: float, hcy: float, hcz: float, sigma_hc: float,
               hits_by_det: List[List[Tuple[float, float, float]]],
               gem_z: Sequence[float], seed_dets: Sequence[int],
               candidate_dets: Sequence[int], params: TrackingParams,
               min_match: int) -> Optional[TrackResult]:
    """Try every (seed_d, seed_idx) pair.  Return the lowest-χ²/dof
    TrackResult that passes the chi²_max gate, or None."""
    best: Optional[TrackResult] = None
    for seed_d in seed_dets:
        if seed_d >= len(hits_by_det):
            continue
        n_seeds = min(len(hits_by_det[seed_d]), params.max_hits_per_det)
        for seed_idx in range(n_seeds):
            r = _try_seed(seed_d, seed_idx, hcx, hcy, hcz, sigma_hc,
                          hits_by_det, gem_z, candidate_dets, params,
                          min_match)
            if r is None:
                continue
            if best is None or r.chi2_per_dof < best.chi2_per_dof:
                best = r
    return best


def best_track_target_seed(hcx: float, hcy: float, hcz: float,
                           target_x: float, target_y: float, target_z: float,
                           sigma_hc: float,
                           hits_by_det: List[List[Tuple[float, float, float]]],
                           gem_z: Sequence[float],
                           params: TrackingParams,
                           min_match: int = 3) -> Optional[TrackResult]:
    """Target-seeded tracker — no GEM in the seed line.

    Seed line goes from (target_x, target_y, target_z) to the HyCal
    cluster.  Project to each GEM plane, take the closest hit within
    `match_nsigma · σ_total`, then fit HyCal + matched GEMs and apply the
    χ²/dof gate.  Need ≥`min_match` GEMs to call it a "good track".

    Per-detector efficiency built on top of this is unbiased *and*
    selects only target-pointing tracks — i.e., it rejects upstream
    halo / vacuum-interaction tracks that the combinatorial unbiased
    mode would still accept.  Comparing the two modes per detector is
    a clean way to estimate the halo fraction."""
    n_dets = len(hits_by_det)
    ax_s, bx_s, ay_s, by_s = seed_line(target_x, target_y, target_z,
                                       hcx, hcy, hcz)

    # σ_HC at target = 0 (target is geometric origin), σ_HC at HyCal =
    # sigma_hc.  Linear interpolation: at z_gem the projected positional
    # uncertainty scales with the lever arm from target.
    lever_hc = (hcz - target_z) if hcz != target_z else 1.0

    matched = [False] * n_dets
    cand: List[Optional[Tuple[float, float, float]]] = [None] * n_dets
    for d in range(n_dets):
        zd = gem_z[d]
        pred_x = ax_s + bx_s * zd
        pred_y = ay_s + by_s * zd
        s_hc_at_gem = sigma_hc * abs((zd - target_z) / lever_hc)
        s_gem = params.gem_pos_res[d] if d < len(params.gem_pos_res) else 0.1
        s_total = math.sqrt(s_hc_at_gem * s_hc_at_gem + s_gem * s_gem)
        cut = params.match_nsigma * s_total
        idx, _ = find_closest(hits_by_det[d], pred_x, pred_y, cut)
        if idx >= 0:
            matched[d] = True
            cand[d] = hits_by_det[d][idx]

    if sum(matched) < min_match:
        return None

    z_arr: List[float] = [hcz]
    x_arr: List[float] = [hcx]
    y_arr: List[float] = [hcy]
    w_arr: List[float] = [1.0 / (sigma_hc * sigma_hc)]
    for d in range(n_dets):
        if not matched[d] or cand[d] is None:
            continue
        h = cand[d]
        z_arr.append(h[2])
        x_arr.append(h[0])
        y_arr.append(h[1])
        s_gem = params.gem_pos_res[d] if d < len(params.gem_pos_res) else 0.1
        w_arr.append(1.0 / (s_gem * s_gem))

    fit = fit_weighted_line(z_arr, x_arr, y_arr, w_arr)
    if fit is None:
        return None
    ax, bx, ay, by, chi2_per_dof = fit
    if chi2_per_dof > params.max_chi2:
        return None
    return TrackResult(chi2_per_dof, matched, cand, (ax, bx, ay, by))


def best_track_unbiased(hcx: float, hcy: float, hcz: float, sigma_hc: float,
                        hits_by_det: List[List[Tuple[float, float, float]]],
                        gem_z: Sequence[float], params: TrackingParams,
                        min_gems: int = 3) -> Optional[TrackResult]:
    """Combinatorial best-track finder — no seed, no projection-and-match.
    For every subset of GEM detectors with size in [min_gems, n_dets],
    enumerate one cluster from each detector in the subset (capped at
    params.max_hits_per_det per detector), fit HyCal + that combo with
    fitWeightedLine, and return the lowest-χ²/dof fit that passes the
    chi²_max gate.

    Per-detector efficiency built on top of this is unbiased: every
    detector competes on equal footing for inclusion in the best track,
    so no detector gets the "automatic match" bonus the seed has in
    AppState::runGemEfficiency."""
    n_dets = len(hits_by_det)
    # Per-detector candidate lists, capped.
    capped: List[List[Tuple[float, float, float]]] = [
        hits_by_det[d][:params.max_hits_per_det] for d in range(n_dets)
    ]
    available = [d for d in range(n_dets) if capped[d]]
    if len(available) < min_gems:
        return None

    inv_sigma2_hc = 1.0 / (sigma_hc * sigma_hc)
    best: Optional[TrackResult] = None

    # Subsets of size {min_gems, …, len(available)} drawn from `available`.
    for k in range(min_gems, len(available) + 1):
        for det_subset in combinations(available, k):
            cluster_lists = [capped[d] for d in det_subset]
            # Pre-compute weights for this subset (constant across combos).
            w_gem = [1.0 / (params.gem_pos_res[d] *
                            params.gem_pos_res[d])
                     for d in det_subset]
            for combo in product(*cluster_lists):
                z = [hcz]
                x = [hcx]
                y = [hcy]
                for h in combo:
                    z.append(h[2]); x.append(h[0]); y.append(h[1])
                w = [inv_sigma2_hc] + w_gem
                fit = fit_weighted_line(z, x, y, w)
                if fit is None:
                    continue
                ax, bx, ay, by, chi2 = fit
                if chi2 > params.max_chi2:
                    continue
                if best is not None and chi2 >= best.chi2_per_dof:
                    continue
                matched = [False] * n_dets
                cand: List[Optional[Tuple[float, float, float]]] = (
                    [None] * n_dets)
                for i, d in enumerate(det_subset):
                    matched[d] = True
                    cand[d] = combo[i]
                best = TrackResult(chi2, matched, cand, (ax, bx, ay, by))
    return best


# ---------------------------------------------------------------------------
# Per-event accumulator
# ---------------------------------------------------------------------------

@dataclass
class ModeStats:
    name:           str
    n_tracks:       int = 0
    chi2_list:      List[float] = field(default_factory=list)
    num:            List[int]   = field(default_factory=lambda: [0]*4)
    den:            List[int]   = field(default_factory=lambda: [0]*4)


@dataclass
class LooStats:
    n_attempted:    List[int]    = field(default_factory=lambda: [0]*4)
    n_matched:      List[int]    = field(default_factory=lambda: [0]*4)
    residuals_x:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    residuals_y:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    C.add_common_args(ap)
    ap.add_argument("--match-nsigma",     type=float, default=3.0,
                    help="Matching window in σ_total (mirrors monitor "
                         "default).")
    ap.add_argument("--max-chi2",         type=float, default=10.0,
                    help="χ²/dof gate for accepting a fit.  The script "
                         "also sweeps a range below this in the "
                         "efficiency-vs-cut plot.")
    ap.add_argument("--max-hits-per-det", type=int,   default=3,
                    help="Cap seed/match candidates per detector (mirrors "
                         "AppState::gem_eff_max_hits_per_det = 3).")
    ap.add_argument("--n-dets",           type=int,   default=4,
                    help="Number of GEM detectors to consider (default 4).")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    p = C.setup_pipeline(
        evio_path     = args.evio_path,
        max_events    = args.max_events,
        run_num       = args.run_num,
        gem_ped_file  = args.gem_ped_file,
        gem_cm_file   = args.gem_cm_file,
        hc_calib_file = args.hc_calib_file,
        daq_config    = args.daq_config,
        gem_map_file  = args.gem_map_file,
        hc_map_file   = args.hc_map_file,
    )

    (pr_A, pr_B, pr_C), gem_pos_res = C.load_matching_config()
    print(f"[setup] HC σ(E) = sqrt(({pr_A:.3f}/√E_GeV)²"
          f"+({pr_B:.3f}/E_GeV)²+{pr_C:.3f}²) mm")
    print(f"[setup] σ_GEM   = {gem_pos_res} mm")
    print(f"[setup] match_nsigma={args.match_nsigma}  "
          f"max_chi2={args.max_chi2}  max_hits_per_det={args.max_hits_per_det}")

    n_dets = min(args.n_dets, 4, p.gem_sys.get_n_detectors())
    gem_z = list(p.geo.gem_z[:n_dets])
    print(f"[setup] GEM lab z = {gem_z}")
    print(f"[setup] target    = ({p.geo.target_x:g}, {p.geo.target_y:g}, "
          f"{p.geo.target_z:g}) mm  (from runinfo)")

    params = TrackingParams(
        match_nsigma     = args.match_nsigma,
        max_chi2         = args.max_chi2,
        max_hits_per_det = args.max_hits_per_det,
        gem_pos_res      = gem_pos_res,
    )

    mode_current  = ModeStats("current (seed=0,1)")
    mode_allseeds = ModeStats("all-seeds (0,1,2,3)")
    mode_unbiased = ModeStats("unbiased (HyCal + ≥3 GEMs, combinatorial)")
    mode_target   = ModeStats(
        f"target-seed (target=({p.geo.target_x:g},{p.geo.target_y:g},"
        f"{p.geo.target_z:g}) → HyCal)")
    loo           = LooStats()
    n_events_filter_pass = 0

    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    t0 = time.monotonic()
    n_phys = n_used = 0

    for fpath in p.evio_files:
        if ch.open_auto(fpath) != dec.Status.success:
            print(f"[WARN] skip (cannot open): {fpath}")
            continue
        print(f"[file] {fpath}")
        done = False
        while ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                decoded = ch.decode_event(i, with_ssp=True)
                if not decoded["ok"]:
                    continue
                n_phys += 1
                if args.max_events and n_phys >= args.max_events:
                    done = True

                fadc_evt = decoded["event"]
                ssp_evt  = decoded["ssp"]
                if fadc_evt.info.trigger_bits != C.PHYSICS_TRIGGER_BITS:
                    continue

                # HyCal clusters → lab.
                hc_raw = C.reconstruct_hycal(p, fadc_evt)
                if not hc_raw:
                    continue
                hc_lab: List[Tuple[float, float, float, float]] = []
                for h in hc_raw:
                    z_local = det.shower_depth(h.center_id, h.energy)
                    x, y, z = p.hycal_xform.to_lab(h.x, h.y, z_local)
                    hc_lab.append((x, y, z, float(h.energy)))

                # GEM hits → lab, per detector.
                C.reconstruct_gem(p, ssp_evt)
                hits_by_det: List[List[Tuple[float, float, float]]] = [
                    [] for _ in range(n_dets)
                ]
                for d in range(n_dets):
                    xform = p.gem_xforms[d]
                    for g in p.gem_sys.get_hits(d):
                        x, y, z = xform.to_lab(g.x, g.y)
                        hits_by_det[d].append((x, y, z))

                # Event filter — at least 3 GEM detectors must each have
                # ≥1 cluster.  Mirrors the user's gating spec; a track
                # needs 4 points (HyCal + ≥3 GEMs) to be a non-trivial
                # 4-parameter fit (dof = 2N − 4 ≥ 2).
                n_dets_with_hits = sum(1 for h in hits_by_det if h)
                if n_dets_with_hits < 3:
                    continue
                n_events_filter_pass += 1

                all_dets = list(range(n_dets))
                for hcx, hcy, hcz, energy in hc_lab:
                    if hcz <= 0:
                        continue
                    sigma_hc = C.hycal_pos_resolution(pr_A, pr_B, pr_C, energy)
                    n_used += 1

                    # current = seeds from {0, 1} only
                    cur = best_track(hcx, hcy, hcz, sigma_hc, hits_by_det,
                                     gem_z, seed_dets=[0, 1],
                                     candidate_dets=all_dets,
                                     params=params, min_match=3)
                    if cur is not None:
                        mode_current.n_tracks += 1
                        mode_current.chi2_list.append(cur.chi2_per_dof)
                        for d in range(n_dets):
                            mode_current.den[d] += 1
                            if cur.matched[d]:
                                mode_current.num[d] += 1

                    # all-seeds = seeds from every detector with hits
                    seeds_all = [d for d in range(n_dets)
                                 if hits_by_det[d]]
                    al = best_track(hcx, hcy, hcz, sigma_hc, hits_by_det,
                                    gem_z, seed_dets=seeds_all,
                                    candidate_dets=all_dets,
                                    params=params, min_match=3)
                    if al is not None:
                        mode_allseeds.n_tracks += 1
                        mode_allseeds.chi2_list.append(al.chi2_per_dof)
                        for d in range(n_dets):
                            mode_allseeds.den[d] += 1
                            if al.matched[d]:
                                mode_allseeds.num[d] += 1

                    # unbiased = HyCal + every k-subset of GEMs (k≥3).
                    # No seed, no projection-and-match.  At most one track
                    # per HyCal cluster; per-detector efficiency is
                    # tracks-with-this-detector / total-good-tracks.
                    ub = best_track_unbiased(hcx, hcy, hcz, sigma_hc,
                                             hits_by_det, gem_z, params,
                                             min_gems=3)
                    if ub is not None:
                        mode_unbiased.n_tracks += 1
                        mode_unbiased.chi2_list.append(ub.chi2_per_dof)
                        for d in range(n_dets):
                            mode_unbiased.den[d] += 1
                            if ub.matched[d]:
                                mode_unbiased.num[d] += 1

                    # target-seed = (target → HyCal) line, no GEM in seed.
                    # Per-detector counting is the same as `unbiased` —
                    # numerator is matched-in-good-track, denominator is
                    # good tracks.  Selects target-pointing tracks only.
                    ts = best_track_target_seed(hcx, hcy, hcz,
                            p.geo.target_x, p.geo.target_y, p.geo.target_z,
                            sigma_hc, hits_by_det, gem_z, params,
                            min_match=3)
                    if ts is not None:
                        mode_target.n_tracks += 1
                        mode_target.chi2_list.append(ts.chi2_per_dof)
                        for d in range(n_dets):
                            mode_target.den[d] += 1
                            if ts.matched[d]:
                                mode_target.num[d] += 1

                    # leave-one-out — for each test detector D, fit using
                    # the other 3 only (≥2 of them matched), then probe D.
                    for test_d in range(n_dets):
                        others = [d for d in range(n_dets) if d != test_d]
                        seeds_for_loo = [d for d in others if hits_by_det[d]]
                        if not seeds_for_loo:
                            continue
                        track_loo = best_track(hcx, hcy, hcz, sigma_hc,
                            hits_by_det, gem_z,
                            seed_dets=seeds_for_loo,
                            candidate_dets=others,
                            params=params, min_match=2)
                        if track_loo is None:
                            continue
                        ax, bx, ay, by = track_loo.fit
                        zd = gem_z[test_d]
                        pred_x = ax + bx * zd
                        pred_y = ay + by * zd
                        s_hc_at_gem = sigma_hc * abs(zd / hcz)
                        s_gem = params.gem_pos_res[test_d]
                        s_total = math.sqrt(s_hc_at_gem ** 2 + s_gem ** 2)
                        cut = params.match_nsigma * s_total
                        idx, _ = find_closest(hits_by_det[test_d],
                                              pred_x, pred_y, cut)
                        loo.n_attempted[test_d] += 1
                        if idx >= 0:
                            loo.n_matched[test_d] += 1
                            h = hits_by_det[test_d][idx]
                            loo.residuals_x[test_d].append(h[0] - pred_x)
                            loo.residuals_y[test_d].append(h[1] - pred_y)

            if done:
                break
        ch.close()
        if done:
            break

    elapsed = time.monotonic() - t0
    print(f"[done] {n_phys} physics events, "
          f"{n_events_filter_pass} pass filter (≥1 HyCal + ≥3 GEMs), "
          f"{n_used} HyCal clusters used, {elapsed:.1f}s")

    write_outputs(out_dir,
                  mode_current, mode_allseeds, mode_unbiased, mode_target,
                  loo, params, n_phys, n_events_filter_pass, n_used)
    return 0


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _eff(num: int, den: int) -> str:
    if den == 0:
        return "—"
    return f"{100.0 * num / den:5.1f}%  ({num}/{den})"


def write_outputs(out_dir: Path,
                  cur: ModeStats, allseeds: ModeStats, unbiased: ModeStats,
                  target: ModeStats,
                  loo: LooStats, params: TrackingParams,
                  n_phys: int, n_events_filter_pass: int, n_used: int
                  ) -> None:

    # ---- text summary (stdout) --------------------------------------------
    print()
    print("GEM tracking-efficiency audit")
    print("=============================")
    print()
    print(f"physics events processed : {n_phys}")
    print(f"events passing filter    : {n_events_filter_pass}  "
          f"(≥1 HyCal cluster AND ≥3 GEM detectors with hits)")
    print(f"HyCal clusters considered: {n_used}")
    print(f"match_nsigma             : {params.match_nsigma}")
    print(f"max_chi2_per_dof         : {params.max_chi2}")
    print(f"max_hits_per_det         : {params.max_hits_per_det}")
    print(f"σ_GEM (mm)               : {params.gem_pos_res}")
    print()

    for mode in (cur, allseeds, unbiased, target):
        print(f"--- {mode.name} ---")
        print(f"  tracks accepted: {mode.n_tracks}")
        for d in range(4):
            print(f"  GEM{d}: {_eff(mode.num[d], mode.den[d])}")
        if mode.chi2_list:
            arr = sorted(mode.chi2_list)
            med  = arr[len(arr)//2]
            p90  = arr[int(0.9 * len(arr))]
            p99  = arr[int(0.99 * len(arr))]
            print(f"  χ²/dof: median={med:.3f}  "
                  f"90%={p90:.3f}  99%={p99:.3f}")
        print()

    print(f"--- leave-one-out (unbiased per-detector) ---")
    print(f"  Each row: track built from the other 3 GEMs (≥2 matched),")
    print(f"  projected onto the test detector, then matched within "
          f"{params.match_nsigma}σ.")
    for d in range(4):
        print(f"  GEM{d}: {_eff(loo.n_matched[d], loo.n_attempted[d])}")

    # ---- plots (matplotlib optional) --------------------------------------
    plt = _import_pyplot()
    if plt is None:
        print("[plot] matplotlib not available; skipping PNGs")
        return

    _plot_chi2(plt, [cur, allseeds, unbiased, target],
               params, out_dir / "chi2_per_dof.png")
    _plot_eff_vs_chi2(plt, [cur, allseeds, unbiased, target],
                      params, out_dir / "efficiency_vs_chi2.png")
    _plot_residuals(plt, loo, out_dir / "residuals_loo.png")


def _import_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None


def _plot_chi2(plt, modes: List[ModeStats], params: TrackingParams,
               out: Path) -> None:
    if not any(m.chi2_list for m in modes):
        return
    import numpy as np
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, max(params.max_chi2, 1.0) * 1.1, 60)
    for i, m in enumerate(modes):
        if not m.chi2_list:
            continue
        ax.hist(m.chi2_list, bins=bins, color=f"C{i}", alpha=0.55,
                label=f"{m.name}  (n={m.n_tracks})")
    ax.axvline(params.max_chi2, color="k", ls="--", lw=1,
               label=f"current cut = {params.max_chi2}")
    ax.set_xlabel("χ²/dof")
    ax.set_ylabel("tracks")
    ax.set_title("GEM line-fit χ²/dof distribution")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _plot_eff_vs_chi2(plt, modes: List[ModeStats], params: TrackingParams,
                      out: Path) -> None:
    """Acceptance fraction vs χ²/dof cut, per mode.  Lets you see where
    tightening the gate starts to lose real tracks (curve flattens) vs
    where it's still cutting noise (curve still climbs)."""
    if not any(m.chi2_list for m in modes):
        return
    import numpy as np
    fig, ax = plt.subplots(figsize=(9, 5))
    cuts = np.linspace(0.5, max(params.max_chi2, 1.0) * 1.5, 60)
    for i, m in enumerate(modes):
        if not m.chi2_list:
            continue
        arr = np.asarray(m.chi2_list)
        accepted_at_cut = np.array(
            [(arr <= c).sum() for c in cuts], dtype=float)
        total = float(len(arr))
        ax.plot(cuts, 100.0 * accepted_at_cut / total, color=f"C{i}",
                lw=1.6, label=m.name)
    ax.axvline(params.max_chi2, color="k", ls="--", lw=1,
               label=f"current cut = {params.max_chi2}")
    ax.set_xlabel("χ²/dof cut")
    ax.set_ylabel("fraction of tracks kept (%)")
    ax.set_title("Track acceptance vs χ²/dof gate")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _plot_residuals(plt, loo: LooStats, out: Path) -> None:
    if not any(loo.residuals_x):
        return
    import numpy as np
    fig, axes = plt.subplots(2, 4, figsize=(16, 6), sharex='row')
    for d in range(4):
        for row, (label, arr) in enumerate(
            (("Δx (mm)", loo.residuals_x[d]),
             ("Δy (mm)", loo.residuals_y[d]))):
            ax = axes[row, d]
            if not arr:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes)
                ax.set_xlabel(label)
                ax.set_title(f"GEM{d}")
                continue
            data = np.asarray(arr)
            lo, hi = np.percentile(data, [1, 99])
            pad = max(1.0, (hi - lo) * 0.1)
            bins = np.linspace(lo - pad, hi + pad, 60)
            ax.hist(data, bins=bins, color=f"C{d}", alpha=0.8)
            med = float(np.median(data))
            std = float(np.std(data))
            ax.axvline(med, color="k", ls="--", lw=0.8,
                       label=f"med={med:.2f}\nσ={std:.2f}")
            ax.set_xlabel(label)
            ax.set_title(f"GEM{d}  (n={len(data)})")
            ax.legend(loc="best", fontsize=8)
            ax.grid(True, alpha=0.3)
    fig.suptitle("Leave-one-out residuals: hit − projected (test detector)")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


if __name__ == "__main__":
    sys.exit(main())
