#!/usr/bin/env python3
"""
gem_eff_audit.py — offline LOO audit of the GEM tracking-efficiency monitor.

Purpose
-------
For each test detector D in {0,1,2,3}, the OTHER 3 GEMs + HyCal define a
straight-line anchor track; we then project the anchor to D and count
whether D recorded a hit at the predicted position.  D itself contributes
nothing to the anchor (matching candidate set excludes D, fit excludes D),
so the test is genuinely unbiased.  Three LOO variants are run side-by-side
on every event so they can be compared:

  * loo              GEM-seeded anchor (HyCal + each OTHER GEM hit drawn as a
                     seed line; lowest-χ² survivor wins).  Fit through
                     HyCal + 3 OTHER GEMs.  Vertex-agnostic — accepts any
                     straight track that lights up the 3 OTHER GEMs.
  * loo-target-in    Same GEM-seeded matching, but the fit additionally
                     includes (target_x, target_y, target_z) as a weighted
                     measurement with σ = `--sigma-target` (≈ beam spot).
                     Pulls the line toward the target → kills upstream halo
                     and beam-gas tracks; the χ² gate then rejects anything
                     that doesn't actually point back to (0,0,0).
  * loo-target-seed  Single-seed variant: the seed line is
                     (target → HyCal cluster); no GEM-pair seeding.
                     Fit through HyCal + 3 OTHER GEMs (target only seeds,
                     never enters the fit).  Cheapest of the three.

Match definition (applied to the anchor in every variant):
  - all 3 OTHER detectors found a hit within match_nsigma · σ_total of the
    seed-line projection (σ_total = sqrt(σ_HC@gem² + σ_GEM²)),
  - weighted line fit's χ²/dof ≤ max_chi2,
  - per-detector residual against the FIT line within match_nsigma · σ_GEM[d].

Test (applied at D after the anchor passes):
  - D has a hit within match_nsigma · σ_GEM[D] of the projection.
  - Increment denominator[D] (anchor was good); increment numerator[D]
    if the test hit was found.

Event filter: ≥1 HyCal cluster AND ≥3 GEM detectors with at least one
cluster — minimum to make any anchor possible.

Usage
-----
    python analysis/pyscripts/gem_eff_audit.py <evio_path> <out_dir> \\
        [--match-nsigma 3.0] [--max-chi2 10.0] [--sigma-gem 0.5] \\
        [--sigma-target 1.0] [--min-cluster-energy 500] \\
        [--max-hits-per-det 3] [--max-events N]

`<evio_path>` accepts a glob (quote it), directory, or single split.
`<out_dir>` is created if missing; the script writes `anchor_chi2.png`
(two-row anchor-quality plot — χ²/dof distribution + cumulative
acceptance, NOT detector efficiency) and one `residuals_<variant>.png`
per LOO mode.  Detector efficiency numbers go to stdout in the text
summary.
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
                      y: Sequence[float],
                      w_x: Sequence[float],
                      w_y: Optional[Sequence[float]] = None
                      ) -> Optional[Tuple[float, float, float, float, float]]:
    """Independent (z,x) and (z,y) weighted LSQ — 4-parameter line.
    `w_y` defaults to `w_x` when both axes share the same per-point σ
    (the common case: GEM/HyCal hits, where σ_x = σ_y).  Pass distinct
    `w_y` to handle anisotropic uncertainties (e.g. the target point in
    `loo-target-in`, where σ_x and σ_y inherit different
    slope×σ_target_z couplings).  Returns (ax, bx, ay, by, chi2_per_dof)
    or None if the normal-equations determinant is degenerate.
    dof = 2N − 4."""
    N = len(z)
    if N < 2:
        return None
    if w_y is None:
        w_y = w_x
    # x-fit
    Swx = Szx = Szzx = Sx = Sxz = 0.0
    for wi, zi, xi in zip(w_x, z, x):
        Swx  += wi
        Szx  += wi * zi
        Szzx += wi * zi * zi
        Sx   += wi * xi
        Sxz  += wi * xi * zi
    delta_x = Swx * Szzx - Szx * Szx
    if abs(delta_x) < 1e-9:
        return None
    bx = (Swx * Sxz - Szx * Sx) / delta_x
    ax = (Sx - bx * Szx) / Swx
    # y-fit
    Swy = Szy = Szzy = Sy = Syz = 0.0
    for wi, zi, yi in zip(w_y, z, y):
        Swy  += wi
        Szy  += wi * zi
        Szzy += wi * zi * zi
        Sy   += wi * yi
        Syz  += wi * yi * zi
    delta_y = Swy * Szzy - Szy * Szy
    if abs(delta_y) < 1e-9:
        return None
    by = (Swy * Syz - Szy * Sy) / delta_y
    ay = (Sy - by * Szy) / Swy
    # chi2 (both axes contribute, possibly with different weights)
    dof = 2 * N - 4
    if dof > 0:
        chi2 = 0.0
        for wxi, wyi, zi, xi, yi in zip(w_x, w_y, z, x, y):
            dxp = (ax + bx * zi) - xi
            dyp = (ay + by * zi) - yi
            chi2 += wxi * dxp * dxp + wyi * dyp * dyp
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
# Match definition: a track is accepted when (a) the weighted line fit's
# χ²/dof passes the gate, AND (b) every matched detector's hit lies within
# `match_nsigma · σ_GEM[d]` of the fit line.  The seed-line projection
# window stays at `σ_total = sqrt(σ_HC@gem² + σ_GEM²)` because that's the
# uncertainty before the fit refines the line.  The post-fit gate uses just
# σ_GEM[d] because by then HyCal is already in the fit and the projection
# at any plane is dominated by the GEM resolution.
# ---------------------------------------------------------------------------

def fit_residuals_within_window(
        fit_params: Tuple[float, float, float, float],
        matched: Sequence[bool],
        cand: Sequence[Optional[Tuple[float, float, float]]],
        params: "TrackingParams",
        ) -> bool:
    ax, bx, ay, by = fit_params
    for d, m in enumerate(matched):
        if not m or cand[d] is None:
            continue
        h = cand[d]
        pred_x = ax + bx * h[2]
        pred_y = ay + by * h[2]
        s_gem = (params.gem_pos_res[d]
                 if d < len(params.gem_pos_res) else 0.1)
        cut = params.match_nsigma * s_gem
        dx = h[0] - pred_x
        dy = h[1] - pred_y
        if dx * dx + dy * dy > cut * cut:
            return False
    return True


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
              min_match: int,
              target_in_fit: bool = False,
              target_x: float = 0.0, target_y: float = 0.0,
              target_z: float = 0.0,
              sigma_target_x: float = 1.0,
              sigma_target_y: float = 1.0,
              sigma_target_z: float = 20.0,
              ) -> Optional[TrackResult]:
    """Build a track using seed_d/seed_idx.  Match candidates only on the
    detectors listed in `candidate_dets` (always include seed_d).  Need at
    least `min_match` matched detectors for a valid fit (counts seed).

    If `target_in_fit`, append (target_x, target_y, target_z) to the
    weighted fit as a soft "track originated at target" constraint.
    σ_target_z (the target's longitudinal extent) couples to the
    transverse measurement at z = target_z via the track slope:
        σ_x_eff² = σ_target_x² + (bx_est · σ_target_z)²
        σ_y_eff² = σ_target_y² + (by_est · σ_target_z)²
    where the slope estimate comes from the (target → HyCal cluster)
    line; for central tracks this leaves the σ_x,y contributions
    unchanged, while peripheral tracks (large |HyCal x,y|) get the
    expected σ_z lever-arm widening."""
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

    # Weighted fit of HyCal + every matched GEM (+ target if requested).
    # x and y use the same per-point σ except at the target: σ_target_z
    # couples differently into σ_x_eff and σ_y_eff via the slope.
    z_arr: List[float]  = [hcz]
    x_arr: List[float]  = [hcx]
    y_arr: List[float]  = [hcy]
    w_x_arr: List[float] = [1.0 / (sigma_hc * sigma_hc)]
    w_y_arr: List[float] = [1.0 / (sigma_hc * sigma_hc)]
    if target_in_fit:
        # Slope estimate from (target → HyCal cluster) — used to widen the
        # target's transverse σ by the σ_z lever arm.  Falls back to zero
        # slope when HyCal sits on the target plane.
        if hcz != target_z:
            bx_est = (hcx - target_x) / (hcz - target_z)
            by_est = (hcy - target_y) / (hcz - target_z)
        else:
            bx_est = by_est = 0.0
        sx_eff_sq = (sigma_target_x * sigma_target_x
                     + (bx_est * sigma_target_z) ** 2)
        sy_eff_sq = (sigma_target_y * sigma_target_y
                     + (by_est * sigma_target_z) ** 2)
        z_arr.append(target_z)
        x_arr.append(target_x)
        y_arr.append(target_y)
        w_x_arr.append(1.0 / sx_eff_sq)
        w_y_arr.append(1.0 / sy_eff_sq)
    for d in range(n_dets):
        if not matched[d] or cand[d] is None:
            continue
        h = cand[d]
        z_arr.append(h[2])
        x_arr.append(h[0])
        y_arr.append(h[1])
        s_gem = params.gem_pos_res[d] if d < len(params.gem_pos_res) else 0.1
        w_x_arr.append(1.0 / (s_gem * s_gem))
        w_y_arr.append(1.0 / (s_gem * s_gem))

    fit = fit_weighted_line(z_arr, x_arr, y_arr, w_x_arr, w_y_arr)
    if fit is None:
        return None
    ax, bx, ay, by, chi2_per_dof = fit
    if chi2_per_dof > params.max_chi2:
        return None
    if not fit_residuals_within_window((ax, bx, ay, by), matched, cand, params):
        return None
    return TrackResult(chi2_per_dof, matched, cand, (ax, bx, ay, by))


def best_track(hcx: float, hcy: float, hcz: float, sigma_hc: float,
               hits_by_det: List[List[Tuple[float, float, float]]],
               gem_z: Sequence[float], seed_dets: Sequence[int],
               candidate_dets: Sequence[int], params: TrackingParams,
               min_match: int,
               target_in_fit: bool = False,
               target_x: float = 0.0, target_y: float = 0.0,
               target_z: float = 0.0,
               sigma_target_x: float = 1.0,
               sigma_target_y: float = 1.0,
               sigma_target_z: float = 20.0,
               ) -> Optional[TrackResult]:
    """Try every (seed_d, seed_idx) pair.  Return the lowest-χ²/dof
    TrackResult that passes the chi²_max gate, or None.  All target_*
    args are forwarded to _try_seed unchanged."""
    best: Optional[TrackResult] = None
    for seed_d in seed_dets:
        if seed_d >= len(hits_by_det):
            continue
        n_seeds = min(len(hits_by_det[seed_d]), params.max_hits_per_det)
        for seed_idx in range(n_seeds):
            r = _try_seed(seed_d, seed_idx, hcx, hcy, hcz, sigma_hc,
                          hits_by_det, gem_z, candidate_dets, params,
                          min_match,
                          target_in_fit=target_in_fit,
                          target_x=target_x, target_y=target_y,
                          target_z=target_z,
                          sigma_target_x=sigma_target_x,
                          sigma_target_y=sigma_target_y,
                          sigma_target_z=sigma_target_z)
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
                           min_match: int = 3,
                           candidate_dets: Optional[Sequence[int]] = None
                           ) -> Optional[TrackResult]:
    """Target-seeded tracker — no GEM in the seed line.

    Seed line goes from (target_x, target_y, target_z) to the HyCal
    cluster.  Project to each candidate GEM plane, take the closest hit
    within `match_nsigma · σ_total`, then fit HyCal + matched GEMs and
    apply the χ²/dof gate.  Need ≥`min_match` GEMs to call it a "good
    track".  `candidate_dets` defaults to all GEM planes; the LOO mode
    passes the OTHER 3 detectors so the test detector is not in the fit.

    Per-detector efficiency built on top of this is unbiased *and*
    selects only target-pointing tracks — i.e., it rejects upstream
    halo / vacuum-interaction tracks that the combinatorial unbiased
    mode would still accept.  Comparing the two modes per detector is
    a clean way to estimate the halo fraction."""
    n_dets = len(hits_by_det)
    if candidate_dets is None:
        candidate_dets = range(n_dets)
    ax_s, bx_s, ay_s, by_s = seed_line(target_x, target_y, target_z,
                                       hcx, hcy, hcz)

    # σ_HC at target = 0 (target is geometric origin), σ_HC at HyCal =
    # sigma_hc.  Linear interpolation: at z_gem the projected positional
    # uncertainty scales with the lever arm from target.
    lever_hc = (hcz - target_z) if hcz != target_z else 1.0

    matched = [False] * n_dets
    cand: List[Optional[Tuple[float, float, float]]] = [None] * n_dets
    for d in candidate_dets:
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
    if not fit_residuals_within_window((ax, bx, ay, by), matched, cand, params):
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
                if not fit_residuals_within_window(
                        (ax, bx, ay, by), matched, cand, params):
                    continue
                best = TrackResult(chi2, matched, cand, (ax, bx, ay, by))
    return best


# ---------------------------------------------------------------------------
# Per-event accumulator — one set per LOO variant.  Each variant runs N times
# per HyCal cluster (one per test detector), so the denominator is per-detector
# (the "anchor exists" counter for that detector excluded), and the numerator
# is per-detector too.  chi2_list aggregates anchor χ²/dof across all test
# detectors.
# ---------------------------------------------------------------------------

@dataclass
class LooStats:
    name:           str
    n_attempted:    List[int]         = field(default_factory=lambda: [0]*4)
    n_matched:      List[int]         = field(default_factory=lambda: [0]*4)
    chi2_list:      List[float]       = field(default_factory=list)
    residuals_x:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    residuals_y:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])


def _record_loo(stats: "LooStats", test_d: int,
                track: Optional[TrackResult],
                hits_by_det: List[List[Tuple[float, float, float]]],
                gem_z: Sequence[float],
                params: TrackingParams) -> None:
    """Bookkeeping for a single LOO test attempt.  If `track` is None the
    anchor failed (3-GEM fit didn't pass χ²/residual gates); we don't
    increment any counter — denominator only counts events where the
    OTHER 3 GEMs delivered a clean anchor."""
    if track is None:
        return
    stats.chi2_list.append(track.chi2_per_dof)
    ax, bx, ay, by = track.fit
    zd = gem_z[test_d]
    pred_x = ax + bx * zd
    pred_y = ay + by * zd
    s_gem = (params.gem_pos_res[test_d]
             if test_d < len(params.gem_pos_res) else 0.1)
    cut = params.match_nsigma * s_gem
    idx, _ = find_closest(hits_by_det[test_d], pred_x, pred_y, cut)
    stats.n_attempted[test_d] += 1
    if idx >= 0:
        stats.n_matched[test_d] += 1
        h = hits_by_det[test_d][idx]
        stats.residuals_x[test_d].append(h[0] - pred_x)
        stats.residuals_y[test_d].append(h[1] - pred_y)


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
    ap.add_argument("--max-chi2",         type=float, default=3.5,
                    help="χ²/dof gate for accepting an anchor fit.  Default "
                         "3.5 — picks up real tracks while cutting the "
                         "anchor χ² tail of pile-up combinations.  The "
                         "anchor_chi2.png plot sweeps a range around this.")
    ap.add_argument("--max-hits-per-det", type=int,   default=3,
                    help="Cap seed/match candidates per detector (mirrors "
                         "AppState::gem_eff_max_hits_per_det = 3).")
    ap.add_argument("--min-cluster-energy", type=float, default=500.0,
                    help="Minimum HyCal cluster energy (MeV) to enter the "
                         "denominator.  Mirrors AppState::"
                         "gem_eff_min_cluster_energy.  Default 500 MeV.")
    ap.add_argument("--sigma-gem",        type=float, default=None,
                    help="Override σ_GEM (mm) for ALL detectors — useful for "
                         "absorbing residual alignment into a 'working' σ "
                         "before per-detector calibration is applied.  "
                         "Default: use reconstruction_config.json:matching:"
                         "gem_pos_res.")
    ap.add_argument("--sigma-target-x",   type=float, default=None,
                    help="σ_target_x (mm) — transverse beam-spot size in x. "
                         "Used by loo-target-in only.  Default: from "
                         "reconstruction_config.json:matching:target_pos_res.")
    ap.add_argument("--sigma-target-y",   type=float, default=None,
                    help="σ_target_y (mm) — transverse beam-spot size in y. "
                         "Default: from reconstruction_config.json.")
    ap.add_argument("--sigma-target-z",   type=float, default=None,
                    help="σ_target_z (mm) — target longitudinal extent.  "
                         "Couples to σ_x_eff and σ_y_eff via the track "
                         "slope at the target plane.  Default: from "
                         "reconstruction_config.json.")
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

    (pr_A, pr_B, pr_C), gem_pos_res, target_pos_res = C.load_matching_config()
    if args.sigma_gem is not None:
        gem_pos_res = [args.sigma_gem] * 4
    cfg_tgt_x, cfg_tgt_y, cfg_tgt_z = target_pos_res
    sigma_target_x = (args.sigma_target_x
                      if args.sigma_target_x is not None else cfg_tgt_x)
    sigma_target_y = (args.sigma_target_y
                      if args.sigma_target_y is not None else cfg_tgt_y)
    sigma_target_z = (args.sigma_target_z
                      if args.sigma_target_z is not None else cfg_tgt_z)
    print(f"[setup] HC σ(E) = sqrt(({pr_A:.3f}/√E_GeV)²"
          f"+({pr_B:.3f}/E_GeV)²+{pr_C:.3f}²) mm")
    print(f"[setup] σ_GEM   = {gem_pos_res} mm"
          + ("  (overridden)" if args.sigma_gem is not None else ""))
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
    print(f"[setup] min HyCal cluster E = {args.min_cluster_energy:.0f} MeV")

    overridden_target = any(a is not None for a in (
        args.sigma_target_x, args.sigma_target_y, args.sigma_target_z))
    print(f"[setup] σ_target = (x={sigma_target_x:.2f}, "
          f"y={sigma_target_y:.2f}, z={sigma_target_z:.2f}) mm  "
          f"(only used by loo-target-in"
          f"{', overridden' if overridden_target else ''})")

    loo            = LooStats("loo (GEM-seeded, HyCal + 3 OTHER GEMs in fit)")
    loo_target_in  = LooStats(
        f"loo-target-in (GEM-seeded, target+HyCal+3 OTHER GEMs in fit, "
        f"σ_target=(x={sigma_target_x:.2f},y={sigma_target_y:.2f},"
        f"z={sigma_target_z:.2f}) mm)")
    loo_target_seed = LooStats(
        f"loo-target-seed (target=({p.geo.target_x:g},{p.geo.target_y:g},"
        f"{p.geo.target_z:g})→HyCal seed, HyCal + 3 OTHER GEMs in fit)")
    loo_modes = (loo, loo_target_in, loo_target_seed)
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

                for hcx, hcy, hcz, energy in hc_lab:
                    if hcz <= 0:
                        continue
                    if energy < args.min_cluster_energy:
                        continue
                    sigma_hc = C.hycal_pos_resolution(pr_A, pr_B, pr_C, energy)
                    n_used += 1

                    # Three LOO variants run per HyCal cluster, per test
                    # detector D.  The OTHER 3 GEMs anchor the fit; D is
                    # excluded from both the candidate matching and the fit,
                    # so the prediction at D is genuinely unbiased.
                    #
                    #   loo              GEM-seeded  · fit = HyCal+3 GEMs
                    #   loo-target-in    GEM-seeded  · fit = target+HyCal+3 GEMs
                    #   loo-target-seed  target→HyCal seed · fit = HyCal+3 GEMs
                    for test_d in range(n_dets):
                        others = [d for d in range(n_dets) if d != test_d]
                        seeds_for_loo = [d for d in others if hits_by_det[d]]

                        # --- loo: GEM-seeded, no target ---
                        if seeds_for_loo:
                            track = best_track(hcx, hcy, hcz, sigma_hc,
                                hits_by_det, gem_z,
                                seed_dets=seeds_for_loo,
                                candidate_dets=others,
                                params=params, min_match=3)
                            _record_loo(loo, test_d, track,
                                        hits_by_det, gem_z, params)

                        # --- loo-target-in: GEM-seeded, target in fit ---
                        if seeds_for_loo:
                            track = best_track(hcx, hcy, hcz, sigma_hc,
                                hits_by_det, gem_z,
                                seed_dets=seeds_for_loo,
                                candidate_dets=others,
                                params=params, min_match=3,
                                target_in_fit=True,
                                target_x=p.geo.target_x,
                                target_y=p.geo.target_y,
                                target_z=p.geo.target_z,
                                sigma_target_x=sigma_target_x,
                                sigma_target_y=sigma_target_y,
                                sigma_target_z=sigma_target_z)
                            _record_loo(loo_target_in, test_d, track,
                                        hits_by_det, gem_z, params)

                        # --- loo-target-seed: target→HyCal seed, no target in fit ---
                        track = best_track_target_seed(hcx, hcy, hcz,
                            p.geo.target_x, p.geo.target_y, p.geo.target_z,
                            sigma_hc, hits_by_det, gem_z, params,
                            min_match=3, candidate_dets=others)
                        _record_loo(loo_target_seed, test_d, track,
                                    hits_by_det, gem_z, params)

            if done:
                break
        ch.close()
        if done:
            break

    elapsed = time.monotonic() - t0
    print(f"[done] {n_phys} physics events, "
          f"{n_events_filter_pass} pass filter (≥1 HyCal + ≥3 GEMs), "
          f"{n_used} HyCal clusters used, {elapsed:.1f}s")

    write_outputs(out_dir, loo_modes, params,
                  n_phys, n_events_filter_pass, n_used,
                  args.min_cluster_energy)
    return 0


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _eff(num: int, den: int) -> str:
    if den == 0:
        return "—"
    return f"{100.0 * num / den:5.1f}%  ({num}/{den})"


def write_outputs(out_dir: Path,
                  loo_modes: Sequence["LooStats"],
                  params: TrackingParams,
                  n_phys: int, n_events_filter_pass: int, n_used: int,
                  min_cluster_energy: float
                  ) -> None:

    # ---- text summary (stdout) --------------------------------------------
    print()
    print("GEM tracking-efficiency audit (leave-one-out)")
    print("=============================================")
    print()
    print(f"physics events processed : {n_phys}")
    print(f"events passing filter    : {n_events_filter_pass}  "
          f"(≥1 HyCal cluster AND ≥3 GEM detectors with hits)")
    print(f"HyCal clusters considered: {n_used}  "
          f"(after E ≥ {min_cluster_energy:.0f} MeV cut)")
    print(f"match_nsigma             : {params.match_nsigma}")
    print(f"max_chi2_per_dof         : {params.max_chi2}")
    print(f"max_hits_per_det         : {params.max_hits_per_det}")
    print(f"σ_GEM (mm)               : {params.gem_pos_res}")
    print()
    print("Each LOO variant runs once per (HyCal cluster, test detector D):")
    print("  - the OTHER 3 GEMs anchor a line fit (all 3 matched in seed")
    print("    window, χ²/dof ≤ max gate, per-det residual within "
          f"{params.match_nsigma}·σ_GEM),")
    print("  - the fit is projected to D and a hit is searched within "
          f"{params.match_nsigma}·σ_GEM[D].")
    print("Denominator = anchors that succeeded; numerator = anchors with a "
          "hit at D.")
    print()

    for mode in loo_modes:
        print(f"--- {mode.name} ---")
        for d in range(4):
            print(f"  GEM{d}: {_eff(mode.n_matched[d], mode.n_attempted[d])}")
        if mode.chi2_list:
            arr = sorted(mode.chi2_list)
            med  = arr[len(arr)//2]
            p90  = arr[int(0.9 * len(arr))]
            p99  = arr[int(0.99 * len(arr))]
            print(f"  anchor χ²/dof: median={med:.3f}  "
                  f"90%={p90:.3f}  99%={p99:.3f}  (n={len(arr)})")
        print()

    # ---- plots (matplotlib optional) --------------------------------------
    plt = _import_pyplot()
    if plt is None:
        print("[plot] matplotlib not available; skipping PNGs")
        return

    _plot_efficiency_bars(plt, loo_modes, out_dir / "efficiency.png")
    _plot_anchor_chi2(plt, loo_modes, params,
                      out_dir / "anchor_chi2.png")
    for mode in loo_modes:
        slug = mode.name.split()[0].replace("/", "_")
        _plot_residuals(plt, mode, out_dir / f"residuals_{slug}.png")


def _import_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None


def _plot_efficiency_bars(plt, modes: Sequence[LooStats], out: Path) -> None:
    """One panel per GEM (2x2 grid).  Each panel stacks one horizontal
    progress bar per LOO variant: full width = 100%, filled portion =
    detector efficiency, annotation = "XX.X%  (num/den)".  The bar IS the
    efficiency — denominator goes into the count text on the right."""
    n_dets = 4
    n_modes = len(modes)
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5))
    axes_flat = axes.flatten()
    bar_h = 0.65
    for d in range(n_dets):
        ax = axes_flat[d]
        # Plot bottom-up so loo (mode 0) ends up on top
        ys = list(range(n_modes - 1, -1, -1))
        for i, m in enumerate(modes):
            n_match = m.n_matched[d]
            n_att   = m.n_attempted[d]
            eff     = (100.0 * n_match / n_att) if n_att > 0 else 0.0
            y       = ys[i]
            # background "track" (full 100%)
            ax.barh(y, 100.0, height=bar_h, color="lightgray", alpha=0.4,
                    edgecolor="gray", linewidth=0.5, zorder=1)
            # filled portion (efficiency)
            ax.barh(y, eff, height=bar_h, color=f"C{i}", alpha=0.85,
                    zorder=2)
            # variant label on the left
            ax.text(-2, y, m.name.split()[0], va="center", ha="right",
                    fontsize=9, color=f"C{i}", fontweight="bold")
            # efficiency + counts on the right
            ann = (f"{eff:5.1f}%  ({n_match}/{n_att})"
                   if n_att > 0 else "no anchors")
            ax.text(102, y, ann, va="center", ha="left", fontsize=9)
        ax.set_xlim(0, 140)
        ax.set_ylim(-0.6, n_modes - 0.4)
        ax.set_yticks([])
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("efficiency (%)" if d >= 2 else "")
        ax.set_title(f"GEM{d}", fontsize=11, fontweight="bold")
        ax.grid(True, axis="x", alpha=0.3, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("LOO per-detector efficiency", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _plot_anchor_chi2(plt, modes: Sequence[LooStats],
                      params: TrackingParams, out: Path) -> None:
    """Two-row anchor-quality plot for the LOO line fits (HyCal + 3 OTHER
    GEMs).  Top row: per-variant χ²/dof histograms.  Bottom row: cumulative
    fraction of anchors retained as the χ²/dof cut sweeps.  Both rows show
    *anchor* statistics — they describe the quality of the 3-GEM-anchor
    sample, NOT detector efficiency (which is the GEM0..GEM3 numbers in the
    text summary)."""
    if not any(m.chi2_list for m in modes):
        return
    import numpy as np
    fig, (ax_hist, ax_acc) = plt.subplots(2, 1, figsize=(9, 8),
                                          sharex=True)

    # ---- top: χ²/dof histograms ------------------------------------------
    bins = np.linspace(0, max(params.max_chi2, 1.0) * 1.1, 60)
    for i, m in enumerate(modes):
        if not m.chi2_list:
            continue
        label = f"{m.name.split()[0]}  (n={len(m.chi2_list)})"
        ax_hist.hist(m.chi2_list, bins=bins, color=f"C{i}", alpha=0.55,
                     label=label)
    ax_hist.axvline(params.max_chi2, color="k", ls="--", lw=1,
                    label=f"current cut = {params.max_chi2}")
    ax_hist.set_ylabel("LOO anchors accepted")
    ax_hist.set_title("LOO anchor χ²/dof distribution "
                      "(HyCal + 3 OTHER GEMs)")
    ax_hist.grid(True, alpha=0.3)
    ax_hist.legend(loc="best", fontsize=9)

    # ---- bottom: cumulative anchor acceptance ----------------------------
    cuts = np.linspace(0.5, max(params.max_chi2, 1.0) * 1.5, 60)
    for i, m in enumerate(modes):
        if not m.chi2_list:
            continue
        arr = np.asarray(m.chi2_list)
        accepted_at_cut = np.array(
            [(arr <= c).sum() for c in cuts], dtype=float)
        total = float(len(arr))
        ax_acc.plot(cuts, 100.0 * accepted_at_cut / total, color=f"C{i}",
                    lw=1.6, label=m.name.split()[0])
    ax_acc.axvline(params.max_chi2, color="k", ls="--", lw=1,
                   label=f"current cut = {params.max_chi2}")
    ax_acc.set_xlabel("χ²/dof cut")
    ax_acc.set_ylabel("anchors retained (%)  [NOT detector efficiency]")
    ax_acc.set_title("Cumulative anchor acceptance vs χ²/dof gate")
    ax_acc.grid(True, alpha=0.3)
    ax_acc.legend(loc="best", fontsize=9)

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
    fig.suptitle(f"LOO residuals: hit − projected at test detector  "
                 f"[{loo.name.split()[0]}]")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


if __name__ == "__main__":
    sys.exit(main())
