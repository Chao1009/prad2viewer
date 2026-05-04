// TrackMatcher.cpp — single-source-of-truth straight-line track matcher
//=============================================================================
// Algorithm (per call):
//   1. Build seed line(s) according to `seed_mode`.
//   2. For each seed: project to every fit-eligible plane, take the
//      closest hit within match_nsigma · σ_total of the projection.
//      σ_total = √(σ_HC@plane² + σ_plane²),
//      σ_HC@plane = σ_HC · |(z_plane - z_target) / (z_HC - z_target)|.
//   3. If ≥ require_planes matched, build a weighted-LSQ fit through
//      HyCal + matched planes (+ target if requested) and gate on
//      χ²/dof ≤ max_chi2.
//   4. Per-plane post-fit residual within match_nsigma · σ_plane.
//   5. Keep the lowest-χ² track that passes; return nullopt otherwise.
//
// Frame: lab.  All inputs carry their lab-frame (x, y, z); the matcher
// only consults Plane::xform when a caller asks for LOO predictions in
// the test plane's local frame (handled via projectLineToLocal).
//=============================================================================

#include "TrackMatcher.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace prad2::trk {

namespace {

// Effective scalar σ for an isotropic-window match at a plane that has
// (potentially) anisotropic σ_x / σ_y.  Matches what runGemEfficiency
// historically did (σ_x = σ_y = sigma_gem[d]); for the anisotropic case
// we use the geometric mean so a circular cut is close to the equivalent
// area of the σ_x × σ_y ellipse.
inline float planeSigma(const Plane &p)
{
    if (p.sigma_x == p.sigma_y) return p.sigma_x;
    return std::sqrt(p.sigma_x * p.sigma_y);
}

inline bool bit(uint32_t mask, int d) { return (mask >> d) & 1u; }

inline int popcount32(uint32_t v)
{
    int c = 0;
    for (; v; v >>= 1) c += (int)(v & 1u);
    return c;
}

// HyCal σ at a given z, scaled along the line target→HyCal.
//   σ_HC@z = σ_HC · |(z - z_target) / (z_HC - z_target)|
// Falls back to σ_HC when the lever arm is degenerate.
inline float hcSigmaAt(float sigma_hc, float z_plane, float z_hc, float z_tgt)
{
    float lever = z_hc - z_tgt;
    if (std::abs(lever) < 1e-6f) return sigma_hc;
    return sigma_hc * std::abs((z_plane - z_tgt) / lever);
}

// Distance² between (x, y) and (px, py).
inline float dist2(float x, float y, float px, float py)
{
    float dx = x - px, dy = y - py;
    return dx*dx + dy*dy;
}

}  // namespace

TrackMatcher::TrackMatcher(MatcherConfig cfg) : cfg_(std::move(cfg)) {}

float TrackMatcher::hcSigma(float energy_MeV) const
{
    float E = energy_MeV * 0.001f;        // → GeV
    if (E <= 0.f) return cfg_.hc_sigma_C;
    float a = cfg_.hc_sigma_A / std::sqrt(E);
    float b = cfg_.hc_sigma_B / E;
    float c = cfg_.hc_sigma_C;
    return std::sqrt(a*a + b*b + c*c);
}

// ---------------------------------------------------------------------------
// Internal anchor-builder.  Given a seed line, for each fit-eligible
// plane (cand_mask) find the closest hit within match_nsigma · σ_total.
// If ≥ require_planes match, build the fit and gate on χ²/dof + per-
// plane residual.  Returns the populated Track or nullopt.
// ---------------------------------------------------------------------------
namespace {

struct AnchorScratch {
    bool                            matched[kMaxPlanes] = {};
    PlaneHit                        cand[kMaxPlanes]    = {};
    int                             n_matched = 0;
};

bool runFromSeed(const MatcherConfig &cfg,
                 const Line3D &seed,
                 const ClusterHit &hc,
                 float sigma_hc,
                 const std::vector<std::vector<PlaneHit>> &hits_by_plane,
                 uint32_t fit_mask,
                 bool target_in_fit,
                 int require_planes,
                 int seed_plane, int seed_hit_idx,   // for output bookkeeping; seed plane is auto-matched if ≥0
                 Stats *diag,
                 Track &out)
{
    const int n_planes = (int)cfg.planes.size();
    AnchorScratch sc{};

    // Seed plane (if any) is automatically a matched candidate.
    if (seed_plane >= 0 && seed_plane < n_planes && seed_hit_idx >= 0) {
        const auto &hits = hits_by_plane[seed_plane];
        if (seed_hit_idx < (int)hits.size()) {
            sc.matched[seed_plane] = true;
            sc.cand[seed_plane]    = hits[seed_hit_idx];
            ++sc.n_matched;
        }
    }

    // Project seed to every fit-eligible plane (skipping the seed plane,
    // which is already matched), pick closest within window.  Matching is
    // done in the plane's local frame because σ_plane is defined locally
    // (strip pitch / resolution are local quantities) — for a tilted plane
    // the lab-frame distance differs from the local one and would drift
    // the cut.
    for (int d = 0; d < n_planes; ++d) {
        if (!bit(fit_mask, d)) continue;
        if (sc.matched[d])     continue;
        const auto &plane = cfg.planes[d];
        const auto &hits  = hits_by_plane[d];
        if (hits.empty())      continue;

        float pred_lx, pred_ly;
        projectLineToLocal(plane.xform, seed, pred_lx, pred_ly);

        const float s_hc_at = hcSigmaAt(sigma_hc, plane.z, hc.z, cfg.target_z);
        const float s_pl    = planeSigma(plane);
        const float s_total = std::sqrt(s_hc_at*s_hc_at + s_pl*s_pl);
        const float cut     = cfg.match_nsigma * s_total;
        const float cut2    = cut * cut;

        const int max_n = std::min((int)hits.size(), cfg.max_hits_per_plane);
        int   best_i = -1;
        float best_d2 = cut2;
        for (int i = 0; i < max_n; ++i) {
            float lx, ly, lz;
            plane.xform.labToLocal(hits[i].x, hits[i].y, hits[i].z, lx, ly, lz);
            float d2 = dist2(lx, ly, pred_lx, pred_ly);
            if (d2 < best_d2) { best_d2 = d2; best_i = i; }
        }
        if (best_i >= 0) {
            sc.matched[d] = true;
            sc.cand[d]    = hits[best_i];
            ++sc.n_matched;
        }
    }

    if (sc.n_matched < require_planes) return false;
    if (diag) ++diag->n_min_match;

    // Build fit arrays.
    constexpr int CAP = kMaxPlanes + 2;   // + HyCal + target
    float zarr[CAP], xarr[CAP], yarr[CAP], wxarr[CAP], wyarr[CAP];
    int   N = 0;
    const float w_h = 1.f / (sigma_hc * sigma_hc);
    zarr[N] = hc.z;  xarr[N] = hc.x;  yarr[N] = hc.y;
    wxarr[N] = wyarr[N] = w_h;
    ++N;

    if (target_in_fit) {
        const float lever = (hc.z != cfg.target_z) ? (hc.z - cfg.target_z) : 1.f;
        const float bx_est = (hc.x - cfg.target_x) / lever;
        const float by_est = (hc.y - cfg.target_y) / lever;
        const float sx2 = cfg.target_sigma_x * cfg.target_sigma_x
                          + (bx_est * cfg.target_sigma_z)
                          * (bx_est * cfg.target_sigma_z);
        const float sy2 = cfg.target_sigma_y * cfg.target_sigma_y
                          + (by_est * cfg.target_sigma_z)
                          * (by_est * cfg.target_sigma_z);
        zarr[N]  = cfg.target_z;
        xarr[N]  = cfg.target_x;
        yarr[N]  = cfg.target_y;
        wxarr[N] = 1.f / sx2;
        wyarr[N] = 1.f / sy2;
        ++N;
    }

    for (int d = 0; d < n_planes; ++d) {
        if (!sc.matched[d]) continue;
        const auto &plane = cfg.planes[d];
        zarr[N]  = sc.cand[d].z;
        xarr[N]  = sc.cand[d].x;
        yarr[N]  = sc.cand[d].y;
        wxarr[N] = 1.f / (plane.sigma_x * plane.sigma_x);
        wyarr[N] = 1.f / (plane.sigma_y * plane.sigma_y);
        ++N;
    }

    Line3D fit;
    if (!fitWeightedLine(N, zarr, xarr, yarr, wxarr, wyarr, fit)) return false;
    if (fit.chi2_per_dof > cfg.max_chi2) return false;
    if (diag) ++diag->n_pass_chi2;

    // Per-plane post-fit residual gate, again in plane-local frame so the
    // σ_plane interpretation stays consistent with the seed-window cut.
    for (int d = 0; d < n_planes; ++d) {
        if (!sc.matched[d]) continue;
        const auto &plane = cfg.planes[d];
        float plx, ply;
        projectLineToLocal(plane.xform, fit, plx, ply);
        float hlx, hly, hlz;
        plane.xform.labToLocal(sc.cand[d].x, sc.cand[d].y, sc.cand[d].z,
                               hlx, hly, hlz);
        const float s = planeSigma(plane);
        const float c = cfg.match_nsigma * s;
        if (dist2(hlx, hly, plx, ply) > c * c) return false;
    }
    if (diag) ++diag->n_pass_resid;

    // Populate Track.
    out = Track{};
    out.fit          = fit;
    out.seed_plane   = seed_plane;
    out.seed_hit_idx = seed_hit_idx;
    out.target_in_fit = target_in_fit;
    for (int d = 0; d < n_planes; ++d) {
        out.matched[d] = sc.matched[d];
        if (sc.matched[d]) {
            out.hit[d] = sc.cand[d];
            ++out.n_matched;
        }
    }
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public: findBestTrack (full-control overload).
// ---------------------------------------------------------------------------

std::optional<Track>
TrackMatcher::findBestTrack(const ClusterHit &hc,
                            const std::vector<std::vector<PlaneHit>> &hits_by_plane,
                            Seed seed_mode,
                            uint32_t seed_planes_mask,
                            uint32_t fit_planes_mask,
                            bool require_target_in_fit,
                            int  require_planes,
                            Stats *diag) const
{
    if (diag) ++diag->n_call;
    const int n_planes = (int)cfg_.planes.size();
    if (n_planes <= 0 || (int)hits_by_plane.size() < n_planes) return std::nullopt;

    const float sigma_hc = hcSigma(hc.energy);

    Track  best{};
    bool   have_best = false;
    float  best_chi2 = std::numeric_limits<float>::infinity();

    auto try_seed = [&](const Line3D &seed, int seed_plane, int seed_hit_idx) {
        if (diag) ++diag->n_seed_tried;
        Track cand;
        if (!runFromSeed(cfg_, seed, hc, sigma_hc, hits_by_plane,
                         fit_planes_mask, require_target_in_fit,
                         require_planes, seed_plane, seed_hit_idx,
                         diag, cand))
            return;
        if (cand.fit.chi2_per_dof < best_chi2) {
            best_chi2 = cand.fit.chi2_per_dof;
            best      = std::move(cand);
            have_best = true;
        }
    };

    switch (seed_mode) {
    case Seed::TargetToHC: {
        Line3D seed = seedLine(cfg_.target_x, cfg_.target_y, cfg_.target_z,
                               hc.x, hc.y, hc.z);
        try_seed(seed, /*seed_plane=*/-1, /*seed_hit_idx=*/-1);
        break;
    }
    case Seed::HCAndPlaneHit: {
        for (int s = 0; s < n_planes; ++s) {
            if (!bit(seed_planes_mask, s))    continue;
            if (!bit(fit_planes_mask,   s))   continue;  // seed plane must be in fit
            const auto &hits = hits_by_plane[s];
            const int max_n  = std::min((int)hits.size(), cfg_.max_hits_per_plane);
            for (int i = 0; i < max_n; ++i) {
                const auto &g = hits[i];
                Line3D seed = seedLine(hc.x, hc.y, hc.z, g.x, g.y, g.z);
                try_seed(seed, s, i);
            }
        }
        break;
    }
    case Seed::FreeCombinatorial: {
        // No seed line; enumerate every fit-eligible-plane subset of size
        // ≥ require_planes, then every per-plane hit cross-product.  The
        // "fit" is the seed.  Bounded by max_hits_per_plane^k * 2^n_planes
        // — small in practice (≤ 4 planes, ≤ 50 hits).
        std::vector<int> fit_planes;
        fit_planes.reserve(n_planes);
        for (int d = 0; d < n_planes; ++d)
            if (bit(fit_planes_mask, d) && !hits_by_plane[d].empty())
                fit_planes.push_back(d);
        const int K = (int)fit_planes.size();
        if (K < require_planes) break;

        // Iterate subsets of fit_planes via bitmask of length K.
        for (uint32_t sub = 1; sub < (1u << K); ++sub) {
            if (popcount32(sub) < require_planes) continue;

            // Subset-of-planes → mask in plane-index space.
            uint32_t sub_mask = 0;
            for (int b = 0; b < K; ++b)
                if (sub & (1u << b)) sub_mask |= (1u << fit_planes[b]);

            // Enumerate hit cross-product across the subset.  Use indices
            // into capped per-plane lists.
            int idx_max[kMaxPlanes] = {};
            int planes_in_subset[kMaxPlanes] = {};
            int K_sub = 0;
            for (int b = 0; b < K; ++b) {
                if (sub & (1u << b)) {
                    int d = fit_planes[b];
                    planes_in_subset[K_sub] = d;
                    idx_max[K_sub] = std::min((int)hits_by_plane[d].size(),
                                              cfg_.max_hits_per_plane);
                    ++K_sub;
                }
            }

            int idx[kMaxPlanes] = {};
            while (true) {
                // Build a fit through HC + selected hits.  No projection
                // step — the fit IS the line, χ² gates it.
                constexpr int CAP = kMaxPlanes + 1;
                float zarr[CAP], xarr[CAP], yarr[CAP], warr[CAP];
                int N = 0;
                const float w_h = 1.f / (sigma_hc * sigma_hc);
                zarr[N] = hc.z; xarr[N] = hc.x; yarr[N] = hc.y;
                warr[N] = w_h;
                ++N;
                for (int b = 0; b < K_sub; ++b) {
                    const int d = planes_in_subset[b];
                    const auto &h = hits_by_plane[d][idx[b]];
                    zarr[N] = h.z; xarr[N] = h.x; yarr[N] = h.y;
                    warr[N] = 1.f / (planeSigma(cfg_.planes[d])
                                     * planeSigma(cfg_.planes[d]));
                    ++N;
                }
                Line3D fit;
                if (diag) ++diag->n_seed_tried;
                if (fitWeightedLine(N, zarr, xarr, yarr, warr, fit)
                    && fit.chi2_per_dof <= cfg_.max_chi2)
                {
                    // Wrap as a Track via the same residual gate path —
                    // call runFromSeed with the fit as the "seed" line and
                    // the subset as fit_mask.  The matched-set will be
                    // identical because the fit goes through these hits.
                    if (cfg_.match_nsigma > 0) {
                        Track cand;
                        if (runFromSeed(cfg_, fit, hc, sigma_hc, hits_by_plane,
                                        sub_mask, require_target_in_fit,
                                        require_planes, /*seed_plane=*/-1,
                                        /*seed_hit_idx=*/-1, diag, cand)
                            && cand.fit.chi2_per_dof < best_chi2)
                        {
                            best_chi2 = cand.fit.chi2_per_dof;
                            best      = std::move(cand);
                            have_best = true;
                        }
                    }
                }

                // Increment the K_sub-digit index counter.
                int b = 0;
                while (b < K_sub) {
                    if (++idx[b] < idx_max[b]) break;
                    idx[b++] = 0;
                }
                if (b == K_sub) break;   // odometer wrapped
            }
        }
        break;
    }
    }

    if (!have_best) return std::nullopt;
    return best;
}

// ---------------------------------------------------------------------------
// Public: findBestTrack (convenience overload — physics default).
// ---------------------------------------------------------------------------

std::optional<Track>
TrackMatcher::findBestTrack(const ClusterHit &hc,
                            const std::vector<std::vector<PlaneHit>> &hits_by_plane,
                            bool free_seed,
                            int  require_planes,
                            Stats *diag) const
{
    const int n_planes = (int)cfg_.planes.size();
    uint32_t all_mask = 0;
    for (int d = 0; d < n_planes; ++d) all_mask |= (1u << d);
    return findBestTrack(hc, hits_by_plane,
                         free_seed ? Seed::FreeCombinatorial : Seed::TargetToHC,
                         /*seed_planes_mask=*/all_mask,
                         /*fit_planes_mask=*/all_mask,
                         /*require_target_in_fit=*/false,
                         require_planes,
                         diag);
}

// ---------------------------------------------------------------------------
// Public: runLoo.
// ---------------------------------------------------------------------------

std::optional<TrackMatcher::LooResult>
TrackMatcher::runLoo(int test_plane,
                     const ClusterHit &hc,
                     const std::vector<std::vector<PlaneHit>> &hits_by_plane,
                     Seed seed_mode,
                     bool target_in_fit,
                     Stats *diag) const
{
    const int n_planes = (int)cfg_.planes.size();
    if (test_plane < 0 || test_plane >= n_planes) return std::nullopt;

    uint32_t all_mask = 0;
    for (int d = 0; d < n_planes; ++d) all_mask |= (1u << d);
    const uint32_t fit_mask  = all_mask & ~(1u << test_plane);
    const uint32_t seed_mask = fit_mask;     // seed never uses the test plane

    // LOO requires every OTHER plane to have a matched hit — that's the
    // "clean basis" rule the audit script and AppState::runGemEfficiency
    // both use (n_planes - 1 for one held-out plane).
    auto anchor = findBestTrack(hc, hits_by_plane, seed_mode,
                                seed_mask, fit_mask, target_in_fit,
                                /*require_planes=*/n_planes - 1,
                                diag);
    if (!anchor) return std::nullopt;

    // Search for closest hit at the test plane within match_nsigma · σ_plane
    // (no σ_HC term — HC is already in the fit and the projection here is
    // dominated by σ_plane).  Match in local frame for the same reason as
    // anchor matching above.
    const auto &tp = cfg_.planes[test_plane];
    float pred_lx, pred_ly;
    projectLineToLocal(tp.xform, anchor->fit, pred_lx, pred_ly);

    const float s_pl = planeSigma(tp);
    const float cut  = cfg_.match_nsigma * s_pl;
    const float cut2 = cut * cut;

    const auto &hits = hits_by_plane[test_plane];
    const int max_n = std::min((int)hits.size(), cfg_.max_hits_per_plane);
    int   best_i  = -1;
    float best_d2 = cut2;
    for (int i = 0; i < max_n; ++i) {
        float lx, ly, lz;
        tp.xform.labToLocal(hits[i].x, hits[i].y, hits[i].z, lx, ly, lz);
        float d2 = dist2(lx, ly, pred_lx, pred_ly);
        if (d2 < best_d2) { best_d2 = d2; best_i = i; }
    }

    // Pred reported back to caller is naive lab-frame at z = z_plane (where
    // the line crosses the plane's nominal z), not the on-surface point —
    // callers that want the on-surface lab pred should re-project via
    // projectLineToLocal + xform.toLab themselves (see
    // AppState::runGemEfficiency snapshot population).
    LooResult r;
    r.anchor     = std::move(*anchor);
    r.test_plane = test_plane;
    r.pred_x     = r.anchor.fit.ax + r.anchor.fit.bx * tp.z;
    r.pred_y     = r.anchor.fit.ay + r.anchor.fit.by * tp.z;
    if (best_i >= 0) r.hit_at_test = hits[best_i];
    return r;
}

// ---------------------------------------------------------------------------
// Public: findPerPlaneMatches.
// ---------------------------------------------------------------------------

PerPlaneMatch
TrackMatcher::findPerPlaneMatches(
    const ClusterHit &hc,
    const std::vector<std::vector<PlaneHit>> &hits_by_plane) const
{
    PerPlaneMatch out{};
    out.seed = seedLine(cfg_.target_x, cfg_.target_y, cfg_.target_z,
                        hc.x, hc.y, hc.z);

    const float sigma_hc = hcSigma(hc.energy);
    const int n_planes = (int)cfg_.planes.size();
    if ((int)hits_by_plane.size() < n_planes) return out;

    for (int d = 0; d < n_planes; ++d) {
        const auto &plane = cfg_.planes[d];
        const auto &hits  = hits_by_plane[d];
        if (hits.empty()) continue;

        float pred_lx, pred_ly;
        projectLineToLocal(plane.xform, out.seed, pred_lx, pred_ly);

        const float s_hc_at = hcSigmaAt(sigma_hc, plane.z, hc.z, cfg_.target_z);
        const float s_pl    = planeSigma(plane);
        const float s_total = std::sqrt(s_hc_at*s_hc_at + s_pl*s_pl);
        const float cut     = cfg_.match_nsigma * s_total;
        const float cut2    = cut * cut;

        const int max_n = std::min((int)hits.size(), cfg_.max_hits_per_plane);
        int   best_i  = -1;
        float best_d2 = cut2;
        for (int i = 0; i < max_n; ++i) {
            float lx, ly, lz;
            plane.xform.labToLocal(hits[i].x, hits[i].y, hits[i].z, lx, ly, lz);
            float d2 = dist2(lx, ly, pred_lx, pred_ly);
            if (d2 < best_d2) { best_d2 = d2; best_i = i; }
        }
        if (best_i >= 0) {
            out.matched[d]  = true;
            out.hit[d]      = hits[best_i];
            out.residual[d] = std::sqrt(best_d2);
        }
    }
    return out;
}

} // namespace prad2::trk
