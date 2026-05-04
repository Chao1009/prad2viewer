#pragma once
//=============================================================================
// TrackMatcher.h — single-source-of-truth straight-line track matcher
//
// Replaces the four scattered implementations of "find a line through one
// HyCal cluster and a bag of GEM hits" that existed before:
//
//   - analysis/{include,src}/MatchingTools.{h,cpp}  (square/circular cut, no fit)
//   - src/app_state.cpp::runGemEfficiency           (full LOO + χ² + funnel)
//   - analysis/pyscripts/gem_eff_audit.py           (Python LOO re-impl)
//   - analysis/pyscripts/gem_hycal_matching.py      (simple σ-cut, target seed)
//
// Use cases:
//   1. Physics: best track per HyCal cluster (target seed by default;
//      free-combinatorial fallback for halo / non-target tracks).
//   2. Efficiency: leave-one-out — fit excludes a test plane, project to
//      it, count whether a hit falls in the per-plane σ window.
//
// Single-HC matching only.  Multi-cluster contention (greedy-by-E or
// Hungarian) lives at the caller — the matcher is stateless and pure
// per call so callers can compose policies on top of it.
//
// Frame contract: all positions in the lab frame (target at the origin
// in the canonical PRad-II setup, but `target_x/y/z` in MatcherConfig
// support arbitrary placement).  Per-plane σ_x, σ_y are at the plane's
// own z; σ_HC is evaluated at the HyCal face.
//=============================================================================

#include "DetectorTransform.h"
#include "TrackGeometry.h"

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace prad2::trk {

inline constexpr int kMaxPlanes = 8;

// One detector layer.  `xform` is only consulted when the matcher needs
// to project a line into the plane's local frame (LOO inside-fiducial
// checks, residual snapshots).  Lab-frame-only callers can leave it
// default-constructed.
struct Plane {
    int    id      = -1;     // arbitrary tag (e.g. 0..3 for GEM1..GEM4)
    float  z       = 0.f;    // lab z (mm)
    float  sigma_x = 0.1f;   // mm — position resolution on this plane
    float  sigma_y = 0.1f;   // mm — anisotropic σ supported from day one
    DetectorTransform xform = {};
};

// One reconstructed hit on a plane, in lab coords.  `hit_idx` is the
// caller's index into its own per-plane hit list — auxiliary payload
// (charge, size, timing…) stays caller-side, looked up via (plane_idx,
// hit_idx).  The matcher never inspects payload.
struct PlaneHit {
    int    plane_idx = -1;
    int    hit_idx   = -1;
    float  x = 0.f, y = 0.f, z = 0.f;
};

// One HyCal cluster — anchor point of every track.
struct ClusterHit {
    float    x = 0.f, y = 0.f, z = 0.f;
    float    energy = 0.f;          // MeV
    uint16_t center_id = 0;
    uint32_t flag = 0;
};

struct MatcherConfig {
    std::vector<Plane> planes;

    // HyCal position resolution at the HyCal face:
    //   σ_HC(E[GeV]) = √((A/√E)² + (B/E)² + C²)   [mm]
    // The matcher scales this to any z via the seed-line lever arm.
    float  hc_sigma_A = 2.5f;
    float  hc_sigma_B = 0.f;
    float  hc_sigma_C = 0.f;

    // Target geometry + transverse beam-spot σ.  σ_target_z is the
    // longitudinal target half-extent; couples to σ_x_eff / σ_y_eff via
    // the slope at z = target_z when the target is in the fit.
    float  target_x = 0.f, target_y = 0.f, target_z = 0.f;
    float  target_sigma_x = 1.f;
    float  target_sigma_y = 1.f;
    float  target_sigma_z = 20.f;

    int    max_hits_per_plane = 50;   // cap on seed enumeration / candidate scan
    float  match_nsigma       = 3.0f; // window for both seed projection and post-fit residual
    float  max_chi2           = 3.5f; // χ²/dof gate for accepting a track
};

enum class Seed : uint8_t {
    TargetToHC,         // line goes (target_x,target_y,target_z) → HyCal cluster
    HCAndPlaneHit,      // try every (HC, plane-d hit) pair on planes set in seed_planes_mask
    FreeCombinatorial,  // no seed: enumerate every fit-plane subset × hit cross-product
};

// Per-call funnel counters — useful for tuning σ / nsigma / max_chi2
// without recompiling.  Owned + reset by the caller; the matcher only
// increments.
struct Stats {
    int n_call         = 0;   // entries to findBestTrack
    int n_seed_tried   = 0;   // seed lines built (HCAndPlaneHit) or subsets enumerated (Free)
    int n_min_match    = 0;   // ≥require_planes matched in the seed window
    int n_pass_chi2    = 0;   // χ²/dof ≤ max_chi2
    int n_pass_resid   = 0;   // every per-plane residual within match_nsigma · σ_GEM
};

// Result of a successful match.  Index into `matched` / `hit` is the
// plane index in MatcherConfig::planes (= outer index of hits_by_plane).
struct Track {
    Line3D                            fit;
    std::array<bool, kMaxPlanes>      matched{};
    std::array<PlaneHit, kMaxPlanes>  hit{};
    int   n_matched     = 0;
    int   seed_plane    = -1;   // -1 when target-seeded or free
    int   seed_hit_idx  = -1;
    bool  target_in_fit = false;
};

// Per-plane independent matching — for each plane, the closest hit (in
// plane-local coords) within match_nsigma · σ_total of the target→HC
// seed-line projection.  No fit, no χ²/dof, no per-plane residual gate;
// each plane is an independent decision.  Used by replay tools and the
// gem_hycal_matching.py per-cluster-per-plane TSV writer.  The legacy
// MatchingTools::MatchPerChamber covered the same intent with a fixed
// 15 mm square cut; this is the σ-scaled successor.
struct PerPlaneMatch {
    Line3D                            seed;        // target → HC, lab frame
    std::array<bool, kMaxPlanes>      matched{};
    std::array<PlaneHit, kMaxPlanes>  hit{};
    std::array<float, kMaxPlanes>     residual{};  // local-frame distance
};

class TrackMatcher {
public:
    explicit TrackMatcher(MatcherConfig cfg);

    const MatcherConfig &config() const { return cfg_; }

    // Core entry point.  Iterates over the seeds implied by `seed_mode`
    // and `seed_planes_mask`, builds a fit through HC + matched planes
    // (+ optional target), gates on χ²/dof and per-plane residual, and
    // returns the lowest-χ² track that passes.
    //
    //   seed_planes_mask : bit d ⇒ plane d may donate seed hits (Seed::HCAndPlaneHit)
    //   fit_planes_mask  : bit d ⇒ plane d eligible to enter the fit
    //   require_target_in_fit : add target as a soft fit point (anisotropic σ)
    //   require_planes  : minimum # of planes (excluding HC/target) in the fit
    //
    // Returns nullopt on no acceptable track.  `diag`, when non-null,
    // accumulates funnel counters for this call (caller owns reset).
    std::optional<Track>
    findBestTrack(const ClusterHit &hc,
                  const std::vector<std::vector<PlaneHit>> &hits_by_plane,
                  Seed seed_mode,
                  uint32_t seed_planes_mask,
                  uint32_t fit_planes_mask,
                  bool require_target_in_fit,
                  int  require_planes,
                  Stats *diag = nullptr) const;

    // Convenience: physics "best track or none".
    //   free_seed = false → Seed::TargetToHC, all planes eligible for fit, no target in fit
    //   free_seed = true  → Seed::FreeCombinatorial, all planes eligible, no target in fit
    // require_planes defaults to 3 (matches current LOO acceptance).
    std::optional<Track>
    findBestTrack(const ClusterHit &hc,
                  const std::vector<std::vector<PlaneHit>> &hits_by_plane,
                  bool free_seed = false,
                  int  require_planes = 3,
                  Stats *diag = nullptr) const;

    // LOO: fit excludes test_plane; project the resulting line to
    // test_plane, search for the closest hit within match_nsigma · σ_y
    // (≈σ_x; matcher uses √(σ_x²+σ_y²)/√2 for the effective scalar
    // window).  `seed_mode` is HCAndPlaneHit or TargetToHC; the seed
    // never uses test_plane.  `target_in_fit` controls whether the
    // target enters the anchor fit.
    struct LooResult {
        Track    anchor;
        int      test_plane = -1;
        float    pred_x = 0.f, pred_y = 0.f;     // lab @ test plane z
        std::optional<PlaneHit> hit_at_test;     // nullopt ⇒ inefficient
    };
    std::optional<LooResult>
    runLoo(int test_plane,
           const ClusterHit &hc,
           const std::vector<std::vector<PlaneHit>> &hits_by_plane,
           Seed seed_mode,
           bool target_in_fit,
           Stats *diag = nullptr) const;

    // Per-plane independent best-match against the target → HC seed line.
    // No fit / χ² / residual gate — purely the seed-window cut at each
    // plane.  Use this for tools that want one row per (HC, plane) pair
    // (replay's per-chamber output, gem_hycal_matching.py TSV) without
    // requiring the matched hits to form a coherent track.
    PerPlaneMatch
    findPerPlaneMatches(const ClusterHit &hc,
                        const std::vector<std::vector<PlaneHit>> &hits_by_plane) const;

    // σ_HC at HyCal face, evaluated from cfg_.hc_sigma_{A,B,C}.
    // energy is in MeV (matches HyCalCluster output).
    float hcSigma(float energy_MeV) const;

private:
    MatcherConfig cfg_;
};

} // namespace prad2::trk
