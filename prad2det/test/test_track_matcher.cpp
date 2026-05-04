// test_track_matcher.cpp — smoke tests for prad2::trk::TrackMatcher
//=============================================================================
// Self-contained: no external deps, returns 0 on pass.  Builds only when
// PRAD2DET_BUILD_TESTS=ON via the prad2det CMakeLists.txt.
//
// Coverage:
//   - exact recovery of a perfect 4-plane track (TargetToHC, HCAndPlaneHit,
//     FreeCombinatorial all find the same line within 1e-3)
//   - target-in-fit doesn't break a target-pointing track
//   - LOO at each plane finds the held-out hit
//   - LOO with a missing held-out hit reports inefficient
//   - σ filter rejects an outlier hit
//=============================================================================

#include "TrackMatcher.h"
#include "TrackGeometry.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace prad2::trk;

#define CHECK(cond) do {                                                      \
    if (!(cond)) {                                                            \
        std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);  \
        std::exit(1);                                                         \
    }                                                                         \
} while (0)

#define CHECK_NEAR(a, b, tol) do {                                            \
    float _a = (a), _b = (b);                                                 \
    if (std::fabs(_a - _b) > (tol)) {                                         \
        std::fprintf(stderr, "FAIL %s:%d: |%g - %g| > %g (%s vs %s)\n",       \
                     __FILE__, __LINE__, _a, _b, (float)(tol), #a, #b);      \
        std::exit(1);                                                         \
    }                                                                         \
} while (0)

namespace {

MatcherConfig make_config()
{
    MatcherConfig cfg;
    cfg.planes.resize(4);
    // GEM-ish layout: two pairs at z ≈ 1.5 m and z ≈ 5.2 m.
    const float zs[4] = {1500.f, 1505.f, 5200.f, 5205.f};
    for (int d = 0; d < 4; ++d) {
        cfg.planes[d].id      = d;
        cfg.planes[d].z       = zs[d];
        cfg.planes[d].sigma_x = 0.5f;
        cfg.planes[d].sigma_y = 0.5f;
        cfg.planes[d].xform.set(0.f, 0.f, zs[d], 0.f, 0.f, 0.f);
    }
    cfg.hc_sigma_A = 2.5f;
    cfg.hc_sigma_B = 0.f;
    cfg.hc_sigma_C = 0.f;
    cfg.target_x   = 0.f;
    cfg.target_y   = 0.f;
    cfg.target_z   = 0.f;
    cfg.target_sigma_x = 1.f;
    cfg.target_sigma_y = 1.f;
    cfg.target_sigma_z = 20.f;
    cfg.max_hits_per_plane = 50;
    cfg.match_nsigma       = 3.0f;
    cfg.max_chi2           = 5.0f;
    return cfg;
}

// Build a perfect track through (0,0,0) at angles (slope_x, slope_y).
ClusterHit make_hc(float slope_x, float slope_y, float z_hc, float energy_MeV)
{
    ClusterHit hc;
    hc.x = slope_x * z_hc;
    hc.y = slope_y * z_hc;
    hc.z = z_hc;
    hc.energy = energy_MeV;
    return hc;
}

std::vector<std::vector<PlaneHit>>
make_perfect_hits(const MatcherConfig &cfg, float slope_x, float slope_y)
{
    std::vector<std::vector<PlaneHit>> hits(cfg.planes.size());
    for (int d = 0; d < (int)cfg.planes.size(); ++d) {
        PlaneHit h;
        h.plane_idx = d;
        h.hit_idx   = 0;
        h.z = cfg.planes[d].z;
        h.x = slope_x * h.z;
        h.y = slope_y * h.z;
        hits[d].push_back(h);
    }
    return hits;
}

void test_perfect_track_target_seed()
{
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    const float z_hc = 5800.f;
    auto hc   = make_hc(sx, sy, z_hc, /*E=*/1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    Stats diag{};
    auto track = m.findBestTrack(hc, hits, /*free_seed=*/false,
                                 /*require_planes=*/3, &diag);
    CHECK(track.has_value());
    CHECK(track->n_matched == 4);
    for (int d = 0; d < 4; ++d) CHECK(track->matched[d]);
    CHECK_NEAR(track->fit.bx, sx, 1e-4);
    CHECK_NEAR(track->fit.by, sy, 1e-4);
    CHECK(track->fit.chi2_per_dof < 1e-3f);
    CHECK(track->seed_plane == -1);     // target-seeded
    CHECK(diag.n_call == 1);
    CHECK(diag.n_pass_resid == 1);
    std::printf("  perfect_track_target_seed OK (chi2/dof=%g)\n",
                track->fit.chi2_per_dof);
}

void test_perfect_track_hc_plane_seed()
{
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    Stats diag{};
    uint32_t all_mask = 0xF;       // 4 planes
    auto track = m.findBestTrack(hc, hits,
                                 Seed::HCAndPlaneHit,
                                 /*seed_planes_mask=*/all_mask,
                                 /*fit_planes_mask=*/all_mask,
                                 /*require_target_in_fit=*/false,
                                 /*require_planes=*/3,
                                 &diag);
    CHECK(track.has_value());
    CHECK(track->n_matched == 4);
    CHECK_NEAR(track->fit.bx, sx, 1e-4);
    CHECK_NEAR(track->fit.by, sy, 1e-4);
    CHECK(track->seed_plane >= 0 && track->seed_plane < 4);
    CHECK(diag.n_seed_tried == 4);   // one per plane
    std::printf("  perfect_track_hc_plane_seed OK (seed plane=%d)\n",
                track->seed_plane);
}

void test_perfect_track_free()
{
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    Stats diag{};
    auto track = m.findBestTrack(hc, hits, /*free_seed=*/true,
                                 /*require_planes=*/3, &diag);
    CHECK(track.has_value());
    CHECK(track->n_matched == 4);
    CHECK_NEAR(track->fit.bx, sx, 1e-4);
    CHECK_NEAR(track->fit.by, sy, 1e-4);
    std::printf("  perfect_track_free OK\n");
}

void test_target_in_fit()
{
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    uint32_t all_mask = 0xF;
    auto track = m.findBestTrack(hc, hits,
                                 Seed::HCAndPlaneHit, all_mask, all_mask,
                                 /*require_target_in_fit=*/true,
                                 /*require_planes=*/3, nullptr);
    CHECK(track.has_value());
    CHECK(track->target_in_fit);
    CHECK_NEAR(track->fit.bx, sx, 1e-4);
    CHECK_NEAR(track->fit.by, sy, 1e-4);
    std::printf("  target_in_fit OK\n");
}

void test_loo_each_plane()
{
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    for (int d = 0; d < 4; ++d) {
        auto loo = m.runLoo(d, hc, hits, Seed::HCAndPlaneHit,
                            /*target_in_fit=*/false, nullptr);
        CHECK(loo.has_value());
        CHECK(loo->test_plane == d);
        CHECK(loo->hit_at_test.has_value());
        CHECK_NEAR(loo->pred_x, sx * cfg.planes[d].z, 1e-3);
        CHECK_NEAR(loo->pred_y, sy * cfg.planes[d].z, 1e-3);
        // Anchor must NOT have used the test plane.
        CHECK(!loo->anchor.matched[d]);
    }
    std::printf("  loo_each_plane OK\n");
}

void test_loo_inefficient()
{
    // Held-out hit absent → LOO returns the anchor but no hit_at_test.
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    const int test_d = 2;
    hits[test_d].clear();   // simulate detector inefficiency

    auto loo = m.runLoo(test_d, hc, hits, Seed::HCAndPlaneHit,
                        /*target_in_fit=*/false, nullptr);
    CHECK(loo.has_value());
    CHECK(!loo->hit_at_test.has_value());
    CHECK(!loo->anchor.matched[test_d]);
    std::printf("  loo_inefficient OK\n");
}

void test_outlier_rejected()
{
    // Add a wildly off-axis hit on plane 0; matcher should still find the
    // straight-line track because the σ window kills the outlier.
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    PlaneHit bad;
    bad.plane_idx = 0; bad.hit_idx = 1;
    bad.z = cfg.planes[0].z;
    bad.x = sx * bad.z + 200.f;     // 200 mm off, way past 3σ
    bad.y = sy * bad.z;
    hits[0].push_back(bad);

    auto track = m.findBestTrack(hc, hits, /*free_seed=*/false,
                                 /*require_planes=*/3, nullptr);
    CHECK(track.has_value());
    CHECK(track->matched[0]);
    CHECK(track->hit[0].hit_idx == 0);   // good hit, not the outlier
    std::printf("  outlier_rejected OK\n");
}

void test_find_per_plane_matches()
{
    // Per-plane independent matching — should populate matched[d] / hit[d]
    // / residual[d] for every plane that has a hit in its seed window, with
    // NO χ²/dof gate or post-fit residual gate.
    auto cfg = make_config();
    TrackMatcher m(cfg);
    const float sx = 0.02f, sy = -0.015f;
    auto hc   = make_hc(sx, sy, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, sx, sy);

    auto pp = m.findPerPlaneMatches(hc, hits);
    for (int d = 0; d < 4; ++d) {
        CHECK(pp.matched[d]);
        CHECK(pp.hit[d].plane_idx == d);
        CHECK_NEAR(pp.residual[d], 0.f, 1e-3);
    }
    std::printf("  find_per_plane_matches OK\n");

    // Drop plane 1's hits — only planes 0/2/3 should match, no fit involved.
    hits[1].clear();
    pp = m.findPerPlaneMatches(hc, hits);
    CHECK(pp.matched[0]);
    CHECK(!pp.matched[1]);
    CHECK(pp.matched[2]);
    CHECK(pp.matched[3]);
    std::printf("  find_per_plane_matches_partial OK\n");
}

void test_low_multiplicity_no_track()
{
    // Only 2 planes have hits — below require_planes=3 → no track.
    auto cfg = make_config();
    TrackMatcher m(cfg);
    auto hc   = make_hc(0.02f, -0.015f, 5800.f, 1500.f);
    auto hits = make_perfect_hits(cfg, 0.02f, -0.015f);
    hits[2].clear();
    hits[3].clear();

    auto track = m.findBestTrack(hc, hits, /*free_seed=*/false,
                                 /*require_planes=*/3, nullptr);
    CHECK(!track.has_value());
    std::printf("  low_multiplicity_no_track OK\n");
}

}  // namespace

int main()
{
    std::printf("[prad2det] TrackMatcher smoke tests\n");
    test_perfect_track_target_seed();
    test_perfect_track_hc_plane_seed();
    test_perfect_track_free();
    test_target_in_fit();
    test_loo_each_plane();
    test_loo_inefficient();
    test_outlier_rejected();
    test_find_per_plane_matches();
    test_low_multiplicity_no_track();
    std::printf("[prad2det] all TrackMatcher tests passed\n");
    return 0;
}
