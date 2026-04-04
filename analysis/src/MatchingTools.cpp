// MatchingTools.cpp — tools for matching HyCal clusters to GEM hits
//=============================================================================
// Adapted from PRadAnalyzer/PRadDetMatch.cpp.
// Usage:
//   MatchingTools matcher;
//   auto matches = matcher.Match(hycalHits, gem1Hits, gem2Hits, gem3Hits, gem4Hits);
//   for (auto &m : matches) { 
//       m.hycal_hit is the HyCal cluster
//       m.gem is the best-matched GEM hit (if any)
//       m.gem1_hits, m.gem2_hits, m.gem3_hits, m.gem4_hits are the candidates in each plane (sorted by distance)
//       m.mflag has bits set for each plane with a match
//       m.hycal_idx is the index of the cluster in the original vector
//   }
//=============================================================================


#include "MatchingTools.h"
#include <algorithm>
#include <set>

namespace analysis {

// ============================================================================
// Projection
// ============================================================================

ProjectHit GetProjectionHits(float x, float y, float z,
                                            float projection_z)
{   
    // simple linear projection from (x,y,z) to (x_proj, y_proj, projection_z)
    // in target and beam center coordinates
    float scale = projection_z / z;
    return ProjectHit(x * scale, y * scale, projection_z);
}

void GetProjection(HCHit &hc, float projection_z)
{
    ProjectHit proj = GetProjectionHits(hc.x, hc.y, hc.z, projection_z);
    hc.x = proj.x_proj;
    hc.y = proj.y_proj;
    hc.z = proj.z_proj;
}

void GetProjection(GEMHit &gem, float projection_z)
{
    ProjectHit proj = GetProjectionHits(gem.x, gem.y, gem.z, projection_z);
    gem.x = proj.x_proj;
    gem.y = proj.y_proj;
    gem.z = proj.z_proj;
}

// Distance between HyCal cluster and GEM hit after projecting GEM to HyCal z
float MatchingTools::ProjectionDistance(const analysis::HCHit &h,
                                       const analysis::GEMHit &g) const
{
    ProjectHit proj = GetProjectionHits(g.x, g.y, g.z, h.z);
    float dx = h.x - proj.x_proj;
    float dy = h.y - proj.y_proj;
    return std::sqrt(dx * dx + dy * dy);
}

// Distance between two GEM hits projected to a common reference z
float MatchingTools::ProjectionDistance(const analysis::GEMHit &g1,
                                       const analysis::GEMHit &g2,
                                       float ref_z) const
{
    ProjectHit p1 = GetProjectionHits(g1.x, g1.y, g1.z, ref_z);
    ProjectHit p2 = GetProjectionHits(g2.x, g2.y, g2.z, ref_z);
    float dx = p1.x_proj - p2.x_proj;
    float dy = p1.y_proj - p2.y_proj;
    return std::sqrt(dx * dx + dy * dy);
}

// ============================================================================
// Pre-match: check if a GEM hit falls within the matching window of a cluster
// ============================================================================

bool MatchingTools::PreMatch(const analysis::HCHit &hycal,
                             const analysis::GEMHit &gem) const
{
    ProjectHit proj = GetProjectionHits(gem.x, gem.y, gem.z, hycal_z_);
    float dx = std::fabs(hycal.x - proj.x_proj);
    float dy = std::fabs(hycal.y - proj.y_proj);

    if (squareSel_) {
        return (dx <= matchRange_) && (dy <= matchRange_);
    } else {
        return (dx * dx + dy * dy) <= matchRange_ * matchRange_;
    }
}

// ============================================================================
// Post-match: sort candidates per plane, set flags, pick best GEM hit
// ============================================================================

void MatchingTools::PostMatch(MatchHit &h) const
{
    if( (h.gem1_hits.empty() && h.gem2_hits.empty() ) ||
        (h.gem3_hits.empty() && h.gem4_hits.empty() ) )
        return; // require at least one match in both upstream and downstream pairs

    // sort each plane's candidates by projection distance (closest first)
    auto by_dist = [this, &h](const analysis::GEMHit &a, const analysis::GEMHit &b) {
        return ProjectionDistance(h.hycal_hit, a) < ProjectionDistance(h.hycal_hit, b);
    };
    std::sort(h.gem1_hits.begin(), h.gem1_hits.end(), by_dist);
    std::sort(h.gem2_hits.begin(), h.gem2_hits.end(), by_dist);
    std::sort(h.gem3_hits.begin(), h.gem3_hits.end(), by_dist);
    std::sort(h.gem4_hits.begin(), h.gem4_hits.end(), by_dist);

    // set match flag for each plane that has candidates
    if (!h.gem1_hits.empty()) fdec::set_bit(h.mflag, kGEM1Match);
    if (!h.gem2_hits.empty()) fdec::set_bit(h.mflag, kGEM2Match);
    if (!h.gem3_hits.empty()) fdec::set_bit(h.mflag, kGEM3Match);
    if (!h.gem4_hits.empty()) fdec::set_bit(h.mflag, kGEM4Match);

    // pick the best match from upstream pair (GEM1/GEM2)
    float best_up = 1e9f;
    analysis::GEMHit best_gem_up{};
    auto check_up = [&](const std::vector<analysis::GEMHit> &plane) {
        if (!plane.empty()) {
            float d = ProjectionDistance(h.hycal_hit, plane.front());
            if (d < best_up) {
                best_up = d;
                best_gem_up = plane.front();
            }
        }
    };
    check_up(h.gem1_hits);
    check_up(h.gem2_hits);
    h.gem[0] = best_gem_up;

    // pick the best match from downstream pair (GEM3/GEM4)
    float best_down = 1e9f;
    analysis::GEMHit best_gem_down{};
    auto check_down = [&](const std::vector<analysis::GEMHit> &plane) {
        if (!plane.empty()) {
            float d = ProjectionDistance(h.hycal_hit, plane.front());
            if (d < best_down) {
                best_down = d;
                best_gem_down = plane.front();
            }
        }
    };
    check_down(h.gem3_hits);
    check_down(h.gem4_hits);
    h.gem[1] = best_gem_down;
}

// ============================================================================
// Main matching — adapted from PRadDetMatch::Match for 4 GEM planes
// ============================================================================

// comparator so GEMHit can be stored in std::set (identity by position)
struct GEMHitCmp {
    bool operator()(const analysis::GEMHit &a, const analysis::GEMHit &b) const
    {
        if (a.z != b.z) return a.z < b.z;
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y;
    }
};

std::vector<MatchHit> MatchingTools::Match(
    std::vector<analysis::HCHit> &hycalHits,
    const std::vector<analysis::GEMHit> &gem1,
    const std::vector<analysis::GEMHit> &gem2,
    const std::vector<analysis::GEMHit> &gem3,
    const std::vector<analysis::GEMHit> &gem4) const
{
    std::vector<MatchHit> result;

    // keep track of GEM hits already claimed (higher-E cluster gets priority)
    std::set<analysis::GEMHit, GEMHitCmp> used1, used2, used3, used4;

    // sort HyCal clusters by energy descending — highest energy matched first
    std::sort(hycalHits.begin(), hycalHits.end(),
              [](const analysis::HCHit &a, const analysis::HCHit &b) {
                  return b.energy < a.energy;
              });

    for (size_t i = 0; i < hycalHits.size(); ++i) {
        const auto &hit = hycalHits[i];

        // collect candidates in each GEM plane
        std::vector<analysis::GEMHit> cand1, cand2, cand3, cand4;

        for (const auto &g : gem1)
            if (PreMatch(hit, g) && used1.find(g) == used1.end())
                cand1.push_back(g);
        for (const auto &g : gem2)
            if (PreMatch(hit, g) && used2.find(g) == used2.end())
                cand2.push_back(g);
        for (const auto &g : gem3)
            if (PreMatch(hit, g) && used3.find(g) == used3.end())
                cand3.push_back(g);
        for (const auto &g : gem4)
            if (PreMatch(hit, g) && used4.find(g) == used4.end())
                cand4.push_back(g);

        // skip if no candidates in any plane
        if (cand1.empty() && cand2.empty() && cand3.empty() && cand4.empty())
            continue;

        result.emplace_back(hit, cand1, cand2, cand3, cand4);
        MatchHit &mhit = result.back();
        mhit.hycal_idx = i;

        // resolve best match and set flags
        PostMatch(mhit);

        // mark the winning GEM hit in each flagged plane as used
        if (fdec::test_bit(mhit.mflag, kGEM1Match))
            used1.insert(mhit.gem1_hits.front());
        if (fdec::test_bit(mhit.mflag, kGEM2Match))
            used2.insert(mhit.gem2_hits.front());
        if (fdec::test_bit(mhit.mflag, kGEM3Match))
            used3.insert(mhit.gem3_hits.front());
        if (fdec::test_bit(mhit.mflag, kGEM4Match))
            used4.insert(mhit.gem4_hits.front());
    }

    return result;
}

} // namespace analysis