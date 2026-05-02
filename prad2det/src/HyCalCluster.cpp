//=============================================================================
// HyCalCluster.cpp — Island clustering for HyCal
//
// Ported from PRadIslandCluster / PRadHyCalReconstructor (PRadAnalyzer).
// Uses pre-computed neighbor lists from HyCalSystem for fast adjacency checks.
//
// Chao Peng (original PRadAnalyzer), adapted for prad2decoder.
//=============================================================================

#include "HyCalCluster.h"
#include <algorithm>
#include <cmath>

namespace fdec
{

static constexpr int ISLAND_GROUP_RESERVE = 50;
static constexpr int POS_RECON_HITS       = 15;

//=============================================================================
// shower_depth (free function — declared in HyCalCluster.h)
//
// Maximum-shower-development depth into the calorimeter face for an EM
// shower of energy E, t = X0 · (ln(E/Ec) − Cf), with Cf = 0.5 for photons.
// Constants verbatim from the legacy analysis::PhysicsTools::GetShowerDepth
// (PRadAnalyzer lineage) and not expected to change — they're physical
// properties of the W and G modules.
//=============================================================================
float shower_depth(int center_id, float energy_mev)
{
    if (energy_mev <= 0.f) return 0.f;
    if (center_id >= PWO_ID0)            // PbWO4 (W-modules)
        return 8.6f  * (std::log(energy_mev / 1.1f)  - 0.5f);
    return            26.7f * (std::log(energy_mev / 2.84f) - 0.5f);  // PbGlass
}

//=============================================================================
// Construction / setup
//=============================================================================

HyCalCluster::HyCalCluster(const HyCalSystem &sys)
    : sys_(sys)
    , profile_(new SimpleProfile())
    , owns_profile_(true)
{
}

HyCalCluster::~HyCalCluster()
{
    if (owns_profile_) delete profile_;
}

void HyCalCluster::SetProfile(IClusterProfile *prof)
{
    if (owns_profile_) delete profile_;
    profile_ = prof;
    owns_profile_ = false;
}

//=============================================================================
// Per-event interface
//=============================================================================

void HyCalCluster::Clear()
{
    hits_.clear();
    groups_.clear();
    clusters_.clear();
}

void HyCalCluster::AddHit(int module_index, float energy, float time)
{
    if (module_index < 0 || module_index >= sys_.module_count()) return;
    if (energy > config_.min_module_energy)
        hits_.push_back({module_index, energy, time});
}

void HyCalCluster::FormClusters()
{
    clusters_.clear();
    groups_.clear();

    // step 1: group adjacent hits using DFS
    group_hits();

    // step 2: find maxima and split each group
    for (auto &group : groups_)
        split_cluster(group);
}

void HyCalCluster::ReconstructHits(std::vector<ClusterHit> &out) const
{
    out.clear();
    out.reserve(clusters_.size());

    for (auto &cl : clusters_) {
        if (cl.energy < config_.min_cluster_energy) continue;
        if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) continue;
        out.push_back(reconstruct_pos(cl));
    }
}

void HyCalCluster::ReconstructMatched(std::vector<RecoResult> &out) const
{
    out.clear();
    out.reserve(clusters_.size());

    for (auto &cl : clusters_) {
        if (cl.energy < config_.min_cluster_energy) continue;
        if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) continue;
        out.push_back({&cl, reconstruct_pos(cl)});
    }
}

//=============================================================================
// DFS grouping — form connected components of adjacent hits
//=============================================================================

void HyCalCluster::group_hits()
{
    // Build reverse map: module index → hit index
    mod_to_hit_.assign(sys_.module_count(), -1);
    for (size_t i = 0; i < hits_.size(); ++i)
        mod_to_hit_[hits_[i].index] = static_cast<int>(i);

    hit_group_id_.assign(hits_.size(), -1);
    std::vector<bool> visited(hits_.size(), false);

    for (size_t i = 0; i < hits_.size(); ++i) {
        if (visited[i]) continue;
        int gid = static_cast<int>(groups_.size());
        groups_.emplace_back();
        groups_.back().reserve(ISLAND_GROUP_RESERVE);
        dfs_group(groups_.back(), gid, static_cast<int>(i), visited);
    }
}

void HyCalCluster::dfs_group(std::vector<int> &group, int group_id, int hit_idx,
                              std::vector<bool> &visited)
{
    group.push_back(hit_idx);
    visited[hit_idx] = true;
    hit_group_id_[hit_idx] = group_id;

    sys_.for_each_neighbor(hits_[hit_idx].index, config_.corner_conn, [&](int ni) {
        int hi = mod_to_hit_[ni];
        if (hi >= 0 && !visited[hi])
            dfs_group(group, group_id, hi, visited);
    });
}

//=============================================================================
// Split cluster — find local maxima, distribute hits
//=============================================================================

void HyCalCluster::split_cluster(const std::vector<int> &group)
{
    auto maxima = find_maxima(group);
    if (maxima.empty()) return;

    if (maxima.size() == 1 ||
        group.size() >= SPLIT_MAX_HITS ||
        maxima.size() >= SPLIT_MAX_MAXIMA)
    {
        // single cluster from this group
        auto &seed = hits_[maxima[0]];
        clusters_.emplace_back();
        auto &cl = clusters_.back();
        cl.center = seed;
        cl.flag   = sys_.module(seed.index).flag;

        for (int hi : group)
            cl.add_hit(hits_[hi]);
    }
    else {
        split_hits(maxima, group);
    }
}

std::vector<int> HyCalCluster::find_maxima(const std::vector<int> &group) const
{
    std::vector<int> local_max;
    local_max.reserve(20);
    if (group.empty()) return local_max;

    const int gid = hit_group_id_[group[0]];

    for (int hi : group) {
        auto &hit = hits_[hi];
        if (hit.energy < config_.min_center_energy)
            continue;

        bool is_max = true;
        // include corners when checking for maxima (same as old code)
        sys_.for_each_neighbor(hit.index, true, [&](int ni) {
            if (!is_max) return;
            int hj = mod_to_hit_[ni];
            if (hj >= 0 && hit_group_id_[hj] == gid && hits_[hj].energy > hit.energy)
                is_max = false;
        });

        if (is_max)
            local_max.push_back(hi);
    }

    return local_max;
}

//=============================================================================
// Hit splitting — distribute shared hits among multiple maxima
//=============================================================================

void HyCalCluster::split_hits(const std::vector<int> &maxima,
                               const std::vector<int> &group)
{
    SplitContainer split;  // ~4KB on stack, safe per-call

    int nmax  = static_cast<int>(maxima.size());
    int nhits = static_cast<int>(group.size());

    // initialize fractions from profile
    for (int i = 0; i < nmax; ++i) {
        auto &center = hits_[maxima[i]];
        for (int j = 0; j < nhits; ++j) {
            auto &hit = hits_[group[j]];
            split.frac[j][i] = get_profile_frac(center, hit) * center.energy;
        }
    }

    // iterative refinement
    eval_fraction(maxima, group, split);

    // create clusters from final fractions
    for (int i = 0; i < nmax; ++i) {
        clusters_.emplace_back();
        auto &cl = clusters_.back();
        cl.center = hits_[maxima[i]];
        cl.flag   = sys_.module(cl.center.index).flag;

        for (int j = 0; j < nhits; ++j) {
            if (split.frac[j][i] == 0.f) continue;

            float nf = split.norm_frac(i, j);
            if (nf < config_.least_split) {
                split.total[j] -= split.frac[j][i];
                continue;
            }

            ModuleHit new_hit = hits_[group[j]];
            new_hit.energy *= nf;
            cl.add_hit(new_hit);

            if (new_hit.index == cl.center.index)
                cl.center.energy = new_hit.energy;

            set_bit(cl.flag, kSplit);
        }
    }
}

void HyCalCluster::eval_fraction(const std::vector<int> &maxima,
                                  const std::vector<int> &group,
                                  SplitContainer &split) const
{
    int nmax  = static_cast<int>(maxima.size());
    int nhits = static_cast<int>(group.size());

    struct BaseHit { float x, y, E; };
    BaseHit temp[POS_RECON_HITS];

    int iters = config_.split_iter;
    while (iters-- > 0) {
        split.sum_frac(nhits, nmax);

        for (int i = 0; i < nmax; ++i) {
            auto &center_hit = hits_[maxima[i]];
            const auto &center_mod = sys_.module(center_hit.index);

            // gather 3x3 neighbors for position reconstruction
            float tot_E = center_hit.energy;
            int count = 0;

            for (int j = 0; j < nhits; ++j) {
                auto &hit = hits_[group[j]];
                if (hit.index == center_hit.index || split.frac[j][i] == 0.f)
                    continue;

                // check if within 3x3 using pre-computed neighbors
                const auto &hit_mod = sys_.module(hit.index);
                double dx, dy;
                sys_.qdist(center_mod, hit_mod, dx, dy);

                if (std::abs(dx) < 1.01 && std::abs(dy) < 1.01 && count < POS_RECON_HITS) {
                    float frac_E = hit.energy * split.norm_frac(i, j);
                    temp[count] = {static_cast<float>(dx), static_cast<float>(dy), frac_E};
                    tot_E += frac_E;
                    count++;
                }
            }

            // reconstruct position (log-weighted)
            float wx = 0.f, wy = 0.f;
            float wtot = get_weight(center_hit.energy, tot_E);

            for (int k = 0; k < count; ++k) {
                float w = get_weight(temp[k].E, tot_E);
                if (w > 0.f) {
                    wx += temp[k].x * w;
                    wy += temp[k].y * w;
                    wtot += w;
                }
            }

            float cx = center_mod.x, cy = center_mod.y;
            if (wtot > 0.f) {
                cx += (wx / wtot) * center_mod.size_x;
                cy += (wy / wtot) * center_mod.size_y;
            }

            // update fractions with new center position
            for (int j = 0; j < nhits; ++j) {
                auto &hit = hits_[group[j]];
                split.frac[j][i] = get_profile_frac_at(cx, cy, tot_E, hit) * tot_E;
            }
        }
    }
    split.sum_frac(nhits, nmax);
}

//=============================================================================
// Position reconstruction — log-weighted centroid
//=============================================================================

ClusterHit HyCalCluster::reconstruct_pos(const ModuleCluster &cl) const
{
    const auto &center_mod = sys_.module(cl.center.index);

    struct BaseHit { float x, y, E; };
    BaseHit temp[POS_RECON_HITS];
    int count = 0;

    // gather 3x3 neighbors
    for (auto &hit : cl.hits) {
        if (hit.index == cl.center.index) continue;
        if (count >= POS_RECON_HITS) break;

        double dx, dy;
        sys_.qdist(center_mod, sys_.module(hit.index), dx, dy);
        if (std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
            temp[count++] = {static_cast<float>(dx), static_cast<float>(dy), hit.energy};
        }
    }

    // total energy
    float tot_E = cl.energy;

    // weighted position
    float wx = 0.f, wy = 0.f;
    float wtot = get_weight(cl.center.energy, tot_E);
    int npos = (wtot > 0.f) ? 1 : 0;

    for (int i = 0; i < count; ++i) {
        float w = get_weight(temp[i].E, tot_E);
        if (w > 0.f) {
            wx += temp[i].x * w;
            wy += temp[i].y * w;
            wtot += w;
            npos++;
        }
    }

    ClusterHit result;
    result.center_id = center_mod.id;
    result.energy    = cl.energy;
    result.nblocks   = static_cast<int>(cl.hits.size());
    result.flag      = cl.flag;

    if (wtot > 0.f) {
        result.x = center_mod.x + (wx / wtot) * center_mod.size_x;
        result.y = center_mod.y + (wy / wtot) * center_mod.size_y;
    } else {
        result.x = center_mod.x;
        result.y = center_mod.y;
    }
    result.npos = npos;

    return result;
}

float HyCalCluster::get_weight(float E, float E_total) const
{
    if (E_total <= 0.f) return 0.f;
    float w = config_.log_weight_thres + std::log(E / E_total);
    return (w > 0.f) ? w : 0.f;
}

//=============================================================================
// Profile helpers
//=============================================================================

float HyCalCluster::get_profile_frac(const ModuleHit &center, const ModuleHit &hit) const
{
    const auto &m1 = sys_.module(center.index);
    const auto &m2 = sys_.module(hit.index);
    double dx, dy;
    sys_.qdist(m1, m2, dx, dy);
    float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    // center module contains ~78% of energy, scale up to estimate total
    return profile_->GetFraction(m1.type, dist, center.energy / 0.78f);
}

float HyCalCluster::get_profile_frac_at(float cx, float cy, float cE,
                                          const ModuleHit &hit) const
{
    const auto &m = sys_.module(hit.index);
    int sid = sys_.get_sector_id(cx, cy);
    double dx, dy;
    sys_.qdist(cx, cy, sid, m.x, m.y, m.sector, dx, dy);
    float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    ModuleType type = sys_.sector_info(sid).mtype;
    return profile_->GetFraction(type, dist, cE);
}

} // namespace fdec
