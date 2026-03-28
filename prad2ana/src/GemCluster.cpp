//============================================================================//
// GEM Cluster Reconstruction                                                 //
//                                                                            //
// Ported from mpd_gem_view_ssp GEMCluster                                    //
// Original authors: Xinzhan Bai, Kondo Gnanvo, Chao Peng                     //
//============================================================================//

#include "GemCluster.h"
#include <algorithm>
#include <cmath>

using namespace gem;

//=============================================================================
// Construction / destruction
//=============================================================================

GemCluster::GemCluster()
{
    // Default characteristic distances for cross-talk identification (mm)
    cfg_.charac_dists = {6.4f, 17.6f, 24.4f, 24.8f, 25.2f, 25.6f,
                         26.0f, 26.4f, 26.8f, 33.6f, 44.8f};
}

GemCluster::~GemCluster() = default;

//=============================================================================
// FormClusters — main entry point
//=============================================================================

void GemCluster::FormClusters(std::vector<StripHit> &hits,
                              std::vector<StripCluster> &clusters) const
{
    clusters.clear();
    if (hits.empty()) return;

    // group consecutive hits → preliminary clusters (with splitting)
    groupHits(hits, clusters);

    // reconstruct cluster position
    for (auto &cluster : clusters)
        reconstructCluster(cluster);

    // mark cross-talk clusters
    setCrossTalk(clusters);

    // filter out bad clusters
    filterClusters(clusters);
}

//=============================================================================
// groupHits — sort by strip, then cluster consecutive strips
//=============================================================================

void GemCluster::groupHits(std::vector<StripHit> &hits,
                           std::vector<StripCluster> &clusters) const
{
    // sort by strip number
    std::sort(hits.begin(), hits.end(),
              [](const StripHit &a, const StripHit &b) {
                  return a.strip < b.strip;
              });

    // cluster consecutive hits
    auto cbeg = hits.begin();
    for (auto it = hits.begin(); it != hits.end(); ++it) {
        auto it_n = it + 1;
        if (it_n == hits.end() ||
            it_n->strip - it->strip > cfg_.consecutive_thres)
        {
            splitCluster(cbeg, it_n, cfg_.split_thres, clusters);
            cbeg = it_n;
        }
    }
}

//=============================================================================
// splitCluster — recursively split at local charge minima (valleys)
//=============================================================================

void GemCluster::splitCluster(std::vector<StripHit>::iterator beg,
                              std::vector<StripHit>::iterator end,
                              float thres,
                              std::vector<StripCluster> &clusters) const
{
    auto size = end - beg;
    if (size <= 0) return;

    // Don't split clusters smaller than 3 strips
    if (size < 3) {
        StripCluster cl;
        cl.hits.assign(beg, end);
        clusters.push_back(std::move(cl));
        return;
    }

    // Find the first local minimum (valley)
    bool descending = false, extremum = false;
    auto minimum = beg;

    for (auto it = beg, it_n = beg + 1; it_n != end; ++it, ++it_n) {
        if (descending) {
            if (it->charge < minimum->charge)
                minimum = it;
            // ascending trend confirms valley
            if (it_n->charge - it->charge > thres) {
                extremum = true;
                break;
            }
        } else {
            // descending trend — potential valley ahead
            if (it->charge - it_n->charge > thres) {
                descending = true;
                minimum = it_n;
            }
        }
    }

    if (extremum) {
        // halve the charge of the overlap strip
        minimum->charge /= 2.f;

        // left sub-cluster
        StripCluster cl;
        cl.hits.assign(beg, minimum);
        clusters.push_back(std::move(cl));

        // recurse on right portion
        splitCluster(minimum, end, thres, clusters);
    } else {
        StripCluster cl;
        cl.hits.assign(beg, end);
        clusters.push_back(std::move(cl));
    }
}

//=============================================================================
// reconstructCluster — charge-weighted position
//=============================================================================

void GemCluster::reconstructCluster(StripCluster &cluster) const
{
    if (cluster.hits.empty()) return;

    cluster.total_charge = 0.f;
    cluster.peak_charge  = 0.f;
    cluster.max_timebin  = -1;
    float weight_pos = 0.f;

    for (auto &hit : cluster.hits) {
        if (hit.charge > cluster.peak_charge) {
            cluster.peak_charge = hit.charge;
            cluster.max_timebin = hit.max_timebin;
        }
        cluster.total_charge += hit.charge;
        weight_pos += hit.position * hit.charge;
    }

    if (cluster.total_charge > 0.f)
        cluster.position = weight_pos / cluster.total_charge;
}

//=============================================================================
// setCrossTalk — mark clusters at characteristic cross-talk distances
//=============================================================================

namespace {

// Check if all hits in a cluster are cross-talk strips
inline bool isPureCrossTalk(const StripCluster &cl)
{
    for (auto &hit : cl.hits)
        if (!hit.cross_talk) return false;
    return true;
}

// Check if cluster is at a characteristic CT distance from any later cluster
inline bool atCTDistance(std::vector<StripCluster>::iterator it,
                        std::vector<StripCluster>::iterator end,
                        float width,
                        const std::vector<float> &charac)
{
    for (auto itn = it + 1; itn != end; ++itn) {
        float delta = std::abs(it->position - itn->position);
        for (float dist : charac) {
            if (delta > dist - width && delta < dist + width)
                return true;
        }
    }
    return false;
}

} // anonymous namespace

void GemCluster::setCrossTalk(std::vector<StripCluster> &clusters) const
{
    if (cfg_.charac_dists.empty()) return;

    // sort by peak charge ascending (check weakest clusters first)
    std::sort(clusters.begin(), clusters.end(),
              [](const StripCluster &a, const StripCluster &b) {
                  return a.peak_charge < b.peak_charge;
              });

    for (auto it = clusters.begin(); it != clusters.end(); ++it) {
        if (!isPureCrossTalk(*it)) continue;
        it->cross_talk = atCTDistance(it, clusters.end(),
                                     cfg_.cross_talk_width, cfg_.charac_dists);
    }
}

//=============================================================================
// filterClusters — remove bad clusters
//=============================================================================

void GemCluster::filterClusters(std::vector<StripCluster> &clusters) const
{
    clusters.erase(
        std::remove_if(clusters.begin(), clusters.end(),
            [this](const StripCluster &cl) {
                // bad size
                int sz = static_cast<int>(cl.hits.size());
                if (sz < cfg_.min_cluster_hits || sz > cfg_.max_cluster_hits)
                    return true;
                // cross-talk
                if (cl.cross_talk)
                    return true;
                return false;
            }),
        clusters.end());
}

//=============================================================================
// seedMeanTime — ADC-weighted mean time of the seed (max-charge) strip
//=============================================================================

static float seedMeanTime(const StripCluster &cl, float ts_period)
{
    // find seed strip (highest charge)
    const StripHit *seed = nullptr;
    for (auto &h : cl.hits)
        if (!seed || h.charge > seed->charge) seed = &h;
    if (!seed || seed->ts_adc.empty()) return -1.f;

    float sum_wt = 0.f, sum_w = 0.f;
    for (size_t i = 0; i < seed->ts_adc.size(); ++i) {
        float w = seed->ts_adc[i];
        if (w > 0.f) {
            sum_wt += w * static_cast<float>(i + 1) * ts_period;
            sum_w  += w;
        }
    }
    return (sum_w > 0.f) ? sum_wt / sum_w : -1.f;
}

//=============================================================================
// CartesianReconstruct — match X and Y clusters to form 2D hits
//
// Full Cartesian product with optional cuts from mpd_gem_view_ssp:
//   - ADC asymmetry: |Qx_peak - Qy_peak| / (Qx_peak + Qy_peak) <= threshold
//   - Timing:        |mean_time_x_seed - mean_time_y_seed| <= threshold
//=============================================================================

void GemCluster::CartesianReconstruct(
    const std::vector<StripCluster> &x_clusters,
    const std::vector<StripCluster> &y_clusters,
    std::vector<GEMHit> &container,
    int det_id, float resolution) const
{
    container.clear();

    const float adc_asym_cut = cfg_.match_adc_asymmetry;
    const float time_cut     = cfg_.match_time_diff;
    const float ts_ns        = cfg_.ts_period;

    for (auto &xc : x_clusters) {
        float x_mean_t = (time_cut >= 0.f) ? seedMeanTime(xc, ts_ns) : 0.f;

        for (auto &yc : y_clusters) {
            // ADC asymmetry cut
            if (adc_asym_cut >= 0.f) {
                float sum = xc.peak_charge + yc.peak_charge;
                if (sum > 0.f) {
                    float asym = std::abs(xc.peak_charge - yc.peak_charge) / sum;
                    if (asym > adc_asym_cut) continue;
                }
            }

            // timing asymmetry cut
            if (time_cut >= 0.f && x_mean_t >= 0.f) {
                float y_mean_t = seedMeanTime(yc, ts_ns);
                if (y_mean_t >= 0.f && std::abs(x_mean_t - y_mean_t) > time_cut)
                    continue;
            }

            GEMHit hit;
            hit.x = xc.position;
            hit.y = yc.position;
            hit.z = 0.f;
            hit.det_id = det_id;
            hit.x_charge = xc.total_charge;
            hit.y_charge = yc.total_charge;
            hit.x_peak   = xc.peak_charge;
            hit.y_peak   = yc.peak_charge;
            hit.x_max_timebin = xc.max_timebin;
            hit.y_max_timebin = yc.max_timebin;
            hit.x_size = static_cast<int>(xc.hits.size());
            hit.y_size = static_cast<int>(yc.hits.size());
            hit.sig_pos = resolution;
            container.push_back(hit);
        }
    }
}
