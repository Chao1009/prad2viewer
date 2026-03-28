#pragma once
//=============================================================================
// GemCluster.h — GEM strip clustering and 2D hit reconstruction
//
// Algorithms ported from mpd_gem_view_ssp GEMCluster:
//   1. Group consecutive strip hits
//   2. Split multi-peak clusters at local minima
//   3. Charge-weighted position reconstruction
//   4. Cross-talk identification by characteristic distance
//   5. Cartesian X-Y cluster matching → 2D hits
//=============================================================================

#include "GemSystem.h"
#include <vector>

namespace gem
{

struct ClusterConfig {
    int   min_cluster_hits  = 1;
    int   max_cluster_hits  = 20;
    int   consecutive_thres = 1;    // max gap between consecutive strips
    float split_thres       = 14.f; // charge valley depth for splitting
    float cross_talk_width  = 2.f;  // mm
    float position_res      = 0.08f;// mm
    std::vector<float> charac_dists;// cross-talk characteristic distances

    // XY matching cuts (Cartesian reconstruction)
    float match_adc_asymmetry = 0.8f;  // max |Qx-Qy|/(Qx+Qy), <0 to disable
    float match_time_diff     = 50.f;  // max |mean_t_x - mean_t_y| in ns, <0 to disable
    float ts_period           = 25.f;  // ns per time sample
};

class GemCluster
{
public:
    GemCluster();
    ~GemCluster();

    void SetConfig(const ClusterConfig &cfg) { cfg_ = cfg; }
    const ClusterConfig &GetConfig() const   { return cfg_; }

    // Form 1D strip clusters from a list of strip hits.
    // hits will be sorted by strip number.
    void FormClusters(std::vector<StripHit> &hits,
                      std::vector<StripCluster> &clusters) const;

    // Match X and Y clusters to form 2D hits via Cartesian product.
    void CartesianReconstruct(const std::vector<StripCluster> &x_clusters,
                              const std::vector<StripCluster> &y_clusters,
                              std::vector<GEMHit> &hits,
                              int det_id, float resolution) const;

private:
    // Group consecutive hits into preliminary clusters
    void groupHits(std::vector<StripHit> &hits,
                   std::vector<StripCluster> &clusters) const;

    // Recursively split clusters at local charge minima
    void splitCluster(std::vector<StripHit>::iterator beg,
                      std::vector<StripHit>::iterator end,
                      float thres,
                      std::vector<StripCluster> &clusters) const;

    // Compute charge-weighted position for a cluster
    void reconstructCluster(StripCluster &cluster) const;

    // Mark cross-talk clusters by characteristic distance
    void setCrossTalk(std::vector<StripCluster> &clusters) const;

    // Filter out bad clusters
    void filterClusters(std::vector<StripCluster> &clusters) const;

    ClusterConfig cfg_;
};

} // namespace gem
