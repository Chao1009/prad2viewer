#pragma once
//=============================================================================
// HyCalCluster.h — Island clustering algorithm for HyCal
//
// Ported from PRadIslandCluster / PRadHyCalReconstructor (PRadAnalyzer).
// Operates on index-based hits referencing HyCalSystem modules.
//
// Usage:
//   HyCalSystem sys;
//   sys.Init("hycal_modules.json", "hycal_daq_map.json");
//
//   HyCalCluster clusterer(sys);
//   // per-event:
//   clusterer.Clear();
//   clusterer.AddHit(module_index, energy, time);
//   clusterer.FormClusters();
//   for (auto &cl : clusterer.GetClusters()) { ... }
//=============================================================================

#include "HyCalSystem.h"
#include <vector>
#include <cmath>

namespace fdec
{

// --- configuration ----------------------------------------------------------
struct ClusterConfig {
    // hit thresholds
    float min_module_energy  = 0.f;       // MeV, minimum hit energy
    float min_center_energy  = 10.f;      // MeV, minimum seed energy
    float min_cluster_energy = 50.f;      // MeV, minimum total cluster energy
    int   min_cluster_size   = 1;         // minimum number of hits in cluster

    // island algorithm
    bool  corner_conn        = false;     // include diagonal neighbors in grouping
    int   split_iter         = 6;         // iterations for fraction refinement
    float least_split        = 0.01f;     // minimum fraction to keep a split hit

    // position reconstruction
    float log_weight_thres   = 3.6f;      // W = max(0, thres + ln(E_i/E_tot))
};

// --- per-event hit ----------------------------------------------------------
struct ModuleHit {
    int   index;        // module index in HyCalSystem
    float energy;       // calibrated energy (MeV)
    float time;         // ADC peaking time (ns)
};

// --- cluster result ---------------------------------------------------------
struct ModuleCluster {
    ModuleHit              center;     // seed module (highest energy local max)
    std::vector<ModuleHit> hits;       // all hits (energy may be split)
    float                  energy = 0.f;
    float                  time   = 0.f; // peaking time of center modules (ns)
    uint32_t               flag   = 0;

    void add_hit(const ModuleHit &h)
    {
        hits.push_back(h);
        energy += h.energy;
        time = h.time;
    }
};

// --- reconstructed hit position ---------------------------------------------
struct ClusterHit {
    int   center_id;    // PrimEx ID of center module
    float x, y;         // reconstructed position (mm)
    float energy;       // total cluster energy (MeV)
    float time;         // ADC peaking time (ns) of center module
    int   nblocks;      // number of modules in cluster
    int   npos;         // number of modules used in position reconstruction
    uint32_t flag;      // cluster flags
};

// Shower-max depth into the calorimeter face for an EM shower of energy `E`
// (MeV) in the module with PrimEx id `center_id`.  Discriminates W/G by the
// PWO_ID0 boundary.  Returns 0 for E ≤ 0.  Units: mm.
//
//   t = X0 · (ln(E/Ec) − Cf)
//   Cf  = 0.5  for photon-induced showers
//   PWO4:    X0 = 8.6 mm, Ec = 1.1  MeV
//   PbGlass: X0 = 26.7 mm, Ec = 2.84 MeV
//
// (Same formula as the legacy analysis::PhysicsTools::GetShowerDepth, moved
// here so prad2det owns it and the python binding can expose it.)
float shower_depth(int center_id, float energy_mev);

// --- cluster profile (energy sharing lookup) --------------------------------
// Abstract interface — users can plug in their own profile data.
// Default implementation uses a simple analytical approximation.
struct IClusterProfile {
    virtual ~IClusterProfile() = default;
    // Returns the fraction of energy at quantized distance `dist` for a cluster
    // of total energy `energy` (MeV) on a module of given type.
    virtual float GetFraction(ModuleType type, float dist, float energy) const = 0;
};

// Simple analytical profile (exponential falloff in Moliere radius units)
struct SimpleProfile : public IClusterProfile {
    float GetFraction(ModuleType type, float dist, float /*energy*/) const override
    {
        // approximate transverse shower profile in quantized distance units
        // PbWO4 Moliere radius ~20mm ≈ module size, PbGlass ~38mm ≈ module size
        // so quantized distance ~1 corresponds to ~1 Moliere radius
        if (dist < 0.01f) return 0.78f;    // center module: ~78% of energy
        float sigma = (type == ModuleType::PbWO4) ? 0.36f : 0.40f;
        return 0.78f * std::exp(-dist * dist / (2.f * sigma * sigma));
    }
};

// --- split container (static, reused across calls) --------------------------
static constexpr int SPLIT_MAX_HITS   = 100;
static constexpr int SPLIT_MAX_MAXIMA = 10;

struct SplitContainer {
    float frac[SPLIT_MAX_HITS][SPLIT_MAX_MAXIMA];
    float total[SPLIT_MAX_HITS];

    void sum_frac(int nhits, int nmax)
    {
        for (int i = 0; i < nhits; ++i) {
            total[i] = 0.f;
            for (int j = 0; j < nmax; ++j)
                total[i] += frac[i][j];
        }
    }

    float norm_frac(int imax, int ihit) const
    {
        return (total[ihit] > 0.f) ? frac[ihit][imax] / total[ihit] : 0.f;
    }
};

// --- main clustering class --------------------------------------------------
class HyCalCluster
{
public:
    explicit HyCalCluster(const HyCalSystem &sys);
    ~HyCalCluster();

    // non-copyable (owns profile pointer)
    HyCalCluster(const HyCalCluster &) = delete;
    HyCalCluster &operator=(const HyCalCluster &) = delete;

    // set configuration
    void SetConfig(const ClusterConfig &cfg) { config_ = cfg; }
    const ClusterConfig &GetConfig() const   { return config_; }

    // set cluster profile (takes ownership)
    void SetProfile(IClusterProfile *prof);

    // --- per-event interface ------------------------------------------------
    void Clear();
    void AddHit(int module_index, float energy, float time);
    void FormClusters();
    void ReconstructHits(std::vector<ClusterHit> &out) const;

    // Reconstruct and return paired (cluster, hit) for clusters passing thresholds.
    // This avoids fragile parallel iteration between GetClusters() and ReconstructHits().
    struct RecoResult {
        const ModuleCluster *cluster;
        ClusterHit hit;
    };
    void ReconstructMatched(std::vector<RecoResult> &out) const;

    // access results
    const std::vector<ModuleCluster> &GetClusters() const { return clusters_; }

private:
    // island algorithm steps
    void group_hits();
    void dfs_group(std::vector<int> &group, int group_id, int hit_idx,
                   std::vector<bool> &visited);
    void split_cluster(const std::vector<int> &group);
    std::vector<int> find_maxima(const std::vector<int> &group) const;
    void split_hits(const std::vector<int> &maxima,
                    const std::vector<int> &group);
    void eval_fraction(const std::vector<int> &maxima,
                       const std::vector<int> &group,
                       SplitContainer &split) const;

    // position reconstruction
    ClusterHit reconstruct_pos(const ModuleCluster &cl) const;
    float get_weight(float E, float E_total) const;

    // profile helper
    float get_profile_frac(const ModuleHit &center, const ModuleHit &hit) const;
    float get_profile_frac_at(float cx, float cy, float cE,
                              const ModuleHit &hit) const;

    const HyCalSystem     &sys_;
    ClusterConfig          config_;
    IClusterProfile       *profile_;
    bool                   owns_profile_;

    // per-event data
    std::vector<ModuleHit>              hits_;
    std::vector<std::vector<int>>       groups_;        // groups of hit indices
    std::vector<ModuleCluster>          clusters_;
    std::vector<int>                    mod_to_hit_;    // module_index → hit_index (-1 if not hit)
    std::vector<int>                    hit_group_id_;  // hit_index → group_id
};

} // namespace fdec
