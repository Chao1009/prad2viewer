#pragma once
//=============================================================================
// GemSystem.h — GEM detector system for PRad-II
//
// Manages the GEM detector hierarchy, DAQ↔detector channel mapping,
// pedestal subtraction, common mode correction, zero suppression,
// and strip hit collection.
//
// Usage:
//   gem::GemSystem sys;
//   sys.Init("gem_map.json");
//   sys.LoadPedestals("gem_ped.dat");
//
//   // per-event:
//   sys.Clear();
//   sys.ProcessEvent(ssp_evt);        // decoded SSP data
//   sys.Reconstruct(cluster);         // clustering + 2D hits
//   auto &hits = sys.GetHits(det_id); // reconstructed GEM hits
//=============================================================================

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

// Forward-declare SSP data types (from prad2dec)
namespace ssp { struct SspEventData; struct ApvData; }

namespace gem
{

// --- data structures --------------------------------------------------------

struct StripHit {
    int32_t strip;          // plane-wise strip number
    float   charge;         // max charge across time samples
    short   max_timebin;    // time sample with max charge
    float   position;       // physical position in mm
    bool    cross_talk;
    std::vector<float> ts_adc;  // all time sample ADC values
};

struct StripCluster {
    float   position;       // charge-weighted position (mm)
    float   peak_charge;    // highest strip charge in cluster
    float   total_charge;   // sum of all strip charges
    short   max_timebin;
    bool    cross_talk;
    std::vector<StripHit> hits;
};

struct GEMHit {
    float x, y, z;
    int   det_id;
    float x_charge, y_charge;
    float x_peak,   y_peak;
    short x_max_timebin, y_max_timebin;
    int   x_size, y_size;
    float sig_pos;          // position resolution
};

// --- APV pedestal -----------------------------------------------------------

struct ApvPedestal {
    float offset = 0.f;
    float noise  = 5000.f;     // large default → no hits until calibrated
};

// --- configuration ----------------------------------------------------------

struct ApvConfig {
    // DAQ address
    int crate_id    = -1;
    int mpd_id      = -1;
    int adc_ch      = -1;

    // Detector mapping
    int det_id      = -1;      // detector index
    int plane_type  = -1;      // 0=X, 1=Y
    int orient      = 0;       // 0 or 1 (strip reversal)
    int plane_index = -1;      // APV position on plane
    int det_pos     = 0;       // detector position in layer

    // Strip mapping parameters (APV channel → physical strip)
    // readout_offset: center point for APV25-to-strip odd/even mapping.
    //   32 = normal, 48 for special APVs. 0 = skip readout mapping.
    int  readout_offset = 32;
    // strip_offset: total offset added to plane-wide strip number (default 0).
    //   Accounts for all wiring differences (e.g. -144 for pos-11:
    //   -128 shares X column with pos-10, -16 for disconnected strips).
    int  strip_offset   = 0;
    bool hybrid_board   = true; // apply hybrid board pin conversion (MPD electronics)
    // match: half-strip intersection constraint for beam hole region.
    //   "" = full strip (intersects all cross-plane strips).
    //   "+Y" = upper half only (above beam hole).
    //   "-Y" = lower half only (below beam hole).
    std::string match;

    // Pedestals (per-strip)
    ApvPedestal pedestal[128];

    // Common mode range (for Danning algorithm)
    float cm_range_min = 0.f;
    float cm_range_max = 5000.f;
};

struct PlaneConfig {
    int   type      = -1;       // 0=X, 1=Y
    float size      = 0.f;      // mm
    int   n_apvs    = 0;        // number of APVs on this plane
    float pitch     = 0.4f;     // strip pitch (mm)
};

struct DetectorConfig {
    std::string name;
    int    id         = -1;
    std::string type  = "PRADGEM";  // detector type for strip mapping
    PlaneConfig planes[2];          // [0]=X, [1]=Y
};

// --- GemSystem class --------------------------------------------------------

class GemCluster;   // forward declaration

class GemSystem
{
public:
    GemSystem();
    ~GemSystem();

    // --- initialization -----------------------------------------------------
    void Init(const std::string &map_file);
    void LoadPedestals(const std::string &ped_file);
    void LoadCommonModeRange(const std::string &cm_file);

    // --- per-event processing -----------------------------------------------
    void Clear();
    void ProcessEvent(const ssp::SspEventData &evt);
    void Reconstruct(GemCluster &clusterer);

    // --- accessors ----------------------------------------------------------
    int GetNDetectors() const { return static_cast<int>(detectors_.size()); }
    const std::vector<DetectorConfig>& GetDetectors() const { return detectors_; }

    const std::vector<StripHit>& GetPlaneHits(int det, int plane) const;
    const std::vector<StripCluster>& GetPlaneClusters(int det, int plane) const;
    const std::vector<GEMHit>& GetHits(int det) const;
    const std::vector<GEMHit>& GetAllHits() const { return all_hits_; }

    // DAQ→APV index lookup (O(1))
    int FindApvIndex(int crate, int mpd, int adc) const;

    // Configuration
    float GetCommonModeThreshold() const { return common_thres_; }
    float GetZeroSupThreshold()    const { return zerosup_thres_; }
    void  SetCommonModeThreshold(float v) { common_thres_ = v; }
    void  SetZeroSupThreshold(float v)    { zerosup_thres_ = v; }

private:
    // --- per-APV processing -------------------------------------------------
    void processApv(int apv_idx, const ssp::ApvData &data);
    float commonModeSorting(float *buf, int size, int apv_idx);
    float commonModeDanning(float *buf, int size, int apv_idx);
    void collectHits(int apv_idx);

    // --- strip mapping ------------------------------------------------------
    void buildStripMap(int apv_idx);

    // --- detector hierarchy -------------------------------------------------
    std::vector<DetectorConfig> detectors_;
    std::vector<ApvConfig> apvs_;
    std::unordered_map<uint64_t, int> apv_map_;  // packed(crate,mpd,adc) → apv index

    static uint64_t packApvKey(int crate, int mpd, int adc)
    {
        return (static_cast<uint64_t>(static_cast<uint16_t>(crate)) << 32) |
               (static_cast<uint64_t>(static_cast<uint16_t>(mpd))  << 16) |
               static_cast<uint64_t>(static_cast<uint16_t>(adc));
    }

    // --- per-APV working data (pre-allocated) -------------------------------
    static constexpr int APV_STRIP_SIZE   = 128;
    static constexpr int SSP_TIME_SAMPLES = 6;
    static constexpr int NUM_HIGH_STRIPS  = 20;  // for sorting CM algorithm

    struct ApvWorkData {
        float raw[APV_STRIP_SIZE * SSP_TIME_SAMPLES];
        bool  hit_pos[APV_STRIP_SIZE];
        int   strip_map[APV_STRIP_SIZE];    // APV channel → plane strip
    };
    std::vector<ApvWorkData> apv_work_;

    // --- per-plane data (hits + clusters) -----------------------------------
    struct PlaneData {
        std::vector<StripHit> hits;
        std::vector<StripCluster> clusters;
    };
    // plane_data_[det_id][plane_type]
    std::vector<std::array<PlaneData, 2>> plane_data_;

    // --- per-detector reconstructed hits ------------------------------------
    std::vector<std::vector<GEMHit>> det_hits_;
    std::vector<GEMHit> all_hits_;

    // --- thresholds ---------------------------------------------------------
    float common_thres_     = 20.f;
    float zerosup_thres_    = 5.f;
    float crosstalk_thres_  = 8.f;
    bool  online_zero_sup_  = false;
    float position_res_     = 0.08f;   // mm
};

} // namespace gem
