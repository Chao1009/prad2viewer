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
    int  pin_rotate  = 0;    // rotated connector pins (e.g. 16 for pos-11 near beam hole)
    int  shared_pos  = -1;   // effective plane position (-1 = use plane_index)
    bool hybrid_board = true; // hybrid board pin conversion (MPD electronics)
    // match: half-strip intersection constraint for beam hole region.
    //   "" = full strip, "+Y" = above hole, "-Y" = below hole.
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

// --- strip-mapping pipeline (pure, stateless) -------------------------------
//
// Maps one APV25 channel index to the plane-wide strip number using the
// 6-step pipeline shared with the original PRadAnalyzer / mpd_gem_view_ssp
// code.  Parameters mirror ApvConfig so off-line analyses (including the
// Python bindings) can reproduce the layout without instantiating a
// GemSystem.  GemSystem::buildStripMap() delegates to this function —
// single source of truth.
//
//   ch             0..apv_channels-1
//   plane_index    APV position on the plane (0..n_apvs-1)
//   orient         strip-reversal flag (0 or 1)
//   pin_rotate     rotated connector pins (e.g. 16 for pos-11 near beam hole)
//   shared_pos     effective plane position (-1 = use plane_index)
//   hybrid_board   true → apply the MPD hybrid-board pin conversion (step 2)
//   apv_channels   N = channels per APV chip (default 128)
//   readout_center default 32; combined with pin_rotate to form the readout
//                  offset; readout_off<=0 skips step 3
int MapStrip(int ch, int plane_index, int orient,
             int  pin_rotate    = 0,
             int  shared_pos    = -1,
             bool hybrid_board  = true,
             int  apv_channels  = 128,
             int  readout_center = 32);

// Convenience: compute plane-wide strip numbers for every channel of an APV.
// Equivalent to calling MapStrip(ch, ...) for ch in [0, apv_channels).
std::vector<int> MapApvStrips(int plane_index, int orient,
                              int  pin_rotate    = 0,
                              int  shared_pos    = -1,
                              bool hybrid_board  = true,
                              int  apv_channels  = 128,
                              int  readout_center = 32);

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

    // APV config access (for diagnostics/serialization)
    int GetNApvs() const { return static_cast<int>(apvs_.size()); }
    const ApvConfig& GetApvConfig(int idx) const { return apvs_[idx]; }

    // Compute beam hole X offset from detector center (mm), using match APV strip positions.
    // Returns 0 if no match APVs found. Uses detector 0 as reference (all identical).
    float GetHoleXOffset() const;

    // Per-APV zero-suppression results (valid after ProcessEvent)
    bool  IsChannelHit(int apv_idx, int ch) const { return apv_work_[apv_idx].hit_pos[ch]; }
    bool  HasApvZsHits(int apv_idx) const {
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch)
            if (apv_work_[apv_idx].hit_pos[ch]) return true;
        return false;
    }
    float GetProcessedAdc(int apv_idx, int ch, int ts) const {
        return apv_work_[apv_idx].raw[ts * APV_STRIP_SIZE + ch];
    }

    // Configuration
    float GetCommonModeThreshold() const { return common_thres_; }
    float GetZeroSupThreshold()    const { return zerosup_thres_; }
    float GetCrossTalkThreshold()  const { return crosstalk_thres_; }
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

    // --- global APV parameters -----------------------------------------------
    int   apv_channels_     = 128;     // channels per APV chip
    int   readout_center_   = 32;      // default readout mapping center

    // --- thresholds ---------------------------------------------------------
    float common_thres_     = 20.f;
    float zerosup_thres_    = 5.f;
    float crosstalk_thres_  = 8.f;
    bool  online_zero_sup_  = false;
    float position_res_     = 0.08f;   // mm

    // --- strip-level cuts (from mpd_gem_view_ssp) -------------------------
    bool  reject_first_timebin_ = true;   // reject if peak at first time bin
    bool  reject_last_timebin_  = true;   // reject if peak at last time bin
    float min_peak_adc_         = 0.f;    // min peak ADC per strip (0=disabled)
    float min_sum_adc_          = 0.f;    // min sum ADC per strip (0=disabled)

    // --- XY matching (passed to GemCluster at Reconstruct time) -----------
    int   match_mode_           = 1;     // 0=sorted, 1=Cartesian+cuts
    float match_adc_asymmetry_  = 0.8f;
    float match_time_diff_      = 50.f;
    float match_ts_period_      = 25.f;
};

} // namespace gem
