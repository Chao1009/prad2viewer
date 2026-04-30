#pragma once
//=============================================================================
// EventData.h — Shared data structures for ROOT replay trees
//
// Used by:
//   - analysis/Replay (writer: EVIO → ROOT)
//   - root_data_source (reader: ROOT → viewer)
//   - analysis tools (reader: ROOT → physics analysis)
//
// These structs define the branch layout of ROOT TTrees produced by
// replay_rawdata ("events" tree) and replay_recon ("recon" tree).
// Changing a struct here automatically updates all readers and writers.
//
// NOTE: No ROOT headers needed — uses standard C++ types only.
//       TTree branch setup uses these as plain arrays.
//=============================================================================

#include "Fadc250Data.h"   // MAX_SAMPLES, MAX_PEAKS, MAX_ROCS, MAX_SLOTS
#include "SspData.h"       // MAX_MPDS, MAX_APVS_PER_MPD, APV_STRIP_SIZE, SSP_TIME_SAMPLES

#include <cstdint>
#include <vector>

namespace prad2 {

// ── Capacity constants ───────────────────────────────────────────────────

static constexpr int kMaxChannels  = fdec::MAX_ROCS * fdec::MAX_SLOTS * 16;
static constexpr int kMaxGemStrips = ssp::MAX_MPDS * ssp::MAX_APVS_PER_MPD * ssp::APV_STRIP_SIZE;
static constexpr int kMaxClusters  = 100;
static constexpr int kMaxGemHits   = 400;

// ── Module type categorisation ────────────────────────────────────────────
//
// Single source of truth at the data-tree level.  Values come from the "t"
// field of hycal_modules.json (PbGlass / PbWO4 / SCINT / LMS), parsed at
// load time and stored per-channel in RawEventData::module_type.  Numeric
// values are arbitrary but stable — kept as a uint8_t so the TTree branch
// stays compact (1 byte per channel).
//
// Consumers categorise channels by this field; they should not parse the
// module name.  An unrecognised type falls through to MOD_UNKNOWN and the
// channel is still written to the tree (so nothing is silently dropped).
enum ModuleType : uint8_t {
    MOD_UNKNOWN = 0,
    MOD_PbGlass = 1,
    MOD_PbWO4   = 2,
    MOD_SCINT   = 3,   // Veto scintillators (V1..V4)
    MOD_LMS     = 4,   // LMS reference PMTs (LMSPin, LMS1..3)
};

// ── Raw replay ("events" tree) ───────────────────────────────────────────
//
// One flat FADC250 channel array.  HyCal, Veto, and LMS channels all live
// in the same arrays — distinguish via `module_type[i]`.  Branch prefix
// stays "hycal.*" for backwards compatibility with existing consumers; they
// do `hycal.module_by_id(...)` then `if (!mod || !mod->is_hycal()) continue;`,
// which transparently skips Veto/LMS entries because their module_id values
// (3001..3004 / 3100..3103) are not registered in HyCalSystem.
//
// module_id encoding (globally unique across types):
//   MOD_PbGlass : 1..1156      (matches HyCalSystem G-module IDs)
//   MOD_PbWO4   : 1001..2152   (HyCal W-module IDs + 1000)
//   MOD_SCINT   : 3001..3004   (V1..V4)
//   MOD_LMS     : 3100..3103   (LMSPin=3100, LMS1..3 = 3101..3103)
struct RawEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits      = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    // FADC250 per-channel data (HyCal + Veto + LMS, distinguished by module_type)
    int          nch = 0;
    uint16_t     module_id[kMaxChannels]   = {};
    uint8_t      module_type[kMaxChannels] = {};   // ModuleType enum value
    int          nsamples[kMaxChannels]    = {};
    uint16_t     samples[kMaxChannels][fdec::MAX_SAMPLES] = {};
    float        ped_mean[kMaxChannels]    = {};
    float        ped_rms[kMaxChannels]     = {};
    float        gain_factor[kMaxChannels] = {};   // 1.0 for non-HyCal types

    // Optional soft-analyzer peak data (gated on -p flag in replay_rawdata)
    int     npeaks[kMaxChannels] = {};
    float   peak_height[kMaxChannels][fdec::MAX_PEAKS]   = {};
    float   peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    float   peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};

    // Optional firmware-mode (FADC250 Modes 1/2/3) peak data — also gated on -p.
    // Produced by Fadc250FwAnalyzer using the soft pedestal mean as PED.
    //   daq_npeaks       — number of pulses kept (≤ Fadc250FwConfig.MAX_PULSES)
    //   daq_peak_vp      — Vpeak (pedestal-subtracted ADC counts)
    //   daq_peak_integral— Σ over [cross−NSB, cross+NSA] (Mode 2 integral)
    //   daq_peak_time    — interpolated mid-amplitude time (ns)
    //   daq_peak_cross   — Tcross sample index (Mode 1 "first sample number")
    //   daq_peak_pos     — sample index of Vp itself (different from Tcross
    //                      whenever the leading edge spans multiple samples)
    //   daq_peak_coarse  — 4-ns clock index of Vba (10-bit firmware field)
    //   daq_peak_fine    — sub-sample fine bits, 0..63 (62.5 ps LSB)
    //   daq_peak_quality — bitmask: Q_PEAK_AT_BOUNDARY|Q_NSB_TRUNCATED|
    //                      Q_NSA_TRUNCATED|Q_VA_OUT_OF_RANGE (see Fadc250Data.h)
    int     daq_npeaks[kMaxChannels] = {};
    float   daq_peak_vp[kMaxChannels][fdec::MAX_PEAKS]       = {};
    float   daq_peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};
    float   daq_peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    int     daq_peak_cross[kMaxChannels][fdec::MAX_PEAKS]    = {};
    int     daq_peak_pos[kMaxChannels][fdec::MAX_PEAKS]      = {};
    int     daq_peak_coarse[kMaxChannels][fdec::MAX_PEAKS]   = {};
    int     daq_peak_fine[kMaxChannels][fdec::MAX_PEAKS]     = {};
    uint8_t daq_peak_quality[kMaxChannels][fdec::MAX_PEAKS]  = {};

    // GEM per-strip data
    int        gem_nch = 0;
    uint8_t mpd_crate[kMaxGemStrips]  = {};
    uint8_t mpd_fiber[kMaxGemStrips]  = {};
    uint8_t apv[kMaxGemStrips]        = {};
    uint8_t strip[kMaxGemStrips]      = {};
    int16_t ssp_samples[kMaxGemStrips][ssp::SSP_TIME_SAMPLES] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;
};

// ── Reconstructed replay ("recon" tree) ──────────────────────────────────

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    // HyCal clusters
    float total_energy = 0.f;
    int     n_clusters = 0;
    float cl_x[kMaxClusters]       = {};
    float cl_y[kMaxClusters]       = {};
    float cl_z[kMaxClusters]       = {};
    float cl_energy[kMaxClusters]  = {};
    uint8_t cl_nblocks[kMaxClusters] = {};
    uint16_t cl_center[kMaxClusters]  = {};
    uint32_t cl_flag[kMaxClusters]    = {};
    // Matching results
    uint32_t matchFlag[kMaxClusters] = {};
    float    matchGEMx[kMaxClusters][4] = {};
    float    matchGEMy[kMaxClusters][4] = {};
    float    matchGEMz[kMaxClusters][4] = {};
    int      matchNum = 0; // number of clusters with matches (for quick access, can be derived from matchFlag)
    //for quick simple access to each matched hit on HC and GEM planes
    // HC_Energy, HC_x/y/z, GEM_x/y/z (in mm, beam center and target center coordinate)
    float    mHit_E[kMaxClusters] = {};
    float    mHit_x[kMaxClusters] = {};
    float    mHit_y[kMaxClusters] = {};
    float    mHit_z[kMaxClusters] = {};
    float    mHit_gx[kMaxClusters][2] = {};
    float    mHit_gy[kMaxClusters][2] = {};
    float    mHit_gz[kMaxClusters][2] = {};
    float    mHit_gid[kMaxClusters][2] = {}; //det_id for matched GEM hits

    // GEM reconstructed hits
    int        n_gem_hits = 0;
    uint8_t det_id[kMaxGemHits]       = {};
    float   gem_x[kMaxGemHits]        = {};
    float   gem_y[kMaxGemHits]        = {};
    float   gem_z[kMaxGemHits]        = {};
    float   gem_x_charge[kMaxGemHits] = {};
    float   gem_y_charge[kMaxGemHits] = {};
    float   gem_x_peak[kMaxGemHits]   = {};
    float   gem_y_peak[kMaxGemHits]   = {};
    uint8_t gem_x_size[kMaxGemHits]   = {};
    uint8_t gem_y_size[kMaxGemHits]   = {};
    uint8_t gem_x_mTbin[kMaxGemHits]   = {};
    uint8_t gem_y_mTbin[kMaxGemHits]   = {};

    //veto information
    int      veto_nch = 0;
    uint8_t veto_id[4]   = {}; // 0,1,2,3 for veto1-4
    int veto_npeaks[4] = {};
    float veto_peak_time[4][fdec::MAX_PEAKS]     = {};
    float veto_peak_integral[4][fdec::MAX_PEAKS] = {};

    //LMS reference PMT information
    int      lms_nch = 0;
    uint8_t lms_id[4]   = {}; // 0,1,2,3 for lms1-4
    int lms_npeaks[4] = {};
    float lms_peak_time[4][fdec::MAX_PEAKS]     = {};
    float lms_peak_integral[4][fdec::MAX_PEAKS] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;
};

} // namespace prad2
