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

namespace prad2 {

// ── Capacity constants ───────────────────────────────────────────────────

static constexpr int kMaxChannels  = fdec::MAX_ROCS * fdec::MAX_SLOTS * 16;
static constexpr int kMaxGemStrips = ssp::MAX_MPDS * ssp::MAX_APVS_PER_MPD * ssp::APV_STRIP_SIZE;
static constexpr int kMaxClusters  = 100;
static constexpr int kMaxGemHits   = 400;

// ── Raw replay ("events" tree) ───────────────────────────────────────────

struct RawEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger      = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    int64_t  timestamp    = 0;

    // HyCal per-channel data
    int          nch = 0;
    uint8_t crate[kMaxChannels]   = {};
    uint8_t slot[kMaxChannels]    = {};
    uint8_t channel[kMaxChannels] = {};
    uint16_t     module_id[kMaxChannels] = {};
    uint8_t nsamples[kMaxChannels] = {};
    uint16_t     samples[kMaxChannels][fdec::MAX_SAMPLES] = {};
    float   ped_mean[kMaxChannels] = {};
    float   ped_rms[kMaxChannels]  = {};
    float   integral[kMaxChannels] = {};

    // Optional peak data
    uint8_t npeaks[kMaxChannels] = {};
    float   peak_height[kMaxChannels][fdec::MAX_PEAKS]   = {};
    float   peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    float   peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};

    // GEM per-strip data
    int        gem_nch = 0;
    uint8_t mpd_crate[kMaxGemStrips]  = {};
    uint8_t mpd_fiber[kMaxGemStrips]  = {};
    uint8_t apv[kMaxGemStrips]        = {};
    uint8_t strip[kMaxGemStrips]      = {};
    int16_t ssp_samples[kMaxGemStrips][ssp::SSP_TIME_SAMPLES] = {};

    // ssp trigger bank tag
    uint8_t n_ssp_triggers = 0;
    uint32_t ssp_trigger_tags[ssp::SSP_TIME_SAMPLES] = {};
};

// ── Reconstructed replay ("recon" tree) ──────────────────────────────────

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    int64_t  timestamp    = 0;

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

    // GEM reconstructed hits
    int        n_gem_hits = 0;
    uint8_t det_id[kMaxGemHits]       = {};
    float   gem_x[kMaxGemHits]        = {};
    float   gem_y[kMaxGemHits]        = {};
    float   gem_x_charge[kMaxGemHits] = {};
    float   gem_y_charge[kMaxGemHits] = {};
    float   gem_x_peak[kMaxGemHits]   = {};
    float   gem_y_peak[kMaxGemHits]   = {};
    uint8_t gem_x_size[kMaxGemHits]   = {};
    uint8_t gem_y_size[kMaxGemHits]   = {};

    // Detector matching results
    int      match_num = 0;
    float matchHC_x[kMaxClusters] = {};
    float matchHC_y[kMaxClusters] = {};
    float matchHC_z[kMaxClusters] = {};
    float matchHC_energy[kMaxClusters] = {};
    uint16_t matchHC_center[kMaxClusters] = {};
    uint32_t matchHC_flag[kMaxClusters] = {};
    float matchG_x[kMaxClusters][2] = {}; // up/down GEM for each cluster
    float matchG_y[kMaxClusters][2] = {};
    float matchG_z[kMaxClusters][2] = {};
    uint8_t matchG_det_id[kMaxClusters][2] = {};

    // ssp trigger bank tag
    uint8_t n_ssp_triggers = 0;
    uint32_t ssp_trigger_tags[ssp::SSP_TIME_SAMPLES] = {};
};

} // namespace prad2
