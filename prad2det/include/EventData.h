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
    int     nch = 0;
    uint8_t crate[kMaxChannels]   = {};
    uint8_t slot[kMaxChannels]    = {};
    uint8_t channel[kMaxChannels] = {};
    int     module_id[kMaxChannels] = {};
    uint8_t nsamples[kMaxChannels] = {};
    int     samples[kMaxChannels][fdec::MAX_SAMPLES] = {};
    float   ped_mean[kMaxChannels] = {};
    float   ped_rms[kMaxChannels]  = {};
    float   integral[kMaxChannels] = {};

    // Optional peak data
    uint8_t npeaks[kMaxChannels] = {};
    float   peak_height[kMaxChannels][fdec::MAX_PEAKS]   = {};
    float   peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    float   peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};

    // GEM per-strip data
    int     gem_nch = 0;
    uint8_t mpd_crate[kMaxGemStrips]  = {};
    uint8_t mpd_fiber[kMaxGemStrips]  = {};
    uint8_t apv[kMaxGemStrips]        = {};
    uint8_t strip[kMaxGemStrips]      = {};
    float   ssp_samples[kMaxGemStrips][ssp::SSP_TIME_SAMPLES] = {};
};

// ── Reconstructed replay ("recon" tree) ──────────────────────────────────

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    int64_t  timestamp    = 0;

    // HyCal clusters
    int   n_clusters = 0;
    float cl_x[kMaxClusters]       = {};
    float cl_y[kMaxClusters]       = {};
    float cl_energy[kMaxClusters]  = {};
    int   cl_nblocks[kMaxClusters] = {};
    int   cl_center[kMaxClusters]  = {};

    // GEM reconstructed hits
    int     n_gem_hits = 0;
    uint8_t det_id[kMaxGemHits]       = {};
    float   gem_x[kMaxGemHits]        = {};
    float   gem_y[kMaxGemHits]        = {};
    float   gem_x_charge[kMaxGemHits] = {};
    float   gem_y_charge[kMaxGemHits] = {};
    float   gem_x_peak[kMaxGemHits]   = {};
    float   gem_y_peak[kMaxGemHits]   = {};
    int     gem_x_size[kMaxGemHits]   = {};
    int     gem_y_size[kMaxGemHits]   = {};
};

} // namespace prad2
