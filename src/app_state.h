#pragma once
//=============================================================================
// app_state.h — shared application state for evc_viewer and evc_monitor
//
// Owns all configuration, accumulated data (histograms, LMS), and HyCal system.
// Both viewer and monitor create a single AppState instance and delegate to it.
//=============================================================================

#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"
#include "viewer_utils.h"

#include <nlohmann/json.hpp>

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <mutex>
#include <atomic>

struct AppState {
    // ---- Configuration (set once at startup, then read-only) ---------------
    HistConfig hist_cfg;
    int hist_nbins = 0;
    int pos_nbins  = 0;

    evc::DaqConfig daq_cfg;
    fdec::HyCalSystem hycal;
    fdec::ClusterConfig cluster_cfg;

    std::unordered_map<int, int> roc_to_crate;  // ROC tag → crate index
    nlohmann::json crate_roc_json;              // crate→ROC tag JSON
    nlohmann::json base_config;                 // modules, daq, crate_roc for /api/config

    // LMS config
    int      lms_trigger_bit  = 16;
    float    lms_warn_thresh  = 0.1f;
    float    lms_warn_min_mean = 100.f;  // warn if mean below this
    int      lms_max_history  = 5000;
    uint32_t lms_trigger_mask = 0;

    // LMS reference channels (for normalization)
    struct LmsRefChannel {
        std::string name;
        int module_index = -1;
    };
    std::vector<LmsRefChannel> lms_ref_channels;

    // online refresh rates (ms), served to frontend
    int refresh_event_ms = 200;
    int refresh_ring_ms  = 500;
    int refresh_hist_ms  = 2000;
    int refresh_lms_ms   = 2000;

    // color range defaults: key "tab:metric" → [min, max]
    std::map<std::string, std::pair<float, float>> color_range_defaults;

    // cluster config
    uint32_t cluster_skip_mask = 0;
    float    adc_to_mev        = 1.0f;
    float    cl_hist_min       = 0.f;
    float    cl_hist_max       = 3000.f;
    float    cl_hist_step      = 10.f;
    int      nclusters_hist_min  = 0;
    int      nclusters_hist_max  = 20;
    int      nclusters_hist_step = 1;
    int      nblocks_hist_min    = 0;
    int      nblocks_hist_max    = 40;
    int      nblocks_hist_step   = 1;

    // ---- Accumulated data (guarded by data_mtx) ----------------------------
    mutable std::mutex data_mtx;
    std::map<std::string, Histogram> histograms;
    std::map<std::string, Histogram> pos_histograms;
    std::map<std::string, int>       occupancy;
    std::map<std::string, int>       occupancy_tcut;
    std::atomic<int>                 events_processed{0};

    Histogram cluster_energy_hist;
    int       cluster_events_processed = 0;

    // ---- LMS data (guarded by lms_mtx) -------------------------------------
    mutable std::mutex lms_mtx;
    std::map<int, std::vector<LmsEntry>> lms_history;
    std::atomic<int> lms_events{0};
    uint64_t lms_first_ts = 0;

    // Sync reference point for absolute time display
    // sync_unix = absolute time, sync_rel_sec = relative time on LMS axis
    uint32_t sync_unix    = 0;
    double   sync_rel_sec = 0.;

    // ---- Initialization (call once at startup) -----------------------------

    // Load all configs from db_dir. daq_config_file may be empty (PRad-II defaults).
    // config_file: main config (config.json or -c override). Empty = auto-find.
    void init(const std::string &db_dir,
              const std::string &daq_config_file,
              const std::string &config_file = "");

    // ---- Per-event processing ----------------------------------------------

    // Fill DQ histograms + occupancy for one event.
    void fillHist(fdec::EventData &event,
                  fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Run clustering on one event, fill cluster_energy_hist.
    void clusterEvent(fdec::EventData &event,
                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Process LMS data for one event (checks trigger mask internally).
    void processLms(fdec::EventData &event,
                    fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Process one fully-decoded event: histograms + clustering + LMS.
    // Thread-safe (locks internally).
    void processEvent(fdec::EventData &event,
                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Encode one decoded event as JSON (channels with waveforms, peaks, pedestal).
    nlohmann::json encodeEventJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Compute clusters for one decoded event, return JSON response.
    nlohmann::json computeClustersJson(fdec::EventData &event, int ev_id,
                                       fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Record a sync event's absolute time. Call when a Sync event is scanned.
    // last_ti_ts is the TI timestamp of the most recent physics event.
    void recordSyncTime(uint32_t unix_time, uint64_t last_ti_ts);

    // ---- Clearing ----------------------------------------------------------
    void clearHistograms();   // locks data_mtx
    void clearLms();          // locks lms_mtx

    // ---- API response builders (thread-safe) -------------------------------
    nlohmann::json apiHist(bool integral, const std::string &key) const;
    nlohmann::json apiClusterHist() const;
    nlohmann::json apiOccupancy() const;
    nlohmann::json apiColorRanges() const;
    nlohmann::json apiLmsSummary(int ref_index = -1) const;
    nlohmann::json apiLmsModule(int module_index, int ref_index = -1) const;
    nlohmann::json apiLmsRefChannels() const;
};
