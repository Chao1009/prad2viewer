#pragma once
//=============================================================================
// app_state.h — shared application state for the event viewer/monitor
//
// Owns all configuration, accumulated data (histograms, LMS), and HyCal system.
// The ViewerServer maintains two AppState instances (file and online) with
// separate accumulators but identical configuration.
//=============================================================================

#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "EpicsStore.h"
#include "DaqConfig.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "GemSystem.h"
#include "GemCluster.h"
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
    uint32_t waveform_trigger_mask = 0;  // 0 = accept all
    int hist_nbins = 0;
    int pos_nbins  = 0;

    evc::DaqConfig daq_cfg;
    fdec::HyCalSystem hycal;
    fdec::ClusterConfig cluster_cfg;

    // GEM system
    gem::GemSystem gem_sys;
    gem::GemCluster gem_clusterer;
    bool gem_enabled = false;       // true if gem_map.json loaded successfully

    // GEM per-detector lab-frame transform (same type as HyCal)
    std::vector<DetectorTransform> gem_transforms;  // indexed by detector id

    // GEM occupancy (accumulated per-detector 2D histograms)
    static constexpr int GEM_OCC_NX = 50;
    static constexpr int GEM_OCC_NY = 30;
    std::vector<Histogram2D> gem_occupancy;  // one per detector

    // GEM accumulated histograms
    int   gem_ncl_min=0, gem_ncl_max=50, gem_ncl_step=1;
    float gem_theta_min=0.f, gem_theta_max=10.f, gem_theta_step=0.2f;
    Histogram gem_nclusters_hist;
    Histogram gem_theta_hist;

    std::unordered_map<int, int> roc_to_crate;  // ROC tag → crate index
    nlohmann::json crate_roc_json;              // crate→ROC tag JSON
    nlohmann::json base_config;                 // modules, daq, crate_roc for /api/config

    // LMS config
    uint32_t lms_trigger_mask = 0;       // 0 = accept all
    float    lms_warn_thresh  = 0.1f;
    float    lms_warn_min_mean = 100.f;  // warn if mean below this
    int      lms_max_history  = 5000;

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

    // Physics / coordinate config
    float target_x=0, target_y=0, target_z=0;  // target position in lab frame (mm)
    DetectorTransform hycal_transform;           // HyCal position + tilting
    float ea_angle_min=0.f, ea_angle_max=8.f, ea_angle_step=0.2f;   // degrees
    float ea_energy_min=0.f, ea_energy_max=3000.f, ea_energy_step=100.f; // MeV
    float beam_energy = 2200.f;  // MeV (for elastic line overlay)
    uint32_t physics_trigger_mask = 0;  // 0 = accept all

    // Møller selection config
    float moller_energy_tol = 0.1f;     // energy sum within this fraction of beam_energy
    float moller_angle_min  = 1.0f;     // deg — require one cluster in this range
    float moller_angle_max  = 1.2f;     // deg
    // Møller XY histogram
    float moller_xy_x_min=-600.f, moller_xy_x_max=600.f, moller_xy_x_step=5.f;  // mm
    float moller_xy_y_min=-600.f, moller_xy_y_max=600.f, moller_xy_y_step=5.f;  // mm
    // Møller energy histogram
    float moller_e_min=0.f, moller_e_max=2500.f, moller_e_step=10.f; // MeV

    // EPICS config
    int   epics_max_history = 5000;
    float epics_warn_thresh  = 0.1f;
    float epics_alert_thresh = 0.2f;
    int   epics_min_avg_pts  = 10;
    int   epics_mean_window  = 20;   // compute mean from most recent N snapshots
    std::vector<std::vector<std::string>> epics_default_slots;  // per-slot channel lists

    // Elog config
    std::string elog_url;
    std::string elog_logbook;
    std::string elog_author;
    std::vector<std::string> elog_tags;
    std::string elog_cert;         // SSL client certificate path
    std::string elog_key;          // SSL client key path

    // color range defaults: key "tab:metric" → [min, max]
    std::map<std::string, std::pair<float, float>> color_range_defaults;

    // cluster config
    uint32_t cluster_trigger_mask = 0;   // 0 = accept all
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
    Histogram nclusters_hist;
    Histogram nblocks_hist;
    Histogram2D energy_angle_hist;
    Histogram2D moller_xy_hist;
    Histogram   moller_energy_hist;
    int         moller_events = 0;
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

    // ---- EPICS data (guarded by epics_mtx) ----------------------------------
    mutable std::mutex epics_mtx;
    fdec::EpicsStore epics;
    std::atomic<int> epics_events{0};

    // ---- Initialization (call once at startup) -----------------------------

    // Load all configs from db_dir. daq_config_file may be empty (uses daq_config.json from db_dir).
    // config_file: main config (config.json or -c override). Empty = auto-find.
    void init(const std::string &db_dir,
              const std::string &daq_config_file,
              const std::string &config_file = "");

    // ---- Per-event processing ----------------------------------------------

    // Process GEM SSP data for one event. Call after DecodeEvent with ssp_evt.
    void processGemEvent(const ssp::SspEventData &ssp_evt);

    // Process one fully-decoded event: histograms + clustering + LMS.
    // Single pass over all channels (analyzes each channel once).
    // Thread-safe (acquires data_mtx + lms_mtx internally).
    void processEvent(fdec::EventData &event,
                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Encode one decoded event as JSON.
    // include_samples=false: summaries only (peaks + pedestal, ~20KB).
    // include_samples=true:  full waveforms included (~800KB, used for ring buffer).
    nlohmann::json encodeEventJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres,
                                   bool include_samples = false);

    // Encode raw waveform for a single channel (key = "roc_slot_ch").
    nlohmann::json encodeWaveformJson(fdec::EventData &event, const std::string &chan_key,
                                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Compute clusters for one decoded event, return JSON response.
    nlohmann::json computeClustersJson(fdec::EventData &event, int ev_id,
                                       fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Record a sync event's absolute time. Call when a Sync event is scanned.
    // last_ti_ts is the TI timestamp of the most recent physics event.
    void recordSyncTime(uint32_t unix_time, uint64_t last_ti_ts);

    // ---- EPICS processing ---------------------------------------------------
    void processEpics(const std::string &text, int32_t event_number, uint64_t timestamp);
    void clearEpics();        // locks epics_mtx

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
    nlohmann::json apiEnergyAngle() const;
    nlohmann::json apiMoller() const;
    nlohmann::json apiEpicsChannels() const;
    nlohmann::json apiEpicsChannel(const std::string &name) const;
    nlohmann::json apiEpicsBatch(const std::vector<std::string> &names) const;
    nlohmann::json apiEpicsLatest() const;
    nlohmann::json apiGemHits() const;
    nlohmann::json apiGemConfig() const;
    nlohmann::json apiGemOccupancy() const;
    nlohmann::json apiGemHist() const;

    // Fill common config fields into a JSON object (used by both viewer and monitor).
    void fillConfigJson(nlohmann::json &cfg) const;

    // Handle a read-only API route. Returns {handled, response_json}.
    // Does NOT handle /api/config, clear endpoints, or mode-specific routes.
    struct ApiResult { bool handled; std::string body; };
    ApiResult handleReadApi(const std::string &uri) const;
};
