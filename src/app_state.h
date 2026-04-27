#pragma once
//=============================================================================
// app_state.h — shared application state for the event viewer/monitor
//
// Owns all configuration, accumulated data (histograms, LMS), and HyCal system.
// The ViewerServer maintains two AppState instances (file and online) with
// separate accumulators but identical configuration.
//=============================================================================

// forward declaration (full definition in data_source.h)
struct ReconEventData;

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
#include <unordered_set>
#include <vector>
#include <string>
#include <mutex>
#include <atomic>
#include <iostream>

// Trigger filter: reject overrides accept. accept==0 means accept all.
struct TriggerFilter {
    uint32_t accept = 0;  // 0 = accept all
    uint32_t reject = 0;  // 0 = reject none

    bool operator()(uint32_t bits) const {
        if (reject && (bits & reject)) return false;
        return accept == 0 || (bits & accept);
    }

    // parse from JSON section containing accept_trigger_bits / reject_trigger_bits
    // values can be: bit numbers (8, 24), or names ("LMS", "Pulser") resolved
    // against the trigger_bits_def lookup table from trigger_bits.json
    void parse(const nlohmann::json &section,
               const nlohmann::json &bits_defs = nlohmann::json::array()) {
        accept = maskFrom(section, "accept_trigger_bits", bits_defs);
        reject = maskFrom(section, "reject_trigger_bits", bits_defs);
    }

    nlohmann::json toJson() const {
        return {{"trigger_accept", accept}, {"trigger_reject", reject}};
    }

    friend std::ostream &operator<<(std::ostream &os, const TriggerFilter &f) {
        return os << "trigger accept=0x" << std::hex << f.accept
                  << " reject=0x" << f.reject << std::dec;
    }

private:
    static uint32_t maskFrom(const nlohmann::json &section, const char *key,
                              const nlohmann::json &bits_defs) {
        if (!section.contains(key)) return 0;
        auto &arr = section[key];
        if (!arr.is_array() || arr.empty()) return 0;
        uint32_t m = 0;
        for (auto &item : arr) {
            if (item.is_number()) {
                m |= (1u << item.get<int>());
            } else if (item.is_string()) {
                auto s = item.get<std::string>();
                bool found = false;
                for (auto &def : bits_defs) {
                    if (def.value("name", "") == s || def.value("label", "") == s) {
                        m |= (1u << def.value("bit", 0));
                        found = true;
                        break;
                    }
                }
                if (!found)
                    std::cerr << "TriggerFilter: unknown trigger bit name '" << s << "'\n";
            }
        }
        return m;
    }
};

struct AppState {
    // ---- Configuration (set once at startup, then read-only) ---------------
    HistConfig hist_cfg;
    TriggerFilter waveform_trigger;
    int hist_nbins   = 0;
    int pos_nbins    = 0;
    int height_nbins = 0;

    // reference lines for plots (pass-through to frontend)
    // key = plot name, value = array of {axis, pos, lw, ls, color}
    nlohmann::json ref_lines = nlohmann::json::object();

    // trigger definitions (pass-through to frontend)
    nlohmann::json trigger_bits_def = nlohmann::json::array();  // FP bits: [{bit, name, label}, ...]
    nlohmann::json trigger_type_def = nlohmann::json::array();  // main trigger: [{type, tag, name, label, primary_bit}, ...]

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
    TriggerFilter lms_trigger;
    TriggerFilter alpha_trigger;          // Am-241 alpha-source trigger (LMS ref channels)
    float    lms_warn_thresh  = 0.1f;
    float    lms_warn_min_mean = 100.f;  // warn if mean below this
    int      lms_max_history  = 5000;

    // LMS reference channels (for normalization)
    struct LmsRefChannel {
        std::string name;
        int module_index = -1;
    };
    std::vector<LmsRefChannel> lms_ref_channels;

    // Latest per-channel integrals for ref correction (LMS_signal / Alpha_signal).
    // Updated unconditionally regardless of lms_max_history saturation, so the
    // correction factor always reflects the most recent readings.  Only the ref
    // module entries are read by the correction code; others are bookkeeping.
    std::map<int, float> latest_lms_integral;     // module_index → latest LMS-trigger integral
    std::map<int, float> latest_alpha_integral;   // module_index → latest Alpha-trigger integral

    // DAQ livetime polling — optional shell command (typical: "caget -t <channel>").
    // Empty disables the poller; ViewerServer reads these after init() and only
    // launches the poll thread when livetime_cmd is non-empty.  Avoids a build-
    // time EPICS dependency by shelling out to whatever tool the host provides.
    // healthy/warning are percent thresholds for the frontend display color
    // (≥ healthy → green, ≥ warning → orange, otherwise red).
    std::string livetime_cmd;
    int         livetime_poll_sec   = 30;
    float       livetime_healthy    = 90.f;
    float       livetime_warning    = 80.f;

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
    // Single-source beam energy: MBSY2C_energy from EPICS overrides; runinfo is fallback.
    // Read in physics paths (cluster filling, plots); written by init() (runinfo) and
    // processEpics() (EPICS). Atomic so processEpics can update without holding data_mtx.
    std::atomic<float> beam_energy{2200.f};                  // MeV
    float       beam_energy_runinfo = 0.f;                   // fallback (loaded from runinfo)
    std::string beam_energy_epics_channel = "MBSY2C_energy"; // EPICS channel name
    float       beam_energy_min_valid = 100.f;               // ignore stale/zero EPICS reads
    TriggerFilter physics_trigger;

    // Møller selection config
    float moller_energy_tol = 0.1f;     // energy sum within this fraction of beam_energy
    float moller_angle_min  = 1.0f;     // deg — require one cluster in this range
    float moller_angle_max  = 1.2f;     // deg
    // Møller XY histogram
    float moller_xy_x_min=-600.f, moller_xy_x_max=600.f, moller_xy_x_step=5.f;  // mm
    float moller_xy_y_min=-600.f, moller_xy_y_max=600.f, moller_xy_y_step=5.f;  // mm
    // Møller energy histogram
    float moller_e_min=0.f, moller_e_max=2500.f, moller_e_step=10.f; // MeV

    // HyCal cluster-hit XY (single-cluster ep-elastic candidates) — cuts + hist
    int   hxy_n_clusters      = 1;        // require Ncl == this
    float hxy_energy_frac_min = 0.9f;     // require E_cl >= frac * beam_energy
    int   hxy_nblocks_min     = 5;
    int   hxy_nblocks_max     = 20;
    float hxy_x_min=-600.f, hxy_x_max=600.f, hxy_x_step=5.f;  // mm
    float hxy_y_min=-600.f, hxy_y_max=600.f, hxy_y_step=5.f;  // mm

    // GEM↔HyCal matching: per-detector residuals filled when ep candidate fires
    bool  gem_match_require_ep = true;    // gate on hxy_* selection (clean track)
    float gem_match_window_mm  = 10.f;    // only fill if sqrt(dx²+dy²) < this (mm)
    float gem_resid_min = -50.f, gem_resid_max = 50.f, gem_resid_step = 0.5f;  // mm

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
    TriggerFilter cluster_trigger;
    float    adc_to_mev        = 1.0f;
    float    cl_hist_min       = 0.f;
    float    cl_hist_max       = 3000.f;
    float    cl_hist_step      = 10.f;
    // nclusters_hist range is float so the user can shift bin edges by half
    // a step (default 0.5 .. 10.5 / 1 → bin centers land on 1, 2, …, 10).
    // The bucket index that a given Ncl event falls in is reused as the
    // index into cluster_energy_hist_by_ncl / nblocks_hist_by_ncl, so the
    // dependent histograms can be filtered to "events with this many
    // clusters" by clicking a bar in the GUI.
    float    nclusters_hist_min  = 0.5f;
    float    nclusters_hist_max  = 10.5f;
    float    nclusters_hist_step = 1.0f;
    int      nblocks_hist_min    = 0;
    int      nblocks_hist_max    = 40;
    int      nblocks_hist_step   = 1;

    // ---- Event filters (loaded from external JSON via loadFilter) -----------
    // trigger_type filter: if enabled, only events with trigger_type in accept pass
    struct TriggerTypeFilter {
        bool enable = false;
        std::vector<uint8_t> accept;  // empty when disabled
        bool operator()(uint8_t tt) const {
            if (!enable || accept.empty()) return true;
            for (auto t : accept) if (t == tt) return true;
            return false;
        }
    } trigger_type_filter;

    WaveformFilter waveform_filter;
    ClusterFilter  cluster_filter;
    // resolved indices for fast filter evaluation
    std::unordered_set<std::string> filter_wf_keys;   // DAQ keys ("roc_slot_ch")
    std::unordered_set<int> filter_cl_includes;        // module indices for includes check
    std::unordered_set<int> filter_cl_centers;          // module IDs (PrimEx) for center check

    // ---- Accumulated data (guarded by data_mtx) ----------------------------
    mutable std::mutex data_mtx;
    std::map<std::string, Histogram> histograms;
    std::map<std::string, Histogram> pos_histograms;
    std::map<std::string, Histogram> height_histograms;
    std::map<std::string, int>       occupancy;
    std::map<std::string, int>       occupancy_tcut;
    std::atomic<int>                 events_processed{0};

    Histogram cluster_energy_hist;
    Histogram nclusters_hist;
    Histogram nblocks_hist;
    // Per-Ncl-bucket dependent histograms.  Index = the bucket of the
    // event's Ncl in nclusters_hist (i.e. floor((Ncl - min) / step)),
    // so size matches nclusters_hist.bins.size() once init() runs.
    // Filled alongside the unfiltered hists; the GUI selects which one
    // to plot when the user clicks an Ncl bar.
    std::vector<Histogram> cluster_energy_hist_by_ncl;
    std::vector<Histogram> nblocks_hist_by_ncl;
    Histogram2D energy_angle_hist;
    Histogram2D moller_xy_hist;
    Histogram   moller_energy_hist;
    Histogram2D hycal_xy_hist;       // single-cluster ep-elastic candidates (cuts in hxy_*)
    int         hycal_xy_events = 0; // events passing hycal_cluster_hit cuts

    // GEM↔HyCal residuals (one Histogram per GEM detector for each axis).
    // gem_match_events = events that contributed; gem_match_hits[d] = matched hits per det.
    std::vector<Histogram> gem_dx_hist;   // size = nDetectors
    std::vector<Histogram> gem_dy_hist;
    int                    gem_match_events = 0;
    std::vector<int>       gem_match_hits;  // per-det count of in-window hits
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
    uint32_t pending_sync_unix = 0;  // stashed until first LMS event arrives
    uint64_t pending_sync_ti   = 0;

    // ---- EPICS data (guarded by epics_mtx) ----------------------------------
    mutable std::mutex epics_mtx;
    fdec::EpicsStore epics;
    std::atomic<int> epics_events{0};

    // ---- GEM calibration revision -------------------------------------------
    // Bumped on any change to per-APV pedestal noise (e.g. LoadPedestals).
    // The frontend caches /api/gem/calib (which carries this rev) and
    // refetches when it sees a different value embedded in /api/gem/apv/<n>.
    // Threshold (zs_sigma) changes do NOT bump this — the frontend reads
    // zs_sigma per event so band tracks the encoded hits.
    std::atomic<int> gem_calib_rev{0};

    // ---- Initialization (call once at startup) -----------------------------

    // Load all configs from db_dir. daq_config_file may be empty (uses daq_config.json from db_dir).
    // config_file: main config (config.json or -c override). Empty = auto-find.
    void init(const std::string &db_dir,
              const std::string &daq_config_file,
              const std::string &config_file = "");

    // ---- Per-event processing ----------------------------------------------

    // Run the GEM clear + ProcessEvent + Reconstruct pipeline without
    // touching accumulators.  Used by the on-demand APV endpoint so a
    // re-request for an event already accumulated still leaves gem_sys
    // populated with that event's per-APV working buffers.
    void prepareGemForView(const ssp::SspEventData &ssp_evt);

    // Process GEM SSP data for one event. Call after DecodeEvent with ssp_evt.
    // Calls prepareGemForView, then accumulates occupancy + histograms.
    void processGemEvent(const ssp::SspEventData &ssp_evt);

    // Process one fully-decoded event: histograms + clustering + LMS.
    // Single pass over all channels (analyzes each channel once).
    // Thread-safe (acquires data_mtx + lms_mtx internally).
    void processEvent(fdec::EventData &event,
                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Process a pre-computed recon event (from ROOT recon files).
    // Fills cluster/physics histograms from pre-computed clusters.
    void processReconEvent(const struct ReconEventData &recon);

    // Project a lab-frame point through the target onto the HyCal local
    // plane (z_local = 0). Returns HyCal-local (px, py).
    // Used for GEM hit overlays and matching residuals.
    void projectToHyCalLocal(float Gx, float Gy, float Gz,
                             float &px, float &py) const;

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

    // Encode pre-computed recon clusters (from ROOT recon files) as JSON.
    nlohmann::json encodeReconClustersJson(const struct ReconEventData &recon, int ev_id);

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
    // type: 0=integral, 1=position, 2=height
    nlohmann::json apiHist(int type, const std::string &key) const;
    nlohmann::json apiClusterHist() const;
    nlohmann::json apiOccupancy() const;
    nlohmann::json apiColorRanges() const;
    nlohmann::json apiLmsSummary(int ref_index = -1) const;
    nlohmann::json apiLmsModule(int module_index, int ref_index = -1) const;
    nlohmann::json apiLmsRefChannels() const;
    nlohmann::json apiEnergyAngle() const;
    nlohmann::json apiMoller() const;
    nlohmann::json apiHycalXY() const;
    nlohmann::json apiGemResiduals() const;
    nlohmann::json apiEpicsChannels() const;
    nlohmann::json apiEpicsChannel(const std::string &name) const;
    nlohmann::json apiEpicsBatch(const std::vector<std::string> &names) const;
    nlohmann::json apiEpicsLatest() const;
    nlohmann::json apiGemHits() const;
    nlohmann::json apiGemConfig() const;
    nlohmann::json apiGemOccupancy() const;
    nlohmann::json apiGemHist() const;
    // Per-event APV waveform dump (for the GEM APV monitor tab).
    // Caller must have just populated gem_sys with this event (e.g. via
    // processGemEvent or accumulate); this method only reads.  Raw ADC
    // samples come from ssp_evt; processed values + ZS hit mask come
    // from gem_sys's per-APV working buffer.
    nlohmann::json apiGemApv(const ssp::SspEventData &ssp_evt, int evnum) const;

    // One-shot calibration payload for the GEM APV tab: returns
    //   {rev, zs_sigma, apvs:[{id, noise:[128]}, ...]}
    // The frontend caches this and refetches only when the calib_rev
    // embedded in /api/gem/apv/<n> diverges from the cached value.
    nlohmann::json apiGemCalib() const;

    // Update the software N-sigma cut on this AppState's gem_sys.  Does
    // NOT bump gem_calib_rev — noise is unchanged, and frontend uses the
    // per-event zs_sigma so band always tracks encoded hits.
    void setGemZsSigma(float v);

    // ---- Filters ---------------------------------------------------------------
    std::string loadFilter(const nlohmann::json &j);
    void unloadFilter();
    nlohmann::json filterToJson() const;
    bool filterActive() const { return trigger_type_filter.enable || waveform_filter.enable || cluster_filter.enable; }
    bool evaluateFilter(fdec::EventData &event, ssp::SspEventData *ssp) const;

    // Fill common config fields into a JSON object (used by both viewer and monitor).
    void fillConfigJson(nlohmann::json &cfg) const;

    // Handle a read-only API route. Returns {handled, response_json}.
    // Does NOT handle /api/config, clear endpoints, or mode-specific routes.
    struct ApiResult { bool handled; std::string body; };
    ApiResult handleReadApi(const std::string &uri) const;

private:
    void resolveFilterKeys();
};
