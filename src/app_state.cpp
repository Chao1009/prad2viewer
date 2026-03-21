#include "app_state.h"
#include "load_daq_config.h"

#include <fstream>
#include <iostream>
#include <cmath>

using json = nlohmann::json;

//=============================================================================
// Initialization
//=============================================================================

void AppState::init(const std::string &db_dir,
                    const std::string &daq_config_file,
                    const std::string &config_file)
{
    // resolve main config file: -c override > config.json > reconstruction.json
    std::string main_config = config_file;
    if (main_config.empty()) main_config = findFile("config.json", db_dir);
    if (main_config.empty()) main_config = findFile("reconstruction.json", db_dir);

    // --- DAQ config + pedestals ---
    if (!daq_config_file.empty()) {
        if (evc::load_daq_config(daq_config_file, daq_cfg)) {
            std::cerr << "DAQ config: " << daq_config_file
                      << " (adc_format=" << daq_cfg.adc_format << ")\n";
            std::ifstream dcf(daq_config_file);
            if (dcf.is_open()) {
                auto dcj = json::parse(dcf, nullptr, false, true);
                if (dcj.contains("pedestal_file")) {
                    std::string ped_file = findFile(dcj["pedestal_file"].get<std::string>(), db_dir);
                    if (evc::load_pedestals(ped_file, daq_cfg))
                        std::cerr << "Pedestals : " << ped_file
                                  << " (" << daq_cfg.pedestals.size() << " channels)\n";
                }
            }
        } else {
            std::cerr << "Warning: failed to load DAQ config: " << daq_config_file << "\n";
        }
    }

    // --- Waveform / histogram config ---
    // Try "waveform" section from config.json (loaded later in reco_file),
    // then fall back to separate hist_config.json, then -H override.
    auto loadWaveformConfig = [&](const json &w) {
        if (w.contains("time_cut")) {
            auto &tc = w["time_cut"];
            if (tc.contains("min")) hist_cfg.time_min = tc["min"];
            if (tc.contains("max")) hist_cfg.time_max = tc["max"];
        }
        if (w.contains("integral_hist")) {
            auto &ih = w["integral_hist"];
            if (ih.contains("min"))  hist_cfg.bin_min  = ih["min"];
            if (ih.contains("max"))  hist_cfg.bin_max  = ih["max"];
            if (ih.contains("step")) hist_cfg.bin_step = ih["step"];
        }
        if (w.contains("time_hist")) {
            auto &th = w["time_hist"];
            if (th.contains("min"))  hist_cfg.pos_min  = th["min"];
            if (th.contains("max"))  hist_cfg.pos_max  = th["max"];
            if (th.contains("step")) hist_cfg.pos_step = th["step"];
        }
        if (w.contains("thresholds")) {
            auto &t = w["thresholds"];
            if (t.contains("min_peak_height"))          hist_cfg.threshold      = t["min_peak_height"];
            if (t.contains("min_secondary_peak_ratio")) hist_cfg.min_peak_ratio = t["min_secondary_peak_ratio"];
        }
    };

    // legacy: load from separate hist_config.json or -H override
    auto loadLegacyHistConfig = [&](const json &hcfg) {
        if (hcfg.contains("hist")) {
            auto &h = hcfg["hist"];
            if (h.contains("time_min"))  hist_cfg.time_min  = h["time_min"];
            if (h.contains("time_max"))  hist_cfg.time_max  = h["time_max"];
            if (h.contains("bin_min"))   hist_cfg.bin_min   = h["bin_min"];
            if (h.contains("bin_max"))   hist_cfg.bin_max   = h["bin_max"];
            if (h.contains("bin_step"))  hist_cfg.bin_step  = h["bin_step"];
            if (h.contains("threshold")) hist_cfg.threshold = h["threshold"];
            if (h.contains("pos_min"))   hist_cfg.pos_min   = h["pos_min"];
            if (h.contains("pos_max"))   hist_cfg.pos_max   = h["pos_max"];
            if (h.contains("pos_step"))  hist_cfg.pos_step  = h["pos_step"];
            if (h.contains("min_peak_ratio")) hist_cfg.min_peak_ratio = h["min_peak_ratio"];
        }
    };

    // load waveform config from main config, or legacy hist_config.json
    bool waveform_loaded = false;
    if (!main_config.empty()) {
        std::string s = readFile(main_config);
        if (!s.empty()) {
            auto j = json::parse(s, nullptr, false);
            if (j.contains("waveform")) { loadWaveformConfig(j["waveform"]); waveform_loaded = true; }
            // legacy "hist" section
            else if (j.contains("hist")) { loadLegacyHistConfig(j); waveform_loaded = true; }
        }
    }
    if (!waveform_loaded) {
        std::string hcfg_path = findFile("hist_config.json", db_dir);
        std::string hcfg_str = readFile(hcfg_path);
        if (!hcfg_str.empty())
            loadLegacyHistConfig(json::parse(hcfg_str, nullptr, false));
    }

    hist_nbins = std::max(1, (int)std::ceil(
        (hist_cfg.bin_max - hist_cfg.bin_min) / hist_cfg.bin_step));
    pos_nbins = std::max(1, (int)std::ceil(
        (hist_cfg.pos_max - hist_cfg.pos_min) / hist_cfg.pos_step));
    std::cerr << "Waveform  : time_cut=[" << hist_cfg.time_min << "," << hist_cfg.time_max
              << "] threshold=" << hist_cfg.threshold << "\n";

    // --- HyCal system ---
    std::string modules_filename = "hycal_modules.json";
    std::string daq_filename     = "daq_map.json";
    if (!daq_config_file.empty()) {
        std::ifstream dcf2(daq_config_file);
        if (dcf2.is_open()) {
            auto dcj2 = json::parse(dcf2, nullptr, false, true);
            if (dcj2.contains("modules_file")) modules_filename = dcj2["modules_file"].get<std::string>();
            if (dcj2.contains("daq_map_file")) daq_filename = dcj2["daq_map_file"].get<std::string>();
        }
    }
    std::string mod_file = findFile(modules_filename, db_dir);
    std::string daq_file = findFile(daq_filename, db_dir);

    if (!mod_file.empty() && !daq_file.empty()) {
        if (hycal.Init(mod_file, daq_file))
            std::cerr << "HyCal     : " << hycal.module_count() << " modules\n";
        else
            std::cerr << "Warning: HyCal system initialization failed\n";
    }

    // --- crate_roc map ---
    crate_roc_json = json::object();
    if (!daq_config_file.empty()) {
        std::ifstream dcf3(daq_config_file);
        if (dcf3.is_open()) {
            auto dcj3 = json::parse(dcf3, nullptr, false, true);
            if (dcj3.contains("roc_tags")) {
                for (auto &entry : dcj3["roc_tags"]) {
                    if (entry.contains("crate") && entry.contains("tag")) {
                        int crate = entry["crate"].get<int>();
                        uint32_t tag = evc::parse_hex(entry["tag"]);
                        crate_roc_json[std::to_string(crate)] = tag;
                    }
                }
            }
        }
    }
    if (crate_roc_json.empty())
        crate_roc_json = {{"0",0x80},{"1",0x82},{"2",0x84},{"3",0x86},{"4",0x88},{"5",0x8a},{"6",0x8c}};

    roc_to_crate.clear();
    for (auto &[k, v] : crate_roc_json.items())
        roc_to_crate[v.get<int>()] = std::stoi(k);

    // --- Reconstruction config (from same main config file) ---
    std::string reco_str = readFile(main_config);
    if (!reco_str.empty()) {
        auto rcfg = json::parse(reco_str, nullptr, false);

        if (rcfg.contains("clustering")) {
            auto &cc = rcfg["clustering"];
            auto loadCfg = [](const json &j, fdec::ClusterConfig &cfg) {
                if (j.contains("min_module_energy"))  cfg.min_module_energy  = j["min_module_energy"];
                if (j.contains("min_center_energy"))  cfg.min_center_energy  = j["min_center_energy"];
                if (j.contains("min_cluster_energy")) cfg.min_cluster_energy = j["min_cluster_energy"];
                if (j.contains("min_cluster_size"))   cfg.min_cluster_size   = j["min_cluster_size"];
                if (j.contains("corner_conn"))        cfg.corner_conn        = j["corner_conn"];
                if (j.contains("split_iter"))         cfg.split_iter         = j["split_iter"];
                if (j.contains("least_split"))        cfg.least_split        = j["least_split"];
                if (j.contains("log_weight_thres"))   cfg.log_weight_thres   = j["log_weight_thres"];
            };
            loadCfg(cc, cluster_cfg);
            if (cc.contains("skip_trigger_bits")) {
                cluster_skip_mask = 0;
                for (auto &b : cc["skip_trigger_bits"])
                    cluster_skip_mask |= (1u << b.get<int>());
            }
            if (cc.contains("energy_hist")) {
                auto &eh = cc["energy_hist"];
                if (eh.contains("min"))  cl_hist_min  = eh["min"];
                if (eh.contains("max"))  cl_hist_max  = eh["max"];
                if (eh.contains("step")) cl_hist_step = eh["step"];
            }
            std::cerr << "Clustering: min_mod=" << cluster_cfg.min_module_energy
                      << " min_center=" << cluster_cfg.min_center_energy
                      << " min_cluster=" << cluster_cfg.min_cluster_energy
                      << " skip_mask=0x" << std::hex << cluster_skip_mask << std::dec
                      << " hist=[" << cl_hist_min << "," << cl_hist_max
                      << "]/" << cl_hist_step << "\n";
        }

        if (rcfg.contains("lms_monitor")) {
            auto &lm = rcfg["lms_monitor"];
            if (lm.contains("trigger_bit"))    lms_trigger_bit   = lm["trigger_bit"];
            if (lm.contains("warn_threshold")) lms_warn_thresh     = lm["warn_threshold"];
            if (lm.contains("warn_min_mean"))  lms_warn_min_mean  = lm["warn_min_mean"];
            if (lm.contains("max_history"))    lms_max_history    = lm["max_history"];
            lms_trigger_mask = (1u << lms_trigger_bit);
            if (lm.contains("reference_channels")) {
                for (auto &name : lm["reference_channels"]) {
                    std::string n = name.get<std::string>();
                    const auto *mod = hycal.module_by_name(n);
                    lms_ref_channels.push_back({n, mod ? mod->index : -1});
                }
            }
            std::cerr << "LMS       : trigger_bit=" << lms_trigger_bit
                      << " mask=0x" << std::hex << lms_trigger_mask << std::dec
                      << " warn=" << lms_warn_thresh
                      << " refs=" << lms_ref_channels.size() << "\n";
        }

        if (rcfg.contains("color_ranges")) {
            for (auto &[key, val] : rcfg["color_ranges"].items()) {
                if (val.is_array() && val.size() == 2)
                    color_range_defaults[key] = {val[0].get<float>(), val[1].get<float>()};
            }
            std::cerr << "Color ranges: " << color_range_defaults.size() << " entries\n";
        }

        if (rcfg.contains("calibration")) {
            auto &cal = rcfg["calibration"];
            if (cal.contains("adc_to_mev")) adc_to_mev = cal["adc_to_mev"];
            if (cal.contains("calibration_file")) {
                std::string calib_file = findFile(cal["calibration_file"].get<std::string>(), db_dir);
                int nmatched = hycal.LoadCalibration(calib_file);
                if (nmatched >= 0)
                    std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";
            }
        }
        std::cerr << "Reco      : " << main_config
                  << " (adc_to_mev=" << adc_to_mev << ")\n";
    }

    // init cluster energy histogram
    int cl_nbins = std::max(1, (int)std::ceil((cl_hist_max - cl_hist_min) / cl_hist_step));
    cluster_energy_hist.init(cl_nbins);
}

//=============================================================================
// Per-event processing
//=============================================================================

void AppState::fillHist(fdec::EventData &event,
                        fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.Analyze(cd.samples, cd.nsamples, wres);

                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                bool has_peak = false, has_peak_tcut = false;
                float best = -1;
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    if (pk.height < hist_cfg.threshold) continue;
                    has_peak = true;
                    if (pk.time >= hist_cfg.time_min && pk.time <= hist_cfg.time_max) {
                        has_peak_tcut = true;
                        if (pk.integral > best) best = pk.integral;
                    }
                }
                if (best >= 0) {
                    auto &h = histograms[key];
                    if (h.bins.empty()) h.init(hist_nbins);
                    h.fill(best, hist_cfg.bin_min, hist_cfg.bin_step);
                }
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    if (pk.height < hist_cfg.threshold) continue;
                    auto &ph = pos_histograms[key];
                    if (ph.bins.empty()) ph.init(pos_nbins);
                    ph.fill(pk.time, hist_cfg.pos_min, hist_cfg.pos_step);
                }
                if (has_peak)      occupancy[key]++;
                if (has_peak_tcut) occupancy_tcut[key]++;
            }
        }
    }
}

void AppState::clusterEvent(fdec::EventData &event,
                            fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    if (cluster_skip_mask != 0 &&
        (event.info.trigger_bits & cluster_skip_mask)) return;

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);

    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        auto cit = roc_to_crate.find(roc.tag);
        if (cit == roc_to_crate.end()) continue;
        int crate = cit->second;

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                const auto *mod = hycal.module_by_daq(crate, s, c);
                if (!mod || !mod->is_hycal()) continue;

                float adc_val = 0;
                if (is_adc1881m) {
                    adc_val = cd.samples[0];
                } else {
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    adc_val = bestPeakInWindow(wres, hist_cfg.threshold,
                                               hist_cfg.time_min, hist_cfg.time_max);
                }
                if (adc_val <= 0) continue;

                float energy = (mod->cal_factor > 0.)
                    ? static_cast<float>(mod->energize(adc_val))
                    : adc_val * adc_to_mev;
                clusterer.AddHit(mod->index, energy);
            }
        }
    }

    clusterer.FormClusters();
    std::vector<fdec::ClusterHit> reco_hits;
    clusterer.ReconstructHits(reco_hits);
    for (auto &rh : reco_hits)
        cluster_energy_hist.fill(rh.energy, cl_hist_min, cl_hist_step);
    cluster_events_processed++;
}

void AppState::processLms(fdec::EventData &event,
                          fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // debug: log first few trigger bits seen
    static int dbg_count = 0;
    if (dbg_count < 20) {
        std::cerr << "LMS check: trig_bits=0x" << std::hex << event.info.trigger_bits
                  << " mask=0x" << lms_trigger_mask << std::dec
                  << " match=" << ((event.info.trigger_bits & lms_trigger_mask) ? "YES" : "NO")
                  << "\n";
        dbg_count++;
    }

    if (lms_trigger_mask == 0 ||
        !(event.info.trigger_bits & lms_trigger_mask)) return;

    if (lms_first_ts == 0)
        lms_first_ts = event.info.timestamp;
    double time_sec = static_cast<double>(event.info.timestamp - lms_first_ts) * TI_TICK_SEC;

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");

    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        auto cit = roc_to_crate.find(roc.tag);
        if (cit == roc_to_crate.end()) continue;
        int crate = cit->second;

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                const auto *mod = hycal.module_by_daq(crate, s, c);
                if (!mod) continue;  // include LMS modules

                float val = 0;
                if (is_adc1881m) {
                    val = cd.samples[0];
                } else {
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    val = bestPeakInWindow(wres, hist_cfg.threshold,
                                           hist_cfg.time_min, hist_cfg.time_max);
                }
                if (val <= 0) continue;

                auto &hist = lms_history[mod->index];
                if (static_cast<int>(hist.size()) < lms_max_history)
                    hist.push_back({time_sec, val});
            }
        }
    }
    lms_events++;
}

json AppState::encodeEventJson(fdec::EventData &event, int ev_id,
                               fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    json channels = json::object();
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.Analyze(cd.samples, cd.nsamples, wres);
                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                json sarr = json::array();
                for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);

                json parr = json::array();
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    parr.push_back({
                        {"p", pk.pos}, {"t", std::round(pk.time * 10) / 10},
                        {"h", std::round(pk.height * 10) / 10},
                        {"i", std::round(pk.integral * 10) / 10},
                        {"l", pk.left}, {"r", pk.right},
                        {"o", pk.overflow ? 1 : 0},
                    });
                }
                channels[key] = {
                    {"s", sarr},
                    {"pm", std::round(wres.ped.mean * 10) / 10},
                    {"pr", std::round(wres.ped.rms * 10) / 10},
                    {"pk", parr},
                };
            }
        }
    }
    return {{"event", ev_id}, {"channels", channels},
            {"event_number", event.info.event_number},
            {"trigger_bits", event.info.trigger_bits}};
}

json AppState::computeClustersJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    if (cluster_skip_mask != 0 &&
        (event.info.trigger_bits & cluster_skip_mask))
        return {{"event", ev_id}, {"hits", json::object()}, {"clusters", json::array()},
                {"info", "trigger filtered"}};

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);

    int nmod = hycal.module_count();
    std::vector<float> mod_energy(nmod, 0.f);

    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        auto cit = roc_to_crate.find(roc.tag);
        if (cit == roc_to_crate.end()) continue;
        int crate = cit->second;

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                const auto *mod = hycal.module_by_daq(crate, s, c);
                if (!mod || !mod->is_hycal()) continue;

                float adc_val = 0;
                if (is_adc1881m) {
                    adc_val = cd.samples[0];
                } else {
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    adc_val = bestPeakInWindow(wres, hist_cfg.threshold,
                                               hist_cfg.time_min, hist_cfg.time_max);
                }
                if (adc_val <= 0) continue;

                float energy = (mod->cal_factor > 0.)
                    ? static_cast<float>(mod->energize(adc_val))
                    : adc_val * adc_to_mev;
                mod_energy[mod->index] = energy;
                clusterer.AddHit(mod->index, energy);
            }
        }
    }

    clusterer.FormClusters();

    json hits_j = json::object();
    for (int i = 0; i < nmod; ++i)
        if (mod_energy[i] > 0.f)
            hits_j[std::to_string(i)] = std::round(mod_energy[i] * 100) / 100;

    std::vector<fdec::HyCalCluster::RecoResult> reco;
    clusterer.ReconstructMatched(reco);

    json cl_arr = json::array();
    for (auto &r : reco) {
        auto &cmod = hycal.module(r.cluster->center.index);
        json indices = json::array();
        for (auto &h : r.cluster->hits) indices.push_back(h.index);
        cl_arr.push_back({
            {"id", static_cast<int>(cl_arr.size())},
            {"center", cmod.name}, {"center_id", cmod.id},
            {"x", std::round(r.hit.x * 10) / 10},
            {"y", std::round(r.hit.y * 10) / 10},
            {"energy", std::round(r.hit.energy * 10) / 10},
            {"nblocks", r.hit.nblocks}, {"npos", r.hit.npos},
            {"modules", indices},
        });
    }

    return {{"event", ev_id}, {"hits", hits_j}, {"clusters", cl_arr}};
}

void AppState::recordSyncTime(uint32_t unix_time, uint64_t last_ti_ts)
{
    if (unix_time == 0) return;
    std::lock_guard<std::mutex> lk(lms_mtx);
    if (sync_unix != 0) return;   // only record first sync after LMS starts
    if (lms_first_ts == 0) return; // no LMS events yet — ignore this sync
    if (last_ti_ts == 0) return;   // no TI reference available
    sync_unix = unix_time;
    sync_rel_sec = static_cast<double>(last_ti_ts - lms_first_ts) * TI_TICK_SEC;
}

void AppState::processEvent(fdec::EventData &event,
                            fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    {
        std::lock_guard<std::mutex> lk(data_mtx);
        fillHist(event, ana, wres);
        clusterEvent(event, ana, wres);
        events_processed++;
    }
    {
        std::lock_guard<std::mutex> lk(lms_mtx);
        processLms(event, ana, wres);
    }
}

//=============================================================================
// Clearing
//=============================================================================

void AppState::clearHistograms()
{
    std::lock_guard<std::mutex> lk(data_mtx);
    for (auto &[k, h] : histograms)     h.clear();
    for (auto &[k, h] : pos_histograms) h.clear();
    occupancy.clear();
    occupancy_tcut.clear();
    events_processed = 0;
    cluster_energy_hist.clear();
    cluster_events_processed = 0;
}

void AppState::clearLms()
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    lms_history.clear();
    lms_events = 0;
    lms_first_ts = 0;
    sync_unix = 0;
    sync_rel_sec = 0.;
}

//=============================================================================
// API response builders
//=============================================================================

json AppState::apiColorRanges() const
{
    json obj = json::object();
    for (auto &[k, v] : color_range_defaults)
        obj[k] = {v.first, v.second};
    return obj;
}

json AppState::apiHist(bool integral, const std::string &key) const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    auto &hmap = integral ? histograms : pos_histograms;
    int nbins = integral ? hist_nbins : pos_nbins;
    auto it = hmap.find(key);
    if (it == hmap.end())
        return {{"bins", std::vector<int>(nbins, 0)}, {"underflow", 0}, {"overflow", 0},
                {"events", events_processed.load()}};
    auto &h = it->second;
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"events", events_processed.load()}};
}

json AppState::apiClusterHist() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    if (cluster_energy_hist.bins.empty())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"events", 0}, {"min", cl_hist_min}, {"max", cl_hist_max},
                {"step", cl_hist_step}};
    auto &h = cluster_energy_hist;
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"events", cluster_events_processed},
            {"min", cl_hist_min}, {"max", cl_hist_max}, {"step", cl_hist_step}};
}

json AppState::apiOccupancy() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    json jocc = json::object(), jtcut = json::object();
    for (auto &[k,v] : occupancy) jocc[k] = v;
    for (auto &[k,v] : occupancy_tcut) jtcut[k] = v;
    return {{"occ", jocc}, {"occ_tcut", jtcut}, {"total", events_processed.load()}};
}

// Reference correction: builds time→value map and computes mean for correction factor.
// Correction: corrected = signal * (ref_mean / ref_signal_at_time)
// This removes LMS-own fluctuation while keeping values in original units.
struct RefCorrection {
    std::map<double, float> ref_map;  // time → ref signal
    float ref_mean = 0.f;             // mean of all ref signals
    bool active = false;
};

static RefCorrection buildRefCorrection(
    const std::map<int, std::vector<LmsEntry>> &lms_history,
    const std::vector<AppState::LmsRefChannel> &refs, int ref_index)
{
    RefCorrection rc;
    if (ref_index < 0 || ref_index >= static_cast<int>(refs.size())) return rc;
    int ri = refs[ref_index].module_index;
    if (ri < 0) return rc;
    auto it = lms_history.find(ri);
    if (it == lms_history.end() || it->second.empty()) return rc;

    double sum = 0;
    for (auto &e : it->second) {
        rc.ref_map[e.time_sec] = e.integral;
        sum += e.integral;
    }
    rc.ref_mean = static_cast<float>(sum / it->second.size());
    rc.active = (rc.ref_mean > 0);
    return rc;
}

// Apply correction: returns signal * (ref_mean / ref_at_time), or -1 if ref missing.
static float applyRefCorrection(float val, double time_sec, const RefCorrection &rc)
{
    if (!rc.active) return val;
    auto it = rc.ref_map.find(time_sec);
    if (it == rc.ref_map.end() || it->second <= 0) return -1.f;
    return val * (rc.ref_mean / it->second);
}

json AppState::apiLmsSummary(int ref_index) const
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    auto rc = buildRefCorrection(lms_history, lms_ref_channels, ref_index);

    json mods = json::object();
    for (auto &[idx, hist] : lms_history) {
        if (hist.empty()) continue;
        double sum = 0, sum2 = 0;
        int count = 0;
        for (auto &e : hist) {
            float v = applyRefCorrection(e.integral, e.time_sec, rc);
            if (v < 0) continue;
            sum += v; sum2 += v * v;
            count++;
        }
        if (count == 0) continue;
        double mean = sum / count;
        double var = sum2 / count - mean * mean;
        double rms = var > 0 ? std::sqrt(var) : 0;
        bool warn = (mean > 0 && rms / mean > lms_warn_thresh) ||
                    (mean < lms_warn_min_mean);
        if (idx >= 0 && idx < hycal.module_count()) {
            auto &mod = hycal.module(idx);
            mods[std::to_string(idx)] = {
                {"name", mod.name}, {"mean", std::round(mean * 10) / 10},
                {"rms", std::round(rms * 100) / 100},
                {"count", count}, {"warn", warn}};
        }
    }
    return {{"modules", mods}, {"events", lms_events.load()},
            {"trigger_bit", lms_trigger_bit},
            {"ref_index", ref_index},
            {"ref_mean", rc.ref_mean},
            {"sync_unix", sync_unix}, {"sync_rel_sec", sync_rel_sec}};
}

json AppState::apiLmsModule(int mod_idx, int ref_index) const
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    auto it = lms_history.find(mod_idx);
    if (it == lms_history.end() || it->second.empty())
        return {{"time", json::array()}, {"integral", json::array()}, {"events", 0}};

    auto rc = buildRefCorrection(lms_history, lms_ref_channels, ref_index);

    auto &hist = it->second;
    json t_arr = json::array(), v_arr = json::array();
    for (auto &e : hist) {
        float v = applyRefCorrection(e.integral, e.time_sec, rc);
        if (v < 0) continue;
        t_arr.push_back(std::round(e.time_sec * 100) / 100);
        v_arr.push_back(std::round(v * 10) / 10);
    }
    std::string name = (mod_idx >= 0 && mod_idx < hycal.module_count())
        ? hycal.module(mod_idx).name : "";
    return {{"name", name}, {"time", t_arr}, {"integral", v_arr},
            {"events", (int)t_arr.size()},
            {"ref_index", ref_index},
            {"sync_unix", sync_unix}, {"sync_rel_sec", sync_rel_sec}};
}

json AppState::apiLmsRefChannels() const
{
    json arr = json::array();
    for (size_t i = 0; i < lms_ref_channels.size(); ++i) {
        arr.push_back({
            {"index", (int)i},
            {"name", lms_ref_channels[i].name},
            {"module_index", lms_ref_channels[i].module_index},
        });
    }
    return arr;
}
