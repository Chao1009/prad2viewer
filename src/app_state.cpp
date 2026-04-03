#include "app_state.h"
#include "data_source.h"
#include "load_daq_config.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

using json = nlohmann::json;

// Serialize a Histogram to JSON.
static json histToJson(const Histogram &h, float mn, float mx, float st)
{
    if (h.bins.empty())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"min", mn}, {"max", mx}, {"step", st}};
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"min", mn}, {"max", mx}, {"step", st}};
}

// Load position/tilting from a JSON object into a DetectorTransform.
// Calls prepare() so the rotation matrix is cached before multithreaded use.
static void loadTransform(DetectorTransform &t, const json &j)
{
    if (j.contains("position") && j["position"].is_array() && j["position"].size()>=3) {
        t.x = j["position"][0]; t.y = j["position"][1]; t.z = j["position"][2];
    }
    if (j.contains("tilting") && j["tilting"].is_array() && j["tilting"].size()>=3) {
        t.rx = j["tilting"][0]; t.ry = j["tilting"][1]; t.rz = j["tilting"][2];
    }
    t.prepare();
}

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

    // --- DAQ config + pedestals (required) ---
    std::string daq_cfg_path = daq_config_file;
    if (daq_cfg_path.empty())
        daq_cfg_path = findFile("daq_config.json", db_dir);

    if (daq_cfg_path.empty() || !evc::load_daq_config(daq_cfg_path, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_cfg_path.empty() ? " (not found)" : ": " + daq_cfg_path)
                  << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::cerr << "DAQ config: " << daq_cfg_path
              << " (adc_format=" << daq_cfg.adc_format << ")\n";
    {
        std::ifstream dcf(daq_cfg_path);
        if (dcf.is_open()) {
            auto dcj = json::parse(dcf, nullptr, false, true);
            if (dcj.contains("pedestal_file")) {
                std::string ped_file = findFile(dcj["pedestal_file"].get<std::string>(), db_dir);
                if (evc::load_pedestals(ped_file, daq_cfg))
                    std::cerr << "Pedestals : " << ped_file
                              << " (" << daq_cfg.pedestals.size() << " channels)\n";
            }
        }
    }


    // --- Waveform / histogram config ---
    // Try "waveform" section from config.json (loaded later in reco_file),
    // then fall back to separate hist_config.json, then -H override.
    auto loadWaveformConfig = [&](const json &w) {
        waveform_trigger.parse(w, trigger_bits_def);
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
        if (w.contains("height_hist")) {
            auto &hh = w["height_hist"];
            if (hh.contains("min"))  hist_cfg.height_min  = hh["min"];
            if (hh.contains("max"))  hist_cfg.height_max  = hh["max"];
            if (hh.contains("step")) hist_cfg.height_step = hh["step"];
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
            if (h.contains("height_min"))  hist_cfg.height_min  = h["height_min"];
            if (h.contains("height_max"))  hist_cfg.height_max  = h["height_max"];
            if (h.contains("height_step")) hist_cfg.height_step = h["height_step"];
            if (h.contains("min_peak_ratio")) hist_cfg.min_peak_ratio = h["min_peak_ratio"];
        }
    };

    // load trigger definitions first — needed for name resolution in trigger filters
    {
        std::string tbpath = findFile("trigger_bits.json", db_dir);
        std::string tbs = readFile(tbpath);
        if (!tbs.empty()) {
            auto tb = json::parse(tbs, nullptr, false);
            if (tb.is_array()) {
                trigger_bits_def = tb;
            } else if (tb.is_object()) {
                if (tb.contains("trigger_bits"))
                    trigger_bits_def = tb["trigger_bits"];
                if (tb.contains("trigger_type"))
                    trigger_type_def = tb["trigger_type"];
            }
        }
    }

    // load waveform config from main config, or legacy hist_config.json
    bool waveform_loaded = false;
    if (!main_config.empty()) {
        std::string s = readFile(main_config);
        if (!s.empty()) {
            auto j = json::parse(s, nullptr, false);
            if (j.contains("waveform")) { loadWaveformConfig(j["waveform"]); waveform_loaded = true; }
            // legacy "hist" section
            else if (j.contains("hist")) { loadLegacyHistConfig(j); waveform_loaded = true; }
            if (j.contains("ref_lines") && j["ref_lines"].is_object())
                ref_lines = j["ref_lines"];
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
    height_nbins = std::max(1, (int)std::ceil(
        (hist_cfg.height_max - hist_cfg.height_min) / hist_cfg.height_step));
    std::cerr << "Waveform  : time_cut=[" << hist_cfg.time_min << "," << hist_cfg.time_max
              << "] threshold=" << hist_cfg.threshold
              << " " << waveform_trigger << "\n";

    // --- HyCal system ---
    std::string modules_filename = "hycal_modules.json";
    std::string daq_filename     = "daq_map.json";
    {
        std::ifstream dcf2(daq_cfg_path);
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

    // --- GEM system (optional) ---
    {
        std::string gem_map_filename = "gem_map.json";
        std::string gem_ped_filename;
        std::ifstream dcf_gem(daq_cfg_path);
        if (dcf_gem.is_open()) {
            auto dcj_gem = json::parse(dcf_gem, nullptr, false, true);
            if (dcj_gem.contains("gem_map_file"))
                gem_map_filename = dcj_gem["gem_map_file"].get<std::string>();
            if (dcj_gem.contains("gem_pedestal_file"))
                gem_ped_filename = dcj_gem["gem_pedestal_file"].get<std::string>();
        }
        std::string gem_map_file = findFile(gem_map_filename, db_dir);
        if (!gem_map_file.empty()) {
            gem_sys.Init(gem_map_file);
            gem_enabled = (gem_sys.GetNDetectors() > 0);
            if (gem_enabled) {
                std::cerr << "GEM       : " << gem_sys.GetNDetectors() << " detectors\n";
                // init per-detector data (identity transform by default, pre-prepared)
                int ndet = gem_sys.GetNDetectors();
                gem_transforms.resize(ndet);
                for (auto &t : gem_transforms) t.prepare();
                gem_occupancy.resize(ndet);
                for (auto &h : gem_occupancy) h.init(GEM_OCC_NX, GEM_OCC_NY);
                if (!gem_ped_filename.empty()) {
                    std::string gem_ped_file = findFile(gem_ped_filename, db_dir);
                    if (!gem_ped_file.empty()) {
                        gem_sys.LoadPedestals(gem_ped_file);
                        std::cerr << "GEM peds  : " << gem_ped_file << "\n";
                    }
                }
            }
        }
    }

    // --- crate_roc map ---
    crate_roc_json = json::object();
    {
        std::ifstream dcf3(daq_cfg_path);
        if (dcf3.is_open()) {
            auto dcj3 = json::parse(dcf3, nullptr, false, true);
            if (dcj3.contains("roc_tags")) {
                for (auto &entry : dcj3["roc_tags"]) {
                    if (entry.contains("crate") && entry.contains("tag")) {
                        // only use data ROCs (type "roc") for crate→tag mapping;
                        // ti_slaves share crate numbers but have different tags
                        std::string rtype = entry.value("type", "");
                        if (!rtype.empty() && rtype != "roc" && rtype != "gem") continue;
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

        if (rcfg.contains("online")) {
            auto &on = rcfg["online"];
            if (on.contains("refresh_ms")) {
                auto &r = on["refresh_ms"];
                if (r.contains("event"))     refresh_event_ms = r["event"];
                if (r.contains("ring"))      refresh_ring_ms  = r["ring"];
                if (r.contains("histogram")) refresh_hist_ms  = r["histogram"];
                if (r.contains("lms"))       refresh_lms_ms   = r["lms"];
            }
        }

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
            cluster_trigger.parse(cc, trigger_bits_def);
            if (cc.contains("energy_hist")) {
                auto &eh = cc["energy_hist"];
                if (eh.contains("min"))  cl_hist_min  = eh["min"];
                if (eh.contains("max"))  cl_hist_max  = eh["max"];
                if (eh.contains("step")) cl_hist_step = eh["step"];
            }
            if (cc.contains("nclusters_hist")) {
                auto &nh = cc["nclusters_hist"];
                if (nh.contains("min"))  nclusters_hist_min  = nh["min"];
                if (nh.contains("max"))  nclusters_hist_max  = nh["max"];
                if (nh.contains("step")) nclusters_hist_step = nh["step"];
            }
            if (cc.contains("nblocks_hist")) {
                auto &bh = cc["nblocks_hist"];
                if (bh.contains("min"))  nblocks_hist_min  = bh["min"];
                if (bh.contains("max"))  nblocks_hist_max  = bh["max"];
                if (bh.contains("step")) nblocks_hist_step = bh["step"];
            }
            std::cerr << "Clustering: min_mod=" << cluster_cfg.min_module_energy
                      << " min_center=" << cluster_cfg.min_center_energy
                      << " min_cluster=" << cluster_cfg.min_cluster_energy
                      << " " << cluster_trigger
                      << " hist=[" << cl_hist_min << "," << cl_hist_max
                      << "]/" << cl_hist_step << "\n";
        }

        if (rcfg.contains("lms_monitor")) {
            auto &lm = rcfg["lms_monitor"];
            lms_trigger.parse(lm, trigger_bits_def);
            if (lm.contains("warn_threshold")) lms_warn_thresh     = lm["warn_threshold"];
            if (lm.contains("warn_min_mean"))  lms_warn_min_mean  = lm["warn_min_mean"];
            if (lm.contains("max_history"))    lms_max_history    = lm["max_history"];
            if (lm.contains("reference_channels")) {
                for (auto &name : lm["reference_channels"]) {
                    std::string n = name.get<std::string>();
                    const auto *mod = hycal.module_by_name(n);
                    lms_ref_channels.push_back({n, mod ? mod->index : -1});
                }
            }
            std::cerr << "LMS       : " << lms_trigger
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

        if (rcfg.contains("runinfo")) {
            // runinfo can be inline object or a string path to an external file
            json ri;
            if (rcfg["runinfo"].is_string()) {
                std::string ri_file = findFile(rcfg["runinfo"].get<std::string>(), db_dir);
                if (!ri_file.empty()) {
                    std::ifstream rif(ri_file);
                    if (rif.is_open()) {
                        ri = json::parse(rif, nullptr, false, true);
                        std::cerr << "RunInfo   : loaded from " << ri_file << "\n";
                    }
                }
            } else {
                ri = rcfg["runinfo"];
            }
            if (ri.contains("beam_energy")) beam_energy = ri["beam_energy"];
            if (ri.contains("target") && ri["target"].is_array() && ri["target"].size()>=3) {
                target_x=ri["target"][0]; target_y=ri["target"][1]; target_z=ri["target"][2];
            }
            if (ri.contains("hycal"))
                loadTransform(hycal_transform, ri["hycal"]);
            if (ri.contains("calibration")) {
                auto &cal = ri["calibration"];
                if (cal.contains("default_adc2mev")) adc_to_mev = cal["default_adc2mev"];
                if (cal.contains("file")) {
                    std::string calib_file = findFile(cal["file"].get<std::string>(), db_dir);
                    int nmatched = hycal.LoadCalibration(calib_file);
                    if (nmatched >= 0)
                        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";
                }
            }
            std::cerr << "RunInfo   : beam=" << beam_energy << "MeV default_adc2mev=" << adc_to_mev
                      << " target=(" << target_x << "," << target_y << "," << target_z
                      << ") HyCal=(" << hycal_transform.x << "," << hycal_transform.y << ","
                      << hycal_transform.z << ")\n";

            // GEM per-detector transforms (same position/tilting format as HyCal)
            if (gem_enabled && ri.contains("gem") && ri["gem"].is_array()) {
                for (auto &entry : ri["gem"]) {
                    int id = entry.value("id", -1);
                    if (id < 0 || id >= (int)gem_transforms.size()) continue;
                    loadTransform(gem_transforms[id], entry);
                }
                std::cerr << "GEM geom  : " << gem_transforms.size() << " detectors configured\n";
            }
        }
        if (rcfg.contains("elog")) {
            auto &el = rcfg["elog"];
            if (el.contains("url"))      elog_url      = el["url"];
            if (el.contains("logbook"))  elog_logbook   = el["logbook"];
            if (el.contains("author"))   elog_author    = el["author"];
            if (el.contains("tags"))
                for (auto &t : el["tags"]) elog_tags.push_back(t);
            if (el.contains("cert")) elog_cert = el["cert"];
            if (el.contains("key"))  elog_key  = el["key"];
            std::cerr << "Elog      : " << elog_url
                      << " logbook=" << elog_logbook
                      << (elog_cert.empty() ? "" : " cert=" + elog_cert)
                      << "\n";
        }

        if (rcfg.contains("physics")) {
            auto &ph = rcfg["physics"];
            physics_trigger.parse(ph, trigger_bits_def);
            if (ph.contains("energy_angle_hist")) {
                auto &ea = ph["energy_angle_hist"];
                if (ea.contains("angle_min"))   ea_angle_min   = ea["angle_min"];
                if (ea.contains("angle_max"))   ea_angle_max   = ea["angle_max"];
                if (ea.contains("angle_step"))  ea_angle_step  = ea["angle_step"];
                if (ea.contains("energy_min"))  ea_energy_min  = ea["energy_min"];
                if (ea.contains("energy_max"))  ea_energy_max  = ea["energy_max"];
                if (ea.contains("energy_step")) ea_energy_step = ea["energy_step"];
            }
            if (ph.contains("moller")) {
                auto &ml = ph["moller"];
                if (ml.contains("energy_tolerance")) moller_energy_tol = ml["energy_tolerance"];
                if (ml.contains("angle_min"))        moller_angle_min  = ml["angle_min"];
                if (ml.contains("angle_max"))        moller_angle_max  = ml["angle_max"];
                if (ml.contains("xy_hist")) {
                    auto &xy = ml["xy_hist"];
                    if (xy.contains("x_min"))  moller_xy_x_min  = xy["x_min"];
                    if (xy.contains("x_max"))  moller_xy_x_max  = xy["x_max"];
                    if (xy.contains("x_step")) moller_xy_x_step = xy["x_step"];
                    if (xy.contains("y_min"))  moller_xy_y_min  = xy["y_min"];
                    if (xy.contains("y_max"))  moller_xy_y_max  = xy["y_max"];
                    if (xy.contains("y_step")) moller_xy_y_step = xy["y_step"];
                }
                if (ml.contains("energy_hist")) {
                    auto &eh = ml["energy_hist"];
                    if (eh.contains("min"))  moller_e_min  = eh["min"];
                    if (eh.contains("max"))  moller_e_max  = eh["max"];
                    if (eh.contains("step")) moller_e_step = eh["step"];
                }
            }
            std::cerr << "Physics   : " << physics_trigger
                      << " Moller: tol=" << moller_energy_tol
                      << " angle=[" << moller_angle_min << "," << moller_angle_max << "]\n";
        }

        if (rcfg.contains("epics")) {
            auto &ep = rcfg["epics"];
            if (ep.contains("max_history"))     epics_max_history  = ep["max_history"];
            if (ep.contains("warn_threshold"))  epics_warn_thresh  = ep["warn_threshold"];
            if (ep.contains("alert_threshold")) epics_alert_thresh = ep["alert_threshold"];
            if (ep.contains("min_avg_points"))  epics_min_avg_pts  = ep["min_avg_points"];
            if (ep.contains("mean_window"))     epics_mean_window  = ep["mean_window"];
            if (ep.contains("slots")) {
                for (auto &slot : ep["slots"]) {
                    std::vector<std::string> names;
                    for (auto &ch : slot) names.push_back(ch);
                    epics_default_slots.push_back(std::move(names));
                }
            }
            std::cerr << "EPICS     : max_history=" << epics_max_history
                      << " slots=" << epics_default_slots.size() << "\n";
        }

        // GEM histogram config (optional section in config.json)
        if (rcfg.contains("gem_histograms")) {
            auto &gh = rcfg["gem_histograms"];
            if (gh.contains("nclusters")) {
                auto &nc = gh["nclusters"];
                if (nc.contains("min"))  gem_ncl_min  = nc["min"];
                if (nc.contains("max"))  gem_ncl_max  = nc["max"];
                if (nc.contains("step")) gem_ncl_step = nc["step"];
            }
            if (gh.contains("theta")) {
                auto &th = gh["theta"];
                if (th.contains("min"))  gem_theta_min  = th["min"];
                if (th.contains("max"))  gem_theta_max  = th["max"];
                if (th.contains("step")) gem_theta_step = th["step"];
            }
        }

        std::cerr << "Reco      : " << main_config
                  << " (adc_to_mev=" << adc_to_mev << ")\n";
    }

    // init cluster histograms
    int cl_nbins = std::max(1, (int)std::ceil((cl_hist_max - cl_hist_min) / cl_hist_step));
    cluster_energy_hist.init(cl_nbins);
    nclusters_hist.init(std::max(1, (nclusters_hist_max - nclusters_hist_min) / nclusters_hist_step));
    nblocks_hist.init(std::max(1, (nblocks_hist_max - nblocks_hist_min) / nblocks_hist_step));
    int ea_nx = std::max(1, (int)std::ceil((ea_angle_max - ea_angle_min) / ea_angle_step));
    int ea_ny = std::max(1, (int)std::ceil((ea_energy_max - ea_energy_min) / ea_energy_step));
    energy_angle_hist.init(ea_nx, ea_ny);
    int ml_nx = std::max(1, (int)std::ceil((moller_xy_x_max - moller_xy_x_min) / moller_xy_x_step));
    int ml_ny = std::max(1, (int)std::ceil((moller_xy_y_max - moller_xy_y_min) / moller_xy_y_step));
    moller_xy_hist.init(ml_nx, ml_ny);
    moller_energy_hist.init(std::max(1, (int)std::ceil((moller_e_max - moller_e_min) / moller_e_step)));
    gem_nclusters_hist.init(std::max(1, (gem_ncl_max - gem_ncl_min) / gem_ncl_step));
    gem_theta_hist.init(std::max(1, (int)std::ceil((gem_theta_max - gem_theta_min) / gem_theta_step)));

    // ensure all transforms are pre-prepared before multithreaded use
    hycal_transform.prepare();
}

//=============================================================================
// Per-event processing
//=============================================================================

// Encode peak array for one channel.
static json encodePeaks(const fdec::WaveResult &wres)
{
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
    return parr;
}

json AppState::encodeEventJson(fdec::EventData &event, int ev_id,
                               fdec::WaveAnalyzer &ana, fdec::WaveResult &wres,
                               bool include_samples)
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

                json ch_j = {
                    {"pm", std::round(wres.ped.mean * 10) / 10},
                    {"pr", std::round(wres.ped.rms * 10) / 10},
                    {"pk", encodePeaks(wres)},
                };
                if (include_samples) {
                    json sarr = json::array();
                    for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);
                    ch_j["s"] = std::move(sarr);
                }
                channels[key] = std::move(ch_j);
            }
        }
    }
    return {{"event", ev_id}, {"channels", channels},
            {"event_number", event.info.event_number},
            {"trigger_type", event.info.trigger_type},
            {"trigger_bits", event.info.trigger_bits}};
}

json AppState::encodeWaveformJson(fdec::EventData &event, const std::string &chan_key,
                                  fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // parse "roc_slot_ch" key
    int roc_tag = 0, sl = 0, ch = 0;
    if (std::sscanf(chan_key.c_str(), "%d_%d_%d", &roc_tag, &sl, &ch) != 3)
        return {{"error", "invalid channel key"}};

    // find the channel in the event
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present || roc.tag != roc_tag) continue;
        if (!roc.slots[sl].present) break;
        if (!(roc.slots[sl].channel_mask & (1ull << ch))) break;
        auto &cd = roc.slots[sl].channels[ch];
        if (cd.nsamples <= 0) break;

        ana.Analyze(cd.samples, cd.nsamples, wres);

        json sarr = json::array();
        for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);

        return {{"key", chan_key}, {"s", sarr},
                {"pm", std::round(wres.ped.mean * 10) / 10},
                {"pr", std::round(wres.ped.rms * 10) / 10},
                {"pk", encodePeaks(wres)}};
    }
    return {{"error", "channel not found"}};
}

json AppState::computeClustersJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    if (!cluster_trigger(event.info.trigger_bits))
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
    if (sync_unix != 0) return;   // already have a sync reference

    if (lms_first_ts == 0) {
        // No LMS events yet — stash for later.
        // Will be applied when the first LMS event sets lms_first_ts.
        pending_sync_unix = unix_time;
        pending_sync_ti = last_ti_ts;
        return;
    }

    sync_unix = unix_time;
    sync_rel_sec = (last_ti_ts != 0)
        ? static_cast<double>(last_ti_ts - lms_first_ts) * TI_TICK_SEC
        : 0.;
}

void AppState::processEvent(fdec::EventData &event,
                            fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // --- check which consumers need this event ---
    uint32_t tb = event.info.trigger_bits;
    bool do_hist    = waveform_trigger(tb);
    bool do_cluster = cluster_trigger(tb);
    bool do_lms     = lms_trigger.accept != 0 && lms_trigger(tb);

    if (!do_hist && !do_cluster && !do_lms) {
        std::lock_guard<std::mutex> lk(data_mtx);
        events_processed++;
        return;
    }

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");

    // clustering setup (stack-allocated, per-event)
    fdec::HyCalCluster clusterer(hycal);
    if (do_cluster) clusterer.SetConfig(cluster_cfg);

    // LMS timing
    double lms_time = 0;

    // acquire both locks for the merged pass
    std::unique_lock<std::mutex> lk1(data_mtx, std::defer_lock);
    std::unique_lock<std::mutex> lk2(lms_mtx, std::defer_lock);
    std::lock(lk1, lk2);

    if (do_lms) {
        if (lms_first_ts == 0) {
            lms_first_ts = event.info.timestamp;
            // apply stashed sync time from a control event that arrived before LMS data
            if (pending_sync_unix != 0 && sync_unix == 0) {
                sync_unix = pending_sync_unix;
                // PRESTART/GO arrives before physics events, so pending_sync_ti is
                // typically 0. In that case sync_rel_sec = 0 (run start = LMS start).
                sync_rel_sec = (pending_sync_ti != 0)
                    ? static_cast<double>(pending_sync_ti - lms_first_ts) * TI_TICK_SEC
                    : 0.;
            }
        }
        lms_time = static_cast<double>(event.info.timestamp - lms_first_ts) * TI_TICK_SEC;
    }

    // --- single pass: analyze once per channel, feed all consumers ---
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;

        // crate lookup (needed by cluster + LMS consumers)
        int crate = -1;
        if (do_cluster || do_lms) {
            auto cit = roc_to_crate.find(roc.tag);
            if (cit != roc_to_crate.end()) crate = cit->second;
        }

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                // ── analyze ONCE ──
                float peak_in_window = -1;
                if (!is_adc1881m) {
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    peak_in_window = bestPeakInWindow(wres, hist_cfg.threshold,
                                                       hist_cfg.time_min, hist_cfg.time_max);
                } else {
                    wres.npeaks = 0;
                    peak_in_window = cd.samples[0];
                }

                // ── histogram consumer ──
                if (do_hist && !is_adc1881m) {
                    std::string key = std::to_string(roc.tag) + "_"
                                   + std::to_string(s) + "_" + std::to_string(c);
                    bool has_peak = false, has_peak_tcut = false;
                    float best = -1, best_height = -1;
                    for (int p = 0; p < wres.npeaks; ++p) {
                        auto &pk = wres.peaks[p];
                        if (pk.height < hist_cfg.threshold) continue;
                        has_peak = true;
                        if (pk.time >= hist_cfg.time_min && pk.time <= hist_cfg.time_max) {
                            has_peak_tcut = true;
                            if (pk.integral > best) { best = pk.integral; best_height = pk.height; }
                        }
                    }
                    if (best >= 0) {
                        auto &h = histograms[key];
                        if (h.bins.empty()) h.init(hist_nbins);
                        h.fill(best, hist_cfg.bin_min, hist_cfg.bin_step);
                        auto &hh = height_histograms[key];
                        if (hh.bins.empty()) hh.init(height_nbins);
                        hh.fill(best_height, hist_cfg.height_min, hist_cfg.height_step);
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

                // ── cluster consumer ──
                if (do_cluster && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod && mod->is_hycal()) {
                        float adc_val = is_adc1881m ? (float)cd.samples[0] : peak_in_window;
                        if (adc_val > 0) {
                            float energy = (mod->cal_factor > 0.)
                                ? static_cast<float>(mod->energize(adc_val))
                                : adc_val * adc_to_mev;
                            clusterer.AddHit(mod->index, energy);
                        }
                    }
                }

                // ── LMS consumer ──
                if (do_lms && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod) {
                        float val = is_adc1881m ? (float)cd.samples[0] : peak_in_window;
                        if (val > 0) {
                            auto &hist = lms_history[mod->index];
                            if (static_cast<int>(hist.size()) < lms_max_history)
                                hist.push_back({lms_time, val});
                        }
                    }
                }
            }
        }
    }

    // --- post-loop: clustering + physics histograms ---
    if (do_cluster) {
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> reco_hits;
        clusterer.ReconstructHits(reco_hits);

        struct ClusterInfo { float lx, ly, lz, theta; };
        std::vector<ClusterInfo> cinfo(reco_hits.size());
        for (size_t i = 0; i < reco_hits.size(); ++i) {
            auto &rh = reco_hits[i];
            auto &ci = cinfo[i];
            hycal_transform.toLab(rh.x, rh.y, ci.lx, ci.ly, ci.lz);
            float dx = ci.lx - target_x, dy = ci.ly - target_y, dz = ci.lz - target_z;
            float rv = std::sqrt(dx*dx + dy*dy);
            ci.theta = std::atan2(rv, dz) * (180.f / 3.14159265f);
        }

        for (size_t i = 0; i < reco_hits.size(); ++i) {
            cluster_energy_hist.fill(reco_hits[i].energy, cl_hist_min, cl_hist_step);
            nblocks_hist.fill(reco_hits[i].nblocks, nblocks_hist_min, nblocks_hist_step);
        }
        nclusters_hist.fill(reco_hits.size(), nclusters_hist_min, nclusters_hist_step);
        cluster_events_processed++;

        bool physics_accept = physics_trigger(tb);
        if (physics_accept) {
            for (size_t i = 0; i < reco_hits.size(); ++i) {
                energy_angle_hist.fill(cinfo[i].theta, reco_hits[i].energy,
                    ea_angle_min, ea_angle_step, ea_energy_min, ea_energy_step);
            }
            if (reco_hits.size() == 2 && beam_energy > 0) {
                float esum = reco_hits[0].energy + reco_hits[1].energy;
                bool energy_ok = std::abs(esum - beam_energy) < moller_energy_tol * beam_energy;
                bool angle_ok = false;
                for (int j = 0; j < 2; ++j)
                    if (cinfo[j].theta >= moller_angle_min && cinfo[j].theta <= moller_angle_max)
                        angle_ok = true;
                if (energy_ok && angle_ok) {
                    moller_events++;
                    for (int j = 0; j < 2; ++j) {
                        moller_xy_hist.fill(cinfo[j].lx, cinfo[j].ly,
                            moller_xy_x_min, moller_xy_x_step, moller_xy_y_min, moller_xy_y_step);
                        moller_energy_hist.fill(reco_hits[j].energy, moller_e_min, moller_e_step);
                    }
                }
            }
        }
    }

    events_processed++;
    if (do_lms) lms_events++;
}

void AppState::processReconEvent(const ReconEventData &recon)
{
    uint32_t tb = recon.trigger_bits;
    bool do_cluster = cluster_trigger(tb);
    bool do_physics = physics_trigger(tb);

    std::lock_guard<std::mutex> lk(data_mtx);
    events_processed++;

    if (do_cluster && !recon.clusters.empty()) {
        for (auto &cl : recon.clusters) {
            cluster_energy_hist.fill(cl.energy, cl_hist_min, cl_hist_step);
            nblocks_hist.fill(cl.nblocks, nblocks_hist_min, nblocks_hist_step);
        }
        nclusters_hist.fill(recon.clusters.size(), nclusters_hist_min, nclusters_hist_step);
        cluster_events_processed++;
    }

    if (do_physics && !recon.clusters.empty()) {
        struct CI { float lx, ly, lz, theta; };
        std::vector<CI> cinfo(recon.clusters.size());
        for (size_t i = 0; i < recon.clusters.size(); ++i) {
            auto &cl = recon.clusters[i];
            auto &ci = cinfo[i];
            hycal_transform.toLab(cl.x, cl.y, ci.lx, ci.ly, ci.lz);
            float dx = ci.lx - target_x, dy = ci.ly - target_y, dz = ci.lz - target_z;
            float r = std::sqrt(dx*dx + dy*dy);
            ci.theta = std::atan2(r, dz) * (180.f / 3.14159265f);
        }
        for (size_t i = 0; i < recon.clusters.size(); ++i)
            energy_angle_hist.fill(cinfo[i].theta, recon.clusters[i].energy,
                ea_angle_min, ea_angle_step, ea_energy_min, ea_energy_step);

        if (recon.clusters.size() == 2 && beam_energy > 0) {
            float esum = recon.clusters[0].energy + recon.clusters[1].energy;
            bool energy_ok = std::abs(esum - beam_energy) < moller_energy_tol * beam_energy;
            bool angle_ok = false;
            for (int j = 0; j < 2; ++j)
                if (cinfo[j].theta >= moller_angle_min && cinfo[j].theta <= moller_angle_max)
                    angle_ok = true;
            if (energy_ok && angle_ok) {
                moller_events++;
                for (int j = 0; j < 2; ++j) {
                    moller_xy_hist.fill(cinfo[j].lx, cinfo[j].ly,
                        moller_xy_x_min, moller_xy_x_step, moller_xy_y_min, moller_xy_y_step);
                    moller_energy_hist.fill(recon.clusters[j].energy, moller_e_min, moller_e_step);
                }
            }
        }
    }
}

json AppState::encodeReconClustersJson(const ReconEventData &recon, int ev_id)
{
    json hits_j = json::object();
    json cl_arr = json::array();

    for (size_t i = 0; i < recon.clusters.size(); ++i) {
        auto &cl = recon.clusters[i];
        std::string center_name;
        if (cl.center_id >= 0 && cl.center_id < hycal.module_count())
            center_name = hycal.module(cl.center_id).name;
        hits_j[std::to_string(cl.center_id)] =
            std::round(cl.energy * 100) / 100;
        cl_arr.push_back({
            {"id", (int)i}, {"center", center_name},
            {"center_id", cl.center_id},
            {"x", std::round(cl.x * 10) / 10},
            {"y", std::round(cl.y * 10) / 10},
            {"energy", std::round(cl.energy * 10) / 10},
            {"nblocks", cl.nblocks}, {"npos", 0},
            {"modules", json::array({cl.center_id})},
        });
    }
    return {{"event", ev_id}, {"hits", hits_j}, {"clusters", cl_arr}};
}

void AppState::processGemEvent(const ssp::SspEventData &ssp_evt)
{
    if (!gem_enabled || ssp_evt.nmpds == 0) return;
    gem_sys.Clear();
    gem_sys.ProcessEvent(ssp_evt);
    gem_sys.Reconstruct(gem_clusterer);

    // accumulate occupancy + histograms in a single pass
    std::lock_guard<std::mutex> lk(data_mtx);
    int total_clusters = 0;
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        float xSize = det.planes[0].size;
        float ySize = det.planes[1].size;
        float xStep = xSize / GEM_OCC_NX;
        float yStep = ySize / GEM_OCC_NY;
        auto &xform = gem_transforms[d];
        auto &hits = gem_sys.GetHits(d);
        total_clusters += static_cast<int>(hits.size());
        for (auto &h : hits) {
            // rotation only for occupancy (local detector coords)
            float ox, oy;
            xform.rotate(h.x, h.y, ox, oy);
            gem_occupancy[d].fill(ox, oy, -xSize/2, xStep, -ySize/2, yStep);
            // full transform for theta (lab frame)
            float lx, ly, lz;
            xform.toLab(h.x, h.y, lx, ly, lz);
            float r = std::sqrt(lx*lx + ly*ly);
            float theta = std::atan2(r, lz) * (180.f / 3.14159265f);
            gem_theta_hist.fill(theta, gem_theta_min, gem_theta_step);
        }
    }
    gem_nclusters_hist.fill(static_cast<float>(total_clusters),
                            static_cast<float>(gem_ncl_min),
                            static_cast<float>(gem_ncl_step));
}

//=============================================================================
// GEM API builders
//=============================================================================

nlohmann::json AppState::apiGemHits() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    result["n_detectors"] = gem_sys.GetNDetectors();
    json detectors = json::array();
    json all_hits = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json dj;
        dj["name"] = det.name;
        dj["id"] = det.id;

        // 1D clusters per plane
        for (int p = 0; p < 2; ++p) {
            std::string pname = (p == 0) ? "x_clusters" : "y_clusters";
            json clusters = json::array();
            for (auto &cl : gem_sys.GetPlaneClusters(d, p)) {
                clusters.push_back({
                    {"position", cl.position},
                    {"peak_charge", cl.peak_charge},
                    {"total_charge", cl.total_charge},
                    {"size", (int)cl.hits.size()},
                    {"max_timebin", cl.max_timebin}
                });
            }
            dj[pname] = clusters;
        }

        // 2D hits (transformed to lab frame) — build per-det and all_hits in one pass
        auto &xform = gem_transforms[d];
        json hits = json::array();
        for (auto &h : gem_sys.GetHits(d)) {
            float lx, ly, lz;
            xform.toLab(h.x, h.y, lx, ly, lz);
            hits.push_back({
                {"x", lx}, {"y", ly},
                {"x_charge", h.x_charge}, {"y_charge", h.y_charge},
                {"x_size", h.x_size}, {"y_size", h.y_size}
            });
            all_hits.push_back({
                {"x", lx}, {"y", ly}, {"det", d},
                {"x_charge", h.x_charge}, {"y_charge", h.y_charge}
            });
        }
        dj["hits_2d"] = hits;
        detectors.push_back(dj);
    }
    result["detectors"] = detectors;
    result["all_hits"] = all_hits;
    return result;
}

nlohmann::json AppState::apiGemConfig() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    result["n_detectors"] = gem_sys.GetNDetectors();
    json layers = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json lj = {
            {"id", det.id},
            {"name", det.name},
            {"type", det.type},
            {"x_pitch", det.planes[0].pitch},
            {"y_pitch", det.planes[1].pitch},
            {"x_apvs", det.planes[0].n_apvs},
            {"y_apvs", det.planes[1].n_apvs},
            {"x_size", det.planes[0].size},
            {"y_size", det.planes[1].size}
        };
        auto &t = gem_transforms[d];
        lj["position"] = {t.x, t.y, t.z};
        lj["tilting"]  = {t.rx, t.ry, t.rz};
        layers.push_back(lj);
    }
    result["layers"] = layers;
    result["occ_nx"] = GEM_OCC_NX;
    result["occ_ny"] = GEM_OCC_NY;
    return result;
}

nlohmann::json AppState::apiGemOccupancy() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    std::lock_guard<std::mutex> lk(data_mtx);
    json dets = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json dj;
        dj["id"] = det.id;
        dj["name"] = det.name;
        dj["x_size"] = det.planes[0].size;
        dj["y_size"] = det.planes[1].size;
        dj["nx"] = GEM_OCC_NX;
        dj["ny"] = GEM_OCC_NY;
        dj["bins"] = gem_occupancy[d].bins;
        dets.push_back(dj);
    }
    result["detectors"] = dets;
    result["total"] = events_processed.load();
    return result;
}

nlohmann::json AppState::apiGemHist() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {
        {"nclusters", histToJson(gem_nclusters_hist, (float)gem_ncl_min, (float)gem_ncl_max, (float)gem_ncl_step)},
        {"theta",     histToJson(gem_theta_hist, gem_theta_min, gem_theta_max, gem_theta_step)}
    };
}

//=============================================================================
// Clearing
//=============================================================================

void AppState::clearHistograms()
{
    std::lock_guard<std::mutex> lk(data_mtx);
    for (auto &[k, h] : histograms)        h.clear();
    for (auto &[k, h] : pos_histograms)   h.clear();
    for (auto &[k, h] : height_histograms) h.clear();
    occupancy.clear();
    occupancy_tcut.clear();
    events_processed = 0;
    cluster_energy_hist.clear();
    nclusters_hist.clear();
    nblocks_hist.clear();
    energy_angle_hist.clear();
    moller_xy_hist.clear();
    moller_energy_hist.clear();
    moller_events = 0;
    cluster_events_processed = 0;
    for (auto &h : gem_occupancy) h.clear();
    gem_nclusters_hist.clear();
    gem_theta_hist.clear();
}

void AppState::clearLms()
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    lms_history.clear();
    lms_events = 0;
    lms_first_ts = 0;
    sync_unix = 0;
    sync_rel_sec = 0.;
    pending_sync_unix = 0;
    pending_sync_ti = 0;
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

json AppState::apiHist(int type, const std::string &key) const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    auto &hmap = (type == 0) ? histograms : (type == 1) ? pos_histograms : height_histograms;
    int nbins  = (type == 0) ? hist_nbins : (type == 1) ? pos_nbins : height_nbins;
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
    json r = histToJson(cluster_energy_hist, cl_hist_min, cl_hist_max, cl_hist_step);
    r["events"] = cluster_events_processed;
    r["nclusters"] = histToJson(nclusters_hist,
        (float)nclusters_hist_min, (float)nclusters_hist_max, (float)nclusters_hist_step);
    r["nblocks"] = histToJson(nblocks_hist,
        (float)nblocks_hist_min, (float)nblocks_hist_max, (float)nblocks_hist_step);
    return r;
}

json AppState::apiEnergyAngle() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {{"bins", energy_angle_hist.bins},
            {"nx", energy_angle_hist.nx}, {"ny", energy_angle_hist.ny},
            {"angle_min", ea_angle_min}, {"angle_max", ea_angle_max}, {"angle_step", ea_angle_step},
            {"energy_min", ea_energy_min}, {"energy_max", ea_energy_max}, {"energy_step", ea_energy_step},
            {"target", {target_x, target_y, target_z}},
            {"hycal_z", hycal_transform.z},
            {"beam_energy", beam_energy},
            {"events", cluster_events_processed}};
}

json AppState::apiMoller() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {{"xy_bins", moller_xy_hist.bins},
            {"xy_nx", moller_xy_hist.nx}, {"xy_ny", moller_xy_hist.ny},
            {"xy_x_min", moller_xy_x_min}, {"xy_x_max", moller_xy_x_max}, {"xy_x_step", moller_xy_x_step},
            {"xy_y_min", moller_xy_y_min}, {"xy_y_max", moller_xy_y_max}, {"xy_y_step", moller_xy_y_step},
            {"energy_hist", histToJson(moller_energy_hist, moller_e_min, moller_e_max, moller_e_step)},
            {"moller_events", moller_events},
            {"total_events", cluster_events_processed},
            {"cuts", {{"energy_tolerance", moller_energy_tol},
                      {"angle_min", moller_angle_min}, {"angle_max", moller_angle_max}}}};
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
            {"trigger", lms_trigger.toJson()},
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

//=============================================================================
// EPICS
//=============================================================================

void AppState::processEpics(const std::string &text, int32_t event_number, uint64_t timestamp)
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    epics.Feed(event_number, timestamp, text);
    epics.Trim(epics_max_history);
    epics_events++;
}

void AppState::clearEpics()
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    epics.Clear();
    epics_events = 0;
}

json AppState::apiEpicsChannels() const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    json names = json::array();
    for (auto &n : epics.GetChannelNames()) names.push_back(n);
    json slots = json::array();
    for (auto &s : epics_default_slots) slots.push_back(s);
    return {{"channels", names}, {"slots", slots},
            {"events", epics_events.load()}};
}

json AppState::apiEpicsChannel(const std::string &name) const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    int id = epics.GetChannelId(name);
    if (id < 0)
        return {{"name", name}, {"time", json::array()}, {"value", json::array()}, {"count", 0}};

    int nsnap = epics.GetSnapshotCount();
    json t_arr = json::array(), v_arr = json::array();

    // time relative to first snapshot's timestamp
    uint64_t t0 = (nsnap > 0) ? epics.GetSnapshot(0).timestamp : 0;
    for (int i = 0; i < nsnap; ++i) {
        auto &snap = epics.GetSnapshot(i);
        double t_sec = static_cast<double>(snap.timestamp - t0) * TI_TICK_SEC;
        float val = (id < (int)snap.values.size()) ? snap.values[id] : 0.f;
        t_arr.push_back(std::round(t_sec * 100) / 100);
        v_arr.push_back(val);
    }
    return {{"name", name}, {"time", t_arr}, {"value", v_arr}, {"count", nsnap}};
}

json AppState::apiEpicsBatch(const std::vector<std::string> &names) const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    int nsnap = epics.GetSnapshotCount();
    uint64_t t0 = (nsnap > 0) ? epics.GetSnapshot(0).timestamp : 0;

    // build shared time array once
    json t_arr = json::array();
    for (int i = 0; i < nsnap; ++i) {
        double t_sec = static_cast<double>(epics.GetSnapshot(i).timestamp - t0) * TI_TICK_SEC;
        t_arr.push_back(std::round(t_sec * 100) / 100);
    }

    json channels = json::array();
    for (auto &name : names) {
        int id = epics.GetChannelId(name);
        if (id < 0) {
            channels.push_back({{"name", name}, {"value", json::array()}, {"count", 0}});
            continue;
        }
        json v_arr = json::array();
        for (int i = 0; i < nsnap; ++i) {
            auto &snap = epics.GetSnapshot(i);
            v_arr.push_back((id < (int)snap.values.size()) ? snap.values[id] : 0.f);
        }
        channels.push_back({{"name", name}, {"value", v_arr}, {"count", nsnap}});
    }
    return {{"time", t_arr}, {"channels", channels}};
}

json AppState::apiEpicsLatest() const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    json channels = json::array();
    int nsnap = epics.GetSnapshotCount();
    int nch = epics.GetChannelCount();
    if (nsnap == 0 || nch == 0)
        return {{"channels", channels}, {"events", epics_events.load()}};

    auto &latest = epics.GetSnapshot(nsnap - 1);

    // compute per-channel mean from most recent mean_window snapshots
    int win_start = std::max(0, nsnap - epics_mean_window);
    std::vector<double> sums(nch, 0.0);
    std::vector<int> counts(nch, 0);
    for (int i = win_start; i < nsnap; ++i) {
        auto &snap = epics.GetSnapshot(i);
        for (int ch = 0; ch < std::min(nch, (int)snap.values.size()); ++ch) {
            sums[ch] += snap.values[ch];
            counts[ch]++;
        }
    }

    for (int ch = 0; ch < nch; ++ch) {
        float val = (ch < (int)latest.values.size()) ? latest.values[ch] : 0.f;
        float mean = (counts[ch] > 0) ? static_cast<float>(sums[ch] / counts[ch]) : val;
        channels.push_back({
            {"name", epics.GetChannelName(ch)},
            {"value", std::round(val * 1000) / 1000},
            {"mean", std::round(mean * 1000) / 1000},
            {"count", counts[ch]},
        });
    }
    return {{"channels", channels}, {"events", epics_events.load()}};
}

//=============================================================================
// Shared config + API routing (used by both viewer and monitor)
//=============================================================================

void AppState::fillConfigJson(json &cfg) const
{
    cfg["hist"] = {
        {"time_min", hist_cfg.time_min}, {"time_max", hist_cfg.time_max},
        {"bin_min", hist_cfg.bin_min}, {"bin_max", hist_cfg.bin_max},
        {"bin_step", hist_cfg.bin_step}, {"threshold", hist_cfg.threshold},
        {"pos_min", hist_cfg.pos_min}, {"pos_max", hist_cfg.pos_max},
        {"pos_step", hist_cfg.pos_step},
        {"height_min", hist_cfg.height_min}, {"height_max", hist_cfg.height_max},
        {"height_step", hist_cfg.height_step},
    };
    cfg["ref_lines"] = ref_lines;
    cfg["trigger_bits"] = trigger_bits_def;
    cfg["trigger_type"] = trigger_type_def;
    cfg["trigger_filter"] = {
        {"dq",      waveform_trigger.toJson()},
        {"cluster", cluster_trigger.toJson()},
        {"lms",     lms_trigger.toJson()},
        {"physics", physics_trigger.toJson()},
    };
    cfg["cluster_hist"] = {{"min", cl_hist_min}, {"max", cl_hist_max}, {"step", cl_hist_step}};
    cfg["nclusters_hist"] = {{"min", nclusters_hist_min}, {"max", nclusters_hist_max}, {"step", nclusters_hist_step}};
    cfg["nblocks_hist"] = {{"min", nblocks_hist_min}, {"max", nblocks_hist_max}, {"step", nblocks_hist_step}};
    cfg["color_ranges"] = apiColorRanges();
    cfg["refresh_ms"] = {{"event", refresh_event_ms}, {"ring", refresh_ring_ms},
                         {"histogram", refresh_hist_ms}, {"lms", refresh_lms_ms}};
    cfg["lms"] = {
        {"trigger", lms_trigger.toJson()},
        {"warn_threshold", lms_warn_thresh},
        {"events", lms_events.load()}, {"ref_channels", apiLmsRefChannels()},
    };
    cfg["runinfo"] = {
        {"beam_energy", beam_energy},
        {"calibration", {{"default_adc2mev", adc_to_mev}}},
        {"target", {target_x, target_y, target_z}},
        {"hycal", {
            {"position", {hycal_transform.x, hycal_transform.y, hycal_transform.z}},
            {"tilting", {hycal_transform.rx, hycal_transform.ry, hycal_transform.rz}},
        }},
    };
    cfg["physics"] = {
        {"trigger", physics_trigger.toJson()},
        {"energy_angle_hist", {
            {"angle_min", ea_angle_min}, {"angle_max", ea_angle_max}, {"angle_step", ea_angle_step},
            {"energy_min", ea_energy_min}, {"energy_max", ea_energy_max}, {"energy_step", ea_energy_step},
        }},
        {"moller", {
            {"energy_tolerance", moller_energy_tol},
            {"angle_min", moller_angle_min}, {"angle_max", moller_angle_max},
        }},
    };
    cfg["elog"] = {
        {"url", elog_url}, {"logbook", elog_logbook},
        {"author", elog_author}, {"tags", elog_tags},
    };
    cfg["epics"] = {
        {"max_history", epics_max_history},
        {"warn_threshold", epics_warn_thresh}, {"alert_threshold", epics_alert_thresh},
        {"min_avg_points", epics_min_avg_pts}, {"mean_window", epics_mean_window},
        {"slots", epics_default_slots},
    };
    cfg["gem"] = apiGemConfig();
}

AppState::ApiResult AppState::handleReadApi(const std::string &uri) const
{
    if (uri == "/api/occupancy")
        return {true, apiOccupancy().dump()};
    if (uri == "/api/physics/energy_angle")
        return {true, apiEnergyAngle().dump()};
    if (uri == "/api/physics/moller")
        return {true, apiMoller().dump()};
    if (uri == "/api/cluster_hist")
        return {true, apiClusterHist().dump()};
    if (uri.rfind("/api/hist/", 0) == 0)
        return {true, apiHist(0, uri.substr(10)).dump()};
    if (uri.rfind("/api/poshist/", 0) == 0)
        return {true, apiHist(1, uri.substr(13)).dump()};
    if (uri.rfind("/api/heighthist/", 0) == 0)
        return {true, apiHist(2, uri.substr(16)).dump()};
    if (uri == "/api/lms/refs")
        return {true, apiLmsRefChannels().dump()};
    if (uri.rfind("/api/lms/", 0) == 0) {
        int ref = -1;
        auto qpos = uri.find('?');
        std::string path = (qpos != std::string::npos) ? uri.substr(9, qpos - 9) : uri.substr(9);
        if (qpos != std::string::npos) {
            std::string q = uri.substr(qpos + 1);
            if (q.rfind("ref=", 0) == 0) ref = std::atoi(q.c_str() + 4);
        }
        if (path == "summary") return {true, apiLmsSummary(ref).dump()};
        if (path == "clear")   return {false, ""};  // clear handled by caller
        return {true, apiLmsModule(std::atoi(path.c_str()), ref).dump()};
    }
    if (uri.rfind("/api/epics/", 0) == 0) {
        std::string path = uri.substr(11);
        if (path == "channels") return {true, apiEpicsChannels().dump()};
        if (path == "latest")   return {true, apiEpicsLatest().dump()};
        if (path == "clear")    return {false, ""};  // clear handled by caller
        if (path.rfind("batch?", 0) == 0) {
            // /api/epics/batch?ch=name1&ch=name2&...
            std::string query = path.substr(6);
            std::vector<std::string> names;
            for (size_t pos = 0; pos < query.size();) {
                size_t amp = query.find('&', pos);
                if (amp == std::string::npos) amp = query.size();
                std::string kv = query.substr(pos, amp - pos);
                if (kv.rfind("ch=", 0) == 0) {
                    // URL-decode
                    std::string raw = kv.substr(3), name;
                    for (size_t i = 0; i < raw.size(); ++i) {
                        if (raw[i] == '%' && i + 2 < raw.size()) {
                            int hi = 0, lo = 0;
                            if (std::sscanf(raw.c_str() + i + 1, "%1x%1x", &hi, &lo) == 2) {
                                name += static_cast<char>((hi << 4) | lo);
                                i += 2; continue;
                            }
                        }
                        if (raw[i] == '+') name += ' ';
                        else name += raw[i];
                    }
                    names.push_back(name);
                }
                pos = amp + 1;
            }
            return {true, apiEpicsBatch(names).dump()};
        }
        if (path.rfind("channel/", 0) == 0) {
            // URL-decode the channel name (e.g. %3A → :)
            std::string raw = path.substr(8), name;
            for (size_t i = 0; i < raw.size(); ++i) {
                if (raw[i] == '%' && i + 2 < raw.size()) {
                    int hi = 0, lo = 0;
                    if (std::sscanf(raw.c_str() + i + 1, "%1x%1x", &hi, &lo) == 2) {
                        name += static_cast<char>((hi << 4) | lo);
                        i += 2;
                        continue;
                    }
                }
                name += raw[i];
            }
            return {true, apiEpicsChannel(name).dump()};
        }
    }
    if (uri == "/api/gem/hits")
        return {true, apiGemHits().dump()};
    if (uri == "/api/gem/config")
        return {true, apiGemConfig().dump()};
    if (uri == "/api/gem/occupancy")
        return {true, apiGemOccupancy().dump()};
    if (uri == "/api/gem/hist")
        return {true, apiGemHist().dump()};
    return {false, ""};
}

//=============================================================================
// Filters
//=============================================================================

void AppState::resolveFilterKeys()
{
    filter_wf_keys.clear();
    for (auto &name : waveform_filter.modules) {
        const auto *mod = hycal.module_by_name(name);
        if (!mod || mod->daq.crate < 0) continue;
        for (auto &[roc_tag, crate_id] : roc_to_crate) {
            if (crate_id == mod->daq.crate) {
                filter_wf_keys.insert(std::to_string(roc_tag) + "_"
                    + std::to_string(mod->daq.slot) + "_"
                    + std::to_string(mod->daq.channel));
                break;
            }
        }
    }

    filter_cl_includes.clear();
    for (auto &name : cluster_filter.includes_modules) {
        const auto *mod = hycal.module_by_name(name);
        if (mod) filter_cl_includes.insert(mod->index);
    }

    filter_cl_centers.clear();
    for (auto &name : cluster_filter.center_modules) {
        const auto *mod = hycal.module_by_name(name);
        if (mod) filter_cl_centers.insert(mod->id);
    }
}

std::string AppState::loadFilter(const json &j)
{
    waveform_filter = {};
    cluster_filter = {};

    if (j.contains("waveform")) {
        auto &w = j["waveform"];
        auto &f = waveform_filter;
        if (w.contains("enable"))       f.enable       = w["enable"];
        if (w.contains("n_peaks_min"))  f.n_peaks_min  = w["n_peaks_min"];
        if (w.contains("n_peaks_max"))  f.n_peaks_max  = w["n_peaks_max"];
        if (w.contains("time_min"))     f.time_min     = w["time_min"];
        if (w.contains("time_max"))     f.time_max     = w["time_max"];
        if (w.contains("integral_min")) f.integral_min = w["integral_min"];
        if (w.contains("integral_max")) f.integral_max = w["integral_max"];
        if (w.contains("height_min"))   f.height_min   = w["height_min"];
        if (w.contains("height_max"))   f.height_max   = w["height_max"];
        if (w.contains("modules"))
            for (auto &m : w["modules"])
                f.modules.push_back(m.get<std::string>());
    }

    if (j.contains("clustering")) {
        auto &c = j["clustering"];
        auto &f = cluster_filter;
        if (c.contains("enable"))           f.enable       = c["enable"];
        if (c.contains("n_min"))            f.n_min        = c["n_min"];
        if (c.contains("n_max"))            f.n_max        = c["n_max"];
        if (c.contains("energy_min"))       f.energy_min   = c["energy_min"];
        if (c.contains("energy_max"))       f.energy_max   = c["energy_max"];
        if (c.contains("size_min"))         f.size_min     = c["size_min"];
        if (c.contains("size_max"))         f.size_max     = c["size_max"];
        if (c.contains("includes_min"))     f.includes_min = c["includes_min"];
        if (c.contains("includes_modules"))
            for (auto &m : c["includes_modules"])
                f.includes_modules.push_back(m.get<std::string>());
        if (c.contains("center_modules"))
            for (auto &m : c["center_modules"])
                f.center_modules.push_back(m.get<std::string>());
    }

    resolveFilterKeys();

    std::cerr << "Filter loaded: waveform "
              << (waveform_filter.enable ? "ON" : "off")
              << " (" << filter_wf_keys.size() << " modules)"
              << ", cluster "
              << (cluster_filter.enable ? "ON" : "off") << "\n";
    return "";
}

void AppState::unloadFilter()
{
    waveform_filter = {};
    cluster_filter = {};
    filter_wf_keys.clear();
    filter_cl_includes.clear();
    filter_cl_centers.clear();
    std::cerr << "Filters unloaded\n";
}

json AppState::filterToJson() const
{
    json r;
    r["active"] = filterActive();
    {
        auto &f = waveform_filter;
        json w;
        w["enable"] = f.enable;
        if (!f.modules.empty()) w["modules"] = f.modules;
        w["n_peaks_min"] = f.n_peaks_min;
        w["n_peaks_max"] = f.n_peaks_max;
        if (f.time_min > -1e20f)     w["time_min"]     = f.time_min;
        if (f.time_max <  1e20f)     w["time_max"]     = f.time_max;
        if (f.integral_min > -1e20f) w["integral_min"] = f.integral_min;
        if (f.integral_max <  1e20f) w["integral_max"] = f.integral_max;
        if (f.height_min > -1e20f)   w["height_min"]   = f.height_min;
        if (f.height_max <  1e20f)   w["height_max"]   = f.height_max;
        r["waveform"] = w;
    }
    {
        auto &f = cluster_filter;
        json c;
        c["enable"] = f.enable;
        c["n_min"] = f.n_min;
        c["n_max"] = f.n_max;
        if (f.energy_min > 0)    c["energy_min"] = f.energy_min;
        if (f.energy_max < 1e20f) c["energy_max"] = f.energy_max;
        c["size_min"] = f.size_min;
        c["size_max"] = f.size_max;
        if (!f.includes_modules.empty()) {
            c["includes_modules"] = f.includes_modules;
            c["includes_min"] = f.includes_min;
        }
        if (!f.center_modules.empty()) c["center_modules"] = f.center_modules;
        r["clustering"] = c;
    }
    return r;
}

bool AppState::evaluateFilter(fdec::EventData &event,
                              ssp::SspEventData *ssp) const
{
    if (!filterActive()) return true;

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");

    // --- waveform filter ---
    if (waveform_filter.enable) {
        if (filter_wf_keys.empty()) return false;  // enabled but no modules resolved

        bool any_module_pass = false;
        fdec::WaveAnalyzer ana;
        ana.cfg.min_peak_ratio = hist_cfg.min_peak_ratio;
        fdec::WaveResult wres;

        for (int r = 0; r < event.nrocs && !any_module_pass; ++r) {
            auto &roc = event.rocs[r];
            if (!roc.present) continue;
            for (int s = 0; s < fdec::MAX_SLOTS && !any_module_pass; ++s) {
                if (!roc.slots[s].present) continue;
                auto &slot = roc.slots[s];
                for (int c = 0; c < fdec::MAX_CHANNELS && !any_module_pass; ++c) {
                    if (!(slot.channel_mask & (1ull << c))) continue;
                    auto &cd = slot.channels[c];
                    if (cd.nsamples <= 0) continue;

                    std::string key = std::to_string(roc.tag) + "_"
                                    + std::to_string(s) + "_" + std::to_string(c);
                    if (filter_wf_keys.find(key) == filter_wf_keys.end()) continue;

                    if (is_adc1881m) continue;  // no peak analysis for ADC1881M

                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    int n_qual = 0;
                    for (int p = 0; p < wres.npeaks; ++p) {
                        auto &pk = wres.peaks[p];
                        if (pk.height < waveform_filter.height_min ||
                            pk.height > waveform_filter.height_max) continue;
                        if (pk.time < waveform_filter.time_min ||
                            pk.time > waveform_filter.time_max) continue;
                        if (pk.integral < waveform_filter.integral_min ||
                            pk.integral > waveform_filter.integral_max) continue;
                        n_qual++;
                    }
                    if (n_qual >= waveform_filter.n_peaks_min &&
                        n_qual <= waveform_filter.n_peaks_max)
                        any_module_pass = true;
                }
            }
        }
        if (!any_module_pass) return false;
    }

    // --- cluster filter ---
    if (cluster_filter.enable) {
        // build clusters from this event (local clusterer, no side effects)
        fdec::HyCalCluster clusterer(hycal);
        clusterer.SetConfig(cluster_cfg);

        fdec::WaveAnalyzer ana;
        ana.cfg.min_peak_ratio = hist_cfg.min_peak_ratio;
        fdec::WaveResult wres;

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

                    float adc_val;
                    if (is_adc1881m) {
                        adc_val = (float)cd.samples[0];
                    } else {
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        adc_val = bestPeakInWindow(wres, hist_cfg.threshold,
                                                   hist_cfg.time_min, hist_cfg.time_max);
                    }
                    if (adc_val > 0) {
                        float energy = (mod->cal_factor > 0.)
                            ? static_cast<float>(mod->energize(adc_val))
                            : adc_val * adc_to_mev;
                        clusterer.AddHit(mod->index, energy);
                    }
                }
            }
        }

        clusterer.FormClusters();
        std::vector<fdec::HyCalCluster::RecoResult> results;
        clusterer.ReconstructMatched(results);

        int n_qual = 0;
        for (auto &rr : results) {
            auto &hit = rr.hit;
            auto *cl = rr.cluster;

            // energy range
            if (hit.energy < cluster_filter.energy_min ||
                hit.energy > cluster_filter.energy_max) continue;
            // size range
            if (hit.nblocks < cluster_filter.size_min ||
                hit.nblocks > cluster_filter.size_max) continue;
            // includes check
            if (!filter_cl_includes.empty()) {
                int count = 0;
                for (auto &h : cl->hits)
                    if (filter_cl_includes.count(h.index)) count++;
                if (count < cluster_filter.includes_min) continue;
            }
            // center check
            if (!filter_cl_centers.empty()) {
                if (!filter_cl_centers.count(hit.center_id)) continue;
            }
            n_qual++;
        }
        if (n_qual < cluster_filter.n_min || n_qual > cluster_filter.n_max)
            return false;
    }

    return true;
}
