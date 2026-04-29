#include "app_state.h"
#include "data_source.h"
#include "load_daq_config.h"
#include "RunInfoConfig.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

using json = nlohmann::json;

static void setTransform(DetectorTransform &t,
                         float x, float y, float z,
                         float rx, float ry, float rz)
{
    t.x = x; t.y = y; t.z = z;
    t.rx = rx; t.ry = ry; t.rz = rz;
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
        std::ifstream dcf_gem(daq_cfg_path);
        if (dcf_gem.is_open()) {
            auto dcj_gem = json::parse(dcf_gem, nullptr, false, true);
            if (dcj_gem.contains("gem_map_file"))
                gem_map_filename = dcj_gem["gem_map_file"].get<std::string>();
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
                // Pedestals + common-mode ranges are loaded later from the
                // runinfo block (per-run calibration data).
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
            // Alpha-source trigger lives under lms_monitor.alpha because
            // Alpha events are always fed to the same LMS ref channels.
            if (lm.contains("alpha"))
                alpha_trigger.parse(lm["alpha"], trigger_bits_def);
            std::cerr << "LMS       : " << lms_trigger
                      << " warn=" << lms_warn_thresh
                      << " refs=" << lms_ref_channels.size()
                      << " alpha=" << alpha_trigger << "\n";
        }

        if (rcfg.contains("livetime")) {
            auto &lt = rcfg["livetime"];
            if (lt.contains("command"))  livetime_cmd      = lt["command"].get<std::string>();
            if (lt.contains("poll_sec")) livetime_poll_sec = std::max(1, (int)lt["poll_sec"]);
            if (lt.contains("healthy"))  livetime_healthy  = lt["healthy"];
            if (lt.contains("warning"))  livetime_warning  = lt["warning"];
            // DSC2-scaler config (bank/slot/source/channel) lives in
            // daq_config.json under `dsc_scaler` — see DaqConfig::DscScaler.
            const auto &ds = daq_cfg.dsc_scaler;
            using DSrc = evc::DaqConfig::DscScaler::Source;
            const char *src_name = (ds.source == DSrc::Ref) ? "ref"
                                 : (ds.source == DSrc::Trg) ? "trg" : "tdc";
            std::cerr << "Livetime  : "
                      << (livetime_cmd.empty() ? "disabled"
                                               : ("'" + livetime_cmd + "' every "
                                                  + std::to_string(livetime_poll_sec) + "s"))
                      << " healthy>=" << livetime_healthy
                      << " warn>=" << livetime_warning;
            if (ds.enabled()) {
                std::cerr << " | DSC2 bank=0x" << std::hex << ds.bank_tag << std::dec
                          << " slot=" << ds.slot << " src=" << src_name;
                if (ds.source != DSrc::Ref)
                    std::cerr << " ch=" << ds.channel;
            } else {
                std::cerr << " | DSC2 disabled";
            }
            std::cerr << "\n";
        }

        if (rcfg.contains("color_ranges")) {
            for (auto &[key, val] : rcfg["color_ranges"].items()) {
                if (val.is_array() && val.size() == 2)
                    color_range_defaults[key] = {val[0].get<float>(), val[1].get<float>()};
            }
            std::cerr << "Color ranges: " << color_range_defaults.size() << " entries\n";
        }

        if (rcfg.contains("runinfo") && rcfg["runinfo"].is_string()) {
            // Run number isn't known at init time (no event yet) — pick the
            // entry with the largest run_number ("latest"). The shared loader
            // logs which entry it selected.
            std::string ri_file = findFile(rcfg["runinfo"].get<std::string>(), db_dir);
            if (ri_file.empty()) {
                std::cerr << "Warning: runinfo file '"
                          << rcfg["runinfo"].get<std::string>()
                          << "' not found in " << db_dir << "\n";
            } else {
                prad2::RunConfig rc = prad2::LoadRunConfig(ri_file, /*run_num=*/-1);

                // Single-source beam energy: runinfo is the fallback. EPICS
                // (MBSY2C_energy) overrides at runtime in processEpics().
                beam_energy_runinfo = rc.Ebeam;
                beam_energy.store(rc.Ebeam);
                target_x = rc.target_x;
                target_y = rc.target_y;
                target_z = rc.target_z;
                adc_to_mev = rc.default_adc2mev;
                setTransform(hycal_transform,
                             rc.hycal_x, rc.hycal_y, rc.hycal_z,
                             rc.hycal_tilt_x, rc.hycal_tilt_y, rc.hycal_tilt_z);

                if (!rc.energy_calib_file.empty()) {
                    std::string calib_file = findFile(rc.energy_calib_file, db_dir);
                    if (calib_file.empty()) {
                        std::cerr << "Warning: calibration file '"
                                  << rc.energy_calib_file
                                  << "' not found in " << db_dir << "\n";
                    } else {
                        int nmatched = hycal.LoadCalibration(calib_file);
                        if (nmatched >= 0)
                            std::cerr << "Calibration: " << calib_file
                                      << " (" << nmatched << " modules)\n";
                    }
                }
                std::cerr << "RunInfo   : beam=" << beam_energy.load()
                          << "MeV default_adc2mev=" << adc_to_mev
                          << " target=(" << target_x << "," << target_y << "," << target_z
                          << ") HyCal=(" << hycal_transform.x << ","
                          << hycal_transform.y << "," << hycal_transform.z << ")\n";

                if (gem_enabled) {
                    int n = std::min<int>(4, (int)gem_transforms.size());
                    for (int id = 0; id < n; ++id) {
                        setTransform(gem_transforms[id],
                                     rc.gem_x[id], rc.gem_y[id], rc.gem_z[id],
                                     rc.gem_tilt_x[id], rc.gem_tilt_y[id], rc.gem_tilt_z[id]);
                    }
                    std::cerr << "GEM geom  : " << gem_transforms.size()
                              << " detectors configured\n";

                    // Build hardware-crate -> logical-crate remap from
                    // daq_config.roc_tags so the upstream pedestal/CM files
                    // (which key by EVIO bank tag = decimal hardware ID,
                    // e.g. 146/147) match our gem_map.json (logical 1/2).
                    std::map<int, int> gem_crate_remap;
                    for (const auto &re : daq_cfg.roc_tags) {
                        if (re.type == "gem")
                            gem_crate_remap[(int)re.tag] = re.crate;
                    }

                    if (!rc.gem_pedestal_file.empty()) {
                        std::string ped = findFile(rc.gem_pedestal_file, db_dir);
                        if (ped.empty())
                            std::cerr << "Warning: gem pedestal file '"
                                      << rc.gem_pedestal_file
                                      << "' not found in " << db_dir << "\n";
                        else
                            gem_sys.LoadPedestals(ped, gem_crate_remap);
                    }
                    if (!rc.gem_common_mode_file.empty()) {
                        std::string cm = findFile(rc.gem_common_mode_file, db_dir);
                        if (cm.empty())
                            std::cerr << "Warning: gem common-mode file '"
                                      << rc.gem_common_mode_file
                                      << "' not found in " << db_dir << "\n";
                        else
                            gem_sys.LoadCommonModeRange(cm, gem_crate_remap);
                    }
                }
            }
        } else if (rcfg.contains("runinfo")) {
            std::cerr << "Warning: 'runinfo' must be a path string to a "
                         "configurations-format JSON file\n";
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
            if (ph.contains("beam_energy")) {
                auto &be = ph["beam_energy"];
                if (be.contains("epics_channel")) beam_energy_epics_channel = be["epics_channel"];
                if (be.contains("min_valid"))     beam_energy_min_valid     = be["min_valid"];
            }
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
            }
            if (ph.contains("hycal_cluster_hit")) {
                auto &hc = ph["hycal_cluster_hit"];
                if (hc.contains("n_clusters"))      hxy_n_clusters      = hc["n_clusters"];
                if (hc.contains("energy_frac_min")) hxy_energy_frac_min = hc["energy_frac_min"];
                if (hc.contains("nblocks_min"))     hxy_nblocks_min     = hc["nblocks_min"];
                if (hc.contains("nblocks_max"))     hxy_nblocks_max     = hc["nblocks_max"];
                if (hc.contains("xy_hist")) {
                    auto &xy = hc["xy_hist"];
                    if (xy.contains("x_min"))  hxy_x_min  = xy["x_min"];
                    if (xy.contains("x_max"))  hxy_x_max  = xy["x_max"];
                    if (xy.contains("x_step")) hxy_x_step = xy["x_step"];
                    if (xy.contains("y_min"))  hxy_y_min  = xy["y_min"];
                    if (xy.contains("y_max"))  hxy_y_max  = xy["y_max"];
                    if (xy.contains("y_step")) hxy_y_step = xy["y_step"];
                }
            }
            std::cerr << "Physics   : " << physics_trigger
                      << " Moller: tol=" << moller_energy_tol
                      << " angle=[" << moller_angle_min << "," << moller_angle_max << "]"
                      << " HyCalXY: Ncl=" << hxy_n_clusters
                      << " E>=" << hxy_energy_frac_min << "*Eb"
                      << " blocks=[" << hxy_nblocks_min << "," << hxy_nblocks_max << "]"
                      << " beam_src='" << beam_energy_epics_channel << "'\n";
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

        // Top-level "gem" block — diagnostics for the GEM tab (HyCal-anchored
        // matching residuals + tracking-efficiency monitor).  Pulled out of
        // "physics" so the GEM tab owns its own configuration.
        if (rcfg.contains("gem")) {
            auto &gemcfg = rcfg["gem"];
            if (gemcfg.contains("hycal_match")) {
                auto &gm = gemcfg["hycal_match"];
                if (gm.contains("require_ep_candidate")) gem_match_require_ep = gm["require_ep_candidate"];
                if (gm.contains("window_mm"))            gem_match_window_mm  = gm["window_mm"];
                if (gm.contains("residual_hist")) {
                    auto &rh = gm["residual_hist"];
                    if (rh.contains("min"))  gem_resid_min  = rh["min"];
                    if (rh.contains("max"))  gem_resid_max  = rh["max"];
                    if (rh.contains("step")) gem_resid_step = rh["step"];
                }
            }
            if (gemcfg.contains("efficiency")) {
                auto &ge = gemcfg["efficiency"];
                if (ge.contains("min_cluster_energy"))    gem_eff_min_cluster_energy = ge["min_cluster_energy"];
                if (ge.contains("match_window_mm"))       gem_eff_match_window_mm    = ge["match_window_mm"];
                if (ge.contains("max_chi2_per_dof"))      gem_eff_max_chi2           = ge["max_chi2_per_dof"];
                if (ge.contains("max_hits_per_detector")) gem_eff_max_hits_per_det   = ge["max_hits_per_detector"];
                if (ge.contains("min_denom_for_eff"))     gem_eff_min_denom          = ge["min_denom_for_eff"];
                if (ge.contains("healthy"))               gem_eff_healthy            = ge["healthy"];
                if (ge.contains("warning"))               gem_eff_warning            = ge["warning"];
            }
            std::cerr << "GEM cfg   : match window=" << gem_match_window_mm
                      << "mm  efficiency: window=" << gem_eff_match_window_mm
                      << "mm chi2/dof<=" << gem_eff_max_chi2 << "\n";
        }

        std::cerr << "Reco      : " << main_config
                  << " (adc_to_mev=" << adc_to_mev << ")\n";
    }

    // init cluster histograms
    int cl_nbins = std::max(1, (int)std::ceil((cl_hist_max - cl_hist_min) / cl_hist_step));
    cluster_energy_hist.init(cl_nbins);
    int nb_nclusters = std::max(1, (int)std::ceil(
        (nclusters_hist_max - nclusters_hist_min) / nclusters_hist_step));
    nclusters_hist.init(nb_nclusters);
    int nb_blocks = std::max(1, (nblocks_hist_max - nblocks_hist_min) / nblocks_hist_step);
    nblocks_hist.init(nb_blocks);
    // One dependent histogram per Ncl bucket — sized to match the
    // unfiltered ones so the frontend can swap them in 1:1.
    cluster_energy_hist_by_ncl.assign(nb_nclusters, Histogram{});
    nblocks_hist_by_ncl.assign(nb_nclusters, Histogram{});
    for (auto &h : cluster_energy_hist_by_ncl) h.init(cl_nbins);
    for (auto &h : nblocks_hist_by_ncl)        h.init(nb_blocks);
    int ea_nx = std::max(1, (int)std::ceil((ea_angle_max - ea_angle_min) / ea_angle_step));
    int ea_ny = std::max(1, (int)std::ceil((ea_energy_max - ea_energy_min) / ea_energy_step));
    energy_angle_hist.init(ea_nx, ea_ny);
    int ml_nx = std::max(1, (int)std::ceil((moller_xy_x_max - moller_xy_x_min) / moller_xy_x_step));
    int ml_ny = std::max(1, (int)std::ceil((moller_xy_y_max - moller_xy_y_min) / moller_xy_y_step));
    moller_xy_hist.init(ml_nx, ml_ny);
    int hxy_nx = std::max(1, (int)std::ceil((hxy_x_max - hxy_x_min) / hxy_x_step));
    int hxy_ny = std::max(1, (int)std::ceil((hxy_y_max - hxy_y_min) / hxy_y_step));
    hycal_xy_hist.init(hxy_nx, hxy_ny);
    // GEM↔HyCal residual histograms — one per GEM detector for each axis.
    {
        int n_gem = gem_enabled ? (int)gem_transforms.size() : 0;
        int resid_nbins = std::max(1, (int)std::ceil((gem_resid_max - gem_resid_min) / gem_resid_step));
        gem_dx_hist.assign(n_gem, Histogram{});
        gem_dy_hist.assign(n_gem, Histogram{});
        gem_match_hits.assign(n_gem, 0);
        for (auto &h : gem_dx_hist) h.init(resid_nbins);
        for (auto &h : gem_dy_hist) h.init(resid_nbins);
    }
    initGemEfficiency();
    // ensure all transforms are pre-prepared before multithreaded use
    hycal_transform.prepare();
}
