#include "app_state.h"

#include <cmath>

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
        nclusters_hist_min, nclusters_hist_max, nclusters_hist_step);
    r["nblocks"] = histToJson(nblocks_hist,
        (float)nblocks_hist_min, (float)nblocks_hist_max, (float)nblocks_hist_step);
    // Per-Ncl bucket dependent histograms — bins_by_ncl[i] is the
    // bins array of the i-th bucket (same indexing as nclusters_hist).
    // Frontend uses these to redraw the energy / blocks histos when the
    // user clicks a particular Ncl bar.
    json energy_by_ncl = json::array();
    json blocks_by_ncl = json::array();
    for (auto &h : cluster_energy_hist_by_ncl) energy_by_ncl.push_back(h.bins);
    for (auto &h : nblocks_hist_by_ncl)        blocks_by_ncl.push_back(h.bins);
    r["bins_by_ncl"]            = energy_by_ncl;
    r["nblocks"]["bins_by_ncl"] = blocks_by_ncl;
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
            {"beam_energy", beam_energy.load()},
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

json AppState::apiHycalXY() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {{"xy_bins", hycal_xy_hist.bins},
            {"xy_nx", hycal_xy_hist.nx}, {"xy_ny", hycal_xy_hist.ny},
            {"xy_x_min", hxy_x_min}, {"xy_x_max", hxy_x_max}, {"xy_x_step", hxy_x_step},
            {"xy_y_min", hxy_y_min}, {"xy_y_max", hxy_y_max}, {"xy_y_step", hxy_y_step},
            {"events", hycal_xy_events},
            {"total_events", cluster_events_processed},
            {"beam_energy", beam_energy.load()},
            {"cuts", {{"n_clusters", hxy_n_clusters},
                      {"energy_frac_min", hxy_energy_frac_min},
                      {"nblocks_min", hxy_nblocks_min},
                      {"nblocks_max", hxy_nblocks_max}}}};
}

json AppState::apiGemResiduals() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    json dets = json::array();
    int n = (int)gem_dx_hist.size();
    int n_dets_runtime = std::min(n, gem_sys.GetNDetectors());
    for (int d = 0; d < n; ++d) {
        std::string name = (d < n_dets_runtime)
            ? gem_sys.GetDetectors()[d].name
            : ("GEM" + std::to_string(d));
        dets.push_back({
            {"id", d},
            {"name", name},
            {"dx_hist", histToJson(gem_dx_hist[d], gem_resid_min, gem_resid_max, gem_resid_step)},
            {"dy_hist", histToJson(gem_dy_hist[d], gem_resid_min, gem_resid_max, gem_resid_step)},
            {"matched_hits", gem_match_hits[d]},
        });
    }
    return {{"enabled", gem_enabled},
            {"detectors", dets},
            {"events", gem_match_events},
            {"cuts", {{"window_mm", gem_match_window_mm},
                      {"require_ep_candidate", gem_match_require_ep}}}};
}

json AppState::apiOccupancy() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    json jocc = json::object(), jtcut = json::object();
    for (auto &[k,v] : occupancy) jocc[k] = v;
    for (auto &[k,v] : occupancy_tcut) jtcut[k] = v;
    return {{"occ", jocc}, {"occ_tcut", jtcut}, {"total", events_processed.load()}};
}

// Reference correction: scalar factor = Alpha_latest / LMS_latest for the chosen
// ref channel.  Both readings are the most recent of their respective triggers
// (Am-241 alpha source vs. LMS pulser).  Applied uniformly to all history
// entries: corrected = signal * factor.  The Alpha source is stable, so this
// removes the LMS pulser's own drift while leaving real channel-gain changes.
struct RefCorrection {
    float factor = 1.f;
    float lms = 0.f;          // latest LMS integral on the ref module (for telemetry)
    float alpha = 0.f;        // latest Alpha integral on the ref module
    bool  active = false;
};

static RefCorrection buildRefCorrection(
    const std::map<int, float> &latest_lms,
    const std::map<int, float> &latest_alpha,
    const std::vector<AppState::LmsRefChannel> &refs, int ref_index)
{
    RefCorrection rc;
    if (ref_index < 0 || ref_index >= static_cast<int>(refs.size())) return rc;
    int ri = refs[ref_index].module_index;
    if (ri < 0) return rc;
    auto lit = latest_lms.find(ri);
    auto ait = latest_alpha.find(ri);
    if (lit == latest_lms.end() || ait == latest_alpha.end()) return rc;
    if (lit->second <= 0 || ait->second <= 0) return rc;
    rc.lms    = lit->second;
    rc.alpha  = ait->second;
    rc.factor = rc.alpha / rc.lms;
    rc.active = true;
    return rc;
}

// Apply correction: returns signal * factor (uniform across history).
static float applyRefCorrection(float val, const RefCorrection &rc)
{
    if (!rc.active) return val;
    return val * rc.factor;
}

json AppState::apiLmsSummary(int ref_index) const
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    auto rc = buildRefCorrection(latest_lms_integral, latest_alpha_integral,
                                  lms_ref_channels, ref_index);

    json mods = json::object();
    for (auto &[idx, hist] : lms_history) {
        if (hist.empty()) continue;
        double sum = 0, sum2 = 0;
        int count = 0;
        for (auto &e : hist) {
            float v = applyRefCorrection(e.integral, rc);
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
            {"ref_factor", rc.factor},
            {"ref_lms", rc.lms},
            {"ref_alpha", rc.alpha},
            {"sync_unix", sync_unix}, {"sync_rel_sec", sync_rel_sec}};
}

json AppState::apiLmsModule(int mod_idx, int ref_index) const
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    auto it = lms_history.find(mod_idx);
    if (it == lms_history.end() || it->second.empty())
        return {{"time", json::array()}, {"integral", json::array()}, {"events", 0}};

    auto rc = buildRefCorrection(latest_lms_integral, latest_alpha_integral,
                                  lms_ref_channels, ref_index);

    auto &hist = it->second;
    json t_arr = json::array(), v_arr = json::array();
    for (auto &e : hist) {
        float v = applyRefCorrection(e.integral, rc);
        t_arr.push_back(std::round(e.time_sec * 100) / 100);
        v_arr.push_back(std::round(v * 10) / 10);
    }
    std::string name = (mod_idx >= 0 && mod_idx < hycal.module_count())
        ? hycal.module(mod_idx).name : "";
    return {{"name", name}, {"time", t_arr}, {"integral", v_arr},
            {"events", (int)t_arr.size()},
            {"ref_index", ref_index},
            {"ref_factor", rc.factor},
            {"ref_lms", rc.lms},
            {"ref_alpha", rc.alpha},
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

    // Single-source beam energy: latest valid MBSY2C_energy reading overrides the
    // runinfo fallback. Skip values below min_valid (zero/garbage during beam trips).
    if (!beam_energy_epics_channel.empty()) {
        int id = epics.GetChannelId(beam_energy_epics_channel);
        if (id >= 0) {
            int n = epics.GetSnapshotCount();
            if (n > 0) {
                const auto &snap = epics.GetSnapshot(n - 1);
                if (id < (int)snap.values.size()) {
                    float v = snap.values[id];
                    if (v > beam_energy_min_valid) beam_energy.store(v);
                }
            }
        }
    }
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
    cfg["livetime"] = {
        {"enabled", !livetime_cmd.empty()},
        {"poll_sec", livetime_poll_sec},
        {"healthy", livetime_healthy},
        {"warning", livetime_warning},
    };
    cfg["runinfo"] = {
        {"beam_energy", beam_energy.load()},
        {"beam_energy_runinfo", beam_energy_runinfo},
        {"calibration", {{"default_adc2mev", adc_to_mev}}},
        {"target", {target_x, target_y, target_z}},
        {"hycal", {
            {"position", {hycal_transform.x, hycal_transform.y, hycal_transform.z}},
            {"tilting", {hycal_transform.rx, hycal_transform.ry, hycal_transform.rz}},
        }},
    };
    cfg["physics"] = {
        {"trigger", physics_trigger.toJson()},
        {"beam_energy", {
            {"epics_channel", beam_energy_epics_channel},
            {"min_valid", beam_energy_min_valid},
        }},
        {"energy_angle_hist", {
            {"angle_min", ea_angle_min}, {"angle_max", ea_angle_max}, {"angle_step", ea_angle_step},
            {"energy_min", ea_energy_min}, {"energy_max", ea_energy_max}, {"energy_step", ea_energy_step},
        }},
        {"moller", {
            {"energy_tolerance", moller_energy_tol},
            {"angle_min", moller_angle_min}, {"angle_max", moller_angle_max},
        }},
        {"hycal_cluster_hit", {
            {"n_clusters", hxy_n_clusters},
            {"energy_frac_min", hxy_energy_frac_min},
            {"nblocks_min", hxy_nblocks_min},
            {"nblocks_max", hxy_nblocks_max},
        }},
        {"gem_hycal_match", {
            {"require_ep_candidate", gem_match_require_ep},
            {"window_mm", gem_match_window_mm},
            {"residual_hist", {
                {"min", gem_resid_min}, {"max", gem_resid_max}, {"step", gem_resid_step},
            }},
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
    if (uri == "/api/physics/hycal_xy")
        return {true, apiHycalXY().dump()};
    if (uri == "/api/gem/residuals")
        return {true, apiGemResiduals().dump()};
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

