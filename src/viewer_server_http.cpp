#include "viewer_server.h"
#include "http_compress.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

namespace fs = std::filesystem;
using json = nlohmann::json;

// =========================================================================
// Resource serving
// =========================================================================

bool ViewerServer::serveResource(const std::string &uri,
                                 WsServer::connection_ptr con)
{
    if (res_dir_.empty()) return false;
    std::string relpath = (uri == "/") ? "viewer.html" : uri.substr(1);
    if (relpath.find("..") != std::string::npos || relpath[0] == '/')
        return false;
    std::string fullpath = res_dir_ + "/" + relpath;
    std::string content = readFile(fullpath);
    if (content.empty()) return false;
    con->set_status(websocketpp::http::status_code::ok);
    con->set_body(content);
    con->append_header("Content-Type", contentType(fullpath));
    return true;
}

json ViewerServer::listFiles(const std::string &subdir)
{
    json entries = json::array();
    if (cfg_.data_dir.empty()) return entries;
    try {
        fs::path root(cfg_.data_dir);
        fs::path dir = subdir.empty() ? root : root / subdir;
        // security: ensure dir is under root
        auto canon_root = fs::canonical(root);
        auto canon_dir  = fs::canonical(dir);
        if (canon_dir.string().rfind(canon_root.string(), 0) != 0)
            return entries;

        for (auto &entry : fs::directory_iterator(
                 canon_dir, fs::directory_options::skip_permission_denied)) {
            auto rel = fs::relative(entry.path(), root).string();
            if (entry.is_directory()) {
                // count data files inside (non-recursive quick scan)
                int count = 0;
                try {
                    for (auto &child : fs::recursive_directory_iterator(
                             entry.path(), fs::directory_options::skip_permission_denied)) {
                        if (!child.is_regular_file()) continue;
                        auto fn = child.path().filename().string();
                        if (fn.find(".evio") != std::string::npos ||
                            fn.find(".root") != std::string::npos)
                            count++;
                    }
                } catch (...) {}
                if (count > 0)
                    entries.push_back(json{{"type", "dir"}, {"name", rel}, {"count", count}});
            } else if (entry.is_regular_file()) {
                auto fn = entry.path().filename().string();
                if (fn.find(".evio") == std::string::npos &&
                    fn.find(".root") == std::string::npos)
                    continue;
                auto sz = entry.file_size();
                entries.push_back(json{{"type", "file"}, {"name", rel},
                                       {"size", sz},
                                       {"size_mb", std::round(sz / 1048576.0 * 10) / 10}});
            }
        }
    } catch (...) {}
    std::sort(entries.begin(), entries.end(),
              [](const json &a, const json &b) {
                  // dirs first, then by name
                  bool da = a["type"] == "dir", db = b["type"] == "dir";
                  if (da != db) return da > db;
                  return a["name"] < b["name"];
              });
    return entries;
}

std::string ViewerServer::resolveDataFile(const std::string &relpath)
{
    if (cfg_.data_dir.empty()) return "";
    try {
        fs::path full = fs::canonical(fs::path(cfg_.data_dir) / relpath);
        fs::path root = fs::canonical(fs::path(cfg_.data_dir));
        if (full.string().rfind(root.string(), 0) != 0) return "";
        if (!fs::is_regular_file(full)) return "";
        return full.string();
    } catch (...) { return ""; }
}

// =========================================================================
// Config JSON
// =========================================================================

json ViewerServer::buildConfig()
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }

    auto &app = activeApp();
    json cfg = app.base_config;
    cfg["mode"] = mode();
#ifdef WITH_ET
    cfg["et_available"] = true;
    cfg["et_connected"] = et_connected_.load();
    cfg["ring_buffer_size"] = ring_size_;
    cfg["et_config"] = {
        {"host", et_cfg_.host}, {"port", et_cfg_.port},
        {"et_file", et_cfg_.et_file}, {"station", et_cfg_.station},
    };
#else
    cfg["et_available"] = false;
    cfg["et_connected"] = false;
    cfg["ring_buffer_size"] = 0;
    cfg["et_config"] = json::object();
#endif
    cfg["file_available"] = !cfg_.data_dir.empty() || (data != nullptr);
    cfg["total_events"] = data ? data->event_count : 0;
    cfg["current_file"] = data ? data->filepath : "";
    cfg["data_dir_enabled"] = !cfg_.data_dir.empty();
    cfg["data_dir"] = cfg_.data_dir;
    cfg["hist_enabled"] = (mode_.load() == Mode::Online) ? true : hist_enabled_.load();
    cfg["filter_active"] = app_file_.filterActive();
    cfg["filtered_count"] = filtered_indices_.empty() ? (data ? data->event_count : 0)
                                                       : (int)filtered_indices_.size();

    // data source capabilities
    // In online mode without a file, report EVIO-native capabilities
    DataSourceCaps caps;
    if (data) {
        caps = data->caps;
    } else if (mode_.load() == Mode::Online) {
        caps.source_type   = "evio";
        caps.has_waveforms = true;
        caps.has_peaks     = true;
        caps.has_pedestals = true;
        caps.has_epics     = true;
        caps.has_sync      = true;
    }
    cfg["source"] = {
        {"type", caps.source_type},
        {"has_waveforms", caps.has_waveforms},
        {"has_peaks", caps.has_peaks},
        {"has_pedestals", caps.has_pedestals},
        {"has_clusters", caps.has_clusters},
        {"has_gem_raw", caps.has_gem_raw},
        {"has_gem_hits", caps.has_gem_hits},
        {"has_epics", caps.has_epics},
        {"has_sync", caps.has_sync},
    };

    app.fillConfigJson(cfg);
    return cfg;
}

// =========================================================================
// Elog post
// =========================================================================

json ViewerServer::handleElogPost(const std::string &body)
{
    if (body.empty())
        return {{"ok", false}, {"error", "Empty body"}};
    if (activeApp().elog_url.empty())
        return {{"ok", false}, {"error", "No elog URL configured"}};

    auto req = json::parse(body, nullptr, false);
    if (req.is_discarded() || !req.contains("xml"))
        return {{"ok", false}, {"error", "Invalid request"}};

    std::string xml_body = req["xml"].get<std::string>();
    std::string tmp = "/tmp/prad2_elog_" + std::to_string(std::time(nullptr)) + ".xml";
    { std::ofstream f(tmp); f << xml_body; }

    std::string cert_flag;
    auto &app = activeApp();
    if (!app.elog_cert.empty())
        cert_flag = " --cert '" + app.elog_cert + "' --key '" + app.elog_key + "'";
    std::string cmd = "curl -s -o /dev/null -w '%{http_code}'" + cert_flag
                    + " --upload-file '" + tmp + "' '"
                    + app.elog_url + "/incoming/prad2_report.xml' 2>/dev/null";
    std::string http_code;
    FILE *p = popen(cmd.c_str(), "r");
    if (p) {
        char buf[256] = {};
        if (fgets(buf, sizeof(buf), p)) http_code = buf;
        while (!http_code.empty() &&
               (http_code.back() == '\n' || http_code.back() == '\r'))
            http_code.pop_back();
        pclose(p);
    }
    std::remove(tmp.c_str());
    bool ok = (http_code.find("200") != std::string::npos ||
               http_code.find("201") != std::string::npos);
    std::cerr << "Elog post: " << app.elog_url << " -> HTTP " << http_code
              << (ok ? " OK" : " FAIL") << "\n";
    return {{"ok", ok}, {"status", http_code}};
}

// =========================================================================
// HTTP handler
// =========================================================================

void ViewerServer::onHttp(WsServer *srv, websocketpp::connection_hdl hdl)
{
    auto con = srv->get_con_from_hdl(hdl);
    std::string uri = con->get_resource();

    // gzip if the client advertised it (browsers always do; the urllib-
    // based Python tools don't and keep getting plain bytes).  Read once
    // per request — cheap header lookup, but no point doing it twice.
    bool wants_gzip = prad2::client_accepts_gzip(
        con->get_request_header("Accept-Encoding"));

    // reply(body, [content_type], [pre_gz]) — when pre_gz is non-null and
    // the client accepts gzip, the cached compressed bytes are served
    // verbatim (saves re-deflating the same payload for each viewer).
    // Otherwise, large bodies are compressed on demand if the client
    // accepts gzip; small bodies and gzip-disabled clients get plain.
    auto reply = [&](const std::string &body,
                     const std::string &ct = "application/json",
                     const std::string *pre_gz = nullptr) {
        con->set_status(websocketpp::http::status_code::ok);
        con->append_header("Content-Type", ct);
        if (wants_gzip && pre_gz && !pre_gz->empty()) {
            con->set_body(*pre_gz);
            con->append_header("Content-Encoding", "gzip");
            return;
        }
        if (wants_gzip && body.size() >= prad2::kGzipMinBytes) {
            try {
                con->set_body(prad2::gzip_compress(body));
                con->append_header("Content-Encoding", "gzip");
                return;
            } catch (...) {
                // Fall through to plain body on any zlib failure.
            }
        }
        con->set_body(body);
    };

    // Any exception thrown below (e.g. std::stoi on a malformed %xx in the
    // URL, json type_error, std::bad_alloc) would otherwise unwind into
    // websocketpp / asio and terminate the server.  Convert to a 400/500
    // response and keep the io_context alive.
    try {

    // --- static resources ---
    if (serveResource(uri, con)) return;

    // --- config ---
    if (uri == "/api/config") { reply(buildConfig().dump()); return; }

    // --- runtime hist_cfg updates (time_min/time_max/threshold) ---
    // POST {time_min, time_max, threshold}.  Updates both file & online
    // AppStates, then re-clusters every cached online ring entry from the
    // stored raw EventData so the online tab refreshes instantly.  File
    // mode has no cache and recomputes per request.
    if (uri == "/api/hist_config") {
        std::string body = con->get_request_body();
        auto j = json::parse(body, nullptr, false);
        if (j.is_discarded() || !j.is_object()) {
            reply("{\"error\":\"invalid JSON\"}"); return;
        }
        auto applyTo = [&](AppState &app) {
            if (j.contains("time_min")  && j["time_min"].is_number())
                app.hist_cfg.time_min  = j["time_min"].get<float>();
            if (j.contains("time_max")  && j["time_max"].is_number())
                app.hist_cfg.time_max  = j["time_max"].get<float>();
            if (j.contains("threshold") && j["threshold"].is_number())
                app.hist_cfg.threshold = j["threshold"].get<float>();
        };
        applyTo(app_file_);
        applyTo(app_online_);
#ifdef WITH_ET
        // Re-cluster every cached ring entry under the new window so the
        // online tab updates instantly instead of waiting for fresh events.
        // json_str (raw peaks) is unaffected by hist_cfg, so we leave it.
        // Snapshot shared_ptrs under the lock, recompute outside, write back
        // under the lock — keeps the ET reader thread unblocked during the
        // heavy WaveAnalyzer pass.
        std::vector<std::pair<int, std::shared_ptr<fdec::EventData>>> snap;
        {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            snap.reserve(ring_.size());
            for (auto &e : ring_) snap.emplace_back(e.seq, e.event_data);
        }
        std::vector<std::pair<int, std::string>> updated;
        updated.reserve(snap.size());
        {
            fdec::WaveAnalyzer ana(app_online_.daq_cfg.wave_cfg);
            ana.cfg.min_peak_ratio = app_online_.hist_cfg.min_peak_ratio;
            ana.SetTemplateStore(&app_online_.template_store);
            fdec::WaveResult wres;
            for (auto &p : snap) {
                if (!p.second) continue;
                updated.emplace_back(p.first,
                    app_online_.computeClustersJson(
                        *p.second, p.first, ana, wres).dump());
            }
        }
        {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &u : updated) {
                for (auto &e : ring_) {
                    if (e.seq == u.first) { e.cluster_str = std::move(u.second); break; }
                }
            }
        }
#endif
        wsBroadcast(json({{"type", "hist_config_updated"},
                          {"time_min",  app_file_.hist_cfg.time_min},
                          {"time_max",  app_file_.hist_cfg.time_max},
                          {"threshold", app_file_.hist_cfg.threshold}}).dump());
        reply(json({{"ok", true},
                    {"time_min",  app_file_.hist_cfg.time_min},
                    {"time_max",  app_file_.hist_cfg.time_max},
                    {"threshold", app_file_.hist_cfg.threshold}}).dump());
        return;
    }

    // --- mode switching ---
    if (uri == "/api/mode/online") {
#ifdef WITH_ET
        {
            std::lock_guard<std::mutex> lk(mode_mtx_);
            // apply optional ET config overrides (serialised by mode_mtx_)
            std::string body = con->get_request_body();
            if (!body.empty()) {
                auto j = json::parse(body, nullptr, false);
                if (!j.is_discarded()) {
                    if (j.contains("host"))    et_cfg_.host    = j["host"];
                    if (j.contains("port"))    et_cfg_.port    = j["port"];
                    if (j.contains("et_file")) et_cfg_.et_file = j["et_file"];
                    if (j.contains("station")) et_cfg_.station = j["station"];
                }
            }
            // bump generation so ET reader reconnects with new config
            if (mode_.load() == Mode::Online)
                et_generation_++;
            et_active_ = true;
            setMode(Mode::Online);
        }
        reply(json({{"mode", "online"}}).dump());
#else
        reply(json({{"error", "ET support not compiled"}}).dump());
#endif
        return;
    }
    if (uri == "/api/mode/file") {
        std::lock_guard<std::mutex> lk(mode_mtx_);
        if (mode_.load() == Mode::Online) {
#ifdef WITH_ET
            et_active_ = false;
#endif
            std::shared_ptr<FileData> data;
            { std::lock_guard<std::mutex> lk2(file_data_mtx_); data = file_data_; }
            setMode(data ? Mode::File : Mode::Idle);
        }
        reply(json({{"mode", mode()}}).dump());
        return;
    }

    // --- event/latest (online only) ---
    if (uri == "/api/event/latest") {
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            if (!ring_.empty()) { reply(ring_.back().json_str); return; }
            reply("{\"error\":\"no events yet\"}"); return;
        }
#endif
        reply("{\"error\":\"not in online mode\"}"); return;
    }

    // --- event/<n> (mode-dependent) ---
    if (uri.rfind("/api/event/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 11);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &e : ring_) {
                if (e.seq == evnum) { reply(e.json_str); return; }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        reply(decodeEvent(evnum).dump()); return;
    }

    // --- waveform/<n>/<key> (file mode only — on-demand single-channel samples) ---
    if (uri.rfind("/api/waveform/", 0) == 0) {
        // parse /api/waveform/<evnum>/<roc_slot_ch>
        std::string rest = uri.substr(14);
        auto slash = rest.find('/');
        if (slash == std::string::npos) {
            reply("{\"error\":\"usage: /api/waveform/<event>/<roc_slot_ch>\"}"); return;
        }
        int evnum = std::atoi(rest.substr(0, slash).c_str());
        std::string chan_key = rest.substr(slash + 1);

        auto event_ptr = std::make_unique<fdec::EventData>();
        auto ssp_ptr = std::make_unique<ssp::SspEventData>();
        std::string err = decodeRawEvent(evnum, *event_ptr, ssp_ptr.get());
        if (!err.empty()) { reply(json({{"error", err}}).dump()); return; }

        accumulate(evnum, *event_ptr, ssp_ptr.get());

        fdec::WaveAnalyzer ana(activeApp().daq_cfg.wave_cfg);
        ana.cfg.min_peak_ratio = activeApp().hist_cfg.min_peak_ratio;
        ana.SetTemplateStore(&activeApp().template_store);
        fdec::WaveResult wres;
        reply(activeApp().encodeWaveformJson(*event_ptr, chan_key, ana, wres).dump());
        return;
    }

    // --- clusters/<n> (mode-dependent) ---
    if (uri.rfind("/api/clusters/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 14);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &e : ring_) {
                if (e.seq == evnum) {
                    if (e.cluster_str.empty()) {
                        // Cache was invalidated by /api/hist_config; we can't
                        // recompute here without the raw EventData, so report
                        // pending and let the next live event refill.
                        reply("{\"error\":\"clusters pending — config changed, "
                              "wait for next event\"}");
                    } else {
                        reply(e.cluster_str);
                    }
                    return;
                }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        reply(computeClusters(evnum).dump()); return;
    }

    // --- gem/calib (one-shot per-APV pedestal noise + global zs_sigma) ---
    // The frontend caches this and refetches only when the calib_rev
    // embedded in /api/gem/apv/<n> diverges from the cached value.
    if (uri == "/api/gem/calib") {
        reply(activeApp().apiGemCalib().dump());
        return;
    }

    // --- gem/threshold (POST {zs_sigma:N} — applies to all consumers) ---
    // Updates both file & online AppStates so the new threshold takes
    // effect for the active mode immediately and for the other mode on
    // the next event it processes.  Old ring entries keep their encoded
    // hits[] (cosmetic lag of ~ring_size events); each new event reflects
    // the new threshold via its zs_sigma field.  Broadcasts a WS notice
    // so other open viewers can refresh their toolbar input.
    if (uri == "/api/gem/threshold") {
        std::string body = con->get_request_body();
        auto j = json::parse(body, nullptr, false);
        if (j.is_discarded() || !j.is_object() ||
            !j.contains("zs_sigma") || !j["zs_sigma"].is_number()) {
            reply("{\"error\":\"expected {\\\"zs_sigma\\\":N}\"}"); return;
        }
        float new_sigma = j["zs_sigma"].get<float>();
        app_file_.setGemZsSigma(new_sigma);
        app_online_.setGemZsSigma(new_sigma);
        wsBroadcast(json({{"type", "gem_threshold_updated"},
                          {"zs_sigma", new_sigma}}).dump());
        reply(json({{"ok", true},
                    {"zs_sigma", new_sigma}}).dump());
        return;
    }

    // --- gem/apv/<n> (per-event GEM APV waveforms, mode-dependent) ---
    // Online: served from a per-ring-entry pre-encoded string so older
    // events don't disturb the live gem_sys state.
    // File: decode + accumulate (which fills gem_sys), then build JSON.
    if (uri.rfind("/api/gem/apv/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 13);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &e : ring_) {
                if (e.seq == evnum) {
                    if (e.gem_apv_str.empty())
                        reply("{\"error\":\"gem apv pending\"}");
                    else
                        reply(e.gem_apv_str, "application/json", &e.gem_apv_gz);
                    return;
                }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        auto event_ptr = std::make_unique<fdec::EventData>();
        auto ssp_ptr   = std::make_unique<ssp::SspEventData>();
        std::string err = decodeRawEvent(evnum, *event_ptr, ssp_ptr.get());
        if (!err.empty()) { reply(json({{"error", err}}).dump()); return; }
        accumulate(evnum, *event_ptr, ssp_ptr.get());
        // accumulate() dedupes by event id, so on a re-request gem_sys
        // may still hold a different event's working buffers.  Force a
        // re-process (no histogram side effects) so apiGemApv reads the
        // requested event regardless of cache state.
        activeApp().prepareGemForView(*ssp_ptr);
        reply(activeApp().apiGemApv(*ssp_ptr, evnum).dump());
        return;
    }

    // --- ring buffer ---
    if (uri == "/api/ring") {
#ifdef WITH_ET
        std::lock_guard<std::mutex> lk(ring_mtx_);
        json arr = json::array();
        for (auto &e : ring_) arr.push_back(e.seq);
        reply(json({{"ring", arr},
                     {"latest", ring_.empty() ? 0 : ring_.back().seq}}).dump());
#else
        reply(json({{"ring", json::array()}, {"latest", 0}}).dump());
#endif
        return;
    }

    // --- progress ---
    if (uri == "/api/progress") { reply(progress_.toJson().dump()); return; }

    // LIVETIME — TS (caget poll) and Meas. (DSC2 from data stream).  Either
    // may be <0 to mean "not available"; the frontend hides each value
    // independently.  The DSC2 path is read from activeApp() so it follows
    // file/online mode without extra state.
    if (uri == "/api/livetime") {
        double ts = -1.0;
        double meas = activeApp().measured_livetime.load();
#ifdef WITH_ET
        ts = livetime_.load();
#endif
        reply(json({{"livetime", ts}, {"measured", meas}}).dump());
        return;
    }

    // --- clear endpoints (always available, clears active mode's data) ---
    // --- filter endpoints ---
    if (uri == "/api/filter") {
        reply(activeApp().filterToJson().dump()); return;
    }
    if (uri == "/api/filter/load") {
        std::string body = con->get_request_body();
        auto fj = json::parse(body, nullptr, false);
        if (fj.is_discarded()) { reply("{\"error\":\"invalid JSON\"}"); return; }
        std::string err = applyFilter(fj);
        if (!err.empty()) { reply(json({{"error", err}}).dump()); return; }
        reply(json({{"status", "ok"}, {"filter", activeApp().filterToJson()},
                     {"filtered_count", (int)filtered_indices_.size()}}).dump());
        return;
    }
    if (uri == "/api/filter/unload") {
        clearFilter();
        reply(json({{"status", "ok"}, {"filter_active", false}}).dump());
        return;
    }
    if (uri == "/api/filter/indices") {
        reply(json(filtered_indices_).dump()); return;
    }

    // --- clear endpoints (always available, clears active mode's data) ---
    if (uri == "/api/hist/clear") {
        activeApp().clearHistograms();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"hist_cleared\"}");
        return;
    }
    if (uri == "/api/lms/clear") {
        activeApp().clearLms();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"lms_cleared\"}");
        return;
    }
    if (uri == "/api/epics/clear") {
        activeApp().clearEpics();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"epics_cleared\"}");
        return;
    }

    // --- shared read-only API routes ---
    auto result = activeApp().handleReadApi(uri);
    if (result.handled) { reply(result.body); return; }

    // --- file browser ---
    if (uri == "/api/files" || uri.rfind("/api/files?", 0) == 0) {
        std::string subdir;
        auto qpos = uri.find('?');
        if (qpos != std::string::npos) {
            std::string q = uri.substr(qpos + 1);
            if (q.rfind("dir=", 0) == 0) {
                subdir = q.substr(4);
                // URL-decode
                std::string dec;
                for (size_t i = 0; i < subdir.size(); ++i) {
                    if (subdir[i] == '%' && i + 2 < subdir.size()) {
                        dec += (char)std::stoi(subdir.substr(i + 1, 2), nullptr, 16);
                        i += 2;
                    } else if (subdir[i] == '+') dec += ' ';
                    else dec += subdir[i];
                }
                subdir = dec;
            }
        }
        reply(json({{"entries", listFiles(subdir)}}).dump()); return;
    }

    // --- load file (relative path from data_dir) ---
    if (uri.rfind("/api/load?", 0) == 0) {
        auto qpos = uri.find('?');
        std::string query = uri.substr(qpos + 1);
        std::string relpath;
        bool do_hist = false;
        for (size_t pos = 0; pos < query.size();) {
            size_t amp = query.find('&', pos);
            if (amp == std::string::npos) amp = query.size();
            std::string kv = query.substr(pos, amp - pos);
            auto eq = kv.find('=');
            if (eq != std::string::npos) {
                std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
                if (k == "file") relpath = v;
                if (k == "hist") do_hist = (v == "1");
            }
            pos = amp + 1;
        }
        // URL-decode %xx
        std::string decoded;
        for (size_t i = 0; i < relpath.size(); ++i) {
            if (relpath[i] == '%' && i + 2 < relpath.size()) {
                decoded += (char)std::stoi(relpath.substr(i + 1, 2), nullptr, 16);
                i += 2;
            } else if (relpath[i] == '+') decoded += ' ';
            else decoded += relpath[i];
        }
        relpath = decoded;

        std::string fullpath = resolveDataFile(relpath);
        if (fullpath.empty()) { reply("{\"error\":\"invalid path\"}"); return; }

        loadFile(fullpath, do_hist);
        reply(json({{"status", "loading"}, {"file", relpath},
                     {"hist_enabled", do_hist}}).dump());
        return;
    }

    // --- elog post ---
    if (uri == "/api/elog/post") {
        std::string body = con->get_request_body();
        json result = handleElogPost(body);
        if (!result.value("ok", true)) {
            con->set_status(result.contains("error") &&
                            result["error"] == "Empty body"
                            ? websocketpp::http::status_code::bad_request
                            : websocketpp::http::status_code::ok);
        } else {
            con->set_status(websocketpp::http::status_code::ok);
        }
        con->set_body(result.dump());
        con->append_header("Content-Type", "application/json");
        return;
    }

    // --- 404 ---
    con->set_status(websocketpp::http::status_code::not_found);
    con->set_body("404 Not Found");

    } catch (const std::exception &e) {
        con->set_status(websocketpp::http::status_code::bad_request);
        con->set_body(json({{"error", e.what()}, {"uri", uri}}).dump());
        con->append_header("Content-Type", "application/json");
        std::cerr << "[http] " << uri << " → " << e.what() << "\n";
    } catch (...) {
        con->set_status(websocketpp::http::status_code::internal_server_error);
        con->set_body("{\"error\":\"unknown exception\"}");
        con->append_header("Content-Type", "application/json");
        std::cerr << "[http] " << uri << " → unknown exception\n";
    }
}
