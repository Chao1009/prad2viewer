#include "viewer_server.h"

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

    auto reply = [&](const std::string &body,
                     const std::string &ct = "application/json") {
        con->set_status(websocketpp::http::status_code::ok);
        con->set_body(body);
        con->append_header("Content-Type", ct);
    };

    // --- static resources ---
    if (serveResource(uri, con)) return;

    // --- config ---
    if (uri == "/api/config") { reply(buildConfig().dump()); return; }

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

        fdec::WaveAnalyzer ana;
        ana.cfg.min_peak_ratio = activeApp().hist_cfg.min_peak_ratio;
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
                if (e.seq == evnum) { reply(e.cluster_str); return; }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        reply(computeClusters(evnum).dump()); return;
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

    // LIVETIME — temporary
    if (uri == "/api/livetime") {
#ifdef WITH_ET
        double v = livetime_.load();
        reply(json({{"livetime", v}}).dump());
#else
        reply(json({{"livetime", -1}}).dump());
#endif
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
}
