// src/evc_viewer.cpp — HyCal event viewer
//
// Usage:
//   evc_viewer [evio_file] [-p port] [-H] [-c config.json] [-d /path/to/data] [-D daq_config.json]
//   evc_viewer -d /data/stage6 -H              # browse and pick from GUI
//   evc_viewer data.evio -H                    # open file directly
//   evc_viewer data.evio -H -d /data/stage6    # open file + enable browsing
//   evc_viewer                                 # empty viewer, no file browser
//   evc_viewer prad.evio -D prad_daq_config.json  # open PRad file with PRad DAQ config
//
// -d enables file browsing: the viewer shows a file picker limited to
// .evio files under that directory tree. Selecting a new file triggers
// background re-indexing + histogram building with progress updates.

#include "EvChannel.h"
#include "HyCalCluster.h"
#include "app_state.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <csignal>
#include <getopt.h>

using json = nlohmann::json;
using WsServer = websocketpp::server<websocketpp::config::asio>;
namespace fs = std::filesystem;
using namespace evc;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif
#ifndef RESOURCE_DIR
#define RESOURCE_DIR "."
#endif

// -------------------------------------------------------------------------
// Viewer-specific types
// -------------------------------------------------------------------------
struct EventIndex { int buffer_num, sub_event; };

struct FileData {
    std::string filepath;
    std::vector<EventIndex> index;
};

struct Progress {
    std::atomic<bool> loading{false};
    std::atomic<int>  phase{0};
    std::atomic<int>  current{0};
    std::atomic<int>  total{0};
    std::string       target_file;
    std::mutex        mtx;

    json toJson() {
        std::lock_guard<std::mutex> lk(mtx);
        return {{"loading", loading.load()},
                {"phase", phase == 1 ? "indexing" : phase == 2 ? "histograms" : "idle"},
                {"current", current.load()}, {"total", total.load()}, {"file", target_file}};
    }
    void setFile(const std::string &f) {
        std::lock_guard<std::mutex> lk(mtx);
        target_file = f;
    }
};

// -------------------------------------------------------------------------
// Globals
// -------------------------------------------------------------------------
static AppState g_app;                        // shared state
static std::shared_ptr<FileData> g_data;
static std::mutex g_data_mtx;
static std::thread g_load_thread;
static std::mutex g_load_mtx;
static bool g_hist_enabled = false;
static std::string g_data_dir;
static std::string g_res_dir;
static Progress g_progress;

// cached file reader
static struct CachedReader {
    EvChannel ch;
    std::string filepath;
    int current_buf = 0;
    std::mutex mtx;

    std::string seekTo(const std::string &path, int buf_num) {
        if (path != filepath || buf_num < current_buf) {
            ch.Close();
            ch.SetConfig(g_app.daq_cfg);
            if (ch.Open(path) != status::success) {
                filepath.clear(); current_buf = 0;
                return "cannot open file";
            }
            filepath = path;
            current_buf = 0;
        }
        while (current_buf < buf_num) {
            if (ch.Read() != status::success) {
                ch.Close(); filepath.clear(); current_buf = 0;
                return "read error";
            }
            current_buf++;
        }
        return "";
    }
    void invalidate() { ch.Close(); filepath.clear(); current_buf = 0; }
} g_reader;

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------
static bool serveResource(const std::string &uri, WsServer::connection_ptr con) {
    if (g_res_dir.empty()) return false;
    std::string relpath = (uri == "/") ? "viewer.html" : uri.substr(1);
    if (relpath.find("..") != std::string::npos || relpath[0] == '/') return false;
    std::string fullpath = g_res_dir + "/" + relpath;
    std::string content = readFile(fullpath);
    if (content.empty()) return false;
    con->set_status(websocketpp::http::status_code::ok);
    con->set_body(content);
    con->append_header("Content-Type", contentType(fullpath));
    return true;
}

static json listFiles(const std::string &data_dir) {
    json files = json::array();
    if (data_dir.empty()) return files;
    try {
        fs::path root(data_dir);
        for (auto &entry : fs::recursive_directory_iterator(root, fs::directory_options::skip_permission_denied)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().filename().string().find(".evio") == std::string::npos) continue;
            auto rel = fs::relative(entry.path(), root).string();
            auto sz = entry.file_size();
            files.push_back({{"path", rel}, {"size", sz}, {"size_mb", std::round(sz / 1048576.0 * 10) / 10}});
        }
    } catch (...) {}
    std::sort(files.begin(), files.end(), [](const json &a, const json &b) { return a["path"] < b["path"]; });
    return files;
}

static std::string resolveDataFile(const std::string &relpath) {
    if (g_data_dir.empty()) return "";
    try {
        fs::path full = fs::canonical(fs::path(g_data_dir) / relpath);
        fs::path root = fs::canonical(fs::path(g_data_dir));
        if (full.string().rfind(root.string(), 0) != 0) return "";
        if (!fs::is_regular_file(full)) return "";
        return full.string();
    } catch (...) { return ""; }
}

// -------------------------------------------------------------------------
// Build index
// -------------------------------------------------------------------------
static void buildIndex(const std::string &path, std::vector<EventIndex> &index, Progress &prog) {
    index.clear();
    EvChannel ch;
    ch.SetConfig(g_app.daq_cfg);
    if (ch.Open(path) != status::success) return;
    prog.phase = 1; prog.current = 0;
    int buf = 0;
    while (ch.Read() == status::success) {
        ++buf; prog.current = buf;
        if (!ch.Scan()) continue;
        for (int i = 0; i < ch.GetNEvents(); ++i)
            index.push_back({buf, i});
    }
    ch.Close();
    prog.total = (int)index.size();
}

// -------------------------------------------------------------------------
// Build histograms (delegates to AppState)
// -------------------------------------------------------------------------
static void buildHistograms(const std::string &path, Progress &prog) {
    g_app.clearHistograms();
    g_app.clearLms();
    g_app.clearEpics();

    EvChannel ch;
    ch.SetConfig(g_app.daq_cfg);
    if (ch.Open(path) != status::success) return;

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_app.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;

    prog.phase = 2; prog.current = 0;
    int buf = 0;
    uint64_t last_ti_ts = 0;
    while (ch.Read() == status::success) {
        ++buf;
        if (!ch.Scan()) continue;
        // capture sync time reference
        if (g_app.sync_unix == 0) {
            uint32_t ct = ch.GetControlTime();
            if (ct != 0) g_app.recordSyncTime(ct, last_ti_ts);
        }
        // EPICS slow control event
        if (ch.GetEventType() == evc::EventType::Epics) {
            std::string text = ch.ExtractEpicsText();
            if (!text.empty())
                g_app.processEpics(text, g_app.events_processed.load(), last_ti_ts);
        }
        for (int i = 0; i < ch.GetNEvents(); ++i) {
            if (!ch.DecodeEvent(i, event)) continue;
            prog.current = g_app.events_processed.load() + 1;
            last_ti_ts = event.info.timestamp;

            // process: histograms + clustering + LMS (single-threaded, no locks needed)
            g_app.fillHist(event, ana, wres);
            g_app.clusterEvent(event, ana, wres);
            g_app.processLms(event, ana, wres);
            g_app.events_processed++;
        }
    }
    ch.Close();
}

// -------------------------------------------------------------------------
// Load file async
// -------------------------------------------------------------------------
static void loadFileAsync(const std::string &filepath) {
    g_progress.loading = true;
    g_progress.setFile(filepath);
    g_progress.phase = 0; g_progress.current = 0; g_progress.total = 0;

    auto data = std::make_shared<FileData>();
    data->filepath = filepath;

    std::cerr << "Loading: " << filepath << "\n";
    buildIndex(filepath, data->index, g_progress);
    std::cerr << "  Indexed " << data->index.size() << " events\n";

    if (g_hist_enabled) {
        g_progress.total = (int)data->index.size();
        buildHistograms(filepath, g_progress);
        std::cerr << "  Histograms: " << g_app.events_processed.load() << " events"
                  << ", clusters: " << g_app.cluster_events_processed
                  << ", LMS: " << g_app.lms_events.load() << "\n";
    }

    { std::lock_guard<std::mutex> lk(g_data_mtx); g_data = data; }
    { std::lock_guard<std::mutex> lk(g_reader.mtx); g_reader.invalidate(); }

    g_progress.loading = false;
    g_progress.phase = 0;
    std::cerr << "  Ready\n";
}

// -------------------------------------------------------------------------
// Decode raw event from file
// -------------------------------------------------------------------------
static std::string decodeRawEvent(int ev1, fdec::EventData &event) {
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(g_data_mtx); data = g_data; }
    if (!data) return "no file loaded";
    int idx = ev1 - 1;
    if (idx < 0 || idx >= (int)data->index.size()) return "event out of range";

    auto &ei = data->index[idx];
    std::lock_guard<std::mutex> lk(g_reader.mtx);
    std::string err = g_reader.seekTo(data->filepath, ei.buffer_num);
    if (!err.empty()) return err;
    if (!g_reader.ch.Scan()) return "scan error";
    if (!g_reader.ch.DecodeEvent(ei.sub_event, event)) return "decode error";
    return "";
}

// -------------------------------------------------------------------------
// Decode one event → JSON (delegates to AppState)
// -------------------------------------------------------------------------
static json decodeEvent(int ev1) {
    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    std::string err = decodeRawEvent(ev1, event);
    if (!err.empty()) return {{"error", err}};

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_app.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    return g_app.encodeEventJson(event, ev1, ana, wres);
}

// -------------------------------------------------------------------------
// Compute clusters for one event (delegates to AppState)
// -------------------------------------------------------------------------
static json computeClusters(int ev1) {
    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    std::string err = decodeRawEvent(ev1, event);
    if (!err.empty()) return {{"error", err}};

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_app.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    return g_app.computeClustersJson(event, ev1, ana, wres);
}

// -------------------------------------------------------------------------
// Build config JSON
// -------------------------------------------------------------------------
static json buildConfig() {
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(g_data_mtx); data = g_data; }

    json cfg = g_app.base_config;
    cfg["total_events"] = data ? (int)data->index.size() : 0;
    cfg["current_file"] = data ? data->filepath : "";
    cfg["data_dir_enabled"] = !g_data_dir.empty();
    cfg["data_dir"] = g_data_dir;
    cfg["hist_enabled"] = g_hist_enabled;
    cfg["mode"] = "file";
    g_app.fillConfigJson(cfg);
    return cfg;
}

// -------------------------------------------------------------------------
// HTTP handler
// -------------------------------------------------------------------------
static void onHttp(WsServer *srv, websocketpp::connection_hdl hdl) {
    auto con = srv->get_con_from_hdl(hdl);
    std::string uri = con->get_resource();

    auto reply = [&](const std::string &body, const std::string &ct = "application/json") {
        con->set_status(websocketpp::http::status_code::ok);
        con->set_body(body);
        con->append_header("Content-Type", ct);
    };

    if (serveResource(uri, con)) return;

    if (uri == "/api/config") { reply(buildConfig().dump()); return; }

    // /api/event/<num>
    if (uri.rfind("/api/event/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 11);
        reply(decodeEvent(evnum).dump()); return;
    }

    // /api/clusters/<num>
    if (uri.rfind("/api/clusters/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 14);
        reply(computeClusters(evnum).dump()); return;
    }

    // /api/progress
    if (uri == "/api/progress") { reply(g_progress.toJson().dump()); return; }

    // shared read-only API routes
    auto result = g_app.handleReadApi(uri);
    if (result.handled) { reply(result.body); return; }

    // /api/files
    if (uri == "/api/files") { reply(json({{"files", listFiles(g_data_dir)}}).dump()); return; }

    // /api/load?file=<relpath>&hist=0|1
    if (uri.rfind("/api/load?", 0) == 0) {
        auto qpos = uri.find('?');
        std::string query = uri.substr(qpos + 1);
        // parse file= and hist= from query
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

        g_hist_enabled = do_hist;
        {
            std::lock_guard<std::mutex> lk(g_load_mtx);
            if (g_load_thread.joinable()) g_load_thread.join();
            g_load_thread = std::thread([fullpath]() { loadFileAsync(fullpath); });
        }
        reply(json({{"status", "loading"}, {"file", relpath},
                     {"hist_enabled", g_hist_enabled}}).dump());
        return;
    }

    con->set_status(websocketpp::http::status_code::not_found);
    con->set_body("404 Not Found");
}

// -------------------------------------------------------------------------
// Main
// -------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    std::string evio_file;
    int port = 5050;
    std::string config_file;
    std::string daq_config_file;

    static struct option long_opts[] = {
        {"port",        required_argument, nullptr, 'p'},
        {"hist",        no_argument,       nullptr, 'H'},
        {"config",      required_argument, nullptr, 'c'},
        {"data-dir",    required_argument, nullptr, 'd'},
        {"daq-config",  required_argument, nullptr, 'D'},
        {"help",        no_argument,       nullptr, '?'},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "p:Hc:d:D:", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'p': port = std::atoi(optarg); break;
        case 'H': g_hist_enabled = true; break;
        case 'c': config_file = optarg; break;
        case 'd': g_data_dir = optarg; break;
        case 'D': daq_config_file = optarg; break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " [evio_file] [-p port] [-H] [-c config.json]"
                      << " [-d data_dir] [-D daq_config.json]\n";
            return 1;
        }
    }
    if (optind < argc) evio_file = argv[optind];

    std::string db_dir  = DATABASE_DIR;
    std::string res_dir = RESOURCE_DIR;

    // initialize shared state
    g_app.init(db_dir, daq_config_file, config_file);

    // build base_config for /api/config
    std::string mod_file = findFile("hycal_modules.json", db_dir);
    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(mod_file); if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    {
        std::string daq_fn = "daq_map.json";
        if (!daq_config_file.empty()) {
            std::ifstream dcf(daq_config_file);
            if (dcf.is_open()) {
                auto dcj = json::parse(dcf, nullptr, false, true);
                if (dcj.contains("daq_map_file")) daq_fn = dcj["daq_map_file"].get<std::string>();
            }
        }
        std::string s = readFile(findFile(daq_fn, db_dir));
        if (!s.empty()) daq_j = json::parse(s, nullptr, false);
    }
    g_app.base_config = {
        {"modules", modules_j}, {"daq", daq_j}, {"crate_roc", g_app.crate_roc_json},
    };

    // resources
    g_res_dir = res_dir;
    if (readFile(g_res_dir + "/viewer.html").empty())
        std::cerr << "Warning: viewer.html not found in " << g_res_dir << "\n";

    std::cerr << "Database  : " << db_dir << "\n"
              << "Resources : " << res_dir << "\n";
    if (!g_data_dir.empty())
        std::cerr << "Data dir  : " << g_data_dir << "\n";

    // load initial file
    if (!evio_file.empty())
        loadFileAsync(evio_file);

    // start server
    static WsServer *g_server_ptr = nullptr;
    WsServer server;
    g_server_ptr = &server;

    std::signal(SIGINT, [](int) {
        if (g_server_ptr) {
            try { g_server_ptr->stop_listening(); g_server_ptr->stop(); } catch (...) {}
        }
    });

    server.set_access_channels(websocketpp::log::alevel::none);
    server.set_error_channels(websocketpp::log::elevel::warn | websocketpp::log::elevel::rerror);
    server.init_asio();
    server.set_reuse_addr(true);
    server.set_http_handler([&server](websocketpp::connection_hdl hdl) { onHttp(&server, hdl); });
    server.listen(port);
    server.start_accept();

    {
        auto data = g_data;
        std::cout << "Viewer at http://localhost:" << port << "\n"
                  << "  " << (data ? data->index.size() : 0) << " events"
                  << (g_hist_enabled ? ", histograms enabled" : "")
                  << (g_data_dir.empty() ? "" : ", file browser enabled") << "\n"
                  << "  Ctrl+C to stop\n";
    }

    server.run();

    { std::lock_guard<std::mutex> lk(g_load_mtx); if (g_load_thread.joinable()) g_load_thread.join(); }
    return 0;
}
