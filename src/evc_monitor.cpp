// src/evc_monitor.cpp — HyCal online event monitor
//
// Connects to ET system, reads events in a background thread,
// accumulates histograms + LMS data, stores latest N events in a ring buffer,
// pushes WebSocket notifications to the viewer on new events.
//
// Usage: evc_monitor [-p port] [-c config.json] [-D daq_config.json]

#include "EtChannel.h"
#include "app_state.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <set>
#include <deque>
#include <mutex>
#include <thread>
#include <atomic>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <csignal>
#include <getopt.h>

using json = nlohmann::json;
using WsServer = websocketpp::server<websocketpp::config::asio>;
using namespace evc;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif
#ifndef RESOURCE_DIR
#define RESOURCE_DIR "."
#endif

// -------------------------------------------------------------------------
// Monitor-specific config
// -------------------------------------------------------------------------
struct EtConfig {
    std::string host    = "localhost";
    int         port    = 11111;
    std::string et_file = "/tmp/et_sys_prad2";
    std::string station = "prad2_monitor";
};

// -------------------------------------------------------------------------
// Globals
// -------------------------------------------------------------------------
static AppState g_app;          // shared state: config, histograms, LMS, HyCal
static EtConfig g_et_cfg;
static int      g_ring_size = 20;

// ring buffer: decoded event JSON strings, newest at back
struct RingEntry {
    int seq;
    std::string json_str;      // encoded event (channels + waveforms)
    std::string cluster_str;   // clustering result JSON
};
static std::deque<RingEntry> g_ring;
static std::mutex g_ring_mtx;

// WebSocket connections
static std::set<websocketpp::connection_hdl,
                std::owner_less<websocketpp::connection_hdl>> g_ws_clients;
static std::mutex g_ws_mtx;
static WsServer *g_server_ptr = nullptr;

// state
static std::atomic<bool> g_running{true};
static std::atomic<bool> g_et_connected{false};
static std::string g_res_dir;

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------
static bool serveResource(const std::string &uri, WsServer::connection_ptr con)
{
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

static void wsBroadcast(const std::string &msg)
{
    std::lock_guard<std::mutex> lk(g_ws_mtx);
    for (auto &hdl : g_ws_clients) {
        try { g_server_ptr->send(hdl, msg, websocketpp::frame::opcode::text); }
        catch (...) {}
    }
}

static void sleepMs(int ms) {
    for (int elapsed = 0; elapsed < ms && g_running; elapsed += 100)
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
}


// -------------------------------------------------------------------------
// ET reader thread
// -------------------------------------------------------------------------
static void etReaderThread()
{
    EtChannel ch;
    ch.SetConfig(g_app.daq_cfg);
    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_app.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    uint64_t last_ti_ts = 0;

    int retry_ms = 3000;
    const int max_retry = 30000;
    int retry_count = 0;
    auto retry_start = std::chrono::steady_clock::now();

    while (g_running) {
        if (retry_count == 0) {
            std::cerr << "ET: connecting to " << g_et_cfg.host << ":" << g_et_cfg.port
                      << "  " << g_et_cfg.et_file << " ...\n";
            retry_start = std::chrono::steady_clock::now();
        }

        if (ch.Connect(g_et_cfg.host, g_et_cfg.port, g_et_cfg.et_file) != status::success) {
            retry_count++;
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - retry_start).count();
            std::cerr << "\rET: waiting for ET system... "
                      << retry_count << " attempts, " << elapsed << "s elapsed   " << std::flush;
            wsBroadcast("{\"type\":\"status\",\"connected\":false,\"waiting\":true,\"retries\":"
                        + std::to_string(retry_count) + "}");
            sleepMs(retry_ms);
            retry_ms = std::min(retry_ms * 2, max_retry);
            continue;
        }

        if (retry_count > 0) std::cerr << "\n";

        if (ch.Open(g_et_cfg.station) != status::success) {
            std::cerr << "ET: station open failed, retrying...\n";
            ch.Disconnect();
            sleepMs(retry_ms);
            retry_ms = std::min(retry_ms * 2, max_retry);
            continue;
        }

        retry_ms = 3000;
        retry_count = 0;
        g_et_connected = true;
        wsBroadcast("{\"type\":\"status\",\"connected\":true}");
        std::cerr << "ET: connected, reading events\n";

        while (g_running) {
            auto st = ch.Read();
            if (st == status::empty) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                continue;
            }
            if (st != status::success) {
                std::cerr << "ET: read error, reconnecting\n";
                break;
            }
            if (!ch.Scan()) continue;

            // capture sync time reference
            if (g_app.sync_unix == 0) {
                uint32_t ct = ch.GetControlTime();
                if (ct != 0) g_app.recordSyncTime(ct, last_ti_ts);
            }

            // EPICS slow control event
            if (ch.GetEventType() == evc::EventType::Epics) {
                std::string text = ch.ExtractEpicsText();
                if (!text.empty()) {
                    int seq = g_app.events_processed.load();
                    g_app.processEpics(text, seq, last_ti_ts);
                    wsBroadcast("{\"type\":\"epics_event\",\"count\":" +
                                std::to_string(g_app.epics_events.load()) + "}");
                }
            }

            for (int i = 0; i < ch.GetNEvents(); ++i) {
                if (!ch.DecodeEvent(i, event)) continue;
                last_ti_ts = event.info.timestamp;

                int seq = g_app.events_processed.load() + 1;

                // encode event + clusters for ring buffer
                std::string evjson = g_app.encodeEventJson(event, seq, ana, wres).dump();
                std::string cljson = g_app.computeClustersJson(event, seq, ana, wres).dump();

                // process: histograms + clustering + LMS (thread-safe)
                g_app.processEvent(event, ana, wres);

                // LMS WebSocket notification
                if (g_app.lms_trigger_mask != 0 &&
                    (event.info.trigger_bits & g_app.lms_trigger_mask))
                    wsBroadcast("{\"type\":\"lms_event\",\"count\":" +
                                std::to_string(g_app.lms_events.load()) + "}");

                // push to ring buffer
                {
                    std::lock_guard<std::mutex> lk(g_ring_mtx);
                    g_ring.push_back({seq, std::move(evjson), std::move(cljson)});
                    while ((int)g_ring.size() > g_ring_size) g_ring.pop_front();
                }

                wsBroadcast("{\"type\":\"new_event\",\"seq\":" + std::to_string(seq) + "}");
            }
        }

        g_et_connected = false;
        ch.Close();
        ch.Disconnect();
        wsBroadcast("{\"type\":\"status\",\"connected\":false}");

        if (g_running) {
            std::cerr << "ET: disconnected, retrying in " << retry_ms/1000 << "s\n";
            sleepMs(retry_ms);
        }
    }
}

// -------------------------------------------------------------------------
// HTTP handler
// -------------------------------------------------------------------------
static void onHttp(WsServer *srv, websocketpp::connection_hdl hdl)
{
    auto con = srv->get_con_from_hdl(hdl);
    std::string uri = con->get_resource();

    auto reply = [&](const std::string &body, const std::string &ct = "application/json") {
        con->set_status(websocketpp::http::status_code::ok);
        con->set_body(body);
        con->append_header("Content-Type", ct);
    };

    if (serveResource(uri, con)) return;

    // /api/config
    if (uri == "/api/config") {
        json cfg = g_app.base_config;
        cfg["total_events"] = g_app.events_processed.load();
        cfg["et_connected"] = g_et_connected.load();
        cfg["mode"] = "online";
        cfg["hist_enabled"] = true;
        cfg["ring_buffer_size"] = g_ring_size;
        cfg["hist"] = {
            {"time_min", g_app.hist_cfg.time_min}, {"time_max", g_app.hist_cfg.time_max},
            {"bin_min", g_app.hist_cfg.bin_min}, {"bin_max", g_app.hist_cfg.bin_max},
            {"bin_step", g_app.hist_cfg.bin_step}, {"threshold", g_app.hist_cfg.threshold},
            {"pos_min", g_app.hist_cfg.pos_min}, {"pos_max", g_app.hist_cfg.pos_max},
            {"pos_step", g_app.hist_cfg.pos_step},
        };
        cfg["cluster_hist"] = {
            {"min", g_app.cl_hist_min}, {"max", g_app.cl_hist_max}, {"step", g_app.cl_hist_step},
        };
        cfg["nclusters_hist"] = {
            {"min", g_app.nclusters_hist_min}, {"max", g_app.nclusters_hist_max}, {"step", g_app.nclusters_hist_step},
        };
        cfg["nblocks_hist"] = {
            {"min", g_app.nblocks_hist_min}, {"max", g_app.nblocks_hist_max}, {"step", g_app.nblocks_hist_step},
        };
        cfg["color_ranges"] = g_app.apiColorRanges();
        cfg["refresh_ms"] = {{"event", g_app.refresh_event_ms}, {"ring", g_app.refresh_ring_ms},
                             {"histogram", g_app.refresh_hist_ms}, {"lms", g_app.refresh_lms_ms}};
        cfg["lms"] = {
            {"trigger_bit", g_app.lms_trigger_bit},
            {"warn_threshold", g_app.lms_warn_thresh},
            {"events", g_app.lms_events.load()},
            {"ref_channels", g_app.apiLmsRefChannels()},
        };
        cfg["elog"] = {
            {"url", g_app.elog_url}, {"logbook", g_app.elog_logbook},
            {"author", g_app.elog_author}, {"tags", g_app.elog_tags},
        };
        cfg["epics"] = {
            {"max_history", g_app.epics_max_history},
            {"warn_threshold", g_app.epics_warn_thresh},
            {"alert_threshold", g_app.epics_alert_thresh},
            {"min_avg_points", g_app.epics_min_avg_pts},
            {"mean_window", g_app.epics_mean_window},
            {"slots", g_app.epics_default_slots},
        };
        reply(cfg.dump()); return;
    }

    // /api/ring
    if (uri == "/api/ring") {
        std::lock_guard<std::mutex> lk(g_ring_mtx);
        json arr = json::array();
        for (auto &e : g_ring) arr.push_back(e.seq);
        reply(json({{"ring", arr}, {"latest", g_ring.empty() ? 0 : g_ring.back().seq}}).dump());
        return;
    }

    // /api/event/latest
    if (uri == "/api/event/latest") {
        std::lock_guard<std::mutex> lk(g_ring_mtx);
        if (g_ring.empty()) { reply("{\"error\":\"no events yet\"}"); return; }
        reply(g_ring.back().json_str); return;
    }

    // /api/event/<seq>
    if (uri.rfind("/api/event/", 0) == 0) {
        int seq = std::atoi(uri.c_str() + 11);
        std::lock_guard<std::mutex> lk(g_ring_mtx);
        for (auto &e : g_ring) {
            if (e.seq == seq) { reply(e.json_str); return; }
        }
        reply("{\"error\":\"event not in ring buffer\"}"); return;
    }

    // /api/clusters/<seq> — fetch pre-computed clusters from ring buffer
    if (uri.rfind("/api/clusters/", 0) == 0) {
        int seq = std::atoi(uri.c_str() + 14);
        std::lock_guard<std::mutex> lk(g_ring_mtx);
        for (auto &e : g_ring) {
            if (e.seq == seq) { reply(e.cluster_str); return; }
        }
        reply("{\"error\":\"event not in ring buffer\"}"); return;
    }

    // /api/hist/clear
    if (uri == "/api/hist/clear") {
        g_app.clearHistograms();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"hist_cleared\"}");
        return;
    }

    // /api/occupancy
    if (uri == "/api/occupancy") { reply(g_app.apiOccupancy().dump()); return; }

    // /api/cluster_hist
    if (uri == "/api/cluster_hist") { reply(g_app.apiClusterHist().dump()); return; }

    // /api/hist/<key>
    if (uri.rfind("/api/hist/", 0) == 0) {
        reply(g_app.apiHist(true, uri.substr(10)).dump()); return;
    }

    // /api/poshist/<key>
    if (uri.rfind("/api/poshist/", 0) == 0) {
        reply(g_app.apiHist(false, uri.substr(13)).dump()); return;
    }

    // /api/lms/refs
    if (uri == "/api/lms/refs") { reply(g_app.apiLmsRefChannels().dump()); return; }

    // /api/lms/*
    if (uri.rfind("/api/lms/", 0) == 0) {
        // parse ref= query param
        int ref = -1;
        auto qpos = uri.find('?');
        std::string path_part = (qpos != std::string::npos) ? uri.substr(9, qpos - 9) : uri.substr(9);
        if (qpos != std::string::npos) {
            std::string q = uri.substr(qpos + 1);
            if (q.rfind("ref=", 0) == 0) ref = std::atoi(q.c_str() + 4);
        }
        if (path_part == "clear") {
            g_app.clearLms();
            reply("{\"cleared\":true}");
            wsBroadcast("{\"type\":\"lms_cleared\"}");
            return;
        }
        if (path_part == "summary") { reply(g_app.apiLmsSummary(ref).dump()); return; }
        reply(g_app.apiLmsModule(std::atoi(path_part.c_str()), ref).dump()); return;
    }

    // /api/epics/*
    if (uri.rfind("/api/epics/", 0) == 0) {
        std::string path_part = uri.substr(11);
        if (path_part == "channels") { reply(g_app.apiEpicsChannels().dump()); return; }
        if (path_part == "latest")   { reply(g_app.apiEpicsLatest().dump()); return; }
        if (path_part == "clear") {
            g_app.clearEpics();
            reply("{\"cleared\":true}");
            wsBroadcast("{\"type\":\"epics_cleared\"}");
            return;
        }
        // /api/epics/channel/<name>
        if (path_part.rfind("channel/", 0) == 0) {
            std::string name = path_part.substr(8);
            reply(g_app.apiEpicsChannel(name).dump()); return;
        }
    }

    // /api/elog/post — proxy elog submission via curl
    if (uri == "/api/elog/post") {
        std::string body = con->get_request_body();
        if (body.empty()) {
            con->set_status(websocketpp::http::status_code::bad_request);
            con->set_body("{\"ok\":false,\"error\":\"Empty body\"}");
            con->append_header("Content-Type", "application/json");
            return;
        }
        if (g_app.elog_url.empty()) {
            reply("{\"ok\":false,\"error\":\"No elog URL configured\"}");
            return;
        }
        // parse JSON: {xml}
        auto req = json::parse(body, nullptr, false);
        if (req.is_discarded() || !req.contains("xml")) {
            reply("{\"ok\":false,\"error\":\"Invalid request\"}");
            return;
        }
        std::string xml_body = req["xml"].get<std::string>();
        // write XML to temp file
        std::string tmp = "/tmp/prad2_elog_" + std::to_string(std::time(nullptr)) + ".xml";
        { std::ofstream f(tmp); f << xml_body; }
        // curl it to elog server with SSL client certificate
        std::string cert_flag;
        if (!g_app.elog_cert.empty())
            cert_flag = " --cert '" + g_app.elog_cert + "' --key '" + g_app.elog_key + "'";
        std::string cmd = "curl -s -o /dev/null -w '%{http_code}'" + cert_flag
                        + " --upload-file '" + tmp + "' '"
                        + g_app.elog_url + "/incoming/prad2_report.xml' 2>/dev/null";
        std::string http_code;
        FILE *p = popen(cmd.c_str(), "r");
        if (p) {
            char buf[256] = {};
            if (fgets(buf, sizeof(buf), p)) http_code = buf;
            // trim whitespace/newlines
            while (!http_code.empty() && (http_code.back()=='\n'||http_code.back()=='\r'))
                http_code.pop_back();
            pclose(p);
        }
        std::remove(tmp.c_str());
        bool ok = (http_code.find("200") != std::string::npos ||
                   http_code.find("201") != std::string::npos);
        std::cerr << "Elog post: " << g_app.elog_url
                  << " -> HTTP " << http_code
                  << (ok ? " OK" : " FAIL") << "\n";
        reply(json({{"ok", ok}, {"status", http_code}}).dump());
        return;
    }

    con->set_status(websocketpp::http::status_code::not_found);
    con->set_body("404 Not Found");
}

// -------------------------------------------------------------------------
// Main
// -------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    int port = 5051;
    std::string config_file;
    std::string daq_config_file;

    static struct option long_opts[] = {
        {"port",        required_argument, nullptr, 'p'},
        {"config",      required_argument, nullptr, 'c'},
        {"daq-config",  required_argument, nullptr, 'D'},
        {"help",        no_argument,       nullptr, '?'},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "p:c:D:", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'p': port = std::atoi(optarg); break;
        case 'c': config_file = optarg; break;
        case 'D': daq_config_file = optarg; break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " [-p port] [-c config.json] [-D daq_config.json]\n";
            return 1;
        }
    }

    std::string db_dir  = DATABASE_DIR;
    std::string res_dir = RESOURCE_DIR;

    // load ET config from config.json "online" section, or legacy online_config.json
    if (config_file.empty()) {
        config_file = findFile("config.json", db_dir);
        if (config_file.empty()) config_file = findFile("online_config.json", db_dir);
    }
    std::string cfg_str = readFile(config_file);
    if (!cfg_str.empty()) {
        auto cfg = json::parse(cfg_str, nullptr, false);
        // new format: "online" section in config.json
        if (cfg.contains("online")) {
            auto &e = cfg["online"];
            if (e.contains("et_host"))    g_et_cfg.host    = e["et_host"];
            if (e.contains("et_port"))    g_et_cfg.port    = e["et_port"];
            if (e.contains("et_file"))    g_et_cfg.et_file = e["et_file"];
            if (e.contains("et_station")) g_et_cfg.station = e["et_station"];
            if (e.contains("ring_buffer_size")) g_ring_size = e["ring_buffer_size"];
        }
        // legacy format: "et" section in online_config.json
        else if (cfg.contains("et")) {
            auto &e = cfg["et"];
            if (e.contains("host"))    g_et_cfg.host    = e["host"];
            if (e.contains("port"))    g_et_cfg.port    = e["port"];
            if (e.contains("et_file")) g_et_cfg.et_file = e["et_file"];
            if (e.contains("station")) g_et_cfg.station = e["station"];
            if (cfg.contains("ring_buffer_size")) g_ring_size = cfg["ring_buffer_size"];
        }
        std::cerr << "Config    : " << config_file << "\n";
    } else {
        std::cerr << "Config    : using defaults\n";
    }

    // initialize shared state (DAQ config, HyCal, histograms, clustering, LMS)
    g_app.init(db_dir, daq_config_file, config_file);

    // build base_config JSON for /api/config
    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(findFile("hycal_modules.json", db_dir));
      if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    { // use the same daq file that AppState resolved
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
        {"modules", modules_j},
        {"daq", daq_j},
        {"crate_roc", g_app.crate_roc_json},
    };

    // resources
    g_res_dir = res_dir;
    if (readFile(g_res_dir + "/viewer.html").empty())
        std::cerr << "Warning: viewer.html not found in " << g_res_dir << "\n";

    std::cerr << "ET target : " << g_et_cfg.host << ":" << g_et_cfg.port << "\n"
              << "Ring buf  : " << g_ring_size << " events\n"
              << "Hist bins : " << g_app.hist_nbins << " integral, "
              << g_app.pos_nbins << " position\n";

    // signal handler
    std::signal(SIGINT, [](int) {
        g_running = false;
        if (g_server_ptr) {
            try { g_server_ptr->stop_listening(); g_server_ptr->stop(); } catch (...) {}
        }
    });

    // start server
    WsServer server;
    g_server_ptr = &server;

    server.set_access_channels(websocketpp::log::alevel::none);
    server.set_error_channels(websocketpp::log::elevel::warn | websocketpp::log::elevel::rerror);
    server.init_asio();
    server.set_reuse_addr(true);

    server.set_http_handler([&server](websocketpp::connection_hdl hdl) { onHttp(&server, hdl); });
    server.set_open_handler([](websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(g_ws_mtx);
        g_ws_clients.insert(hdl);
    });
    server.set_close_handler([](websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(g_ws_mtx);
        g_ws_clients.erase(hdl);
    });

    server.listen(port);
    server.start_accept();

    std::thread reader(etReaderThread);

    std::cout << "Monitor at http://localhost:" << port << "\n"
              << "  ET: " << g_et_cfg.host << ":" << g_et_cfg.port << "\n"
              << "  Ctrl+C to stop\n";

    server.run();

    std::cerr << "\nShutting down...\n";
    g_running = false;
    if (reader.joinable()) reader.join();
    std::cerr << "Done.\n";
    return 0;
}
