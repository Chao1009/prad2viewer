// src/evc_monitor.cpp — HyCal online event monitor
//
// Connects to ET system, reads events in a background thread,
// accumulates histograms, stores latest N events in a ring buffer,
// pushes WebSocket notifications to the viewer on new events.
//
// Usage: evc_monitor [port] [--config online_config.json]

#include "EtChannel.h"
#include "Fadc250Data.h"
#include "Fadc250Decoder.h"
#include "WaveAnalyzer.h"

#include <nlohmann/json.hpp>

#define ASIO_STANDALONE
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <mutex>
#include <thread>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <csignal>

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
// Configuration
// -------------------------------------------------------------------------
struct EtConfig {
    std::string host    = "localhost";
    int         port    = 11111;
    std::string et_file = "/tmp/et_sys_prad2";
    std::string station = "prad2_monitor";
};

struct HistConfig {
    float time_min  = 170;
    float time_max  = 190;
    float bin_min   = 0;
    float bin_max   = 20000;
    float bin_step  = 100;
    float threshold = 3.0;
    float pos_min   = 0;
    float pos_max   = 400;
    float pos_step  = 4;
};

struct Histogram {
    int underflow = 0, overflow = 0;
    std::vector<int> bins;
    void init(int n) { bins.assign(n, 0); underflow = overflow = 0; }
    void fill(float v, float bmin, float bstep) {
        if (v < bmin) { ++underflow; return; }
        int b = (int)((v - bmin) / bstep);
        if (b >= (int)bins.size()) { ++overflow; return; }
        ++bins[b];
    }
    void clear() { std::fill(bins.begin(), bins.end(), 0); underflow = overflow = 0; }
};

// -------------------------------------------------------------------------
// Globals
// -------------------------------------------------------------------------
static EtConfig  g_et_cfg;
static HistConfig g_hist_cfg;
static int g_ring_size = 20;
static int g_hist_nbins = 0, g_pos_nbins = 0;

// ring buffer: decoded event JSON strings, newest at back
struct RingEntry {
    int         seq;
    std::string json_str;
};
static std::deque<RingEntry> g_ring;
static std::mutex g_ring_mtx;

// histograms
static std::map<std::string, Histogram> g_histograms;
static std::map<std::string, Histogram> g_pos_histograms;
static std::atomic<int> g_events_processed{0};
static std::mutex g_hist_mtx;

// WebSocket connections
static std::set<websocketpp::connection_hdl,
                std::owner_less<websocketpp::connection_hdl>> g_ws_clients;
static std::mutex g_ws_mtx;
static WsServer *g_server_ptr = nullptr;

// state
static std::atomic<bool> g_running{true};
static std::atomic<bool> g_et_connected{false};
static std::string g_viewer_html;
static json g_config;

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------
static std::string readFile(const std::string &path) {
    std::ifstream f(path);
    if (!f) return "";
    return {std::istreambuf_iterator<char>(f), {}};
}
static std::string findFile(const std::string &name, const std::string &base) {
    { std::ifstream f(name); if (f.good()) return name; }
    std::string p = base + "/" + name;
    { std::ifstream f(p); if (f.good()) return p; }
    return "";
}

// -------------------------------------------------------------------------
// Encode one decoded event as JSON string
// -------------------------------------------------------------------------
static std::string encodeEvent(fdec::EventData &event, int seq,
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
                if (!(slot.channel_mask & (1u << c))) continue;
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
    return json({{"event", seq}, {"channels", channels}}).dump();
}

// -------------------------------------------------------------------------
// Fill histograms for one event
// -------------------------------------------------------------------------
static void fillHist(fdec::EventData &event, fdec::WaveAnalyzer &ana,
                     fdec::WaveResult &wres)
{
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1u << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.Analyze(cd.samples, cd.nsamples, wres);

                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                // integral histogram: largest peak within time cut
                float best = -1;
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    if (pk.height < g_hist_cfg.threshold) continue;
                    if (pk.time >= g_hist_cfg.time_min && pk.time <= g_hist_cfg.time_max)
                        if (pk.integral > best) best = pk.integral;
                }
                if (best >= 0) {
                    auto &h = g_histograms[key];
                    if (h.bins.empty()) h.init(g_hist_nbins);
                    h.fill(best, g_hist_cfg.bin_min, g_hist_cfg.bin_step);
                }

                // position histogram: all peaks above threshold, no time cut
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    if (pk.height < g_hist_cfg.threshold) continue;
                    auto &ph = g_pos_histograms[key];
                    if (ph.bins.empty()) ph.init(g_pos_nbins);
                    ph.fill(pk.time, g_hist_cfg.pos_min, g_hist_cfg.pos_step);
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
// Notify all connected WebSocket clients
// -------------------------------------------------------------------------
static void wsBroadcast(const std::string &msg)
{
    std::lock_guard<std::mutex> lk(g_ws_mtx);
    for (auto &hdl : g_ws_clients) {
        try { g_server_ptr->send(hdl, msg, websocketpp::frame::opcode::text); }
        catch (...) {}
    }
}

// -------------------------------------------------------------------------
// ET reader thread
// -------------------------------------------------------------------------
static void etReaderThread()
{
    EtChannel ch;
    fdec::EventData event;
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;

    while (g_running) {
        // connect
        std::cerr << "ET: connecting to " << g_et_cfg.host << ":" << g_et_cfg.port
                  << "  " << g_et_cfg.et_file << " ...\n";

        if (ch.Connect(g_et_cfg.host, g_et_cfg.port, g_et_cfg.et_file) != status::success) {
            std::cerr << "ET: connect failed, retrying in 3s\n";
            std::this_thread::sleep_for(std::chrono::seconds(3));
            continue;
        }

        if (ch.Open(g_et_cfg.station) != status::success) {
            std::cerr << "ET: station open failed, retrying in 3s\n";
            ch.Disconnect();
            std::this_thread::sleep_for(std::chrono::seconds(3));
            continue;
        }

        g_et_connected = true;
        wsBroadcast("{\"type\":\"status\",\"connected\":true}");
        std::cerr << "ET: connected, reading events\n";

        // read loop
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

            for (int i = 0; i < ch.GetNEvents(); ++i) {
                if (!ch.DecodeEvent(i, event)) continue;

                int seq = ++g_events_processed;

                // encode
                std::string evjson = encodeEvent(event, seq, ana, wres);

                // fill histograms
                {
                    std::lock_guard<std::mutex> lk(g_hist_mtx);
                    fillHist(event, ana, wres);
                }

                // push to ring buffer
                {
                    std::lock_guard<std::mutex> lk(g_ring_mtx);
                    g_ring.push_back({seq, std::move(evjson)});
                    while ((int)g_ring.size() > g_ring_size)
                        g_ring.pop_front();
                }

                // notify viewers
                wsBroadcast("{\"type\":\"new_event\",\"seq\":" + std::to_string(seq) + "}");
            }
        }

        g_et_connected = false;
        ch.Close();
        ch.Disconnect();
        wsBroadcast("{\"type\":\"status\",\"connected\":false}");

        if (g_running)
            std::this_thread::sleep_for(std::chrono::seconds(3));
    }
}

// -------------------------------------------------------------------------
// Histogram JSON helper
// -------------------------------------------------------------------------
static json getHist(const std::map<std::string, Histogram> &hmap, const std::string &key)
{
    std::lock_guard<std::mutex> lk(g_hist_mtx);
    auto it = hmap.find(key);
    if (it == hmap.end())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"events", g_events_processed.load()}};
    auto &h = it->second;
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"events", g_events_processed.load()}};
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

    if (uri == "/") {
        reply(g_viewer_html, "text/html; charset=utf-8"); return;
    }

    if (uri == "/api/config") {
        // update live fields before sending
        g_config["total_events"] = g_events_processed.load();
        g_config["et_connected"] = g_et_connected.load();
        reply(g_config.dump()); return;
    }

    // /api/ring — list of available event seq numbers
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

    // /api/event/<seq_number>  — fetch by sequence number from ring
    if (uri.rfind("/api/event/", 0) == 0) {
        int seq = std::atoi(uri.c_str() + 11);
        std::lock_guard<std::mutex> lk(g_ring_mtx);
        for (auto &e : g_ring) {
            if (e.seq == seq) { reply(e.json_str); return; }
        }
        reply("{\"error\":\"event not in ring buffer\"}"); return;
    }

    // /api/hist/clear — clear all histograms
    if (uri == "/api/hist/clear") {
        {
            std::lock_guard<std::mutex> lk(g_hist_mtx);
            for (auto &[k, h] : g_histograms)     h.clear();
            for (auto &[k, h] : g_pos_histograms)  h.clear();
            g_events_processed = 0;
        }
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"hist_cleared\"}");
        return;
    }

    // /api/hist/<key>
    if (uri.rfind("/api/hist/", 0) == 0) {
        reply(getHist(g_histograms, uri.substr(10)).dump()); return;
    }

    // /api/poshist/<key>
    if (uri.rfind("/api/poshist/", 0) == 0) {
        reply(getHist(g_pos_histograms, uri.substr(13)).dump()); return;
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

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) config_file = argv[++i];
        else if (arg[0] != '-') port = std::atoi(argv[i]);
    }

    std::string db_dir  = DATABASE_DIR;
    std::string res_dir = RESOURCE_DIR;

    // load config
    if (config_file.empty())
        config_file = findFile("online_config.json", db_dir);

    std::string cfg_str = readFile(config_file);
    if (!cfg_str.empty()) {
        auto cfg = json::parse(cfg_str, nullptr, false);
        if (cfg.contains("et")) {
            auto &e = cfg["et"];
            if (e.contains("host"))    g_et_cfg.host    = e["host"];
            if (e.contains("port"))    g_et_cfg.port    = e["port"];
            if (e.contains("et_file")) g_et_cfg.et_file = e["et_file"];
            if (e.contains("station")) g_et_cfg.station = e["station"];
        }
        if (cfg.contains("ring_buffer_size"))
            g_ring_size = cfg["ring_buffer_size"];
        if (cfg.contains("hist")) {
            auto &h = cfg["hist"];
            if (h.contains("time_min"))  g_hist_cfg.time_min  = h["time_min"];
            if (h.contains("time_max"))  g_hist_cfg.time_max  = h["time_max"];
            if (h.contains("bin_min"))   g_hist_cfg.bin_min   = h["bin_min"];
            if (h.contains("bin_max"))   g_hist_cfg.bin_max   = h["bin_max"];
            if (h.contains("bin_step"))  g_hist_cfg.bin_step  = h["bin_step"];
            if (h.contains("threshold")) g_hist_cfg.threshold = h["threshold"];
            if (h.contains("pos_min"))   g_hist_cfg.pos_min   = h["pos_min"];
            if (h.contains("pos_max"))   g_hist_cfg.pos_max   = h["pos_max"];
            if (h.contains("pos_step"))  g_hist_cfg.pos_step  = h["pos_step"];
        }
        std::cerr << "Config: " << config_file << "\n";
    } else {
        std::cerr << "Config: using defaults\n";
    }

    g_hist_nbins = std::max(1, (int)std::ceil(
        (g_hist_cfg.bin_max - g_hist_cfg.bin_min) / g_hist_cfg.bin_step));
    g_pos_nbins = std::max(1, (int)std::ceil(
        (g_hist_cfg.pos_max - g_hist_cfg.pos_min) / g_hist_cfg.pos_step));

    // load viewer and database
    std::string html_file = findFile("viewer.html", res_dir);
    std::string mod_file  = findFile("hycal_modules.json", db_dir);
    std::string daq_file  = findFile("daq_map.json", db_dir);

    g_viewer_html = readFile(html_file);
    if (g_viewer_html.empty())
        std::cerr << "Warning: viewer.html not found\n";

    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(mod_file); if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    { std::string s = readFile(daq_file); if (!s.empty()) daq_j     = json::parse(s, nullptr, false); }

    g_config = {
        {"modules", modules_j},
        {"daq", daq_j},
        {"crate_roc", {{"0",0x80},{"1",0x82},{"2",0x84},{"3",0x86},{"4",0x88},{"5",0x8a},{"6",0x8c}}},
        {"total_events", 0},
        {"mode", "online"},
        {"hist_enabled", true},
        {"et_connected", false},
        {"ring_buffer_size", g_ring_size},
        {"hist", {
            {"time_min", g_hist_cfg.time_min}, {"time_max", g_hist_cfg.time_max},
            {"bin_min", g_hist_cfg.bin_min}, {"bin_max", g_hist_cfg.bin_max},
            {"bin_step", g_hist_cfg.bin_step}, {"threshold", g_hist_cfg.threshold},
            {"pos_min", g_hist_cfg.pos_min}, {"pos_max", g_hist_cfg.pos_max},
            {"pos_step", g_hist_cfg.pos_step},
        }},
    };

    std::cerr << "Database  : " << db_dir << " ("
              << modules_j.size() << " modules, " << daq_j.size() << " DAQ channels)\n"
              << "ET target : " << g_et_cfg.host << ":" << g_et_cfg.port << "\n"
              << "Ring buf  : " << g_ring_size << " events\n"
              << "Hist bins : " << g_hist_nbins << " integral, " << g_pos_nbins << " position\n";

    // signal handler
    std::signal(SIGINT, [](int) { g_running = false; });

    // start WebSocket/HTTP server
    WsServer server;
    g_server_ptr = &server;

    server.set_access_channels(websocketpp::log::alevel::none);
    server.set_error_channels(websocketpp::log::elevel::warn | websocketpp::log::elevel::rerror);
    server.init_asio();
    server.set_reuse_addr(true);

    server.set_http_handler([&server](websocketpp::connection_hdl hdl) {
        onHttp(&server, hdl);
    });
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

    // start ET reader thread
    std::thread reader(etReaderThread);

    std::cout << "Monitor at http://localhost:" << port << "\n"
              << "  ET: " << g_et_cfg.host << ":" << g_et_cfg.port << "\n"
              << "  Ctrl+C to stop\n";

    server.run();

    // cleanup
    g_running = false;
    if (reader.joinable()) reader.join();
    return 0;
}
