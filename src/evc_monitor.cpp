// src/evc_monitor.cpp — HyCal online event monitor
//
// Connects to ET system, reads events in a background thread,
// accumulates histograms, stores latest N events in a ring buffer,
// pushes WebSocket notifications to the viewer on new events.
//
// Usage: evc_monitor [-p port] [-c online_config.json] [-H hist_config.json]

#include "EtChannel.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"
#include "HyCalSystem.h"
#include "load_daq_config.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <unordered_map>
#include <set>
#include <deque>
#include <mutex>
#include <thread>
#include <atomic>
#include <cmath>
#include <cstdlib>
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
    float min_peak_ratio = 0.3f;
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
static std::map<std::string, int> g_occupancy;
static std::map<std::string, int> g_occupancy_tcut;
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
static std::string g_res_dir;
static json g_config;

// HyCal system for module lookup (LMS + clustering)
static fdec::HyCalSystem g_hycal;
static std::unordered_map<int, int> g_roc_to_crate;
static evc::DaqConfig g_daq_cfg;

// LMS monitoring
struct LmsEntry {
    double time_sec;
    float  integral;
};
static std::map<int, std::vector<LmsEntry>> g_lms_history;  // module_index → entries
static std::atomic<int> g_lms_events{0};
static uint64_t g_lms_first_ts = 0;
static std::mutex g_lms_mtx;

// LMS config
static int      g_lms_trigger_bit = 16;
static float    g_lms_warn_thresh = 0.1f;
static int      g_lms_max_history = 5000;
static uint32_t g_lms_trigger_mask = 0;

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

static std::string contentType(const std::string &path) {
    if (path.size() >= 5 && path.substr(path.size()-5) == ".html") return "text/html; charset=utf-8";
    if (path.size() >= 4 && path.substr(path.size()-4) == ".css")  return "text/css; charset=utf-8";
    if (path.size() >= 3 && path.substr(path.size()-3) == ".js")   return "application/javascript; charset=utf-8";
    return "application/octet-stream";
}

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
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.Analyze(cd.samples, cd.nsamples, wres);

                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                bool has_peak = false, has_peak_tcut = false;

                // integral histogram: all peaks within time cut
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    if (pk.height < g_hist_cfg.threshold) continue;
                    has_peak = true;
                    if (pk.time >= g_hist_cfg.time_min && pk.time <= g_hist_cfg.time_max) {
                        has_peak_tcut = true;
                        auto &h = g_histograms[key];
                        if (h.bins.empty()) h.init(g_hist_nbins);
                        h.fill(pk.integral, g_hist_cfg.bin_min, g_hist_cfg.bin_step);
                    }
                }

                // position histogram
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    if (pk.height < g_hist_cfg.threshold) continue;
                    auto &ph = g_pos_histograms[key];
                    if (ph.bins.empty()) ph.init(g_pos_nbins);
                    ph.fill(pk.time, g_hist_cfg.pos_min, g_hist_cfg.pos_step);
                }

                // occupancy
                if (has_peak)      g_occupancy[key]++;
                if (has_peak_tcut) g_occupancy_tcut[key]++;
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

// Interruptible sleep: returns early if g_running becomes false
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
    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;

    int retry_ms = 3000;       // start at 3s
    const int max_retry = 30000; // cap at 30s
    int retry_count = 0;
    auto retry_start = std::chrono::steady_clock::now();

    while (g_running) {
        // connect
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

        if (retry_count > 0) std::cerr << "\n";  // newline after in-place updates

        if (ch.Open(g_et_cfg.station) != status::success) {
            std::cerr << "ET: station open failed, retrying...\n";
            ch.Disconnect();
            sleepMs(retry_ms);
            retry_ms = std::min(retry_ms * 2, max_retry);
            continue;
        }

        // connected — reset backoff
        retry_ms = 3000;
        retry_count = 0;
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

                // LMS monitoring
                if (g_lms_trigger_mask != 0 &&
                    (event.info.trigger_bits & g_lms_trigger_mask))
                {
                    std::lock_guard<std::mutex> lk(g_lms_mtx);
                    if (g_lms_first_ts == 0) g_lms_first_ts = event.info.timestamp;
                    double time_sec = static_cast<double>(event.info.timestamp - g_lms_first_ts) * 4e-9;

                    bool is_adc1881m = (g_daq_cfg.adc_format == "adc1881m");

                    for (int r = 0; r < event.nrocs; ++r) {
                        auto &roc = event.rocs[r];
                        if (!roc.present) continue;
                        auto cit = g_roc_to_crate.find(roc.tag);
                        if (cit == g_roc_to_crate.end()) continue;
                        int crate = cit->second;

                        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                            if (!roc.slots[s].present) continue;
                            auto &slot = roc.slots[s];
                            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                                if (!(slot.channel_mask & (1ull << c))) continue;
                                auto &cd = slot.channels[c];
                                if (cd.nsamples <= 0) continue;

                                const auto *mod = g_hycal.module_by_daq(crate, s, c);
                                if (!mod || !mod->is_hycal()) continue;

                                float val = 0;
                                if (is_adc1881m) {
                                    val = cd.samples[0];
                                } else {
                                    // use wres from encodeEvent above
                                    ana.Analyze(cd.samples, cd.nsamples, wres);
                                    float best = -1;
                                    for (int p = 0; p < wres.npeaks; ++p) {
                                        auto &pk = wres.peaks[p];
                                        if (pk.height < g_hist_cfg.threshold) continue;
                                        if (pk.time >= g_hist_cfg.time_min && pk.time <= g_hist_cfg.time_max)
                                            if (pk.integral > best) best = pk.integral;
                                    }
                                    val = best;
                                }
                                if (val <= 0) continue;

                                auto &hist = g_lms_history[mod->index];
                                if (static_cast<int>(hist.size()) < g_lms_max_history)
                                    hist.push_back({time_sec, val});
                            }
                        }
                    }
                    g_lms_events++;
                    wsBroadcast("{\"type\":\"lms_event\",\"count\":" + std::to_string(g_lms_events.load()) + "}");
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

        if (g_running) {
            std::cerr << "ET: disconnected, retrying in " << retry_ms/1000 << "s\n";
            sleepMs(retry_ms);
        }
    }
}

// -------------------------------------------------------------------------
// Histogram JSON helper
// -------------------------------------------------------------------------
static json getHist(const std::map<std::string, Histogram> &hmap, int nbins,
                    const std::string &key)
{
    std::lock_guard<std::mutex> lk(g_hist_mtx);
    auto it = hmap.find(key);
    if (it == hmap.end()) {
        // return zero-filled histogram with correct bin count for proper axis range
        return {{"bins", std::vector<int>(nbins, 0)}, {"underflow", 0}, {"overflow", 0},
                {"events", g_events_processed.load()}};
    }
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

    if (serveResource(uri, con)) return;

    if (uri == "/api/config") {
        // update live fields before sending
        g_config["total_events"] = g_events_processed.load();
        g_config["et_connected"] = g_et_connected.load();
        g_config["lms"] = {
            {"trigger_bit", g_lms_trigger_bit},
            {"warn_threshold", g_lms_warn_thresh},
            {"events", g_lms_events.load()},
        };
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

    // /api/hist/clear — clear all histograms and occupancy
    if (uri == "/api/hist/clear") {
        {
            std::lock_guard<std::mutex> lk(g_hist_mtx);
            for (auto &[k, h] : g_histograms)     h.clear();
            for (auto &[k, h] : g_pos_histograms)  h.clear();
            g_occupancy.clear();
            g_occupancy_tcut.clear();
            g_events_processed = 0;
        }
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"hist_cleared\"}");
        return;
    }

    // /api/occupancy
    if (uri == "/api/occupancy") {
        std::lock_guard<std::mutex> lk(g_hist_mtx);
        json jocc = json::object(), jtcut = json::object();
        for (auto &[k,v] : g_occupancy) jocc[k] = v;
        for (auto &[k,v] : g_occupancy_tcut) jtcut[k] = v;
        reply(json({{"occ", jocc}, {"occ_tcut", jtcut},
                     {"total", g_events_processed.load()}}).dump());
        return;
    }

    // /api/hist/<key>
    if (uri.rfind("/api/hist/", 0) == 0) {
        reply(getHist(g_histograms, g_hist_nbins, uri.substr(10)).dump()); return;
    }

    // /api/poshist/<key>
    if (uri.rfind("/api/poshist/", 0) == 0) {
        reply(getHist(g_pos_histograms, g_pos_nbins, uri.substr(13)).dump()); return;
    }

    // /api/lms/<module_index> or /api/lms/summary or /api/lms/clear
    if (uri.rfind("/api/lms/", 0) == 0) {
        std::string sub = uri.substr(9);

        if (sub == "clear") {
            {
                std::lock_guard<std::mutex> lk(g_lms_mtx);
                g_lms_history.clear();
                g_lms_events = 0;
                g_lms_first_ts = 0;
            }
            reply("{\"cleared\":true}");
            wsBroadcast("{\"type\":\"lms_cleared\"}");
            return;
        }

        if (sub == "summary") {
            std::lock_guard<std::mutex> lk(g_lms_mtx);
            json mods = json::object();
            for (auto &[idx, hist] : g_lms_history) {
                if (hist.empty()) continue;
                double sum = 0, sum2 = 0;
                for (auto &e : hist) { sum += e.integral; sum2 += e.integral * e.integral; }
                double mean = sum / hist.size();
                double var = sum2 / hist.size() - mean * mean;
                double rms = var > 0 ? std::sqrt(var) : 0;
                bool warn = (mean > 0 && rms / mean > g_lms_warn_thresh);
                if (idx >= 0 && idx < g_hycal.module_count()) {
                    auto &mod = g_hycal.module(idx);
                    mods[std::to_string(idx)] = {
                        {"name", mod.name}, {"mean", std::round(mean * 10) / 10},
                        {"rms", std::round(rms * 100) / 100},
                        {"count", (int)hist.size()}, {"warn", warn}};
                }
            }
            reply(json({{"modules", mods}, {"events", g_lms_events.load()},
                         {"trigger_bit", g_lms_trigger_bit}}).dump());
            return;
        }

        // /api/lms/<module_index>
        int mod_idx = std::atoi(sub.c_str());
        std::lock_guard<std::mutex> lk(g_lms_mtx);
        auto it = g_lms_history.find(mod_idx);
        if (it == g_lms_history.end() || it->second.empty()) {
            reply(json({{"time", json::array()}, {"integral", json::array()}, {"events", 0}}).dump());
            return;
        }
        auto &hist = it->second;
        json t_arr = json::array(), v_arr = json::array();
        for (auto &e : hist) {
            t_arr.push_back(std::round(e.time_sec * 100) / 100);
            v_arr.push_back(std::round(e.integral * 10) / 10);
        }
        std::string name = (mod_idx >= 0 && mod_idx < g_hycal.module_count())
            ? g_hycal.module(mod_idx).name : "";
        reply(json({{"name", name}, {"time", t_arr}, {"integral", v_arr},
                     {"events", (int)hist.size()}}).dump());
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
    std::string hist_config_file;

    static struct option long_opts[] = {
        {"port",        required_argument, nullptr, 'p'},
        {"config",      required_argument, nullptr, 'c'},
        {"hist-config", required_argument, nullptr, 'H'},
        {"help",        no_argument,       nullptr, '?'},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "p:c:H:", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'p': port = std::atoi(optarg); break;
        case 'c': config_file = optarg; break;
        case 'H': hist_config_file = optarg; break;
        default:
            std::cerr << "Usage: " << argv[0] << " [-p port] [-c online_config.json] [-H hist_config.json]\n";
            return 1;
        }
    }

    std::string db_dir  = DATABASE_DIR;
    std::string res_dir = RESOURCE_DIR;

    // load histogram config (shared with evc_viewer)
    if (hist_config_file.empty())
        hist_config_file = findFile("hist_config.json", db_dir);

    std::string hcfg_str = readFile(hist_config_file);
    if (!hcfg_str.empty()) {
        auto hcfg = json::parse(hcfg_str, nullptr, false);
        if (hcfg.contains("hist")) {
            auto &h = hcfg["hist"];
            if (h.contains("time_min"))  g_hist_cfg.time_min  = h["time_min"];
            if (h.contains("time_max"))  g_hist_cfg.time_max  = h["time_max"];
            if (h.contains("bin_min"))   g_hist_cfg.bin_min   = h["bin_min"];
            if (h.contains("bin_max"))   g_hist_cfg.bin_max   = h["bin_max"];
            if (h.contains("bin_step"))  g_hist_cfg.bin_step  = h["bin_step"];
            if (h.contains("threshold")) g_hist_cfg.threshold = h["threshold"];
            if (h.contains("pos_min"))   g_hist_cfg.pos_min   = h["pos_min"];
            if (h.contains("pos_max"))   g_hist_cfg.pos_max   = h["pos_max"];
            if (h.contains("pos_step"))  g_hist_cfg.pos_step  = h["pos_step"];
            if (h.contains("min_peak_ratio")) g_hist_cfg.min_peak_ratio = h["min_peak_ratio"];
        }
        std::cerr << "Hist config: " << hist_config_file << "\n";
    }

    // load ET / online config
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
        std::cerr << "ET config: " << config_file << "\n";
    } else {
        std::cerr << "ET config: using defaults\n";
    }

    g_hist_nbins = std::max(1, (int)std::ceil(
        (g_hist_cfg.bin_max - g_hist_cfg.bin_min) / g_hist_cfg.bin_step));
    g_pos_nbins = std::max(1, (int)std::ceil(
        (g_hist_cfg.pos_max - g_hist_cfg.pos_min) / g_hist_cfg.pos_step));

    // load reconstruction config (LMS monitoring)
    std::string reco_file = findFile("reconstruction.json", db_dir);
    std::string reco_str = readFile(reco_file);
    if (!reco_str.empty()) {
        auto rcfg = json::parse(reco_str, nullptr, false);
        if (rcfg.contains("lms_monitor")) {
            auto &lm = rcfg["lms_monitor"];
            if (lm.contains("trigger_bit"))    g_lms_trigger_bit   = lm["trigger_bit"];
            if (lm.contains("warn_threshold")) g_lms_warn_thresh   = lm["warn_threshold"];
            if (lm.contains("max_history"))    g_lms_max_history   = lm["max_history"];
            g_lms_trigger_mask = (1u << g_lms_trigger_bit);
        }
        std::cerr << "Reco config: " << reco_file << "\n";
        if (g_lms_trigger_mask)
            std::cerr << "LMS       : trigger_bit=" << g_lms_trigger_bit
                      << " mask=0x" << std::hex << g_lms_trigger_mask << std::dec
                      << " warn=" << g_lms_warn_thresh << "\n";
    }

    // resources directory (viewer.html, viewer.css, viewer.js)
    g_res_dir = res_dir;
    if (readFile(g_res_dir + "/viewer.html").empty())
        std::cerr << "Warning: viewer.html not found in " << g_res_dir << "\n";

    // load database
    std::string mod_file  = findFile("hycal_modules.json", db_dir);
    std::string daq_file  = findFile("daq_map.json", db_dir);

    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(mod_file); if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    { std::string s = readFile(daq_file); if (!s.empty()) daq_j     = json::parse(s, nullptr, false); }

    // initialize HyCal system for module lookup (LMS monitoring)
    if (!mod_file.empty() && !daq_file.empty()) {
        if (g_hycal.Init(mod_file, daq_file))
            std::cerr << "HyCal     : " << g_hycal.module_count() << " modules\n";
    }

    // build ROC tag → crate map from DAQ config
    json crate_roc_j = {{"0",0x80},{"1",0x82},{"2",0x84},{"3",0x86},{"4",0x88},{"5",0x8a},{"6",0x8c}};
    for (auto &[k, v] : crate_roc_j.items())
        g_roc_to_crate[v.get<int>()] = std::stoi(k);

    g_config = {
        {"modules", modules_j},
        {"daq", daq_j},
        {"crate_roc", crate_roc_j},
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

    // signal handler — stop server event loop
    std::signal(SIGINT, [](int) {
        g_running = false;
        if (g_server_ptr) {
            try {
                g_server_ptr->stop_listening();
                g_server_ptr->stop();
            } catch (...) {}
        }
    });

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
    std::cerr << "\nShutting down...\n";
    g_running = false;
    if (reader.joinable()) reader.join();
    std::cerr << "Done.\n";
    return 0;
}
