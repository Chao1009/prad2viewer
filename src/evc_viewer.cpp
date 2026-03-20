// src/evc_viewer.cpp — HyCal event viewer
//
// Usage:
//   evc_viewer [evio_file] [-p port] [-H] [-c config.json] [-d /path/to/data]
//   evc_viewer -d /data/stage6 -H              # browse and pick from GUI
//   evc_viewer data.evio -H                    # open file directly
//   evc_viewer data.evio -H -d /data/stage6    # open file + enable browsing
//   evc_viewer                                 # empty viewer, no file browser
//
// -d enables file browsing: the viewer shows a file picker limited to
// .evio files under that directory tree. Selecting a new file triggers
// background re-indexing + histogram building with progress updates.

#include "EvChannel.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
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
// Forward declarations
// -------------------------------------------------------------------------
struct HistConfig;
struct Histogram;
struct EventIndex { int buffer_num, sub_event; };

// -------------------------------------------------------------------------
// Histogram types
// -------------------------------------------------------------------------
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
};

// -------------------------------------------------------------------------
// Data container: holds everything for one loaded file
// -------------------------------------------------------------------------
struct FileData {
    std::string filepath;
    std::vector<EventIndex> index;
    std::map<std::string, Histogram> histograms;
    std::map<std::string, Histogram> pos_histograms;
    std::map<std::string, int> occupancy;       // events with ≥1 peak above threshold
    std::map<std::string, int> occupancy_tcut;   // same, but peak must be in time window
    int hist_events_processed = 0;
};

// -------------------------------------------------------------------------
// Loading progress
// -------------------------------------------------------------------------
struct Progress {
    std::atomic<bool> loading{false};
    std::atomic<int>  phase{0};         // 0=idle, 1=indexing, 2=histograms
    std::atomic<int>  current{0};       // buffers or events processed so far
    std::atomic<int>  total{0};         // estimated total (0 if unknown)
    std::string       target_file;      // file being loaded
    std::mutex        mtx;              // protects target_file

    json toJson() {
        std::lock_guard<std::mutex> lk(mtx);
        return {
            {"loading", loading.load()},
            {"phase", phase == 1 ? "indexing" : phase == 2 ? "histograms" : "idle"},
            {"current", current.load()},
            {"total", total.load()},
            {"file", target_file},
        };
    }

    void setFile(const std::string &f) {
        std::lock_guard<std::mutex> lk(mtx);
        target_file = f;
    }
};

// -------------------------------------------------------------------------
// Globals
// -------------------------------------------------------------------------
static std::shared_ptr<FileData> g_data;  // current active data (read by HTTP threads)
static std::mutex g_data_mtx;             // protects swap of g_data pointer
static std::thread g_load_thread;         // background loading thread
static std::mutex g_load_mtx;             // protects g_load_thread join

static HistConfig g_hist_cfg;
static bool g_hist_enabled = false;
static int  g_hist_nbins = 0, g_pos_nbins = 0;

static std::string g_data_dir;            // sandboxed data directory (empty = disabled)
static std::string g_res_dir;             // resources directory (viewer.html, .css, .js)
static json g_base_config;                // config without per-file fields
static Progress g_progress;

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

// Serve a file from the resources directory (no directory traversal)
static bool serveResource(const std::string &uri, websocketpp::server<websocketpp::config::asio>::connection_ptr con)
{
    if (g_res_dir.empty()) return false;

    // map "/" to "/viewer.html"
    std::string relpath = (uri == "/") ? "viewer.html" : uri.substr(1);

    // reject anything with .. or absolute paths
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
// File listing (sandboxed)
// -------------------------------------------------------------------------
static json listFiles(const std::string &data_dir)
{
    json files = json::array();
    if (data_dir.empty()) return files;

    try {
        fs::path root(data_dir);
        for (auto &entry : fs::recursive_directory_iterator(root,
                fs::directory_options::skip_permission_denied)) {
            if (!entry.is_regular_file()) continue;
            // match: *.evio, *.evio.0, *.evio.00000, etc.
            std::string fname = entry.path().filename().string();
            if (fname.find(".evio") == std::string::npos) continue;
            {
                // relative path from data_dir
                auto rel = fs::relative(entry.path(), root).string();
                auto sz = entry.file_size();
                files.push_back({
                    {"path", rel},
                    {"size", sz},
                    {"size_mb", std::round(sz / 1048576.0 * 10) / 10},
                });
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "listFiles error: " << e.what() << "\n";
    }

    // sort by path
    std::sort(files.begin(), files.end(),
        [](const json &a, const json &b) { return a["path"] < b["path"]; });
    return files;
}

// -------------------------------------------------------------------------
// Validate a requested file path (prevent directory traversal)
// -------------------------------------------------------------------------
static std::string resolveDataFile(const std::string &relpath)
{
    if (g_data_dir.empty()) return "";
    try {
        fs::path full = fs::canonical(fs::path(g_data_dir) / relpath);
        fs::path root = fs::canonical(fs::path(g_data_dir));
        // must be under root
        auto rootStr = root.string();
        if (full.string().rfind(rootStr, 0) != 0) return "";
        if (!fs::is_regular_file(full)) return "";
        return full.string();
    } catch (...) { return ""; }
}

// -------------------------------------------------------------------------
// Build index for a file
// -------------------------------------------------------------------------
static void buildIndex(const std::string &path, std::vector<EventIndex> &index,
                       Progress &prog)
{
    index.clear();
    EvChannel ch;
    if (ch.Open(path) != status::success) return;

    prog.phase = 1;
    prog.current = 0;
    int buf = 0;
    while (ch.Read() == status::success) {
        ++buf;
        prog.current = buf;
        if (!ch.Scan()) continue;
        for (int i = 0; i < ch.GetNEvents(); ++i)
            index.push_back({buf, i});
    }
    ch.Close();
    prog.total = (int)index.size();
}

// -------------------------------------------------------------------------
// Build histograms
// -------------------------------------------------------------------------
static void buildHistograms(const std::string &path,
                            std::map<std::string, Histogram> &hists,
                            std::map<std::string, Histogram> &pos_hists,
                            std::map<std::string, int> &occ,
                            std::map<std::string, int> &occ_tcut,
                            int &events_out, Progress &prog)
{
    hists.clear();
    pos_hists.clear();
    occ.clear();
    occ_tcut.clear();
    events_out = 0;

    EvChannel ch;
    if (ch.Open(path) != status::success) return;

    fdec::EventData event;
    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;

    prog.phase = 2;
    prog.current = 0;

    int buf = 0, total = 0;
    while (ch.Read() == status::success) {
        ++buf;
        if (!ch.Scan()) continue;
        for (int i = 0; i < ch.GetNEvents(); ++i) {
            if (!ch.DecodeEvent(i, event)) continue;
            ++total;
            prog.current = total;

            // fill histograms
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

                        bool has_peak = false, has_peak_tcut = false;

                        // integral histogram + occupancy
                        float best = -1;
                        for (int p = 0; p < wres.npeaks; ++p) {
                            auto &pk = wres.peaks[p];
                            if (pk.height < g_hist_cfg.threshold) continue;
                            has_peak = true;
                            if (pk.time >= g_hist_cfg.time_min && pk.time <= g_hist_cfg.time_max) {
                                has_peak_tcut = true;
                                if (pk.integral > best) best = pk.integral;
                            }
                        }
                        if (best >= 0) {
                            auto &h = hists[key];
                            if (h.bins.empty()) h.init(g_hist_nbins);
                            h.fill(best, g_hist_cfg.bin_min, g_hist_cfg.bin_step);
                        }

                        // position histogram
                        for (int p = 0; p < wres.npeaks; ++p) {
                            auto &pk = wres.peaks[p];
                            if (pk.height < g_hist_cfg.threshold) continue;
                            auto &ph = pos_hists[key];
                            if (ph.bins.empty()) ph.init(g_pos_nbins);
                            ph.fill(pk.time, g_hist_cfg.pos_min, g_hist_cfg.pos_step);
                        }

                        // occupancy
                        if (has_peak)      occ[key]++;
                        if (has_peak_tcut) occ_tcut[key]++;
                    }
                }
            }
        }
    }
    ch.Close();
    events_out = total;
}

// -------------------------------------------------------------------------
// Load a file: index + optional histograms, then swap into g_data
// -------------------------------------------------------------------------
static void loadFileAsync(const std::string &filepath)
{
    g_progress.loading = true;
    g_progress.setFile(filepath);
    g_progress.phase = 0;
    g_progress.current = 0;
    g_progress.total = 0;

    auto data = std::make_shared<FileData>();
    data->filepath = filepath;

    std::cerr << "Loading: " << filepath << "\n";

    // index
    buildIndex(filepath, data->index, g_progress);
    std::cerr << "  Indexed " << data->index.size() << " events\n";

    // histograms
    if (g_hist_enabled) {
        g_progress.total = (int)data->index.size();
        buildHistograms(filepath, data->histograms, data->pos_histograms,
                        data->occupancy, data->occupancy_tcut,
                        data->hist_events_processed, g_progress);
        std::cerr << "  Histograms: " << data->hist_events_processed << " events, "
                  << data->histograms.size() << " channels\n";
    }

    // atomic swap
    {
        std::lock_guard<std::mutex> lk(g_data_mtx);
        g_data = data;
    }

    g_progress.loading = false;
    g_progress.phase = 0;
    std::cerr << "  Ready\n";
}

// -------------------------------------------------------------------------
// Decode one event → JSON (from current g_data)
// -------------------------------------------------------------------------
static json decodeEvent(int ev1)
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(g_data_mtx); data = g_data; }
    if (!data) return {{"error", "no file loaded"}};

    int idx = ev1 - 1;
    if (idx < 0 || idx >= (int)data->index.size())
        return {{"error", "event out of range"}};

    auto &ei = data->index[idx];
    EvChannel ch;
    if (ch.Open(data->filepath) != status::success)
        return {{"error", "cannot open file"}};

    for (int b = 0; b < ei.buffer_num; ++b)
        if (ch.Read() != status::success) { ch.Close(); return {{"error", "read error"}}; }
    if (!ch.Scan()) { ch.Close(); return {{"error", "scan error"}}; }

    fdec::EventData event;
    if (!ch.DecodeEvent(ei.sub_event, event)) { ch.Close(); return {{"error", "decode error"}}; }
    ch.Close();

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = g_hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
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
    return {{"event", ev1}, {"channels", channels}};
}

// -------------------------------------------------------------------------
// Histogram API
// -------------------------------------------------------------------------
static json getHist(bool integral, const std::string &key)
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(g_data_mtx); data = g_data; }
    if (!data || !g_hist_enabled)
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0}, {"events", 0}};

    auto &hmap = integral ? data->histograms : data->pos_histograms;
    auto it = hmap.find(key);
    if (it == hmap.end())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"events", data->hist_events_processed}};
    auto &h = it->second;
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"events", data->hist_events_processed}};
}

// -------------------------------------------------------------------------
// Build config JSON for current state
// -------------------------------------------------------------------------
static json buildConfig()
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(g_data_mtx); data = g_data; }

    json cfg = g_base_config;
    cfg["total_events"] = data ? (int)data->index.size() : 0;
    cfg["current_file"] = data ? data->filepath : "";
    cfg["data_dir_enabled"] = !g_data_dir.empty();
    cfg["data_dir"] = g_data_dir;
    cfg["hist_enabled"] = g_hist_enabled;
    // always include hist config so the client can use ranges for coloring
    cfg["hist"] = {
        {"time_min", g_hist_cfg.time_min}, {"time_max", g_hist_cfg.time_max},
        {"bin_min", g_hist_cfg.bin_min}, {"bin_max", g_hist_cfg.bin_max},
        {"bin_step", g_hist_cfg.bin_step}, {"threshold", g_hist_cfg.threshold},
        {"pos_min", g_hist_cfg.pos_min}, {"pos_max", g_hist_cfg.pos_max},
        {"pos_step", g_hist_cfg.pos_step},
    };
    return cfg;
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

    if (uri == "/api/config") { reply(buildConfig().dump()); return; }

    // /api/event/<num>
    if (uri.rfind("/api/event/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 11);
        reply(decodeEvent(evnum).dump()); return;
    }

    // /api/hist/<key>
    if (uri.rfind("/api/hist/", 0) == 0) {
        reply(getHist(true, uri.substr(10)).dump()); return;
    }

    // /api/poshist/<key>
    if (uri.rfind("/api/poshist/", 0) == 0) {
        reply(getHist(false, uri.substr(13)).dump()); return;
    }

    // /api/occupancy — per-channel event counts (for geo view)
    if (uri == "/api/occupancy") {
        std::shared_ptr<FileData> data;
        { std::lock_guard<std::mutex> lk(g_data_mtx); data = g_data; }
        if (!data) { reply("{\"occ\":{},\"occ_tcut\":{},\"total\":0}"); return; }
        json jocc = json::object(), jtcut = json::object();
        for (auto &[k,v] : data->occupancy) jocc[k] = v;
        for (auto &[k,v] : data->occupancy_tcut) jtcut[k] = v;
        reply(json({{"occ", jocc}, {"occ_tcut", jtcut},
                     {"total", data->hist_events_processed}}).dump());
        return;
    }

    // /api/files — list .evio files under data-dir
    if (uri == "/api/files") {
        reply(json({{"files", listFiles(g_data_dir)}}).dump()); return;
    }

    // /api/progress — loading progress
    if (uri == "/api/progress") {
        reply(g_progress.toJson().dump()); return;
    }

    // /api/load?file=relative/path.evio&hist=1 — switch to a new file
    if (uri.rfind("/api/load?", 0) == 0) {
        if (g_data_dir.empty()) {
            reply("{\"error\":\"file browsing not enabled (use --data-dir)\"}"); return;
        }
        bool expected = false;
        if (!g_progress.loading.compare_exchange_strong(expected, true)) {
            reply("{\"error\":\"already loading a file\"}"); return;
        }

        // parse query string
        std::string qs = uri.substr(10);  // after "/api/load?"
        auto urlDecode = [](const std::string &raw) {
            std::string out;
            for (size_t i = 0; i < raw.size(); ++i) {
                if (raw[i] == '+') out += ' ';
                else if (raw[i] == '%' && i + 2 < raw.size()) {
                    out += (char)std::stoi(raw.substr(i+1, 2), nullptr, 16);
                    i += 2;
                } else out += raw[i];
            }
            return out;
        };
        auto getParam = [&](const std::string &key) -> std::string {
            std::string prefix = key + "=";
            size_t pos = qs.find(prefix);
            if (pos == std::string::npos) return "";
            size_t start = pos + prefix.size();
            size_t end = qs.find('&', start);
            return urlDecode(qs.substr(start, end == std::string::npos ? end : end - start));
        };

        std::string relpath = getParam("file");
        std::string hist_param = getParam("hist");

        if (relpath.empty()) {
            g_progress.loading = false;
            reply("{\"error\":\"missing file parameter\"}"); return;
        }

        std::string fullpath = resolveDataFile(relpath);
        if (fullpath.empty()) {
            g_progress.loading = false;
            reply("{\"error\":\"file not found or access denied\"}"); return;
        }

        // update histogram enabled flag
        if (!hist_param.empty())
            g_hist_enabled = (hist_param == "1" || hist_param == "true");

        // start background loading (join any previous load first)
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
int main(int argc, char *argv[])
{
    std::string evio_file;
    int port = 5050;
    std::string hist_config_file;

    static struct option long_opts[] = {
        {"port",        required_argument, nullptr, 'p'},
        {"hist",        no_argument,       nullptr, 'H'},
        {"hist-config", required_argument, nullptr, 'c'},
        {"data-dir",    required_argument, nullptr, 'd'},
        {"help",        no_argument,       nullptr, '?'},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "p:Hc:d:", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'p': port = std::atoi(optarg); break;
        case 'H': g_hist_enabled = true; break;
        case 'c': hist_config_file = optarg; g_hist_enabled = true; break;
        case 'd': g_data_dir = optarg; break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " [evio_file] [-p port] [-H] [-c hist_config.json] [-d data_dir]\n";
            return 1;
        }
    }
    if (optind < argc)
        evio_file = argv[optind];

    std::string db_dir  = DATABASE_DIR;
    std::string res_dir = RESOURCE_DIR;

    // always load histogram config (needed if user enables hist via GUI later)
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

    g_hist_nbins = std::max(1, (int)std::ceil(
        (g_hist_cfg.bin_max - g_hist_cfg.bin_min) / g_hist_cfg.bin_step));
    g_pos_nbins = std::max(1, (int)std::ceil(
        (g_hist_cfg.pos_max - g_hist_cfg.pos_min) / g_hist_cfg.pos_step));

    // resources directory (viewer.html, viewer.css, viewer.js)
    g_res_dir = res_dir;
    if (readFile(g_res_dir + "/viewer.html").empty())
        std::cerr << "Warning: viewer.html not found in " << g_res_dir << "\n";

    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(findFile("hycal_modules.json", db_dir));
      if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    { std::string s = readFile(findFile("daq_map.json", db_dir));
      if (!s.empty()) daq_j = json::parse(s, nullptr, false); }

    // base config (file-independent fields)
    g_base_config = {
        {"modules", modules_j},
        {"daq", daq_j},
        {"crate_roc", {{"0",0x80},{"1",0x82},{"2",0x84},{"3",0x86},{"4",0x88},{"5",0x8a},{"6",0x8c}}},
        {"mode", "file"},
    };

    std::cerr << "Database  : " << db_dir << "\n"
              << "Resources : " << res_dir << "\n";
    if (!g_data_dir.empty())
        std::cerr << "Data dir  : " << g_data_dir << "\n";

    // load initial file if provided (blocking)
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

    // clean shutdown: join any in-progress load thread
    {
        std::lock_guard<std::mutex> lk(g_load_mtx);
        if (g_load_thread.joinable()) g_load_thread.join();
    }
    return 0;
}
