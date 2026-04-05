// =========================================================================
// viewer_server.cpp — Unified server implementation
//
// Merges file-based event viewing and online ET monitoring into one server.
// =========================================================================

#include "viewer_server.h"

#ifdef WITH_ET
#include "EtChannel.h"
#endif

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <chrono>

using json = nlohmann::json;
namespace fs = std::filesystem;
using namespace evc;

// =========================================================================
// Progress
// =========================================================================

json Progress::toJson() const
{
    std::lock_guard<std::mutex> lk(mtx);
    return {{"loading", loading.load()},
            {"phase", phase == 1 ? "indexing" : phase == 2 ? "histograms" : "idle"},
            {"current", current.load()}, {"total", total.load()}, {"file", target_file}};
}

void Progress::setFile(const std::string &f)
{
    std::lock_guard<std::mutex> lk(mtx);
    target_file = f;
}


// =========================================================================
// ViewerServer — lifecycle
// =========================================================================

ViewerServer::ViewerServer() = default;

ViewerServer::~ViewerServer()
{
    stop();
    joinAll();
}

void ViewerServer::joinAll()
{
#ifdef WITH_ET
    if (et_thread_.joinable()) et_thread_.join();
#endif
    { std::lock_guard<std::mutex> lk(load_mtx_);
      if (load_thread_.joinable()) load_thread_.join(); }
    if (server_thread_.joinable()) server_thread_.join();
}

void ViewerServer::init(const Config &cfg)
{
    cfg_ = cfg;
    res_dir_ = cfg.resource_dir;
    hist_enabled_ = cfg.hist_enabled;

    const auto &db_dir = cfg.database_dir;

    // --- resolve DAQ config path ---
    std::string daq_cfg_file = cfg.daq_config_file;
    if (daq_cfg_file.empty())
        daq_cfg_file = findFile("daq_config.json", db_dir);

    // --- load ET config from config.json "online" section ---
#ifdef WITH_ET
    {
        std::string cf = cfg.config_file;
        if (cf.empty()) cf = findFile("config.json", db_dir);
        if (cf.empty()) cf = findFile("online_config.json", db_dir);
        std::string s = readFile(cf);
        if (!s.empty()) {
            auto j = json::parse(s, nullptr, false);
            if (j.contains("online")) {
                auto &e = j["online"];
                if (e.contains("et_host"))    et_cfg_.host    = e["et_host"];
                if (e.contains("et_port"))    et_cfg_.port    = e["et_port"];
                if (e.contains("et_file"))    et_cfg_.et_file = e["et_file"];
                if (e.contains("et_station")) et_cfg_.station = e["et_station"];
                if (e.contains("ring_buffer_size")) ring_size_ = e["ring_buffer_size"];
            } else if (j.contains("et")) {
                auto &e = j["et"];
                if (e.contains("host"))    et_cfg_.host    = e["host"];
                if (e.contains("port"))    et_cfg_.port    = e["port"];
                if (e.contains("et_file")) et_cfg_.et_file = e["et_file"];
                if (e.contains("station")) et_cfg_.station = e["station"];
                if (j.contains("ring_buffer_size")) ring_size_ = j["ring_buffer_size"];
            }
        }
    }
#endif

    // --- initialize both AppState instances (same config, separate accumulators) ---
    app_file_.init(db_dir, daq_cfg_file, cfg.config_file);
    app_online_.init(db_dir, daq_cfg_file, cfg.config_file);

    // --- build base_config JSON for /api/config ---
    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(findFile("hycal_modules.json", db_dir));
      if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    {
        std::string daq_fn = "daq_map.json";
        if (!daq_cfg_file.empty()) {
            std::ifstream dcf(daq_cfg_file);
            if (dcf.is_open()) {
                auto dcj = json::parse(dcf, nullptr, false, true);
                if (dcj.contains("daq_map_file"))
                    daq_fn = dcj["daq_map_file"].get<std::string>();
            }
        }
        std::string s = readFile(findFile(daq_fn, db_dir));
        if (!s.empty()) daq_j = json::parse(s, nullptr, false);
    }
    json base_cfg = {
        {"modules", modules_j}, {"daq", daq_j}, {"crate_roc", app_file_.crate_roc_json},
    };
    app_file_.base_config = std::move(base_cfg);
    app_online_.base_config = app_file_.base_config;

    // --- build crate→ROC tag map for ROOT data sources ---
    for (auto &[k, v] : app_file_.crate_roc_json.items())
        crate_to_roc_[std::stoi(k)] = v.get<uint32_t>();

    // --- validate resources ---
    if (readFile(res_dir_ + "/viewer.html").empty())
        std::cerr << "Warning: viewer.html not found in " << res_dir_ << "\n";

    std::cerr << "Database  : " << db_dir << "\n"
              << "Resources : " << res_dir_ << "\n";
    if (!cfg.data_dir.empty())
        std::cerr << "Data dir  : " << cfg.data_dir << "\n";
#ifdef WITH_ET
    std::cerr << "ET target : " << et_cfg_.host << ":" << et_cfg_.port << "\n"
              << "Ring buf  : " << ring_size_ << " events\n";
#endif
}

// =========================================================================
// Server setup & run
// =========================================================================

void ViewerServer::setupServer(int port)
{
    server_ = std::make_unique<WsServer>();
    server_->set_access_channels(websocketpp::log::alevel::none);
    server_->set_error_channels(websocketpp::log::elevel::warn |
                                websocketpp::log::elevel::rerror);
    server_->init_asio();
    server_->set_reuse_addr(true);

    server_->set_http_handler([this](websocketpp::connection_hdl hdl) {
        onHttp(server_.get(), hdl);
    });
    server_->set_open_handler([this](websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(ws_mtx_);
        ws_clients_.insert(hdl);
    });
    server_->set_close_handler([this](websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(ws_mtx_);
        ws_clients_.erase(hdl);
    });

    // bind — retry a few ports if port==0
    if (port == 0) {
        for (int p = 15050; p < 15150; ++p) {
            try { server_->listen(p); port_ = p; break; }
            catch (...) { continue; }
        }
        if (port_ == 0) { server_->listen(15050); port_ = 15050; }
    } else {
        server_->listen(port);
        port_ = port;
    }

    server_->start_accept();
}

void ViewerServer::run()
{
    setupServer(cfg_.port);

    // start ET reader thread (sleeps until et_active_)
#ifdef WITH_ET
    et_thread_ = std::thread([this]() { etReaderThread(); });
#endif

    // load initial file (blocking)
    if (!cfg_.initial_file.empty()) {
        mode_ = Mode::File;
        hist_enabled_ = cfg_.hist_enabled;
        loadFileInternal(cfg_.initial_file);
    }

    // if --et, switch to online mode
    if (cfg_.start_online) {
        mode_ = Mode::Online;
#ifdef WITH_ET
        et_active_ = true;
#endif
    }

    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    std::cout << "Server at http://localhost:" << port_ << "\n"
              << "  Mode: " << mode() << "\n"
              << "  " << (data ? data->event_count : 0) << " events"
              << (hist_enabled_ ? ", histograms enabled" : "")
              << (!cfg_.data_dir.empty() ? ", file browser enabled" : "")
              << (cfg_.interactive ? "\n  Type 'help' for commands, " : "\n  ")
              << "Ctrl+C to stop\n";

    // interactive stdin command loop (detached — will be killed on exit)
    if (cfg_.interactive)
        std::thread([this]() { commandLoop(); }).detach();

    server_->run();

    // cleanup — stop() already set flags; joinAll() waits for threads
    stop();
    joinAll();
}

int ViewerServer::startAsync(int port)
{
    if (port != 0) cfg_.port = port;
    setupServer(cfg_.port);

    // start ET reader thread
#ifdef WITH_ET
    et_thread_ = std::thread([this]() { etReaderThread(); });
#endif

    // load initial file (async)
    if (!cfg_.initial_file.empty())
        loadFile(cfg_.initial_file, cfg_.hist_enabled);

    // start online if requested
    if (cfg_.start_online) {
        mode_ = Mode::Online;
#ifdef WITH_ET
        et_active_ = true;
#endif
    }

    // run server in background
    server_thread_ = std::thread([this]() { server_->run(); });
    return port_;
}

void ViewerServer::stop()
{
    running_ = false;
#ifdef WITH_ET
    et_active_ = false;
#endif

    if (server_) {
        try { server_->stop_listening(); server_->stop(); } catch (...) {}
    }

    // Thread joins are handled by run()/startAsync() cleanup — not here.
    // stop() may be called from a signal handler where blocking on join()
    // would hang if a thread is stuck in a library call (e.g. et_open).
}

// =========================================================================
// Mode switching
// =========================================================================

std::string ViewerServer::mode() const
{
    switch (mode_.load()) {
    case Mode::File:   return "file";
    case Mode::Online: return "online";
    default:           return "idle";
    }
}

AppState &ViewerServer::activeApp()
{
    return mode_.load() == Mode::Online ? app_online_ : app_file_;
}

void ViewerServer::setMode(Mode m)
{
    mode_ = m;
    wsBroadcast(json({{"type", "mode_changed"}, {"mode", mode()}}).dump());
}

// =========================================================================
// WebSocket broadcast
// =========================================================================

void ViewerServer::wsBroadcast(const std::string &msg)
{
    std::lock_guard<std::mutex> lk(ws_mtx_);
    for (auto &hdl : ws_clients_) {
        try { server_->send(hdl, msg, websocketpp::frame::opcode::text); }
        catch (...) {}
    }
}

// =========================================================================
// File mode — loading
// =========================================================================

void ViewerServer::loadFile(const std::string &path, bool hist)
{
    {
        std::lock_guard<std::mutex> lk(mode_mtx_);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online)
            et_active_ = false;
#endif
        hist_enabled_ = hist;
        setMode(Mode::File);
    }

    std::lock_guard<std::mutex> lk(load_mtx_);
    if (load_thread_.joinable()) load_thread_.join();
    load_thread_ = std::thread([this, path]() { loadFileInternal(path); });
}

void ViewerServer::loadFileInternal(const std::string &filepath)
{
    progress_.loading = true;
    progress_.setFile(filepath);
    progress_.phase = 0; progress_.current = 0; progress_.total = 0;

    auto data = std::make_shared<FileData>();
    data->filepath = filepath;

    // reset on-demand accumulation tracking
    { std::lock_guard<std::mutex> lk(ondemand_mtx_); ondemand_processed_.clear(); }

    std::cerr << "Loading: " << filepath << "\n";

    // create and open the appropriate data source
    progress_.phase = 1;
    auto source = createDataSource(filepath, app_file_.daq_cfg, crate_to_roc_);
    if (!source) {
        std::cerr << "  Error: unsupported file type\n";
        progress_.loading = false; return;
    }
    std::string err = source->open(filepath);
    if (!err.empty()) {
        std::cerr << "  Error: " << err << "\n";
        progress_.loading = false; return;
    }

    data->event_count = source->eventCount();
    data->caps = source->capabilities();
    std::cerr << "  Indexed " << data->event_count << " events"
              << " (source: " << data->caps.source_type << ")\n";

    // install the new data source before building histograms
    // (buildHistograms reads data_source_; loadFileInternal runs on the load
    // thread, and HTTP threads acquire data_source_mtx_ for random access)
    { std::lock_guard<std::mutex> lk(data_source_mtx_); data_source_ = std::move(source); }

    // always reset accumulators when opening a new file
    app_file_.clearHistograms();
    app_file_.clearLms();
    app_file_.clearEpics();

    if (hist_enabled_) {
        progress_.total = data->event_count;
        buildHistograms();
    }

    { std::lock_guard<std::mutex> lk(file_data_mtx_); file_data_ = data; }

    progress_.loading = false;
    progress_.phase = 0;
    std::cerr << "  Ready\n";
}

void ViewerServer::buildHistograms()
{
    app_file_.clearHistograms();
    app_file_.clearLms();
    app_file_.clearEpics();

    if (!data_source_) return;

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = app_file_.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;

    progress_.phase = 2; progress_.current = 0;

    data_source_->iterateAll(
        // physics events (EVIO / ROOT raw)
        [&](int idx, fdec::EventData &event, ssp::SspEventData *ssp) {
            progress_.current = app_file_.events_processed.load() + 1;
            app_file_.processEvent(event, ana, wres);
            if (ssp) app_file_.processGemEvent(*ssp);
        },
        // recon events (ROOT recon)
        [&](int idx, const ReconEventData &recon) {
            progress_.current = app_file_.events_processed.load() + 1;
            app_file_.processReconEvent(recon);
        },
        // control events (sync/prestart/go)
        [&](uint32_t unix_time, uint64_t last_ti_ts) {
            if (app_file_.sync_unix == 0)
                app_file_.recordSyncTime(unix_time, last_ti_ts);
        },
        // EPICS events
        [&](const std::string &text, int32_t ev_num, uint64_t ts) {
            app_file_.processEpics(text, app_file_.events_processed.load(), ts);
        }
    );

    std::cerr << "  Histograms: " << app_file_.events_processed.load() << " events"
              << ", clusters: " << app_file_.cluster_events_processed
              << ", LMS: " << app_file_.lms_events.load() << "\n";
}

// =========================================================================
// File mode — event decoding
// =========================================================================

std::string ViewerServer::decodeRawEvent(int ev1, fdec::EventData &event,
                                         ssp::SspEventData *ssp_evt)
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    if (!data) return "no file loaded";
    int idx = ev1 - 1;
    if (idx < 0 || idx >= data->event_count) return "event out of range";

    std::lock_guard<std::mutex> lk(data_source_mtx_);
    if (!data_source_) return "no data source";
    return data_source_->decodeEvent(idx, event, ssp_evt);
}

// Central accumulation gate for file mode.
// Mirrors online-mode logic: processEvent + processGemEvent are called once
// per event.  Preprocessed files skip this (buildHistograms already did it).
// Any new accumulation added to processEvent() works automatically.
void ViewerServer::accumulate(int ev1, fdec::EventData &event,
                              ssp::SspEventData *ssp)
{
    if (hist_enabled_) return;  // preprocessed — already done

    std::lock_guard<std::mutex> lk(ondemand_mtx_);
    if (!ondemand_processed_.insert(ev1).second) return;  // already seen

    if (ssp) app_file_.processGemEvent(*ssp);

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = app_file_.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    app_file_.processEvent(event, ana, wres);
}

json ViewerServer::decodeEvent(int ev1)
{
    // check if this is a recon source (no per-channel data)
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    if (data && data->caps.source_type == "root_recon") {
        // return minimal event info + empty channels
        ReconEventData recon;
        { std::lock_guard<std::mutex> lk(data_source_mtx_);
          if (!data_source_ || !data_source_->decodeReconEvent(ev1 - 1, recon))
              return {{"error", "decode error"}}; }
        return {{"event", ev1}, {"channels", json::object()},
                {"event_number", recon.event_num},
                {"trigger_bits", recon.trigger_bits}};
    }

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    std::string err = decodeRawEvent(ev1, event, ssp_ptr.get());
    if (!err.empty()) return {{"error", err}};

    accumulate(ev1, event, ssp_ptr.get());

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = app_file_.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    return app_file_.encodeEventJson(event, ev1, ana, wres);
}

json ViewerServer::computeClusters(int ev1)
{
    // recon source: return pre-computed clusters
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    if (data && data->caps.source_type == "root_recon") {
        ReconEventData recon;
        { std::lock_guard<std::mutex> lk(data_source_mtx_);
          if (!data_source_ || !data_source_->decodeReconEvent(ev1 - 1, recon))
              return {{"error", "decode error"}}; }
        return app_file_.encodeReconClustersJson(recon, ev1);
    }

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    std::string err = decodeRawEvent(ev1, event, ssp_ptr.get());
    if (!err.empty()) return {{"error", err}};

    accumulate(ev1, event, ssp_ptr.get());

    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = app_file_.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    return app_file_.computeClustersJson(event, ev1, ana, wres);
}

// =========================================================================
// Interactive command loop (stdin)
// =========================================================================

void ViewerServer::commandLoop()
{
    std::string line;
    while (running_ && std::getline(std::cin, line)) {
        // trim whitespace
        auto s = line.find_first_not_of(" \t");
        if (s == std::string::npos) continue;
        line = line.substr(s, line.find_last_not_of(" \t") - s + 1);
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string cmd;
        iss >> cmd;

        if (cmd == "help" || cmd == "?") {
            std::cout << "Commands:\n"
                      << "  status            — show server status\n"
                      << "  load <file> [H]   — load an evio file (H=1 for histograms)\n"
#ifdef WITH_ET
                      << "  online            — switch to ET/online mode\n"
#endif
                      << "  offline           — switch to file mode\n"
                      << "  clear <what>      — clear hist, lms, or epics\n"
                      << "  quit / exit       — stop the server\n";
        }
        else if (cmd == "status") {
            std::shared_ptr<FileData> data;
            { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
            std::cout << "Mode: " << mode() << "\n";
            if (data)
                std::cout << "File: " << data->filepath << " (" << data->event_count << " events)\n";
#ifdef WITH_ET
            std::cout << "ET: " << (et_connected_.load() ? "connected" : "not connected")
                      << " (" << et_cfg_.host << ":" << et_cfg_.port << ")\n";
            std::cout << "Events processed: " << app_online_.events_processed.load() << "\n";
#endif
        }
        else if (cmd == "load") {
            std::string path;
            int hist = 0;
            iss >> path >> hist;
            if (path.empty()) {
                std::cout << "Usage: load <filepath> [hist=0|1]\n";
            } else {
                loadFile(path, hist != 0);
                std::cout << "Loading " << path << (hist ? " (with histograms)" : "") << "\n";
            }
        }
#ifdef WITH_ET
        else if (cmd == "online") {
            std::lock_guard<std::mutex> lk(mode_mtx_);
            if (mode_.load() == Mode::Online)
                et_generation_++;
            et_active_ = true;
            setMode(Mode::Online);
            std::cout << "Switched to online mode\n";
        }
#endif
        else if (cmd == "offline") {
            std::lock_guard<std::mutex> lk(mode_mtx_);
            if (mode_.load() == Mode::Online) {
#ifdef WITH_ET
                et_active_ = false;
#endif
                std::shared_ptr<FileData> data;
                { std::lock_guard<std::mutex> lk2(file_data_mtx_); data = file_data_; }
                setMode(data ? Mode::File : Mode::Idle);
            }
            std::cout << "Switched to " << mode() << " mode\n";
        }
        else if (cmd == "clear") {
            std::string what;
            iss >> what;
            if (what == "hist") {
                activeApp().clearHistograms();
                wsBroadcast("{\"type\":\"hist_cleared\"}");
                std::cout << "Histograms cleared\n";
            } else if (what == "lms") {
                activeApp().clearLms();
                wsBroadcast("{\"type\":\"lms_cleared\"}");
                std::cout << "LMS cleared\n";
            } else if (what == "epics") {
                activeApp().clearEpics();
                wsBroadcast("{\"type\":\"epics_cleared\"}");
                std::cout << "EPICS cleared\n";
            } else {
                std::cout << "Usage: clear hist|lms|epics\n";
            }
        }
        else if (cmd == "quit" || cmd == "exit") {
            std::cout << "Stopping server...\n";
            stop();
            break;
        }
        else {
            std::cout << "Unknown command: " << cmd << " (type 'help')\n";
        }
    }
}

// =========================================================================
// Online mode — ET reader thread
// =========================================================================

#ifdef WITH_ET

void ViewerServer::sleepMs(int ms)
{
    for (int elapsed = 0; elapsed < ms && running_ && et_active_; elapsed += 100)
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void ViewerServer::etReaderThread()
{
    EtChannel ch;
    ch.SetConfig(app_online_.daq_cfg);
    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    auto &ssp_evt = *ssp_ptr;
    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = app_online_.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    uint64_t last_ti_ts = 0;

    while (running_) {
        // sleep until activated
        while (running_ && !et_active_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        if (!running_) break;

        int retry_ms = 3000;
        const int max_retry = 30000;
        int retry_count = 0;
        auto retry_start = std::chrono::steady_clock::now();
        int gen = et_generation_.load();

        while (running_ && et_active_ && et_generation_.load() == gen) {
            if (retry_count == 0) {
                std::cerr << "ET: connecting to " << et_cfg_.host << ":" << et_cfg_.port
                          << "  " << et_cfg_.et_file << " ...\n";
                retry_start = std::chrono::steady_clock::now();
            }

            if (ch.Connect(et_cfg_.host, et_cfg_.port, et_cfg_.et_file)
                    != status::success) {
                retry_count++;
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::steady_clock::now() - retry_start).count();
                std::cerr << "\rET: waiting for ET system... "
                          << retry_count << " attempts, " << elapsed << "s elapsed   "
                          << std::flush;
                wsBroadcast("{\"type\":\"status\",\"connected\":false,\"waiting\":true,"
                            "\"retries\":" + std::to_string(retry_count) + "}");
                sleepMs(retry_ms);
                retry_ms = std::min(retry_ms * 2, max_retry);
                continue;
            }

            if (retry_count > 0) std::cerr << "\n";

            if (ch.Open(et_cfg_.station) != status::success) {
                std::cerr << "ET: station open failed, retrying...\n";
                ch.Disconnect();
                sleepMs(retry_ms);
                retry_ms = std::min(retry_ms * 2, max_retry);
                continue;
            }

            retry_ms = 3000;
            retry_count = 0;
            et_connected_ = true;
            wsBroadcast("{\"type\":\"status\",\"connected\":true}");
            std::cerr << "ET: connected, reading events\n";

            int gen = et_generation_.load();
            auto last_ring_push = std::chrono::steady_clock::now();
            constexpr auto ring_interval = std::chrono::milliseconds(50);
            auto last_lms_notify = last_ring_push;
            constexpr auto lms_notify_interval = std::chrono::milliseconds(200);

            while (running_ && et_active_ && et_generation_.load() == gen) {
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

                if (app_online_.sync_unix == 0) {
                    uint32_t ct = ch.GetControlTime();
                    if (ct != 0) app_online_.recordSyncTime(ct, last_ti_ts);
                }

                if (ch.GetEventType() == EventType::Epics) {
                    std::string text = ch.ExtractEpicsText();
                    if (!text.empty()) {
                        int seq = app_online_.events_processed.load();
                        app_online_.processEpics(text, seq, last_ti_ts);
                        wsBroadcast("{\"type\":\"epics_event\",\"count\":" +
                                    std::to_string(app_online_.epics_events.load()) + "}");
                    }
                }

                for (int i = 0; i < ch.GetNEvents(); ++i) {
                    ssp_evt.clear();
                    if (!ch.DecodeEvent(i, event, &ssp_evt)) continue;
                    last_ti_ts = event.info.timestamp;

                    app_online_.processGemEvent(ssp_evt);
                    app_online_.processEvent(event, ana, wres);

                    int seq = app_online_.events_processed.load();

                    if (app_online_.lms_trigger.accept != 0 &&
                        app_online_.lms_trigger(event.info.trigger_bits)) {
                        auto now = std::chrono::steady_clock::now();
                        if (now - last_lms_notify >= lms_notify_interval) {
                            last_lms_notify = now;
                            wsBroadcast("{\"type\":\"lms_event\",\"count\":" +
                                        std::to_string(app_online_.lms_events.load()) + "}");
                        }
                    }

                    auto now = std::chrono::steady_clock::now();
                    if (now - last_ring_push >= ring_interval) {
                        last_ring_push = now;

                        std::string evjson = app_online_.encodeEventJson(
                            event, seq, ana, wres, true).dump();
                        std::string cljson = app_online_.computeClustersJson(
                            event, seq, ana, wres).dump();

                        {
                            std::lock_guard<std::mutex> lk(ring_mtx_);
                            ring_.push_back({seq, std::move(evjson),
                                             std::move(cljson)});
                            while ((int)ring_.size() > ring_size_)
                                ring_.pop_front();
                        }

                        wsBroadcast("{\"type\":\"new_event\",\"seq\":" +
                                    std::to_string(seq) + "}");
                    }
                }
            }

            et_connected_ = false;
            ch.Close();
            ch.Disconnect();
            wsBroadcast("{\"type\":\"status\",\"connected\":false}");

            if (running_ && et_active_) {
                std::cerr << "ET: disconnected, retrying in "
                          << retry_ms / 1000 << "s\n";
                sleepMs(retry_ms);
            }
        }
    }
}

#endif // WITH_ET

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
                    entries.push_back({{"type", "dir"}, {"name", rel}, {"count", count}});
            } else if (entry.is_regular_file()) {
                auto fn = entry.path().filename().string();
                if (fn.find(".evio") == std::string::npos &&
                    fn.find(".root") == std::string::npos)
                    continue;
                auto sz = entry.file_size();
                entries.push_back({{"type", "file"}, {"name", rel},
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
