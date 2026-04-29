#pragma once
// =========================================================================
// viewer_server.h — Unified HTTP/WebSocket server for PRad-II event viewer
//
// Combines file-based viewing and online ET monitoring into a single server.
// Mode switching between "idle", "file", and "online" via API or user actions.
// =========================================================================

#include "data_source.h"
#include "app_state.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "TdcData.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <atomic>
#include <deque>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

// ── Shared types ─────────────────────────────────────────────────────────

struct FileData {
    std::string filepath;
    int event_count = 0;
    DataSourceCaps caps;
};

struct Progress {
    std::atomic<bool> loading{false};
    std::atomic<int>  phase{0};      // 0=idle, 1=indexing, 2=histograms
    std::atomic<int>  current{0};
    std::atomic<int>  total{0};
    std::string       target_file;
    mutable std::mutex mtx;

    nlohmann::json toJson() const;
    void setFile(const std::string &f);
};

struct RingEntry {
    int seq;
    std::string json_str;       // pre-encoded event JSON
    std::string cluster_str;    // pre-encoded cluster JSON
    std::string gem_apv_str;    // pre-encoded GEM per-APV waveform JSON
    // Pre-compressed gzip bytes for the big payload — deflated once per
    // event in the ET reader thread so each viewer's HTTP request just
    // copies the cached blob instead of re-running zlib.  Empty when
    // gem_apv_str is empty (e.g. GEM disabled).
    std::string gem_apv_gz;

    // Raw event copies kept so /api/hist_config can recompute cluster_str
    // (and re-encode json_str) under a new time/threshold window without
    // waiting for new events to arrive.  ~8MB per entry — bounded by ring_size_.
    std::shared_ptr<fdec::EventData>    event_data;
    std::shared_ptr<ssp::SspEventData>  ssp_data;
};

// ── ViewerServer ─────────────────────────────────────────────────────────

class ViewerServer {
public:
    using WsServer = websocketpp::server<websocketpp::config::asio>;

    struct Config {
        std::string database_dir;
        std::string resource_dir;
        std::string data_dir;           // file browsing root (empty = disabled)
        std::string daq_config_file;    // empty = auto-find in database_dir
        std::string config_file;        // empty = auto-find in database_dir
        std::string initial_file;       // .evio file to open on startup
        int    port         = 5051;
        bool   hist_enabled = false;
        bool   start_online = false;    // connect ET on startup
        bool   interactive  = false;    // enable stdin command loop
        std::string filter_file;        // external filter JSON (-f)
    };

    ViewerServer();
    ~ViewerServer();

    // Initialize application state. Must be called before run/startAsync.
    void init(const Config &cfg);

    // Run the server (blocking). Loads initial file, then serves.
    void run();

    // Start server in a background thread. Returns the actual port.
    int startAsync(int port = 0);

    // Stop the server and all background threads.
    void stop();

    // Load a file by absolute path. Switches to file mode.
    // Non-blocking: spawns a background load thread.
    void loadFile(const std::string &path, bool hist);

    int  port() const { return port_; }
    std::string mode() const;
    bool isLoading() const { return progress_.loading.load(); }
    nlohmann::json getProgress() const { return progress_.toJson(); }

    // Active AppState for the current mode.
    AppState &activeApp();

private:
    // ── Mode ─────────────────────────────────────────────────────────────
    enum class Mode { Idle, File, Online };
    std::atomic<Mode> mode_{Mode::Idle};
    std::mutex mode_mtx_;       // serialises mode transitions

    void setMode(Mode m);       // set + broadcast

    // ── Dual AppState (file vs online, never mixed) ──────────────────────
    AppState    app_file_;
    AppState    app_online_;
    Config      cfg_;
    std::string res_dir_;
    int         port_ = 0;
    std::atomic<bool> running_{true};

    // ── WebSocket ────────────────────────────────────────────────────────
    std::unique_ptr<WsServer> server_;
    std::thread server_thread_;
    std::set<websocketpp::connection_hdl,
             std::owner_less<websocketpp::connection_hdl>> ws_clients_;
    std::mutex ws_mtx_;

    void wsBroadcast(const std::string &msg);
    void handleWsMessage(websocketpp::connection_hdl hdl,
                         const std::string &payload);

    // ── Tagger live stream (binary broadcast to subscribed WebSocket clients) ─
    // Zero-cost when no one is subscribed: the ET reader only runs the TDC
    // decoder when tagger_subs_count_ > 0.  Frames are a 24-byte header
    // followed by N × 16B packed BinHit records — see viewer_server_et.cpp
    // for the exact layout (the Python client in scripts/tagger_viewer.py
    // mirrors it).
    std::set<websocketpp::connection_hdl,
             std::owner_less<websocketpp::connection_hdl>> tagger_subs_;
    std::mutex                 tagger_subs_mtx_;
    std::atomic<int>           tagger_subs_count_{0};
    std::atomic<uint64_t>      tagger_dropped_frames_{0};  // incremented on per-subscriber send failure

    void taggerSubscribe(websocketpp::connection_hdl hdl);
    void taggerUnsubscribe(websocketpp::connection_hdl hdl);
    void taggerBroadcastBinary(const void *data, size_t nbytes);

    // ── File mode ────────────────────────────────────────────────────────
    std::shared_ptr<FileData> file_data_;
    std::mutex file_data_mtx_;

    std::unique_ptr<DataSource> data_source_;
    std::mutex data_source_mtx_;
    std::unordered_map<int, uint32_t> crate_to_roc_;  // for ROOT data sources

    mutable Progress progress_;
    std::atomic<bool> hist_enabled_{false};
    std::thread load_thread_;
    std::mutex load_mtx_;

    // On-demand accumulation: mirrors the online-mode logic for file browsing.
    // - Preprocessed (hist_enabled_): all events already processed by
    //   buildHistograms(); further calls are no-ops.
    // - Not preprocessed: processEvent/processGemEvent are called once per
    //   event as the user browses (deduped by ondemand_processed_).
    // Any new accumulation added to processEvent() automatically follows
    // this pattern — no per-endpoint code is needed.
    std::unordered_set<int> ondemand_processed_;
    std::mutex ondemand_mtx_;
    void accumulate(int ev1, fdec::EventData &event, ssp::SspEventData *ssp);

    // ── Filters ──────────────────────────────────────────────────────────
    std::vector<int> filtered_indices_;   // 1-based event indices passing filter
    void buildFilteredIndex();
    std::string applyFilter(const nlohmann::json &fj);
    void clearFilter();

    void buildHistograms();
    void loadFileInternal(const std::string &filepath);

    std::string decodeRawEvent(int ev1, fdec::EventData &event,
                               ssp::SspEventData *ssp_evt = nullptr);
    nlohmann::json decodeEvent(int ev1);
    nlohmann::json computeClusters(int ev1);

    // ── Online mode (ET) ─────────────────────────────────────────────────
#ifdef WITH_ET
    struct EtCfg {
        std::string host    = "localhost";
        int         port    = 11111;
        std::string et_file = "/tmp/et_sys_prad2";
        std::string station = "prad2_monitor";
    } et_cfg_;

    int ring_size_ = 20;
    std::deque<RingEntry> ring_;
    std::mutex ring_mtx_;

    std::atomic<bool> et_active_{false};
    std::atomic<bool> et_connected_{false};
    std::atomic<int>  et_generation_{0};    // bumped to trigger reconnect
    std::thread et_thread_;

    void etReaderThread();
    void sleepMs(int ms);

    // DAQ livetime (TS, percent).  <0 = not available.
    // Written by livetimePollThread() — an optional shell-command poll
    // (typical: "caget -t <channel>") configured via AppState::livetime_cmd.
    // Disabled when livetime_cmd is empty.  Avoids a build-time EPICS
    // dependency by shelling out to whatever tool the host provides.
    // The "measured" companion lives on AppState::measured_livetime and is
    // populated from the DSC2 scaler bank in the EVIO stream.
    std::atomic<double> livetime_{-1.0};
    std::thread         livetime_thread_;
    void                livetimePollThread();
#endif

    void joinAll();      // join all background threads (safe to call multiple times)
    void commandLoop();  // interactive stdin command loop

    // ── HTTP / resource handling ─────────────────────────────────────────
    void setupServer(int port);
    bool serveResource(const std::string &uri, WsServer::connection_ptr con);
    nlohmann::json listFiles(const std::string &subdir = "");
    std::string resolveDataFile(const std::string &relpath);
    nlohmann::json buildConfig();
    nlohmann::json handleElogPost(const std::string &body);
    void onHttp(WsServer *srv, websocketpp::connection_hdl hdl);
};
