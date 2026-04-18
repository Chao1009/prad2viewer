#include "viewer_server.h"

#ifdef WITH_ET
#include "EtChannel.h"
#endif

#include <cstring>

using namespace evc;
using json = nlohmann::json;

#ifdef WITH_ET

// =========================================================================
// TDC live-stream frame format
//
// Header (little-endian, 24 bytes, matches the dtype expected by
// scripts/tdc_viewer.py):
//
//   char     magic[4];       // "TDC1"
//   uint32_t flags;           // bit 0 set when dropped_count > 0
//   uint32_t n_hits;          // number of 16-byte BinHit records that follow
//   uint32_t first_seq;       // event seq of the first hit in the frame
//   uint32_t last_seq;        // event seq of the last hit
//   uint32_t dropped;         // total dropped frames (since server start)
//
// Hits use a 16-byte packed layout (same as scripts/tdc_viewer.py RAW_DTYPE):
//
//   uint32_t event_num;
//   uint32_t trigger_bits;
//   uint16_t roc_tag;
//   uint8_t  slot;
//   uint8_t  channel_edge;    // bit 7 = edge, bits 6:0 = channel
//   uint32_t tdc;
// =========================================================================
namespace {

constexpr size_t TDC_HIT_SIZE     = 16;
constexpr size_t TDC_HDR_SIZE     = 24;
constexpr uint32_t TDC_BATCH_MAX  = 256;      // flush at this many hits
constexpr auto     TDC_BATCH_MS   = std::chrono::milliseconds(10);

#pragma pack(push, 1)
struct TdcBinHit {
    uint32_t event_num;
    uint32_t trigger_bits;
    uint16_t roc_tag;
    uint8_t  slot;
    uint8_t  channel_edge;
    uint32_t tdc;
};
#pragma pack(pop)
static_assert(sizeof(TdcBinHit) == TDC_HIT_SIZE, "TdcBinHit layout");

} // namespace

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
    auto tdc_ptr = std::make_unique<tdc::TdcEventData>();
    auto &tdc_evt = *tdc_ptr;
    fdec::WaveAnalyzer ana;
    ana.cfg.min_peak_ratio = app_online_.hist_cfg.min_peak_ratio;
    fdec::WaveResult wres;
    uint64_t last_ti_ts = 0;

    // TDC batch state (local to the thread — only this thread writes it).
    // Pre-allocate enough for TDC_BATCH_MAX hits + header.
    std::vector<uint8_t> tdc_batch;
    tdc_batch.reserve(TDC_HDR_SIZE + TDC_BATCH_MAX * TDC_HIT_SIZE);
    uint32_t tdc_batch_hits = 0;
    uint32_t tdc_batch_first_seq = 0;
    uint32_t tdc_batch_last_seq = 0;
    auto tdc_batch_last_flush = std::chrono::steady_clock::now();

    auto tdc_flush = [&]() {
        if (tdc_batch_hits == 0) { tdc_batch_last_flush = std::chrono::steady_clock::now(); return; }
        // Fill header in place.
        uint8_t *p = tdc_batch.data();
        std::memcpy(p + 0,  "TDC1", 4);
        uint32_t drops = static_cast<uint32_t>(tdc_dropped_frames_.load());
        uint32_t flags = (drops > 0) ? 1u : 0u;
        std::memcpy(p + 4,  &flags,               4);
        std::memcpy(p + 8,  &tdc_batch_hits,      4);
        std::memcpy(p + 12, &tdc_batch_first_seq, 4);
        std::memcpy(p + 16, &tdc_batch_last_seq,  4);
        std::memcpy(p + 20, &drops,               4);
        tdcBroadcastBinary(tdc_batch.data(),
                           TDC_HDR_SIZE + tdc_batch_hits * TDC_HIT_SIZE);
        // Reset batch (keep allocation).
        tdc_batch.resize(TDC_HDR_SIZE);
        tdc_batch_hits = 0;
        tdc_batch_first_seq = 0;
        tdc_batch_last_seq = 0;
        tdc_batch_last_flush = std::chrono::steady_clock::now();
    };

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

                // skip monitoring events (TI only, no waveforms) — same
                // as file viewer path in evio_data_source.cpp
                if (app_online_.daq_cfg.is_monitoring(ch.GetEvHeader().tag))
                    continue;

                if (app_online_.sync_unix == 0) {
                    auto et = ch.GetEventType();
                    if (et == EventType::Prestart || et == EventType::Go ||
                        et == EventType::End      || et == EventType::Sync)
                    {
                        const auto &s = ch.Sync();
                        if (s.unix_time != 0)
                            app_online_.recordSyncTime(s.unix_time, last_ti_ts);
                    }
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
                    const bool want_tdc =
                        tdc_subs_count_.load(std::memory_order_relaxed) > 0;
                    tdc::TdcEventData *tdc_arg = want_tdc ? &tdc_evt : nullptr;
                    if (tdc_arg) tdc_arg->clear();
                    if (!ch.DecodeEvent(i, event, &ssp_evt, nullptr, tdc_arg)) continue;
                    last_ti_ts = event.info.timestamp;

                    app_online_.processGemEvent(ssp_evt);
                    app_online_.processEvent(event, ana, wres);

                    int seq = app_online_.events_processed.load();

                    // --- live TDC stream: batch + broadcast ------------------
                    if (want_tdc && tdc_evt.n_hits > 0) {
                        if (tdc_batch_hits == 0) {
                            // Grow to header size on first hit of the batch.
                            tdc_batch.resize(TDC_HDR_SIZE);
                            tdc_batch_first_seq = static_cast<uint32_t>(seq);
                        }
                        const uint32_t evnum = static_cast<uint32_t>(event.info.event_number);
                        const uint32_t tbits = event.info.trigger_bits;
                        for (int h = 0; h < tdc_evt.n_hits; ++h) {
                            const auto &src = tdc_evt.hits[h];
                            TdcBinHit bh;
                            bh.event_num    = evnum;
                            bh.trigger_bits = tbits;
                            bh.roc_tag      = static_cast<uint16_t>(src.roc_tag);
                            bh.slot         = src.slot;
                            bh.channel_edge =
                                static_cast<uint8_t>(((src.edge & 0x1) << 7) |
                                                      (src.channel & 0x7F));
                            bh.tdc          = src.value;
                            size_t off = tdc_batch.size();
                            tdc_batch.resize(off + TDC_HIT_SIZE);
                            std::memcpy(tdc_batch.data() + off, &bh, TDC_HIT_SIZE);
                            ++tdc_batch_hits;
                            if (tdc_batch_hits >= TDC_BATCH_MAX) { tdc_flush(); break; }
                        }
                        tdc_batch_last_seq = static_cast<uint32_t>(seq);
                    }
                    // Time-based flush (covers sparse streams).
                    if (tdc_batch_hits > 0 &&
                        std::chrono::steady_clock::now() - tdc_batch_last_flush >= TDC_BATCH_MS)
                    {
                        tdc_flush();
                    }

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

                        // Snapshot raw event data so /api/hist_config can
                        // recompute clusters under a new window without
                        // waiting for the next live event.
                        auto ev_copy  = std::make_shared<fdec::EventData>(event);
                        auto ssp_copy = std::make_shared<ssp::SspEventData>(ssp_evt);

                        {
                            std::lock_guard<std::mutex> lk(ring_mtx_);
                            ring_.push_back({seq, std::move(evjson),
                                             std::move(cljson),
                                             std::move(ev_copy),
                                             std::move(ssp_copy)});
                            while ((int)ring_.size() > ring_size_)
                                ring_.pop_front();
                        }

                        wsBroadcast("{\"type\":\"new_event\",\"seq\":" +
                                    std::to_string(seq) + "}");
                    }
                }
            }

            // Drain any partial TDC batch before losing the connection context.
            tdc_flush();

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

// LIVETIME — temporary: poll external command every 30s while ET is active
void ViewerServer::livetimePollThread()
{
    while (running_) {
        while (running_ && !et_active_)
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        if (!running_) break;

        double val = -1.0;
        FILE *p = popen(livetime_cmd_.c_str(), "r");
        if (p) {
            char buf[512];
            std::string out;
            while (fgets(buf, sizeof(buf), p)) out += buf;
            pclose(p);
            // parse "Livetime = XX.X%" from output
            auto pos = out.find("Livetime");
            if (pos != std::string::npos) {
                auto eq = out.find('=', pos);
                if (eq != std::string::npos) {
                    try { val = std::stod(out.substr(eq + 1)); }
                    catch (...) {}
                }
            }
        }
        livetime_.store(val);

        // sleep livetime_poll_sec_, but wake early if shutting down or ET deactivated
        for (int i = 0; i < livetime_poll_sec_ * 10 && running_ && et_active_; ++i)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

#endif // WITH_ET
