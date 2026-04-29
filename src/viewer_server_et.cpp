#include "viewer_server.h"
#include "http_compress.h"

#ifdef WITH_ET
#include "EtChannel.h"
#endif

#include <cstdio>
#include <cstring>

using namespace evc;
using json = nlohmann::json;

#ifdef WITH_ET

// =========================================================================
// Tagger live-stream frame format
//
// Header (little-endian, 24 bytes, matches the dtype expected by
// scripts/tagger_viewer.py):
//
//   char     magic[4];       // "TGR1"
//   uint32_t flags;           // bit 0 set when dropped_count > 0
//   uint32_t n_hits;          // number of 16-byte BinHit records that follow
//   uint32_t first_seq;       // event seq of the first hit in the frame
//   uint32_t last_seq;        // event seq of the last hit
//   uint32_t dropped;         // total dropped frames (since server start)
//
// Hits use a 16-byte packed layout (same as scripts/tagger_viewer.py RAW_DTYPE):
//
//   uint32_t event_num;
//   uint32_t trigger_bits;
//   uint16_t roc_tag;
//   uint8_t  slot;
//   uint8_t  channel_edge;    // bit 7 = edge, bits 6:0 = channel
//   uint32_t tdc;             // raw V1190 TDC value (hardware-level name)
// =========================================================================
namespace {

constexpr size_t TAGGER_HIT_SIZE    = 16;
constexpr size_t TAGGER_HDR_SIZE    = 24;
constexpr uint32_t TAGGER_BATCH_MAX = 256;    // flush at this many hits
constexpr auto     TAGGER_BATCH_MS  = std::chrono::milliseconds(10);

#pragma pack(push, 1)
struct TaggerBinHit {
    uint32_t event_num;
    uint32_t trigger_bits;
    uint16_t roc_tag;
    uint8_t  slot;
    uint8_t  channel_edge;
    uint32_t tdc;
};
#pragma pack(pop)
static_assert(sizeof(TaggerBinHit) == TAGGER_HIT_SIZE, "TaggerBinHit layout");

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

    // Tagger batch state (local to the thread — only this thread writes it).
    // Pre-allocate enough for TAGGER_BATCH_MAX hits + header.
    std::vector<uint8_t> tagger_batch;
    tagger_batch.reserve(TAGGER_HDR_SIZE + TAGGER_BATCH_MAX * TAGGER_HIT_SIZE);
    uint32_t tagger_batch_hits = 0;
    uint32_t tagger_batch_first_seq = 0;
    uint32_t tagger_batch_last_seq = 0;
    auto tagger_batch_last_flush = std::chrono::steady_clock::now();

    auto tagger_flush = [&]() {
        if (tagger_batch_hits == 0) { tagger_batch_last_flush = std::chrono::steady_clock::now(); return; }
        // Fill header in place.
        uint8_t *p = tagger_batch.data();
        std::memcpy(p + 0,  "TGR1", 4);
        uint32_t drops = static_cast<uint32_t>(tagger_dropped_frames_.load());
        uint32_t flags = (drops > 0) ? 1u : 0u;
        std::memcpy(p + 4,  &flags,                  4);
        std::memcpy(p + 8,  &tagger_batch_hits,      4);
        std::memcpy(p + 12, &tagger_batch_first_seq, 4);
        std::memcpy(p + 16, &tagger_batch_last_seq,  4);
        std::memcpy(p + 20, &drops,                  4);
        taggerBroadcastBinary(tagger_batch.data(),
                              TAGGER_HDR_SIZE + tagger_batch_hits * TAGGER_HIT_SIZE);
        // Reset batch (keep allocation).
        tagger_batch.resize(TAGGER_HDR_SIZE);
        tagger_batch_hits = 0;
        tagger_batch_first_seq = 0;
        tagger_batch_last_seq = 0;
        tagger_batch_last_flush = std::chrono::steady_clock::now();
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

                // DSC2 scaler bank → measured livetime (Sync events typically;
                // some sites embed it in physics events too, so check both).
                if (app_online_.daq_cfg.dsc_scaler.enabled()) {
                    auto et = ch.GetEventType();
                    if (et == EventType::Sync || et == EventType::Physics) {
                        const auto *node = ch.FindFirstByTag(
                            (uint32_t)app_online_.daq_cfg.dsc_scaler.bank_tag);
                        if (node && node->data_words > 0)
                            app_online_.processDscBank(ch.GetData(*node), node->data_words);
                    }
                }

                for (int i = 0; i < ch.GetNEvents(); ++i) {
                    ssp_evt.clear();
                    const bool want_tagger =
                        tagger_subs_count_.load(std::memory_order_relaxed) > 0;
                    tdc::TdcEventData *tagger_arg = want_tagger ? &tdc_evt : nullptr;
                    if (tagger_arg) tagger_arg->clear();
                    if (!ch.DecodeEvent(i, event, &ssp_evt, nullptr, tagger_arg)) continue;
                    last_ti_ts = event.info.timestamp;

                    app_online_.processGemEvent(ssp_evt);
                    app_online_.processEvent(event, ana, wres);

                    int seq = app_online_.events_processed.load();

                    // --- live tagger stream: batch + broadcast ---------------
                    if (want_tagger && tdc_evt.n_hits > 0) {
                        if (tagger_batch_hits == 0) {
                            // Grow to header size on first hit of the batch.
                            tagger_batch.resize(TAGGER_HDR_SIZE);
                            tagger_batch_first_seq = static_cast<uint32_t>(seq);
                        }
                        const uint32_t evnum = static_cast<uint32_t>(event.info.event_number);
                        const uint32_t tbits = event.info.trigger_bits;
                        for (int h = 0; h < tdc_evt.n_hits; ++h) {
                            const auto &src = tdc_evt.hits[h];
                            TaggerBinHit bh;
                            bh.event_num    = evnum;
                            bh.trigger_bits = tbits;
                            bh.roc_tag      = static_cast<uint16_t>(src.roc_tag);
                            bh.slot         = src.slot;
                            bh.channel_edge =
                                static_cast<uint8_t>(((src.edge & 0x1) << 7) |
                                                      (src.channel & 0x7F));
                            bh.tdc          = src.value;
                            size_t off = tagger_batch.size();
                            tagger_batch.resize(off + TAGGER_HIT_SIZE);
                            std::memcpy(tagger_batch.data() + off, &bh, TAGGER_HIT_SIZE);
                            ++tagger_batch_hits;
                            if (tagger_batch_hits >= TAGGER_BATCH_MAX) { tagger_flush(); break; }
                        }
                        tagger_batch_last_seq = static_cast<uint32_t>(seq);
                    }
                    // Time-based flush (covers sparse streams).
                    if (tagger_batch_hits > 0 &&
                        std::chrono::steady_clock::now() - tagger_batch_last_flush >= TAGGER_BATCH_MS)
                    {
                        tagger_flush();
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
                        // GEM per-APV waveforms — encoded here so the API
                        // for older ring events doesn't need to re-process
                        // gem_sys (which would clobber the live state used
                        // by /api/gem/hits etc.).  gem_sys was just filled
                        // by processGemEvent above for this event.
                        std::string gemapvjson = app_online_.gem_enabled
                            ? app_online_.apiGemApv(ssp_evt, seq).dump()
                            : std::string("{\"enabled\":false}");
                        // Pre-compress the gem_apv payload once here so
                        // every viewer's HTTP fetch serves the cached
                        // bytes (vs deflating the same ~1.3 MB JSON for
                        // each viewer × 5 Hz refresh).  Skip below the
                        // gzip threshold — the disabled stub is tiny.
                        std::string gemapvgz;
                        if (gemapvjson.size() >= prad2::kGzipMinBytes) {
                            try {
                                gemapvgz = prad2::gzip_compress(gemapvjson);
                            } catch (...) {
                                gemapvgz.clear();   // serve plain on failure
                            }
                        }

                        // Snapshot raw event data so /api/hist_config can
                        // recompute clusters under a new window without
                        // waiting for the next live event.
                        auto ev_copy  = std::make_shared<fdec::EventData>(event);
                        auto ssp_copy = std::make_shared<ssp::SspEventData>(ssp_evt);

                        {
                            std::lock_guard<std::mutex> lk(ring_mtx_);
                            ring_.push_back({seq, std::move(evjson),
                                             std::move(cljson),
                                             std::move(gemapvjson),
                                             std::move(gemapvgz),
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

            // Drain any partial tagger batch before losing the connection context.
            tagger_flush();

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

// Optional DAQ-livetime poller: shells out to the user-configured command
// (typical: "caget -t <epics_channel>") and parses the first floating-point
// number from stdout.  Empty command disables the thread entirely; this is
// gated in run()/startAsync() so we never spawn a useless poller.
//
// Polls only while ET is active.  On bad output the value goes back to <0
// (frontend hides the display) — the reading is consumed as a snapshot, not
// integrated, so transient parse failures self-heal on the next poll.
void ViewerServer::livetimePollThread()
{
    const std::string cmd = app_file_.livetime_cmd;
    const int poll_sec = std::max(1, app_file_.livetime_poll_sec);
    if (cmd.empty()) return;

    int consecutive_failures = 0;
    while (running_) {
        while (running_ && !et_active_)
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        if (!running_) break;

        double val = -1.0;
        FILE *p = popen(cmd.c_str(), "r");
        if (p) {
            char buf[512];
            std::string out;
            while (fgets(buf, sizeof(buf), p)) out += buf;
            pclose(p);

            // Extract the first floating-point number found in stdout.
            // Works with `caget -t CHAN` (just "<num>"), bare `caget CHAN`
            // ("<chan>  <num>"), or any tool that prints "Livetime = X%".
            size_t i = 0;
            while (i < out.size()) {
                char c = out[i];
                if (c == '-' || c == '+' || c == '.' || (c >= '0' && c <= '9')) {
                    try {
                        size_t consumed = 0;
                        double v = std::stod(out.substr(i), &consumed);
                        if (consumed > 0) { val = v; break; }
                    } catch (...) {}
                }
                ++i;
            }
        }
        livetime_.store(val);
        if (val < 0) {
            if (consecutive_failures++ == 0)
                std::cerr << "Livetime  : poll command produced no number "
                          << "(`" << cmd << "`) — leaving as 'not available'\n";
        } else {
            consecutive_failures = 0;
        }

        for (int i = 0; i < poll_sec * 10 && running_ && et_active_; ++i)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

#endif // WITH_ET
