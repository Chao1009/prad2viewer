#include "viewer_server.h"

#ifdef WITH_ET
#include "EtChannel.h"
#endif

using namespace evc;
using json = nlohmann::json;

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
