// =========================================================================
// evio_data_source.cpp — EVIO file data source implementation
// =========================================================================

#include "evio_data_source.h"

#include <iostream>
#include <memory>

using namespace evc;

// =========================================================================
// Open / Close
// =========================================================================

std::string EvioDataSource::open(const std::string &path)
{
    close();
    filepath_ = path;

    // scan file to build event index
    EvChannel ch;
    ch.SetConfig(cfg_);
    if (ch.Open(path) != status::success)
        return "cannot open file";

    int buf = 0;
    while (ch.Read() == status::success) {
        ++buf;
        if (!ch.Scan()) continue;
        for (int i = 0; i < ch.GetNEvents(); ++i)
            index_.push_back({buf, i});
    }
    ch.Close();
    return "";
}

void EvioDataSource::close()
{
    index_.clear();
    invalidateReader();
    filepath_.clear();
}

// =========================================================================
// Capabilities
// =========================================================================

DataSourceCaps EvioDataSource::capabilities() const
{
    return {
        true,   // has_waveforms
        true,   // has_peaks (computed by WaveAnalyzer)
        true,   // has_pedestals
        true,   // has_clusters (computed by HyCalCluster)
        true,   // has_gem_raw
        true,   // has_gem_hits (computed by GemCluster)
        true,   // has_epics
        true,   // has_sync
        "evio"  // source_type
    };
}

// =========================================================================
// Random-access event decoding
// =========================================================================

std::string EvioDataSource::seekTo(int buf_num)
{
    if (filepath_ != reader_path_ || buf_num < reader_buf_) {
        reader_.Close();
        reader_.SetConfig(cfg_);
        if (reader_.Open(filepath_) != status::success) {
            invalidateReader();
            return "cannot open file";
        }
        reader_path_ = filepath_;
        reader_buf_ = 0;
    }
    while (reader_buf_ < buf_num) {
        if (reader_.Read() != status::success) {
            invalidateReader();
            return "read error";
        }
        reader_buf_++;
    }
    return "";
}

void EvioDataSource::invalidateReader()
{
    reader_.Close();
    reader_path_.clear();
    reader_buf_ = 0;
}

std::string EvioDataSource::decodeEvent(int index, fdec::EventData &evt,
                                         ssp::SspEventData *ssp)
{
    if (index < 0 || index >= (int)index_.size())
        return "event out of range";

    auto &ei = index_[index];
    std::lock_guard<std::mutex> lk(reader_mtx_);
    std::string err = seekTo(ei.buffer_num);
    if (!err.empty()) return err;
    if (!reader_.Scan()) return "scan error";
    if (ssp) ssp->clear();
    // DecodeEvent populates evt.info (event#, timestamp, trigger_type) regardless
    // of return value. Returns false for events without detector data (e.g. monitoring
    // triggers) — this is normal, not an error.
    reader_.DecodeEvent(ei.sub_event, evt, ssp);
    return "";
}

// =========================================================================
// Full iteration (for histogram building)
// =========================================================================

void EvioDataSource::iterateAll(EventCallback ev_cb, ReconCallback /*recon_cb*/,
                                ControlCallback ctrl_cb, EpicsCallback epics_cb)
{
    EvChannel ch;
    ch.SetConfig(cfg_);
    if (ch.Open(filepath_) != status::success) return;

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    auto &ssp_evt = *ssp_ptr;
    uint64_t last_ti_ts = 0;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;

        // control events (sync/prestart/go/end)
        if (ctrl_cb) {
            uint32_t ct = ch.GetControlTime();
            if (ct != 0) ctrl_cb(ct, last_ti_ts);
        }

        // EPICS events
        if (epics_cb && ch.GetEventType() == EventType::Epics) {
            std::string text = ch.ExtractEpicsText();
            if (!text.empty())
                epics_cb(text, 0, last_ti_ts);
        }

        // physics events
        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ssp_evt.clear();
            if (!ch.DecodeEvent(i, event, &ssp_evt)) continue;
            last_ti_ts = event.info.timestamp;
            if (ev_cb) ev_cb(i, event, &ssp_evt);
        }
    }
    ch.Close();
}
