#pragma once
// =========================================================================
// evio_data_source.h — EVIO file data source for the event viewer
// =========================================================================

#include "data_source.h"
#include "EvChannel.h"
#include "DaqConfig.h"

#include <mutex>
#include <string>
#include <vector>

class EvioDataSource : public DataSource {
public:
    explicit EvioDataSource(const evc::DaqConfig &cfg) : cfg_(cfg) {}

    std::string open(const std::string &path) override;
    void close() override;
    DataSourceCaps capabilities() const override;
    int eventCount() const override { return (int)index_.size(); }

    std::string decodeEvent(int index, fdec::EventData &evt,
                             ssp::SspEventData *ssp = nullptr) override;

    evc::EventType eventTypeAt(int index) const override
    {
        if (index < 0 || index >= (int)index_.size())
            return evc::EventType::Unknown;
        return index_[index].type;
    }

    void iterateAll(EventCallback ev_cb, ReconCallback recon_cb,
                    ControlCallback ctrl_cb, EpicsCallback epics_cb,
                    DscCallback dsc_cb, int dsc_bank_tag) override;

private:
    evc::DaqConfig cfg_;
    std::string filepath_;

    // 0-based evio event index (for EvChannel::ReadEventByIndex), the
    // sub-event index within that event's built-trigger block, and the
    // classified event type so we can label non-Physics samples in the
    // viewer without re-scanning.
    struct EvioIndex { int evio_event, sub_event; evc::EventType type; };
    std::vector<EvioIndex> index_;

    // Persistent reader.  Opened during open() via EvChannel::OpenAuto so
    // we get random-access mode when the file supports it, sequential
    // otherwise.  decodeEvent() dispatches on reader_.IsRandomAccess():
    //   * RA mode  → reader_.ReadEventByIndex(ei.evio_event) — O(1).
    //   * seq mode → close/reopen on backward jumps, walk forward to target.
    evc::EvChannel reader_;
    std::string reader_path_;
    int reader_pos_ = -1;   // sequential-mode read cursor (-1 = pre-first)
    // Index of the event currently decoded in reader_'s lazy cache, or -1 if
    // the cache is invalid.  When decodeEvent() is called with the same index
    // a second time (common for viewer_server's decodeEvent + computeClusters
    // pair on the same click), we skip Scan+decode and copy straight out of
    // reader_.Fadc()/Gem().
    int last_decoded_index_ = -1;
    std::mutex reader_mtx_;

    // Walk the sequential reader forward (or close/reopen on backward jumps)
    // until it's positioned on evio_target.  Returns an error string on
    // failure, empty on success.
    std::string seekTo(int evio_target);
    void invalidateReader();
};
