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

    void iterateAll(EventCallback ev_cb, ReconCallback recon_cb,
                    ControlCallback ctrl_cb, EpicsCallback epics_cb) override;

private:
    evc::DaqConfig cfg_;
    std::string filepath_;

    struct EvioIndex { int buffer_num, sub_event; };
    std::vector<EvioIndex> index_;

    // cached sequential reader for random access
    evc::EvChannel reader_;
    std::string reader_path_;
    int reader_buf_ = 0;
    // Index of the event currently decoded in reader_'s lazy cache, or -1 if
    // the cache is invalid.  When decodeEvent() is called with the same index
    // a second time (common for viewer_server's decodeEvent + computeClusters
    // pair on the same click), we skip Scan+decode and copy straight out of
    // reader_.Fadc()/Gem().
    int last_decoded_index_ = -1;
    std::mutex reader_mtx_;

    std::string seekTo(int buf_num);
    void invalidateReader();
};
