#pragma once
// =========================================================================
// data_source.h — Abstract data source interface for the event viewer
//
// Allows the viewer to read events from EVIO files, ROOT raw replay files,
// or ROOT recon files through a uniform interface.
// =========================================================================

#include "Fadc250Data.h"
#include "SspData.h"
#include "DscData.h"
#include "DaqConfig.h"   // evc::EventType (small POD enum, cheap to include)

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fdec { class HyCalSystem; }

// ── Capabilities ─────────────────────────────────────────────────────────

struct DataSourceCaps {
    bool has_waveforms  = false;   // raw ADC samples per channel
    bool has_peaks      = false;   // per-channel peak info
    bool has_pedestals  = false;   // per-channel pedestal mean/rms
    bool has_clusters   = false;   // cluster data (computed or pre-computed)
    bool has_gem_raw    = false;   // GEM raw strip data
    bool has_gem_hits   = false;   // GEM reconstructed hits
    bool has_epics      = false;   // EPICS slow control events
    bool has_sync       = false;   // sync/control events (absolute time)
    std::string source_type;       // "evio", "root_raw", "root_recon"
};

// ── Pre-computed reconstruction data (ROOT recon files) ──────────────────

struct ReconCluster {
    float x, y, energy;
    int nblocks, center_id;
};

struct ReconGemHit {
    int det_id;
    float x, y;
    float x_charge, y_charge;
    float x_peak, y_peak;
    int x_size, y_size;
};

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;
    uint32_t trigger_bits = 0;
    uint64_t timestamp    = 0;
    std::vector<ReconCluster> clusters;
    std::vector<ReconGemHit>  gem_hits;
};

// ── DataSource interface ─────────────────────────────────────────────────

class DataSource {
public:
    virtual ~DataSource() = default;

    // Open a file. Returns empty string on success, error message on failure.
    virtual std::string open(const std::string &path) = 0;
    virtual void close() = 0;

    // Capabilities of this data source.
    virtual DataSourceCaps capabilities() const = 0;

    // Total number of physics events (available after open).
    virtual int eventCount() const = 0;

    // Decode event by 0-based index into EventData.
    // For recon sources, fills only EventInfo fields (nrocs=0).
    // Returns empty string on success, error message on failure.
    virtual std::string decodeEvent(int index, fdec::EventData &evt,
                                     ssp::SspEventData *ssp = nullptr) = 0;

    // Classified event type for the given 0-based index (Physics / Sync /
    // Epics / control / Unknown).  Used by the viewer to label non-Physics
    // samples in the status bar — those events are kept in the index so the
    // EPICS/control bookkeeping sees them, but they decode to empty FADC.
    // Default returns Physics for sources that don't track the distinction
    // (e.g. ROOT recon files where every entry is a real readout).
    virtual evc::EventType eventTypeAt(int /*index*/) const
    {
        return evc::EventType::Physics;
    }

    // Decode pre-computed cluster/GEM data (recon sources only).
    // Returns false if not supported or index out of range.
    virtual bool decodeReconEvent(int index, ReconEventData &recon) { return false; }

    // Iterate all events for histogram/LMS accumulation.
    // EVIO/ROOT raw sources call ev_cb for each physics event.
    // ROOT recon sources call recon_cb instead.
    // EVIO sources also call ctrl_cb (sync/control) and epics_cb.
    using EventCallback   = std::function<void(int idx, fdec::EventData &evt,
                                                ssp::SspEventData *ssp)>;
    using ReconCallback   = std::function<void(int idx, const ReconEventData &recon)>;
    using ControlCallback = std::function<void(uint32_t unix_time, uint64_t last_ti_ts)>;
    using EpicsCallback   = std::function<void(const std::string &text,
                                                int32_t ev_num, uint64_t timestamp)>;
    // Decoded DSC2 scaler record.  Fired on every Sync/Physics event whose
    // configured DSC2 bank decodes to a record with `present == true`
    // (Sync events typically; some sites also embed the bank in physics
    // events).  EVIO-only — ROOT sources ignore the parameter.
    using DscCallback     = std::function<void(const dsc::DscEventData &dsc)>;

    virtual void iterateAll(EventCallback ev_cb,
                            ReconCallback recon_cb = nullptr,
                            ControlCallback ctrl_cb = nullptr,
                            EpicsCallback epics_cb = nullptr,
                            DscCallback dsc_cb = nullptr,
                            int dsc_bank_tag = -1) = 0;
};

// ── Factory ──────────────────────────────────────────────────────────────

// Create the appropriate DataSource for a file path.
// Auto-detects by extension (.evio → EVIO, .root → ROOT) and tree name.
// crate_to_roc maps crate IDs (0,1,...) to ROC tags (0x80,0x82,...) for
// ROOT files where the replay stores crate IDs.
// hycal is required for ROOT raw files to reverse module_id back to DAQ
// (crate, slot, channel) addressing. Pass nullptr for EVIO-only usage.
// Returns nullptr if the file type is unrecognized or support not compiled in.
std::unique_ptr<DataSource> createDataSource(
    const std::string &path,
    const evc::DaqConfig &daq_cfg,
    const std::unordered_map<int, uint32_t> &crate_to_roc,
    const fdec::HyCalSystem *hycal = nullptr);
