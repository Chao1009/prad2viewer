#pragma once
//=============================================================================
// Replay.h — convert raw DAQ data (EVIO) to ROOT trees
//
// Decodes EVIO events and writes per-channel waveform/peak data to a TTree.
// Depends on prad2dec (decoder) and ROOT (TFile/TTree).
//=============================================================================

#include "EvChannel.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"
#include "DaqConfig.h"
#include "load_daq_config.h"

#include <TFile.h>
#include <TTree.h>

#include <string>
#include <unordered_map>

namespace analysis {

class Replay
{
public:
    Replay() = default;

    // Load DAQ configuration (event tags, ADC format, etc.).
    void LoadDaqConfig(const std::string &json_path) { evc::load_daq_config(json_path, daq_cfg_); }

    // Load DAQ map (module name lookup by crate/slot/channel).
    void LoadDaqMap(const std::string &json_path);

    // Convert an EVIO file to a ROOT file with a TTree.
    // max_events <= 0 means process all. peaks=true adds peak branches.
    bool Process(const std::string &input_evio, const std::string &output_root,
                 int max_events = -1, bool write_peaks = false);

private:
    // per-event data (sized to worst case, reused)
    static constexpr int kMaxCh = fdec::MAX_ROCS * fdec::MAX_SLOTS * 16;

    struct EventVars {
        int     event_num = 0;
        int     trigger = 0;
        Long64_t timestamp = 0;
        int     nch = 0;
        int     crate[kMaxCh] = {};
        int     slot[kMaxCh] = {};
        int     channel[kMaxCh] = {};
        int     module_id[kMaxCh] = {};
        int     nsamples[kMaxCh] = {};
        int     samples[kMaxCh][fdec::MAX_SAMPLES] = {};
        float   ped_mean[kMaxCh] = {};
        float   ped_rms[kMaxCh] = {};
        float   integral[kMaxCh] = {};
        int     npeaks[kMaxCh] = {};
        float   peak_height[kMaxCh][fdec::MAX_PEAKS] = {};
        float   peak_time[kMaxCh][fdec::MAX_PEAKS] = {};
        float   peak_integral[kMaxCh][fdec::MAX_PEAKS] = {};
    };

    void setupBranches(TTree *tree, EventVars &ev, bool write_peaks);
    void clearEvent(EventVars &ev);

    std::string moduleName(int roc, int slot, int ch) const;
    int moduleID(int roc, int slot, int ch) const;
    float computeIntegral(const fdec::ChannelData &cd, float pedestal) const;

    using DaqMap = std::unordered_map<std::string, std::string>;  // "roc_slot_ch" -> name
    DaqMap daq_map_;
    evc::DaqConfig daq_cfg_;
};

} // namespace analysis
