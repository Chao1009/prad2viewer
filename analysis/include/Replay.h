#pragma once
//=============================================================================
// Replay.h — convert raw DAQ data (EVIO) to ROOT trees
//
// Decodes EVIO events and writes per-channel waveform/peak data to a TTree.
// Depends on prad2dec (decoder) and ROOT (TFile/TTree).
//=============================================================================

#include "EvChannel.h"
#include "EventData.h"
#include "WaveAnalyzer.h"
#include "DaqConfig.h"
#include "ConfigSetup.h"
#include "load_daq_config.h"

#include <TFile.h>
#include <TTree.h>

#include <string>
#include <unordered_map>

namespace analysis {

// Aliases for the shared replay data structures
using EventVars       = prad2::RawEventData;
using EventVars_Recon = prad2::ReconEventData;

class Replay
{
public:
    Replay() = default;

    // Load DAQ configuration (event tags, ADC format, etc.).
    void LoadDaqConfig(const std::string &json_path) { evc::load_daq_config(json_path, daq_cfg_); }

    // Load DAQ map (module name lookup by crate/slot/channel).
    void LoadDaqMap(const std::string &json_path);

    std::string moduleName(int roc, int slot, int ch) const;
    int moduleID(int roc, int slot, int ch) const;

    // Convert an EVIO file to a ROOT file with a TTree.
    // max_events <= 0 means process all. peaks=true adds peak branches.
    bool Process(const std::string &input_evio, const std::string &output_root, const std::string &db_dir,
                 int max_events = -1, bool write_peaks = false, const std::string &daq_config_file = "");

    bool ProcessWithRecon(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                            const std::string &db_dir,
                            const std::string &daq_config_file = "",
                            const std::string &gem_ped_file = "", float zerosup_override = 0.f,
                            bool prad1 = false);

private:
    void setupBranches(TTree *tree, EventVars &ev, bool write_peaks);
    void clearEvent(EventVars &ev);

    void setupReconBranches(TTree *tree, EventVars_Recon &ev);
    void clearReconEvent(EventVars_Recon &ev);

    float computeIntegral(const fdec::ChannelData &cd, float pedestal) const;

    using DaqMap = std::unordered_map<std::string, std::string>;  // "roc_slot_ch" -> name
    DaqMap daq_map_;
    evc::DaqConfig daq_cfg_;
};

} // namespace analysis
