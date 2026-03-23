#pragma once

#include "EvChannel.h"
#include "EvStruct.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <nlohmann/json.hpp>

#define PROGRESS_COUNT 100
#define FADC_CHANNELS 16

using namespace fdec;
using DaqMap = std::unordered_map<std::string, std::string>;

//variables for root tree
static const int max_chs = MAX_ROCS * MAX_SLOTS * FADC_CHANNELS;
static int event_num, ftrigger;
static Long64_t ftimestamp;
static int Nchs, sample_num[max_chs];
static int crate_id[max_chs], slot_id[max_chs], channel_id[max_chs], module_id[max_chs];
static int samples[max_chs][MAX_SAMPLES];
static float ped[max_chs], ped_err[max_chs], integral[max_chs];
static int Npeaks[max_chs];
static float peak[max_chs][MAX_PEAKS], peak_time[max_chs][MAX_PEAKS], peak_inte[max_chs][MAX_PEAKS];

//

class replay
{
    public:
        replay( ) {}

        bool rawdata2root( const std::string &input_path, const std::string &output_path,
            const int max_events = -1, const bool mode_peaks = false);

        void setupTTree(TTree *tree, bool mode_peaks = false);
        void setupTTreeRead(TTree *tree, bool mode_peaks = false);
        
        void clear();
    private:
        DaqMap daq_map;

        DaqMap loadDaqMap(const std::string &json_path);
        std::string GetModuleName(int roc, int slot, int channel);
        float GetIntegral(const fdec::ChannelData &cd, const float pedestal);
};