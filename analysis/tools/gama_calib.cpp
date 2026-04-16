// epCalib.cpp: tool to get calibration constants from photon beam in calibration runs
// on each module, read the rawdata.root(peak mode) file, reconstruction and fit the peak,
// get the ratio of expected/measured peak position, and write to a database file. 
//=============================================================================
//
// Usage: gama_calib <input_raw.root|dir> input_calib_file output_calib_file
//                         --option       [-o output_root_file] [-n max_events]
//
// Reads rawdata(adc level).root (peak mode), runs HyCal clustering, fills per-module energy histograms
//=============================================================================

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "EventData.h"
#include "load_daq_config.h"

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars       = prad2::RawEventData;
void SetReadBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->SetBranchAddress("event_num", &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp", &ev.timestamp);
    tree->SetBranchAddress("hycal.nch", &ev.nch);
    tree->SetBranchAddress("hycal.module_id", ev.module_id);
    tree->SetBranchAddress("hycal.nsamples", ev.nsamples);
    tree->SetBranchAddress("hycal.samples", ev.samples);
    tree->SetBranchAddress("hycal.ped_mean", ev.ped_mean);
    tree->SetBranchAddress("hycal.ped_rms", ev.ped_rms);
    tree->SetBranchAddress("hycal.integral", ev.integral);
    if (write_peaks) {
        tree->SetBranchAddress("hycal.npeaks", &ev.npeaks);
        tree->SetBranchAddress("hycal.peak_height", ev.peak_height);
        tree->SetBranchAddress("hycal.peak_time", ev.peak_time);
        tree->SetBranchAddress("hycal.peak_integral", ev.peak_integral);
    }
}

static std::vector<std::string> collectRootFiles(const std::string &path);

int main(int argc, char *argv[])
{
    std::string in_calib_file, out_calib_file, output_root_file;
    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    int max_events = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:")) != -1) {
        switch (opt) {
            case 'o': output_root_file = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }

    if(optind + 2 != argc) {
        std::cerr << "Usage: " << argv[0] << " <input_raw.root|dir> <input_calib_file> <output_calib_file>\n"
                  << "Options:\n"
                  << "  -o output_root_file   Output ROOT file for reconstructed variables (optional)\n"
                  << "  -n max_events         Maximum number of events to process (default: all)\n";
        return 1;
    }
    std::string input = argv[optind];
    in_calib_file = argv[optind + 1];
    out_calib_file = argv[optind + 2];

    // collect input files (can be files, directories)
    std::vector<std::string> input_files = collectRootFiles(input);
    if (input_files.empty()) {
        std::cerr << "No input files found in " << input << "\n";
        return 1;
    }

    // --- setup TChain and branches ---
    TChain *chain = new TChain("events");
    for (const auto &f : input_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }
    auto ev = std::make_unique<EventVars>();
    SetReadBranches(tree, *ev, true);

    //setup for reconstruction
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    std::string daq_config_file = db_dir + "/daq_config.json"; // default DAQ config for PRad2
    evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json");

    int nmatched = hycal.LoadCalibration(in_calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << in_calib_file << " (" << nmatched << " modules)\n";

    analysis::PhysicsTools physics(hycal);
    fdec::ClusterConfig cl_cfg;

    //loop over events
    
}


// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_raw.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}
