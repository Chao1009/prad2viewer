// epCalib.cpp: tool to get calibration constants from elastic e-pevents
// on each module, read the rawdata.root(peak mode) file, fit the elastic peak,
// get the ratio of expected/measured peak position, and write to a database file. 
//=============================================================================
//
// Usage: epCalib <input.evio> [-o output_calib_file] [-D daq_config.json] [-n max_events]
//
// Reads rawdata.root (peak mode), runs HyCal clustering, fills per-module energy histograms
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "WaveAnalyzer.h"

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

float hycal_z = 6225.f; // distance from target to HyCal front face in mm, used for kinematics

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

void SetReadBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->SetBranchAddress("event_num", &ev.event_num);
    tree->SetBranchAddress("trigger",   &ev.trigger);
    tree->SetBranchAddress("timestamp", &ev.timestamp);
    tree->SetBranchAddress("nch",       &ev.nch);
    tree->SetBranchAddress("crate",     ev.crate);
    tree->SetBranchAddress("slot",      ev.slot);
    tree->SetBranchAddress("channel",   ev.channel);
    tree->SetBranchAddress("module_id", ev.module_id);
    tree->SetBranchAddress("nsamples",  ev.nsamples);
    tree->SetBranchAddress("ped_mean",  ev.ped_mean);
    tree->SetBranchAddress("ped_rms",   ev.ped_rms);
    tree->SetBranchAddress("integral",  ev.integral);
    if (write_peaks) {
        tree->SetBranchAddress("npeaks",       ev.npeaks);
        tree->SetBranchAddress("peak_height",  ev.peak_height);
        tree->SetBranchAddress("peak_time",    ev.peak_time);
        tree->SetBranchAddress("peak_integral",ev.peak_integral);
    }
}

int main(int argc, char *argv[])
{
    std::string input, output_calib_file, config_file, daq_config_file;
    std::string db_dir = DATABASE_DIR;
    int max_events = -1;

    //hardcoded beam energy for yield histograms, can be made configurable if needed
    float Ebeam = 3500.f; // MeV

    int opt;
    while ((opt = getopt(argc, argv, "o:D:n:")) != -1) {
        switch (opt) {
            case 'o': output_calib_file = optarg; break;
            case 'D': daq_config_file = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty()) {
        std::cerr << "Usage: epCalib <input.evio> [-o output_calib_file] "
                  << " [-D daq_config.json] [-n max_events]\n";
        return 1;
    }

    if (output_calib_file.empty()) {
        output_calib_file = db_dir + "/fast_ep_calibration/calib.txt";
    }

    TFile *infile = TFile::Open(input.c_str(), "READ");
    if (!infile || !infile->IsOpen()) {
        std::cerr << "Cannot open " << input << "\n";
        return 1;
    }
    TTree *tree = (TTree *)infile->Get("events");
    if (!tree) {
        std::cerr << "Cannot find TTree 'events' in " << input << "\n";
        return 1;
    }

    TFile outfile("ep_calib.root", "RECREATE");
    if (!outfile.IsOpen()) {
        std::cerr << "Cannot create ep_calib.root\n";
        return 1;
    }

    auto ev = std::make_unique<EventVars>();
    SetReadBranches(tree, *ev, true);

    //setup for reconstruction
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    if (!daq_config_file.empty()) evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json");
    analysis::PhysicsTools physics(hycal);
    fdec::ClusterConfig cl_cfg;

    //loop over events, fill histograms
    int nentries = tree->GetEntries();
    nentries = (max_events > 0) ? std::min(nentries, max_events) : nentries;
    for(int i = 0; i < nentries; i++){
        tree->GetEntry(i);
        //reconstruct clusters, fill histograms
        fdec::HyCalCluster clusterer(hycal);
        clusterer.SetConfig(cl_cfg);
        for(int j = 0; j < ev->nch; j++){
            const auto *mod = hycal.module_by_daq(ev->crate[j], ev->slot[j], ev->channel[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev->npeaks[j] <= 0) continue;
            float adc = ev->peak_integral[j][0];
            float energy = (mod->cal_factor > 0.) ?
                static_cast<float>(mod->energize(adc)) : adc * 0.078f;
            clusterer.AddHit(mod->index, energy);
        }

        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        if(hits.size() == 1) physics.FillModuleEnergy(hits[0].center_id, hits[0].energy);
    }

    //after loop, fit peaks, get calibration constants, write to database    
    int nmod = hycal.module_count();
    //ratio of expected/measured peak position for each module
    TH1F *ratio_module_all = new TH1F("ratio_all", "Ratio of Expected/Measured Peak Position for All Modules;Ratio;Modules", 100, 0, 2);
    std::vector<float> ratio_values(nmod, 0.f);
    for (int m = 0; m < nmod; m++) {
        auto [peak, resolution] = physics.FitPeakResolution(m);
        if (peak > 0 && resolution > 0) {
            std::string name = hycal.module(m).name;
            float theta_deg = atan(sqrt(pow(hycal.module(m).x, 2) + pow(hycal.module(m).y, 2)) / hycal_z) * 180.f / 3.14159265f;
            float expected_peak = physics.ExpectedEnergy(theta_deg, Ebeam, "ep");
            float ratio = expected_peak / peak;
            ratio_module_all->Fill(ratio);
            ratio_values[m] = ratio;

            double current_factor = hycal.GetCalibConstant(hycal.module(m).id);
            double new_factor = current_factor * ratio;
            hycal.SetCalibConstant(hycal.module(m).id, new_factor);
        }
    }
    //write the new calibration constants to database file
    hycal.PrintCalibConstants(output_calib_file);
    outfile.mkdir("module_energy");
    outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); ++i) {
        TH1F *h = physics.GetModuleEnergyHist(i);
        if (h && h->GetEntries() > 0) h->Write();
    }
    outfile.cd();
    ratio_module_all->Write();
    outfile.Close();
    infile->Close();

    return 0;
}