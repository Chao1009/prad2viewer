// epCalib.cpp: tool to get calibration constants from elastic e-pevents
// on each module, read the rawdata.root(peak mode) file, fit the elastic peak,
// get the ratio of expected/measured peak position, and write to a database file. 
//=============================================================================
//
// Usage: epCalib <input_raw.root|dir> [-o output_calib_file] [-r output_root_file] 
//                                      [-D daq_config.json] [-n max_events]
//
// Reads rawdata(adc level).root (peak mode), runs HyCal clustering, fills per-module energy histograms
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "WaveAnalyzer.h"
#include "EventData.h"

#include <TFile.h>
#include <TTree.h>
#include <TLatex.h>
#include <TCanvas.h>
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
    std::string output_calib_file, output_root_file, config_file, daq_config_file;
    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    int max_events = -1;

    //hardcoded beam energy for yield histograms, can be made configurable if needed
    float Ebeam = 3500.f; // MeV
    float hycal_z = 6225.f; // distance from target to HyCal front face in mm, used for kinematics

    int opt;
    while ((opt = getopt(argc, argv, "o:r:D:n:")) != -1) {
        switch (opt) {
            case 'o': output_calib_file = optarg; break;
            case 'r': output_root_file = optarg; break;
            case 'D': daq_config_file = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }

    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_raw.root|dir> [more files...] [-o out_calib.txt] [-r out_root.root] [-n max_events]\n";
        return 1;
    }

    if (output_calib_file.empty()) {
        output_calib_file = db_dir + "/fast_ep_calibration/calib.txt";
    }

    // --- setup TChain and branches ---
    TChain *chain = new TChain("events");
    for (const auto &f : root_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }

    // --- output file ---
    TString outName = output_root_file;
    if (outName.IsNull()) {
        outName = root_files[0];
        outName.ReplaceAll("_recon.root", "_epCalibResult.root");
    }
    TFile outfile(outName, "RECREATE");

    auto ev = std::make_unique<EventVars>();
    SetReadBranches(tree, *ev, true);

    //setup for reconstruction
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    if (!daq_config_file.empty()) evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json");

    std::string calib_file = db_dir + "/prad1/prad_calibration.json";
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    analysis::PhysicsTools physics(hycal);
    fdec::ClusterConfig cl_cfg;

    //loop over events, fill histograms
    int nentries = tree->GetEntries();
    nentries = (max_events > 0) ? std::min(nentries, max_events) : nentries;
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cl_cfg);
    for(int i = 0; i < nentries; i++){
        tree->GetEntry(i);
        //reconstruct clusters, fill histograms
        clusterer.Clear();
        for(int j = 0; j < ev->nch; j++){
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev->npeaks[j] <= 0) continue;
            float adc = ev->peak_integral[j][0];
            float energy = (mod->cal_factor > 0.2) ?
                static_cast<float>(mod->energize(adc)) : adc;
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
    TH2F *module_ratio = new TH2F("#cbar#bar{E_{recon}} - E_{expect}#cbar #/ E_{expect}",
                                  "#cbar#bar{E_{recon}} - E_{expect}#cbar #/ E_{expect}",
                                  34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);
    std::vector<float> ratio_values(nmod, 0.f);
    TLatex t;
    t.SetTextSize(0.01);
    t.SetTextColor(kBlack);
    TCanvas *c = new TCanvas("c", "Calibration", 1200, 1200);
    c->cd();
    for (int m = 0; m < nmod; m++) {
        auto [peak, resolution] = physics.FitPeakResolution(m);
        if (peak > 0 && resolution > 0) {
            std::string name = hycal.module(m).name;
            if(name[0] != 'W') continue; 
            float theta_deg = atan(sqrt(pow(hycal.module(m).x, 2) + pow(hycal.module(m).y, 2)) / hycal_z) * 180.f / 3.14159265f;
            float expected_peak = physics.ExpectedEnergy(theta_deg, Ebeam, "ep");
            float ratio = expected_peak / peak;
            ratio_module_all->Fill(ratio);
            ratio_values[m] = ratio;

            double current_factor = hycal.GetCalibConstant(hycal.module(m).id);
            double new_factor = current_factor * ratio;
            hycal.SetCalibConstant(hycal.module(m).id, new_factor);

            module_ratio->Fill(hycal.module(m).x, hycal.module(m).y, abs(1.f-1.f/ratio));
            t.DrawLatex(hycal.module(m).x, hycal.module(m).y, name.c_str());
        }
    }
    module_ratio->SetStats(0);
    module_ratio->Draw("COLZ");
    module_ratio->Write();
    c->Write();
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

    return 0;
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