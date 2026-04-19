// epCalib.cpp: tool to get calibration constants from elastic e-pevents
// on each module, read the rawdata.root(peak mode) file, fit the elastic peak,
// get the ratio of expected/measured peak position, and write to a database file. 
//=============================================================================
//
// Usage: epCalib <input_raw.root|dir> [-i iteration] [-o output_root_file]
//                                     [-E Ebeam] [-D daq_config.json] [-n max_events]
//   - input_raw.root|dir: input ROOT file or directory containing ROOT files with raw data (peak mode) to be analyzed. The tool will look for TTree named "events"
//   - iteration: iteration number for calibration, used for bookkeeping and output file naming (default: 0)
//   - output_root_file: output ROOT file to write calibration constants (default: <db_dir>/calibration/<run_number>/ep_calib.root)
//   - Ebeam: beam energy in MeV, used for calculating expected elastic peak position (default: 2100 MeV)
//   - daq_config.json: DAQ config file for mapping ROC tags to crate indices, needed for decoding raw data (default: none, but recommended to provide)
//   - max_events: maximum number of events to process (default: -1 for all events)
//
// Reads rawdata(adc level).root (peak mode), runs HyCal clustering, fills per-module energy histograms
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "InstallPaths.h"
#include "load_daq_config.h"

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
void SetReadBranches(TTree *tree, EventVars &ev)
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
    tree->SetBranchAddress("hycal.npeaks", &ev.npeaks);
    tree->SetBranchAddress("hycal.peak_height", ev.peak_height);
    tree->SetBranchAddress("hycal.peak_time", ev.peak_time);
    tree->SetBranchAddress("hycal.peak_integral", ev.peak_integral);
}

static std::vector<std::string> collectRootFiles(const std::string &path);

int main(int argc, char *argv[])
{
    std::string input_calib_file, output_calib_file;
    std::string output_root_file, daq_config_file;
    int iteration = 1; // default to iteration 1, which uses cosmic calibration as input. User can specify -i 0 to start from scratch (no input calib file, just use adc*0.1 as energy)
    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    int max_events = -1;

    //hardcoded beam energy for yield histograms, can be made configurable if needed
    float Ebeam = 2100.f; // MeV
    float hycal_z = 6225.f; // distance from target to HyCal front face in mm, used for kinematics

    int opt;
    while ((opt = getopt(argc, argv, "i:o:E:D:n:")) != -1) {
        switch (opt) {
            case 'i': iteration = std::atoi(optarg); break;
            case 'o': output_root_file = optarg; break;
            case 'E': Ebeam = std::atof(optarg); break;
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
        std::cerr << "Usage: epCalib <input_raw.root|dir> [more files...] -i iteration -o output_root_file -E Ebeam -D daq_config.json -n max_events]\n";
        return 1;
    }

    // extract run number from first input file name (e.g. prad_023626.00000_raw.root -> 23626)
    std::string run_str = "unknown";
    {
        std::string fname = fs::path(root_files[0]).filename().string();
        auto ppos = fname.find("prad_");
        if (ppos != std::string::npos) {
            size_t s = ppos + 5;
            size_t e = s;
            while (e < fname.size() && std::isdigit((unsigned char)fname[e])) e++;
            if (e > s) run_str = std::to_string(std::stoul(fname.substr(s, e - s)));
        }
    }
    // make output directory: db_dir/calibration/<run_number>/
    std::string run_out_dir = db_dir + "/calibration/" + run_str;
    fs::create_directories(run_out_dir);
    std::cerr << "Output directory: " << run_out_dir << "\n";
    if( iteration == 1 )
        input_calib_file = db_dir + "/calibration/adc_to_mev_factors_cosmic.json";
    else if( iteration > 1 )
        input_calib_file = run_out_dir + Form("/calib_iter%d.json",   iteration-1);
    else{
        std::cerr << "Invalid iteration number: " << iteration << ". Must be >= 1.\n";
        return 1;
    }
    output_calib_file = run_out_dir + Form("/calib_iter%d.json", iteration);
    

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
    if (output_root_file.empty()) {
        output_root_file = run_out_dir + Form("/CalibResult_iter%d.root", iteration);
    }
    TString outName = output_root_file;
    TFile outfile(outName, "RECREATE");

    auto ev = std::make_unique<EventVars>();
    SetReadBranches(tree, *ev);

    //setup for reconstruction
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    if (daq_config_file.empty()) daq_config_file = db_dir + "/daq_config.json"; // default DAQ config for PRad2
    evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json");

    std::string calib_file = input_calib_file;
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    analysis::PhysicsTools physics(hycal);
    fdec::ClusterConfig cl_cfg;

    // --- histograms ---
    TH2F *hit_pos = new TH2F("hit_pos",
        "One Cluster Event Hit positions;hycal X (mm);hycal Y (mm)", 250, -500, 500, 250, -500, 500);
    TH1F *h_E_1cl = new TH1F("one_cluster_energy",
        "Single-cluster Event energy;E (MeV);Counts", 1000, 0, 5000);
    TH1F *h_mearured_peak = new TH1F("measured_peak", 
        "Measured Peak Position;Energy (MeV);Counts", 1000, 0, 5000);
    TH1F *h_recon_sigma = new TH1F("recon_sigma", 
        "Reconstructed Cluster Energy Resolution;Sigma (MeV);Counts", 100, 0, 200);
    TH1F *h_recon_chi2 = new TH1F("recon_chi2/ndf", 
        "Reconstructed E hist Fit Chi2/ndf;Chi2/ndf;Counts", 100, 0, 50);
    //ratio of expected/measured peak position for each module
    TH1F *ratio_module_all = new TH1F("ratio_all", 
        "Ratio of Expected/Measured Peak Position for All Modules;Ratio;Modules", 200, 0, 4);
    TH2F *module_ratio = new TH2F("#cbar#bar{E_{recon}} - E_{expect}#cbar #/ E_{expect}",
        "#cbar#bar{E_{recon}} - E_{expect}#cbar #/ E_{expect}",
        34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);
    //loop over events, fill histograms
    int nentries = tree->GetEntries();
    nentries = (max_events > 0) ? std::min(nentries, max_events) : nentries;
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cl_cfg);
    for(int i = 0; i < nentries; i++){
        tree->GetEntry(i);
        if (i % 1000 == 0) 
            std::cout << "Processing event " << i+1 << "/" << nentries << "\r" << std::flush;

        static constexpr uint32_t TBIT_sum = (1u << 8);
        static constexpr uint32_t TBIT_lms  = (1u << 24);
        if (!(ev->trigger_bits & TBIT_sum)) continue;
        if (  ev->trigger_bits & TBIT_lms  ) continue;

        //reconstruct clusters, fill histograms
        clusterer.Clear();
        for(int j = 0; j < ev->nch; j++){
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev->npeaks[j] <= 0) continue;
            float adc = ev->peak_integral[j][0];
            float energy = (mod->cal_factor > 0) ?
                static_cast<float>(mod->energize(adc)) : adc*0.1f;
            clusterer.AddHit(mod->index, energy);
        }

        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        if(hits.size() == 1) {   
            physics.FillModuleEnergy(hits[0].center_id, hits[0].energy);
            hit_pos->Fill(hits[0].x, hits[0].y);
            h_E_1cl->Fill(hits[0].energy);
            float theta = atan(sqrt(pow(hits[0].x, 2) + pow(hits[0].y, 2)) / hycal_z) * 180.f / 3.14159265f;
            physics.FillEnergyVsTheta(theta, hits[0].energy);
        }
    }
    std::cerr << "\nFinished processing " << nentries << " events.\n";
    std::cerr << "Fitting peaks and calculating calibration constants...\n";

    //after loop, fit peaks, get calibration constants, write to database  
    std::string dat_out_path = run_out_dir + Form("/fitting_parameters_iter%d.dat", iteration);
    std::ofstream dat_out(dat_out_path);
    if (!dat_out.is_open()) {
        std::cerr << "Cannot open output file " << dat_out_path << "\n";
        return 1;
    }  
    int nmod = hycal.module_count();
    std::vector<float> ratio_values(nmod, 0.f);

    TCanvas *c = new TCanvas("c", "Calibration", 1200, 1200);
    c->cd();
    dat_out << std::left;
    dat_out << std::setw(8) << "Module" << std::setw(16) << "ExpectedPeak" << std::setw(16) <<
    "MeasuredPeak" << std::setw(16) << "Ratio" << std::setw(16) << "Sigma" << std::setw(16) << "Chi2/ndf" << "\n";
    int n_calibrated = 0;
    for (int m = 0; m < nmod; m++) {
        int mod_id = hycal.module(m).id;
        auto [peak, sigma, chi2] = physics.FitPeakResolution(mod_id);
        if (peak > 0 && sigma > 0 && chi2 < 100.f) {
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

            h_mearured_peak->Fill(peak);
            h_recon_sigma->Fill(sigma);
            h_recon_chi2->Fill(chi2);

            dat_out << std::setw(8) << name << std::setw(16) << expected_peak << std::setw(16) << peak << std::setw(16) << ratio 
            << std::setw(16) << sigma << std::setw(16) << chi2 << "\n";
            n_calibrated++;
        }
    }
    std::cerr << "Calibrated " << n_calibrated << " modules. Results written to " << dat_out_path << "\n";

    module_ratio->SetStats(0);
    module_ratio->Draw("COLZ");
    module_ratio->Write();

    TLatex t;
    t.SetTextSize(0.0122);
    t.SetTextColor(kBlack);
    for(int m = 0; m < nmod; m++) {
        std::string name = hycal.module(m).name;
        if(name[0] != 'W') continue;
        t.DrawLatex(hycal.module(m).x-6., hycal.module(m).y-2., hycal.module(m).name.c_str());
    }
    
    c->Update();
    c->Write();

    TCanvas *c_hit = new TCanvas("c_hit", "Hit Position", 1200, 1200);
    c_hit->cd();
    hit_pos->Draw("COLZ");
    hit_pos->Write();
    outfile.cd();
    c_hit->Write();

    physics.GetEnergyVsThetaHist()->Write();
    h_E_1cl->Write();
    h_mearured_peak->Write();
    h_recon_sigma->Write();
    h_recon_chi2->Write();
    ratio_module_all->Write();

    physics.FillNeventsModuleMap();
    TH2F *h_map = physics.GetNeventsModuleMapHist();
    TCanvas *c_map = new TCanvas("c_map", "Number of Events per Module", 1200, 1200);
    h_map->Draw("COLZ");
    c_map->Write();
    h_map->Write();

    outfile.mkdir("module_energy");
    outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); ++i) {
        int mod_id = hycal.module(i).id;
        TH1F *h = physics.GetModuleEnergyHist(mod_id);
        if (h && h->GetEntries() > 0) h->Write();
    }

    //write the new calibration constants to database file
    hycal.PrintCalibConstants(output_calib_file);

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