// quick_check.C — ROOT script version 
//
// Reads reconstructed ROOT tree (output of replay_recon), runs physics
// analysis using PhysicsTools and MatchingTools from prad2det, and saves
// histograms to an output ROOT file.
// Usage:
//   quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]
//   -o  output ROOT file (default: input filename with _quick_check.root suffix)
//   -n  max events to process (default: all)
// Example:
//   quick_check recon.root -o recon_check.root -n 10000
//   quick_check recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// ── Tree branch struct ───────────────────────────────────────────────────

static const int kMaxCl = 100;
static const int kMaxGEMHits = 400;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;
void setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/i");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    tree->Branch("total_energy", &ev.total_energy, "total_energy/F");
    // HyCal cluster branches
    tree->Branch("n_clusters",   &ev.n_clusters,   "n_clusters/I");
    tree->Branch("cl_x",         ev.cl_x,          "cl_x[n_clusters]/F");
    tree->Branch("cl_y",         ev.cl_y,          "cl_y[n_clusters]/F");
    tree->Branch("cl_z",         ev.cl_z,          "cl_z[n_clusters]/F");
    tree->Branch("cl_energy",    ev.cl_energy,     "cl_energy[n_clusters]/F");
    tree->Branch("cl_nblocks",   ev.cl_nblocks,    "cl_nblocks[n_clusters]/b");
    tree->Branch("cl_center",    ev.cl_center,     "cl_center[n_clusters]/s");
    tree->Branch("cl_flag",      ev.cl_flag,       "cl_flag[n_clusters]/i");
    // GEM part
    tree->Branch("n_gem_hits",   &ev.n_gem_hits,   "n_gem_hits/I");
    tree->Branch("det_id",       ev.det_id,        "det_id[n_gem_hits]/b");
    tree->Branch("gem_x",        ev.gem_x,         "gem_x[n_gem_hits]/F");
    tree->Branch("gem_y",        ev.gem_y,         "gem_y[n_gem_hits]/F");
    tree->Branch("gem_x_charge", ev.gem_x_charge,  "gem_x_charge[n_gem_hits]/F");
    tree->Branch("gem_y_charge", ev.gem_y_charge,  "gem_y_charge[n_gem_hits]/F");
    tree->Branch("gem_x_peak",   ev.gem_x_peak,    "gem_x_peak[n_gem_hits]/F");
    tree->Branch("gem_y_peak",   ev.gem_y_peak,    "gem_y_peak[n_gem_hits]/F");
    tree->Branch("gem_x_size",   ev.gem_x_size,    "gem_x_size[n_gem_hits]/b");
    tree->Branch("gem_y_size",   ev.gem_y_size,    "gem_y_size[n_gem_hits]/b");
    // Matching results
    tree->Branch("match_num",       &ev.match_num,       "match_num/I");
    tree->Branch("matchHC_x",       ev.matchHC_x,        "matchHC_x[match_num]/F");
    tree->Branch("matchHC_y",       ev.matchHC_y,        "matchHC_y[match_num]/F");
    tree->Branch("matchHC_z",       ev.matchHC_z,        "matchHC_z[match_num]/F");
    tree->Branch("matchHC_energy",  ev.matchHC_energy,   "matchHC_energy[match_num]/F");
    tree->Branch("matchHC_center",  ev.matchHC_center,   "matchHC_center[match_num]/s");
    tree->Branch("matchHC_flag",    ev.matchHC_flag,     "matchHC_flag[match_num]/i");
    tree->Branch("matchG_x",        ev.matchG_x,         Form("matchG_x[match_num][2]/F"));
    tree->Branch("matchG_y",        ev.matchG_y,         Form("matchG_y[match_num][2]/F"));
    tree->Branch("matchG_z",        ev.matchG_z,         Form("matchG_z[match_num][2]/F"));
    tree->Branch("matchG_det_id",   ev.matchG_det_id,    Form("matchG_det_id[match_num][2]/b"));

    tree->Branch("gem.n_ssp_triggers", &ev.n_ssp_triggers, "n_ssp_triggers/b");
    tree->Branch("gem.ssp_trigger_tags", ev.ssp_trigger_tags, Form("ssp_trigger_tags[n_ssp_triggers][%d]/i", ssp::SSP_TIME_SAMPLES));
};

static std::vector<std::string> collectRootFiles(const std::string &path);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 1100.f;
    int run_id = 12345;

    // --- geometry constants (can be made configurable) ---
    const float hycal_z = 6225.f; //5646 mm for prad1, 6225mm for prad2
    const float gem_z[4] = {5407.f + 39.71f/2, 5407.f - 39.71f/2,
                            5807.f + 39.71f/2, 5807.f - 39.71f/2};
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
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
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = DATABASE_DIR;

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_modules.json",
               dbDir + "/daq_map.json");
    PhysicsTools physics(hycal);
    MatchingTools matching;
    
    // --- setup TChain and branches ---
    TChain *chain = new TChain("recon");
    for (const auto &f : root_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev;
    setupReconBranches(tree, ev);

    // --- output file ---
    TString outName = output;
    if (outName.IsNull()) {
        outName = root_files[0];
        outName.ReplaceAll("_recon.root", "_quick_check.root");
    }
    TFile outfile(outName, "RECREATE");

    // --- histograms ---
    TH2F *hit_pos = new TH2F("hit_pos",
        "Hit positions;X (mm);Y (mm)", 250, -500, 500, 250, -500, 500);
    TH1F *h_1cl = new TH1F("one_cluster_energy",
        "Single-cluster energy;E (MeV);Counts", 1000, 0, 4000);
    TH1F *h_2cl = new TH1F("two_cluster_energy",
        "Two-cluster energy;E (MeV);Counts", 1000, 0, 4000);
    TH1F *h_all = new TH1F("clusters_energy",
        "All clusters;E (MeV);Counts", 1000, 0, 4000);
    TH1F *h_tot = new TH1F("total_energy",
        "Total energy per event;E (MeV);Counts", 1000, 0, 4000);

    PhysicsTools::MollerData hycal_mollers;

    // --- event loop : basic distributions + Moller candidates on HyCal ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;

    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "\rPass 1: " << i << " / " << N << std::flush;

        for (int j = 0; j < ev.n_clusters; j++) {
            float r = std::sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]);
            float theta = std::atan(r / hycal_z) * 180.f / M_PI;

            physics.FillEnergyVsModule(ev.cl_center[j], ev.cl_energy[j]);
            physics.FillEnergyVsTheta(theta, ev.cl_energy[j]);
            hit_pos->Fill(ev.cl_x[j], ev.cl_y[j]);
            h_all->Fill(ev.cl_energy[j]);
        }
        h_tot->Fill(ev.total_energy);

        if (ev.n_clusters == 1) {
            physics.FillModuleEnergy(ev.cl_center[0], ev.cl_energy[0]);
            h_1cl->Fill(ev.cl_energy[0]);
        }

        if (ev.n_clusters == 2) {
            h_2cl->Fill(ev.cl_energy[0]);
            h_2cl->Fill(ev.cl_energy[1]);

            float Epair = ev.cl_energy[0] + ev.cl_energy[1];
            float sigma = Ebeam * 0.025f / std::sqrt(Ebeam / 1000.f);
            if (std::abs(Epair - Ebeam) < 3. * sigma) {
                PhysicsTools::MollerEvent mp(
                    {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                    {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
                hycal_mollers.push_back(mp);
                physics.FillMollerPhiDiff(physics.GetMollerPhiDiff(mp));
                physics.Fill2armMollerPosHist(mp.first.x, mp.first.y);
                physics.Fill2armMollerPosHist(mp.second.x, mp.second.y);
            }
        }
    }

    // --- Moller vertex analysis ---
    for (size_t i = 0; i < hycal_mollers.size(); i++) {
        physics.FillMollerZ(physics.GetMollerZdistance(hycal_mollers[i], Ebeam));
        if (i >= 1) {
            auto c = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            physics.FillMollerXY(c[0], c[1]);
        }
    }

    // --- write output ---
    outfile.cd();
    hit_pos->Write();

    outfile.mkdir("energy_plots"); outfile.cd("energy_plots");
    if (physics.GetEnergyVsModuleHist()) physics.GetEnergyVsModuleHist()->Write();
    if (physics.GetEnergyVsThetaHist())  physics.GetEnergyVsThetaHist()->Write();
    h_1cl->Write(); h_2cl->Write(); h_all->Write(); h_tot->Write();

    outfile.cd();
    outfile.mkdir("physics_yields"); outfile.cd("physics_yields");
    auto h_ep = physics.GetEpYieldHist(physics.GetEnergyVsThetaHist(), Ebeam);
    auto h_ee = physics.GetEeYieldHist(physics.GetEnergyVsThetaHist(), Ebeam);
    auto h_ratio = physics.GetYieldRatioHist(h_ep.get(), h_ee.get());
    if (h_ep) h_ep->Write();
    if (h_ee) h_ee->Write();
    if (h_ratio) h_ratio->Write();

    outfile.cd();
    outfile.mkdir("moller_analysis"); outfile.cd("moller_analysis");
    if (physics.Get2armMollerPosHist()) physics.Get2armMollerPosHist()->Write();
    if (physics.GetMollerPhiDiffHist()) physics.GetMollerPhiDiffHist()->Write();
    if (physics.GetMollerXHist()) physics.GetMollerXHist()->Write();
    if (physics.GetMollerYHist()) physics.GetMollerYHist()->Write();
    if (physics.GetMollerZHist()) physics.GetMollerZHist()->Write();

    outfile.mkdir("module_energy"); outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); i++) {
        TH1F *h = physics.GetModuleEnergyHist(i);
        if (h && h->GetEntries() > 0) h->Write();
    }

    outfile.Close();
    
    physics.Resolution2Database(run_id); // example run ID

    std::cerr << "Result saved -> " << outName.Data() << "\n";
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_recon.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}