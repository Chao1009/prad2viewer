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
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TCanvas.h>

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

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;
// ── Tree branch struct ───────────────────────────────────────────────────
void setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    tree->SetBranchAddress("event_num",    &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",    &ev.timestamp);
    tree->SetBranchAddress("total_energy", &ev.total_energy);
    // HyCal cluster branches
    tree->SetBranchAddress("n_clusters",   &ev.n_clusters);
    tree->SetBranchAddress("cl_x",         ev.cl_x);
    tree->SetBranchAddress("cl_y",         ev.cl_y);
    tree->SetBranchAddress("cl_z",         ev.cl_z);
    tree->SetBranchAddress("cl_energy",    ev.cl_energy);
    tree->SetBranchAddress("cl_nblocks",   ev.cl_nblocks);
    tree->SetBranchAddress("cl_center",    ev.cl_center);
    tree->SetBranchAddress("cl_flag",      ev.cl_flag);
    // Matching results
    tree->SetBranchAddress("matchFlag",    ev.matchFlag);
    tree->SetBranchAddress("matchGEMx",    ev.matchGEMx);
    tree->SetBranchAddress("matchGEMy",    ev.matchGEMy);
    tree->SetBranchAddress("matchGEMz",    ev.matchGEMz);
    tree->SetBranchAddress("matchNum",     &ev.matchNum);
    //quick and simple matching results for quick check
    tree->SetBranchAddress("mHit_E", ev.mHit_E);
    tree->SetBranchAddress("mHit_x", ev.mHit_x);
    tree->SetBranchAddress("mHit_y", ev.mHit_y);
    tree->SetBranchAddress("mHit_z", ev.mHit_z);
    tree->SetBranchAddress("mHit_gx", ev.mHit_gx);
    tree->SetBranchAddress("mHit_gy", ev.mHit_gy);
    tree->SetBranchAddress("mHit_gz", ev.mHit_gz);
    tree->SetBranchAddress("mHit_gid", ev.mHit_gid);
    // GEM part
    tree->SetBranchAddress("n_gem_hits",   &ev.n_gem_hits);
    tree->SetBranchAddress("det_id",       ev.det_id);
    tree->SetBranchAddress("gem_x",        ev.gem_x);
    tree->SetBranchAddress("gem_y",        ev.gem_y);
    tree->SetBranchAddress("gem_x_charge", ev.gem_x_charge);
    tree->SetBranchAddress("gem_y_charge", ev.gem_y_charge);
    tree->SetBranchAddress("gem_x_peak",   ev.gem_x_peak);
    tree->SetBranchAddress("gem_y_peak",   ev.gem_y_peak);
    tree->SetBranchAddress("gem_x_size",   ev.gem_x_size);
    tree->SetBranchAddress("gem_y_size",   ev.gem_y_size);
};

static std::vector<std::string> collectRootFiles(const std::string &path);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 1100.f;
    int run_id = 12345;

    // --- geometry constants (can be made configurable) ---
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
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

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

    MollerData hycal_mollers;

    // --- event loop : basic distributions + Moller candidates on HyCal ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;

    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "\rPass 1: " << i << " / " << N << std::flush;

        for (int j = 0; j < ev.n_clusters; j++) {
            float r = std::sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]);
            float theta = std::atan(r / ev.cl_z[j]) * 180.f / M_PI;

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
                MollerEvent mp(
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

    physics.FillNeventsModuleMap();
    TH2F *h_map = physics.GetNeventsModuleMapHist();
    TCanvas *c_map = new TCanvas("c_map", "Number of Events per Module", 1200, 1200);
    h_map->Draw("COLZ");
    c_map->Write();
    h_map->Write();

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