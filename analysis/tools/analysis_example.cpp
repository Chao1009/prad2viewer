//=============================================================================
// analysis_example.cpp some examples for offline physics analysis
//
// Usage: analysis_example <input_recon.root> [-o output.root] [-n max_events]
//
// Reads the reconstructed root files, call help functions from PhysicsTools, fills per-module energy histograms
// and moller event analysis histograms, and saves to output ROOT file.
//=============================================================================

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "Replay.h"

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

//hardcoded beam energy and position for yield histograms, can be made configurable if needed
float hycal_z = 5646.f; // distance from target to HyCal front face in mm
float Ebeam = 1100.f; // MeV
//Todo: get run ID from filename or config
int run_id = 12345;

const int kMaxCl = 100;
const int kMaxGEMHits = 400;

struct EventVars_Recon {
    uint32_t event_num = 0;
    uint32_t trigger_bits = 0;
    Long64_t timestamp = 0;
    int n_clusters = 0;
    float cl_x[kMaxCl];
    float cl_y[kMaxCl];
    float cl_energy[kMaxCl];
    int cl_nblocks[kMaxCl];
    int cl_center[kMaxCl];
    // GEM part
    int n_gem_hits = 0;
    uint8_t det_id[kMaxGEMHits];
    float gem_x[kMaxGEMHits];
    float gem_y[kMaxGEMHits];
    float gem_x_charge[kMaxGEMHits];
    float gem_y_charge[kMaxGEMHits];
    float gem_x_peak[kMaxGEMHits];
    float gem_y_peak[kMaxGEMHits];
    int gem_x_size[kMaxGEMHits];
    int gem_y_size[kMaxGEMHits];
};

void setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    tree->SetBranchAddress("event_num",    &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",    &ev.timestamp);
    tree->SetBranchAddress("n_clusters",   &ev.n_clusters);
    tree->SetBranchAddress("cl_x",         ev.cl_x);
    tree->SetBranchAddress("cl_y",         ev.cl_y);
    tree->SetBranchAddress("cl_energy",    ev.cl_energy);
    tree->SetBranchAddress("cl_nblocks",   ev.cl_nblocks);
    tree->SetBranchAddress("cl_center",    ev.cl_center);
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
}

//analysis histograms and variables
TH2F *hit_pos = new TH2F("hit_pos", "Hit positions;X (mm);Y (mm)", 250, -500, 500, 250, -500, 500);
TH1F *one_cluster_energy = new TH1F("one_cluster_energy", "Energy of single-cluster events;Energy (MeV);Counts", 1000, 0, 4000);
TH1F *two_cluster_energy = new TH1F("two_cluster_energy", "Energy of 2-cluster events;Energy (MeV);Counts", 1000, 0, 4000);
TH1F *clusters_energy = new TH1F("clusters_energy", "Energy of all clusters;Energy (MeV);Counts", 1000, 0, 4000);
TH1F *total_energy = new TH1F("total_energy", "Total energy per event;Energy (MeV);Counts", 1000, 0, 4000);
analysis::PhysicsTools::MollerEvent MollerPair;
analysis::PhysicsTools::MollerData hycal_mollers;

int main(int argc, char *argv[])
{
    std::string input_file, output;
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }
    if (optind < argc) input_file = argv[optind];

    if (input_file.empty()) {
        std::cerr << "Usage: replay_hycalRecon <input.evio> [-o out.root] [-n max_events]\n";
        return 1;
    }

    if (output.empty()) {
        output = input_file;
        auto pos = output.find(".root");
        if (pos != std::string::npos) output = output.substr(0, pos);
        output += "_result.root";
    }

    // --- setup ---
    fdec::HyCalSystem hycal;
    hycal.Init(std::string(DATABASE_DIR) + "/hycal_modules.json",
               std::string(DATABASE_DIR) + "/daq_map.json");
    analysis::PhysicsTools physics(hycal);

    //setup input ROOT file and tree
    TFile *infile = TFile::Open(input_file.c_str(), "READ");
    if (!infile || !infile->IsOpen()) {
        std::cerr << "Cannot open " << input_file << "\n";
        return 1;
    }
    TTree *tree = (TTree *)infile->Get("recon");
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in " << input_file << "\n";
        return 1;
    }

    EventVars_Recon ev;
    setupReconBranches(tree, ev);

    //setup output ROOT file
    TFile outfile(output.c_str(), "RECREATE");

    int Nentries = tree->GetEntries();
    Nentries = (max_events > 0) ? std::min(Nentries, max_events) : Nentries;
    
    for(int i = 0; i < Nentries; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "Reading " << i << " events / " << Nentries << " total events\r" << std::flush;
        
        //loop over all the clusters on HyCal
        int nHits = ev.n_clusters;
        float sum_energy = 0.f;
        for (int j = 0; j < nHits; j++) {
            float theta = std::atan(std::sqrt(ev.cl_x[j] * ev.cl_x[j] + ev.cl_y[j] * ev.cl_y[j]) /hycal_z) * 180.f / 3.14159265f;
            //fill histograms for energy vs moduleID
            physics.FillEnergyVsModule(ev.cl_center[j], ev.cl_energy[j]);
            //fill histograms for clusters energy vs theta
            physics.FillEnergyVsTheta(theta, ev.cl_energy[j]);
            //fill histograms for clusters hit positions and energy distribution
            hit_pos->Fill(ev.cl_x[j], ev.cl_y[j]);
            clusters_energy->Fill(ev.cl_energy[j]);

            sum_energy += ev.cl_energy[j];
        }
        total_energy->Fill(sum_energy);

        // select events with only 1 cluster on HyCal(mostly elastic ep events)
        if (nHits == 1){
            physics.FillModuleEnergy(ev.cl_center[0], ev.cl_energy[0]);
            one_cluster_energy->Fill(ev.cl_energy[0]);
        }

        //select events with 2 clusters on HyCal (potential Moller events)
        // no GEM matching for simple quick start
        if (nHits == 2){
            two_cluster_energy->Fill(ev.cl_energy[0]);
            two_cluster_energy->Fill(ev.cl_energy[1]);
            //try to find good Moller with energy cut
            if(std::abs(ev.cl_energy[0] + ev.cl_energy[1] - Ebeam) < 3.*Ebeam*0.025/sqrt(Ebeam/1000.f)){
                //save these 2-cluster events for Moller analysis
                MollerPair = analysis::PhysicsTools::MollerEvent(
                                {ev.cl_x[0], ev.cl_y[0], 0.f, ev.cl_energy[0]},
                                {ev.cl_x[1], ev.cl_y[1], 0.f, ev.cl_energy[1]});
                hycal_mollers.push_back(MollerPair);
                //fill Moller phi difference histogram
                float phi_diff = physics.GetMollerPhiDiff(MollerPair);
                physics.FillMollerPhiDiff(phi_diff);
                //fill 2-arm Moller position histogram
                physics.Fill2armMollerPosHist(MollerPair.first.x, MollerPair.first.y);
                physics.Fill2armMollerPosHist(MollerPair.second.x, MollerPair.second.y);
            }
        }
    }

    //analyze Moller events, get ditector position information and z-vertex distribution, fill histograms
    for (int i = 0; i < hycal_mollers.size(); i++) {
        float z = physics.GetMollerZdistance(hycal_mollers[i], Ebeam);
        physics.FillMollerZ(z);
        if(i >= 1) {
            auto center = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            physics.FillMollerXY(center[0], center[1]);
        } 
    }

    // write histograms into output ROOT file
    outfile.cd();
    hit_pos->Write();

    outfile.mkdir("energy_plots");
    outfile.cd("energy_plots");
    if (physics.GetEnergyVsModuleHist()) physics.GetEnergyVsModuleHist()->Write();
    if (physics.GetEnergyVsThetaHist())  physics.GetEnergyVsThetaHist()->Write();
    one_cluster_energy->Write();
    two_cluster_energy->Write();
    clusters_energy->Write();
    total_energy->Write();

    outfile.cd();
    outfile.mkdir("physics_yields");
    outfile.cd("physics_yields");
    auto h_ep = physics.GetEpYieldHist(physics.GetEnergyVsThetaHist(), Ebeam);
    auto h_ee = physics.GetEeYieldHist(physics.GetEnergyVsThetaHist(), Ebeam);
    auto h_ratio = physics.GetYieldRatioHist(h_ep.get(), h_ee.get());
    if (h_ep) h_ep->Write();
    if (h_ee) h_ee->Write();
    if (h_ratio) h_ratio->Write();

    outfile.cd();
    outfile.mkdir("moller_analysis");
    outfile.cd("moller_analysis");
    if (physics.Get2armMollerPosHist()) physics.Get2armMollerPosHist()->Write();
    else std::cerr << "No 2-arm Moller position histogram filled.\n";
    if (physics.GetMollerPhiDiffHist()) physics.GetMollerPhiDiffHist()->Write();
    else std::cerr << "No Moller phi difference histogram filled.\n";
    
    if(physics.GetMollerXHist()) physics.GetMollerXHist()->Write();
    if(physics.GetMollerYHist()) physics.GetMollerYHist()->Write();
    if(physics.GetMollerZHist()) physics.GetMollerZHist()->Write();

    outfile.mkdir("module_energy");
    outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); i++) {
        TH1F *h = physics.GetModuleEnergyHist(i);
        if (h && h->GetEntries() > 0) h->Write();
    }

    outfile.Close();
    physics.Resolution2Database(run_id); // example run ID

    std::cerr << "\r" << Nentries << " events analyzed, the result saved -> " << output << "\n";

    return 0;
}
