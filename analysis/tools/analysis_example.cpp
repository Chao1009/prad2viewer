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
#include "MatchingTools.h"
#include "EventData.h"

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

using namespace analysis;

//hardcoded beam energy and position for yield histograms, can be made configurable if needed
float hycal_z = 5646.f; // distance from target to HyCal front face in mm
float gem_z[4] = {5407.+39.71/2., 5407.-39.71/2., 5807.+39.71/2., 5807.-39.71/2.}; //mm, the center of the two GEMs, 39.71 is the gap of 2 GEMs
float Ebeam = 1100.f; // MeV
//Todo: get run ID from filename or config
int run_id = 12345;

const int kMaxCl = 100;
const int kMaxGEMHits = 400;

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

//analysis histograms and variables
TH2F *hit_pos = new TH2F("hit_pos", "Hit positions;X (mm);Y (mm)", 250, -500, 500, 250, -500, 500);
TH1F *one_cluster_energy = new TH1F("one_cluster_energy", "Energy of single-cluster events;Energy (MeV);Counts", 1000, 0, 4000);
TH1F *two_cluster_energy = new TH1F("two_cluster_energy", "Energy of 2-cluster events;Energy (MeV);Counts", 1000, 0, 4000);
TH1F *clusters_energy = new TH1F("clusters_energy", "Energy of all clusters;Energy (MeV);Counts", 1000, 0, 4000);
TH1F *total_energy = new TH1F("total_energy", "Total energy per event;Energy (MeV);Counts", 1000, 0, 4000);
PhysicsTools::MollerEvent MollerPair;
PhysicsTools::MollerData hycal_mollers;

PhysicsTools::MollerData good_Hmollers;
PhysicsTools::MollerData good_Gmollers[4];

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
    PhysicsTools physics(hycal);
    MatchingTools matching;

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
    
    // this event loop do not include GEM matching, very rough selection for quick analysis example, more strict selection can be implemented with GEM matching and other kinematic cuts
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
                MollerPair = PhysicsTools::MollerEvent(
                                {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                                {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
                hycal_mollers.push_back(MollerPair);
                //fill Moller phi difference histogram
                float phi_diff = physics.GetMollerPhiDiff(MollerPair);
                physics.FillMollerPhiDiff(phi_diff);
                //fill 2-arm Moller position histogram
                physics.Fill2armMollerPosHist(MollerPair.first.x, MollerPair.first.y);
                physics.Fill2armMollerPosHist(MollerPair.second.x, MollerPair.second.y);
            }
        }
    }//end of event loop

    //analyze Moller events saved in the loop,
    // get detector position information and z-vertex distribution, fill histograms
    for (int i = 0; i < hycal_mollers.size(); i++) {
        float z = physics.GetMollerZdistance(hycal_mollers[i], Ebeam);
        physics.FillMollerZ(z);
        if(i >= 1) {
            auto center = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            physics.FillMollerXY(center[0], center[1]);
        } 
    }
    std::cerr << "\r" << Nentries << " events analyzed\n";

    //this event loop will include GEM matching and coordinate transformation
    // more strict event selection here for better e-p and e-e
    for(int i = 0; i < Nentries; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "Reading " << i << " events / " << Nentries << " total events\r" << std::flush;

        //store all the hits on HyCal and GEMs in this event
        std::vector<HCHit> hc_hits;
        std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
        for( int j = 0; j < ev.n_clusters; j++) {
            hc_hits.push_back(HCHit{ev.cl_x[j], ev.cl_y[j], ev.cl_z[j],
                               ev.cl_energy[j], ev.cl_center[j], ev.cl_flag[j]});
        }
        for (int j = 0; j < ev.n_gem_hits; j++) {
            gem_hits[ev.det_id[j]].push_back(GEMHit{ev.gem_x[j], ev.gem_y[j], 0.f, ev.det_id[j]});
        }

        //transform detector coordinates to target and beam center coordinates
        TransformDetData(hc_hits, 0.f, 0.f, hycal_z); // assuming beamX=beamY=0 for now
        for(int d = 0; d < 4; d++) 
            TransformDetData(gem_hits[d], 0.f, 0.f, gem_z[d]);

        //then matching between GEM hits and HyCal clusters
        //matching.SetMatchRange(5.0f); // matching radius in mm
        //matching.SetSquareSelection(true); // use square cut instead of circular cut
        std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
        //show how to access the matching result
        for (auto &m : matched_hits) { 
            HCHit hycal_hit = m.hycal_hit;  //the HyCal cluster be matched
            GEMHit gem_up_hit = m.gem[0];    //best-matched GEM hit from upstream pair (GEM1/GEM2)
            GEMHit gem_down_hit = m.gem[1]; //best-matched GEM hit from downstream pair (GEM3/GEM4)
            std::vector<GEMHit> gem1_matches = m.gem1_hits;
            std::vector<GEMHit> gem2_matches = m.gem2_hits;
            std::vector<GEMHit> gem3_matches = m.gem3_hits;
            std::vector<GEMHit> gem4_matches = m.gem4_hits;

            int hycal_idx = m.hycal_idx;  //index of the cluster in the original vector
        }
        
        //select good e-e events
        if (matched_hits.size() == 2){
            auto &m = matched_hits;
            float energy_sum = m[0].hycal_hit.energy + m[1].hycal_hit.energy;
            if(std::abs(energy_sum - Ebeam) < 5.*Ebeam*0.025/sqrt(Ebeam/1000.f)){
                //save these 2-cluster events for Moller analysis
                MollerPair = PhysicsTools::MollerEvent(
                                {m[0].hycal_hit.x, m[0].hycal_hit.y, m[0].hycal_hit.z, m[0].hycal_hit.energy},
                                {m[1].hycal_hit.x, m[1].hycal_hit.y, m[1].hycal_hit.z, m[1].hycal_hit.energy});
                good_Hmollers.push_back(MollerPair);
                
                // check upstream GEM pair (GEM1/GEM2)
                if(m[0].gem[0].det_id == m[1].gem[0].det_id && m[0].gem[0].det_id != 0){
                    MollerPair = PhysicsTools::MollerEvent(
                                    {m[0].gem[0].x, m[0].gem[0].y, m[0].gem[0].z, m[0].hycal_hit.energy},
                                    {m[1].gem[0].x, m[1].gem[0].y, m[1].gem[0].z, m[1].hycal_hit.energy});
                    good_Gmollers[m[0].gem[0].det_id].push_back(MollerPair);
                }
                // check downstream GEM pair (GEM3/GEM4)
                if(m[0].gem[1].det_id == m[1].gem[1].det_id && m[0].gem[1].det_id != 0){
                    MollerPair = PhysicsTools::MollerEvent(
                                    {m[0].gem[1].x, m[0].gem[1].y, m[0].gem[1].z, m[0].hycal_hit.energy},
                                    {m[1].gem[1].x, m[1].gem[1].y, m[1].gem[1].z, m[1].hycal_hit.energy});
                    good_Gmollers[m[0].gem[1].det_id].push_back(MollerPair);
                }
            }
        }

        //select good e-p events
        if (matched_hits.size() == 1){
            auto &m = matched_hits[0];
            float energy = m.hycal_hit.energy;
            if(std::abs(energy - Ebeam) < 3.*Ebeam*0.025/sqrt(Ebeam/1000.f)){
                //save these 1-cluster events for elastic e-p analysis
            }
        }
    }//end of event loop

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

    std::cerr << "The result saved -> " << output << "\n";

    return 0;
}
