//=============================================================================
// replay_hycalRecon — HyCal reconstruction replay with physics analysis
//
// Usage: replay_hycalRecon <input.evio> [-o output.root] [-c config.json]
//                          [-D daq_config.json] [-f max_files] [-p(read prad1)]
//
// Reads EVIO, runs HyCal clustering, fills per-module energy histograms
// and writes reconstruction results to a ROOT tree.
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "WaveAnalyzer.h"

#include <nlohmann/json.hpp>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <filesystem>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

float hycal_z = 5646.f; // distance from target to HyCal front face in mm

std::vector<std::string> getFilesInDir(const std::string &dir_path)
{
    std::vector<std::string> files;
    for (auto &entry : std::filesystem::directory_iterator(dir_path)) {
        if (entry.is_regular_file()) {
            if (entry.path().filename().string().find(".evio") != std::string::npos)
                files.push_back(entry.path().string());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

int main(int argc, char *argv[])
{
    std::string input_dir, output, config_file, daq_config_file;
    std::string db_dir = DATABASE_DIR;
    int max_files = -1;
    bool prad1 = false; // set to true for PRad1-specific config

    //hardcoded beam energy for yield histograms, can be made configurable if needed
    float Ebeam = 1100.f; // MeV
    //Todo: get run ID from filename or config
    int run_id = 12345;

    int opt;
    while ((opt = getopt(argc, argv, "o:c:f:p")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'c': config_file = optarg; break;
            case 'f': max_files = std::atoi(optarg); break;
            case 'p': prad1 = true; break;
        }
    }
    if (optind < argc) input_dir = argv[optind];

    if (input_dir.empty()) {
        std::cerr << "Usage: replay_hycalRecon <input.evio> [-o out.root] "
                  << "[-c config.json] [-f M] [-p(read prad1)]\n";
        return 1;
    }
    // Get list of EVIO files in the directory
    std::vector<std::string> evio_files = getFilesInDir(input_dir);
    int num_files = evio_files.size();
    if(num_files == 0){
        std::cerr << "No files found in directory: " << input_dir << "\n";
        return 1;
    }
    if (max_files > 0) {
        num_files = std::min(num_files, max_files);
    }

    if (output.empty()) {
        output = evio_files[0];
        auto pos = output.find(".evio");
        if (pos != std::string::npos) output = output.substr(0, pos);
        output += "_recon.root";
    }

    // --- setup ---
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;

    daq_config_file = db_dir + "/daq_config.json"; // default DAQ config for PRad2
    std::string daq_map_file = db_dir + "/daq_map.json"; // default DAQ map for PRad2
    if(prad1 == true){
        daq_config_file = db_dir + "/prad1/prad_daq_config.json";
        daq_map_file = db_dir + "/prad1/prad_daq_map.json";
    }    
    if (!daq_config_file.empty()) evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", daq_map_file);
    
    if(prad1 == true) evc::load_pedestals(db_dir + "/prad1/adc1881m_pedestals.json", daq_cfg);
    
    std::string calib_file = db_dir + "/prad1/prad_calibration.json";
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";
    analysis::PhysicsTools physics(hycal);

    //analysis histograms and variables
    TH2F *hit_pos = new TH2F("hit_pos", "Hit positions;X (mm);Y (mm)", 250, -500, 500, 250, -500, 500);
    TH1F *one_cluster_energy = new TH1F("one_cluster_energy", "Energy of single-cluster events;Energy (MeV);Counts", 1000, 0, 4000);
    TH1F *two_cluster_energy = new TH1F("two_cluster_energy", "Energy of 2-cluster events;Energy (MeV);Counts", 1000, 0, 4000);
    TH1F *clusters_energy = new TH1F("clusters_energy", "Energy of all clusters;Energy (MeV);Counts", 1000, 0, 4000);
    TH1F *total_energy = new TH1F("total_energy", "Total energy per event;Energy (MeV);Counts", 1000, 0, 4000);
    analysis::PhysicsTools::MollerEvent MollerPair;
    analysis::PhysicsTools::MollerData hycal_mollers;

    //setup output ROOT file and tree
    TFile outfile(output.c_str(), "RECREATE");
    TTree *tree = new TTree("recon", "HyCal reconstruction");
    // tree branches
    int ev_num = 0, n_clusters = 0;
    static const int kMaxCl = 100;
    float cl_x[kMaxCl], cl_y[kMaxCl], cl_energy[kMaxCl];
    int cl_nblocks[kMaxCl], cl_center[kMaxCl];
    tree->Branch("event_num",  &ev_num, "event_num/I");
    tree->Branch("n_clusters", &n_clusters, "n_clusters/I");
    tree->Branch("cl_x",       cl_x,       "cl_x[n_clusters]/F");
    tree->Branch("cl_y",       cl_y,       "cl_y[n_clusters]/F");
    tree->Branch("cl_energy",  cl_energy,  "cl_energy[n_clusters]/F");
    tree->Branch("cl_nblocks", cl_nblocks, "cl_nblocks[n_clusters]/I");
    tree->Branch("cl_center",  cl_center,  "cl_center[n_clusters]/I");

    // build ROC tag → crate index mapping from DAQ config JSON
    std::unordered_map<int, int> roc_to_crate;
    if (!daq_config_file.empty()) {
        std::cout << "Loading DAQ config from " << daq_config_file << "\n";
        std::ifstream dcf(daq_config_file);
        if (dcf.is_open()) {
            auto dcj = nlohmann::json::parse(dcf, nullptr, false, true);
            if (dcj.contains("roc_tags") && dcj["roc_tags"].is_array()) {
            for (auto &entry : dcj["roc_tags"]) {
                int tag   = std::stoi(entry.at("tag").get<std::string>(), nullptr, 16);
                int crate = entry.at("crate").get<int>();
                roc_to_crate[tag] = crate;
            }
}
        }
    }

    //initialize tools for event decoder and cluster reconstruction
    fdec::EventData event;
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    fdec::ClusterConfig cl_cfg;
    int total = 0;

    std::cout << "Processing " << num_files << " evio files, output will be saved to " << output << "\n";
    std::cout << "Files: ";
    for (int f = 0; f < num_files; ++f) std::cout << evio_files[f] << "\n";
    std::cout << "\n";

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);
    for (int f = 0; f < num_files; ++f) {
        if (ch.Open(evio_files[f]) != evc::status::success) {
            std::cerr << "Cannot open " << evio_files[f] << "\n";
            return 1;
        }
        while (ch.Read() == evc::status::success) {
            if (!ch.Scan()) continue;
            for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
                if (!ch.DecodeEvent(ie, event)) continue;

                fdec::HyCalCluster clusterer(hycal);
                clusterer.SetConfig(cl_cfg);
                float energy_sum = 0.f;
                // feed hits
                for (int r = 0; r < event.nrocs; ++r) {
                    auto &roc = event.rocs[r];
                    if (!roc.present) continue;
                    auto cit = roc_to_crate.find(roc.tag);
                    if (cit == roc_to_crate.end()) continue;
                    int crate = cit->second;
                    for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                        if (!roc.slots[s].present) continue;
                        for (int c = 0; c < 160; ++c) {
                            if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                            auto &cd = roc.slots[s].channels[c];
                            if (cd.nsamples <= 0) continue;
                            const auto *mod = hycal.module_by_daq(crate, s, c);
                            if (!mod || !mod->is_hycal()) continue;
                            float adc = 0.f;
                            if(prad1 == true) 
                                adc = cd.samples[0] * 0.543; //0.543 for prad1 run1308,correct to 1.1GeV
                            else{
                                ana.Analyze(cd.samples, cd.nsamples, wres);
                                if (wres.npeaks <= 0) continue;
                                adc = wres.peaks[0].integral;
                            }
                            float energy = (mod->cal_factor > 0.) ?
                                static_cast<float>(mod->energize(adc)) : adc * 0.078f;
                            clusterer.AddHit(mod->index, energy);
                            energy_sum += energy;
                        }
                    }
                }
                total_energy->Fill(energy_sum);
                //if (energy_sum < 0.5 * Ebeam) continue; // skip low-energy events
                clusterer.FormClusters();
                std::vector<fdec::ClusterHit> hits;
                clusterer.ReconstructHits(hits);

                //event reconstrued, fill root tree and histograms
                ev_num = event.info.event_number;
                n_clusters = std::min((int)hits.size(), kMaxCl);
                for (int i = 0; i < n_clusters; ++i) {
                    cl_x[i]       = hits[i].x;
                    cl_y[i]       = hits[i].y;
                    cl_energy[i]  = hits[i].energy;
                    cl_nblocks[i] = hits[i].nblocks;
                    cl_center[i]  = hits[i].center_id;
                    
                    physics.FillEnergyVsModule(hits[i].center_id, hits[i].energy);
                    float theta = std::atan(std::sqrt(hits[i].x * hits[i].x + hits[i].y * hits[i].y) /hycal_z) * 180.f / 3.14159265f;
                    physics.FillEnergyVsTheta(theta, hits[i].energy);

                    hit_pos->Fill(hits[i].x, hits[i].y);
                    clusters_energy->Fill(hits[i].energy);
                }
                tree->Fill();
                total++;

                if(n_clusters == 1) {
                    physics.FillModuleEnergy(hits[0].center_id, hits[0].energy);
                    one_cluster_energy->Fill(hits[0].energy);
                }
                if(n_clusters == 2){
                    two_cluster_energy->Fill(hits[0].energy);
                    two_cluster_energy->Fill(hits[1].energy);
                    //try to find good Moller events with 2 clusters
                    if(std::abs(hits[0].energy + hits[1].energy - Ebeam) < 3.*Ebeam*0.025/sqrt(Ebeam/1000.f)){
                        MollerPair = analysis::PhysicsTools::MollerEvent(
                                        {hits[0].x, hits[0].y, 0.f, hits[0].energy},
                                        {hits[1].x, hits[1].y, 0.f, hits[1].energy});
                        hycal_mollers.push_back(MollerPair);
                        float phi_diff = physics.GetMollerPhiDiff(MollerPair);
                        physics.FillMollerPhiDiff(phi_diff);
                        physics.Fill2armMollerPosHist(MollerPair.first.x, MollerPair.first.y);
                        physics.Fill2armMollerPosHist(MollerPair.second.x, MollerPair.second.y);
                    }
                }

                if (total % 1000 == 0)
                    std::cerr << "\r" << total << " events" << std::flush;
            }
        }
    }

    std::cerr << "\r" << total << " events reconstructed -> " << output << "\n";
    tree->Write();

    //analyze Moller events, fill histograms, save to output ROOT file
    for (int i = 0; i < hycal_mollers.size(); i++) {
        float z = physics.GetMollerZdistance(hycal_mollers[i], Ebeam);
        physics.FillMollerZ(z);
        if(i >= 1) {
            auto center = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            physics.FillMollerXY(center[0], center[1]);
        } 

    }

    // write per-module histograms
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
    return 0;
}
