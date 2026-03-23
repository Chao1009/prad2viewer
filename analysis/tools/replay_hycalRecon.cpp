//=============================================================================
// replay_hycalRecon — HyCal reconstruction replay with physics analysis
//
// Usage: replay_hycalRecon <input.evio> [-o output.root] [-c config.json]
//                          [-D daq_config.json] [-n max_events]
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

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

int main(int argc, char *argv[])
{
    std::string input, output, config_file, daq_config_file;
    std::string db_dir = DATABASE_DIR;
    int max_events = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:c:D:n:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'c': config_file = optarg; break;
            case 'D': daq_config_file = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty()) {
        std::cerr << "Usage: replay_hycalRecon <input.evio> [-o out.root] "
                  << "[-c config.json] [-D daq_config.json] [-n N]\n";
        return 1;
    }

    if (output.empty()) {
        output = input;
        auto dot = output.rfind('.');
        if (dot != std::string::npos) output = output.substr(0, dot);
        output += "_recon.root";
    }

    // --- setup ---
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    if (!daq_config_file.empty()) evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json");

    analysis::PhysicsTools physics(hycal);

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);
    if (ch.Open(input) != evc::status::success) {
        std::cerr << "Cannot open " << input << "\n";
        return 1;
    }

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
        std::ifstream dcf(daq_config_file);
        if (dcf.is_open()) {
            auto dcj = nlohmann::json::parse(dcf, nullptr, false, true);
            if (dcj.contains("crate_roc")) {
                for (auto &[k, v] : dcj["crate_roc"].items())
                    roc_to_crate[v.get<int>()] = std::stoi(k);
            }
        }
    }

    fdec::EventData event;
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    fdec::ClusterConfig cl_cfg;
    int total = 0;

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            if (!ch.DecodeEvent(ie, event)) continue;
            if (max_events > 0 && total >= max_events) break;

            fdec::HyCalCluster clusterer(hycal);
            clusterer.SetConfig(cl_cfg);

            // feed hits
            for (int r = 0; r < event.nrocs; ++r) {
                auto &roc = event.rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 16; ++c) {
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0) continue;
                        const auto *mod = hycal.module_by_daq(crate, s, c);
                        if (!mod || !mod->is_hycal()) continue;
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        if (wres.npeaks <= 0) continue;
                        float adc = wres.peaks[0].integral;
                        float energy = (mod->cal_factor > 0.) ?
                            static_cast<float>(mod->energize(adc)) : adc * 0.078f;
                        clusterer.AddHit(mod->index, energy);
                    }
                }
            }

            clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hits;
            clusterer.ReconstructHits(hits);

            ev_num = event.info.event_number;
            n_clusters = std::min((int)hits.size(), kMaxCl);
            for (int i = 0; i < n_clusters; ++i) {
                cl_x[i]       = hits[i].x;
                cl_y[i]       = hits[i].y;
                cl_energy[i]  = hits[i].energy;
                cl_nblocks[i] = hits[i].nblocks;
                cl_center[i]  = hits[i].center_id;
                physics.FillModuleEnergy(hits[i].center_id, hits[i].energy);
                physics.FillEnergyVsModule(hits[i].center_id, hits[i].energy);
            }
            tree->Fill();
            total++;

            if (total % 1000 == 0)
                std::cerr << "\r" << total << " events" << std::flush;
        }
        if (max_events > 0 && total >= max_events) break;
    }

    std::cerr << "\r" << total << " events reconstructed -> " << output << "\n";
    tree->Write();

    // write per-module histograms
    outfile.mkdir("module_energy");
    outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); ++i) {
        TH1F *h = physics.GetModuleHist(i);
        if (h && h->GetEntries() > 0) h->Write();
    }
    outfile.cd();
    if (physics.GetEnergyVsModuleHist())
        physics.GetEnergyVsModuleHist()->Write();

    outfile.Close();
    return 0;
}
