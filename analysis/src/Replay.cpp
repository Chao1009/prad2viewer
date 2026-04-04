//=============================================================================
// Replay.cpp — EVIO to ROOT tree conversion
//=============================================================================

#include "Replay.h"
#include "DaqConfig.h"
#include "HyCalSystem.h"
#include "GemSystem.h"
#include "HyCalCluster.h"
#include "GemCluster.h"
#include "MatchingTools.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace analysis {

void Replay::LoadDaqMap(const std::string &json_path)
{
    std::ifstream f(json_path);
    if (!f.is_open()) {
        std::cerr << "Replay: cannot open DAQ map: " << json_path << "\n";
        return;
    }
    auto j = json::parse(f, nullptr, false, true);
    if (j.is_array()) {
        for (auto &entry : j) {
            std::string name = entry.value("name", "");
            int crate   = entry.value("crate", -1);
            int slot    = entry.value("slot", -1);
            int channel = entry.value("channel", -1);
            if (!name.empty() && crate >= 0)
                daq_map_[std::to_string(crate) + "_" + std::to_string(slot) +
                         "_" + std::to_string(channel)] = name;
        }
    }
    std::cerr << "Replay: loaded " << daq_map_.size() << " DAQ map entries\n";
}

std::string Replay::moduleName(int roc, int slot, int ch) const
{
    auto it = daq_map_.find(std::to_string(roc) + "_" + std::to_string(slot) +
                            "_" + std::to_string(ch));
    return (it != daq_map_.end()) ? it->second : "";
}

int Replay::moduleID(int roc, int slot, int ch) const
{
    auto name = moduleName(roc, slot, ch);
    if (name.empty()) return -1;
    if(name[0] == 'G') return std::stoi(name.substr(1));
    else if(name[0] == 'W') return std::stoi(name.substr(1))+1000;
    else return -1; // Unknown module type
}

float Replay::computeIntegral(const fdec::ChannelData &cd, float pedestal) const
{
    float sum = 0.f;
    for (int i = 0; i < cd.nsamples; ++i)
        sum += cd.samples[i] - pedestal;
    return sum;
}

void Replay::clearEvent(EventVars &ev)
{
    ev.event_num = 0;
    ev.trigger = 0;
    ev.timestamp = 0;
    ev.nch = 0;
    ev.gem_nch = 0;
}

void Replay::clearReconEvent(EventVars_Recon &ev)
{
    ev.event_num = 0;
    ev.trigger_bits = 0;
    ev.timestamp = 0;
    ev.total_energy = 0.f;
    ev.n_clusters = 0;
    ev.n_gem_hits = 0;
    ev.match_num = 0;
}

void Replay::setupBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->Branch("event_num", &ev.event_num, "event_num/i");
    tree->Branch("trigger",   &ev.trigger,   "trigger/i");
    tree->Branch("timestamp", &ev.timestamp, "timestamp/L");
    tree->Branch("hycal.nch",       &ev.nch,       "nch/I");
    tree->Branch("hycal.crate",     ev.crate,      "crate[nch]/b");
    tree->Branch("hycal.slot",      ev.slot,       "slot[nch]/b");
    tree->Branch("hycal.channel",   ev.channel,    "channel[nch]/b");
    tree->Branch("hycal.module_id", ev.module_id,  "module_id[nch]/s");
    tree->Branch("hycal.nsamples",  ev.nsamples,   "nsamples[nch]/b");
    tree->Branch("hycal.samples",   ev.samples,    Form("samples[nch][%d]/s", fdec::MAX_SAMPLES));
    tree->Branch("hycal.ped_mean",  ev.ped_mean,   "ped_mean[nch]/F");
    tree->Branch("hycal.ped_rms",   ev.ped_rms,    "ped_rms[nch]/F");
    tree->Branch("hycal.integral",  ev.integral,   "integral[nch]/F");
    if (write_peaks) {
        tree->Branch("hycal.npeaks",       &ev.npeaks,       "npeaks[nch]/b");
        tree->Branch("hycal.peak_height",  ev.peak_height,  Form("peak_height[nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_time",    ev.peak_time,    Form("peak_time[nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_integral",ev.peak_integral, Form("peak_integral[nch][%d]/F", fdec::MAX_PEAKS));
    }
    //GEM part
    tree->Branch("gem.nch",        &ev.gem_nch,   "gem_nch/I");
    tree->Branch("gem.mpd_crate",  ev.mpd_crate,  "mpd_crate[gem_nch]/b");
    tree->Branch("gem.mpd_fiber",  ev.mpd_fiber,  "mpd_fiber[gem_nch]/b");
    tree->Branch("gem.apv",        ev.apv,        "apv[gem_nch]/b");
    tree->Branch("gem.strip",        ev.strip,        "strip[gem_nch]/b");
    tree->Branch("gem.ssp_samples",  ev.ssp_samples,  Form("ssp_samples[gem_nch][%d]/S", ssp::SSP_TIME_SAMPLES));
    tree->Branch("n_ssp_triggers", &ev.n_ssp_triggers, "n_ssp_triggers/b");
    tree->Branch("ssp_trigger_tags", ev.ssp_trigger_tags, Form("ssp_trigger_tags[n_ssp_triggers][%d]/i", ssp::SSP_TIME_SAMPLES));
}

void Replay::setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/i");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    tree->Branch("total_energy", &ev.total_energy, "total_energy/F");
    // HyCal cluster branches
    // detector coordinate system (crystal surface)
    tree->Branch("n_clusters",   &ev.n_clusters,   "n_clusters/I");
    tree->Branch("cl_x",         ev.cl_x,          "cl_x[n_clusters]/F");
    tree->Branch("cl_y",         ev.cl_y,          "cl_y[n_clusters]/F");
    tree->Branch("cl_z",         ev.cl_z,          "cl_z[n_clusters]/F");
    tree->Branch("cl_energy",    ev.cl_energy,     "cl_energy[n_clusters]/F");
    tree->Branch("cl_nblocks",   ev.cl_nblocks,    "cl_nblocks[n_clusters]/b");
    tree->Branch("cl_center",    ev.cl_center,     "cl_center[n_clusters]/s");
    tree->Branch("cl_flag",      ev.cl_flag,       "cl_flag[n_clusters]/i");
    // GEM part
    //detector coordinate system (GEM plane)
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
    // the x,y,z positions in this part are in the target and beam center coordinate system
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

    tree->Branch("n_ssp_triggers", &ev.n_ssp_triggers, "n_ssp_triggers/b");
    tree->Branch("ssp_trigger_tags", ev.ssp_trigger_tags, Form("ssp_trigger_tags[n_ssp_triggers][%d]/i", ssp::SSP_TIME_SAMPLES));
}

bool Replay::Process(const std::string &input_evio, const std::string &output_root,
                     int max_events, bool write_peaks , const std::string &daq_config_file)
{
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
    else {
        std::cerr << "No DAQ config file provided, ROC tag to crate mapping will be unavailable.\n";
    }

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.Open(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TTree *tree = new TTree("events", "PRad2 replay data");
    //EventVars ev;
    auto ev = std::make_unique<EventVars>();
    setupBranches(tree, *ev, write_peaks);

    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    int total = 0;

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != evc::EventType::Physics) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            ssp_evt->clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;
            if (max_events > 0 && total >= max_events) break;

            clearEvent(*ev);
            ev->event_num = event->info.event_number;
            ev->trigger   = event->info.trigger_bits;
            ev->timestamp = event->info.timestamp;

            // decode HyCal FADC250 data
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                int crate;
                if (cit == roc_to_crate.end()) crate = roc.tag;
                else crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 16; ++c) {
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0 || nch >= prad2::kMaxChannels) continue;

                        ev->crate[nch]   = crate;
                        ev->slot[nch]    = s;
                        ev->channel[nch] = c;
                        //ev->module_id[nch] = moduleID(roc.tag, s, c);
                        ev->nsamples[nch] = cd.nsamples;
                        for (int i = 0; i < cd.nsamples && i < fdec::MAX_SAMPLES; ++i)
                            ev->samples[nch][i] = cd.samples[i];

                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        ev->ped_mean[nch] = wres.ped.mean;
                        ev->ped_rms[nch]  = wres.ped.rms;
                        ev->integral[nch] = computeIntegral(cd, wres.ped.mean);

                        if (write_peaks) {
                            ev->npeaks[nch] = wres.npeaks;
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                ev->peak_height[nch][p]   = wres.peaks[p].height;
                                ev->peak_time[nch][p]     = wres.peaks[p].time;
                                ev->peak_integral[nch][p] = wres.peaks[p].integral;
                            }
                        }
                        nch++;
                    }
                }
            }
            ev->nch = nch;

            // decode GEM SSP data
            int gem_ch = 0;
            for (int m = 0; m < ssp_evt->nmpds; ++m) {
                auto &mpd = ssp_evt->mpds[m];
                if (!mpd.present) continue;
                for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
                    auto &apv = mpd.apvs[a];
                    if (!apv.present) continue;
                    int idx = -1; // find APV index in GemSystem if needed
                    for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
                        if (!apv.hasStrip(s)) continue;
                        if (gem_ch >= prad2::kMaxGemStrips) continue;
                        
                        ev->mpd_crate[gem_ch] = mpd.crate_id;
                        ev->mpd_fiber[gem_ch] = mpd.mpd_id;
                        ev->apv[gem_ch]       = a;
                        ev->strip[gem_ch]     = s;
                        for (int t = 0; t < ssp::SSP_TIME_SAMPLES; t++)
                            ev->ssp_samples[gem_ch][t] = apv.strips[s][t];

                        gem_ch++;
                    }
                }
            }
            ev->gem_nch = gem_ch; // total channels = HyCal + GEM
            tree->Fill();
            total++;

            if (total % 1000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
        if (max_events > 0 && total >= max_events) break;
    }

    std::cerr << "\rReplay: " << total << " events written to " << output_root << "\n";
    tree->Write();
    delete outfile;
    return true;
}

bool Replay::ProcessWithRecon(const std::string &input_evio, const std::string &output_root,
                                const std::string &daq_config_file, const std::string &gem_ped_file,
                                const float zerosup_override, bool prad1)
{
    // Similar to Process(), but with HyCal reconstruction and GEM hit reconstruction
    // before filling the ROOT tree.
    // The main differences are:
    // - After decoding, we run the HyCal clusterer to reconstruct clusters and hits.
    // - We also run the GemSystem reconstruction to get GEM hits.
    // - We fill a different TTree with reconstructed quantities instead of raw data.

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;

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
    else {
        std::cerr << "No DAQ config file provided, ROC tag to crate mapping will be unavailable.\n";
    }

    // Setup HyCal system and clusterer
    fdec::HyCalSystem hycal;
     std::string daq_map_file = db_dir + "/daq_map.json";
    if(prad1 == true)
        daq_map_file = db_dir + "/prad1/prad_daq_map.json";
    hycal.Init(db_dir + "/hycal_modules.json", daq_map_file);
    if(prad1 == true) evc::load_pedestals(db_dir + "/prad1/adc1881m_pedestals.json", daq_cfg_);
    std::string calib_file = db_dir + "/prad1/prad_calibration.json";
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    fdec::HyCalCluster clusterer(hycal);
    fdec::ClusterConfig cl_cfg;
    clusterer.SetConfig(cl_cfg);

    MatchingTools matching;

    // Initialize GEM system and clusterer
    std::unique_ptr<gem::GemSystem> gem_sys;
    std::unique_ptr<gem::GemCluster> gem_clusterer;
if(!prad1){
    gem_sys = std::make_unique<gem::GemSystem>();
    std::string gem_map_file = db_dir + "/gem_map.json";
    gem_sys->Init(gem_map_file);
    std::cerr << "GEM map  : " << gem_map_file
                << " (" << gem_sys->GetNDetectors() << " detectors)\n";

    if (!gem_ped_file.empty()) {
        gem_sys->LoadPedestals(gem_ped_file);
        std::cerr << "GEM peds : " << gem_ped_file << "\n";
    }

    if (zerosup_override >= 0.f) {
        gem_sys->SetZeroSupThreshold(zerosup_override);
        std::cerr << "Zero-sup : " << zerosup_override << " sigma (override)\n";
    }
    
    gem_clusterer = std::make_unique<gem::GemCluster>();
}
    //open EVIO file and output ROOT file
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.Open(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    // create TTree and branches for reconstructed data
    TTree *tree = new TTree("recon", "PRad2 replay reconstruction");
    auto ev = std::make_unique<EventVars_Recon>();
    setupReconBranches(tree, *ev);

    //initialize tools for event decoder and cluster reconstruction
    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    
    int total = 0;

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != evc::EventType::Physics) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            ssp_evt->clear();
            clusterer.Clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;

            clearReconEvent(*ev);
            ev->event_num = event->info.event_number;
            ev->trigger_bits = event->info.trigger_bits;
            ev->timestamp = event->info.timestamp;

            bool cluster = false, LMS = false;
            if(ev->trigger_bits & (1 << 8) ||
               ev->trigger_bits & (1 << 9) ||
               ev->trigger_bits & (1 << 10)||
               ev->trigger_bits & (1 << 11)) cluster = true;

            if(ev->trigger_bits & (1 << 24)) LMS = true;
            int fire_num = 0;
            // decode and reconstruct HyCal data
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 160; ++c) { //should be 16, a bigger number to adapt PRad1 data
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
                        hycal.SetCalibConstant(mod->id, 0.2);
                        float energy = static_cast<float>(mod->energize(adc));
                        clusterer.AddHit(mod->index, energy);
                        ev->total_energy += energy;
                        fire_num++;
                    }
                }
            }
            if( !(cluster && !LMS && fire_num < 100) ) continue; //apply trigger and fire number condition to select good events for reconstruction
            clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hits;
            clusterer.ReconstructHits(hits);
            //HyCal event reconstrued, fill root tree and histograms
            ev->n_clusters = std::min((int)hits.size(), prad2::kMaxClusters);
            for (int i = 0; i < ev->n_clusters; ++i) {
                ev->cl_x[i]       = hits[i].x;
                ev->cl_y[i]       = hits[i].y;
                ev->cl_z[i]       = PhysicsTools::GetShowerDepth(hits[i].center_id, hits[i].energy);
                ev->cl_energy[i]  = hits[i].energy;
                ev->cl_nblocks[i] = hits[i].nblocks;
                ev->cl_center[i]  = hits[i].center_id;
                ev->cl_flag[i]    = hits[i].flag;
            }

            //decode GEM data and reconstruct GEM hits
        if(!prad1){
            if (gem_sys) {
                gem_sys->Clear();
                gem_sys->ProcessEvent(*ssp_evt);
                if (gem_clusterer)
                    gem_sys->Reconstruct(*gem_clusterer);
            }
            else {
                ev->n_gem_hits = 0;
                std::cerr << "Warning: GEM system not initialized, skipping GEM reconstruction\n";
            }
            auto &all_hits = gem_sys->GetAllHits();
            ev->n_gem_hits = std::min((int)all_hits.size(), prad2::kMaxGemHits);
            for (int i = 0; i < ev->n_gem_hits; i++) {
                auto &h = all_hits[i];
                ev->det_id[i] = h.det_id;
                ev->gem_x[i] = h.x;
                ev->gem_y[i] = h.y;
                ev->gem_x_charge[i] = h.x_charge;
                ev->gem_y_charge[i] = h.y_charge;
                ev->gem_x_peak[i] = h.x_peak;
                ev->gem_y_peak[i] = h.y_peak;
                ev->gem_x_size[i] = h.x_size;
                ev->gem_y_size[i] = h.y_size;
            }

            // Perform matching between HyCal clusters and GEM hits
            //store all the hits on HyCal and GEMs in this event
            std::vector<HCHit> hc_hits;
            std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
            for (int i = 0; i < ev->n_clusters; ++i)
                hc_hits.push_back({ev->cl_x[i], ev->cl_y[i], ev->cl_z[i], ev->cl_energy[i], ev->cl_center[i], ev->cl_flag[i]});
            for (int i = 0; i < ev->n_gem_hits; ++i)
                gem_hits[ev->det_id[i]].push_back(GEMHit{ev->gem_x[i], ev->gem_y[i], 0.f, ev->det_id[i]});
            //transform the coordinates of detector data
            // (now the default position, TO DO: read from database)
            float gem_z[4] = {5407.+39.71/2., 5407.-39.71/2.,
                              5807.+39.71/2., 5807.-39.71/2.}; //mm, the center of the two GEMs
            TransformDetData(hc_hits, 0.f, 0.f, 6225.f);
            for(int i = 0; i < 4; ++i)
                TransformDetData(gem_hits[i], 0.f, 0.f, gem_z[i]);

            //matching.SetMatchRange(5.0f); // matching radius in mm, 15mm default
            matching.SetSquareSelection(true); // use square cut instead of circular cut
            std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);

            ev->match_num = std::min((int)matched_hits.size(), prad2::kMaxClusters);
            for (int i = 0; i < ev->match_num; i++){
                ev->matchHC_x[i] = matched_hits[i].hycal_hit.x;
                ev->matchHC_y[i] = matched_hits[i].hycal_hit.y;
                ev->matchHC_z[i] = matched_hits[i].hycal_hit.z;
                ev->matchHC_energy[i] = matched_hits[i].hycal_hit.energy;
                ev->matchHC_center[i] = matched_hits[i].hycal_hit.center_id;
                ev->matchHC_flag[i] = matched_hits[i].hycal_hit.flag;
                for (int j = 0; j < 2; j++) {
                    ev->matchG_x[i][j] = matched_hits[i].gem[j].x;
                    ev->matchG_y[i][j] = matched_hits[i].gem[j].y;
                    ev->matchG_z[i][j] = matched_hits[i].gem[j].z;
                    ev->matchG_det_id[i][j] = matched_hits[i].gem[j].det_id;
                }
            }

        } //end of if(PRad1)
            tree->Fill();
            total++;
            if (total % 1000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
    }
    std::cerr << "\rReplay: " << total << " events reconstructed -> " << output_root << "\n";
    tree->Write();
    delete outfile;

    return true;
}

} // namespace analysis
