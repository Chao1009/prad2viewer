//=============================================================================
// Replay.cpp — EVIO to ROOT tree conversion
//=============================================================================

#include "Replay.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

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
}

void Replay::setupBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->Branch("event_num", &ev.event_num, "event_num/I");
    tree->Branch("trigger",   &ev.trigger,   "trigger/I");
    tree->Branch("timestamp", &ev.timestamp, "timestamp/L");
    tree->Branch("nch",       &ev.nch,       "nch/I");
    tree->Branch("crate",     ev.crate,      "crate[nch]/I");
    tree->Branch("slot",      ev.slot,       "slot[nch]/I");
    tree->Branch("channel",   ev.channel,    "channel[nch]/I");
    tree->Branch("module_id", ev.module_id,  "module_id[nch]/I");
    tree->Branch("nsamples",  ev.nsamples,   "nsamples[nch]/I");
    tree->Branch("ped_mean",  ev.ped_mean,   "ped_mean[nch]/F");
    tree->Branch("ped_rms",   ev.ped_rms,    "ped_rms[nch]/F");
    tree->Branch("integral",  ev.integral,   "integral[nch]/F");
    if (write_peaks) {
        tree->Branch("npeaks",       ev.npeaks,       "npeaks[nch]/I");
        tree->Branch("peak_height",  ev.peak_height,  Form("peak_height[nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("peak_time",    ev.peak_time,    Form("peak_time[nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("peak_integral",ev.peak_integral, Form("peak_integral[nch][%d]/F", fdec::MAX_PEAKS));
    }
}

bool Replay::Process(const std::string &input_evio, const std::string &output_root,
                     int max_events, bool write_peaks)
{
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

    fdec::EventData event;
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    int total = 0;

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            if (!ch.DecodeEvent(ie, event)) continue;
            if (max_events > 0 && total >= max_events) break;

            clearEvent(*ev);
            ev->event_num = event.info.event_number;
            ev->trigger   = event.info.trigger_bits;
            ev->timestamp = event.info.timestamp;

            int nch = 0;
            for (int r = 0; r < event.nrocs; ++r) {
                auto &roc = event.rocs[r];
                if (!roc.present) continue;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 16; ++c) {
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0 || nch >= kMaxCh) continue;

                        ev->crate[nch]   = roc.tag;
                        ev->slot[nch]    = s;
                        ev->channel[nch] = c;
                        ev->nsamples[nch] = cd.nsamples;

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

} // namespace analysis
