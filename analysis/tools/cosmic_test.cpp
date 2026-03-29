#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "WaveAnalyzer.h"

#include <TFile.h>
#include <TTree.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

// per-event data (sized to worst case, reused)
static constexpr int kMaxCh = fdec::MAX_ROCS * fdec::MAX_SLOTS * 16;

struct EventVars {
    int     event_num = 0;
    int     trigger = 0;
    Long64_t timestamp = 0;
    int     nch = 0;
    int     crate[kMaxCh] = {};
    int     slot[kMaxCh] = {};
    int     channel[kMaxCh] = {};
    int     module_id[kMaxCh] = {};
    int     nsamples[kMaxCh] = {};
    int     samples[kMaxCh][fdec::MAX_SAMPLES] = {};
    float   ped_mean[kMaxCh] = {};
    float   ped_rms[kMaxCh] = {};
    float   integral[kMaxCh] = {};
    int     npeaks[kMaxCh] = {};
    float   peak_height[kMaxCh][fdec::MAX_PEAKS] = {};
    float   peak_time[kMaxCh][fdec::MAX_PEAKS] = {};
    float   peak_integral[kMaxCh][fdec::MAX_PEAKS] = {};
};

void SetReadBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->SetBranchAddress("event_num", &ev.event_num);
    tree->SetBranchAddress("trigger",   &ev.trigger);
    tree->SetBranchAddress("timestamp", &ev.timestamp);
    tree->SetBranchAddress("nch",       &ev.nch);
    tree->SetBranchAddress("crate",     ev.crate);
    tree->SetBranchAddress("slot",      ev.slot);
    tree->SetBranchAddress("channel",   ev.channel);
    tree->SetBranchAddress("module_id", ev.module_id);
    tree->SetBranchAddress("nsamples",  ev.nsamples);
    tree->SetBranchAddress("ped_mean",  ev.ped_mean);
    tree->SetBranchAddress("ped_rms",   ev.ped_rms);
    tree->SetBranchAddress("integral",  ev.integral);
    if (write_peaks) {
        tree->SetBranchAddress("npeaks",       ev.npeaks);
        tree->SetBranchAddress("peak_height",  ev.peak_height);
        tree->SetBranchAddress("peak_time",    ev.peak_time);
        tree->SetBranchAddress("peak_integral",ev.peak_integral);
    }
}

int main(int argc, char *argv[])
{
    std::string input = "cosmic.root";
    if (argc > 1) input = argv[1];

    TFile *infile = TFile::Open(input.c_str(), "READ");
    if (!infile || !infile->IsOpen()) {
        std::cerr << "Cannot open " << input << "\n";
        return 1;
    }
    TTree *tree = (TTree *)infile->Get("events");
    if (!tree) {
        std::cerr << "Cannot find TTree 'events' in " << input << "\n";
        return 1;
    }
    EventVars ev;
    SetReadBranches(tree, ev, true);

    //setup for reconstruction
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    std::string daq_config_file = db_dir + "/daq_config.json"; // default DAQ config for PRad2
    if (!daq_config_file.empty()) evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json");

    TH1F *peak_hist_module[1156];
    for (int i = 0; i < 1156; i++) {
        std::string name = "peak_module_" + std::to_string(i+1);
        peak_hist_module[i] = new TH1F(name.c_str(), name.c_str(), 200, 0, 4000);
    }

    TH2F *cosmic_visual[100];
    for (int i = 0; i < 100; i++) {
        std::string name = "cosmic_visual_" + std::to_string(i);
        cosmic_visual[i] = new TH2F(name.c_str(), name.c_str(), 34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);
    } 

    int nentries = tree->GetEntries();
    for(int i = 0; i < nentries; i++){
        tree->GetEntry(i);
        std::cout << "Event " << ev.event_num << ": nch = " << ev.nch << "\n";
        if (ev.nch > 100) continue; // skip noisy events
        for (int j = 0; j < ev.nch; j++) {
            const auto *mod = hycal.module_by_daq(ev.crate[j], ev.slot[j], ev.channel[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev.npeaks[j] <= 0) continue;
            float peak = ev.peak_integral[j][0];
            peak_hist_module[mod->id-1]->Fill(peak);
            if(i < 100) cosmic_visual[i]->Fill(mod->x, mod->y, peak);
        }
    }

    TFile outfile("cosmic_analysis.root", "RECREATE");
    for (int i = 0; i < 1156; i++) {
        if (peak_hist_module[i]->GetEntries() > 0) peak_hist_module[i]->Write();
    }
    for (int i = 0; i < 100; i++) {
        if (cosmic_visual[i]->GetEntries() > 0) cosmic_visual[i]->Write();
    }
    outfile.Close();
    infile->Close();

    return 0;
}