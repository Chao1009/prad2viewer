#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "WaveAnalyzer.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TChain.h>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

// per-event data (sized to worst case, reused)
static constexpr int kMaxCh = fdec::MAX_ROCS * fdec::MAX_SLOTS * 16;

using EventVars       = prad2::RawEventData;
void SetReadBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->SetBranchAddress("event_num", &ev.event_num);
    tree->SetBranchAddress("trigger",   &ev.trigger);
    tree->SetBranchAddress("timestamp", &ev.timestamp);
    tree->SetBranchAddress("hycal.nch",       &ev.nch);
    tree->SetBranchAddress("hycal.crate",     ev.crate);
    tree->SetBranchAddress("hycal.slot",      ev.slot);
    tree->SetBranchAddress("hycal.channel",   ev.channel);
    tree->SetBranchAddress("hycal.module_id", ev.module_id);
    tree->SetBranchAddress("hycal.nsamples",  ev.nsamples);
    tree->SetBranchAddress("hycal.samples",   ev.samples);
    tree->SetBranchAddress("hycal.ped_mean",  ev.ped_mean);
    tree->SetBranchAddress("hycal.ped_rms",   ev.ped_rms);
    tree->SetBranchAddress("hycal.integral",  ev.integral);
    if (write_peaks) {
        tree->SetBranchAddress("hycal.npeaks",       &ev.npeaks);
        tree->SetBranchAddress("hycal.peak_height",  ev.peak_height);
        tree->SetBranchAddress("hycal.peak_time",    ev.peak_time);
        tree->SetBranchAddress("hycal.peak_integral",ev.peak_integral);
    }
}

// ── Save first N waveforms with signal for a given module ─────────────
static void saveModuleWaveforms(TTree *tree, fdec::HyCalSystem &hycal,
                                EventVars &ev, int target_id,
                                int max_wf, TFile &outfile)
{
    int found = 0;
    TString dirName = Form("waveforms_module_%d", target_id);
    outfile.mkdir(dirName);
    outfile.cd(dirName);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries && found < max_wf; i++) {
        tree->GetEntry(i);
        for (int j = 0; j < ev.nch; j++) {
            const auto *mod = hycal.module_by_daq(ev.crate[j], ev.slot[j], ev.channel[j]);
            if (!mod || mod->id != target_id) continue;
            if (ev.npeaks[j] <= 0) continue;

            int ns = ev.nsamples[j];
            TString hname  = Form("wf_mod%d_ev%d", target_id, ev.event_num);
            TString htitle = Form("Module %d  Event %d;Sample;ADC", target_id, ev.event_num);
            TH1F *hwf = new TH1F(hname, htitle, ns, 0, ns);
            for (int s = 0; s < ns; s++)
                hwf->SetBinContent(s + 1, ev.samples[j][s]);
            hwf->Write();
            delete hwf;
            found++;
            break;
        }
    }
    outfile.cd();
    std::cerr << "Saved " << found << " waveforms for module " << target_id << "\n";
}

int main(int argc, char *argv[])
{
    std::string input = "cosmic.root";
    int view_module_id = -1;
    int opt;
    while ((opt = getopt(argc, argv, "m:")) != -1) {
        switch (opt) {
            case 'm': view_module_id = std::atoi(optarg); break;
        }
    }
    if (optind < argc) input = argv[optind];

    TChain *cosmic_chain = new TChain("events");
    for(int i = 0; i<=29; i++){
        std::string filename = Form("/data/stage6/prad_023419/prad_023419.000%02d.root", i);
        cosmic_chain->Add(filename.c_str());
    }

     TTree *tree = cosmic_chain;
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
        peak_hist_module[i] = new TH1F(name.c_str(), name.c_str(), 100, 0, 500);
    }
    TH1F *peak_hist_LG_module[1000];
    for (int i = 0; i < 1000; i++) {
        std::string name = "peak_LG_module_" + std::to_string(i+1);
        peak_hist_LG_module[i] = new TH1F(name.c_str(), name.c_str(), 100, 0, 500);
    }

    TH2F *cosmic_eventNum = new TH2F("cosmic_eventNum", "Cosmic Event Number", 34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);
    TH2F *cosmic_eventNum_LG = new TH2F("cosmic_eventNum_LG", "Cosmic Event Number for LG Modules", 34, -17.*38.15, 17.*38.15, 34, -17.*38.15, 17.*38.15);

    int event_num_module[3000] = {};

    int nentries = tree->GetEntries();
    for(int i = 0; i < 10000; i++){
        tree->GetEntry(i);
        std::cout << "Event " << ev.event_num << ": nch = " << ev.nch << "\r" << std::flush;
        if (ev.nch > 100) continue; // skip noisy events
        for (int j = 0; j < ev.nch; j++) {
            const auto *mod = hycal.module_by_daq(ev.crate[j], ev.slot[j], ev.channel[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev.npeaks[j] <= 0) continue;
            // Check module ID bounds
            event_num_module[mod->id]++;
            int module_id = mod->id-1000;
            if (module_id >= 1 && module_id <= 1156){
                float peak = ev.peak_integral[j][0];
                peak_hist_module[module_id-1]->Fill(peak);
                cosmic_eventNum->Fill(mod->x, mod->y);
            }
            else if(module_id < 0 && module_id >= -1000){
                int lg_module_id = module_id+1000;
                float peak = ev.peak_integral[j][0];
                peak_hist_LG_module[lg_module_id-1]->Fill(peak);
                cosmic_eventNum_LG->Fill(mod->x, mod->y);
            }
        }
    }

    TFile outfile("cosmic_23419_test_waveforms.root", "RECREATE");
    outfile.cd();
    outfile.mkdir("peak_histograms")->cd();
    for (int i = 0; i < 1156; i++) {
        if (peak_hist_module[i]->GetEntries() > 0) peak_hist_module[i]->Write();
    }
    outfile.cd();
    outfile.mkdir("peak_histograms_LG")->cd();
    for(int i = 0; i < 1000; i++) {
        if (peak_hist_LG_module[i]->GetEntries() > 0) peak_hist_LG_module[i]->Write();
    }

    outfile.cd();
    if (cosmic_eventNum->GetEntries() > 0) cosmic_eventNum->Write();
    if (cosmic_eventNum_LG->GetEntries() > 0) cosmic_eventNum_LG->Write();

    float peak[1156], rms[1156];
    for (int i = 0; i < 1156; i++) {
        if (peak_hist_module[i]->GetEntries() > 0) {
            float max = peak_hist_module[i]->GetBinCenter(peak_hist_module[i]->GetMaximumBin());

            peak_hist_module[i]->Fit("gaus", "Q", "r", 0, max + 40.);
            TF1 *fit = peak_hist_module[i]->GetFunction("gaus");
            if (fit) {
                peak[i] = fit->GetParameter(1); // mean
                rms[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                peak_hist_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas3/fit_module_" + std::to_string(i+1) + ".png").c_str());
                delete c;
            }
            else {
                peak[i] = 0.1;
                rms[i] = 0.1;
            }
        } else {
            peak[i] = 0.1;
            rms[i] = 0.1;
        }
    }

    TH1F *peak_module = new TH1F("peak_module", "Peak Integral by Module", 100, 0, 500);
    TH1F *rms_module = new TH1F("rms_module", "RMS of Peak Integral by Module", 100, 0, 400);
    for (int i = 0; i < 1156; i++) {
        peak_module->Fill(peak[i]);
        rms_module->Fill(rms[i]);
    }
    peak_module->Write();
    rms_module->Write();

    std::ofstream csv_out("cosmic_peak_23419.dat");
    csv_out << "ModuleID  PeakIntegral  RMS\n";
    for (int i = 0; i < 1156; i++) {
        csv_out << "W" << (i+1) << "  " << peak[i] << "  " << rms[i] << "\n";
    }
    for (int i = 0; i < 1000; i++) {
        csv_out << "G" << (i+1) << "  " << peak_hist_LG_module[i]->GetMean() << "  " << peak_hist_LG_module[i]->GetRMS() << "\n";
    }
    csv_out.close();

    std::ofstream rate_out("cosmic_eventNum_23419.dat");
    rate_out << "ModuleID  EventCount\n";
    for (int i = 0; i < 1156; i++) {
        rate_out << "W" << (i+1) << "  " << event_num_module[i+1000+1] << "\n";
    }
    for (int i = 0; i < 1000; i++) {
        rate_out << "G" << (i+1) << "  " << event_num_module[i+1] << "\n";
    }
    rate_out.close();

    outfile.Close();

    return 0;
}