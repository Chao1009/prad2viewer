#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DaqConfig.h"
#include "WaveAnalyzer.h"
#include "EventData.h"

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
    tree->SetBranchAddress("trigger_bits",   &ev.trigger_bits);
    tree->SetBranchAddress("timestamp", &ev.timestamp);
    tree->SetBranchAddress("hycal.nch",       &ev.nch);
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

const int LG_num = 76;
int LG_module_id[LG_num] = 
{156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174,
 186, 216, 246, 276, 306, 336, 366, 396, 426, 456, 486, 516, 546, 576, 606, 636, 666, 696, 726,
 175, 205, 235, 265, 295, 325, 355, 385, 415, 445, 475, 505, 535, 565, 595, 625, 655, 685, 715,
 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745
};

int main(int argc, char *argv[])
{
    std::string input;
    std::string in_json;
    int run_number = -1, file_number = -1;
    int opt;
    while ((opt = getopt(argc, argv, "r:n:j:")) != -1) {
        switch (opt) {
            case 'r': run_number  = std::atoi(optarg); break;
            case 'n': file_number = std::atoi(optarg); break;
            case 'j': in_json    = optarg; break;
            default:
                std::cerr << "Usage: " << argv[0] << " [-r run_number] [-n file_number] [-j existing_json]\n";
                return 1;
        }
    }
    if (optind < argc) input = argv[optind];

    TChain *cosmic_chain = new TChain("events");
    for(int i = 0; i<=file_number-1; i++){
        std::string filename = Form("prad_023%d.000%02d_raw.root", run_number, i);
        cosmic_chain->Add(filename.c_str());
    }

     TTree *tree = cosmic_chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'events' \n";
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
    TH1F *peakHeight_hist_module[1156];
    TH1F *lms_hist_module[1156];
    TH1F *lms_height_hist_module[1156];
    for (int i = 0; i < 1156; i++) {
        std::string name = "peak_module_" + std::to_string(i+1);
        peak_hist_module[i] = new TH1F(name.c_str(), name.c_str(), 80, 0, 800);
        std::string name_height = "peakHeight_module_" + std::to_string(i+1);
        peakHeight_hist_module[i] = new TH1F(name_height.c_str(), name_height.c_str(), 80, 0, 200);
        std::string name_lms = "lms_module_" + std::to_string(i+1);
        lms_hist_module[i] = new TH1F(name_lms.c_str(), name_lms.c_str(), 80*4, 0, 800*4);
        std::string name_lms_height = "lms_height_module_" + std::to_string(i+1);
        lms_height_hist_module[i] = new TH1F(name_lms_height.c_str(), name_lms_height.c_str(), 80*4, 0, 200*4);
    }
    TH1F *peak_hist_LG_module[LG_num];
    TH1F *peakHeight_hist_LG_module[LG_num];
    TH1F *lms_hist_LG_module[LG_num];
    TH1F *lms_height_hist_LG_module[LG_num];
    for (int i = 0; i < LG_num; i++) {
        std::string name = "peak_LG_module_" + std::to_string(LG_module_id[i]);
        peak_hist_LG_module[i] = new TH1F(name.c_str(), name.c_str(), 80, 0, 800);
        std::string name_height = "peakHeight_LG_module_" + std::to_string(LG_module_id[i]);
        peakHeight_hist_LG_module[i] = new TH1F(name_height.c_str(), name_height.c_str(), 80, 0, 200);
        std::string name_lms = "lms_LG_module_" + std::to_string(LG_module_id[i]);
        lms_hist_LG_module[i] = new TH1F(name_lms.c_str(), name_lms.c_str(), 80*4, 0, 800*4);
        std::string name_lms_height = "lms_height_LG_module_" + std::to_string(LG_module_id[i]);
        lms_height_hist_LG_module[i] = new TH1F(name_lms_height.c_str(), name_lms_height.c_str(), 80*4, 0, 200*4);
    }

    TH2F *cosmic_eventNum = new TH2F("cosmic_eventNum", "Cosmic Event Number", 34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);
    TH2F *cosmic_eventNum_LG = new TH2F("cosmic_eventNum_LG", "Cosmic Event Number for LG Modules", 34, -17.*38.15, 17.*38.15, 34, -17.*38.15, 17.*38.15);

    int event_num_module[3000] = {};

    int nentries = tree->GetEntries();
    for(int i = 0; i < nentries; i++){
        tree->GetEntry(i);
        std::cout << "Event " << ev.event_num << ": nch = " << ev.nch << "\r" << std::flush;
        if (ev.nch > 900) {
                for (int j = 0; j < ev.nch; j++) {
                    const auto *mod = hycal.module_by_id(ev.module_id[j]);
                    if (!mod || !mod->is_hycal()) continue;
                    if (ev.npeaks[j] <= 0) continue;
                    int module_id = mod->id-1000;
                    if (module_id >= 1 && module_id <= 1156){
                        if(ev.npeaks[j] != 1) continue;
                        lms_hist_module[module_id-1]->Fill(ev.peak_integral[j][0]);
                        lms_height_hist_module[module_id-1]->Fill(ev.peak_height[j][0]);
                    }
                    for(int k = 0; k < LG_num; k++){
                        if(mod->id == LG_module_id[k]) {
                            if(ev.npeaks[j] != 1) continue;
                            lms_hist_LG_module[k]->Fill(ev.peak_integral[j][0]);
                            lms_height_hist_LG_module[k]->Fill(ev.peak_height[j][0]);
                        }
                    }
                }
        }
        else if (ev.nch > 4 && ev.nch <= 70) {
            for (int j = 0; j < ev.nch; j++) {
                const auto *mod = hycal.module_by_id(ev.module_id[j]);
                if (!mod || !mod->is_hycal()) continue;
                if (ev.npeaks[j] <= 0) continue;
                event_num_module[mod->id]++;
                int module_id = mod->id-1000;
                if (module_id >= 1 && module_id <= 1156){
                    if(ev.npeaks[j] != 1) continue;
                    peak_hist_module[module_id-1]->Fill(ev.peak_integral[j][0]);
                    peakHeight_hist_module[module_id-1]->Fill(ev.peak_height[j][0]);
                    cosmic_eventNum->Fill(mod->x, mod->y);
                }
                for(int k = 0; k < LG_num; k++){
                    if(mod->id == LG_module_id[k]) {
                        if(ev.npeaks[j] != 1) continue;
                        peak_hist_LG_module[k]->Fill(ev.peak_integral[j][0]);
                        peakHeight_hist_LG_module[k]->Fill(ev.peak_height[j][0]);
                        cosmic_eventNum_LG->Fill(mod->x, mod->y);
                    }
                }
            }
        }
    }

    TFile outfile(Form("cosmic_run_%d.root", run_number), "RECREATE");
    outfile.cd();
    outfile.mkdir("peak_histograms")->cd();
    for (int i = 0; i < 1156; i++) {
        if (peak_hist_module[i]->GetEntries() > 0) peak_hist_module[i]->Write();
        if (peakHeight_hist_module[i]->GetEntries() > 0) peakHeight_hist_module[i]->Write();
    }
    outfile.cd();
    outfile.mkdir("peak_histograms_LG")->cd();
    for(int i = 0; i < LG_num; i++) {
        if (peak_hist_LG_module[i]->GetEntries() > 0) peak_hist_LG_module[i]->Write();
        if (peakHeight_hist_LG_module[i]->GetEntries() > 0) peakHeight_hist_LG_module[i]->Write();
    }

    outfile.cd();
    if (cosmic_eventNum->GetEntries() > 0) cosmic_eventNum->Write();
    if (cosmic_eventNum_LG->GetEntries() > 0) cosmic_eventNum_LG->Write();

    float peak[1156], rms[1156];
    for (int i = 0; i < 1156; i++) {
        if (peak_hist_module[i]->GetEntries() > 0) {
            float max_l1 = peak_hist_module[i]->GetBinContent(peak_hist_module[i]->GetMaximumBin()-1);
            float max_l2 = peak_hist_module[i]->GetBinContent(peak_hist_module[i]->GetMaximumBin()-2);
            if(max_l1 < 0.5 * peak_hist_module[i]->GetBinContent(peak_hist_module[i]->GetMaximumBin())
               || max_l2 < 0.4 * peak_hist_module[i]->GetBinContent(peak_hist_module[i]->GetMaximumBin())){
                for(int b = 1; b <= 8; b++){
                    peak_hist_module[i]->SetBinContent(b, 0);
                }
            }
            float max = peak_hist_module[i]->GetBinCenter(peak_hist_module[i]->GetMaximumBin());
            peak_hist_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = peak_hist_module[i]->GetFunction("gaus");
            if (fit) {
                peak[i] = fit->GetParameter(1); // mean
                rms[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                peak_hist_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_module_" + std::to_string(i+1) + ".png").c_str());
                delete c;
            }
            else {
                peak[i] = 0.1;
                rms[i] = 1e5;
            }
        } else {
            peak[i] = 0.1;
            rms[i] = 1e5;
        }
    }

    float peak_height[1156], rms_height[1156];
    for (int i = 0; i < 1156; i++) {
        if (peakHeight_hist_module[i]->GetEntries() > 0) {
            float max_l1 = peakHeight_hist_module[i]->GetBinContent(peakHeight_hist_module[i]->GetMaximumBin()-1);
            float max_l2 = peakHeight_hist_module[i]->GetBinContent(peakHeight_hist_module[i]->GetMaximumBin()-2);
            if(max_l1 < 0.5 * peakHeight_hist_module[i]->GetBinContent(peakHeight_hist_module[i]->GetMaximumBin())
               || max_l2 < 0.4 * peakHeight_hist_module[i]->GetBinContent(peakHeight_hist_module[i]->GetMaximumBin())){
                for(int b = 1; b <= 8; b++){
                    peakHeight_hist_module[i]->SetBinContent(b, 0);
                }
            }
            float max = peakHeight_hist_module[i]->GetBinCenter(peakHeight_hist_module[i]->GetMaximumBin());
            if(max < 20){
                peak_height[i] = peakHeight_hist_module[i]->GetMean();
                rms_height[i] = 1.;
            }
            peakHeight_hist_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = peakHeight_hist_module[i]->GetFunction("gaus");
            if (fit) {
                peak_height[i] = fit->GetParameter(1); // mean
                rms_height[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                peakHeight_hist_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_peakHeight_module_" + std::to_string(i+1) + ".png").c_str());
                delete c;
            }
            else {
                peak_height[i] = 0.1;
                rms_height[i] = 1e5;
            }
        } else {
            peak_height[i] = 0.1;
            rms_height[i] = 1e5;
        }
    }

    float peak_LG[LG_num], rms_LG[LG_num];
    for (int i = 0; i < LG_num; i++) {
        if (peak_hist_LG_module[i]->GetEntries() > 0) {
            float max = peak_hist_LG_module[i]->GetBinCenter(peak_hist_LG_module[i]->GetMaximumBin());
            peak_hist_LG_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = peak_hist_LG_module[i]->GetFunction("gaus");
            if (fit) {
                peak_LG[i] = fit->GetParameter(1); // mean
                rms_LG[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                peak_hist_LG_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_LG_module_" + std::to_string(LG_module_id[i]) + ".png").c_str());
                delete c;
            }
            else {
                peak_LG[i] = 0.1;
                rms_LG[i] = 1e5;
            }
        } else {
            peak_LG[i] = 0.1;
            rms_LG[i] = 1e5;
        }
    }

    float peakHeight_LG[LG_num], rms_height_LG[LG_num];
    for (int i = 0; i < LG_num; i++) {
        if (peakHeight_hist_LG_module[i]->GetEntries() > 0) {
            float max = peakHeight_hist_LG_module[i]->GetBinCenter(peakHeight_hist_LG_module[i]->GetMaximumBin());
            peakHeight_hist_LG_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = peakHeight_hist_LG_module[i]->GetFunction("gaus");
            if (fit) {
                peakHeight_LG[i] = fit->GetParameter(1); // mean
                rms_height_LG[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                peakHeight_hist_LG_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_peakHeight_LG_module_" + std::to_string(LG_module_id[i]) + ".png").c_str());
                delete c;
            }
            else {
                peakHeight_LG[i] = 0.1;
                rms_height_LG[i] = 1e5;
            }
        } else {
            peakHeight_LG[i] = 0.1;
            rms_height_LG[i] = 1e5;
        }
    }

    float lms_peak[1156], lms_rms[1156];
    for (int i = 0; i < 1156; i++) {
        if (lms_hist_module[i]->GetEntries() > 0) {
            float max = lms_hist_module[i]->GetBinCenter(lms_hist_module[i]->GetMaximumBin());
            lms_hist_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = lms_hist_module[i]->GetFunction("gaus");
            if (fit) {
                lms_peak[i] = fit->GetParameter(1); // mean
                lms_rms[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                lms_hist_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_LMS_module_" + std::to_string(i+1) + ".png").c_str());
                delete c;
            }
            else {
                lms_peak[i] = 0.1;
                lms_rms[i] = 1e5;
            }
        } else {
            lms_peak[i] = 0.1;
            lms_rms[i] = 1e5;
        }
    }

     float lms_peak_LG[LG_num], lms_rms_LG[LG_num];
    for (int i = 0; i < LG_num; i++) {
        if (lms_hist_LG_module[i]->GetEntries() > 0) {
            float max = lms_hist_LG_module[i]->GetBinCenter(lms_hist_LG_module[i]->GetMaximumBin());
            lms_hist_LG_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = lms_hist_LG_module[i]->GetFunction("gaus");
            if (fit) {
                lms_peak_LG[i] = fit->GetParameter(1); // mean
                lms_rms_LG[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                lms_hist_LG_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_LMS_LG_module_" + std::to_string(LG_module_id[i]) + ".png").c_str());
                delete c;
            }
            else {
                lms_peak_LG[i] = 0.1;
                lms_rms_LG[i] = 1e5;
            }
        } else {
            lms_peak_LG[i] = 0.1;
            lms_rms_LG[i] = 1e5;
        }
    }

    float lms_height_peak[1156], lms_height_rms[1156];
    for (int i = 0; i < 1156; i++) {
        if (lms_height_hist_module[i]->GetEntries() > 0) {
            float max = lms_height_hist_module[i]->GetBinCenter(lms_height_hist_module[i]->GetMaximumBin());
            lms_height_hist_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = lms_height_hist_module[i]->GetFunction("gaus");
            if (fit) {
                lms_height_peak[i] = fit->GetParameter(1); // mean
                lms_height_rms[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                lms_height_hist_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_LMS_height_module_" + std::to_string(i+1) + ".png").c_str());
                delete c;
            }
            else {
                lms_height_peak[i] = 0.1;
                lms_height_rms[i] = 1e5;
            }
        } else {
            lms_height_peak[i] = 0.1;
            lms_height_rms[i] = 1e5;
        }
    }

    float lms_height_peak_LG[LG_num], lms_height_rms_LG[LG_num];
    for (int i = 0; i < LG_num; i++) {
        if (lms_height_hist_LG_module[i]->GetEntries() > 0) {
            float max = lms_height_hist_LG_module[i]->GetBinCenter(lms_height_hist_LG_module[i]->GetMaximumBin());
            lms_height_hist_LG_module[i]->Fit("gaus", "Q", "r", max*0.7, max*1.5);
            TF1 *fit = lms_height_hist_LG_module[i]->GetFunction("gaus");
            if (fit) {
                lms_height_peak_LG[i] = fit->GetParameter(1); // mean
                lms_height_rms_LG[i] = fit->GetParameter(2);  // sigma
                TCanvas *c = new TCanvas();
                lms_height_hist_LG_module[i]->Draw();
                fit->Draw("same");
                c->SaveAs(("./fit_canvas/fit_LMS_height_LG_module_" + std::to_string(LG_module_id[i]) + ".png").c_str());
                delete c;
            }
            else {
                lms_height_peak_LG[i] = 0.1;
                lms_height_rms_LG[i] = 1e5;
            }
        } else {
            lms_height_peak_LG[i] = 0.1;
            lms_height_rms_LG[i] = 1e5;
        }
    }

    TH1F *peak_module = new TH1F("peak_module", "Peak Integral by Module", 100, 0, 500);
    TH1F *rms_module = new TH1F("rms_module", "RMS of Peak Integral by Module", 100, 0, 400);
    for (int i = 0; i < 1156; i++) {
        peak_module->Fill(peak[i]);
        rms_module->Fill(rms[i]);
    }
    for(int i = 0; i < LG_num; i++) {
        peak_module->Fill(peak_LG[i]);
        rms_module->Fill(rms_LG[i]);
    }
    peak_module->Write();
    rms_module->Write();

    std::ofstream csv_out(Form("cosmic_peak_%d.dat", run_number));
    csv_out << "ModuleID  PeakIntegral  RMS\n";
    for (int i = 0; i < 1156; i++) {
        csv_out << "W" << (i+1) << "  " << peak[i] << "  " << rms[i] << "\n";
    }
    for (int i = 0; i < LG_num; i++) {
        csv_out << "G" << LG_module_id[i] << "  " << peak_LG[i] << "  " << rms_LG[i] << "\n";
    }
    csv_out.close();

    std::ofstream rate_out(Form("cosmic_eventNum_%d.dat", run_number));
    rate_out << "ModuleID  EventCount\n";
    for (int i = 0; i < 1156; i++) {
        rate_out << "W" << (i+1) << "  " << event_num_module[i+1000+1] << "\n";
    }
    for (int i = 0; i < LG_num; i++) {
        rate_out << "G" << LG_module_id[i] << "  " << event_num_module[LG_module_id[i]] << "\n";
    }
    rate_out.close();

    // ── JSON output ───────────────────────────────────────────────────────
    if (run_number > 0) {
        // Build the new entry string for each module
        auto make_entry = [&](int i) -> std::string {
            char buf[512];
            std::snprintf(buf, sizeof(buf),
                "{\"run\": %d, \"peak_height_mean\": %g"
                ", \"peak_height_sigma\": %g"
                ", \"peak_height_diff\": %g"
                ", \"peak_integral_mean\": %g"
                ", \"peak_integral_sigma\": %g"
                ", \"peak_integral_diff\": %g"
                ", \"count\": %d}",
                run_number,
                peak_height[i], rms_height[i], peak_height[i] - 35.,
                peak[i], rms[i], peak[i] - 250.,
                event_num_module[i+1000+1]);
            return std::string(buf);
        };
        auto make_entry_LG = [&](int i) -> std::string {
            char buf[512];
            std::snprintf(buf, sizeof(buf),
                "{\"run\": %d, \"peak_height_mean\": %g"
                ", \"peak_height_sigma\": %g"
                ", \"peak_height_diff\": %g"
                ", \"peak_integral_mean\": %g"
                ", \"peak_integral_sigma\": %g"
                ", \"peak_integral_diff\": %g"
                ", \"count\": %d}",
                run_number,
                peakHeight_LG[i], rms_height_LG[i], peakHeight_LG[i] - 35.,
                peak_LG[i], rms_LG[i], peak_LG[i] - 250.,
                event_num_module[LG_module_id[i]]);
            return std::string(buf);
        };

        if (!in_json.empty()) {
            // ── Append mode: read existing JSON and insert new entry ──────
            std::ifstream fin(in_json);
            if (!fin) {
                std::cerr << "Cannot open input JSON: " << in_json << "\n";
            } else {
                std::vector<std::string> lines;
                std::string line;
                while (std::getline(fin, line)) lines.push_back(line);
                fin.close();

                int mod_idx = 0;
                for (auto &l : lines) {
                    // Lines with module entries contain "}]" near the end
                    auto pos = l.rfind("}");
                    if (pos != std::string::npos && pos + 1 < l.size() && l[pos+1] == ']') {
                        // insert ", {new_entry}" before the closing "]"
                        std::string new_entry = ", " + make_entry(mod_idx);
                        l.insert(pos + 1, new_entry);
                        mod_idx++;
                    }
                    if (mod_idx >= 1156){
                        if (l.rfind("}]") != std::string::npos) {
                            std::string new_entry = ", " + make_entry_LG(mod_idx - 1156);
                            l.insert(l.rfind("}]"), new_entry);
                        }
                        mod_idx++;
                    }
                }

                std::ofstream fout(in_json);
                for (auto &l : lines) fout << l << "\n";
                fout.close();
                std::cerr << "Appended run " << run_number << " to " << in_json << "\n";
            }
        } else {
            // ── Create mode: write new JSON ───────────────────────────────
            std::string out_path = Form("cosmic_modules_run%d.json", run_number);
            std::ofstream json_out(out_path);
            json_out << "{\n";
            for (int i = 0; i < 1156; i++) {
                json_out << "  \"W" << (i+1) << "\": [" << make_entry(i) << "]";
                if (i < 1155) json_out << ",";
                json_out << "\n";
            }
            for (int i = 0; i < LG_num; i++) {
                json_out << "  \"G" << LG_module_id[i] << "\": [" << make_entry_LG(i) << "]";
                if (i < LG_num - 1) json_out << ",";
                json_out << "\n";
            }
            json_out << "}\n";
            json_out.close();
            std::cerr << "JSON written to " << out_path << "\n";
        }
    }

    //json output for lms results
    std::ofstream lms_json_out(Form("lms_run_%d.json", run_number));
    lms_json_out << "{\n";
    for (int i = 0; i < 1156; i++) {
        lms_json_out << "  \"W" << (i+1) << "\": {\"run\": " << run_number << ", \"lms_peak\": " << lms_peak[i] << ", \"lms_rms\": " << lms_rms[i] << "}";
        if (i < 1155) lms_json_out << ",";
        lms_json_out << "\n";
    }
    for (int i = 0; i < LG_num; i++) {
        lms_json_out << "  \"G" << LG_module_id[i] << "\": {\"run\": " << run_number << ", \"lms_peak\": " << lms_peak_LG[i] << ", \"lms_rms\": " << lms_rms_LG[i] << "}";
        if (i < LG_num - 1) lms_json_out << ",";
        lms_json_out << "\n";
    }
    lms_json_out << "}\n";
    lms_json_out.close();

    outfile.Close();

    return 0;
}