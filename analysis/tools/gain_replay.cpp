// gain_replay.cpp — replay tool for gain monitoring data
// Usage: gain_replay <input evio dir> [-o output.dat] [-D daq_config.json] [-n N]

#include "Replay.h"
#include "DaqConfig.h"
#include "HyCalSystem.h"
#include "PhysicsTools.h"
#include "PulseTemplateStore.h"

#include <TF1.h>
#include <TMath.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

using json = nlohmann::json;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

static std::vector<std::string> getFilesInDir(const std::string &dir_path)
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
    std::string input, output, daq_config;
    int max_events = -1;

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:n:D:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'D': daq_config = optarg; break;
        }
    }
    if (optind < argc) input = argv[optind];
    else {
        std::cerr << "Usage: gain_replay <input evio dir> [-o output.dat] [-D daq_config.json] [-n N]\n";
        return 1;
    }

    std::vector<std::string> evio_files = getFilesInDir(input);
    if (evio_files.empty()) {
        std::cerr << "No EVIO files found in: " << input << "\n";
        return 1;
    }

    if (input.empty()) {
        std::cerr << "Usage: gain_replay <input evio dir> [-o output.dat] [-D daq_config.json] [-n N]\n";
        return 1;
    }

    evc::DaqConfig daq_cfg;
    evc::load_daq_config(daq_config, daq_cfg);

    analysis::Replay replay;
    if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
    replay.LoadDaqMap(db_dir + "/hycal_daq_map.json");
    std::cerr << "Using DAQ map: " << db_dir + "/hycal_daq_map.json" << "\n";
    
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_modules.json", db_dir + "/hycal_daq_map.json");
    analysis::PhysicsTools physics(hycal);

    // build ROC tag → crate index mapping from DAQ config JSON
    std::unordered_map<int, int> roc_to_crate;
    if (!daq_config.empty()) {
        std::cout << "Loading DAQ config from " << daq_config << "\n";
        std::ifstream dcf(daq_config);
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

    if(output.empty()) output = "gain_replay_output.root";
    TFile *outfile = TFile::Open(output.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Cannot create output ROOT file\n";
        return 1;
    }

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);

    auto event = std::make_unique<fdec::EventData>();
    fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
    fdec::PulseTemplateStore template_store;
    if (daq_cfg.wave_cfg.nnls_deconv.enabled
        && !daq_cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        template_store.LoadFromFile(
            db_dir + "/" + daq_cfg.wave_cfg.nnls_deconv.template_file,
            daq_cfg.wave_cfg);
    }
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;
    int total = 0;

    for (const auto &input_evio : evio_files) {
        if (ch.OpenAuto(input_evio) != evc::status::success) {
            std::cerr << "Replay: cannot open " << input_evio << "\n";
            return 1;
        }

        while (ch.Read() == evc::status::success) {
            if (!ch.Scan()) continue;
            if (ch.GetEventType() != evc::EventType::Physics) continue;

            for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
                event->clear();
                if (!ch.DecodeEvent(ie, *event)) continue;
                if (max_events > 0 && total >= max_events) goto done;

                uint32_t trigger_bits = event->info.trigger_bits;
                uint64_t timestamp    = event->info.timestamp;

                static constexpr uint32_t TBIT_lms = (1u << 24);
                static constexpr uint32_t TBIT_alpha = (1u << 25);
                if ( !(trigger_bits & TBIT_lms) && !(trigger_bits & TBIT_alpha) ) continue;

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
                            int mod_id = replay.moduleID(crate, s, c);
                            if(mod_id < 0){
                                std::string mod_name = replay.moduleName(crate, s, c);
                                if(mod_name.empty()) continue;
                                if(mod_name[0] == 'L'){
                                    if(mod_name.length() != 4) continue;
                                    int lms_id;
                                    if(mod_name[3] == 'P') lms_id = 0;
                                    else lms_id = mod_name[3] - '0';
                                    ana.SetChannelKey(roc.tag, s, c);
                                    ana.Analyze(cd.samples, cd.nsamples, wres);
                                    //if(wres.npeaks != 1) continue;
                                    int idx = 0;
                                    for (int k = 0; k < wres.npeaks; ++k)
                                        if (wres.peaks[k].height > wres.peaks[idx].height) idx = k;
                                    float peak_height = wres.peaks[idx].height;
                                    float peak_integral = wres.peaks[idx].integral;
                                    if(trigger_bits & TBIT_lms) {
                                        physics.Fill_lmsCH_lmsHeight(lms_id, peak_height);
                                        physics.Fill_lmsCH_lmsIntegral(lms_id, peak_integral);
                                    }
                                    if(trigger_bits & TBIT_alpha) {
                                        physics.Fill_lmsCH_alphaHeight(lms_id, peak_height);
                                        physics.Fill_lmsCH_alphaIntegral(lms_id, peak_integral);
                                    }
                                    total++;
                                }
                            }else {
                                if(!(trigger_bits & TBIT_lms)) continue;

                                int mod_id = replay.moduleID(crate, s, c);
                                if (mod_id < 0) continue;
                                ana.SetChannelKey(roc.tag, s, c);
                                ana.Analyze(cd.samples, cd.nsamples, wres);
                                if(wres.npeaks != 1) continue;;
                                float peak_height = wres.peaks[0].height;
                                float peak_integral = wres.peaks[0].integral;
                                physics.Fill_modCH_lmsHeight(mod_id, peak_height);
                                physics.Fill_modCH_lmsIntegral(mod_id, peak_integral);
                                total++;
                            }
                        }
                    }
                }
                std::cerr << "Processed event " << total << " / " << max_events << "\r" << std::flush;
            }
        }
        ch.Close();
    }

done:
    // analyze and save results
    // fit LMS or alpha signal peaks and output to a .dat file
    //Format:
    //        LMS1  lms1_peak lms1_sigma lms1_kai2/ndf alpha1_peak alpha1_sigma alpha1_kai2/ndf
    //         ...
    //         W1   lms_peak lms_sigma lms_kai2/ndf    g1 g2 g3
    //         ...
    //       gi = lms_peak * alpha_i_peak / lms_i_peak
    std::ofstream out("lms_alpha_peaks.dat");
    if (!out.is_open()) {
        std::cerr << "Cannot open output file lms_alpha_peaks.dat\n";
        return 1;
    }

    //variables to store fit results
    const int nmod = hycal.module_count();
    float lms_ref[4] = {0}, alpha_ref[4] = {0};
    float lms_sigma[4] = {0}, alpha_sigma[4] = {0};
    float lms_chi2[4] = {0}, alpha_chi2[4] = {0};
    float mod_lms[nmod] = {0};
    float mod_lms_sigma[nmod] = {0};
    float mod_lms_chi2[nmod] = {0};
    float g[nmod][4] = {{0}}; // gain factors for each module and alpha peak

    int n_lms = 0;
    for (int i = 0; i < 4; ++i) {
        TH1F *h_lms = physics.Get_lmsCH_lmsHeightHist(i);
        if (h_lms == nullptr || h_lms->GetEntries() < 10) continue;
        {
            double peak0 = h_lms->GetBinCenter(h_lms->GetMaximumBin());
            double rms0  = h_lms->GetRMS();
            double lo = peak0 - 2.0 * rms0, hi = peak0 + 2.0 * rms0;
            TF1 f_gaus("f_lms", "gaus", lo, hi);
            f_gaus.SetParameters(h_lms->GetMaximum(), peak0, rms0);
            h_lms->Fit(&f_gaus, "RQ0");
            lms_ref[i]   = f_gaus.GetParameter(1);
            lms_sigma[i] = f_gaus.GetParameter(2);
            lms_chi2[i]  = (f_gaus.GetNDF() > 0) ? f_gaus.GetChisquare() / f_gaus.GetNDF() : 0;
        }

        TH1F *h_alpha = physics.Get_lmsCH_alphaHeightHist(i);
        if (h_alpha == nullptr || h_alpha->GetEntries() < 10) continue;
        {
            double peak0 = h_alpha->GetBinCenter(h_alpha->GetMaximumBin());
            double rms0  = h_alpha->GetRMS();
            double lo = peak0 - 2.0 * rms0, hi = peak0 + 2.0 * rms0;
            TF1 f_gaus("f_alpha", "gaus", lo, hi);
            f_gaus.SetParameters(h_alpha->GetMaximum(), peak0, rms0);
            h_alpha->Fit(&f_gaus, "RQ0");
            alpha_ref[i]   = f_gaus.GetParameter(1);
            alpha_sigma[i] = f_gaus.GetParameter(2);
            alpha_chi2[i]  = (f_gaus.GetNDF() > 0) ? f_gaus.GetChisquare() / f_gaus.GetNDF() : 0;
        }

        if (lms_ref[i] > 0) n_lms++;
    }
    out << std::left;
    out << std::setw(6) << "Name" << std::setw(12) << "lms_peak" << std::setw(12) << "lms_sigma" << std::setw(12) << "lms_chi2/ndf" 
        << std::setw(16) << "alpha_peak (g1)" << std::setw(16) << "alpha_sigma (g2)" << std::setw(16) << "alpha_chi2/ndf (g3)" << "\n";
    for(int i = 1; i <= 3; i++){
        out << std::setw(6) << ("LMS" + std::to_string(i))
            << std::setw(12) << std::fixed << std::setprecision(3) << lms_ref[i]
            << std::setw(12) << std::fixed << std::setprecision(3) << lms_sigma[i]
            << std::setw(12) << std::fixed << std::setprecision(3) << lms_chi2[i]
            << std::setw(16) << std::fixed << std::setprecision(3) << alpha_ref[i]
            << std::setw(16) << std::fixed << std::setprecision(3) << alpha_sigma[i]
            << std::setw(16) << std::fixed << std::setprecision(3) << alpha_chi2[i]
            << "\n";
    }
    int n_mod = 0;
    for (int i = 0; i < hycal.module_count(); ++i) {
        if (!hycal.module(i).is_hycal()) continue;
        if (hycal.module(i).name[0] != 'W') continue;

        TH1F *h_lms = physics.Get_modCH_lmsHeightHist(hycal.module(i).id);
        if (h_lms == nullptr) continue;

        {
            double peak0 = h_lms->GetBinCenter(h_lms->GetMaximumBin());
            double rms0  = h_lms->GetRMS();
            double lo = peak0 - 2.0 * rms0, hi = peak0 + 2.0 * rms0;
            TF1 f_gaus("f_mod", "gaus", lo, hi);
            f_gaus.SetParameters(h_lms->GetMaximum(), peak0, rms0);
            h_lms->Fit(&f_gaus, "RQ");
            mod_lms[i]       = f_gaus.GetParameter(1);
            mod_lms_sigma[i] = f_gaus.GetParameter(2);
            mod_lms_chi2[i]  = (f_gaus.GetNDF() > 0) ? f_gaus.GetChisquare() / f_gaus.GetNDF() : 0;
        }

        for (int j = 1; j <= 3; ++j) {
            g[i][j] = (lms_ref[j] > 0 && mod_lms[i] > 0 && alpha_ref[j] > 0) ? 
                        mod_lms[i] * alpha_ref[j] / lms_ref[j] : 0;
        }
        out << std::setw(6) << hycal.module(i).name
            << std::setw(12) << std::fixed << std::setprecision(3) << mod_lms[i]
            << std::setw(12) << std::fixed << std::setprecision(3) << mod_lms_sigma[i]
            << std::setw(12) << std::fixed << std::setprecision(3) << mod_lms_chi2[i]
            << std::setw(16) << std::fixed << std::setprecision(3) << g[i][1]
            << std::setw(16) << std::fixed << std::setprecision(3) << g[i][2]
            << std::setw(16) << std::fixed << std::setprecision(3) << g[i][3]
            << "\n";
        n_mod++;
    }
    out.close();
    std::cerr << "LMS and alpha peak analysis completed for " << n_lms << " LMS channels and " << n_mod << " modules. Results saved to lms_alpha_peaks.dat\n";

    outfile->cd();
    outfile->mkdir("lms");
    outfile->cd("lms");
    for (int i = 0; i < 4; ++i) {
        if (physics.Get_lmsCH_lmsHeightHist(i)) physics.Get_lmsCH_lmsHeightHist(i)->Write();
        if (physics.Get_lmsCH_lmsIntegralHist(i)) physics.Get_lmsCH_lmsIntegralHist(i)->Write();
        if (physics.Get_lmsCH_alphaHeightHist(i)) physics.Get_lmsCH_alphaHeightHist(i)->Write();
        if (physics.Get_lmsCH_alphaIntegralHist(i)) physics.Get_lmsCH_alphaIntegralHist(i)->Write();
    }
    outfile->mkdir("modules");
    outfile->cd("modules");
    for (int i = 0; i < hycal.module_count(); ++i) {
        if (physics.Get_modCH_lmsHeightHist(hycal.module(i).id)) physics.Get_modCH_lmsHeightHist(hycal.module(i).id)->Write();
    }
    outfile->Close();
}