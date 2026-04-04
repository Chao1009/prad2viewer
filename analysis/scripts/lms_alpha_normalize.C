//============================================================================
// lms_alpha_normalize.C — LMS / Alpha normalization using prad2dec decoder
//
// Scans a run's EVIO files (prad_{run}.evio.00000-99999), selects events by
// trigger_bits (LMS = bit 24, Alpha = bit 25), and normalizes every HyCal
// module's LMS response using the most recent Alpha reference on LMS 1/2/3.
//
//   norm_i = integral_i * mean( alpha_ref_j / lms_ref_j )     j = 1,2,3
//
// Compile with ACLiC after loading rootlogon:
//   root -l rootlogon.C
//   .x lms_alpha_normalize.C+("/data/prad/run1234", 1234)
//
// Or from the command line (from the build directory):
//   root -l rootlogon.C 'lms_alpha_normalize.C+("/data/prad/run1234", 1234)'
//============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "WaveAnalyzer.h"
#include "Fadc250Data.h"

#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>

// ── trigger bits (database/trigger_bits.json) ───────────────────────────
static constexpr uint32_t TBIT_LMS   = (1u << 24);
static constexpr uint32_t TBIT_ALPHA = (1u << 25);

// ── LMS reference channels ─────────────────────────────────────────────
static constexpr int N_LMS_REF = 3;   // LMS1, LMS2, LMS3

// ── channel bookkeeping ────────────────────────────────────────────────
struct ChInfo {
    std::string name;
    int  idx       = -1;   // 0..N-1 for HyCal modules, -1 otherwise
    bool is_lms    = false;
    int  lms_id    = -1;   // 0/1/2 for LMS1/2/3
};

static std::map<int, ChInfo> gCh;   // key = crate*10000 + slot*100 + ch
static int gNMod = 0;               // number of HyCal modules

static int packAddr(int c, int s, int ch) { return c * 10000 + s * 100 + ch; }

// Load database/daq_map.json — maps every (crate,slot,channel) to a name.
// W* / G* names are HyCal modules; LMS1/2/3 are reference channels.
static void loadDaqMap(const char *path)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open daq_map " << path << "\n";
        return;
    }
    auto j = nlohmann::json::parse(f);

    int idx = 0;
    for (auto &e : j) {
        std::string nm = e["name"].get<std::string>();
        int cr = e["crate"].get<int>();
        int sl = e["slot"].get<int>();
        int ch = e["channel"].get<int>();

        ChInfo ci;
        ci.name = nm;
        if      (nm == "LMS1") { ci.is_lms = true; ci.lms_id = 0; }
        else if (nm == "LMS2") { ci.is_lms = true; ci.lms_id = 1; }
        else if (nm == "LMS3") { ci.is_lms = true; ci.lms_id = 2; }
        else if (nm[0] == 'W' || nm[0] == 'G') { ci.idx = idx++; }

        gCh[packAddr(cr, sl, ch)] = ci;
    }
    gNMod = idx;
}

// ── file discovery ─────────────────────────────────────────────────────
static std::vector<std::string> discoverFiles(const char *dir, int run)
{
    std::vector<std::string> out;
    TString prefix = Form("prad_%d.evio.", run);

    TSystemDirectory d("", dir);
    TList *list = d.GetListOfFiles();
    if (!list) return out;

    TIter next(list);
    TObject *obj;
    while ((obj = next())) {
        TString name = obj->GetName();
        if (!name.BeginsWith(prefix)) continue;
        TString sfx = name(prefix.Length(), name.Length() - prefix.Length());
        bool ok = sfx.Length() > 0;
        for (int i = 0; i < sfx.Length(); i++)
            if (!isdigit(sfx[i])) { ok = false; break; }
        if (!ok) continue;
        out.push_back(std::string(dir) + "/" + name.Data());
    }
    std::sort(out.begin(), out.end());
    return out;
}

// ── main ───────────────────────────────────────────────────────────────
void lms_alpha_normalize(const char *data_dir, int run_number,
                         const char *daq_cfg_path = nullptr,
                         const char *daq_map_path = nullptr)
{
    // --- locate database directory ---
    TString dbDir = gSystem->Getenv("PRAD2_DATABASE_DIR");
    if (dbDir.IsNull()) dbDir = "database";

    TString cfgFile = daq_cfg_path ? daq_cfg_path
                                   : Form("%s/daq_config.json", dbDir.Data());
    TString mapFile = daq_map_path ? daq_map_path
                                   : Form("%s/daq_map.json", dbDir.Data());

    printf("============================================\n");
    printf(" PRad LMS / Alpha Normalization\n");
    printf(" Run       : %d\n", run_number);
    printf(" Data      : %s\n", data_dir);
    printf(" DAQ cfg   : %s\n", cfgFile.Data());
    printf(" DAQ map   : %s\n", mapFile.Data());
    printf("============================================\n");

    // --- load configs ---
    evc::DaqConfig cfg;
    printf("Loading DAQ config: %s\n", cfgFile.Data());
    if (!evc::load_daq_config(cfgFile.Data(), cfg)) {
        std::cerr << "FATAL: failed to load " << cfgFile << "\n";
        std::cerr << "  Without daq_config.json, physics event tags (0x82, 0xFE, ...)\n"
                  << "  are not recognized — GetNEvents() will always return 0.\n"
                  << "  Set PRAD2_DATABASE_DIR or pass the path explicitly.\n";
        return;
    }
    printf("  Physics tags  :");
    for (auto t : cfg.physics_tags) printf(" 0x%X", t);
    printf("\n");
    printf("  ADC format    : %s\n", cfg.adc_format.c_str());
    printf("  ROC crates    : %zu\n", cfg.roc_tags.size());
    printf("  TI master tag : 0x%X\n", cfg.ti_master_tag);

    printf("Loading DAQ map: %s\n", mapFile.Data());
    loadDaqMap(mapFile.Data());
    if (gNMod == 0) {
        std::cerr << "FATAL: no HyCal modules found in daq_map\n";
        return;
    }
    printf("  HyCal modules : %d\n", gNMod);
    printf("  LMS refs      :");
    for (auto &[addr, ci] : gCh)
        if (ci.is_lms) printf(" %s(crate %d)", ci.name.c_str(), addr / 10000);
    printf("\n");

    // build ROC tag -> crate lookup from daq_config
    std::map<uint32_t, int> tagToCrate;
    for (auto &re : cfg.roc_tags) {
        tagToCrate[re.tag] = re.crate;
        printf("  ROC 0x%02X -> crate %d (%s)\n", re.tag, re.crate, re.name.c_str());
    }

    // --- discover files ---
    auto files = discoverFiles(data_dir, run_number);
    if (files.empty()) {
        fprintf(stderr, "No prad_%d.evio.* files in %s\n", run_number, data_dir);
        return;
    }
    printf("EVIO files : %zu\n\n", files.size());

    // --- book histograms ---
    TH1D *hNormDist = new TH1D("hNormDist",
        "Avg Normalized LMS per Module;Normalized LMS (a.u.);Modules",
        400, 0, 4.0);
    TH1D *hNormMap = new TH1D("hNormMap",
        "Normalized LMS vs Module;Module index;Avg Norm LMS",
        gNMod, 0, gNMod);

    TGraph *grAlpha[N_LMS_REF], *grLMS[N_LMS_REF], *grRatio[N_LMS_REF];
    for (int j = 0; j < N_LMS_REF; j++) {
        grAlpha[j] = new TGraph();
        grAlpha[j]->SetName(Form("grAlpha%d", j + 1));
        grAlpha[j]->SetTitle(Form("Alpha Ref LMS%d;Event #;Integral", j + 1));
        grAlpha[j]->SetMarkerStyle(7);
        grAlpha[j]->SetMarkerColor(kRed + j);

        grLMS[j] = new TGraph();
        grLMS[j]->SetName(Form("grLMS%d", j + 1));
        grLMS[j]->SetTitle(Form("LMS Ref LMS%d;Event #;Integral", j + 1));
        grLMS[j]->SetMarkerStyle(7);
        grLMS[j]->SetMarkerColor(kBlue + j);

        grRatio[j] = new TGraph();
        grRatio[j]->SetName(Form("grRatio%d", j + 1));
        grRatio[j]->SetTitle(Form("Alpha/LMS Ratio LMS%d;LMS Event #;Ratio", j + 1));
        grRatio[j]->SetMarkerStyle(7);
        grRatio[j]->SetMarkerColor(kGreen + 2 + j);
    }

    // --- accumulators ---
    std::vector<double> modSum(gNMod, 0.0);
    std::vector<int>    modCnt(gNMod, 0);

    double alphaRef[N_LMS_REF] = {};
    bool   alphaValid = false;
    int    nTotal = 0, nAlpha = 0, nLMS = 0, nSkipped = 0;

    // --- decoder objects ---
    evc::EvChannel reader;
    reader.SetConfig(cfg);

    fdec::EventData evt;
    fdec::WaveAnalyzer wave;

    // --- event loop ---------------------------------------------------------
    for (size_t fi = 0; fi < files.size(); fi++) {
        printf("[%zu/%zu] %s\n", fi + 1, files.size(), files[fi].c_str());

        if (reader.Open(files[fi]) != evc::status::success) {
            printf("  WARNING: cannot open, skipping\n");
            continue;
        }

        while (reader.Read() == evc::status::success) {
            if (!reader.Scan()) continue;
            if (reader.GetEventType() != evc::EventType::Physics) continue;

            for (int ie = 0; ie < reader.GetNEvents(); ie++) {
                evt.clear();
                if (!reader.DecodeEvent(ie, evt)) continue;
                nTotal++;

                bool isLMS   = (evt.info.trigger_bits & TBIT_LMS)   != 0;
                bool isAlpha = (evt.info.trigger_bits & TBIT_ALPHA) != 0;
                if (!isLMS && !isAlpha) continue;

                // -- compute waveform integrals for all channels ----
                // store as:  packed_addr -> integral
                std::map<int, float> integrals;

                for (int r = 0; r < fdec::MAX_ROCS; r++) {
                    auto &roc = evt.rocs[r];
                    if (!roc.present) continue;
                    auto ct = tagToCrate.find(roc.tag);
                    if (ct == tagToCrate.end()) continue;
                    int crate = ct->second;

                    for (int s = 0; s < fdec::MAX_SLOTS; s++) {
                        if (!roc.slots[s].present) continue;
                        for (int c = 0; c < fdec::MAX_CHANNELS; c++) {
                            auto &cd = roc.slots[s].channels[c];
                            if (cd.nsamples == 0) continue;

                            fdec::WaveResult wres;
                            wave.Analyze(cd.samples, cd.nsamples, wres);

                            float integral = 0.f;
                            for (int k = 0; k < cd.nsamples; k++)
                                integral += cd.samples[k] - wres.ped.mean;

                            integrals[packAddr(crate, s, c)] = integral;
                        }
                    }
                }

                // -- extract LMS reference values from this event ----
                float lmsVal[N_LMS_REF] = {};
                for (auto &[addr, integ] : integrals) {
                    auto it = gCh.find(addr);
                    if (it == gCh.end()) continue;
                    if (it->second.is_lms)
                        lmsVal[it->second.lms_id] = integ;
                }

                // -- Alpha trigger -----------------------------------
                if (isAlpha) {
                    bool good = true;
                    for (int j = 0; j < N_LMS_REF; j++) {
                        if (lmsVal[j] > 0) {
                            alphaRef[j] = lmsVal[j];
                            grAlpha[j]->SetPoint(grAlpha[j]->GetN(),
                                                 nTotal, lmsVal[j]);
                        } else {
                            good = false;
                        }
                    }
                    if (good) alphaValid = true;
                    nAlpha++;
                }

                // -- LMS trigger ------------------------------------
                if (isLMS) {
                    if (!alphaValid) { nSkipped++; nLMS++; continue; }

                    bool refOK = true;
                    for (int j = 0; j < N_LMS_REF; j++) {
                        grLMS[j]->SetPoint(grLMS[j]->GetN(),
                                           nTotal, lmsVal[j]);
                        if (lmsVal[j] <= 0) refOK = false;
                    }
                    if (!refOK) { nLMS++; continue; }

                    // average normalization ratio across 3 references
                    double avgRatio = 0.0;
                    int nRef = 0;
                    for (int j = 0; j < N_LMS_REF; j++) {
                        if (alphaRef[j] > 0 && lmsVal[j] > 0) {
                            double r = alphaRef[j] / lmsVal[j];
                            grRatio[j]->SetPoint(grRatio[j]->GetN(),
                                                 nLMS, r);
                            avgRatio += r;
                            nRef++;
                        }
                    }
                    if (nRef == 0) { nLMS++; continue; }
                    avgRatio /= nRef;

                    // normalize every HyCal module
                    for (auto &[addr, integ] : integrals) {
                        auto it = gCh.find(addr);
                        if (it == gCh.end()) continue;
                        int mi = it->second.idx;
                        if (mi < 0 || mi >= gNMod) continue;
                        if (integ <= 0) continue;

                        double norm = integ * avgRatio;
                        modSum[mi] += norm;
                        modCnt[mi]++;
                    }
                    nLMS++;
                }

                if (nTotal % 50000 == 0)
                    printf("  %d events  (Alpha %d  LMS %d)\r",
                           nTotal, nAlpha, nLMS);
            }
        }
        reader.Close();
    }

    // --- summary ------------------------------------------------------------
    printf("\n\n=== Run %d Summary ===\n", run_number);
    printf("Total events  : %d\n", nTotal);
    printf("Alpha events  : %d\n", nAlpha);
    printf("LMS events    : %d  (skipped %d — no alpha ref yet)\n",
           nLMS, nSkipped);

    // --- fill per-module histograms -----------------------------------------
    int nGood = 0;
    for (int i = 0; i < gNMod; i++) {
        if (modCnt[i] == 0) continue;
        double avg = modSum[i] / modCnt[i];
        hNormDist->Fill(avg);
        hNormMap->SetBinContent(i + 1, avg);
        nGood++;
    }
    printf("Modules with data: %d / %d\n\n", nGood, gNMod);

    // --- save ROOT output ---------------------------------------------------
    TString outname = Form("lms_alpha_run%d.root", run_number);
    TFile *fout = new TFile(outname, "RECREATE");
    hNormDist->Write();
    hNormMap->Write();
    for (int j = 0; j < N_LMS_REF; j++) {
        grAlpha[j]->Write();
        grLMS[j]->Write();
        grRatio[j]->Write();
    }
    fout->Close();
    printf("Saved %s\n", outname.Data());

    // --- draw summary canvas ------------------------------------------------
    gStyle->SetOptStat(111);

    TCanvas *c1 = new TCanvas("c1", "LMS/Alpha Overview", 1400, 900);
    c1->Divide(3, 2);

    c1->cd(1);
    hNormDist->SetLineColor(kBlack);
    hNormDist->Draw();

    c1->cd(2);
    hNormMap->SetMarkerStyle(6);
    hNormMap->Draw("P");

    c1->cd(3);
    if (grAlpha[0]->GetN() > 0) {
        grAlpha[0]->Draw("AP");
        for (int j = 1; j < N_LMS_REF; j++) grAlpha[j]->Draw("P SAME");
        auto *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        for (int j = 0; j < N_LMS_REF; j++)
            leg->AddEntry(grAlpha[j], Form("LMS%d", j + 1), "p");
        leg->Draw();
    }

    c1->cd(4);
    if (grLMS[0]->GetN() > 0) {
        grLMS[0]->Draw("AP");
        for (int j = 1; j < N_LMS_REF; j++) grLMS[j]->Draw("P SAME");
        auto *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        for (int j = 0; j < N_LMS_REF; j++)
            leg->AddEntry(grLMS[j], Form("LMS%d", j + 1), "p");
        leg->Draw();
    }

    c1->cd(5);
    if (grRatio[0]->GetN() > 0) {
        grRatio[0]->Draw("AP");
        for (int j = 1; j < N_LMS_REF; j++) grRatio[j]->Draw("P SAME");
        auto *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        for (int j = 0; j < N_LMS_REF; j++)
            leg->AddEntry(grRatio[j], Form("LMS%d", j + 1), "p");
        leg->Draw();
    }

    c1->Update();
    c1->SaveAs(Form("lms_alpha_run%d.png", run_number));

    printf("Done.\n");
}
