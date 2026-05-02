//=============================================================================
// analysis_example.cpp some examples for offline physics analysis
//
// Usage: analysis_example <input_recon.root> [-o output.root] [-n max_events]
//
// Reads the reconstructed root files, call help functions from PhysicsTools, fills per-module energy histograms
// and moller event analysis histograms, and saves to output ROOT file.
//=============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"
#include "Fadc250Data.h"
#include "InstallPaths.h"

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TROOT.h>

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <regex>
#include <string>
#include <stdexcept>
#include <cctype>
#include <iomanip>

namespace fs = std::filesystem;

std::string EVIODIR = ".";
std::string OUTPUTDIR = ".";

// ── trigger bits (database/trigger_bits.json) ───────────────────────────
static constexpr uint32_t TBIT_LMS   = (1u << 24);
static constexpr uint32_t TBIT_ALPHA = (1u << 25);

// ── LMS reference channels ─────────────────────────────────────────────
static constexpr int N_LMS_REF = 3;   // LMS1, LMS2, LMS3
static constexpr int N_HYCAL_MOD = 1728;   // number of hycal modules

// ── channel bookkeeping ────────────────────────────────────────────────
struct ChInfo {
    std::string name;
    int  idx       = -1;   // 0..N-1 for HyCal modules, -1 otherwise
    bool is_lms    = false;
    int  lms_id    = -1;   // 0/1/2 for LMS1/2/3
    int  hist_idx  = -1;   // precomputed index into LMSHist / AlphaHist
};

static std::unordered_map<int, ChInfo> gCh;   // key = crate*10000 + slot*100 + ch
static int gNMod = 0;               // number of HyCal modules
std::vector<TH1F*> ObjectContainer;   // container for saving root histogram/graph

static int packAddr(int c, int s, int ch) { return c * 10000 + s * 100 + ch; }
static void loadDaqMap(const char *path);
std::vector<std::string> discoverFiles(const char *dir, const std::string run, unsigned int startFileNum, unsigned int endFileNum);
TH1F* Init1DHist(const std::string & name,  const std::string & title, const int & nbin, 
                 const double & min, const double & max, 
                 const std::string & xaxis, const std::string & yaxis, const int & color);
void WritePlotToFile(const std::string fileName);
int channelNameToIndex(const std::string &name);
std::string indexToChannelName(int index);

int main(int argc, char *argv[])
{   
    std::string run_number = "";
    unsigned int startFileNum = 0;
    unsigned int endFileNum = 100;
    int opt;
    while ((opt = getopt(argc, argv, "r:s:e:o:i:")) != -1) {
        switch (opt) {
            case 'r': run_number    = optarg; break;
            case 's': startFileNum  = std::atoi(optarg); break;
            case 'e': endFileNum    = std::atoi(optarg); break;
            case 'o': OUTPUTDIR     = optarg; break;
            case 'i': EVIODIR     = optarg; break;
        }
    }
    if (endFileNum < startFileNum) endFileNum = startFileNum;
    
    TString dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR).c_str();

    TString cfgFile = Form("%s/daq_config.json", dbDir.Data());
    TString mapFile = Form("%s/hycal_daq_map.json", dbDir.Data());

    printf("============================================\n");
    printf(" PRad LMS / Alpha Normalization\n");
    printf(" Run        : %s\n", run_number.data());
    printf(" Data Dir   : %s\n", EVIODIR.data());
    printf(" Output Dir : %s\n", OUTPUTDIR.data());
    printf(" DAQ cfg    : %s\n", cfgFile.Data());
    printf(" DAQ map    : %s\n", mapFile.Data());
    printf("============================================\n");
    
    // --- load configs ---
    evc::DaqConfig cfg;
    if (!evc::load_daq_config(cfgFile.Data(), cfg)) {
        std::cerr << "FATAL: failed to load " << cfgFile << "\n";
        return 0;
    }

    loadDaqMap(mapFile.Data());
    if (gNMod == 0) {
        std::cerr << "FATAL: no HyCal modules found in daq_map\n";
        return 0;
    }
    printf("HyCal modules: %d\n", gNMod);
    
    // build ROC tag -> crate lookup from daq_config
    std::map<uint32_t, int> tagToCrate;
    for (auto &re : cfg.roc_tags)
        tagToCrate[re.tag] = re.crate;

    //populate two help arrays for the hycal modules and index
    // index → name (fastest: vector)
    std::vector<std::string> indexToName(N_LMS_REF+N_HYCAL_MOD);
    // name → index (fast lookup: unordered_map)
    std::unordered_map<std::string, int> nameToIndex;
    nameToIndex.reserve(N_LMS_REF+N_HYCAL_MOD);
    // ----------------------------------------
    // Build maps
    // ----------------------------------------
    for (int i = 0; i < N_LMS_REF+N_HYCAL_MOD; ++i) {
        std::string name = indexToChannelName(i);
        indexToName[i] = name;
        nameToIndex[name] = i;
    }

    // Pre-compute hist_idx in each ChInfo to avoid a second map lookup in the hot loop
    for (auto &[addr, ci] : gCh) {
        auto it = nameToIndex.find(ci.name);
        if (it != nameToIndex.end()) ci.hist_idx = it->second;
    }
    
    //initialize root histograms
    std::vector<TH1F*> LMSHist;
    std::vector<TH1F*> AlphaHist;
    
    for (int i=0; i<N_LMS_REF; ++i){
        AlphaHist.push_back(Init1DHist(Form("LMS%d_Alpha", i+1), 
                                       Form("LMS%d_Alpha", i+1), 
                                       900, 0, 4500, "ADC", "count", 1));
    }
    for (int i=0; i<N_LMS_REF+N_HYCAL_MOD; ++i){
        LMSHist.push_back(Init1DHist(Form("%s_LMS", indexToName[i].c_str()), 
                                     Form("%s_LMS", indexToName[i].c_str()), 
                                     900, 0, 4500, "ADC", "count", 1));
    }
    std::vector<std::string> InputFiles = discoverFiles(EVIODIR.data(), run_number, startFileNum, endFileNum);
    printf("found EVIO files : %zu\n", InputFiles.size());
    
    // --- decoder objects ---
    evc::EvChannel reader;
    reader.SetConfig(cfg);

    static fdec::EventData evt;
    fdec::WaveAnalyzer wave(cfg.wave_cfg);
    fdec::PulseTemplateStore template_store;
    if (cfg.wave_cfg.nnls_deconv.enabled
        && !cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        // dbDir was already env-aware-resolved above (PRAD2_DATABASE_DIR
        // → ../share/... → DATABASE_DIR), so reuse it instead of going
        // back to the build-time DATABASE_DIR macro.
        template_store.LoadFromFile(
            std::string(dbDir.Data()) + "/" + cfg.wave_cfg.nnls_deconv.template_file,
            cfg.wave_cfg);
    }
    wave.SetTemplateStore(&template_store);
    
    double alphaRef[N_LMS_REF] = {};
    bool   alphaValid = false;

    // Hoist outside all loops: avoids repeated heap allocation/deallocation per event.
    // unordered_map gives O(1) average inserts vs O(log n) for std::map.
    std::unordered_map<int, float> integrals;
    integrals.reserve(2048);

    for (size_t fi = 0; fi < InputFiles.size(); fi++) {
        printf("[%zu/%zu] %s\n", fi + 1, InputFiles.size(), InputFiles[fi].c_str());

        if (reader.OpenAuto(InputFiles[fi]) != evc::status::success) {
            printf("  WARNING: cannot open, skipping\n");
            continue;
        }
        
        while (reader.Read() == evc::status::success) {
            if (!reader.Scan()) continue;
            if (reader.GetEventType() != evc::EventType::Physics) continue;

            for (int ie = 0; ie < reader.GetNEvents(); ie++) {
                
                evt.clear();
                if (!reader.DecodeEvent(ie, evt)) continue;
                
                bool isLMS = (evt.info.trigger_bits & TBIT_LMS);
                bool isAlpha = (evt.info.trigger_bits & TBIT_ALPHA);

                if (!(isLMS || isAlpha)) continue;

                // -- compute waveform integrals for all channels ----
                // store as:  packed_addr -> integral
                integrals.clear();

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
                            wave.SetChannelKey(roc.tag, s, c);
                            wave.Analyze(cd.samples, cd.nsamples, wres);

                            // Find peak sample index in window [40, 60]
                            int peakSample = 40;
                            float peakVal = cd.samples[40] - wres.ped.mean;
                            int searchEnd = std::min(60, cd.nsamples - 1);
                            for (int k = 41; k <= searchEnd; ++k) {
                                float v = cd.samples[k] - wres.ped.mean;
                                if (v > peakVal) { peakVal = v; peakSample = k; }
                            }

                            // Average ±3 samples around the peak
                            int avgStart = std::max(0, peakSample - 3);
                            int avgEnd   = std::min(cd.nsamples - 1, peakSample + 3);
                            float sum = 0.f;
                            int   cnt = 0;
                            for (int k = avgStart; k <= avgEnd; ++k) {
                                sum += cd.samples[k] - wres.ped.mean;
                                ++cnt;
                            }
                            float integral = (cnt > 0) ? sum / cnt : 0.f;

                            integrals[packAddr(crate, s, c)] = integral;
                        }
                    }
                }

                const size_t nIntegrals = integrals.size();
                for (auto &[addr, integ] : integrals) {
                    auto it = gCh.find(addr);
                    if (it == gCh.end()) continue;
                    const auto &ci = it->second;
                    if (ci.hist_idx < 0) continue;

                    if (ci.is_lms && isAlpha && nIntegrals < 500){
                        assert(ci.hist_idx <= N_LMS_REF);
                        AlphaHist[ci.hist_idx]->Fill(integ);
                    }
                    if (isLMS && nIntegrals > 500){
                        assert(ci.hist_idx <= N_LMS_REF + N_HYCAL_MOD);
                        LMSHist[ci.hist_idx]->Fill(integ);
                    }
                }
            }
        }
        
        reader.Close();
    }

    //fitting the histograms
    //for (unsigned int i=0; i<AlphaHist.size(); i++) AlphaPara[i] = FitHistogramWithGaussian(AlphaHist[i], 0.1);
    //for (unsigned int i=0; i<LMSHist.size(); i++) LMSPara[i] = FitHistogramWithGaussian(LMSHist[i], 0.1);

    WritePlotToFile(Form("%s/prad_%s_LMS_file_%d_%d.root", OUTPUTDIR.c_str(), run_number.c_str(), startFileNum, endFileNum));
    
    //write gain factors to dat file
    //format: first 3 line reference channel name, alpha peak position, alpha sigma, alpha fit chi2/ndf, lms peak position, lms sigma, lms fit chi2/ndf
    //format: the rest: HyCal module name, lms peak, lms sigma, lms fit chi2/ndf, and three gain factors using 3 reference PMT
    
    /*std::ofstream outDatFile;
    outDatFile.open(Form("%s/prad_%s_LMS.dat", OUTPUTDIR.c_str(), run_number.c_str()));
    
    for (unsigned int i = 0; i<N_LMS_REF; i++){
        outDatFile<<std::setw(9)<<Form("LMS%d", i+1)
                  <<std::setw(15)<<AlphaPara[i].mean<<std::setw(15)<<AlphaPara[i].sigma<<std::setw(15)<<AlphaPara[i].chi2pndf
                  <<std::setw(15)<<LMSPara[i].mean<<std::setw(15)<<LMSPara[i].sigma<<std::setw(15)<<LMSPara[i].chi2pndf<<std::endl;
    }
    
    for (unsigned int i = N_LMS_REF; i < N_LMS_REF + N_HYCAL_MOD; i++){
        float factor[N_LMS_REF] = {0., 0., 0.};
        for (int j = 0; j<N_LMS_REF; j++){
            if (AlphaPara[j].mean > 1. && LMSPara[j].mean > 1.) 
            factor[j] = LMSPara[i].mean * AlphaPara[j].mean / LMSPara[j].mean;
        }
        outDatFile<<std::setw(9)<<indexToName[i]
                  <<std::setw(15)<<LMSPara[i].mean<<std::setw(15)<<LMSPara[i].sigma<<std::setw(15)<<LMSPara[i].chi2pndf
                  <<std::setw(15)<<factor[0]<<std::setw(15)<<factor[1]<<std::setw(15)<<factor[2]<<std::endl;
    }
    outDatFile.close();*/
    return 0;
}

// Load database/hycal_daq_map.json — maps every (crate,slot,channel) to a name.
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
std::vector<std::string> discoverFiles(const char *dir, const std::string run, unsigned int startFileNum, unsigned int endFileNum)
{
    std::vector<std::pair<int, std::string>> files_with_index;

    // Construct subdirectory path: dir/prad_<run>
    fs::path subdir = fs::path(dir) / ("prad_" + run);

    if (!fs::exists(subdir) || !fs::is_directory(subdir)) {
        throw std::runtime_error("Subdirectory not found: " + subdir.string());
    }

    // Regex to match: prad_<run>.evio.XXXXX
    std::string pattern = "prad_" + run + R"(\.evio\.(\d+))";
    std::regex re(pattern);

    for (const auto &entry : fs::directory_iterator(subdir)) {
        if (!entry.is_regular_file()) continue;

        std::string filename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_match(filename, match, re)) {
            int subrun = std::stoi(match[1]);  // extract XXXXX
            if ((unsigned int)subrun < startFileNum || (unsigned int)subrun > endFileNum) continue;
            files_with_index.emplace_back(subrun, entry.path().string());
        }
    }

    // Sort by subrun number
    std::sort(files_with_index.begin(), files_with_index.end(),
              [](const auto &a, const auto &b) {
                  return a.first < b.first;
              });

    // Extract sorted file paths
    std::vector<std::string> result;
    result.reserve(files_with_index.size());

    for (const auto &p : files_with_index) {
        result.push_back(p.second);
    }

    return result;
}

TH1F* Init1DHist(const std::string & name,  const std::string & title, const int & nbin, const double & min, const double & max, 
                 const std::string & xaxis, const std::string & yaxis, const int & color)
{
    TH1F* h = new TH1F(name.c_str(), title.c_str(), nbin, min, max);
    h->GetXaxis()->SetTitle(xaxis.c_str());
    h->GetYaxis()->SetTitle(yaxis.c_str());
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    ObjectContainer.push_back(h);
    return h;
}

void WritePlotToFile(const std::string fileName)
{
    TFile* f = new TFile(fileName.c_str(), "RECREATE");
    f->cd();
    
    for (unsigned int i = 0; i < ObjectContainer.size(); i++)
    ObjectContainer[i]->Write();

    f->Close();
    delete f;
}

int channelNameToIndex(const std::string &name)
{
    // --- LMS first ---
    if (name == "LMS1") return 0;
    if (name == "LMS2") return 1;
    if (name == "LMS3") return 2;

    if (name.size() < 2) {
        throw std::out_of_range("Invalid channel name: " + name);
    }

    const char type = name[0];
    const int num = std::stoi(name.substr(1));

    if (type == 'W') {
        if (num < 1 || num > 1156 || num == 561 || num == 562 || num == 595 || num == 596) {
            throw std::out_of_range("Invalid W channel: " + name);
        }

        int idx = num - 1;
        if (num >= 563) idx -= 2; // skip W561, W562
        if (num >= 597) idx -= 2; // skip W595, W596

        return 3 + idx; // shift by 3
    }

    if (type == 'G') {
        int local = -1;

        if (num >= 1 && num <= 186) {
            local = num - 1;
        }
        else if (num >= 205 && num <= 696) {
            int d = num - 205;
            int block = d / 30;
            int off   = d % 30;

            if (off >= 12) {
                throw std::out_of_range("Missing G channel: " + name);
            }

            local = 186 + block * 12 + off;
        }
        else if (num >= 715 && num <= 900) {
            local = 390 + (num - 715);
        }
        else {
            throw std::out_of_range("Invalid G channel: " + name);
        }

        return 3 + 1152 + local; // shift by 3
    }

    throw std::out_of_range("Invalid channel name: " + name);
}

std::string indexToChannelName(int index)
{
    if (index < 0 || index > 1730) {
        throw std::out_of_range("Index out of range");
    }

    // --- LMS first ---
    if (index == 0) return "LMS1";
    if (index == 1) return "LMS2";
    if (index == 2) return "LMS3";

    // --- W channels ---
    if (index >= 3 && index < 3 + 1152) {
        int local = index - 3;

        int num = local + 1;

        if (local >= 560) num += 2; // reinsert W561, W562
        if (local >= 592) num += 2; // reinsert W595, W596

        return "W" + std::to_string(num);
    }

    // --- G channels ---
    int local = index - (3 + 1152);

    int num = -1;

    if (local < 186) {
        num = local + 1;
    }
    else if (local < 390) {
        int d = local - 186;
        int block = d / 12;
        int off   = d % 12;
        num = 205 + block * 30 + off;
    }
    else {
        num = 715 + (local - 390);
    }

    return "G" + std::to_string(num);
}
