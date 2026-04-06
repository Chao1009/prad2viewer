
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TF1.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;
// ── Tree branch struct ───────────────────────────────────────────────────
void setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    tree->SetBranchAddress("event_num",    &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",    &ev.timestamp);
    tree->SetBranchAddress("total_energy", &ev.total_energy);
    // HyCal cluster branches
    tree->SetBranchAddress("n_clusters",   &ev.n_clusters);
    tree->SetBranchAddress("cl_x",         ev.cl_x);
    tree->SetBranchAddress("cl_y",         ev.cl_y);
    tree->SetBranchAddress("cl_z",         ev.cl_z);
    tree->SetBranchAddress("cl_energy",    ev.cl_energy);
    tree->SetBranchAddress("cl_nblocks",   ev.cl_nblocks);
    tree->SetBranchAddress("cl_center",    ev.cl_center);
    tree->SetBranchAddress("cl_flag",      ev.cl_flag);
    // GEM part
    tree->SetBranchAddress("n_gem_hits",   &ev.n_gem_hits);
    tree->SetBranchAddress("det_id",       ev.det_id);
    tree->SetBranchAddress("gem_x",        ev.gem_x);
    tree->SetBranchAddress("gem_y",        ev.gem_y);
    tree->SetBranchAddress("gem_x_charge", ev.gem_x_charge);
    tree->SetBranchAddress("gem_y_charge", ev.gem_y_charge);
    tree->SetBranchAddress("gem_x_peak",   ev.gem_x_peak);
    tree->SetBranchAddress("gem_y_peak",   ev.gem_y_peak);
    tree->SetBranchAddress("gem_x_size",   ev.gem_x_size);
    tree->SetBranchAddress("gem_y_size",   ev.gem_y_size);
    // Matching results
    tree->SetBranchAddress("match_num",       &ev.match_num);
    tree->SetBranchAddress("matchHC_x",       ev.matchHC_x);
    tree->SetBranchAddress("matchHC_y",       ev.matchHC_y);
    tree->SetBranchAddress("matchHC_z",       ev.matchHC_z);
    tree->SetBranchAddress("matchHC_energy",  ev.matchHC_energy);
    tree->SetBranchAddress("matchHC_center",  ev.matchHC_center);
    tree->SetBranchAddress("matchHC_flag",    ev.matchHC_flag);
    tree->SetBranchAddress("matchG_x",        ev.matchG_x);
    tree->SetBranchAddress("matchG_y",        ev.matchG_y);
    tree->SetBranchAddress("matchG_z",        ev.matchG_z);
    tree->SetBranchAddress("matchG_det_id",   ev.matchG_det_id);
}

static std::vector<std::string> collectRootFiles(const std::string &path);
void projectToHyCalSurface(PhysicsTools::MollerData &m_data, float hycal_z);
double fitAndDraw(TH1F* hist, const std::string& out_path, const double fit_range = 4.);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 1100.f;
    float hycal_z = 6225.f; //mm, the default position of HyCal surface, TO DO: read from database
    int run_id = 12345;
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }
    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = DATABASE_DIR;

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_modules.json",
               dbDir + "/daq_map.json");
    PhysicsTools physics(hycal);

    // --- setup TChain and branches ---
    TChain *chain = new TChain("recon");
    for (const auto &f : root_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev;
    setupReconBranches(tree, ev);

    //output histograms for calibration results
    TH1F *vertex_hycal = new TH1F("vertex_hycal", "Moller vertex z distance HyCal;Z (mm);Counts", 2000, 5500, 7500);
    TH2F *center_hycal = new TH2F("center_hycal", "Moller center distribution HyCal;X (mm);Y (mm)", 200, -100, 100, 200, -100, 100);
    TH1F *center_hycal_x = new TH1F("center_hycal_x", "Moller center X distribution HyCal;X (mm);Counts", 200, -50, 50);
    TH1F *center_hycal_y = new TH1F("center_hycal_y", "Moller center Y distribution HyCal;Y (mm);Counts", 200, -50, 50);

    TH1F *vertex_gem[4];
    TH2F *center_gem[4];
    TH1F *center_gem_x[4];
    TH1F *center_gem_y[4];
    for (int d = 0; d < 4; d++) {
        vertex_gem[d] = new TH1F(Form("vertex_gem%d", d), Form("Moller vertex z distance GEM%d;Z (mm);Counts", d), 2000, 5500, 7500);
        center_gem[d] = new TH2F(Form("center_gem%d", d), Form("Moller center distribution GEM%d;X (mm);Y (mm)", d), 200, -100, 100, 200, -100, 100);
        center_gem_x[d] = new TH1F(Form("center_gem_x%d", d), Form("Moller center X distribution GEM%d;X (mm);Counts", d), 200, -50, 50);
        center_gem_y[d] = new TH1F(Form("center_gem_y%d", d), Form("Moller center Y distribution GEM%d;Y (mm);Counts", d), 200, -50, 50);
    }

    // --- output file ---
    TString outName = output;
    if (outName.IsNull()) {
        outName = root_files[0];
        outName.ReplaceAll("_recon.root", "_posCalib.root");
    }
    TFile outfile(outName, "RECREATE");

    PhysicsTools::MollerData hycal_mollers;
    PhysicsTools::MollerData gem_mollers[4];

    // --- event loop : select Moller events on HyCal and each GEM plane ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;

    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "\rPass 1: " << i << " / " << N << std::flush;

        bool good_moller = false;
        if(ev.match_num == 2){
            float Epair = ev.matchHC_energy[0] + ev.matchHC_energy[1];
            if (std::abs(Epair - Ebeam) < 5.f * Ebeam * 0.025f / std::sqrt(Ebeam / 1000.f)) {
                good_moller = true;
            }
        }
        if(!good_moller) continue;

        //have selected good Moller events for further analysis
        PhysicsTools::MollerEvent h_m;
        PhysicsTools::MollerEvent g_m;
        
        h_m = PhysicsTools::MollerEvent(
                {ev.matchHC_x[0], ev.matchHC_y[0], ev.matchHC_z[0], ev.matchHC_energy[0]},
                {ev.matchHC_x[1], ev.matchHC_y[1], ev.matchHC_z[1], ev.matchHC_energy[1]});
        hycal_mollers.push_back(h_m);
        
        // select two moller on one chamber for upstream GEMs 
        if(ev.matchG_det_id[0][0] == ev.matchG_det_id[1][0]){
            g_m = PhysicsTools::MollerEvent(
                {ev.matchG_x[0][0], ev.matchG_y[0][0], ev.matchG_z[0][0], ev.matchHC_energy[0]},
                {ev.matchG_x[1][0], ev.matchG_y[1][0], ev.matchG_z[1][0], ev.matchHC_energy[1]});
            int det_id = ev.matchG_det_id[0][0];
            if(det_id >= 0 && det_id < 4) gem_mollers[det_id].push_back(g_m);
            else std::cerr << "Warning: Invalid GEM det_id " << det_id << " in event " << ev.event_num << "\n";
        }

        // select two moller on one chamber for downstream GEMs
        if(ev.matchG_det_id[0][1] == ev.matchG_det_id[1][1]){
            g_m = PhysicsTools::MollerEvent(
                {ev.matchG_x[0][1], ev.matchG_y[0][1], ev.matchG_z[0][1], ev.matchHC_energy[0]},
                {ev.matchG_x[1][1], ev.matchG_y[1][1], ev.matchG_z[1][1], ev.matchHC_energy[1]});
            int det_id = ev.matchG_det_id[0][1];
            if(det_id >= 0 && det_id < 4) gem_mollers[det_id].push_back(g_m);
            else std::cerr << "Warning: Invalid GEM det_id " << det_id << " in event " << ev.event_num << "\n";
        }
    }

    // After collecting Moller events, analyze them for detector calibration
    //summary of Moller events on each detector plane
    std::cerr << "\nSummary of selected Moller events:\n";
    std::cerr << "HyCal: " << hycal_mollers.size() << " events\n";
    for (int d = 0; d < 4; d++) {
        std::cerr << "GEM " << d << ": " << gem_mollers[d].size() << " events\n";
    }

    //hycal Moller events
    projectToHyCalSurface(hycal_mollers, hycal_z); //project to HyCal surface
    for (int i = 0; i < hycal_mollers.size(); i++) {
        vertex_hycal->Fill(physics.GetMollerZdistance(hycal_mollers[i], Ebeam));
        if (i >= 1) {
            auto c = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            center_hycal->Fill(c[0], c[1]);
            center_hycal_x->Fill(c[0]);
            center_hycal_y->Fill(c[1]);
        }
    }

    //gem Moller events
    for (int d = 0; d < 4; d++) {
        for (int i = 0; i < gem_mollers[d].size(); i++) {
            vertex_gem[d]->Fill(physics.GetMollerZdistance(gem_mollers[d][i], Ebeam));
            if (i >= 1) {
                auto c = physics.GetMollerCenter(gem_mollers[d][i-1], gem_mollers[d][i]);
                center_gem[d]->Fill(c[0], c[1]);
                center_gem_x[d]->Fill(c[0]);
                center_gem_y[d]->Fill(c[1]);
            }
        }
    }

    //fit histograms, and get the beam position and vertex distance for each detector plane
    double hycal_vertex_z = fitAndDraw(vertex_hycal, "calib_result/hycal_vertex_z", 100.);
    double hycal_center_x = fitAndDraw(center_hycal_x, "calib_result/hycal_center_x", 2.);
    double hycal_center_y = fitAndDraw(center_hycal_y, "calib_result/hycal_center_y", 2.);
    double gem_vertex_z[4];
    double gem_center_x[4];
    double gem_center_y[4];
    for (int d = 0; d < 4; d++) {
        gem_vertex_z[d] = fitAndDraw(vertex_gem[d], Form("calib_result/gem%d_vertex_z", d), 50.);
        gem_center_x[d] = fitAndDraw(center_gem_x[d], Form("calib_result/gem%d_center_x", d), 1.);
        gem_center_y[d] = fitAndDraw(center_gem_y[d], Form("calib_result/gem%d_center_y", d), 4.);
    }
    //print summary of calibration results
    std::cerr << "HyCal vertex z distance: " << hycal_vertex_z << " mm (pre-entered number " << hycal_z << " mm)" << "\n";
    std::cerr << "HyCal center x: " << hycal_center_x << " mm\n";
    std::cerr << "HyCal center y: " << hycal_center_y << " mm\n";
    for (int d = 0; d < 4; d++) {
        std::cerr << "GEM " << d << " vertex z distance: " << gem_vertex_z[d] << " mm\n";
        std::cerr << "GEM " << d << " center x: " << gem_center_x[d] << " mm\n";
        std::cerr << "GEM " << d << " center y: " << gem_center_y[d] << " mm\n";
    }

    //save histograms
    outfile.cd();
    vertex_hycal->Write();
    center_hycal->Write();
    center_hycal_x->Write();
    center_hycal_y->Write();
    for (int d = 0; d < 4; d++) {
        vertex_gem[d]->Write();
        center_gem[d]->Write();
        center_gem_x[d]->Write();
        center_gem_y[d]->Write();
    }
    outfile.Close();
    std::cerr << "Calibration histograms saved to " << outName << "\n";

}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_recon.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

void projectToHyCalSurface(PhysicsTools::MollerData &m_data, float hycal_z)
{
    //project the Moller event from target center(z = 0) to the HyCal surface (z = hycal_z)
    for (auto &evt : m_data) {
        for (auto *dp : {&evt.first, &evt.second}) {
            float scale = hycal_z / dp->z;
            dp->x = dp->x * scale;
            dp->y = dp->y * scale;
            dp->z = hycal_z;
        }
    }
}

double fitAndDraw(TH1F* hist, const std::string& out_path, const double fit_range){
    TCanvas *c = new TCanvas("", "", 800, 600);
    double mean = hist->GetBinCenter(hist->GetMaximumBin());
    hist->Fit("gaus", "rq", "", mean-fit_range, mean+fit_range);
    hist->Draw();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.15, 0.85, Form("%.2f mm +- %.2f mm", hist->GetFunction("gaus")->GetParameter(1), hist->GetFunction("gaus")->GetParError(1)));
    c->SaveAs((out_path + ".png").c_str());
    delete c;

    return hist->GetFunction("gaus")->GetParameter(1);
}