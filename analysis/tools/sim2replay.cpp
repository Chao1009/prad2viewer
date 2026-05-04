// A tool to transform the Geant4 ouput into the replay_recon format,
// include matching between HyCal clusters and GEM hits
// can be used to compare with the replay reconstruction result and test your analysis
// scrirt, 
// Usage:
//   sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_luminosity(nb^-1)> <ee_luminosity(nb^-1)> [options]
//   -o  output ROOT file (default: sim_recon.root)
// Example:
//   sim2replay ep.root ee.root 1e6 1e6 -o sim_recon.root -n 100000 

#include "TrackMatcher.h"
#include "TrackGeometry.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TChain.h>

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

// --- geometry constants (can be made configurable) ---
const float gem_z[4] = {5407.f + 39.71f/2, 5407.f - 39.71f/2,
                        5807.f + 39.71f/2, 5807.f - 39.71f/2};
const float hycal_z = 6225.f;

double beamE = 3500.0; // MeV, can be made configurable

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;
//event data structure for the input Geant4 simulation output
struct SimEventData {
    int VD_n = 0; // number of hits in the virtual detector
    double VD_x[500] = {};
    double VD_y[500] = {};
    double VD_z[500] = {};
    double VD_E[500] = {};
    int GEM_n = 0; // number of hits in the GEMs
    int GEM_id[500] = {};
    double GEM_x_in[500] = {};
    double GEM_y_in[500] = {};
    double GEM_z_in[500] = {};
    double GEM_x_out[500] = {};
    double GEM_y_out[500] = {};
    double GEM_z_out[500] = {};
    double GEM_edep[500] = {};
};
void setupSimBranches(TTree *tree, SimEventData &ev)
{   
    tree->SetBranchAddress("VD.N", &ev.VD_n);
    tree->SetBranchAddress("VD.X", ev.VD_x);
    tree->SetBranchAddress("VD.Y", ev.VD_y);
    tree->SetBranchAddress("VD.Z", ev.VD_z);
    tree->SetBranchAddress("VD.P", ev.VD_E);
    tree->SetBranchAddress("GEM.N", &ev.GEM_n);
    tree->SetBranchAddress("GEM.DID", ev.GEM_id);
    tree->SetBranchAddress("GEM.X", ev.GEM_x_in);
    tree->SetBranchAddress("GEM.Y", ev.GEM_y_in);
    tree->SetBranchAddress("GEM.Z", ev.GEM_z_in);
    tree->SetBranchAddress("GEM.Xout", ev.GEM_x_out);
    tree->SetBranchAddress("GEM.Yout", ev.GEM_y_out);
    tree->SetBranchAddress("GEM.Zout", ev.GEM_z_out);
    tree->SetBranchAddress("GEM.Edep", ev.GEM_edep);
}

static std::vector<std::string> collectRootFiles(const std::string &path);
float EResolution(float E);

// Find the HyCal module ID (PrimEx ID) whose area contains (x, y).
// Returns the module's PrimEx ID, or -1 if not found.
static int findModuleID(const fdec::HyCalSystem &hycal, double x, double y)
{
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        if (std::abs(x - m.x) <= m.size_x / 2. &&
            std::abs(y - m.y) <= m.size_y / 2.)
            return m.id;
    }
    return -1;
}

int main (int argc, char *argv[])
{
    std::string ep_file_dir, ee_file_dir;
    double ep_lumi = 0., ee_lumi = 0.;
    std::string outName = "sim_recon.root";
    int opt;
    while ((opt = getopt(argc, argv, "o:")) != -1) {
        switch (opt) {
            case 'o': outName = optarg; break;
        }
    }
    if (optind + 4 > argc) {
        std::cerr << "Usage: sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_lumi(nb^-1)> <ee_lumi(nb^-1)> [options]\n";
        std::cerr << "Options:\n  -o  output ROOT file (default: sim_recon.root)\n";
        return 1;
    }
    ep_file_dir = argv[optind];
    ee_file_dir = argv[optind + 1];
    ep_lumi = std::stod(argv[optind + 2]);
    ee_lumi = std::stod(argv[optind + 3]);

    // collect input files (can be files, directories)
    std::vector<std::string> ep_files, ee_files;
    auto f = collectRootFiles(ep_file_dir);
    ep_files.insert(ep_files.end(), f.begin(), f.end());
    f = collectRootFiles(ee_file_dir);
    ee_files.insert(ee_files.end(), f.begin(), f.end());

    if (ep_files.empty() || ee_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_lumi(nb^-1)> <ee_lumi(nb^-1)> [options]\n";
        return 1;
    }
    if(ep_lumi <= 0. || ee_lumi <= 0.) {
        std::cerr << "Invalid luminosity values.\n";
        std::cerr << "Usage: sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_lumi(nb^-1)> <ee_lumi(nb^-1)> [options]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // --- TrackMatcher (target-seed, χ²-gated) -----------------------------
    namespace trk = prad2::trk;
    trk::MatcherConfig mcfg;
    mcfg.planes.resize(4);
    for (int d = 0; d < 4; ++d) {
        mcfg.planes[d].id      = d;
        mcfg.planes[d].z       = gem_z[d];
        mcfg.planes[d].sigma_x = 0.5f;          // mm — typical PRad-II GEM σ
        mcfg.planes[d].sigma_y = 0.5f;
        mcfg.planes[d].xform.set(0.f, 0.f, gem_z[d], 0.f, 0.f, 0.f);
    }
    mcfg.hc_sigma_A = 2.5f;
    mcfg.hc_sigma_B = 0.f;
    mcfg.hc_sigma_C = 0.f;
    mcfg.target_x = 0.f; mcfg.target_y = 0.f; mcfg.target_z = 0.f;
    mcfg.match_nsigma = 3.0f;
    mcfg.max_chi2     = 5.0f;
    trk::TrackMatcher matcher(std::move(mcfg));

    // --- setup TChain and branches ---
    auto sim = std::make_unique<SimEventData>();
    TChain *chain_ee = new TChain("T");
    for (const auto &f : ee_files) {
        chain_ee->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree_ee = chain_ee;
    if (!tree_ee) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }
    setupSimBranches(tree_ee, *sim);
    TChain *chain_ep = new TChain("T");
    for (const auto &f : ep_files) {
        chain_ep->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree_ep = chain_ep;
    if (!tree_ep) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }
    setupSimBranches(tree_ep, *sim);

    TFile *outfile = TFile::Open(outName.c_str(), "RECREATE");
    // create TTree and branches for reconstructed data
    TTree *tree_out = new TTree("recon", "PRad2 replay reconstruction from G4");
    auto ev = std::make_unique<EventVars_Recon>();
    prad2::SetReconWriteBranches(tree_out, *ev);

    // caculate luminosity and number of events to process for ep and ee
    double lumi = std::min(ep_lumi, ee_lumi);
    int N_ep = lumi / ep_lumi * tree_ep->GetEntries();
    int N_ee = lumi / ee_lumi * tree_ee->GetEntries();

    std::cerr << "Processing " << N_ep << " ep events and " << N_ee << " ee events (luminosity = " << lumi << " nb^-1)\n";

    // --- loop over events and fill output tree ---
    int ep_count = 0, ee_count = 0;
    for(int i = 0; i < N_ep + N_ee; i++){
        if(i + 1 % 4 == 0 && ep_count < N_ep) {
            tree_ep->GetEntry(ep_count);
            ep_count++;
        } else if(ee_count < N_ee) {
            tree_ee->GetEntry(ee_count);
            ee_count++;
        }
        if(ep_count == N_ep && ee_count < N_ee) {
            tree_ee->GetEntry(ee_count);
            ee_count++;
        }
        if(ee_count == N_ee && ep_count < N_ep) {
            tree_ep->GetEntry(ep_count);
            ep_count++;
        }
        if(ee_count == N_ee && ep_count == N_ep) break;
        if((i + 1) % 1000 == 0) {
            std::cerr << "Processing event " << i + 1 << "/" << N_ep + N_ee << "\r" << std::flush;
        }

        // Here you can fill the ev structure with the HyCal virtual plane hit and GEM hit information
        // For example, you can use the matching tools to match HyCal clusters and GEM hits
        // and fill the ev.matchG_x, etc. arrays accordingly
        // This part will depend on how you want to do the reconstruction and matching based on the input simulation data structure
        // You can also use the physics tools to calculate expected energies, angles, etc. for the clusters and hits, and fill the ev structure with those values as well

        //first total energy cut
        *ev = EventVars_Recon{}; // reset all fields for this event
        float total_energy = 0.f;
        for(int j = 0; j < sim->VD_n; j++){
            if(sim->VD_E[j] < 1. / 300.f * beamE) continue; // add some energy threshold to reduce noise
            sim->VD_E[j] += gRandom->Gaus(0, EResolution(sim->VD_E[j])); // add some energy smearing
            total_energy += sim->VD_E[j];
        }
        if(total_energy < 0.5f * beamE) continue; // only keep events with total energy above 0.5*beamE

        ev->event_num = i;
        ev->trigger_type = 8; // primary_bit 8
        ev->trigger_bits = (1 << 8); //ssp raw sum
        ev->timestamp = 0;
        ev->total_energy = total_energy;

        for(int j = 0; j < sim->VD_n; j++){
            if(sim->VD_E[j] < 1. / 300.f * beamE) continue; // add some energy threshold to reduce noise
            ev->cl_x[ev->n_clusters] = float(sim->VD_x[j] + gRandom->Gaus(0, 2.6/sqrt(sim->VD_E[j] / 1000.f)));
            ev->cl_y[ev->n_clusters] = float(sim->VD_y[j] + gRandom->Gaus(0, 2.6/sqrt(sim->VD_E[j] / 1000.f)));
            ev->cl_z[ev->n_clusters] = float(hycal_z);
            ev->cl_energy[ev->n_clusters] = float(sim->VD_E[j]);
            ev->cl_nblocks[ev->n_clusters] = 1;
            ev->cl_center[ev->n_clusters] = findModuleID(hycal, sim->VD_x[j], sim->VD_y[j]);
            ev->cl_flag[ev->n_clusters] = ( 1 << 1); // set flag to PbWO4
            ev->n_clusters++;
        }

        for(int j = 0; j < sim->GEM_n; j++){
            if(sim->GEM_edep[j] < 26.e-6*2.) continue; // add some energy threshold to reduce noise
            ev->det_id[ev->n_gem_hits] = sim->GEM_id[j];
            ev->gem_x[ev->n_gem_hits] = float(0.5*(sim->GEM_x_in[j] + sim->GEM_x_out[j]) + gRandom->Gaus(0, 0.07)); // add some position smearing
            ev->gem_y[ev->n_gem_hits] = float(0.5*(sim->GEM_y_in[j] + sim->GEM_y_out[j]) + gRandom->Gaus(0, 0.07));
            ev->gem_x_charge[ev->n_gem_hits] = 0.f;
            ev->gem_y_charge[ev->n_gem_hits] = 0.f;
            ev->gem_x_peak[ev->n_gem_hits] = 0.f;
            ev->gem_y_peak[ev->n_gem_hits] = 0.f;
            ev->gem_x_size[ev->n_gem_hits] = 3; // set some default cluster size
            ev->gem_y_size[ev->n_gem_hits] = 3;
            ev->n_gem_hits++;
        }

        // Repackage clusters / GEM hits for TrackMatcher.  Sim positions
        // are already lab-frame (no detector pose applied), so PlaneHit
        // and ClusterHit get them directly.
        std::vector<std::vector<trk::PlaneHit>> hits_by_plane(4);
        for (int i = 0; i < ev->n_gem_hits; ++i) {
            int d = ev->det_id[i];
            if (d < 0 || d >= 4) continue;
            trk::PlaneHit ph;
            ph.plane_idx = d;
            ph.hit_idx   = (int)hits_by_plane[d].size();
            ph.x = ev->gem_x[i];
            ph.y = ev->gem_y[i];
            ph.z = gem_z[d];
            hits_by_plane[d].push_back(ph);
        }

        // Reset all cluster match fields (default = no match).
        for (int i = 0; i < ev->n_clusters; ++i) {
            for (int j = 0; j < 2; j++) {
                ev->matchGEMx[i][j] = -999.f;
                ev->matchGEMy[i][j] = -999.f;
                ev->matchGEMz[i][j] = -999.f;
            }
            ev->matchFlag[i] = 0;
        }

        // Greedy-by-energy: highest-E cluster gets first claim on
        // overlapping GEM hits — same exclusivity policy MatchingTools::Match
        // applied.
        struct Idx { int i; float E; };
        std::vector<Idx> order;
        order.reserve(ev->n_clusters);
        for (int i = 0; i < ev->n_clusters; ++i)
            order.push_back({i, ev->cl_energy[i]});
        std::sort(order.begin(), order.end(),
                  [](const Idx &a, const Idx &b) { return a.E > b.E; });

        ev->matchNum = 0;
        for (auto [cl_idx, _] : order) {
            trk::ClusterHit hc;
            hc.x = ev->cl_x[cl_idx]; hc.y = ev->cl_y[cl_idx]; hc.z = ev->cl_z[cl_idx];
            hc.energy = ev->cl_energy[cl_idx];
            hc.center_id = ev->cl_center[cl_idx];
            hc.flag      = ev->cl_flag[cl_idx];

            auto track = matcher.findBestTrack(hc, hits_by_plane,
                                               /*free_seed=*/false,
                                               /*require_planes=*/3);
            if (!track) continue;

            // Pair-best (smallest residual) for the down (0,1) and up (2,3) layers.
            auto pick_pair = [&](int a, int b) -> trk::PlaneHit {
                trk::PlaneHit none{}; none.plane_idx = -1;
                bool ha = track->matched[a], hb = track->matched[b];
                if (!ha && !hb) return none;
                if (!ha) return track->hit[b];
                if (!hb) return track->hit[a];
                auto rs = [&](const trk::PlaneHit &h) {
                    float px = track->fit.ax + track->fit.bx * h.z;
                    float py = track->fit.ay + track->fit.by * h.z;
                    float dx = h.x - px, dy = h.y - py;
                    return dx*dx + dy*dy;
                };
                return rs(track->hit[a]) < rs(track->hit[b])
                       ? track->hit[a] : track->hit[b];
            };
            trk::PlaneHit g_down = pick_pair(0, 1);
            trk::PlaneHit g_up   = pick_pair(2, 3);

            // Set per-cluster match arrays (j=0 down, j=1 up — matches the
            // historical mHit_g* / matchGEM* layout).
            const trk::PlaneHit *gp[2] = {&g_down, &g_up};
            uint32_t mflag = 0;
            for (int j = 0; j < 2; ++j) {
                if (gp[j]->plane_idx < 0) continue;
                ev->matchGEMx[cl_idx][j] = gp[j]->x;
                ev->matchGEMy[cl_idx][j] = gp[j]->y;
                ev->matchGEMz[cl_idx][j] = gp[j]->z;
                mflag |= (1u << gp[j]->plane_idx);
            }
            ev->matchFlag[cl_idx] = mflag;

            int slot = ev->matchNum;
            if (slot < prad2::kMaxClusters) {
                ev->mHit_E[slot] = hc.energy;
                ev->mHit_x[slot] = hc.x;
                ev->mHit_y[slot] = hc.y;
                ev->mHit_z[slot] = hc.z;
                for (int j = 0; j < 2; ++j) {
                    ev->mHit_gx[slot][j]  = (gp[j]->plane_idx >= 0) ? gp[j]->x : -999.f;
                    ev->mHit_gy[slot][j]  = (gp[j]->plane_idx >= 0) ? gp[j]->y : -999.f;
                    ev->mHit_gz[slot][j]  = (gp[j]->plane_idx >= 0) ? gp[j]->z : -999.f;
                    ev->mHit_gid[slot][j] = (gp[j]->plane_idx >= 0) ? gp[j]->plane_idx : 0;
                }
                ev->matchNum = slot + 1;
            }

            // Remove claimed GEM hits to enforce greedy-by-E exclusivity.
            for (int d = 0; d < 4; ++d) {
                if (!track->matched[d]) continue;
                int hi = track->hit[d].hit_idx;
                auto &v = hits_by_plane[d];
                for (auto it = v.begin(); it != v.end(); ++it) {
                    if (it->hit_idx == hi) { v.erase(it); break; }
                }
            }
        }

        tree_out->Fill();
    }
    outfile->Write();
    delete outfile;
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find(".root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

float EResolution(float E)
{
    // Example energy resolution function, can be modified based on actual detector performance
    return E * 0.026f / sqrt(E / 1000.f); // 2.6% at 1 GeV, scaling with sqrt(E)
}

