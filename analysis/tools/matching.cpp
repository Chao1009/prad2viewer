// A quick check tool to test the matching result between HyCal clusters and GEM hits in the replay output
// Usage:
//   matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]
//   -o  output ROOT file (default: input filename with _matching.root suffix)
//   -n  max events to process (default: all)
// Example:
//   matching recon.root -o recon_matching.root -n 10000
//   matching recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "TrackMatcher.h"
#include "TrackGeometry.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
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

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    
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
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- load detector geometry config from JSON ---
    std::string run_str = get_run_str(root_files[0]);
    int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(dbDir + "/runinfo/2p1_general.json", run_num);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // --- TrackMatcher (target-seed, χ²-gated) -----------------------------
    // gRunConfig provides target + GEM placements; per-plane σ_GEM and HC
    // (A, B, C) coefficients aren't in RunConfig (they live in
    // reconstruction_config.json:matching), so we use working defaults
    // here.  Tools that need per-run tuning should load matching params
    // explicitly the way AppState does in src/app_state_init.cpp.
    namespace trk = prad2::trk;
    auto xforms = analysis::BuildLabTransforms(gRunConfig);
    trk::MatcherConfig mcfg;
    mcfg.planes.resize(4);
    for (int d = 0; d < 4; ++d) {
        mcfg.planes[d].id      = d;
        mcfg.planes[d].z       = gRunConfig.gem_z[d];
        mcfg.planes[d].sigma_x = 0.5f;          // mm — typical PRad-II GEM σ
        mcfg.planes[d].sigma_y = 0.5f;
        mcfg.planes[d].xform   = xforms.gem[d];
    }
    mcfg.hc_sigma_A = hycal.GetPositionResolutionA();
    mcfg.hc_sigma_B = hycal.GetPositionResolutionB();
    mcfg.hc_sigma_C = hycal.GetPositionResolutionC();
    mcfg.target_x = gRunConfig.target_x;
    mcfg.target_y = gRunConfig.target_y;
    mcfg.target_z = gRunConfig.target_z;
    mcfg.match_nsigma = 3.0f;
    mcfg.max_chi2     = 5.0f;        // loose; this is a quick-look tool
    trk::TrackMatcher matcher(std::move(mcfg));

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
    prad2::SetReconReadBranches(tree, ev);

    // --- output file ---
    TString outName = output;
    if (outName.IsNull()) {
        outName = root_files[0];
        outName.ReplaceAll("_recon.root", "_matching.root");
    }
    TFile outfile(outName, "RECREATE");
    int hitN = 0;
    float HC_X[100], HC_Y[100], HC_Z[100], HC_Energy[100];
    float Gup_X[100], Gup_Y[100], Gup_Z[100];
    float Gdown_X[100], Gdown_Y[100], Gdown_Z[100];
    float chi2[100];
    uint16_t center_id[100];
    uint32_t flag[100];

    TTree *outTree = new TTree("matching", "HyCal-GEM matching results");
    outTree->Branch("hitN", &hitN, "hitN/I");
    outTree->Branch("HC_X", HC_X, "HC_X[hitN]/F");
    outTree->Branch("HC_Y", HC_Y, "HC_Y[hitN]/F");
    outTree->Branch("HC_Z", HC_Z, "HC_Z[hitN]/F");
    outTree->Branch("HC_Energy", HC_Energy, "HC_Energy[hitN]/F");
    outTree->Branch("Gup_X", Gup_X, "Gup_X[hitN]/F");
    outTree->Branch("Gup_Y", Gup_Y, "Gup_Y[hitN]/F");
    outTree->Branch("Gup_Z", Gup_Z, "Gup_Z[hitN]/F");
    outTree->Branch("Gdown_X", Gdown_X, "Gdown_X[hitN]/F");
    outTree->Branch("Gdown_Y", Gdown_Y, "Gdown_Y[hitN]/F");
    outTree->Branch("Gdown_Z", Gdown_Z, "Gdown_Z[hitN]/F");
    outTree->Branch("chi2_per_dof", chi2, "chi2_per_dof[hitN]/F");
    outTree->Branch("center_id", center_id, "center_id[hitN]/s");
    outTree->Branch("flag", flag, "flag[hitN]/i");

    // --- event loop  ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;
    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "Reading " << i << " events / " << N << " total events\r" << std::flush;

        // Repackage clusters and GEM hits for TrackMatcher.  ev.cl_x /
        // ev.gem_x are already lab-frame (Replay applied DetectorTransform
        // before writing the tree), so no further frame change here.
        std::vector<std::vector<trk::PlaneHit>> hits_by_plane(4);
        for (int j = 0; j < ev.n_gem_hits; j++) {
            int d = ev.det_id[j];
            if (d < 0 || d >= 4) continue;
            trk::PlaneHit ph;
            ph.plane_idx = d;
            ph.hit_idx   = (int)hits_by_plane[d].size();
            ph.x = ev.gem_x[j]; ph.y = ev.gem_y[j]; ph.z = ev.gem_z[j];
            hits_by_plane[d].push_back(ph);
        }

        // One track per HyCal cluster.  Greedy-by-energy keeps the same
        // exclusivity policy the old MatchingTools::Match used: the highest-E
        // cluster gets first claim on overlapping GEM hits.
        struct ClusterIdx { int idx; float energy; };
        std::vector<ClusterIdx> order;
        order.reserve(ev.n_clusters);
        for (int j = 0; j < ev.n_clusters; ++j)
            order.push_back({j, ev.cl_energy[j]});
        std::sort(order.begin(), order.end(),
                  [](const ClusterIdx &a, const ClusterIdx &b) {
                      return a.energy > b.energy;
                  });

        for (auto [j, _] : order) {
            trk::ClusterHit hc;
            hc.x = ev.cl_x[j]; hc.y = ev.cl_y[j]; hc.z = ev.cl_z[j];
            hc.energy = ev.cl_energy[j];
            hc.center_id = ev.cl_center[j];
            hc.flag      = ev.cl_flag[j];

            // Target-seeded, all 4 GEM planes eligible, ≥3 matched required
            // (downstream + upstream pair must both contribute).
            auto track = matcher.findBestTrack(hc, hits_by_plane,
                                               /*free_seed=*/false,
                                               /*require_planes=*/3);
            if (!track || hitN >= 100) continue;

            // Pick best of {GEM1, GEM2} for "down" and {GEM3, GEM4} for
            // "up" — smaller post-fit residual (in lab frame at the
            // plane's z) wins the pair.  Falls back to whichever plane
            // has a hit when only one of the pair contributed.
            auto pick_pair = [&](int a, int b) -> trk::PlaneHit {
                trk::PlaneHit none{}; none.plane_idx = -1;
                bool ha = track->matched[a], hb = track->matched[b];
                if (!ha && !hb) return none;
                if (!ha) return track->hit[b];
                if (!hb) return track->hit[a];
                auto residSq = [&](const trk::PlaneHit &h) {
                    float px = track->fit.ax + track->fit.bx * h.z;
                    float py = track->fit.ay + track->fit.by * h.z;
                    float dx = h.x - px, dy = h.y - py;
                    return dx*dx + dy*dy;
                };
                return residSq(track->hit[a]) < residSq(track->hit[b])
                       ? track->hit[a] : track->hit[b];
            };
            trk::PlaneHit g_down = pick_pair(0, 1);
            trk::PlaneHit g_up   = pick_pair(2, 3);

            // Remove claimed GEM hits from `hits_by_plane` so a lower-E
            // cluster can't reuse them — preserves the greedy-by-E
            // exclusivity of the old MatchingTools::Match.
            auto erase_from = [&](int d, int hit_idx) {
                if (d < 0 || hit_idx < 0) return;
                auto &v = hits_by_plane[d];
                for (auto it = v.begin(); it != v.end(); ++it) {
                    if (it->hit_idx == hit_idx) { v.erase(it); break; }
                }
            };
            for (int d = 0; d < 4; ++d)
                if (track->matched[d])
                    erase_from(d, track->hit[d].hit_idx);

            // Project onto the HyCal face for downstream plotting (same
            // convention as the legacy tool).
            HCHit hc_proj{hc.x, hc.y, hc.z, hc.energy, hc.center_id, hc.flag};
            GEMHit gu{g_up.x,   g_up.y,   g_up.z,   2};
            GEMHit gd{g_down.x, g_down.y, g_down.z, 0};
            GetProjection(hc_proj, gRunConfig.hycal_z);
            if (g_up.plane_idx >= 0)   GetProjection(gu, gRunConfig.hycal_z);
            if (g_down.plane_idx >= 0) GetProjection(gd, gRunConfig.hycal_z);

            HC_X[hitN] = hc_proj.x;
            HC_Y[hitN] = hc_proj.y;
            HC_Z[hitN] = hc_proj.z;
            HC_Energy[hitN] = hc_proj.energy;
            Gup_X[hitN]   = (g_up.plane_idx   >= 0) ? gu.x : -999.f;
            Gup_Y[hitN]   = (g_up.plane_idx   >= 0) ? gu.y : -999.f;
            Gup_Z[hitN]   = (g_up.plane_idx   >= 0) ? gu.z : -999.f;
            Gdown_X[hitN] = (g_down.plane_idx >= 0) ? gd.x : -999.f;
            Gdown_Y[hitN] = (g_down.plane_idx >= 0) ? gd.y : -999.f;
            Gdown_Z[hitN] = (g_down.plane_idx >= 0) ? gd.z : -999.f;
            chi2[hitN]    = track->fit.chi2_per_dof;
            center_id[hitN] = hc_proj.center_id;
            flag[hitN]      = hc_proj.flag;
            hitN++;
        }
        outTree->Fill();
        hitN = 0; // reset for next event
    }
    std::cerr << "Finished processing " << N << " events.\n";
    outTree->Write();
    outfile.Close();
    std::cerr << "Output saved to " << outName << "\n";
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