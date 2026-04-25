//============================================================================
// gem_clusters_to_root.C — full GEM analysis pipeline → ROOT tree.
//
// Pipeline per event:
//   EvChannel.Read() → DecodeEvent() → SspEventData
//                   → GemSystem.ProcessEvent() (pedestal + CM + ZS)
//                   → GemSystem.Reconstruct(GemCluster) (1D + 2D clustering)
//                   → loop over per-plane clusters, copy hits + 6 time samples
//                   → TTree.Fill()
//
// Output tree layout (one entry per physics event):
//
//   event_num         I       physics event sequence number
//   trigger_bits      i       FP trigger bits (32-bit mask)
//
//   ncl               I       number of 1D strip clusters in this event
//   cl_det            vector<int>     [ncl]   detector id (0..3)
//   cl_plane          vector<int>     [ncl]   0 = X, 1 = Y
//   cl_position       vector<float>   [ncl]   charge-weighted strip pos (mm)
//   cl_peak_charge    vector<float>   [ncl]   highest strip charge in cluster
//   cl_total_charge   vector<float>   [ncl]   sum of strip charges
//   cl_max_tb         vector<short>   [ncl]   peak time sample of cluster
//   cl_cross_talk     vector<bool>    [ncl]   firmware-classified cross-talk
//   cl_nhits          vector<int>     [ncl]   #strips in this cluster
//   cl_first          vector<int>     [ncl]   index into hit_* arrays for the
//                                             first strip of this cluster
//
//   nhits             I       sum of cl_nhits across all clusters
//   hit_cl            vector<int>     [nhits] cluster index this strip belongs to
//   hit_strip         vector<int>     [nhits] plane-wise strip number
//   hit_position      vector<float>   [nhits] strip physical position (mm)
//   hit_charge        vector<float>   [nhits] max charge across time samples
//   hit_max_tb        vector<short>   [nhits] time sample with max charge
//   hit_cross_talk    vector<bool>    [nhits]
//   hit_ts0..ts5      vector<float>   [nhits] full 6-sample ADC waveform
//
// Usage
// -----
//   cd build
//   root -l ../analysis/scripts/rootlogon.C
//
//   // simplest — pedestal/CM auto-discovered from database/config.json's
//   // runinfo (latest entry):
//   .x ../analysis/scripts/gem_clusters_to_root.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.00000", \
//       "gem_clusters_023867.root")
//
//   // pin to a specific run's runinfo entry:
//   .x ../analysis/scripts/gem_clusters_to_root.C+( \
//       "/data/.../prad_023867.evio.00000", \
//       "out.root", nullptr, nullptr, 0, /*run_num=*/23867)
//
//   // explicit pedestal / CM override (paths relative to PRAD2_DATABASE_DIR
//   // or absolute):
//   .x ../analysis/scripts/gem_clusters_to_root.C+( \
//       "/data/.../prad_023867.evio.00000", \
//       "out.root", \
//       "gem_peds/peds_23867.txt", \
//       "gem_peds/cm_23867.txt")
//============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"
#include "SspData.h"

#include "GemSystem.h"
#include "GemCluster.h"
#include "RunInfoConfig.h"     // prad2::LoadRunConfig (header-only)

#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace evc;

namespace {

// Resolve a database-relative path against PRAD2_DATABASE_DIR (set by
// rootlogon.C); returns the input unchanged if it's already absolute.
static std::string resolve_db_path(const std::string &p)
{
    if (p.empty()) return p;
    if (p.size() >= 1 && (p[0] == '/' || p[0] == '\\')) return p;
    if (p.size() >= 2 && p[1] == ':') return p;   // Windows drive letter
    const char *db = std::getenv("PRAD2_DATABASE_DIR");
    if (!db) return p;
    return std::string(db) + "/" + p;
}

// Build the GEM crate-remap from daq_config.json's roc_tags table — same
// logic as src/app_state_init.cpp so the analysis tree matches what the
// monitor sees.  Empty result if no GEM ROCs are configured.
static std::map<int, int> build_gem_crate_remap(const DaqConfig &cfg)
{
    std::map<int, int> remap;
    for (const auto &re : cfg.roc_tags)
        if (re.type == "gem") remap[(int)re.tag] = re.crate;
    return remap;
}

// Read the "runinfo" pointer out of database/config.json and return the
// resolved path to the runinfo JSON.  Returns empty string if not found.
static std::string discover_runinfo_path()
{
    const char *db = std::getenv("PRAD2_DATABASE_DIR");
    std::string db_dir = db ? db : "database";

    std::ifstream f(db_dir + "/config.json");
    if (!f) return {};
    auto j = nlohmann::json::parse(f, nullptr, false, true);
    if (j.is_discarded() || !j.contains("runinfo") || !j["runinfo"].is_string())
        return {};

    std::string rel = j["runinfo"].get<std::string>();
    if (rel.empty()) return {};
    if (rel[0] == '/' || rel[0] == '\\') return rel;          // absolute
    if (rel.size() >= 2 && rel[1] == ':') return rel;          // Windows drive
    return db_dir + "/" + rel;
}

// Per-event vectors filled by the loop, bound to the TTree once at setup.
struct EventVars {
    int   event_num    = 0;
    uint32_t trigger_bits = 0;

    int   ncl  = 0;
    std::vector<int>   cl_det;
    std::vector<int>   cl_plane;
    std::vector<float> cl_position;
    std::vector<float> cl_peak_charge;
    std::vector<float> cl_total_charge;
    std::vector<short> cl_max_tb;
    std::vector<bool>  cl_cross_talk;
    std::vector<int>   cl_nhits;
    std::vector<int>   cl_first;

    int   nhits = 0;
    std::vector<int>   hit_cl;
    std::vector<int>   hit_strip;
    std::vector<float> hit_position;
    std::vector<float> hit_charge;
    std::vector<short> hit_max_tb;
    std::vector<bool>  hit_cross_talk;
    std::vector<float> hit_ts[6];

    void clear()
    {
        ncl = 0;
        cl_det.clear(); cl_plane.clear();
        cl_position.clear(); cl_peak_charge.clear(); cl_total_charge.clear();
        cl_max_tb.clear(); cl_cross_talk.clear();
        cl_nhits.clear();   cl_first.clear();

        nhits = 0;
        hit_cl.clear();   hit_strip.clear();
        hit_position.clear(); hit_charge.clear();
        hit_max_tb.clear(); hit_cross_talk.clear();
        for (int t = 0; t < 6; ++t) hit_ts[t].clear();
    }
};

static void bind_branches(TTree *tree, EventVars &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");

    tree->Branch("ncl",             &ev.ncl, "ncl/I");
    tree->Branch("cl_det",          &ev.cl_det);
    tree->Branch("cl_plane",        &ev.cl_plane);
    tree->Branch("cl_position",     &ev.cl_position);
    tree->Branch("cl_peak_charge",  &ev.cl_peak_charge);
    tree->Branch("cl_total_charge", &ev.cl_total_charge);
    tree->Branch("cl_max_tb",       &ev.cl_max_tb);
    tree->Branch("cl_cross_talk",   &ev.cl_cross_talk);
    tree->Branch("cl_nhits",        &ev.cl_nhits);
    tree->Branch("cl_first",        &ev.cl_first);

    tree->Branch("nhits",           &ev.nhits, "nhits/I");
    tree->Branch("hit_cl",          &ev.hit_cl);
    tree->Branch("hit_strip",       &ev.hit_strip);
    tree->Branch("hit_position",    &ev.hit_position);
    tree->Branch("hit_charge",      &ev.hit_charge);
    tree->Branch("hit_max_tb",      &ev.hit_max_tb);
    tree->Branch("hit_cross_talk",  &ev.hit_cross_talk);
    for (int t = 0; t < 6; ++t)
        tree->Branch(TString::Format("hit_ts%d", t), &ev.hit_ts[t]);
}

} // anonymous namespace

//=============================================================================
// Entry point
//=============================================================================
// Pass nullptr for gem_ped_file / gem_cm_file to auto-discover via the
// runinfo pointed to by database/config.json (falls back to the entry
// matching `run_num` — pass -1 for "latest").  Pass an explicit path to
// override.
int gem_clusters_to_root(const char *evio_path,
                         const char *out_path     = "gem_clusters.root",
                         const char *gem_ped_file = nullptr,
                         const char *gem_cm_file  = nullptr,
                         long        max_events   = 0,
                         int         run_num      = -1,
                         const char *daq_config   = nullptr,
                         const char *gem_map_file = nullptr)
{
    //---- DAQ config + GEM map -----------------------------------------------
    std::string daq_path = daq_config ? daq_config : "";
    if (daq_path.empty())
        daq_path = resolve_db_path("daq_config.json");
    DaqConfig cfg;
    if (!load_daq_config(daq_path, cfg)) {
        std::cerr << "ERROR: cannot load DAQ config " << daq_path << "\n";
        return 1;
    }
    std::cout << "DAQ config : " << daq_path << "\n";

    std::string map_path = gem_map_file ? gem_map_file : "";
    if (map_path.empty()) map_path = resolve_db_path("gem_map.json");

    //---- GemSystem (decoder + clustering) -----------------------------------
    gem::GemSystem    gem_sys;
    gem::GemCluster   gem_clusterer;
    gem_sys.Init(map_path);
    std::cout << "GEM map    : " << map_path
              << "  (" << gem_sys.GetNDetectors() << " detectors)\n";

    auto remap = build_gem_crate_remap(cfg);

    // Pedestal + common-mode files: explicit args win, otherwise look up
    // via runinfo (same path the live monitor uses through app_state_init).
    std::string ped_path = gem_ped_file ? resolve_db_path(gem_ped_file) : "";
    std::string cm_path  = gem_cm_file  ? resolve_db_path(gem_cm_file)  : "";
    if (ped_path.empty() || cm_path.empty()) {
        std::string ri_path = discover_runinfo_path();
        if (ri_path.empty()) {
            std::cerr << "WARN: no runinfo pointer in database/config.json — "
                         "GEM pedestal / common-mode auto-discovery skipped.\n";
        } else {
            std::cout << "RunInfo    : " << ri_path
                      << "  (run_num=" << run_num << ")\n";
            prad2::RunConfig rc = prad2::LoadRunConfig(ri_path, run_num);
            if (ped_path.empty() && !rc.gem_pedestal_file.empty())
                ped_path = resolve_db_path(rc.gem_pedestal_file);
            if (cm_path.empty() && !rc.gem_common_mode_file.empty())
                cm_path = resolve_db_path(rc.gem_common_mode_file);
        }
    }

    if (!ped_path.empty()) {
        gem_sys.LoadPedestals(ped_path, remap);
        std::cout << "Pedestals  : " << ped_path << "\n";
    } else {
        std::cerr << "WARN: no pedestal file — full-readout data will reconstruct empty.\n";
    }
    if (!cm_path.empty()) {
        gem_sys.LoadCommonModeRange(cm_path, remap);
        std::cout << "Common mode: " << cm_path << "\n";
    }

    //---- EVIO file ----------------------------------------------------------
    EvChannel ch;
    ch.SetConfig(cfg);
    if (ch.OpenAuto(evio_path) != status::success) {
        std::cerr << "ERROR: cannot open " << evio_path << "\n";
        return 1;
    }
    std::cout << "EVIO       : " << evio_path << "\n";

    //---- ROOT output --------------------------------------------------------
    TFile fout(out_path, "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "ERROR: cannot create " << out_path << "\n";
        return 1;
    }
    TTree *tree = new TTree("gem", "GEM 1D clusters and constituent strip hits");
    EventVars ev;
    bind_branches(tree, ev);

    //---- event loop ---------------------------------------------------------
    auto t0 = std::chrono::steady_clock::now();
    fdec::EventData    fadc_evt;
    ssp::SspEventData  ssp_evt;

    long n_read = 0, n_phys = 0, n_filled = 0;
    long total_clusters = 0, total_hits = 0;

    while (ch.Read() == status::success) {
        ++n_read;
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ssp_evt.clear();
            // DecodeEvent() fills both FADC + SSP — we only need SSP for
            // GEMs, but FADC gives us the trigger bits + event number.
            if (!ch.DecodeEvent(i, fadc_evt, &ssp_evt)) continue;
            ++n_phys;

            ev.clear();
            ev.event_num    = static_cast<int>(fadc_evt.info.event_number);
            ev.trigger_bits = fadc_evt.info.trigger_bits;

            // GEM pipeline.  Clear() is mandatory: ProcessEvent appends
            // into per-APV buffers, so a stale event would leak through.
            gem_sys.Clear();
            gem_sys.ProcessEvent(ssp_evt);
            gem_sys.Reconstruct(gem_clusterer);

            // Walk every (det, plane) pair and copy clusters + hits into
            // the flat per-event arrays.
            for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
                for (int p = 0; p < 2; ++p) {
                    const auto &clusters = gem_sys.GetPlaneClusters(d, p);
                    for (const auto &cl : clusters) {
                        const int icl = ev.ncl++;
                        ev.cl_det.push_back(d);
                        ev.cl_plane.push_back(p);
                        ev.cl_position.push_back(cl.position);
                        ev.cl_peak_charge.push_back(cl.peak_charge);
                        ev.cl_total_charge.push_back(cl.total_charge);
                        ev.cl_max_tb.push_back(cl.max_timebin);
                        ev.cl_cross_talk.push_back(cl.cross_talk);
                        ev.cl_nhits.push_back(static_cast<int>(cl.hits.size()));
                        ev.cl_first.push_back(ev.nhits);

                        for (const auto &h : cl.hits) {
                            ev.hit_cl.push_back(icl);
                            ev.hit_strip.push_back(h.strip);
                            ev.hit_position.push_back(h.position);
                            ev.hit_charge.push_back(h.charge);
                            ev.hit_max_tb.push_back(h.max_timebin);
                            ev.hit_cross_talk.push_back(h.cross_talk);
                            // ts_adc is sized SSP_TIME_SAMPLES (= 6) by
                            // the strip pipeline; clamp defensively.
                            for (int t = 0; t < 6; ++t) {
                                float v = (t < (int)h.ts_adc.size())
                                    ? h.ts_adc[t] : 0.f;
                                ev.hit_ts[t].push_back(v);
                            }
                            ++ev.nhits;
                        }
                    }
                }
            }

            tree->Fill();
            ++n_filled;
            total_clusters += ev.ncl;
            total_hits     += ev.nhits;

            if (max_events > 0 && n_phys >= max_events) goto done;
        }
        if (n_phys % 10000 == 0 && n_phys > 0)
            std::cerr << "  " << n_phys << " physics events...\r" << std::flush;
    }

done:
    auto t1 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    fout.cd();
    tree->Write();
    fout.Close();

    std::cout << "\n";
    std::cout << "EVIO records          : " << n_read         << "\n";
    std::cout << "physics events        : " << n_phys         << "\n";
    std::cout << "tree entries written  : " << n_filled       << "\n";
    std::cout << "total clusters        : " << total_clusters << "\n";
    std::cout << "total strip hits      : " << total_hits     << "\n";
    std::cout << "elapsed (s)           : " << secs           << "\n";
    std::cout << "wrote                 : " << out_path       << "\n";
    return 0;
}
