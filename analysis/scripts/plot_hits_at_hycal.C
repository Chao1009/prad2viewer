//============================================================================
// plot_hits_at_hycal.C — 2D occupancy of GEM hits (projected to HyCal
// surface) and HyCal cluster centroids (already on HyCal surface), drawn
// side-by-side in the lab / target-centered, beam-aligned frame at z =
// hycal_z.
//
// Pipeline per physics event:
//   EvChannel.Read()  → DecodeEvent() → FADC + SSP buffers
//                     → HyCal: WaveAnalyzer → energize → HyCalCluster
//                     → GEM:   GemSystem.ProcessEvent → Reconstruct
//                     → coord transform to lab (per-detector tilt + offset
//                       via prad2det's DetectorTransform)
//                     → GEM hits: GetProjection(hits, hycal_z) — straight
//                       line from target through (x,y,z) to z=hycal_z
//                     → fill the two TH2F occupancy maps
//
// Both plots share x/y range and binning.  GEM hits from all four
// detectors are combined into the left histogram; HyCal cluster centroids
// (one entry per cluster) populate the right histogram.
//
// Trigger filter: only events with `trigger_bits == 0x100` (production
// physics trigger) contribute.  Everything else (LMS / Alpha / cosmic /
// etc.) is skipped.
//
// Multi-file mode is selected by the input path:
//   * `/data/.../prad_023881.evio.*`  → glob: enumerate every sibling
//     `prad_023881.evio.<digits>`, fold them all into the same two
//     histograms, and warn (to stderr) about any gap in the suffix
//     sequence (including missing from .00000).
//   * `/data/prad_023881/`            → directory: same enumeration,
//     run number sniffed from the directory name.
//   * `/data/.../prad_023881.evio.00000` → single specific split file.
//
// Heap-allocate the big POD-ish decoder structs (fdec::EventData,
// ssp::SspEventData) — see the project's `feedback_heap_allocate_decoder
// _structs` memory: stack-allocating them SEGVs at function prologue.
//
// Usage
// -----
//   cd build
//   root -l ../analysis/scripts/rootlogon.C
//
//   # full run (glob — warns about any missing split):
//   .x ../analysis/scripts/plot_hits_at_hycal.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.*", \
//       "hits_at_hycal.pdf")
//
//   # single split (debugging):
//   .x ../analysis/scripts/plot_hits_at_hycal.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.00000", \
//       "hits_at_hycal_seg0.pdf")
//
//   args (full): evio_path, out_path, max_events, run_num,
//                gem_ped_file, gem_cm_file, hc_calib_file,
//                daq_config, gem_map_file, hc_map_file
//   - out_path  : PDF/PNG/etc. for the canvas; an accompanying .root
//                 file alongside it stores both TH2Fs for re-plotting.
//   - max_events: 0 = all
//   - run_num   : -1 = sniff from EVIO basename (prad_NNNNNN.evio.*)
//   - all "_file"/"daq_config" args: "" = auto-discover via runinfo
//============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "WaveAnalyzer.h"

#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "GemSystem.h"
#include "GemCluster.h"
#include "RunInfoConfig.h"

#include "PhysicsTools.h"     // analysis::GetProjection (moved here from MatchingTools.h)
#include "ConfigSetup.h"      // analysis::ApplyToLab
#include "PipelineBuilder.h"  // prad2::PipelineBuilder — one-stop wiring
#include "script_helpers.h"   // discover_split_files, strip_extension

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace evc;

//=============================================================================
// Forward declaration of the full 10-arg version + convenience overloads
// (cling default-arg marshalling is buggy for mixed-type signatures).
//=============================================================================
int plot_hits_at_hycal(const char *evio_path,
                       const char *out_path,
                       long        max_events,
                       int         run_num,
                       const char *gem_ped_file,
                       const char *gem_cm_file,
                       const char *hc_calib_file,
                       const char *daq_config,
                       const char *gem_map_file,
                       const char *hc_map_file);

int plot_hits_at_hycal(const char *evio_path, const char *out_path)
{
    return plot_hits_at_hycal(evio_path, out_path,
                              0L, -1, "", "", "", "", "", "");
}
int plot_hits_at_hycal(const char *evio_path, const char *out_path,
                       long max_events)
{
    return plot_hits_at_hycal(evio_path, out_path,
                              max_events, -1, "", "", "", "", "", "");
}
int plot_hits_at_hycal(const char *evio_path, const char *out_path,
                       long max_events, int run_num)
{
    return plot_hits_at_hycal(evio_path, out_path,
                              max_events, run_num, "", "", "", "", "", "");
}

//=============================================================================
// Entry point — full version
//=============================================================================
int plot_hits_at_hycal(const char *evio_path,
                       const char *out_path,
                       long        max_events,
                       int         run_num,
                       const char *gem_ped_file,
                       const char *gem_cm_file,
                       const char *hc_calib_file,
                       const char *daq_config,
                       const char *gem_map_file,
                       const char *hc_map_file)
{
    auto blank = [](const char *s) -> bool { return !s || !*s; };

    //---- detector pipeline (DAQ config, runinfo, HyCal, GEM) ----------------
    // One call wires daq_config + reconstruction_config + runinfo +
    // HyCalSystem (Init + LoadCalibration + position-resolution) + GemSystem
    // (Init + LoadPedestals + LoadCommonModeRange + per-detector configs)
    // + DetectorTransforms.  Crate-remap derivation, JSON parsing, and
    // [GEMSYS]/[GEMCFG]/[PEDSUM] diagnostic prints all live inside.
    prad2::Pipeline pipeline;
    try {
        pipeline = prad2::PipelineBuilder()
            .set_daq_config(blank(daq_config)        ? "" : daq_config)
            .set_hycal_calib(blank(hc_calib_file)    ? "" : hc_calib_file)
            .set_gem_pedestal(blank(gem_ped_file)    ? "" : gem_ped_file)
            .set_gem_common_mode(blank(gem_cm_file)  ? "" : gem_cm_file)
            .set_hycal_map(blank(hc_map_file)        ? "" : hc_map_file)
            .set_gem_map(blank(gem_map_file)         ? "" : gem_map_file)
            .set_run_number(run_num > 0 ? run_num : -1)
            .set_run_number_from_evio(evio_path ? evio_path : "")
            .build();
    } catch (const std::exception &e) {
        Printf("[ERROR] %s", e.what());
        return 1;
    }

    auto &cfg            = pipeline.daq_cfg;
    auto &geo            = pipeline.run_cfg;
    auto &hycal          = pipeline.hycal;
    auto &gem_sys        = pipeline.gem;
    auto &hycal_xform    = pipeline.hycal_transform;
    auto &gem_xforms     = pipeline.gem_transforms;
    analysis::gRunConfig = geo;

    std::map<int, int> crate_map;
    for (const auto &re : cfg.roc_tags) crate_map[(int)re.tag] = re.crate;

    fdec::HyCalCluster hc_clusterer(hycal);
    hc_clusterer.SetConfig(pipeline.hycal_cluster_cfg);
    gem::GemCluster gem_clusterer;

    //---- EVIO discovery -----------------------------------------------------
    EvChannel ch;
    ch.SetConfig(cfg);
    auto evio_files = discover_split_files(evio_path ? evio_path : "");
    Printf("[setup] EVIO       : %zu split file(s) for input %s",
           evio_files.size(), evio_path ? evio_path : "(null)");
    for (const auto &f : evio_files) Printf("           %s", f.c_str());

    //---- histograms ---------------------------------------------------------
    // Lab frame (target-centered, beam-aligned) at z = hycal_z.
    // Range +/-650 mm covers PRad-II HyCal LG ring (~580 mm to outer edge);
    // 5 mm bins give a clean occupancy map without being too noisy.
    constexpr float kRange = 650.f;
    constexpr int   kBins  = 260;        // 5 mm bins
    auto h_gem = std::make_unique<TH2F>(
        "h_gem_at_hycal",
        TString::Format(
            "GEM hits projected to HyCal surface (z = %.0f mm);x (mm);y (mm)",
            geo.hycal_z),
        kBins, -kRange, kRange, kBins, -kRange, kRange);
    auto h_hc = std::make_unique<TH2F>(
        "h_hycal",
        TString::Format(
            "HyCal cluster centroids on HyCal surface (z = %.0f mm);x (mm);y (mm)",
            geo.hycal_z),
        kBins, -kRange, kRange, kBins, -kRange, kRange);

    //---- event loop ---------------------------------------------------------
    auto t0 = std::chrono::steady_clock::now();
    auto fadc_evt_ptr = std::make_unique<fdec::EventData>();
    auto ssp_evt_ptr  = std::make_unique<ssp::SspEventData>();
    auto &fadc_evt    = *fadc_evt_ptr;
    auto &ssp_evt     = *ssp_evt_ptr;
    fdec::WaveAnalyzer ana;
    fdec::WaveResult   wres;

    long n_read = 0, n_phys = 0, n_kept = 0;
    long n_files_open = 0;
    long n_hc_clusters = 0, n_gem_hits = 0;

    for (const auto &fpath : evio_files) {
        if (ch.OpenAuto(fpath) != status::success) {
            Printf("[WARN] skip (cannot open): %s", fpath.c_str());
            continue;
        }
        ++n_files_open;
        Printf("[file %ld/%zu] %s",
               n_files_open, evio_files.size(), fpath.c_str());

        while (ch.Read() == status::success) {
            ++n_read;
            if (!ch.Scan()) continue;
            if (ch.GetEventType() != EventType::Physics) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ssp_evt.clear();
            hc_clusterer.Clear();
            if (!ch.DecodeEvent(i, fadc_evt, &ssp_evt)) continue;
            ++n_phys;

            // Trigger filter: keep only events with trigger_bits exactly
            // == 0x100.  max_events still gates against n_phys (raw
            // physics count) so file-scan extent stays predictable.
            if (fadc_evt.info.trigger_bits != 0x100u) {
                if (max_events > 0 && n_phys >= max_events) goto done;
                continue;
            }
            ++n_kept;

            // ---------- HyCal: waveform → energy → clusters ----------
            for (int r = 0; r < fadc_evt.nrocs; ++r) {
                auto &roc = fadc_evt.rocs[r];
                if (!roc.present) continue;
                auto cit = crate_map.find(roc.tag);
                if (cit == crate_map.end()) continue;
                const int crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    auto &slot = roc.slots[s];
                    if (!slot.present) continue;
                    for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                        if (!(slot.channel_mask & (1ull << c))) continue;
                        const auto *mod = hycal.module_by_daq(crate, s, c);
                        if (!mod || !mod->is_hycal()) continue;
                        auto &cd = slot.channels[c];
                        if (cd.nsamples <= 0) continue;
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        if (wres.npeaks <= 0) continue;
                        int   best = -1;
                        float best_h = -1.f;
                        for (int p = 0; p < wres.npeaks; ++p) {
                            const auto &pk = wres.peaks[p];
                            if (pk.time > 100.f && pk.time < 200.f
                                && pk.height > best_h) {
                                best_h = pk.height; best = p;
                            }
                        }
                        if (best < 0) continue;
                        float energy = static_cast<float>(
                            mod->energize(wres.peaks[best].integral));
                        hc_clusterer.AddHit(mod->index, energy, 0.f);
                    }
                }
            }
            hc_clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hc_raw;
            hc_clusterer.ReconstructHits(hc_raw);

            // Build HCHit list with z = 0 (no shower depth) so transform
            // lands them at exactly z = hycal_z — i.e. on the HyCal face.
            std::vector<analysis::HCHit> hc_hits;
            hc_hits.reserve(hc_raw.size());
            for (const auto &h : hc_raw) {
                analysis::HCHit hh;
                hh.x = h.x; hh.y = h.y; hh.z = 0.f;
                hh.energy    = h.energy;
                hh.center_id = h.center_id;
                hh.flag      = h.flag;
                hc_hits.push_back(hh);
            }
            for (auto &h : hc_hits) analysis::ApplyToLab(hycal_xform, h);

            for (const auto &h : hc_hits) h_hc->Fill(h.x, h.y);
            n_hc_clusters += hc_hits.size();

            // ---------- GEM: pedestal → CM → ZS → 1D + 2D ----------
            gem_sys.Clear();
            gem_sys.ProcessEvent(ssp_evt);
            gem_sys.Reconstruct(gem_clusterer);

            // Per-detector lab-frame hit lists.  GetHits(d) returns local
            // plane hits (x, y, z=0); rotate + transform per-detector,
            // then project the line target->hit onto z = hycal_z.
            for (int d = 0; d < gem_sys.GetNDetectors() && d < 4; ++d) {
                const auto &raw = gem_sys.GetHits(d);
                if (raw.empty()) continue;
                std::vector<analysis::GEMHit> lab;
                lab.reserve(raw.size());
                for (const auto &h : raw) {
                    analysis::GEMHit gh;
                    gh.x = h.x; gh.y = h.y; gh.z = 0.f;
                    gh.det_id = d;
                    lab.push_back(gh);
                }
                for (auto &h : lab) analysis::ApplyToLab(gem_xforms[d], h);
                analysis::GetProjection(lab, geo.hycal_z);
                for (const auto &g : lab) h_gem->Fill(g.x, g.y);
                n_gem_hits += lab.size();
            }

            if (max_events > 0 && n_phys >= max_events) goto done;
        }
        if (n_phys > 0 && n_phys % 5000 == 0)
            Printf("[progress] %ld physics events", n_phys);
        }
        ch.Close();
    }

done:
    auto t1 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    //---- draw + save --------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);

    TCanvas c("c_hits_at_hycal", "Hits at HyCal surface", 1600, 720);
    c.Divide(2, 1, 0.005, 0.005);

    c.cd(1);
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.10);
    h_gem->Draw("COLZ");

    c.cd(2);
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.10);
    h_hc->Draw("COLZ");

    c.SaveAs(out_path);

    // Sibling .root file with the two histograms for re-plotting.
    std::string root_out = strip_extension(out_path) + ".root";
    TFile fout(root_out.c_str(), "RECREATE");
    if (!fout.IsZombie()) {
        h_gem->Write();
        h_hc->Write();
        c.Write();
        fout.Close();
        Printf("[setup] Saved hists: %s", root_out.c_str());
    }

    Printf("--- summary ---");
    Printf("  EVIO files opened     : %ld / %zu", n_files_open, evio_files.size());
    Printf("  EVIO records          : %ld", n_read);
    Printf("  physics events        : %ld", n_phys);
    Printf("  passed trig cut 0x100 : %ld", n_kept);
    Printf("  HyCal clusters total  : %ld  (avg %.2f / kept event)",
           n_hc_clusters,
           n_kept ? double(n_hc_clusters) / n_kept : 0.0);
    Printf("  GEM hits total (4 det): %ld  (avg %.2f / kept event)",
           n_gem_hits,
           n_kept ? double(n_gem_hits) / n_kept : 0.0);
    Printf("  elapsed (s)           : %.2f", secs);
    Printf("  wrote canvas          : %s", out_path);
    return 0;
}
