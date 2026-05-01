#pragma once
//=============================================================================
// EventData_io.h — TTree branch I/O helpers for RawEventData / ReconEventData
//
// Header-only.  prad2det's library does NOT compile this file; it is included
// only by ROOT-aware consumers (analysis tools, the viewer's root data
// source, sim2replay).  The library still has no link-time ROOT dependency.
//
// All four helpers are inline.  Writers are unconditional: every call sets
// up the same branch list given the `with_peaks` flag.  Readers use
// TTree::GetBranch to skip any branch that's missing on disk, so they
// happily read older files that pre-date the firmware-peak / ssp_raw /
// per-cluster-flag additions.
//
// Usage:
//   #include "EventData_io.h"
//   prad2::SetRawWriteBranches(tree, ev, /*with_peaks=*/true);
//   ...
//   auto status = prad2::SetRawReadBranches(tree, ev);
//   if (status.has_peaks) { ... }
//=============================================================================
//
// Single source of truth for the replay tree schema.  Match this against
// `analysis/REPLAYED_DATA.md` (or vice versa) when the layout changes.

#include "EventData.h"

#include <TTree.h>
#include <TString.h>   // for Form()

namespace prad2 {

// ── Reader status ────────────────────────────────────────────────────────
struct RawReadStatus {
    bool has_peaks      = false;   // soft-analyzer: ped_mean/rms + peak_*
    bool has_daq_peaks  = false;   // firmware-mode: daq_peak_*
    bool has_gem        = false;   // gem.*
    bool has_ssp_raw    = false;   // ssp_raw vector
};

struct ReconReadStatus {
    bool has_match_num  = false;   // mHit_* quick-access arrays
    bool has_per_cl_match = false; // matchFlag / matchGEM*
    bool has_veto       = false;
    bool has_lms        = false;
    bool has_ssp_raw    = false;
};

// ─────────────────────────────────────────────────────────────────────────
// Raw "events" tree — write
// ─────────────────────────────────────────────────────────────────────────
inline void SetRawWriteBranches(TTree *tree, RawEventData &ev, bool with_peaks)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_type", &ev.trigger_type, "trigger_type/b");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");

    // Unified FADC250 channel array (HyCal + Veto + LMS).  Categorisation
    // is via hycal.module_type per channel — HyCal consumers using
    // hycal.module_by_id() naturally skip SCINT/LMS entries because their
    // module_id values (3001+ / 3100+) are not registered in HyCalSystem.
    tree->Branch("hycal.nch",         &ev.nch,         "hycal.nch/I");
    tree->Branch("hycal.module_id",   ev.module_id,    "hycal.module_id[hycal.nch]/s");
    tree->Branch("hycal.module_type", ev.module_type,  "hycal.module_type[hycal.nch]/b");
    tree->Branch("hycal.nsamples",    ev.nsamples,     "hycal.nsamples[hycal.nch]/b");
    tree->Branch("hycal.samples",     ev.samples,
                 Form("hycal.samples[hycal.nch][%d]/s", fdec::MAX_SAMPLES));
    tree->Branch("hycal.gain_factor", ev.gain_factor,  "hycal.gain_factor[hycal.nch]/F");

    if (with_peaks) {
        // ped_mean / ped_rms / ped_nused / ped_quality / ped_slope are
        // products of the soft analyzer — only meaningful when the
        // analyzer ran.  ped_quality is a Q_PED_* bitmask
        // (NOT_CONVERGED / FLOOR_ACTIVE / TOO_FEW_SAMPLES /
        // PULSE_IN_WINDOW / OVERFLOW / TRAILING_WINDOW; see Fadc250Data.h).
        tree->Branch("hycal.ped_mean",      ev.ped_mean,       "hycal.ped_mean[hycal.nch]/F");
        tree->Branch("hycal.ped_rms",       ev.ped_rms,        "hycal.ped_rms[hycal.nch]/F");
        tree->Branch("hycal.ped_nused",     ev.ped_nused,      "hycal.ped_nused[hycal.nch]/b");
        tree->Branch("hycal.ped_quality",   ev.ped_quality,    "hycal.ped_quality[hycal.nch]/b");
        tree->Branch("hycal.ped_slope",     ev.ped_slope,      "hycal.ped_slope[hycal.nch]/F");
        tree->Branch("hycal.npeaks",        ev.npeaks,         "hycal.npeaks[hycal.nch]/b");
        tree->Branch("hycal.peak_height",   ev.peak_height,
                     Form("hycal.peak_height[hycal.nch][%d]/F",   fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_time",     ev.peak_time,
                     Form("hycal.peak_time[hycal.nch][%d]/F",     fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_integral", ev.peak_integral,
                     Form("hycal.peak_integral[hycal.nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_quality",  ev.peak_quality,
                     Form("hycal.peak_quality[hycal.nch][%d]/b",  fdec::MAX_PEAKS));

        // Firmware-mode (FADC250 Modes 1/2/3) emulation peaks.
        // daq_peak_quality is a Q_* bitmask (peak-at-boundary,
        // NSB/NSA truncation, Va out-of-range).
        tree->Branch("hycal.daq_npeaks",        ev.daq_npeaks,    "hycal.daq_npeaks[hycal.nch]/b");
        tree->Branch("hycal.daq_peak_vp",       ev.daq_peak_vp,
                     Form("hycal.daq_peak_vp[hycal.nch][%d]/F",       fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_integral", ev.daq_peak_integral,
                     Form("hycal.daq_peak_integral[hycal.nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_time",     ev.daq_peak_time,
                     Form("hycal.daq_peak_time[hycal.nch][%d]/F",     fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_cross",    ev.daq_peak_cross,
                     Form("hycal.daq_peak_cross[hycal.nch][%d]/I",    fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_pos",      ev.daq_peak_pos,
                     Form("hycal.daq_peak_pos[hycal.nch][%d]/I",      fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_coarse",   ev.daq_peak_coarse,
                     Form("hycal.daq_peak_coarse[hycal.nch][%d]/I",   fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_fine",     ev.daq_peak_fine,
                     Form("hycal.daq_peak_fine[hycal.nch][%d]/I",     fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_quality",  ev.daq_peak_quality,
                     Form("hycal.daq_peak_quality[hycal.nch][%d]/b",  fdec::MAX_PEAKS));
    }

    // GEM strip data.
    tree->Branch("gem.nch",         &ev.gem_nch,     "gem.nch/I");
    tree->Branch("gem.mpd_crate",   ev.mpd_crate,    "gem.mpd_crate[gem.nch]/b");
    tree->Branch("gem.mpd_fiber",   ev.mpd_fiber,    "gem.mpd_fiber[gem.nch]/b");
    tree->Branch("gem.apv",         ev.apv,          "gem.apv[gem.nch]/b");
    tree->Branch("gem.strip",       ev.strip,        "gem.strip[gem.nch]/b");
    tree->Branch("gem.ssp_samples", ev.ssp_samples,
                 Form("gem.ssp_samples[gem.nch][%d]/S", ssp::SSP_TIME_SAMPLES));

    // Raw 0xE10C SSP trigger bank words.
    tree->Branch("ssp_raw", &ev.ssp_raw);
}

// ─────────────────────────────────────────────────────────────────────────
// Raw "events" tree — read
// Binds the addresses of every branch that exists on `tree`; reports back
// which optional groups are present.
// ─────────────────────────────────────────────────────────────────────────
inline RawReadStatus SetRawReadBranches(TTree *tree, RawEventData &ev)
{
    RawReadStatus s;
    auto bind = [&](const char *name, void *addr) {
        if (tree->GetBranch(name)) tree->SetBranchAddress(name, addr);
    };

    bind("event_num",    &ev.event_num);
    bind("trigger_type", &ev.trigger_type);
    bind("trigger_bits", &ev.trigger_bits);
    bind("timestamp",    &ev.timestamp);

    bind("hycal.nch",         &ev.nch);
    bind("hycal.module_id",   ev.module_id);
    bind("hycal.module_type", ev.module_type);
    bind("hycal.nsamples",    ev.nsamples);
    bind("hycal.samples",     ev.samples);
    bind("hycal.gain_factor", ev.gain_factor);

    s.has_peaks = (tree->GetBranch("hycal.npeaks") != nullptr);
    if (s.has_peaks) {
        bind("hycal.ped_mean",      ev.ped_mean);
        bind("hycal.ped_rms",       ev.ped_rms);
        // ped_nused / ped_quality / ped_slope are post-Mar-2026 additions —
        // bind() silently no-ops on older files that pre-date them.
        bind("hycal.ped_nused",     ev.ped_nused);
        bind("hycal.ped_quality",   ev.ped_quality);
        bind("hycal.ped_slope",     ev.ped_slope);
        bind("hycal.npeaks",        ev.npeaks);
        bind("hycal.peak_height",   ev.peak_height);
        bind("hycal.peak_time",     ev.peak_time);
        bind("hycal.peak_integral", ev.peak_integral);
        // peak_quality is a post-Mar-2026 addition — bind() no-ops on
        // older files that pre-date it.
        bind("hycal.peak_quality",  ev.peak_quality);
    }

    s.has_daq_peaks = (tree->GetBranch("hycal.daq_npeaks") != nullptr);
    if (s.has_daq_peaks) {
        bind("hycal.daq_npeaks",        ev.daq_npeaks);
        bind("hycal.daq_peak_vp",       ev.daq_peak_vp);
        bind("hycal.daq_peak_integral", ev.daq_peak_integral);
        bind("hycal.daq_peak_time",     ev.daq_peak_time);
        bind("hycal.daq_peak_cross",    ev.daq_peak_cross);
        bind("hycal.daq_peak_pos",      ev.daq_peak_pos);
        bind("hycal.daq_peak_coarse",   ev.daq_peak_coarse);
        bind("hycal.daq_peak_fine",     ev.daq_peak_fine);
        bind("hycal.daq_peak_quality",  ev.daq_peak_quality);
    }

    s.has_gem = (tree->GetBranch("gem.nch") != nullptr);
    if (s.has_gem) {
        bind("gem.nch",         &ev.gem_nch);
        bind("gem.mpd_crate",   ev.mpd_crate);
        bind("gem.mpd_fiber",   ev.mpd_fiber);
        bind("gem.apv",         ev.apv);
        bind("gem.strip",       ev.strip);
        bind("gem.ssp_samples", ev.ssp_samples);
    }

    // ssp_raw is std::vector<uint32_t>: ROOT needs a stable
    // `vector<uint32_t>**` address.  Consumers that need it must bind it
    // themselves with their own held pointer:
    //   auto *p = &ev.ssp_raw;
    //   tree->SetBranchAddress("ssp_raw", &p);   // p must outlive GetEntry
    s.has_ssp_raw = (tree->GetBranch("ssp_raw") != nullptr);

    return s;
}

// ─────────────────────────────────────────────────────────────────────────
// Recon tree — write
// ─────────────────────────────────────────────────────────────────────────
inline void SetReconWriteBranches(TTree *tree, ReconEventData &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_type", &ev.trigger_type, "trigger_type/b");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    tree->Branch("total_energy", &ev.total_energy, "total_energy/F");

    // HyCal cluster branches (lab frame: target/beam-centred).
    tree->Branch("n_clusters", &ev.n_clusters, "n_clusters/I");
    tree->Branch("cl_x",       ev.cl_x,        "cl_x[n_clusters]/F");
    tree->Branch("cl_y",       ev.cl_y,        "cl_y[n_clusters]/F");
    tree->Branch("cl_z",       ev.cl_z,        "cl_z[n_clusters]/F");
    tree->Branch("cl_energy",  ev.cl_energy,   "cl_energy[n_clusters]/F");
    tree->Branch("cl_nblocks", ev.cl_nblocks,  "cl_nblocks[n_clusters]/b");
    tree->Branch("cl_center",  ev.cl_center,   "cl_center[n_clusters]/s");
    tree->Branch("cl_flag",    ev.cl_flag,     "cl_flag[n_clusters]/i");

    // Per-cluster HyCal↔GEM matches (one row per HyCal cluster, 4 GEMs).
    tree->Branch("matchFlag", ev.matchFlag, "matchFlag[n_clusters]/i");
    tree->Branch("matchGEMx", ev.matchGEMx, "matchGEMx[n_clusters][4]/F");
    tree->Branch("matchGEMy", ev.matchGEMy, "matchGEMy[n_clusters][4]/F");
    tree->Branch("matchGEMz", ev.matchGEMz, "matchGEMz[n_clusters][4]/F");

    // Quick-access matched pairs (clusters with ≥2 GEMs matched).
    tree->Branch("match_num", &ev.matchNum, "match_num/I");
    tree->Branch("mHit_E",  ev.mHit_E,  "mHit_E[match_num]/F");
    tree->Branch("mHit_x",  ev.mHit_x,  "mHit_x[match_num]/F");
    tree->Branch("mHit_y",  ev.mHit_y,  "mHit_y[match_num]/F");
    tree->Branch("mHit_z",  ev.mHit_z,  "mHit_z[match_num]/F");
    tree->Branch("mHit_gx", ev.mHit_gx, "mHit_gx[match_num][2]/F");
    tree->Branch("mHit_gy", ev.mHit_gy, "mHit_gy[match_num][2]/F");
    tree->Branch("mHit_gz", ev.mHit_gz, "mHit_gz[match_num][2]/F");
    tree->Branch("mHit_gid", ev.mHit_gid, "mHit_gid[match_num][2]/F");

    // GEM hits (lab frame, per-detector plane).
    tree->Branch("n_gem_hits",   &ev.n_gem_hits,   "n_gem_hits/I");
    tree->Branch("det_id",       ev.det_id,        "det_id[n_gem_hits]/b");
    tree->Branch("gem_x",        ev.gem_x,         "gem_x[n_gem_hits]/F");
    tree->Branch("gem_y",        ev.gem_y,         "gem_y[n_gem_hits]/F");
    tree->Branch("gem_z",        ev.gem_z,         "gem_z[n_gem_hits]/F");
    tree->Branch("gem_x_charge", ev.gem_x_charge,  "gem_x_charge[n_gem_hits]/F");
    tree->Branch("gem_y_charge", ev.gem_y_charge,  "gem_y_charge[n_gem_hits]/F");
    tree->Branch("gem_x_peak",   ev.gem_x_peak,    "gem_x_peak[n_gem_hits]/F");
    tree->Branch("gem_y_peak",   ev.gem_y_peak,    "gem_y_peak[n_gem_hits]/F");
    tree->Branch("gem_x_size",   ev.gem_x_size,    "gem_x_size[n_gem_hits]/b");
    tree->Branch("gem_y_size",   ev.gem_y_size,    "gem_y_size[n_gem_hits]/b");
    tree->Branch("gem_x_mTbin",  ev.gem_x_mTbin,   "gem_x_mTbin[n_gem_hits]/b");
    tree->Branch("gem_y_mTbin",  ev.gem_y_mTbin,   "gem_y_mTbin[n_gem_hits]/b");

    // Veto + LMS soft-peak summaries.
    tree->Branch("veto_nch",         &ev.veto_nch,         "veto_nch/I");
    tree->Branch("veto_id",          ev.veto_id,           "veto_id[veto_nch]/b");
    tree->Branch("veto_npeaks",      ev.veto_npeaks,       "veto_npeaks[veto_nch]/I");
    tree->Branch("veto_peak_time",   ev.veto_peak_time,
                 Form("veto_peak_time[veto_nch][%d]/F",     fdec::MAX_PEAKS));
    tree->Branch("veto_peak_integral", ev.veto_peak_integral,
                 Form("veto_peak_integral[veto_nch][%d]/F", fdec::MAX_PEAKS));

    tree->Branch("lms_nch",         &ev.lms_nch,         "lms_nch/I");
    tree->Branch("lms_id",          ev.lms_id,           "lms_id[lms_nch]/b");
    tree->Branch("lms_npeaks",      ev.lms_npeaks,       "lms_npeaks[lms_nch]/I");
    tree->Branch("lms_peak_time",   ev.lms_peak_time,
                 Form("lms_peak_time[lms_nch][%d]/F",     fdec::MAX_PEAKS));
    tree->Branch("lms_peak_integral", ev.lms_peak_integral,
                 Form("lms_peak_integral[lms_nch][%d]/F", fdec::MAX_PEAKS));

    // Raw 0xE10C SSP trigger bank words.
    tree->Branch("ssp_raw", &ev.ssp_raw);
}

// ─────────────────────────────────────────────────────────────────────────
// Recon tree — read
// ─────────────────────────────────────────────────────────────────────────
inline ReconReadStatus SetReconReadBranches(TTree *tree, ReconEventData &ev)
{
    ReconReadStatus s;
    auto bind = [&](const char *name, void *addr) {
        if (tree->GetBranch(name)) tree->SetBranchAddress(name, addr);
    };

    bind("event_num",    &ev.event_num);
    bind("trigger_type", &ev.trigger_type);
    bind("trigger_bits", &ev.trigger_bits);
    bind("timestamp",    &ev.timestamp);
    bind("total_energy", &ev.total_energy);

    bind("n_clusters", &ev.n_clusters);
    bind("cl_x",       ev.cl_x);
    bind("cl_y",       ev.cl_y);
    bind("cl_z",       ev.cl_z);
    bind("cl_energy",  ev.cl_energy);
    bind("cl_nblocks", ev.cl_nblocks);
    bind("cl_center",  ev.cl_center);
    bind("cl_flag",    ev.cl_flag);

    s.has_per_cl_match = (tree->GetBranch("matchFlag") != nullptr);
    if (s.has_per_cl_match) {
        bind("matchFlag", ev.matchFlag);
        bind("matchGEMx", ev.matchGEMx);
        bind("matchGEMy", ev.matchGEMy);
        bind("matchGEMz", ev.matchGEMz);
    }

    s.has_match_num = (tree->GetBranch("match_num") != nullptr);
    if (s.has_match_num) {
        bind("match_num", &ev.matchNum);
        bind("mHit_E",  ev.mHit_E);
        bind("mHit_x",  ev.mHit_x);
        bind("mHit_y",  ev.mHit_y);
        bind("mHit_z",  ev.mHit_z);
        bind("mHit_gx", ev.mHit_gx);
        bind("mHit_gy", ev.mHit_gy);
        bind("mHit_gz", ev.mHit_gz);
        bind("mHit_gid", ev.mHit_gid);
    }

    bind("n_gem_hits",   &ev.n_gem_hits);
    bind("det_id",       ev.det_id);
    bind("gem_x",        ev.gem_x);
    bind("gem_y",        ev.gem_y);
    bind("gem_z",        ev.gem_z);
    bind("gem_x_charge", ev.gem_x_charge);
    bind("gem_y_charge", ev.gem_y_charge);
    bind("gem_x_peak",   ev.gem_x_peak);
    bind("gem_y_peak",   ev.gem_y_peak);
    bind("gem_x_size",   ev.gem_x_size);
    bind("gem_y_size",   ev.gem_y_size);
    bind("gem_x_mTbin",  ev.gem_x_mTbin);
    bind("gem_y_mTbin",  ev.gem_y_mTbin);

    s.has_veto = (tree->GetBranch("veto_nch") != nullptr);
    if (s.has_veto) {
        bind("veto_nch",          &ev.veto_nch);
        bind("veto_id",           ev.veto_id);
        bind("veto_npeaks",       ev.veto_npeaks);
        bind("veto_peak_time",    ev.veto_peak_time);
        bind("veto_peak_integral", ev.veto_peak_integral);
    }

    s.has_lms = (tree->GetBranch("lms_nch") != nullptr);
    if (s.has_lms) {
        bind("lms_nch",          &ev.lms_nch);
        bind("lms_id",           ev.lms_id);
        bind("lms_npeaks",       ev.lms_npeaks);
        bind("lms_peak_time",    ev.lms_peak_time);
        bind("lms_peak_integral", ev.lms_peak_integral);
    }

    // ssp_raw — see note in SetRawReadBranches.
    s.has_ssp_raw = (tree->GetBranch("ssp_raw") != nullptr);

    return s;
}

} // namespace prad2
