// =============================================================================
// live_charge.cpp — accumulate live charge from replayed ROOT files OR raw
// split EVIO files (auto-detected per input).
//
// Inputs (mixed paths allowed):
//   - replayed ROOT files (`scalers` + `epics` side trees from
//     replay_rawdata / replay_recon / replay_filter), OR
//   - raw EVIO files (split-file sequences or directories thereof).
//
// In both cases the integrator does the same arithmetic replay_filter
// applies internally to its passing checkpoint pairs:
//
//     Q = Σ live_fraction · Δt · ½(I_i + I_{i+1})
//
// where live_fraction is the slice-local DSC2 livetime (Δgated / Δungated
// since the previous scaler readout), Δt is the gap between adjacent merged
// checkpoints (TI ticks via the row-stamped `ti_ticks` and
// `ti_ticks_at_arrival`), and I is forward-filled from the configured
// EPICS beam-current channel.
//
// `good` filtering: when ROOT inputs carry the per-row `good` bool that
// replay_filter writes, only passing-passing pairs contribute to the
// canonical live charge (`value_nC` / `live_seconds` / `real_seconds`),
// which is the post-cut quantity downstream tools should use.  An
// ungated counterpart (`ungated_value_nC` / `ungated_live_seconds` /
// `ungated_real_seconds`) is also computed over every valid-data pair
// and reported alongside, so consumers can see how much charge the cut
// throws away.  When the `good` column is absent (raw / recon / EVIO),
// every pair counts as good and the gated values equal the ungated
// values.  Mixing replay_filter ROOT outputs with raw EVIO is allowed
// but probably not what you want — the EVIO rows always pass and will
// inflate the "kept" count.
//
// Beam current is assumed to publish in nA (true for the Hall B IPM
// scalers); the resulting charge is therefore in nC.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>

#include <TFile.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include "EventData.h"
#include "EventData_io.h"

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "InstallPaths.h"

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;
using json   = nlohmann::json;

namespace {

constexpr double TI_TICK_SEC = 4.0e-9;
constexpr int    DSC_NCH     = 16;

// ── Per-input dispatch ──────────────────────────────────────────────────────

enum class InputKind { Root, Evio, Unknown };

InputKind classify(const std::string &p)
{
    const std::string fname = fs::path(p).filename().string();
    if (fname.size() >= 5 && fname.substr(fname.size() - 5) == ".root")
        return InputKind::Root;
    if (fname.find(".evio") != std::string::npos)
        return InputKind::Evio;
    return InputKind::Unknown;
}

// Expand a user-supplied path: file paths pass through, directories are
// listed and each child is classified individually so a single directory can
// contribute both ROOT and EVIO files (rare, but cheap to support).
struct InputFile {
    std::string path;
    InputKind   kind;
};

std::vector<InputFile> collect_inputs(const std::string &arg)
{
    std::vector<InputFile> out;
    auto push = [&](const std::string &p) {
        InputKind k = classify(p);
        if (k != InputKind::Unknown) out.push_back({p, k});
    };
    if (fs::is_directory(arg)) {
        std::vector<std::string> children;
        for (auto &entry : fs::directory_iterator(arg))
            if (entry.is_regular_file())
                children.push_back(entry.path().string());
        std::sort(children.begin(), children.end());
        for (auto &c : children) push(c);
    } else {
        push(arg);
    }
    return out;
}

// ── In-memory rows ──────────────────────────────────────────────────────────

struct ScalerRow {
    int32_t  event_number   = 0;
    int64_t  ti_ticks       = 0;
    uint32_t ref_gated      = 0;
    uint32_t ref_ungated    = 0;
    uint32_t trg_gated[DSC_NCH]   = {};
    uint32_t trg_ungated[DSC_NCH] = {};
    uint32_t tdc_gated[DSC_NCH]   = {};
    uint32_t tdc_ungated[DSC_NCH] = {};
    bool     good           = true;        // defaults to "kept" when no col
};

struct EpicsRow {
    int32_t                       event_number = 0;
    int64_t                       ti_ticks     = 0;
    bool                          good         = true;
    std::map<std::string, double> updates;
};

// ── ROOT-tree readers (replayed `scalers` / `epics` side trees) ─────────────

bool load_scalers_root(const std::vector<std::string> &files,
                       std::vector<ScalerRow>          &out,
                       bool                            &any_good_col)
{
    prad2::RawScalerData sc;
    for (const auto &path : files) {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "live_charge: cannot open " << path << "\n";
            return false;
        }
        TTree *t = dynamic_cast<TTree *>(f->Get("scalers"));
        if (!t) continue;
        prad2::SetScalerReadBranches(t, sc);
        bool good = true;
        const bool has_good = (t->GetBranch("good") != nullptr);
        any_good_col = any_good_col || has_good;
        if (has_good) t->SetBranchAddress("good", &good);
        Long64_t n = t->GetEntries();
        out.reserve(out.size() + n);
        for (Long64_t i = 0; i < n; ++i) {
            good = true;
            t->GetEntry(i);
            ScalerRow r;
            r.event_number = sc.event_number;
            r.ti_ticks     = sc.ti_ticks;
            r.ref_gated    = sc.ref_gated;
            r.ref_ungated  = sc.ref_ungated;
            std::memcpy(r.trg_gated,   sc.trg_gated,   DSC_NCH * sizeof(uint32_t));
            std::memcpy(r.trg_ungated, sc.trg_ungated, DSC_NCH * sizeof(uint32_t));
            std::memcpy(r.tdc_gated,   sc.tdc_gated,   DSC_NCH * sizeof(uint32_t));
            std::memcpy(r.tdc_ungated, sc.tdc_ungated, DSC_NCH * sizeof(uint32_t));
            r.good = good;
            out.push_back(r);
        }
    }
    return true;
}

bool load_epics_root(const std::vector<std::string> &files,
                     std::vector<EpicsRow>           &out,
                     bool                            &any_good_col)
{
    for (const auto &path : files) {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "live_charge: cannot open " << path << "\n";
            return false;
        }
        TTree *t = dynamic_cast<TTree *>(f->Get("epics"));
        if (!t) continue;

        prad2::RawEpicsData ep;
        std::vector<std::string> *cp = &ep.channel;
        std::vector<double>      *vp = &ep.value;
        t->SetBranchAddress("event_number_at_arrival", &ep.event_number_at_arrival);
        const bool has_ticks =
            (t->GetBranch("ti_ticks_at_arrival") != nullptr);
        if (has_ticks)
            t->SetBranchAddress("ti_ticks_at_arrival", &ep.ti_ticks_at_arrival);
        t->SetBranchAddress("channel", &cp);
        t->SetBranchAddress("value",   &vp);
        bool good = true;
        const bool has_good = (t->GetBranch("good") != nullptr);
        any_good_col = any_good_col || has_good;
        if (has_good) t->SetBranchAddress("good", &good);

        Long64_t n = t->GetEntries();
        out.reserve(out.size() + n);
        for (Long64_t i = 0; i < n; ++i) {
            ep.ti_ticks_at_arrival = 0;
            good = true;
            t->GetEntry(i);
            EpicsRow r;
            r.event_number = ep.event_number_at_arrival;
            r.ti_ticks     = has_ticks
                ? static_cast<int64_t>(ep.ti_ticks_at_arrival) : 0;
            r.good         = good;
            const size_t k_max = std::min(ep.channel.size(), ep.value.size());
            for (size_t k = 0; k < k_max; ++k)
                r.updates[ep.channel[k]] = ep.value[k];
            out.push_back(std::move(r));
        }
    }
    return true;
}

// ── EVIO walker (raw split-file sequences) ──────────────────────────────────
//
// One EvChannel reused across every input file — its persistent
// last_physics_event_number_ / last_sync_info_ snapshots survive the
// per-file Close()/OpenAuto() cycle, so EPICS rows arriving at the start of
// file N+1 still get stamped with the last physics event_number from file N
// rather than -1.  Calls only the cheap Info()/Dsc()/Epics() accessors —
// no FADC waveform decode.
//
// At end-of-file we re-run the same `integrate()` arithmetic on just this
// file's row slice so the per-file log line carries elapsed wall-clock
// (first→last physics TI tick), live charge, ⟨livetime⟩ and ⟨I⟩ in addition
// to the raw row counts.  These per-file numbers don't sum exactly to the
// global Q because each slice loses its own first-pair anchor (delta_lt is
// undefined for the leading scaler row), but they're a useful per-file
// signal that's right to a few percent.

// ── Delta-livetime + integration ────────────────────────────────────────────

inline std::pair<uint32_t, uint32_t>
select_pair(const ScalerRow &r, const std::string &source, int channel)
{
    if (source == "ref") return {r.ref_gated, r.ref_ungated};
    const int c = std::clamp(channel, 0, DSC_NCH - 1);
    if (source == "trg") return {r.trg_gated[c], r.trg_ungated[c]};
    if (source == "tdc") return {r.tdc_gated[c], r.tdc_ungated[c]};
    return {0, 0};
}

template <class T>
std::vector<size_t> sort_by_event(const std::vector<T> &v)
{
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) { return v[a].event_number < v[b].event_number; });
    return idx;
}

struct ChargeResult {
    // Gated (post-cut): only adjacent pairs where both endpoints have
    // good=true.  When the input has no `good` column (raw EVIO / recon /
    // unfiltered ROOT) every pair is treated as good, so the gated values
    // equal the ungated values below.  These are the canonical "live
    // charge" numbers downstream tools should consume.
    double  value_nC                   = 0.0;
    double  live_seconds               = 0.0;
    double  real_seconds               = 0.0;   // Σ Δt over integrated pairs
    int64_t n_pairs_kept               = 0;     // both endpoints good=true
    int64_t n_pairs_integrated         = 0;     // contributed to Q
    int64_t n_pairs_skipped            = 0;     // kept but missing data
    // Ungated (pre-cut): every adjacent pair with valid data, ignoring the
    // `good` column.  Reported alongside the gated values so users can see
    // how much charge the cut throws away.  Equal to the gated values when
    // the input has no `good` column.
    double  ungated_value_nC           = 0.0;
    double  ungated_live_seconds       = 0.0;
    double  ungated_real_seconds       = 0.0;
    int64_t n_ungated_pairs_integrated = 0;
    int64_t n_ungated_pairs_skipped    = 0;
    int64_t n_pairs_total              = 0;     // adjacent (i, i+1) pairs walked
    bool    any_good_col               = false; // input carried the `good` bool
};

ChargeResult integrate(const std::vector<ScalerRow> &scalers,
                       const std::vector<EpicsRow>  &epics_rows,
                       bool any_good_col,
                       const std::string &source, int channel,
                       const std::string &beam_current_channel)
{
    ChargeResult r;
    r.any_good_col = any_good_col;

    // 1. delta livetime per scaler row, indexed by load order.
    auto sc_order = sort_by_event(scalers);
    std::vector<double> delta_lt(scalers.size(), -1.0);   // fraction in [0, 1]
    {
        uint32_t pg = 0, pu = 0;
        for (size_t k = 0; k < sc_order.size(); ++k) {
            const size_t orig = sc_order[k];
            const auto [g, u] = select_pair(scalers[orig], source, channel);
            if (g < pg || u < pu) { pg = pu = 0; }   // counter rebase
            const uint32_t dg = g - pg;
            const uint32_t du = u - pu;
            if (du > 0 && dg <= du)
                delta_lt[orig] = double(dg) / double(du);
            pg = g; pu = u;
        }
    }

    // 2. merged-timeline pass with forward-fill of livetime + beam current.
    auto ep_order = sort_by_event(epics_rows);
    struct Cp {
        int64_t ticks         = 0;
        double  live_fraction = std::nan("");
        double  beam_current  = std::nan("");
        bool    good          = true;
    };
    std::vector<Cp> tl;
    tl.reserve(scalers.size() + epics_rows.size());

    double cur_lt = std::nan("");
    double cur_I  = std::nan("");
    size_t i = 0, j = 0;
    while (i < sc_order.size() || j < ep_order.size()) {
        const bool take_sc =
            (i < sc_order.size()) &&
            (j >= ep_order.size() ||
             scalers[sc_order[i]].event_number <=
             epics_rows[ep_order[j]].event_number);
        Cp cp;
        if (take_sc) {
            const size_t orig = sc_order[i++];
            const auto &s = scalers[orig];
            cp.ticks = s.ti_ticks;
            const double v = delta_lt[orig];
            if (v >= 0.0) cur_lt = v;
            cp.live_fraction = cur_lt;
            cp.beam_current  = cur_I;
            cp.good          = s.good;
        } else {
            const size_t orig = ep_order[j++];
            const auto &e = epics_rows[orig];
            cp.ticks = e.ti_ticks;
            for (const auto &kv : e.updates)
                if (kv.first == beam_current_channel) cur_I = kv.second;
            cp.live_fraction = cur_lt;
            cp.beam_current  = cur_I;
            cp.good          = e.good;
        }
        tl.push_back(cp);
    }

    // 3. integrate.  Walk every adjacent pair once and accumulate into
    //    both buckets:
    //      * ungated_*  — every pair with valid data, ignoring the `good`
    //        column (always populated; equals the gated values when the
    //        input has no `good` column).
    //      * value_nC / live_seconds / real_seconds (gated) — only pairs
    //        where both endpoints had good=true.  When the input lacks a
    //        `good` column every pair is treated as good.
    for (size_t k = 1; k < tl.size(); ++k) {
        ++r.n_pairs_total;
        const auto &a = tl[k - 1];
        const auto &b = tl[k];
        const bool good_pair = !any_good_col || (a.good && b.good);
        if (good_pair) ++r.n_pairs_kept;
        const bool data_ok = !(a.ticks <= 0 || b.ticks <= 0 || b.ticks <= a.ticks
            || !std::isfinite(b.live_fraction) || b.live_fraction < 0
            || !std::isfinite(a.beam_current)  || !std::isfinite(b.beam_current));
        if (!data_ok) {
            if (good_pair) ++r.n_pairs_skipped;
            ++r.n_ungated_pairs_skipped;
            continue;
        }
        const double dt = (b.ticks - a.ticks) * TI_TICK_SEC;
        const double I  = 0.5 * (a.beam_current + b.beam_current);
        const double dQ = b.live_fraction * dt * I;
        const double dL = b.live_fraction * dt;
        r.ungated_value_nC     += dQ;
        r.ungated_live_seconds += dL;
        r.ungated_real_seconds += dt;
        ++r.n_ungated_pairs_integrated;
        if (good_pair) {
            r.value_nC     += dQ;
            r.live_seconds += dL;
            r.real_seconds += dt;
            ++r.n_pairs_integrated;
        }
    }
    return r;
}

bool load_from_evio(const std::vector<std::string> &files,
                    const evc::DaqConfig            &daq_cfg,
                    std::vector<ScalerRow>          &scalers,
                    std::vector<EpicsRow>           &epics_rows,
                    int                              max_events_per_file,
                    const std::string               &source,
                    int                              channel,
                    const std::string               &beam_current_channel)
{
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);

    for (const auto &path : files) {
        if (ch.OpenAuto(path) != evc::status::success) {
            std::cerr << "live_charge: cannot open " << path << "\n";
            return false;
        }
        const size_t sc_before = scalers.size();
        const size_t ep_before = epics_rows.size();
        int      total          = 0;
        uint64_t first_phys_ts  = 0;
        uint64_t last_phys_ts   = 0;
        bool     any_phys       = false;

        while (ch.Read() == evc::status::success) {
            if (max_events_per_file > 0 && total >= max_events_per_file) break;
            if (!ch.Scan()) continue;

            const auto et = ch.GetEventType();

            if (et == evc::EventType::Epics) {
                const auto &rec = ch.Epics();
                if (!rec.present) continue;
                EpicsRow r;
                r.event_number = rec.event_number_at_arrival;
                r.ti_ticks     = static_cast<int64_t>(rec.timestamp_at_arrival);
                const size_t k_max = std::min(rec.channel.size(), rec.value.size());
                for (size_t k = 0; k < k_max; ++k)
                    r.updates[rec.channel[k]] = rec.value[k];
                epics_rows.push_back(std::move(r));
                continue;
            }

            if (et != evc::EventType::Physics) continue;

            ch.SelectEvent(0);
            const auto &info = ch.Info();
            if (info.timestamp > 0) {
                if (!any_phys) { first_phys_ts = info.timestamp; any_phys = true; }
                last_phys_ts = info.timestamp;
            }
            const auto &dsc  = ch.Dsc();
            if (dsc.present) {
                ScalerRow r;
                r.event_number = info.event_number;
                r.ti_ticks     = static_cast<int64_t>(info.timestamp);
                r.ref_gated    = dsc.ref_gated;
                r.ref_ungated  = dsc.ref_ungated;
                std::memcpy(r.trg_gated,   dsc.trg_gated,   DSC_NCH * sizeof(uint32_t));
                std::memcpy(r.trg_ungated, dsc.trg_ungated, DSC_NCH * sizeof(uint32_t));
                std::memcpy(r.tdc_gated,   dsc.tdc_gated,   DSC_NCH * sizeof(uint32_t));
                std::memcpy(r.tdc_ungated, dsc.tdc_ungated, DSC_NCH * sizeof(uint32_t));
                scalers.push_back(r);
            }
            ++total;
        }

        const size_t sc_added = scalers.size() - sc_before;
        const size_t ep_added = epics_rows.size() - ep_before;
        const double elapsed_s = any_phys
            ? double(last_phys_ts - first_phys_ts) * TI_TICK_SEC : 0.0;

        // Per-file slice integration.  Cheap: tens of rows.
        ChargeResult rf{};
        if (sc_added > 0 && ep_added > 0) {
            std::vector<ScalerRow> sc_slice(scalers.begin() + sc_before,
                                            scalers.end());
            std::vector<EpicsRow>  ep_slice(epics_rows.begin() + ep_before,
                                            epics_rows.end());
            rf = integrate(sc_slice, ep_slice, /*any_good_col=*/false,
                           source, channel, beam_current_channel);
        }
        const double avg_lt = rf.real_seconds > 0
            ? rf.live_seconds / rf.real_seconds : 0.0;
        const double avg_I  = rf.live_seconds > 0
            ? rf.value_nC / rf.live_seconds : 0.0;

        std::cerr << "  [evio] " << path
                  << "  events="  << total
                  << "  +scaler=" << sc_added
                  << "  +epics="  << ep_added
                  << std::fixed
                  << "  elapsed=" << std::setprecision(2) << elapsed_s << "s"
                  << "  Q="       << std::setprecision(3) << rf.value_nC << "nC"
                  << "  <lt>="    << std::setprecision(4) << avg_lt
                  << "  <I>="     << std::setprecision(3) << avg_I << "nA"
                  << std::defaultfloat << "\n";
        ch.Close();
    }
    return true;
}

// ── CLI ─────────────────────────────────────────────────────────────────────

void print_usage(const char *argv0)
{
    std::cerr <<
"usage: " << argv0 << " <input.(root|evio*)|dir> [more inputs ...]\n"
"        [-c|--beam-current CHAN]   EPICS channel (default hallb_IPM2C21A_CUR)\n"
"        [-s|--source ref|trg|tdc]  livetime DSC2 source (default ref)\n"
"        [-n|--channel N]           DSC2 channel for trg/tdc (default 0)\n"
"        [-d|--daq-config PATH]     EVIO DAQ config (default <db>/daq_config.json)\n"
"        [-N|--max-events N]        max events per file (EVIO only)\n"
"        [-j|--json PATH]           also write a JSON summary to PATH\n"
"        [-h|--help]\n"
"\n"
"Each input is dispatched by extension: *.root → reads `scalers` and\n"
"`epics` side trees from a replayed file; *.evio* → walks the raw EVIO\n"
"stream via EvChannel and pulls the same DSC2 / EPICS records directly\n"
"(no replay step needed).  Directories are expanded; mixing ROOT and\n"
"EVIO inputs is allowed but rarely useful.\n"
"\n"
"Integrates Σ live_fraction · Δt · ½(Iₐ + I_b) over adjacent slow-event\n"
"checkpoints.  When the ROOT trees carry the per-row `good` bool that\n"
"replay_filter writes, only passing-passing pairs contribute (post-cut\n"
"live charge); otherwise every adjacent pair contributes (total live\n"
"charge over the run).  Beam current is assumed to publish in nA, so\n"
"the result is reported in nC.\n";
}

} // namespace

int main(int argc, char **argv)
{
    std::string beam_current = "hallb_IPM2C21A_CUR";
    std::string source       = "ref";
    int         channel      = 0;
    std::string json_path;
    std::string daq_config;
    int         max_events   = -1;

    static struct option long_opts[] = {
        {"beam-current", required_argument, nullptr, 'c'},
        {"source",       required_argument, nullptr, 's'},
        {"channel",      required_argument, nullptr, 'n'},
        {"daq-config",   required_argument, nullptr, 'd'},
        {"max-events",   required_argument, nullptr, 'N'},
        {"json",         required_argument, nullptr, 'j'},
        {"help",         no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0},
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "c:s:n:d:N:j:h", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'c': beam_current = optarg;            break;
        case 's': source       = optarg;            break;
        case 'n': channel      = std::atoi(optarg); break;
        case 'd': daq_config   = optarg;            break;
        case 'N': max_events   = std::atoi(optarg); break;
        case 'j': json_path    = optarg;            break;
        case 'h': print_usage(argv[0]); return 0;
        default:  print_usage(argv[0]); return 2;
        }
    }

    std::vector<std::string> root_inputs, evio_inputs, all_inputs;
    for (int i = optind; i < argc; ++i) {
        for (auto &in : collect_inputs(argv[i])) {
            all_inputs.push_back(in.path);
            if (in.kind == InputKind::Root) root_inputs.push_back(in.path);
            else                            evio_inputs.push_back(in.path);
        }
    }
    if (all_inputs.empty()) { print_usage(argv[0]); return 2; }

    std::vector<ScalerRow> scalers;
    std::vector<EpicsRow>  epics_rows;
    bool any_good_col = false;

    if (!root_inputs.empty()) {
        if (!load_scalers_root(root_inputs, scalers,    any_good_col)) return 1;
        if (!load_epics_root  (root_inputs, epics_rows, any_good_col)) return 1;
    }

    if (!evio_inputs.empty()) {
        const std::string db_dir = prad2::resolve_data_dir(
            "PRAD2_DATABASE_DIR",
            {"../share/prad2evviewer/database"},
            DATABASE_DIR);
        if (daq_config.empty()) daq_config = db_dir + "/daq_config.json";

        evc::DaqConfig daq_cfg;
        if (!evc::load_daq_config(daq_config, daq_cfg)) {
            std::cerr << "live_charge: failed to load DAQ config "
                      << daq_config << "\n";
            return 1;
        }
        std::cerr << "live_charge: " << evio_inputs.size()
                  << " EVIO file(s) — DAQ config " << daq_config << "\n";
        if (!load_from_evio(evio_inputs, daq_cfg, scalers, epics_rows,
                            max_events, source, channel, beam_current))
            return 1;
    }

    if (scalers.empty()) {
        std::cerr << "live_charge: no scaler rows in input — cannot compute "
                     "livetime.\n";
        return 1;
    }
    if (epics_rows.empty()) {
        std::cerr << "live_charge: no EPICS rows in input — cannot read "
                     "beam current.\n";
        return 1;
    }

    const auto r = integrate(scalers, epics_rows, any_good_col,
                             source, channel, beam_current);

    std::cout
        << "live_charge: " << r.value_nC << " nC"
        << (any_good_col ? " (gated by `good` column)" : "") << "\n"
        << "  real time              : " << r.real_seconds << " s\n"
        << "  live time              : " << r.live_seconds << " s\n"
        << "  ⟨livetime⟩             : "
        << (r.real_seconds > 0 ? r.live_seconds / r.real_seconds : 0.0) << "\n"
        << "  ⟨I⟩                    : "
        << (r.live_seconds > 0 ? r.value_nC / r.live_seconds : 0.0) << " nA\n"
        << "  ungated live charge    : " << r.ungated_value_nC << " nC"
        << (any_good_col ? "" : " (= gated; no `good` column)") << "\n"
        << "  ungated real time      : " << r.ungated_real_seconds << " s\n"
        << "  ungated live time      : " << r.ungated_live_seconds << " s\n"
        << "  ungated ⟨livetime⟩     : "
        << (r.ungated_real_seconds > 0
            ? r.ungated_live_seconds / r.ungated_real_seconds : 0.0) << "\n"
        << "  ungated ⟨I⟩            : "
        << (r.ungated_live_seconds > 0
            ? r.ungated_value_nC / r.ungated_live_seconds : 0.0) << " nA\n"
        << "  beam current channel   : " << beam_current << "\n"
        << "  livetime source        : " << source
        << (source == "ref" ? "" : (" ch " + std::to_string(channel))) << "\n"
        << "  ROOT inputs            : " << root_inputs.size() << "\n"
        << "  EVIO inputs            : " << evio_inputs.size() << "\n"
        << "  scaler rows            : " << scalers.size() << "\n"
        << "  EPICS rows             : " << epics_rows.size() << "\n"
        << "  adjacent pairs (total) : " << r.n_pairs_total << "\n"
        << "  pairs kept             : " << r.n_pairs_kept
        << (any_good_col ? " (good=true on both)" : " (no `good` column — every pair kept)")
        << "\n"
        << "  pairs integrated       : " << r.n_pairs_integrated << "\n"
        << "  pairs skipped          : " << r.n_pairs_skipped
        << " (missing tick / livetime / current on either side)\n"
        << "  ungated pairs integ.   : " << r.n_ungated_pairs_integrated << "\n"
        << "  ungated pairs skipped  : " << r.n_ungated_pairs_skipped
        << " (missing tick / livetime / current on either side)\n";

    if (!json_path.empty()) {
        // value_nC / live_seconds / real_seconds are gated by the `good`
        // column when present (canonical post-cut numbers).  The
        // ungated_* counterparts are computed over every valid-data pair
        // and are emitted unconditionally so consumers always see both —
        // for raw / recon inputs without a `good` column they are equal.
        json j = {
            {"value_nC",                   r.value_nC},
            {"unit",                       "nC"},
            {"beam_current_channel",       beam_current},
            {"beam_current_unit",          "nA"},
            {"livetime_source",            source},
            {"livetime_channel",           channel},
            {"live_seconds",               r.live_seconds},
            {"real_seconds",               r.real_seconds},
            {"average_livetime",           r.real_seconds > 0
                                           ? r.live_seconds / r.real_seconds : 0.0},
            {"average_current_nA",         r.live_seconds > 0
                                           ? r.value_nC / r.live_seconds : 0.0},
            {"ungated_value_nC",           r.ungated_value_nC},
            {"ungated_live_seconds",       r.ungated_live_seconds},
            {"ungated_real_seconds",       r.ungated_real_seconds},
            {"ungated_average_livetime",   r.ungated_real_seconds > 0
                                           ? r.ungated_live_seconds / r.ungated_real_seconds
                                           : 0.0},
            {"ungated_average_current_nA", r.ungated_live_seconds > 0
                                           ? r.ungated_value_nC / r.ungated_live_seconds
                                           : 0.0},
            {"n_root_inputs",              root_inputs.size()},
            {"n_evio_inputs",              evio_inputs.size()},
            {"n_scaler_rows",              scalers.size()},
            {"n_epics_rows",               epics_rows.size()},
            {"n_pairs_total",              r.n_pairs_total},
            {"n_pairs_kept",               r.n_pairs_kept},
            {"n_pairs_integrated",         r.n_pairs_integrated},
            {"n_pairs_skipped",            r.n_pairs_skipped},
            {"n_ungated_pairs_integrated", r.n_ungated_pairs_integrated},
            {"n_ungated_pairs_skipped",    r.n_ungated_pairs_skipped},
            {"good_column_present",        any_good_col},
            {"input_files",                all_inputs},
        };
        std::ofstream of(json_path);
        if (!of) {
            std::cerr << "live_charge: cannot write " << json_path << "\n";
            return 1;
        }
        of << j.dump(2) << "\n";
    }
    return 0;
}
