//============================================================================
// tagger_hycal_correlation.C — tagger TDC → HyCal coincidence study (10 pairs)
//
// Two-phase analysis over one read of the evio file.
//
//   Phase 1 (cache + ΔT fits)
//     A single pass collects per-event tuples (T10R TDC, per-pair E‑side
//     TDCs, W1156 height/integral) into memory, and fills one ΔT histogram
//     per pair.  A Gaussian is fit to each peak to extract (μ_k, σ_k).
//
//   Phase 2 (event-wise cut)
//     For each cached event, ΔT_k = tdc(T10R) − tdc(E_k) is tested against
//     the cut ``|ΔT_k − μ_k| < Nσ · σ_k`` for every pair independently.
//     If ANY pair passes AND the event has W1156 FADC samples, the event
//     is "good" and contributes ONCE to two global histograms:
//         W1156_height      W1156 peak height
//         W1156_integral    W1156 peak integral
//
//   The summary canvas shows all 10 ΔT spectra with their fitted μ/±Nσ
//   bounds, the two W1156 histograms, and a terminal table breaks down
//   how many events each pair selected.
//
// Pair layout (update if the DAQ cabling changes):
//
//   T10R           slot 18, channel  0
//   E49 … E53      slot 18, channels 11…15
//   E54 … E58      slot 19, channels  0…4
//
// Compile with ACLiC after loading rootlogon:
//
//     cd build
//     root -l ../analysis/scripts/rootlogon.C
//     .x ../analysis/scripts/tagger_hycal_correlation.C+( \
//         "/data/stage6/prad_023686/prad_023686.evio.00000", \
//         "tagger_w1156_corr.root", 500000)
//============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "VtpData.h"
#include "TdcData.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace evc;

//-----------------------------------------------------------------------------
// Configuration — adjust if the DAQ layout changes
//-----------------------------------------------------------------------------
namespace {

// Reference channel (T10R).
constexpr int T10R_SLOT = 18;
constexpr int T10R_CH   = 0;

// 10 coincidence pairs: T10R vs Eₓ.  slot + channel identify the E-side hit.
struct PairCfg {
    const char *name;
    int         slot;
    int         channel;
};
constexpr PairCfg PAIRS[] = {
    { "E49", 18, 11 }, { "E50", 18, 12 }, { "E51", 18, 13 },
    { "E52", 18, 14 }, { "E53", 18, 15 },
    { "E54", 19,  0 }, { "E55", 19,  1 }, { "E56", 19,  2 },
    { "E57", 19,  3 }, { "E58", 19,  4 },
};
constexpr int N_PAIRS = sizeof(PAIRS) / sizeof(PAIRS[0]);

// Timing cut: |ΔT − μ| < NSIGMA_CUT · σ.
constexpr double NSIGMA_CUT = 3.0;

// HyCal W1156 DAQ address (database/daq_map.json:
// {"W1156", crate 6, slot 7, channel 3}; crate 6 → ROC 0x8C).
constexpr uint32_t W1156_ROC  = 0x8C;
constexpr int      W1156_SLOT = 7;
constexpr int      W1156_CH   = 3;

// Simple FADC peak finder parameters.
constexpr int PED_WINDOW    = 10;
constexpr int INT_HALFWIDTH = 8;

// ΔT histogram axis (wide enough to include typical cable-delay offsets).
constexpr int    DT_BINS  = 800;
constexpr double DT_RANGE = 4000.0;     // ±100 ns (LSB ≈ 25 ps)

// Fit window half-width around the peak bin.
constexpr double DT_FIT_HALF = 40.0;

// W1156 output histograms.
constexpr int    H_BINS = 200;
constexpr double H_MIN  = 0.0;
constexpr double H_MAX  = 4000.0;
constexpr int    I_BINS = 200;
constexpr double I_MIN  = 0.0;
constexpr double I_MAX  = 40000.0;

//-----------------------------------------------------------------------------
// Helpers
//-----------------------------------------------------------------------------

// Earliest hit (smallest TDC value) for (slot, ch) in this event, or -1.
static int first_tdc(const tdc::TdcEventData &t, int slot, int ch)
{
    int best = -1;
    for (int i = 0; i < t.n_hits; ++i) {
        const auto &h = t.hits[i];
        if ((int)h.slot != slot || (int)h.channel != ch) continue;
        if (best < 0 || (int)h.value < best) best = (int)h.value;
    }
    return best;
}

// Pedestal-and-max peak finder: returns false if no samples.
static bool hycal_peak(const fdec::ChannelData &c,
                       float &height, float &integral)
{
    if (c.nsamples <= PED_WINDOW) return false;
    double ped = 0.0;
    for (int i = 0; i < PED_WINDOW; ++i) ped += c.samples[i];
    ped /= PED_WINDOW;

    int tmax = 0;
    double maxv = (double)c.samples[0] - ped;
    for (int i = 1; i < c.nsamples; ++i) {
        double v = (double)c.samples[i] - ped;
        if (v > maxv) { maxv = v; tmax = i; }
    }
    height = (float)maxv;

    int lo = std::max(0, tmax - INT_HALFWIDTH);
    int hi = std::min<int>(c.nsamples, tmax + INT_HALFWIDTH + 1);
    double sum = 0.0;
    for (int i = lo; i < hi; ++i) sum += (double)c.samples[i] - ped;
    integral = (float)sum;
    return true;
}

// Extract W1156 peak (height, integral). Returns false if absent.
static bool w1156_peak(const fdec::EventData &evt,
                       float &height, float &integral)
{
    const fdec::RocData *roc = evt.findRoc(W1156_ROC);
    if (!roc) return false;
    const fdec::SlotData &s = roc->slots[W1156_SLOT];
    if (!s.present) return false;
    const fdec::ChannelData &c = s.channels[W1156_CH];
    if (c.nsamples <= 0) return false;
    return hycal_peak(c, height, integral);
}

// Per-event record cached during Phase 1.
struct Row {
    int32_t t10r;                 // always ≥ 0 — events without T10R are skipped
    int32_t te[N_PAIRS];          // -1 if that pair's E-side did not fire
    bool    has_w1156;
    float   height;
    float   integral;
};

} // anonymous namespace

//=============================================================================
// Entry point
//=============================================================================

int tagger_hycal_correlation(const char *evio_path,
                             const char *out_path   = "tagger_w1156_corr.root",
                             Long64_t    max_events = 0,
                             const char *daq_config = nullptr)
{
    //---- load DAQ config ----------------------------------------------------
    std::string cfg_path = daq_config ? daq_config : "";
    if (cfg_path.empty()) {
        const char *db = std::getenv("PRAD2_DATABASE_DIR");
        cfg_path = std::string(db ? db : "database") + "/daq_config.json";
    }
    DaqConfig cfg;
    if (!load_daq_config(cfg_path, cfg)) {
        std::cerr << "ERROR: cannot load " << cfg_path << "\n";
        return 1;
    }

    //---- open evio ----------------------------------------------------------
    EvChannel ch;
    ch.SetConfig(cfg);
    if (ch.Open(evio_path) != status::success) {
        std::cerr << "ERROR: cannot open " << evio_path << "\n";
        return 1;
    }
    std::cout << "reading " << evio_path << std::endl;

    //---- single pass over the evio file: cache rows + fill ΔT histograms ----
    TFile out(out_path, "RECREATE");
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1110);

    std::vector<TH1D*> h_dt(N_PAIRS, nullptr);
    for (int k = 0; k < N_PAIRS; ++k) {
        const auto &p = PAIRS[k];
        h_dt[k] = new TH1D(
            TString::Format("dt_T10R_%s", p.name),
            TString::Format("#DeltaT = T10R - %s;"
                            "tdc(T10R) - tdc(%s) [LSB];events", p.name, p.name),
            DT_BINS, -DT_RANGE, DT_RANGE);
    }

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto tdc_ptr   = std::make_unique<tdc::TdcEventData>();
    auto &event    = *event_ptr;
    auto &tdc_evt  = *tdc_ptr;

    std::vector<Row> rows;
    rows.reserve(1u << 20);

    Long64_t n_physics = 0;
    Long64_t n_w1156_total = 0;

    // One-line overwriting progress report.
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();
    auto t_last  = t_start;
    constexpr auto progress_interval = std::chrono::milliseconds(300);
    auto report_progress = [&](bool final_line) {
        auto now = clock::now();
        double elapsed = std::chrono::duration<double>(now - t_start).count();
        double rate    = elapsed > 0 ? (double)n_physics / elapsed : 0.0;
        std::cout << "\r  " << std::setw(10) << n_physics << " events";
        if (max_events > 0) {
            double pct = 100.0 * (double)n_physics / (double)max_events;
            double eta = rate > 0
                ? ((double)max_events - (double)n_physics) / rate : 0.0;
            std::cout << " / " << max_events
                      << "  (" << std::fixed << std::setprecision(1) << pct << "%)"
                      << std::defaultfloat
                      << "  " << std::setprecision(3) << rate / 1e3 << "k/s"
                      << "  ETA " << std::fixed << std::setprecision(0)
                      << (int)eta << "s" << std::defaultfloat;
        } else {
            std::cout << "  " << std::setprecision(3) << rate / 1e3 << "k/s";
        }
        std::cout << "  rows: " << std::setw(8) << rows.size()
                  << "  W1156: " << std::setw(8) << n_w1156_total
                  << "     " << std::flush;
        if (final_line) std::cout << "\n";
    };

    while (ch.Read() == status::success) {
        if (!ch.Scan() || ch.GetEventType() != EventType::Physics) continue;

        const int nsub = ch.GetNEvents();
        for (int i = 0; i < nsub; ++i) {
            ch.DecodeEvent(i, event, nullptr, nullptr, &tdc_evt);

            int t0 = first_tdc(tdc_evt, T10R_SLOT, T10R_CH);
            if (t0 < 0) { ++n_physics; goto tick; }

            // Collect this event's row.
            {
                Row r{};
                r.t10r = t0;
                for (int k = 0; k < N_PAIRS; ++k)
                    r.te[k] = first_tdc(tdc_evt, PAIRS[k].slot, PAIRS[k].channel);

                r.has_w1156 = w1156_peak(event, r.height, r.integral);
                if (r.has_w1156) ++n_w1156_total;

                // Fill ΔT histograms right away to avoid a second memory loop.
                for (int k = 0; k < N_PAIRS; ++k)
                    if (r.te[k] >= 0)
                        h_dt[k]->Fill((double)r.t10r - (double)r.te[k]);

                rows.push_back(r);
            }

            ++n_physics;
        tick:
            auto now = clock::now();
            if (now - t_last >= progress_interval) {
                t_last = now; report_progress(false);
            }
            if (max_events > 0 && n_physics >= max_events) goto done;
        }
    }
done:
    report_progress(true);
    ch.Close();
    double elapsed = std::chrono::duration<double>(
        clock::now() - t_start).count();
    std::cout << "scan done: " << n_physics << " physics events in "
              << std::fixed << std::setprecision(1) << elapsed << " s, "
              << rows.size() << " with T10R, "
              << n_w1156_total << " with W1156\n" << std::defaultfloat;

    if (rows.empty()) {
        std::cerr << "no T10R-gated events found — nothing to plot\n";
        out.Close();
        return 1;
    }

    //---- Phase 1b: fit each ΔT histogram to extract (μ, σ) ------------------
    struct Fit {
        double mu = 0, sigma = 0;
        double coarse_peak = 0;
        int    dt_min = 0, dt_max = 0;
        Long64_t n_dt = 0;
        bool   ok = false;
    };
    std::vector<Fit> fits(N_PAIRS);

    std::cout << "\n=== Phase 1: ΔT fits ===\n"
              << "   pair   coarse_peak     mu[LSB]   sigma[LSB]     "
                 "n_dt_filled\n";
    for (int k = 0; k < N_PAIRS; ++k) {
        const auto &p = PAIRS[k];
        Fit &f = fits[k];
        f.n_dt = (Long64_t)h_dt[k]->GetEntries();
        if (f.n_dt == 0) {
            std::cout << "  T10R-" << p.name << "   (no coincidences)\n";
            continue;
        }

        double peak_x = h_dt[k]->GetXaxis()->GetBinCenter(
                            h_dt[k]->GetMaximumBin());
        double bw     = h_dt[k]->GetXaxis()->GetBinWidth(1);
        f.coarse_peak = peak_x;

        TF1 *gfit = new TF1(TString::Format("gfit_%s", p.name), "gaus",
                            peak_x - DT_FIT_HALF, peak_x + DT_FIT_HALF);
        gfit->SetParameter(1, peak_x);
        gfit->SetParameter(2, 5.0);
        int status = (int)h_dt[k]->Fit(gfit, "RQ", "",
                                       peak_x - DT_FIT_HALF,
                                       peak_x + DT_FIT_HALF);
        if (status == 0) {
            f.mu    = gfit->GetParameter(1);
            f.sigma = std::max(std::fabs(gfit->GetParameter(2)), bw);
            f.ok    = true;
        } else {
            f.mu = peak_x;
            f.sigma = 10.0 * bw;
            f.ok = false;
            std::cerr << "  T10R-" << p.name
                      << ": gaussian fit did not converge,"
                      << " falling back to bin peak=" << f.mu
                      << " / σ=" << f.sigma << "\n";
        }

        h_dt[k]->Write();

        std::cout << "  T10R-" << p.name
                  << "   " << std::setw(10) << f.coarse_peak
                  << "   " << std::setw(9) << f.mu
                  << "   " << std::setw(9) << f.sigma
                  << "   " << std::setw(11) << f.n_dt << "\n";
    }

    //---- Phase 2: event-wise cut, fill the two global W1156 histograms ------
    TH1F *h_w_height = new TH1F(
        "W1156_height",
        "W1156 peak height (event passes any T10R-Eₓ cut);"
        "height [ADC];events",
        H_BINS, H_MIN, H_MAX);
    TH1F *h_w_integ  = new TH1F(
        "W1156_integral",
        "W1156 peak integral (event passes any T10R-Eₓ cut);"
        "integral [ADC#upoint sample];events",
        I_BINS, I_MIN, I_MAX);

    Long64_t n_good = 0;                   // events filling the W1156 histograms
    std::vector<Long64_t> n_pass(N_PAIRS, 0);

    for (const auto &r : rows) {
        bool any_pass = false;
        for (int k = 0; k < N_PAIRS; ++k) {
            if (r.te[k] < 0) continue;
            if (!fits[k].ok) continue;      // skip pairs where the fit failed
            double dt   = (double)r.t10r - (double)r.te[k];
            double half = NSIGMA_CUT * fits[k].sigma;
            if (std::fabs(dt - fits[k].mu) < half) {
                ++n_pass[k];
                any_pass = true;
                // don't break — we want the per-pair counters
            }
        }
        if (any_pass && r.has_w1156) {
            h_w_height->Fill(r.height);
            h_w_integ ->Fill(r.integral);
            ++n_good;
        }
    }
    h_w_height->Write();
    h_w_integ ->Write();

    //---- Summary canvas -----------------------------------------------------
    // 4 rows × 10 cols isn't great; use 5 columns × 3 rows for the 10 ΔT's
    // (top two rows), and put the two W1156 histos on the bottom row.
    // Layout: rows 1-2 = ΔT (5 panels each), row 3 = W1156 height + integral.
    TCanvas *canvas = new TCanvas("summary", "tagger-W1156 correlations",
                                  1700, 1000);
    canvas->Divide(5, 3);

    for (int k = 0; k < N_PAIRS; ++k) {
        canvas->cd(k + 1);              // pads 1..10
        h_dt[k]->Draw();
        if (fits[k].ok) {
            double half = NSIGMA_CUT * fits[k].sigma;
            double y = h_dt[k]->GetMaximum() * 1.05;
            auto *l_lo = new TLine(fits[k].mu - half, 0, fits[k].mu - half, y);
            auto *l_hi = new TLine(fits[k].mu + half, 0, fits[k].mu + half, y);
            l_lo->SetLineColor(kRed); l_lo->SetLineStyle(2); l_lo->Draw("same");
            l_hi->SetLineColor(kRed); l_hi->SetLineStyle(2); l_hi->Draw("same");
        }
    }
    canvas->cd(11 + 1); h_w_height->Draw();   // pad 12
    canvas->cd(11 + 2); h_w_integ ->Draw();   // pad 13
    canvas->Write();
    out.Close();

    //---- Terminal summary ---------------------------------------------------
    std::cout << "\n=== Phase 2: event-wise cut ===\n"
              << "  n_good = events with W1156 that pass AT LEAST ONE pair's "
              << NSIGMA_CUT << "-σ cut\n\n"
              << "   pair       mu[LSB]   sigma[LSB]   cut-half[LSB]   "
                 "n_pass   pass%\n";
    for (int k = 0; k < N_PAIRS; ++k) {
        const auto &f = fits[k];
        double half = NSIGMA_CUT * f.sigma;
        double frac = f.n_dt ? 100.0 * (double)n_pass[k] / f.n_dt : 0.0;
        std::cout << "  T10R-" << PAIRS[k].name
                  << "   " << std::setw(9) << f.mu
                  << "   " << std::setw(9) << f.sigma
                  << "   " << std::setw(9) << half
                  << "   " << std::setw(8) << n_pass[k]
                  << "   " << std::setw(5) << frac << "%\n";
    }
    std::cout << "\n  good events (ANY pair cut + W1156): " << n_good
              << " / " << n_w1156_total << " W1156-present events\n"
              << "\nhistograms written to " << out_path << "\n";
    return 0;
}
