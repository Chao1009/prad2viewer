#pragma once
//=============================================================================
// WaveAnalyzer.h — fast FADC250 waveform analysis (zero allocation)
//
// Merges two approaches:
//   - Triangular-kernel smoothing + local-maxima search (robust peak finding)
//   - Iterative outlier-rejection pedestal estimation
//   - Integration with baseline-crossing boundaries
//   - Peak position corrected from smoothed back to raw sample maximum
//   - Threshold adapts to pedestal noise (N × RMS with absolute floor)
//   - All scratch buffers on the stack — no heap allocation in the hot path
//
// Author: Chao Peng (original), merged/rewritten 2025
//=============================================================================

#include "Fadc250Data.h"
#include <cmath>
#include <algorithm>

namespace fdec
{

struct WaveResult {
    Pedestal ped;
    int      npeaks;
    Peak     peaks[MAX_PEAKS];
    void clear() { ped = {}; npeaks = 0; }
};

struct WaveConfig {
    // ---- core soft-analyzer knobs --------------------------------------
    int      smooth_order  = 2;       // kernel order: 1 = identity, N gives a 2N-1 tap triangular kernel
    float    threshold     = 5.0f;    // peak threshold in pedestal RMS units
    float    min_threshold = 3.0f;    // absolute floor (ADC counts above pedestal)
    float    min_peak_ratio = 0.3f;   // new peak must be ≥ this fraction of a nearby peak
    float    int_tail_ratio = 0.1f;   // stop integration when signal drops below this fraction of peak height
    int      tail_break_n  = 2;       // require N consecutive sub-threshold samples to terminate integration
    int      peak_pileup_gap = 2;     // peaks with integration bounds within this many samples are flagged Q_PEAK_PILED
    int      ped_nsamples  = 30;      // max samples for pedestal window
    float    ped_flatness  = 1.0f;    // max RMS for a "flat" pedestal region
    int      ped_max_iter  = 3;       // outlier rejection iterations
    uint16_t overflow      = 4095;    // overflow ADC value (12-bit)
    float    clk_mhz       = 250.0f;  // clock frequency for time conversion

    // ---- NNLS pile-up deconvolution ------------------------------------
    //
    // Deconvolve() solves a non-negative linear-least-squares problem
    // against per-channel pulse templates to recover individual peak
    // amplitudes when pulses overlap.  Knobs match the JSON layout in
    // database/daq_config.json under fadc250_waveform.analyzer.nnls_deconv.
    //
    // - enabled                     Master switch (false ⇒ Deconvolve()
    //                                returns Q_DECONV_DISABLED immediately).
    // - fallback_to_global_template Caller-side hint: when true, the
    //                                upstream code is expected to pass the
    //                                global-median template (with
    //                                tmpl.is_global=true) for channels
    //                                whose own fit is unusable.  The
    //                                analyzer itself just records the
    //                                outcome via Q_DECONV_FALLBACK_GLOBAL.
    // - apply_to_all_peaks          Caller-side hint: when false the caller
    //                                only invokes Deconvolve() on events
    //                                with at least one Q_PEAK_PILED peak.
    // - tau_*_min_ns / tau_*_max_ns Hard validation bounds; templates with
    //                                τ outside this range yield
    //                                Q_DECONV_BAD_TEMPLATE.
    // - cond_number_max             Conditioning ceiling on M^T M; above
    //                                this we abort with Q_DECONV_SINGULAR
    //                                rather than amplify noise.
    // - pre_samples / post_samples  Per-peak window for integral_dec.
    struct NnlsDeconvConfig {
        bool   enabled                     = false;
        bool   fallback_to_global_template = false;
        bool   apply_to_all_peaks          = false;
        float  tau_r_min_ns =  0.5f;
        float  tau_r_max_ns = 10.0f;
        float  tau_f_min_ns =  2.0f;
        float  tau_f_max_ns = 100.0f;
        float  cond_number_max = 1.0e6f;
        int    pre_samples  =  8;
        int    post_samples = 40;
    };
    NnlsDeconvConfig nnls_deconv;
};

class WaveAnalyzer
{
public:
    explicit WaveAnalyzer(const WaveConfig &cfg = {}) : cfg(cfg) {}

    // Analyze one channel. Fills result in-place, no heap allocation.
    void Analyze(const uint16_t *samples, int nsamples, WaveResult &result) const;

    // Public so callers (and Python bindings) can access just the
    // smoothed buffer for plotting / debugging without re-running the
    // whole analyzer.  Stateless: writes only to `buf`.
    void smooth(const uint16_t *raw, int n, float *buf) const;

    // NNLS pile-up deconvolution.  Caller has already invoked Analyze()
    // and passes the resulting `wres` here.  `tmpl` is the per-channel
    // template (caller does the lookup).  Output goes into `dec_out`;
    // see DeconvOutput in Fadc250Data.h for the state semantics.
    //
    // No-op (state set, arrays zeroed, returns immediately) when:
    //   - cfg.nnls_deconv.enabled == false           Q_DECONV_DISABLED
    //   - wres.npeaks == 0                            Q_DECONV_NOT_RUN
    //   - tmpl.tau_r_ns / tau_f_ns out of cfg range   Q_DECONV_BAD_TEMPLATE
    //   - design matrix conditioning > cfg.cond_number_max
    //                                                 Q_DECONV_SINGULAR
    // Otherwise sets Q_DECONV_APPLIED (or Q_DECONV_FALLBACK_GLOBAL when
    // tmpl.is_global == true) and fills the parallel arrays.
    //
    // Stack-only.  Scratch up to MAX_PEAKS × MAX_SAMPLES floats.
    void Deconvolve(const uint16_t *samples, int nsamples,
                    const WaveResult &wres,
                    const PulseTemplate &tmpl,
                    DeconvOutput &dec_out) const;

    WaveConfig cfg;

private:

    // Estimate pedestal mean/rms/slope/nused on samples [start, start+nped)
    // of the smoothed buffer.  Median+MAD bootstrap then iterative σ-clip;
    // sets Q_PED_NOT_CONVERGED / Q_PED_FLOOR_ACTIVE / Q_PED_TOO_FEW_SAMPLES
    // per the converged state.  Q_PED_OVERFLOW / Q_PED_PULSE_IN_WINDOW /
    // Q_PED_TRAILING_WINDOW are set by Analyze() (they need raw samples
    // and the peak-finding result).
    void findPedestal(const float *buf, int start, int nped, Pedestal &ped) const;

    // local-maxima search on smoothed data, fill peaks
    void findPeaks(const uint16_t *raw, const float *buf, int n,
                   float ped_mean, float ped_rms, float thr,
                   WaveResult &result) const;
};

} // namespace fdec
