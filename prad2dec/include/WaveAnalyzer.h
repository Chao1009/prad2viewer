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
#include <string>

namespace fdec
{

class PulseTemplateStore;   // forward decl — defined in PulseTemplateStore.h

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
    // - enabled                     Master switch for the auto-deconv path
    //                                inside Analyze().  Has no effect on
    //                                the explicit Deconvolve() API used by
    //                                the Python diagnostic script.
    // - template_file               Path to the per-channel template JSON
    //                                (see PulseTemplateStore).  Empty ⇒
    //                                application skips the store load and
    //                                deconv silently stays off.
    // - fallback_to_global_template Use the synthesised global-median
    //                                template for channels whose own fit
    //                                failed gates.  Auto-deconv path passes
    //                                this to PulseTemplateStore::Lookup;
    //                                outcome is recorded as
    //                                Q_DECONV_FALLBACK_GLOBAL on success.
    // - apply_to_all_peaks          When false, auto-deconv only fires on
    //                                events with at least one Q_PEAK_PILED
    //                                peak.  Setting true runs NNLS on
    //                                EVERY channel-event with peaks, even
    //                                isolated ones — order-of-magnitude
    //                                slowdown of the replay/viewer.  Use
    //                                only for debugging or when you
    //                                explicitly need a deconvolved height
    //                                on every peak.
    // - tau_*_min_ns / tau_*_max_ns Hard validation bounds; templates with
    //                                τ outside this range yield
    //                                Q_DECONV_BAD_TEMPLATE.
    // - cond_number_max             Conditioning ceiling on M^T M; above
    //                                this we abort with Q_DECONV_SINGULAR
    //                                rather than amplify noise.
    // - pre_samples / post_samples  Per-peak window for integral_dec.
    struct NnlsDeconvConfig {
        bool        enabled                     = false;
        // Path to the per-channel template JSON (output of
        // fit_pulse_template.py).  Caller resolves relative paths against
        // the database directory before handing to PulseTemplateStore.
        // Empty ⇒ no template file configured ⇒ deconv silently disabled.
        std::string template_file;
        bool        fallback_to_global_template = false;
        bool        apply_to_all_peaks          = false;
        float       tau_r_min_ns =  0.5f;
        float       tau_r_max_ns = 10.0f;
        float       tau_f_min_ns =  2.0f;
        float       tau_f_max_ns = 100.0f;
        float       cond_number_max = 1.0e6f;
        int         pre_samples  =  8;
        int         post_samples = 40;
    };
    NnlsDeconvConfig nnls_deconv;
};

// Thread safety: WaveAnalyzer holds per-instance mutable state
// (channel-key context + bound template store) that Analyze() reads on
// every call.  Construct one per thread; do NOT share an instance across
// threads.  The bound PulseTemplateStore itself IS safe to share — see
// PulseTemplateStore.h.
class WaveAnalyzer
{
public:
    explicit WaveAnalyzer(const WaveConfig &cfg = {}) : cfg(cfg) {}

    // Analyze one channel.  Fills result in-place, no heap allocation.
    //
    // When all of the following hold, Analyze ALSO runs NNLS pile-up
    // deconvolution at the end of the pipeline:
    //   * cfg.nnls_deconv.enabled is true
    //   * SetTemplateStore() was called with a valid store
    //   * SetChannelKey() was called with a non-negative triple
    //   * the store has a usable template for that channel
    //   * the event has at least one Q_PEAK_PILED peak (or
    //     cfg.nnls_deconv.apply_to_all_peaks is true)
    //
    // On a successful deconv each affected peak's height/integral are
    // overwritten in place and Q_PEAK_DECONVOLVED is OR-ed into its
    // quality flag.  Any failure path (no template, bad template,
    // singular system) silently leaves the original peaks untouched —
    // production code keeps running without templates.
    void Analyze(const uint16_t *samples, int nsamples, WaveResult &result) const;

    // Public so callers (and Python bindings) can access just the
    // smoothed buffer for plotting / debugging without re-running the
    // whole analyzer.  Stateless: writes only to `buf`.
    void smooth(const uint16_t *raw, int n, float *buf) const;

    // ---- Deconv setup -------------------------------------------------
    //
    // Bind a PulseTemplateStore for the duration of subsequent Analyze()
    // calls.  Pointer must outlive the analyzer.  Pass nullptr to detach.
    void SetTemplateStore(const PulseTemplateStore *store) { template_store_ = store; }
    const PulseTemplateStore *GetTemplateStore() const     { return template_store_; }

    // Identify which channel the next Analyze() call should look the
    // template up for.  Negative values clear the context (subsequent
    // Analyze() calls run without deconv).  Caller is expected to call
    // this once per channel inside its loop.
    void SetChannelKey(int roc_tag, int slot, int channel)
    {
        ck_roc_  = roc_tag;
        ck_slot_ = slot;
        ck_chan_ = channel;
    }
    void ClearChannelKey() { SetChannelKey(-1, -1, -1); }

    // ---- Per-pulse shape fit (calibration utility) -------------------
    //
    // Output of FitPulseShape().  `ok=false` means the LM solve never
    // accepted a step (no useful fit returned); fields below are
    // undefined in that case except `peak_amp` which is set as soon as
    // the slice is validated.
    struct PulseFitResult {
        bool   ok;
        float  t0_ns;          // template onset, ns from start of slice
        float  tau_r_ns;
        float  tau_f_ns;
        float  peak_amp;       // raw pedsub peak height (ADC, un-normalised)
        float  chi2_per_dof;   // on the normalised pulse
        int    n_iter;
    };

    // Output of FitPulseShapeGamma() — same conventions, but the shape is
    // the 3-parameter Gamma model from Li et al. 2024 (PLOS ONE
    // 10.1371/journal.pone.0313999):
    //
    //     v(t) = (t-t0)^b · exp(-c·(t-t0))   for t > t0
    //
    // physically motivated as the impulse response of a multi-stage
    // CR-RC shaping network (which is what the FADC250 anti-alias
    // front end actually is).  The peak is at t = t0 + b/c.  `b` is a
    // shape exponent allowed to be a real number (not just integer
    // filter order); `c_inv_ns` = 1/c stored as ns instead of 1/ns so
    // it's directly comparable to the two-tau τ_f / τ_r values.
    struct PulseFitGammaResult {
        bool   ok;
        float  t0_ns;
        float  b;              // shape exponent (≈ filter order)
        float  c_inv_ns;       // 1/c, ns — asymptotic decay timescale
        float  peak_amp;
        float  chi2_per_dof;
        int    n_iter;
    };

    // Three-parameter Levenberg-Marquardt fit of the unit-amplitude
    // two-tau model T(t; t0, τ_r, τ_f) / T_max(τ_r, τ_f) to a waveform
    // slice after pedsub + per-pulse peak-height normalisation.  Used by
    // analysis/pyscripts/fit_pulse_template.py to build the per-channel
    // template store; the script aggregates per-pulse results to a
    // median/MAD entry per channel.
    //
    // Static because the fit has no analyzer state — purely a math
    // operation.  Stack-only: scratch sized for MAX_SAMPLES floats.
    //
    // `slice` is a window of `nslice` samples around the peak;
    // `peak_idx_in_slice` is the position of the maximum within that
    // window (caller picks).  `ped` and `ped_rms` come from the same
    // analyzer pass that produced the peak.  `model_err_floor` floors
    // the per-sample σ on the normalised pulse (default 1% of unit
    // peak) so χ²/dof stays comparable across amplitudes.
    static PulseFitResult FitPulseShape(const uint16_t *slice, int nslice,
                                        int peak_idx_in_slice,
                                        float ped, float ped_rms,
                                        float clk_ns,
                                        float model_err_floor = 0.01f);

    // Same fit machinery as FitPulseShape but with the Gamma model (see
    // PulseFitGammaResult for the formula).  Three shape params (t0, b,
    // c); the per-pulse peak-height normalisation is identical.  Caller
    // and downstream Python script use --model {two_tau,gamma} to switch.
    static PulseFitGammaResult FitPulseShapeGamma(const uint16_t *slice,
                                                  int nslice,
                                                  int peak_idx_in_slice,
                                                  float ped, float ped_rms,
                                                  float clk_ns,
                                                  float model_err_floor = 0.01f);

    // Power-user / diagnostic API: explicit NNLS deconvolution against a
    // caller-supplied template.  Used by the Python `apply_pulse_template`
    // script to plot before/after comparisons; NOT called by production
    // code (Replay / viewer) — they get deconv automatically through
    // Analyze().  See DeconvOutput in Fadc250Data.h for state semantics.
    //
    // Always runs given valid inputs (samples non-null, npeaks>0, tmpl in
    // the cfg τ ranges, M^TM well-conditioned) — does NOT consult
    // cfg.nnls_deconv.enabled, so the diagnostic can compute deconv
    // values without flipping the production master switch.  Conditioning
    // and τ-range gates from cfg.nnls_deconv still apply.
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

    // Auto-deconv path called from Analyze().  Looks up template via
    // template_store_ + channel-key context; runs NNLS into a stack
    // DeconvOutput, then overwrites the matching peaks' height /
    // integral and OR-s in Q_PEAK_DECONVOLVED on success.  Silent
    // no-op on every failure path (no store, no key, no template,
    // singular, bad template).
    void applyAutoDeconv(const uint16_t *samples, int nsamples,
                         WaveResult &result) const;

    // ---- Deconv binding (mutable so we can set without losing const
    // semantics on Analyze) ----------------------------------------------
    const PulseTemplateStore *template_store_ = nullptr;
    int ck_roc_  = -1;
    int ck_slot_ = -1;
    int ck_chan_ = -1;
};

} // namespace fdec
