#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"

#include <limits>

using namespace fdec;

// --- triangular-kernel smoothing (your SmoothSpectrum, zero-alloc) ----------
void WaveAnalyzer::smooth(const uint16_t *raw, int n, float *buf) const
{
    int res = cfg.smooth_order;
    if (res <= 1) {
        for (int i = 0; i < n; ++i) buf[i] = raw[i];
        return;
    }
    for (int i = 0; i < n; ++i) {
        float val = raw[i];
        float wsum = 1.0f;
        for (int j = 1; j < res; ++j) {
            if (j > i || i + j >= n) continue;
            float w = 1.0f - j / static_cast<float>(res + 1);
            val  += w * (raw[i - j] + raw[i + j]);
            wsum += 2.0f * w;
        }
        buf[i] = val / wsum;
    }
}

// --- iterative pedestal with median/MAD bootstrap + outlier rejection -------
//
// Median + MAD (×1.4826) seed is robust against ≤50% contamination — a
// previous-event tail or early ringing in the leading window biases the
// simple-mean seed badly, which then loosens the σ-clip band and the
// iteration can converge on a contaminated baseline.  Median-seeded σ-clip
// recovers the right baseline immediately and matches the simple-mean
// behaviour on clean baselines.
void WaveAnalyzer::findPedestal(const float *buf, int start, int nped,
                                Pedestal &ped) const
{
    ped = {};
    if (nped <= 0) return;

    // Copy the window plus original sample indices (needed for slope below,
    // since the survivor set after σ-clip is a subset of the window).
    float scratch[MAX_SAMPLES];
    int   orig_idx[MAX_SAMPLES];
    for (int i = 0; i < nped; ++i) {
        scratch[i]  = buf[start + i];
        orig_idx[i] = start + i;
    }
    int active = nped;

    // ── Median + MAD bootstrap.
    float sorted[MAX_SAMPLES];
    for (int i = 0; i < nped; ++i) sorted[i] = scratch[i];
    std::sort(sorted, sorted + nped);
    float mean = (nped % 2 == 1)
               ? sorted[nped / 2]
               : 0.5f * (sorted[nped / 2 - 1] + sorted[nped / 2]);
    for (int i = 0; i < nped; ++i) sorted[i] = std::abs(scratch[i] - mean);
    std::sort(sorted, sorted + nped);
    const float mad = (nped % 2 == 1)
                    ? sorted[nped / 2]
                    : 0.5f * (sorted[nped / 2 - 1] + sorted[nped / 2]);
    float rms = mad * 1.4826f;        // MAD → σ for normally-distributed noise

    // ── Iterative σ-clip from the robust seed.  scratch / orig_idx track
    // surviving samples in lock-step so we can compute slope on the actual
    // survivor set (not on samples that pass the final band post-hoc).
    bool converged = false;
    for (int iter = 0; iter < cfg.ped_max_iter; ++iter) {
        const float band = std::max(rms, cfg.ped_flatness);
        int count = 0;
        for (int i = 0; i < active; ++i) {
            if (std::abs(scratch[i] - mean) < band) {
                scratch[count]  = scratch[i];
                orig_idx[count] = orig_idx[i];
                ++count;
            }
        }
        if (count == active) { converged = true; break; }
        if (count < 5) {
            ped.quality |= Q_PED_TOO_FEW_SAMPLES;
            active = count;
            break;     // keep prior mean/rms — too few survivors to refit
        }
        active = count;
        float sum = 0, sum2 = 0;
        for (int i = 0; i < active; ++i) { sum += scratch[i]; sum2 += scratch[i] * scratch[i]; }
        mean = sum / active;
        const float var = sum2 / active - mean * mean;
        rms = (var > 0) ? std::sqrt(var) : 0;
    }
    if (!converged && !(ped.quality & Q_PED_TOO_FEW_SAMPLES))
        ped.quality |= Q_PED_NOT_CONVERGED;
    if (rms < cfg.ped_flatness)
        ped.quality |= Q_PED_FLOOR_ACTIVE;

    // ── Linear least-squares slope on the survivors (ADC/sample).  Catches
    // baseline drift / pulse-tail contamination that the σ-clip alone can
    // hide (e.g. a slow tail tilts every sample similarly so none of them
    // register as outliers).
    float slope = 0.0f;
    if (active >= 2) {
        double sx = 0, sy = 0;
        for (int i = 0; i < active; ++i) { sx += orig_idx[i]; sy += scratch[i]; }
        const double xbar = sx / active, ybar = sy / active;
        double sxy = 0, sxx = 0;
        for (int i = 0; i < active; ++i) {
            const double dx = orig_idx[i] - xbar;
            sxy += dx * (scratch[i] - ybar);
            sxx += dx * dx;
        }
        if (sxx > 0) slope = static_cast<float>(sxy / sxx);
    }

    ped.mean  = mean;
    ped.rms   = rms;
    ped.nused = static_cast<uint8_t>(active < 255 ? active : 255);
    ped.slope = slope;
}

// --- local-maxima peak search (your SearchMaxima approach, zero-alloc) ------
void WaveAnalyzer::findPeaks(const uint16_t *raw, const float *buf, int n,
                             float ped_mean, float ped_rms, float thr,
                             WaveResult &result) const
{
    result.npeaks = 0;
    if (n < 3) return;

    // track peak-finding ranges (left/right) separately from integration bounds
    int pk_range[MAX_PEAKS][2];  // [i][0]=left, [i][1]=right

    // Trend: +1 rising, -1 falling, 0 flat.  The flat-tolerance scales with
    // the pedestal RMS — a hardcoded 0.1 ADC threshold treats noise-level
    // wiggles as "rising/falling" on quiet channels (ped_rms ~ 0.5), which
    // splits real plateaus into spurious mini-peaks.  Floor at 0.1 keeps
    // behaviour reasonable on raw integer ADC.
    const float trend_tol = std::max(0.1f, 0.5f * ped_rms);
    auto trend = [trend_tol](float a, float b) -> int {
        const float d = a - b;
        return (std::abs(d) < trend_tol) ? 0 : (d > 0 ? 1 : -1);
    };

    for (int i = 1; i < n - 1 && result.npeaks < MAX_PEAKS; ++i) {
        int tr1 = trend(buf[i], buf[i - 1]);  // +1 if buf[i] > left
        int tr2 = trend(buf[i], buf[i + 1]);  // +1 if buf[i] > right

        // local maximum: higher than (or equal to) both neighbors, with at least one strict
        if (tr1 * tr2 < 0 || (tr1 == 0 && tr2 == 0)) continue;

        // handle flat plateau: if flat on the right side, walk to end of plateau
        // and use the center as the peak position
        int flat_end = i;
        if (tr2 == 0) {
            while (flat_end < n - 1 && trend(buf[flat_end], buf[flat_end + 1]) == 0)
                ++flat_end;
            // plateau must fall on the right to be a real maximum
            if (flat_end >= n - 1 || trend(buf[flat_end], buf[flat_end + 1]) <= 0)
                continue;
        }
        int peak_pos = (i + flat_end) / 2;

        // expand peak range: walk left while rising, walk right while falling/flat
        int left = i, right = flat_end;
        while (left > 0 && trend(buf[left], buf[left - 1]) > 0)
            --left;
        while (right < n - 1 && trend(buf[right], buf[right + 1]) >= 0)
            ++right;

        // estimate local baseline from edges (handles peaks on a slope)
        int span = right - left;
        if (span <= 0) continue;
        float base = (buf[left] * (right - peak_pos) + buf[right] * (peak_pos - left))
                   / static_cast<float>(span);

        // height above local baseline on smoothed data
        float smooth_height = buf[peak_pos] - base;
        if (smooth_height < thr) { i = right; continue; }

        // Height above pedestal.  thr already = max(threshold·rms, min_threshold)
        // ≥ 5·rms in the default config, so a separate "≥ 3·rms" guard is
        // redundant (it can only fail when thr already does).
        float ped_height = buf[peak_pos] - ped_mean;
        if (ped_height < thr) { i = right; continue; }

        // --- integrate: walk outward from peak, stop at baseline or tail cutoff ---
        // Termination requires N = cfg.tail_break_n consecutive sub-threshold
        // samples — a single noise dip in the tail no longer truncates the
        // integral early.  Below-threshold samples seen during a not-yet-
        // confirmed run are held in `pending` and either committed (on
        // recovery) or discarded (when the run reaches N).
        //
        // int_left / int_right are INCLUSIVE bounds: they're only advanced
        // when a sample is actually added to `integral`, so they always
        // point to the outermost above-threshold sample on each side.
        float integral = buf[peak_pos] - ped_mean;
        const float tail_cut = ped_height * cfg.int_tail_ratio;
        const int   N_break  = std::max(1, cfg.tail_break_n);
        int int_left = peak_pos, int_right = peak_pos;

        auto is_below = [&](float v) {
            return v < tail_cut || v < ped_rms || v * ped_height < 0;
        };

        {
            int   below_run = 0;
            float pending   = 0.0f;
            for (int j = peak_pos - 1; j >= left; --j) {
                const float v = buf[j] - ped_mean;
                if (is_below(v)) {
                    ++below_run;
                    pending += v;
                    if (below_run >= N_break) break;
                } else {
                    integral += pending + v;
                    pending = 0.0f;
                    below_run = 0;
                    int_left = j;
                }
            }
        }
        {
            int   below_run = 0;
            float pending   = 0.0f;
            for (int j = peak_pos + 1; j <= right; ++j) {
                const float v = buf[j] - ped_mean;
                if (is_below(v)) {
                    ++below_run;
                    pending += v;
                    if (below_run >= N_break) break;
                } else {
                    integral += pending + v;
                    pending = 0.0f;
                    below_run = 0;
                    int_right = j;
                }
            }
        }

        // --- correct peak position: find max in raw samples near smoothed peak ---
        int raw_pos = peak_pos;
        float raw_height = raw[peak_pos] - ped_mean;
        int search = std::max(1, cfg.smooth_order) + (flat_end - i) / 2;  // widen for plateaus
        for (int j = 1; j <= search; ++j) {
            if (peak_pos - j >= 0) {
                float h = raw[peak_pos - j] - ped_mean;
                if (h > raw_height) { raw_height = h; raw_pos = peak_pos - j; }
            }
            if (peak_pos + j < n) {
                float h = raw[peak_pos + j] - ped_mean;
                if (h > raw_height) { raw_height = h; raw_pos = peak_pos + j; }
            }
        }

        // --- reject if overlapping a previous peak and local height too small ---
        // Use peak-finding range (left/right), not integration bounds, for overlap test.
        // smooth_height is the height above the line connecting left/right edges,
        // i.e., how much this peak rises above the tail it sits on.
        bool rejected = false;
        for (int k = 0; k < result.npeaks; ++k) {
            if (left <= pk_range[k][1] && right >= pk_range[k][0]) {
                if (smooth_height < result.peaks[k].height * cfg.min_peak_ratio) {
                    rejected = true;
                    break;
                }
            }
        }
        if (rejected) { i = right; continue; }

        // --- quadratic peak-time interpolation ---
        // Fit y = a x² + b x + c through the 3 raw samples around raw_pos;
        // the parabola vertex sits at δ = (h[-1] - h[+1]) / (2·(h[-1] - 2·h[0] + h[+1]))
        // relative to raw_pos.  Lifts the time resolution from 4 ns
        // (sample-quantised) to ≪ 1 ns for clean peaks.  Guarded by:
        //   - raw_pos not at the buffer edge,
        //   - denom < 0 (real concave-down max — flat plateaus and numerical
        //     noise have denom ≥ 0 and skip interpolation),
        //   - δ clamped to ±1 sample for robustness.
        float t_subsample = 0.0f;
        if (raw_pos > 0 && raw_pos < n - 1) {
            const float h_minus = raw[raw_pos - 1];
            const float h_zero  = raw[raw_pos];
            const float h_plus  = raw[raw_pos + 1];
            const float denom = h_minus - 2.0f * h_zero + h_plus;
            if (denom < -1e-3f) {
                const float delta = 0.5f * (h_minus - h_plus) / denom;
                t_subsample = std::max(-1.0f, std::min(1.0f, delta));
            }
        }

        // --- pile-up detection ---
        // Flag this peak (and the matching previously-found peak) when
        // their integration windows touch or overlap within
        // cfg.peak_pileup_gap samples — diagnostic for downstream cuts on
        // isolated vs piled-up pulses.
        uint8_t my_quality = Q_PEAK_GOOD;
        const int gap = std::max(1, cfg.peak_pileup_gap);
        for (int k = 0; k < result.npeaks; ++k) {
            const Peak &prev = result.peaks[k];
            if (int_left  <= prev.right + gap &&
                int_right >= prev.left  - gap) {
                result.peaks[k].quality |= Q_PEAK_PILED;
                my_quality |= Q_PEAK_PILED;
            }
        }

        // --- fill peak ---
        Peak &p = result.peaks[result.npeaks];
        p.pos      = raw_pos;
        p.left     = int_left;
        p.right    = int_right;
        p.height   = raw_height;
        p.integral = integral;
        p.time     = (raw_pos + t_subsample) * 1e3f / cfg.clk_mhz;  // ns
        p.overflow = (raw[raw_pos] >= cfg.overflow);
        p.quality  = my_quality;
        pk_range[result.npeaks][0] = left;
        pk_range[result.npeaks][1] = right;
        result.npeaks++;

        // skip past this peak's range to avoid double-counting
        i = right;
    }
}

// --- main entry point -------------------------------------------------------
void WaveAnalyzer::Analyze(const uint16_t *samples, int nsamples, WaveResult &result) const
{
    result.clear();
    if (!samples || nsamples <= 0 || nsamples > MAX_SAMPLES) return;

    // stack-allocated scratch buffer for smoothed waveform
    float buf[MAX_SAMPLES];
    smooth(samples, nsamples, buf);

    auto window_overflow = [&](int wstart, int wlen) -> bool {
        const uint16_t ovr = cfg.overflow;
        for (int i = wstart; i < wstart + wlen; ++i)
            if (samples[i] >= ovr) return true;
        return false;
    };

    const int nped_window = std::min(cfg.ped_nsamples, nsamples);

    // ── Leading-window pedestal estimate.
    Pedestal P_lead;
    findPedestal(buf, 0, nped_window, P_lead);
    if (window_overflow(0, nped_window))
        P_lead.quality |= Q_PED_OVERFLOW;

    // ── Adaptive: if the leading window looks suspicious (didn't converge,
    // lost > 50% of samples to rejection, or hit overflow), try the
    // trailing window — only if the two don't overlap.  Pick whichever
    // has the lower RMS (with nused as tiebreaker); flag the choice with
    // Q_PED_TRAILING_WINDOW.
    Pedestal P_use         = P_lead;
    int      ped_win_start = 0;
    const bool lead_suspicious =
        (P_lead.quality & (Q_PED_NOT_CONVERGED |
                           Q_PED_TOO_FEW_SAMPLES |
                           Q_PED_OVERFLOW))
        || (P_lead.nused * 2 < nped_window);

    if (lead_suspicious && nsamples >= 2 * nped_window) {
        const int trail_start = nsamples - nped_window;
        Pedestal P_trail;
        findPedestal(buf, trail_start, nped_window, P_trail);
        if (window_overflow(trail_start, nped_window))
            P_trail.quality |= Q_PED_OVERFLOW;
        const bool trail_better =
            (P_trail.rms < P_lead.rms) ||
            (P_trail.rms == P_lead.rms && P_trail.nused > P_lead.nused);
        if (trail_better) {
            P_use         = P_trail;
            P_use.quality |= Q_PED_TRAILING_WINDOW;
            ped_win_start = trail_start;
        }
    }
    result.ped = P_use;

    // ── Peak finding uses the chosen pedestal.
    const float thr = std::max(cfg.threshold * result.ped.rms, cfg.min_threshold);
    findPeaks(samples, buf, nsamples, result.ped.mean, result.ped.rms, thr, result);

    // ── Post-hoc: was a real pulse inside the pedestal window we used?
    // Diagnostic for downstream filters — doesn't influence the estimate
    // (the median+MAD seed already absorbs single-pulse contamination on
    // most channels), but lets analyses optionally cut on clean events.
    const int ped_win_end = ped_win_start + nped_window;
    for (int p = 0; p < result.npeaks; ++p) {
        const int pos = result.peaks[p].pos;
        if (pos >= ped_win_start && pos < ped_win_end) {
            result.ped.quality |= Q_PED_PULSE_IN_WINDOW;
            break;
        }
    }

    // ── NNLS pile-up deconv.  Silent no-op unless the caller bound a
    // template store and the current channel key — production code that
    // doesn't care about deconv keeps its existing behaviour bit-for-bit.
    applyAutoDeconv(samples, nsamples, result);
}

void WaveAnalyzer::applyAutoDeconv(const uint16_t *samples, int nsamples,
                                   WaveResult &result) const
{
    if (!cfg.nnls_deconv.enabled)               return;
    if (template_store_ == nullptr)             return;
    if (ck_roc_ < 0 || ck_slot_ < 0 || ck_chan_ < 0) return;
    if (result.npeaks <= 0)                     return;

    // Cheap gate: skip clean events unless config asks otherwise.  This
    // is the dominant cost saving — most channels see no pile-up and
    // we'd otherwise pay an NNLS solve on every event.
    if (!cfg.nnls_deconv.apply_to_all_peaks) {
        bool any_piled = false;
        for (int k = 0; k < result.npeaks; ++k) {
            if (result.peaks[k].quality & Q_PEAK_PILED) {
                any_piled = true; break;
            }
        }
        if (!any_piled) return;
    }

    const PulseTemplate *tmpl = template_store_->Lookup(
        ck_roc_, ck_slot_, ck_chan_,
        cfg.nnls_deconv.fallback_to_global_template);
    if (tmpl == nullptr) return;

    DeconvOutput out;
    Deconvolve(samples, nsamples, result, *tmpl, out);

    // On success: replace each peak's height/integral with the deconv
    // values and mark Q_PEAK_DECONVOLVED.  Failure paths
    // (Q_DECONV_BAD_TEMPLATE / Q_DECONV_SINGULAR) leave the peaks as-is
    // so downstream code falls back to WaveAnalyzer's tail-cutoff values.
    if (out.state == Q_DECONV_APPLIED || out.state == Q_DECONV_FALLBACK_GLOBAL) {
        const int K = (out.n < result.npeaks) ? out.n : result.npeaks;
        for (int k = 0; k < K; ++k) {
            result.peaks[k].height    = out.height[k];
            result.peaks[k].integral  = out.integral[k];
            result.peaks[k].quality  |= Q_PEAK_DECONVOLVED;
        }
    }
}

//=============================================================================
// Per-pulse-fit pile-up deconvolution
//=============================================================================
//
// Given pedsub samples b[0..n-1] and K peak times τ_1..τ_K from the
// WaveAnalyzer, we fit a 4K-parameter model
//
//   model(t_i) = Σ_k a_k · T(t_i; t0_k, τ_r_k, τ_f_k) / T_max(τ_r_k, τ_f_k)
//
// via Levenberg-Marquardt with the channel template providing the
// initial guess (a_k = WA peak.height / T_max(template), t0_k from
// peak.time, (τ_r, τ_f) = template) and tight bounds around it
// (cfg.nnls_deconv.shape_window_factor, t0_window_ns, amp_max_factor).
// This handles per-pulse shape variation that a fixed-shape NNLS
// can't capture.
//
// Algorithm: classical LM with forward-difference Jacobian.  The
// shape-param FD is per-peak (each FD perturbation only recomputes
// the k-th template column, not the whole model — analytic block
// structure of the Jacobian).

namespace {

// Lower-triangular Cholesky factorisation of a KxK SPD matrix M (row-major,
// L is also row-major).  Returns the smallest pivot squared (= L[k,k]²) so
// callers can do a conditioning check.  Returns -1 on failure (M not SPD).
inline float cholesky(const float *M, float *L, int K)
{
    float min_pivot_sq = std::numeric_limits<float>::infinity();
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j <= i; ++j) {
            float sum = M[i * K + j];
            for (int k = 0; k < j; ++k)
                sum -= L[i * K + k] * L[j * K + k];
            if (i == j) {
                if (sum <= 0.0f) return -1.0f;
                if (sum < min_pivot_sq) min_pivot_sq = sum;
                L[i * K + j] = std::sqrt(sum);
            } else {
                L[i * K + j] = sum / L[j * K + j];
            }
        }
        // Zero the strict upper triangle so chol_solve is well-defined
        for (int j = i + 1; j < K; ++j) L[i * K + j] = 0.0f;
    }
    return min_pivot_sq;
}

// Solve L L^T x = b using a precomputed Cholesky factor.  Sized to
// accommodate the per-pulse-fit deconvolver's 4·MAX_PEAKS-parameter
// Hessian (so K can be up to 32, not just MAX_PEAKS=8).
inline void chol_solve(const float *L, const float *b, float *x, int K)
{
    float y[4 * MAX_PEAKS];   // = 32; covers K up to 4·MAX_PEAKS
    // Forward: L y = b
    for (int i = 0; i < K; ++i) {
        float s = b[i];
        for (int k = 0; k < i; ++k) s -= L[i * K + k] * y[k];
        y[i] = s / L[i * K + i];
    }
    // Back: L^T x = y
    for (int i = K - 1; i >= 0; --i) {
        float s = y[i];
        for (int k = i + 1; k < K; ++k) s -= L[k * K + i] * x[k];
        x[i] = s / L[i * K + i];
    }
}

// Closed-form template peak position offset (ns) and peak value.
//   t_peak = τ_r · ln((τ_r + τ_f) / τ_r)
//   T_max  = (τ_f / (τ_r + τ_f)) · (τ_r / (τ_r + τ_f))^(τ_r/τ_f)
inline void template_peak(float tr, float tf, float &t_off, float &t_max)
{
    const float u = tr / (tr + tf);
    t_off = tr * std::log(1.0f / u);
    t_max = (1.0f - u) * std::pow(u, tr / tf);
}

// Unit-amplitude template at sample times t_i = i·clk_ns, with template
// onset at t0 (ns).  Writes n values into `col`.  Two exp() per sample.
// In the per-pulse LM the shape parameters (τ_r, τ_f) vary per peak per
// iteration, so the precomputed-grid optimisation that the old fixed-
// shape NNLS used doesn't apply — every iteration computes fresh
// templates analytically.
inline void template_column(float *col, int n, float clk_ns,
                            float t0, float tr, float tf)
{
    for (int i = 0; i < n; ++i) {
        const float t = i * clk_ns;
        if (t <= t0) { col[i] = 0.0f; continue; }
        const float dt = t - t0;
        col[i] = (1.0f - std::exp(-dt / tr)) * std::exp(-dt / tf);
    }
}

} // namespace (anonymous)

void WaveAnalyzer::Deconvolve(const uint16_t *samples, int nsamples,
                              const WaveResult &wres,
                              const PulseTemplate &tmpl,
                              DeconvOutput &dec_out) const
{
    dec_out.clear();

    // Explicit API: always runs when given valid inputs and a usable
    // template.  The cfg.nnls_deconv.enabled gate only governs the auto
    // path (applyAutoDeconv inside Analyze) so the Python diagnostic
    // can compute deconv values without flipping the production switch.
    if (!samples || nsamples <= 0 || nsamples > MAX_SAMPLES) {
        dec_out.state = Q_DECONV_NOT_RUN;
        return;
    }
    const int K = wres.npeaks;
    if (K <= 0 || K > MAX_PEAKS) {
        dec_out.state = Q_DECONV_NOT_RUN;
        return;
    }

    const auto &dcfg = cfg.nnls_deconv;
    const float tr_tmpl = tmpl.tau_r_ns;
    const float tf_tmpl = tmpl.tau_f_ns;
    if (!(tr_tmpl >= dcfg.tau_r_min_ns && tr_tmpl <= dcfg.tau_r_max_ns) ||
        !(tf_tmpl >= dcfg.tau_f_min_ns && tf_tmpl <= dcfg.tau_f_max_ns) ||
        !(tr_tmpl > 0.0f) || !(tf_tmpl > 0.0f)) {
        dec_out.state = Q_DECONV_BAD_TEMPLATE;
        return;
    }

    const float clk_ns = (cfg.clk_mhz > 0.0f) ? (1000.0f / cfg.clk_mhz) : 4.0f;
    const float ped    = wres.ped.mean;
    const float sigma  = std::max(wres.ped.rms, 1.0f);
    const float inv_sigma2 = 1.0f / (sigma * sigma);

    // Pedsub data.
    float b[MAX_SAMPLES];
    for (int i = 0; i < nsamples; ++i)
        b[i] = static_cast<float>(samples[i]) - ped;

    // Initial param vector laid out as 4K floats with the per-peak block
    //   [a_k, t0_k, τ_r_k, τ_f_k]
    // at index 4*k.  Bounds in p_lo / p_hi follow the same layout.
    constexpr int NPP = 4;                          // params per peak
    constexpr int NPMAX = NPP * MAX_PEAKS;          // 32
    const int npar = NPP * K;

    float t_off_tmpl, t_max_tmpl;
    template_peak(tr_tmpl, tf_tmpl, t_off_tmpl, t_max_tmpl);
    const float inv_t_max_tmpl = 1.0f / t_max_tmpl;

    const float fac = std::max(dcfg.shape_window_factor, 1.001f);
    const float t0_win = std::max(dcfg.t0_window_ns, 0.0f);

    float p[NPMAX]   = {0};
    float p_lo[NPMAX] = {0};
    float p_hi[NPMAX] = {0};
    float peak_t[MAX_PEAKS];

    for (int k = 0; k < K; ++k) {
        const float t_pk = wres.peaks[k].time;           // ns from sample 0
        peak_t[k] = t_pk;
        const float h_wa = std::max(wres.peaks[k].height, 1.0f);
        const float a_init  = h_wa * inv_t_max_tmpl;     // height / T_max
        const float t0_init = t_pk - t_off_tmpl;

        p[NPP*k + 0] = a_init;
        p[NPP*k + 1] = t0_init;
        p[NPP*k + 2] = tr_tmpl;
        p[NPP*k + 3] = tf_tmpl;

        p_lo[NPP*k + 0] = 0.0f;
        p_hi[NPP*k + 0] = std::max(dcfg.amp_max_factor, 1.001f) * h_wa
                         * inv_t_max_tmpl;
        p_lo[NPP*k + 1] = t0_init - t0_win;
        p_hi[NPP*k + 1] = t0_init + t0_win;
        p_lo[NPP*k + 2] = tr_tmpl / fac;
        p_hi[NPP*k + 2] = tr_tmpl * fac;
        p_lo[NPP*k + 3] = tf_tmpl / fac;
        p_hi[NPP*k + 3] = tf_tmpl * fac;
    }

    // Helpers to (re)compute the K template columns at given params and
    // evaluate the residual chi².
    auto compute_M = [&](const float *params, float *M) {
        for (int k = 0; k < K; ++k) {
            template_column(&M[k * MAX_SAMPLES], nsamples, clk_ns,
                            params[NPP*k + 1],
                            params[NPP*k + 2],
                            params[NPP*k + 3]);
        }
    };
    auto eval_chi2 = [&](const float *params, const float *M) -> float {
        float c = 0.0f;
        for (int i = 0; i < nsamples; ++i) {
            float fit = 0.0f;
            for (int k = 0; k < K; ++k)
                fit += params[NPP*k + 0] * M[k * MAX_SAMPLES + i];
            const float rr = (b[i] - fit);
            c += rr * rr;
        }
        return c * inv_sigma2;
    };

    // Initial M and chi².
    float M[MAX_PEAKS * MAX_SAMPLES];
    compute_M(p, M);
    float chi2 = eval_chi2(p, M);

    // Best-seen tracking (return whatever we found, scipy-style).
    float chi2_best = chi2;
    float p_best[NPMAX];
    for (int j = 0; j < npar; ++j) p_best[j] = p[j];

    // LM hyperparameters.
    constexpr int   MAX_ITER  = 100;
    constexpr float TOL_PARAM = 1e-5f;
    constexpr float LAMBDA0   = 1.0e-3f;
    constexpr float LAMBDA_UP = 10.0f;
    constexpr float LAMBDA_DN = 10.0f;
    constexpr float LAMBDA_MAX = 1.0e10f;

    float lambda = LAMBDA0;
    int iter = 0;
    bool any_accepted = false;

    // Scratch buffers reused across iterations.
    float r_vec[MAX_SAMPLES];
    float J[NPMAX * MAX_SAMPLES];   // J[j * MAX_SAMPLES + i] = ∂r_i/∂p_j
    float A[NPMAX * NPMAX];
    float Aug[NPMAX * NPMAX];
    float L[NPMAX * NPMAX];
    float g[NPMAX], delta[NPMAX];
    float p_new[NPMAX];
    float M_new[MAX_PEAKS * MAX_SAMPLES];
    float col_perturbed[MAX_SAMPLES];

    for (; iter < MAX_ITER; ++iter) {
        // Residuals at current point.
        for (int i = 0; i < nsamples; ++i) {
            float fit = 0.0f;
            for (int k = 0; k < K; ++k)
                fit += p[NPP*k + 0] * M[k * MAX_SAMPLES + i];
            r_vec[i] = b[i] - fit;
        }

        // Jacobian — block structure: only the k-th template column
        // depends on (t0_k, τ_r_k, τ_f_k); all peaks contribute to a_k
        // through the analytic ∂model/∂a_k = -T_k.
        for (int k = 0; k < K; ++k) {
            const float a_k  = p[NPP*k + 0];
            const float t0_k = p[NPP*k + 1];
            const float tr_k = p[NPP*k + 2];
            const float tf_k = p[NPP*k + 3];

            // ∂r/∂a_k = -T_k(t_i)  (analytic — reuse M).
            for (int i = 0; i < nsamples; ++i)
                J[(NPP*k + 0) * MAX_SAMPLES + i] = -M[k * MAX_SAMPLES + i];

            const float h_t0 = std::max(1e-3f * clk_ns, 1e-6f);
            const float h_tr = std::max(1e-3f * tr_k,    1e-6f);
            const float h_tf = std::max(1e-3f * tf_k,    1e-6f);

            // ∂r/∂t0_k = -a_k · ∂T_k/∂t0_k
            template_column(col_perturbed, nsamples, clk_ns,
                            t0_k + h_t0, tr_k, tf_k);
            for (int i = 0; i < nsamples; ++i) {
                const float dT = (col_perturbed[i] - M[k * MAX_SAMPLES + i]) / h_t0;
                J[(NPP*k + 1) * MAX_SAMPLES + i] = -a_k * dT;
            }
            template_column(col_perturbed, nsamples, clk_ns,
                            t0_k, tr_k + h_tr, tf_k);
            for (int i = 0; i < nsamples; ++i) {
                const float dT = (col_perturbed[i] - M[k * MAX_SAMPLES + i]) / h_tr;
                J[(NPP*k + 2) * MAX_SAMPLES + i] = -a_k * dT;
            }
            template_column(col_perturbed, nsamples, clk_ns,
                            t0_k, tr_k, tf_k + h_tf);
            for (int i = 0; i < nsamples; ++i) {
                const float dT = (col_perturbed[i] - M[k * MAX_SAMPLES + i]) / h_tf;
                J[(NPP*k + 3) * MAX_SAMPLES + i] = -a_k * dT;
            }
        }

        // Build A = J^T J (npar × npar) and g = J^T r (npar).
        for (int j = 0; j < npar; ++j) {
            g[j] = 0;
            for (int j2 = 0; j2 < npar; ++j2) A[j * npar + j2] = 0;
        }
        for (int i = 0; i < nsamples; ++i) {
            for (int j1 = 0; j1 < npar; ++j1) {
                const float jj1 = J[j1 * MAX_SAMPLES + i];
                g[j1] += jj1 * r_vec[i];
                for (int j2 = 0; j2 <= j1; ++j2) {
                    const float jj2 = J[j2 * MAX_SAMPLES + i];
                    A[j1 * npar + j2] += jj1 * jj2;
                }
            }
        }
        // Mirror lower → upper.
        for (int j1 = 0; j1 < npar; ++j1)
            for (int j2 = j1 + 1; j2 < npar; ++j2)
                A[j1 * npar + j2] = A[j2 * npar + j1];

        // (A + λI) δ = g
        for (int j1 = 0; j1 < npar; ++j1)
            for (int j2 = 0; j2 < npar; ++j2)
                Aug[j1 * npar + j2] = A[j1 * npar + j2];
        for (int j = 0; j < npar; ++j) Aug[j * npar + j] += lambda;

        const float min_pivot_sq = cholesky(Aug, L, npar);
        if (min_pivot_sq < 0.0f) {
            // Indefinite — back off and try a tighter LM step.
            lambda *= LAMBDA_UP;
            if (lambda > LAMBDA_MAX) break;
            continue;
        }
        chol_solve(L, g, delta, npar);

        // Trial point — same sign convention as FitPulseShape:
        //   p_new = p − δ, then clamp to bounds.
        float dpar = 0.0f;
        for (int j = 0; j < npar; ++j) {
            p_new[j] = std::clamp(p[j] - delta[j], p_lo[j], p_hi[j]);
            const float scale = std::max(std::abs(p[j]), 1e-3f);
            dpar += std::abs(p_new[j] - p[j]) / scale;
        }

        compute_M(p_new, M_new);
        const float chi2_new = eval_chi2(p_new, M_new);

        if (std::isfinite(chi2_new) && chi2_new < chi2_best) {
            chi2_best = chi2_new;
            for (int j = 0; j < npar; ++j) p_best[j] = p_new[j];
        }

        if (chi2_new < chi2) {
            // Accept.
            for (int j = 0; j < npar; ++j) p[j] = p_new[j];
            for (int idx = 0; idx < K * MAX_SAMPLES; ++idx) M[idx] = M_new[idx];
            chi2 = chi2_new;
            lambda = std::max(lambda / LAMBDA_DN, 1e-12f);
            any_accepted = true;
            if (dpar < TOL_PARAM) { ++iter; break; }
        } else {
            lambda *= LAMBDA_UP;
            if (lambda > LAMBDA_MAX) break;
        }
    }

    if (!any_accepted) {
        dec_out.state = Q_DECONV_LM_NOT_CONVERGED;
        return;
    }

    // Use the best-seen params.  Recompute M one more time so the
    // integral window sums use the converged shape.
    for (int j = 0; j < npar; ++j) p[j] = p_best[j];
    compute_M(p, M);

    dec_out.n = K;
    for (int k = 0; k < K; ++k) {
        const float a_k  = p[NPP*k + 0];
        const float t0_k = p[NPP*k + 1];
        const float tr_k = p[NPP*k + 2];
        const float tf_k = p[NPP*k + 3];

        if (!std::isfinite(a_k) || !std::isfinite(t0_k) ||
            !std::isfinite(tr_k) || !std::isfinite(tf_k)) {
            dec_out.state = Q_DECONV_LM_NOT_CONVERGED;
            return;
        }

        float t_off_k, t_max_k;
        template_peak(tr_k, tf_k, t_off_k, t_max_k);

        dec_out.amplitude[k] = a_k;
        dec_out.height[k]    = a_k * t_max_k;
        dec_out.t0_ns[k]     = t0_k;
        dec_out.tau_r_ns[k]  = tr_k;
        dec_out.tau_f_ns[k]  = tf_k;

        // Per-peak integral: window centred on the WaveAnalyzer peak
        // index (consistent with how downstream consumers slice raw data),
        // summed against the per-peak fitted template.
        const int i_pk = static_cast<int>(std::lround(peak_t[k] / clk_ns));
        const int lo   = std::max(0,        i_pk - dcfg.pre_samples);
        const int hi   = std::min(nsamples, i_pk + dcfg.post_samples + 1);
        float sum = 0.0f;
        for (int i = lo; i < hi; ++i) sum += M[k * MAX_SAMPLES + i];
        dec_out.integral[k] = a_k * sum;
    }

    const int dof = std::max(1, nsamples - npar);
    dec_out.chi2_per_dof = chi2_best / static_cast<float>(dof);
    dec_out.state = tmpl.is_global ? Q_DECONV_FALLBACK_GLOBAL
                                   : Q_DECONV_APPLIED;
}

//=============================================================================
// Per-pulse shape fit (Levenberg-Marquardt on normalised two-tau model)
//=============================================================================
//
// Three-parameter LM solve.  Jacobian via central / forward finite
// differences (FD) on the model itself — analytic partials are correct
// but tedious and bring no real speed-up here since the dominant cost
// is the exp() calls inside the model evaluator and FD reuses those.
//
// Replaces the scipy.optimize.curve_fit call in
// analysis/pyscripts/fit_pulse_template.py and gives the script an
// honest ~100× speed-up (5 ev/s → 500 ev/s territory).

namespace {

inline float two_tau_unit_pt(float t, float t0, float tr, float tf,
                              float t_max_inv)
{
    if (t <= t0) return 0.0f;
    const float dt = t - t0;
    return (1.0f - std::exp(-dt / tr)) * std::exp(-dt / tf) * t_max_inv;
}

inline float t_max_value(float tr, float tf)
{
    const float u = tr / (tr + tf);
    return (1.0f - u) * std::pow(u, tr / tf);
}


// 3×3 SPD solve via Cholesky.  Returns false if A is not positive
// definite (caller bumps λ and retries).  A is row-major; b and x are
// 3-vectors.
inline bool chol3_solve(const float A[9], const float b[3], float x[3])
{
    if (A[0] <= 0.0f) return false;
    const float L00 = std::sqrt(A[0]);
    const float L10 = A[3] / L00;
    const float v11 = A[4] - L10 * L10;
    if (v11 <= 0.0f) return false;
    const float L11 = std::sqrt(v11);
    const float L20 = A[6] / L00;
    const float L21 = (A[7] - L20 * L10) / L11;
    const float v22 = A[8] - L20 * L20 - L21 * L21;
    if (v22 <= 0.0f) return false;
    const float L22 = std::sqrt(v22);

    const float y0 = b[0] / L00;
    const float y1 = (b[1] - L10 * y0) / L11;
    const float y2 = (b[2] - L20 * y0 - L21 * y1) / L22;

    x[2] = y2 / L22;
    x[1] = (y1 - L21 * x[2]) / L11;
    x[0] = (y0 - L10 * x[1] - L20 * x[2]) / L00;
    return true;
}

// Sum of squared residuals on the normalised pulse for the given params.
inline float chi2_eval(const float *y, int n, float clk_ns,
                       float t0, float tr, float tf)
{
    const float t_max_inv = 1.0f / t_max_value(tr, tf);
    float s = 0.0f;
    for (int i = 0; i < n; ++i) {
        const float t = i * clk_ns;
        const float r = y[i] - two_tau_unit_pt(t, t0, tr, tf, t_max_inv);
        s += r * r;
    }
    return s;
}

// ---- Two-tau with rise-edge exponent ---------------------------------
//
// T(t; t0, τ_r, τ_f, p) = [1 − exp(−(t−t0)/τ_r)]^p · exp(−(t−t0)/τ_f)
//
// Closed-form peak:
//   u_peak  = p·τ_f / (τ_r + p·τ_f)        (= 1 − exp(−dt_peak/τ_r))
//   dt_peak = τ_r · ln((τ_r + p·τ_f) / τ_r)
//   T_max   = u_peak^p · (τ_r/(τ_r + p·τ_f))^(τ_r/τ_f)
//
// p = 1 reduces to the standard two-tau form (verify: u_peak = τ_f/(τ_r+τ_f),
// matches t_max_value above).
inline float two_tau_p_unit_pt(float t, float t0, float tr, float tf,
                                float p, float t_max_inv)
{
    if (t <= t0) return 0.0f;
    const float dt = t - t0;
    const float u  = 1.0f - std::exp(-dt / tr);
    if (u <= 0.0f) return 0.0f;
    return std::pow(u, p) * std::exp(-dt / tf) * t_max_inv;
}

inline float t_max_value_p(float tr, float tf, float p)
{
    if (!(tr > 0.0f) || !(tf > 0.0f) || !(p > 0.0f)) return 1.0f;  // defensive
    const float denom = tr + p * tf;
    const float u_peak = p * tf / denom;            // (1 - exp(-dt_peak/τ_r))
    return std::pow(u_peak, p) * std::pow(tr / denom, tr / tf);
}

inline float chi2_eval_two_tau_p(const float *y, int n, float clk_ns,
                                  float t0, float tr, float tf, float p)
{
    const float t_max_inv = 1.0f / t_max_value_p(tr, tf, p);
    float s = 0.0f;
    for (int i = 0; i < n; ++i) {
        const float t = i * clk_ns;
        const float r = y[i] - two_tau_p_unit_pt(t, t0, tr, tf, p, t_max_inv);
        s += r * r;
    }
    return s;
}

} // anon

WaveAnalyzer::PulseFitResult
WaveAnalyzer::FitPulseShape(const uint16_t *slice, int nslice,
                            int peak_idx_in_slice,
                            float ped, float ped_rms,
                            float clk_ns,
                            float model_err_floor)
{
    PulseFitResult res{};
    res.ok = false;

    if (!slice || nslice < 8 || nslice > MAX_SAMPLES) return res;
    if (peak_idx_in_slice < 0 || peak_idx_in_slice >= nslice) return res;

    const float peak_amp = static_cast<float>(slice[peak_idx_in_slice]) - ped;
    if (peak_amp <= 0.0f) return res;
    res.peak_amp = peak_amp;

    // Pedsub + normalise to unit peak.
    float y[MAX_SAMPLES];
    const float inv_amp = 1.0f / peak_amp;
    for (int i = 0; i < nslice; ++i)
        y[i] = (static_cast<float>(slice[i]) - ped) * inv_amp;

    // σ on the normalised pulse: relative noise, floored at the model-
    // error scale so χ²/dof stays sane on high-amplitude pulses.
    const float sigma_noise = std::max(ped_rms, 1.0f) * inv_amp;
    const float sigma       = std::max(sigma_noise, model_err_floor);
    const float w_inv2      = 1.0f / (sigma * sigma);

    // Initial guesses — same as the previous Python fit.
    float t0 = peak_idx_in_slice * clk_ns - 2.0f * clk_ns;
    float tr = 1.0f * clk_ns;
    float tf = 5.0f * clk_ns;

    const float t0_lo = -2.0f * clk_ns;
    const float t0_hi = (nslice - 1) * clk_ns;
    const float tr_lo = 0.2f * clk_ns;
    const float tr_hi = 10.0f * clk_ns;
    const float tf_lo = 1.0f * clk_ns;
    const float tf_hi = 80.0f * clk_ns;

    constexpr int   MAX_ITER  = 100;       // ~scipy curve_fit default budget
    constexpr float TOL_PARAM = 1e-5f;     // relative param step in clk_ns units
    constexpr float LAMBDA0   = 1.0e-3f;
    constexpr float LAMBDA_UP = 10.0f;
    constexpr float LAMBDA_DN = 10.0f;
    constexpr float LAMBDA_MAX = 1.0e10f;

    float lambda = LAMBDA0;
    float chi2   = chi2_eval(y, nslice, clk_ns, t0, tr, tf);

    // Best params seen during the loop — we return these unconditionally
    // at the end, matching scipy.optimize.curve_fit semantics ("hand back
    // the best found, let the caller decide via χ²/dof whether the fit
    // is useful").  Without this, ~12% of fits where the initial guess
    // is already near-optimal silently report ok=false because no LM
    // step strictly improves χ².
    float chi2_best = chi2;
    float t0_best = t0, tr_best = tr, tf_best = tf;

    int iter = 0;
    for (; iter < MAX_ITER; ++iter) {
        // Residuals at current point.
        float r[MAX_SAMPLES];
        const float t_max_inv = 1.0f / t_max_value(tr, tf);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            r[i] = y[i] - two_tau_unit_pt(t, t0, tr, tf, t_max_inv);
        }

        // Forward-difference Jacobian columns (n_data × 3, column-major
        // with stride MAX_SAMPLES).
        float J[MAX_SAMPLES * 3];
        const float h_t0 = std::max(1e-3f * clk_ns, 1e-6f);
        const float h_tr = std::max(1e-3f * tr,     1e-6f);
        const float h_tf = std::max(1e-3f * tf,     1e-6f);

        const float tmi_t0 = t_max_inv;            // T_max independent of t0
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_unit_pt(t, t0,        tr, tf, tmi_t0);
            const float fp = two_tau_unit_pt(t, t0 + h_t0, tr, tf, tmi_t0);
            J[i + 0 * MAX_SAMPLES] = -(fp - f0) / h_t0;     // dr/dt0
        }
        const float tmi_trp = 1.0f / t_max_value(tr + h_tr, tf);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_unit_pt(t, t0, tr,        tf, t_max_inv);
            const float fp = two_tau_unit_pt(t, t0, tr + h_tr, tf, tmi_trp);
            J[i + 1 * MAX_SAMPLES] = -(fp - f0) / h_tr;     // dr/dτ_r
        }
        const float tmi_tfp = 1.0f / t_max_value(tr, tf + h_tf);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_unit_pt(t, t0, tr, tf,        t_max_inv);
            const float fp = two_tau_unit_pt(t, t0, tr, tf + h_tf, tmi_tfp);
            J[i + 2 * MAX_SAMPLES] = -(fp - f0) / h_tf;     // dr/dτ_f
        }

        // Normal equations: A = J^T J, g = J^T r  (uniform weight w_inv2
        // factors out — same scaling on A and g, solution unchanged).
        float A[9] = {0};
        float g[3] = {0};
        for (int i = 0; i < nslice; ++i) {
            const float j0 = J[i + 0 * MAX_SAMPLES];
            const float j1 = J[i + 1 * MAX_SAMPLES];
            const float j2 = J[i + 2 * MAX_SAMPLES];
            A[0] += j0 * j0; A[1] += j0 * j1; A[2] += j0 * j2;
                              A[4] += j1 * j1; A[5] += j1 * j2;
                                                A[8] += j2 * j2;
            g[0] += j0 * r[i]; g[1] += j1 * r[i]; g[2] += j2 * r[i];
        }
        A[3] = A[1]; A[6] = A[2]; A[7] = A[5];

        // Augment A with λ on the diagonal and solve for the step.
        float Aug[9] = {A[0]+lambda, A[1], A[2],
                        A[3], A[4]+lambda, A[5],
                        A[6], A[7], A[8]+lambda};
        float delta[3];
        if (!chol3_solve(Aug, g, delta)) {
            lambda *= LAMBDA_UP;
            if (lambda > LAMBDA_MAX) break;
            continue;
        }

        // Trial point, clamped to bounds.
        //
        // Sign: we set up J_ik = ∂r_i/∂p_k = -∂model_i/∂p_k and g = J^T r,
        // then solved (J^T J + λI) δ = g.  The Gauss-Newton step that
        // minimises ||r + J δ||² is δ_GN = -(J^T J + λI)^{-1} J^T r, i.e.
        // the NEGATIVE of what chol3_solve returned — so the parameter
        // update is `p - delta`, not `p + delta`.
        const float new_t0 = std::clamp(t0 - delta[0], t0_lo, t0_hi);
        const float new_tr = std::clamp(tr - delta[1], tr_lo, tr_hi);
        const float new_tf = std::clamp(tf - delta[2], tf_lo, tf_hi);
        const float new_chi2 = chi2_eval(y, nslice, clk_ns, new_t0, new_tr, new_tf);

        // Always remember the best point we've seen, even if the strict-
        // accept test below rejects it.
        if (std::isfinite(new_chi2) && new_chi2 < chi2_best) {
            chi2_best = new_chi2;
            t0_best = new_t0; tr_best = new_tr; tf_best = new_tf;
        }

        if (new_chi2 < chi2) {
            const float dpar = std::abs(delta[0]) / clk_ns
                             + std::abs(delta[1]) / clk_ns
                             + std::abs(delta[2]) / clk_ns;
            t0 = new_t0; tr = new_tr; tf = new_tf;
            chi2 = new_chi2;
            lambda = std::max(lambda / LAMBDA_DN, 1.0e-12f);
            if (dpar < TOL_PARAM) { ++iter; break; }
        } else {
            lambda *= LAMBDA_UP;
            if (lambda > LAMBDA_MAX) break;
        }
    }

    // Return the best params seen.  ok=false only if every candidate
    // produced NaN/Inf — initial guess can't get worse than that.
    if (!std::isfinite(t0_best) || !std::isfinite(tr_best) || !std::isfinite(tf_best))
        return res;

    const int dof = std::max(1, nslice - 3);
    res.t0_ns        = t0_best;
    res.tau_r_ns     = tr_best;
    res.tau_f_ns     = tf_best;
    res.chi2_per_dof = (chi2_best * w_inv2) / static_cast<float>(dof);
    res.n_iter       = iter;
    res.ok           = true;
    return res;
}

WaveAnalyzer::PulseFitTwoTauPResult
WaveAnalyzer::FitPulseShapeTwoTauP(const uint16_t *slice, int nslice,
                                    int peak_idx_in_slice,
                                    float ped, float ped_rms,
                                    float clk_ns,
                                    float model_err_floor)
{
    PulseFitTwoTauPResult res{};
    res.ok = false;

    if (!slice || nslice < 8 || nslice > MAX_SAMPLES) return res;
    if (peak_idx_in_slice < 0 || peak_idx_in_slice >= nslice) return res;

    const float peak_amp = static_cast<float>(slice[peak_idx_in_slice]) - ped;
    if (peak_amp <= 0.0f) return res;
    res.peak_amp = peak_amp;

    float y[MAX_SAMPLES];
    const float inv_amp = 1.0f / peak_amp;
    for (int i = 0; i < nslice; ++i)
        y[i] = (static_cast<float>(slice[i]) - ped) * inv_amp;

    const float sigma_noise = std::max(ped_rms, 1.0f) * inv_amp;
    const float sigma       = std::max(sigma_noise, model_err_floor);
    const float w_inv2      = 1.0f / (sigma * sigma);

    // Initial guesses — match FitPulseShape's, plus p starting at the
    // sigmoidal value 2.0 (mid-range for the expected [1, 3] span; LM
    // can drift up or down within bounds).
    float t0 = peak_idx_in_slice * clk_ns - 2.0f * clk_ns;
    float tr = 1.0f * clk_ns;
    float tf = 5.0f * clk_ns;
    float p  = 2.0f;

    const float t0_lo = -2.0f * clk_ns;
    const float t0_hi = (nslice - 1) * clk_ns;
    const float tr_lo = 0.2f * clk_ns;
    const float tr_hi = 10.0f * clk_ns;
    const float tf_lo = 1.0f * clk_ns;
    const float tf_hi = 80.0f * clk_ns;
    const float p_lo  = 0.3f;     // < 1 = sharp/under-shoot rise (allowed but unusual)
    const float p_hi  = 10.0f;    // > 10 = pathologically smooth

    constexpr int   MAX_ITER  = 100;
    constexpr float TOL_PARAM = 1e-5f;
    constexpr float LAMBDA0   = 1.0e-3f;
    constexpr float LAMBDA_UP = 10.0f;
    constexpr float LAMBDA_DN = 10.0f;
    constexpr float LAMBDA_MAX = 1.0e10f;

    float lambda = LAMBDA0;
    float chi2   = chi2_eval_two_tau_p(y, nslice, clk_ns, t0, tr, tf, p);

    float chi2_best = chi2;
    float t0_best = t0, tr_best = tr, tf_best = tf, p_best = p;

    int iter = 0;
    for (; iter < MAX_ITER; ++iter) {
        // Residuals at current point.
        float r[MAX_SAMPLES];
        const float t_max_inv = 1.0f / t_max_value_p(tr, tf, p);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            r[i] = y[i] - two_tau_p_unit_pt(t, t0, tr, tf, p, t_max_inv);
        }

        // Forward-difference Jacobian (4 columns, n_data × 4, column-
        // major with stride MAX_SAMPLES).
        float J[MAX_SAMPLES * 4];
        const float h_t0 = std::max(1e-3f * clk_ns, 1e-6f);
        const float h_tr = std::max(1e-3f * tr,     1e-6f);
        const float h_tf = std::max(1e-3f * tf,     1e-6f);
        const float h_p  = std::max(1e-3f * p,      1e-6f);

        const float tmi_t0 = t_max_inv;            // T_max independent of t0
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_p_unit_pt(t, t0,        tr, tf, p, tmi_t0);
            const float fp = two_tau_p_unit_pt(t, t0 + h_t0, tr, tf, p, tmi_t0);
            J[i + 0 * MAX_SAMPLES] = -(fp - f0) / h_t0;     // dr/dt0
        }
        const float tmi_trp = 1.0f / t_max_value_p(tr + h_tr, tf, p);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_p_unit_pt(t, t0, tr,        tf, p, t_max_inv);
            const float fp = two_tau_p_unit_pt(t, t0, tr + h_tr, tf, p, tmi_trp);
            J[i + 1 * MAX_SAMPLES] = -(fp - f0) / h_tr;     // dr/dτ_r
        }
        const float tmi_tfp = 1.0f / t_max_value_p(tr, tf + h_tf, p);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_p_unit_pt(t, t0, tr, tf,        p, t_max_inv);
            const float fp = two_tau_p_unit_pt(t, t0, tr, tf + h_tf, p, tmi_tfp);
            J[i + 2 * MAX_SAMPLES] = -(fp - f0) / h_tf;     // dr/dτ_f
        }
        const float tmi_pp = 1.0f / t_max_value_p(tr, tf, p + h_p);
        for (int i = 0; i < nslice; ++i) {
            const float t = i * clk_ns;
            const float f0 = two_tau_p_unit_pt(t, t0, tr, tf, p,       t_max_inv);
            const float fp = two_tau_p_unit_pt(t, t0, tr, tf, p + h_p, tmi_pp);
            J[i + 3 * MAX_SAMPLES] = -(fp - f0) / h_p;      // dr/dp
        }

        // Normal equations.  A is 4×4 row-major.
        float A[16] = {0};
        float g[4]  = {0};
        for (int i = 0; i < nslice; ++i) {
            const float j0 = J[i + 0 * MAX_SAMPLES];
            const float j1 = J[i + 1 * MAX_SAMPLES];
            const float j2 = J[i + 2 * MAX_SAMPLES];
            const float j3 = J[i + 3 * MAX_SAMPLES];
            A[ 0] += j0*j0; A[ 1] += j0*j1; A[ 2] += j0*j2; A[ 3] += j0*j3;
                            A[ 5] += j1*j1; A[ 6] += j1*j2; A[ 7] += j1*j3;
                                            A[10] += j2*j2; A[11] += j2*j3;
                                                            A[15] += j3*j3;
            g[0] += j0 * r[i]; g[1] += j1 * r[i]; g[2] += j2 * r[i]; g[3] += j3 * r[i];
        }
        // Mirror upper to lower (keep both halves filled for cholesky/chol_solve).
        A[ 4] = A[ 1];  A[ 8] = A[ 2];  A[12] = A[ 3];
        A[ 9] = A[ 6];  A[13] = A[ 7];
        A[14] = A[11];

        // (A + λI) δ = g  via the generic K=4 Cholesky.
        float Aug[16];
        for (int i = 0; i < 16; ++i) Aug[i] = A[i];
        Aug[ 0] += lambda; Aug[ 5] += lambda; Aug[10] += lambda; Aug[15] += lambda;

        float L[16];
        const float min_pivot_sq = cholesky(Aug, L, 4);
        if (min_pivot_sq < 0.0f) {
            lambda *= LAMBDA_UP;
            if (lambda > LAMBDA_MAX) break;
            continue;
        }

        float delta[4];
        chol_solve(L, g, delta, 4);

        // Trial point — same sign convention as FitPulseShape.
        const float new_t0 = std::clamp(t0 - delta[0], t0_lo, t0_hi);
        const float new_tr = std::clamp(tr - delta[1], tr_lo, tr_hi);
        const float new_tf = std::clamp(tf - delta[2], tf_lo, tf_hi);
        const float new_p  = std::clamp(p  - delta[3], p_lo,  p_hi);
        const float new_chi2 = chi2_eval_two_tau_p(y, nslice, clk_ns,
                                                   new_t0, new_tr, new_tf, new_p);

        if (std::isfinite(new_chi2) && new_chi2 < chi2_best) {
            chi2_best = new_chi2;
            t0_best = new_t0; tr_best = new_tr; tf_best = new_tf; p_best = new_p;
        }

        if (new_chi2 < chi2) {
            const float dpar = std::abs(delta[0]) / clk_ns
                             + std::abs(delta[1]) / clk_ns
                             + std::abs(delta[2]) / clk_ns
                             + std::abs(delta[3]) / std::max(p, 0.1f);
            t0 = new_t0; tr = new_tr; tf = new_tf; p = new_p;
            chi2 = new_chi2;
            lambda = std::max(lambda / LAMBDA_DN, 1.0e-12f);
            if (dpar < TOL_PARAM) { ++iter; break; }
        } else {
            lambda *= LAMBDA_UP;
            if (lambda > LAMBDA_MAX) break;
        }
    }

    if (!std::isfinite(t0_best) || !std::isfinite(tr_best) ||
        !std::isfinite(tf_best) || !std::isfinite(p_best))
        return res;

    const int dof = std::max(1, nslice - 4);
    res.t0_ns        = t0_best;
    res.tau_r_ns     = tr_best;
    res.tau_f_ns     = tf_best;
    res.p            = p_best;
    res.chi2_per_dof = (chi2_best * w_inv2) / static_cast<float>(dof);
    res.n_iter       = iter;
    res.ok           = true;
    return res;
}

