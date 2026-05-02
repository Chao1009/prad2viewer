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
// NNLS pile-up deconvolution
//=============================================================================
//
// Given pedsub samples b[0..n-1] and K peak times τ_1..τ_K (each placed so
// that the template peak sits at the WaveAnalyzer-reported peak time), we
// fit non-negative amplitudes a_k by minimising ||M a - b||² subject to
// a_k ≥ 0, where M[i,k] = T(t_i − τ_k_onset; τ_r, τ_f).
//
// Algorithm: Lawson-Hanson active-set NNLS.  Hand-coded for K ≤ MAX_PEAKS
// because (a) we already need stack-only scratch, (b) external linear
// algebra deps would be intrusive in the analyzer hot path, and (c) the
// problem is tiny — typical K is 2..3 and Cholesky on K×K is essentially
// free.

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

// Solve L L^T x = b in place using a precomputed Cholesky factor.
inline void chol_solve(const float *L, const float *b, float *x, int K)
{
    float y[MAX_PEAKS];
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

// Active-set NNLS for small K (≤ MAX_PEAKS).  Inputs:
//   AtA  K×K precomputed M^T M (row-major)
//   Atb  K-vector M^T b
// Output:
//   x    K-vector solution, x_k ≥ 0
// Returns:
//   smallest pivot² ever encountered during the inner Cholesky solves —
//   < 0 if the active submatrix went indefinite (caller should treat
//   that as singular).  The conditioning proxy max_diag(AtA) / min_pivot²
//   is left to the caller.
inline float nnls_smallK(const float *AtA, const float *Atb, int K, float *x)
{
    // KKT tolerance — scale to the problem size so it's invariant to the
    // numerical scale of b.
    float Atb_max = 0.0f;
    for (int j = 0; j < K; ++j)
        Atb_max = std::max(Atb_max, std::abs(Atb[j]));
    const float tol_kkt = 1.0e-6f * std::max(Atb_max, 1.0f);
    const float tol_x   = 1.0e-7f * std::max(Atb_max, 1.0f);

    bool inP[MAX_PEAKS];          // false = active set R, true = passive set P
    for (int j = 0; j < K; ++j) { x[j] = 0.0f; inP[j] = false; }

    float min_pivot_sq_seen = std::numeric_limits<float>::infinity();
    constexpr int MAX_OUTER = 3 * MAX_PEAKS + 5;

    for (int outer = 0; outer < MAX_OUTER; ++outer) {
        // w = M^T (b - M x) = Atb - AtA x — KKT residual / gradient
        float w[MAX_PEAKS];
        for (int j = 0; j < K; ++j) {
            float wj = Atb[j];
            for (int k = 0; k < K; ++k) wj -= AtA[j * K + k] * x[k];
            w[j] = wj;
        }
        // Pick max w over the active set; converged if all are ≤ tol.
        int   best = -1;
        float wmax = tol_kkt;
        for (int j = 0; j < K; ++j) {
            if (!inP[j] && w[j] > wmax) { wmax = w[j]; best = j; }
        }
        if (best < 0) return min_pivot_sq_seen;  // optimal
        inP[best] = true;

        // Inner loop: solve LS on current P, back off if any x_p went < 0.
        constexpr int MAX_INNER = 3 * MAX_PEAKS + 5;
        for (int inner = 0; inner < MAX_INNER; ++inner) {
            // Build sub-matrix indexed by current P.
            int  idx[MAX_PEAKS]; int Kp = 0;
            for (int j = 0; j < K; ++j) if (inP[j]) idx[Kp++] = j;
            if (Kp == 0) break;

            float subM[MAX_PEAKS * MAX_PEAKS];
            float subb[MAX_PEAKS];
            for (int p = 0; p < Kp; ++p) {
                subb[p] = Atb[idx[p]];
                for (int q = 0; q < Kp; ++q)
                    subM[p * Kp + q] = AtA[idx[p] * K + idx[q]];
            }

            float L[MAX_PEAKS * MAX_PEAKS];
            float pivot_sq = cholesky(subM, L, Kp);
            if (pivot_sq < 0.0f) {
                // Sub-problem went singular — bail out and return a flag.
                return -1.0f;
            }
            if (pivot_sq < min_pivot_sq_seen) min_pivot_sq_seen = pivot_sq;

            float s_sub[MAX_PEAKS];
            chol_solve(L, subb, s_sub, Kp);

            float s[MAX_PEAKS];
            for (int j = 0; j < K; ++j) s[j] = 0.0f;
            for (int p = 0; p < Kp; ++p) s[idx[p]] = s_sub[p];

            // Find the most-violating component — smallest s_k for k in P.
            float alpha = 1.0f;
            int blocker = -1;
            for (int p = 0; p < Kp; ++p) {
                int j = idx[p];
                if (s[j] <= tol_x) {
                    const float denom = x[j] - s[j];
                    if (denom > tol_x) {
                        const float a = x[j] / denom;
                        if (a < alpha) { alpha = a; blocker = j; }
                    }
                }
            }

            if (blocker < 0) {
                // All s_k > 0 — accept the LS solution.
                for (int j = 0; j < K; ++j) x[j] = s[j];
                break;  // exit inner loop, back to outer
            }

            // Step part-way and drop blockers from P.
            for (int j = 0; j < K; ++j) x[j] = x[j] + alpha * (s[j] - x[j]);
            for (int j = 0; j < K; ++j)
                if (inP[j] && x[j] <= tol_x) { x[j] = 0.0f; inP[j] = false; }
        }
    }
    return min_pivot_sq_seen;
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
// onset at t0 (ns).  Writes n values into `col`.  Slow path — two exp()
// per sample.  Used for Python-constructed templates that don't carry a
// precomputed grid; production code (which always has a grid) goes
// through template_column_grid below.
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

// Fast path: linear-interpolate the precomputed unit-amplitude grid
// (PulseTemplate::grid sampled at tmpl.grid_clk_ns) onto sample times
// t_i = i·clk_ns shifted by t0.  Skips per-sample exp() — about 10×
// faster than template_column on the deconv hot loop.
//
// Caller guarantees tmpl.grid_clk_ns > 0.  Samples whose shifted index
// falls outside the grid's right edge are set to 0 (template has long
// since decayed below noise) or 0 for negative dt (pre-onset).
inline void template_column_grid(float *col, int n, float clk_ns,
                                 float t0, const PulseTemplate &tmpl)
{
    const float inv_grid = 1.0f / tmpl.grid_clk_ns;
    const int   last     = PulseTemplate::GRID_N - 1;
    for (int i = 0; i < n; ++i) {
        const float dt = i * clk_ns - t0;
        if (dt <= 0.0f) { col[i] = 0.0f; continue; }
        const float pos = dt * inv_grid;
        const int   k   = static_cast<int>(pos);
        if (k >= last) { col[i] = 0.0f; continue; }
        const float frac = pos - static_cast<float>(k);
        col[i] = tmpl.grid[k] + frac * (tmpl.grid[k + 1] - tmpl.grid[k]);
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
    const float tr = tmpl.tau_r_ns;
    const float tf = tmpl.tau_f_ns;
    if (!(tr >= dcfg.tau_r_min_ns && tr <= dcfg.tau_r_max_ns) ||
        !(tf >= dcfg.tau_f_min_ns && tf <= dcfg.tau_f_max_ns) ||
        !(tr > 0.0f) || !(tf > 0.0f)) {
        dec_out.state = Q_DECONV_BAD_TEMPLATE;
        return;
    }

    const float clk_ns  = (cfg.clk_mhz > 0.0f) ? (1000.0f / cfg.clk_mhz) : 4.0f;
    const float ped     = wres.ped.mean;
    // sigma is the per-sample noise that drives chi²/dof.  We take it as
    // the pedestal RMS (floored at 1 ADC) — purely a diagnostic value;
    // no filtering decision is made on chi²/dof here, so the
    // model-error-floor trick the Python fitter applies isn't needed
    // (high-amp pulses will report large chi²/dof but still produce
    // valid deconv values).
    const float sigma   = std::max(wres.ped.rms, 1.0f);

    float t_off = 0.0f, t_max = 0.0f;
    template_peak(tr, tf, t_off, t_max);

    // Build pedsub vector b and design matrix M.  M is laid out as K
    // back-to-back fixed-stride columns (column k at offset k*MAX_SAMPLES),
    // so M[k*MAX_SAMPLES + i] is the value of column k (= the k-th peak's
    // unit-amplitude template) at sample i.  We only fill / read the first
    // nsamples of each column; the tail (nsamples..MAX_SAMPLES-1) is
    // untouched by the AtA / Atb / chi² loops below.  Worst case 1600
    // floats on the stack (K=8, MAX_SAMPLES=200).
    float b[MAX_SAMPLES];
    for (int i = 0; i < nsamples; ++i)
        b[i] = static_cast<float>(samples[i]) - ped;

    float M[MAX_SAMPLES * MAX_PEAKS];
    const bool use_grid = (tmpl.grid_clk_ns > 0.0f);
    for (int k = 0; k < K; ++k) {
        const float t_pk = wres.peaks[k].time;       // ns from sample 0
        const float t0   = t_pk - t_off;
        if (use_grid)
            template_column_grid(&M[k * MAX_SAMPLES], nsamples, clk_ns, t0, tmpl);
        else
            template_column(&M[k * MAX_SAMPLES], nsamples, clk_ns, t0, tr, tf);
    }

    // Normal equations.  AtA is K×K row-major, Atb is K-vector.
    float AtA[MAX_PEAKS * MAX_PEAKS];
    float Atb[MAX_PEAKS];
    float diag_max = 0.0f;
    for (int k = 0; k < K; ++k) {
        float s = 0.0f;
        for (int i = 0; i < nsamples; ++i)
            s += M[k * MAX_SAMPLES + i] * b[i];
        Atb[k] = s;
        for (int j = 0; j <= k; ++j) {
            float dot = 0.0f;
            for (int i = 0; i < nsamples; ++i)
                dot += M[k * MAX_SAMPLES + i] * M[j * MAX_SAMPLES + i];
            AtA[k * K + j] = dot;
            AtA[j * K + k] = dot;
        }
        if (AtA[k * K + k] > diag_max) diag_max = AtA[k * K + k];
    }

    // Conditioning gate — if any active-set Cholesky pivot drops below
    // diag_max / cond_number_max, the system is too ill-conditioned to
    // trust per-peak amplitudes.  We surface that as Q_DECONV_SINGULAR.
    float a[MAX_PEAKS];
    const float min_pivot_sq = nnls_smallK(AtA, Atb, K, a);
    const float cond_proxy = (min_pivot_sq > 0.0f && diag_max > 0.0f)
                           ? (diag_max / min_pivot_sq)
                           : std::numeric_limits<float>::infinity();
    dec_out.cond_number = cond_proxy;

    if (min_pivot_sq < 0.0f || cond_proxy > dcfg.cond_number_max ||
        !std::isfinite(cond_proxy)) {
        dec_out.state = Q_DECONV_SINGULAR;
        return;
    }

    // Per-peak height and integral-over-window from the converged a_k.
    for (int k = 0; k < K; ++k) {
        dec_out.amplitude[k] = a[k];
        dec_out.height[k]    = a[k] * t_max;
    }
    for (int k = 0; k < K; ++k) {
        const int i_pk = static_cast<int>(std::lround(wres.peaks[k].time / clk_ns));
        const int lo   = std::max(0,        i_pk - dcfg.pre_samples);
        const int hi   = std::min(nsamples, i_pk + dcfg.post_samples + 1);
        float sum = 0.0f;
        for (int i = lo; i < hi; ++i) sum += M[k * MAX_SAMPLES + i];
        dec_out.integral[k] = a[k] * sum;
    }

    // Global χ²/dof of the reconstructed waveform.
    float chi2 = 0.0f;
    for (int i = 0; i < nsamples; ++i) {
        float fit = 0.0f;
        for (int k = 0; k < K; ++k) fit += a[k] * M[k * MAX_SAMPLES + i];
        const float r = (b[i] - fit) / sigma;
        chi2 += r * r;
    }
    const int dof = std::max(1, nsamples - K);
    dec_out.chi2_per_dof = chi2 / static_cast<float>(dof);

    dec_out.n     = K;
    dec_out.state = tmpl.is_global ? Q_DECONV_FALLBACK_GLOBAL : Q_DECONV_APPLIED;
}
