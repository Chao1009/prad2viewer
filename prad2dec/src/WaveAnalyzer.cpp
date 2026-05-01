#include "WaveAnalyzer.h"

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
}
