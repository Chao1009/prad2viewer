//=============================================================================
// Fadc250FwAnalyzer.cpp — firmware-faithful FADC250 Mode 1/2/3 emulation.
//
// Algorithm references throughout this file point to:
//
//   FADC250 User's Manual / "FIRMWARE for FADC250 Ver2 ADC FPGA":
//     • §Pedestal Subtraction          (manual ~line 1262)
//     • §Mode 1 (Pulse Mode)           (~1321)
//     • §Mode 2 (Integral Mode)        (~1342)
//     • §Mode 3 TDC Algorithm Overview (~1361)
//     • §Requirements for TDC Algorithm (~2365)
//     • §TDC Algorithm for Mode 3      (~2371)  ← seven literal steps
//
// Two ambiguities in the manual are resolved here against the algorithm
// overview (and the firmware HDL block diagram, which contains a
// `Linear_Interpolation` / `Divide_18By12` block for Tfine):
//
//   [R1]  Manual step 1f literal: Vaverage = (Vpeak − Vmin) / 2 (delta).
//         Algorithm overview: Va is the *value* between Vmin and Vpeak.
//         For the bracketing search to be coherent, Va must be in the same
//         coordinate system as the samples → use absolute mid:
//             Va = Vmin + (Vpeak − Vmin) / 2.
//
//   [R2]  Manual step 2b literal: "Read until Vram > Vmin. Vba = Vram."
//         Taken literally this would put Vba at the first sample above noise
//         (= pulse start), which destroys the sub-sample fine time the
//         algorithm is supposed to produce.  Algorithm overview: "fine time
//         is the interpolating value of mid value away from next sample"
//         → bracket Va, not Vmin:
//             walk k from i_start until samples[k] ≥ Va
//             Vba = samples[k−1], Vaa = samples[k], coarse = k − 1
//
// Both resolutions are also documented in docs/clas_fadc/FADC250_algorithms.md.
//=============================================================================

#include "Fadc250FwAnalyzer.h"

#include <algorithm>
#include <cmath>

namespace fdec
{

void Fadc250FwAnalyzer::Analyze(const uint16_t *raw, int n, float PED,
                                DaqWaveResult &result) const
{
    result.clear();
    if (!raw || n <= 0 || n > MAX_SAMPLES) return;

    // ── Step 0 — Pedestal subtraction (MANUAL §Pedestal Subtraction, ~1262)
    //
    //   "programmable pedestal value in the Trigger Data Path. The result is
    //    not allowed to go below zero."   →   s'[i] = max(0, s[i] − PED)
    //
    // The TET comparison and all subsequent windowing operate on s'.
    float s[MAX_SAMPLES];
    for (int i = 0; i < n; ++i) {
        float v = static_cast<float>(raw[i]) - PED;
        s[i] = (v > 0.0f) ? v : 0.0f;
    }

    // MANUAL §Requirements (~2366): "There must be at least 5 samples
    // (background) before pulse. Four of these samples are used to determine
    // the pedestal (Vnoise) floor."  Below 5 samples we can't even compute
    // the noise floor — bail out cleanly.
    if (n < 5) return;

    // ── Step 1b — Vnoise from first 4 samples (MANUAL §TDC step 1b, ~2375).
    const float Vnoise = (s[0] + s[1] + s[2] + s[3]) * 0.25f;
    result.vnoise = Vnoise;

    // ── Step 1c — Vmin = Vnoise (MANUAL §TDC step 1c, ~2379).
    // Constant for the entire window — the manual's "first value greater
    // than Vnoise" reading in §Requirements describes where a pulse *starts*,
    // not the value used in the Va formula (cf. step 1c which is the actual
    // assignment the firmware performs).  See FADC250_algorithms.md §7.1.
    const float Vmin = Vnoise;

    const float TET = cfg.TET;
    const int   nsb = cfg.NSB;
    const int   nsa = cfg.NSA;
    const int   max_pulses = std::min(cfg.MAX_PULSES, MAX_PEAKS);
    const float clk_per_64 = cfg.CLK_NS / 64.0f;  // ns per fine-time LSB

    // Search starts at index 4 — manual §Requirements: "at least 5 samples
    // (background) before pulse."  Four for Vnoise + at least one to settle.
    int i = 4;
    int pulse_idx = 0;

    while (i < n - 1 && pulse_idx < max_pulses) {

        // ── Find pulse start: first sample above Vnoise.
        while (i < n && s[i] <= Vnoise) ++i;
        if (i >= n) break;
        const int i_start = i;

        // ── Step 1d — Walk to peak.  MANUAL: "Read until Vram < Vram_delay.
        // Vpeak = Vram if Vram is greater than TET."  We treat the test as
        // strict less-than: keep walking through equal samples (plateau);
        // i_peak = last sample before strict decrease.
        int   i_peak = i_start;
        float Vpeak  = s[i_start];
        ++i;
        while (i < n) {
            if (s[i] >= s[i - 1]) {
                i_peak = i;
                Vpeak  = s[i];
                ++i;
            } else {
                // s[i] < s[i-1] → previous sample was the peak.
                break;
            }
        }
        // (If the loop exited because i == n, i_peak / Vpeak already point
        //  at s[n-1] from the last iteration.  Quality flag below covers it.)

        // ── Threshold gate (MANUAL §TDC step 1d, "Vpeak > TET").
        // Reject sub-threshold bumps; advance past the trailing edge to
        // baseline before resuming the outer search.
        if (Vpeak <= TET) {
            while (i < n && s[i] > Vnoise) ++i;
            continue;
        }

        // ── Tcross — first leading-edge sample exceeding TET.
        // MANUAL §Mode 1 (~1321): "When an ADC sample has a value that is
        // greater than Programmable Trigger Energy Threshold (TET), the
        // number of samples before (NSB)... and after (NSA) Vp are sent..."
        // §Pulse Raw Word 1 data format reports this as the "first sample
        // number for pulse" field.
        int cross = i_start;
        while (cross <= i_peak && s[cross] <= TET) ++cross;
        if (cross > i_peak) cross = i_peak;     // pathological — clamp.

        // ── Step 1f [R1] — Va = absolute mid amplitude.
        const float Va = Vmin + (Vpeak - Vmin) * 0.5f;

        // ── Step 2 [R2] — Bracket Va on the leading edge.
        //   Walk k forward from i_start until s[k] ≥ Va.
        //   Vba = s[k−1]   (sample below Va)
        //   Vaa = s[k]     (sample at-or-above Va)
        //   coarse = k − 1 (4-ns clock index of Vba, manual's 10-bit field)
        int k = i_start;
        while (k <= i_peak && s[k] < Va) ++k;

        float Vba, Vaa;
        int   coarse;
        uint8_t quality = Q_GOOD;

        if (k <= i_start) {
            // Va sits at or below the very first pulse sample (rise faster
            // than one sample period — Va was crossed between s[i_start-1]
            // and s[i_start]).  Bracket the transition into the pulse.
            quality |= Q_VA_OUT_OF_RANGE;
            if (i_start - 1 >= 0) {
                Vba    = s[i_start - 1];
                Vaa    = s[i_start];
                coarse = i_start - 1;
            } else {
                Vba    = s[i_start];
                Vaa    = s[i_start];
                coarse = i_start;
            }
        } else if (k > i_peak) {
            // Va above the peak — should not happen with Va = (Vmin+Vp)/2,
            // but guard anyway.
            quality |= Q_VA_OUT_OF_RANGE;
            Vba    = s[i_peak];
            Vaa    = s[i_peak];
            coarse = i_peak;
        } else {
            Vba    = s[k - 1];
            Vaa    = s[k];
            coarse = k - 1;
        }

        // ── Tfine (MANUAL §Mode 3 overview ~1364: "fine time is the
        // interpolating value of mid value away from next sample.  fine
        // value is 6 bit").  fine ∈ [0, 63].
        //
        // Edge case: Va == Vaa exactly (mid lands on the sample) gives
        // f = 64.  The firmware's 6-bit field can't hold 64; the natural
        // resolution is to carry into the coarse step (fine = 0, coarse++)
        // so the reported time lands exactly on the sample.  This keeps
        // T_units = (k − 1) · 64 + 64 ≡ k · 64 + 0 — same total, cleaner
        // representation.
        int fine;
        const float denom = Vaa - Vba;
        if (denom <= 0.0f) {
            fine = 0;
        } else {
            float f = (Va - Vba) / denom * 64.0f;
            int   r = static_cast<int>(std::lround(f));
            if (r < 0) r = 0;
            while (r >= 64) { r -= 64; ++coarse; }
            fine = r;
        }

        const int   time_units = coarse * 64 + fine;
        const float time_ns    = time_units * clk_per_64;

        // ── Mode 1 / Mode 2 windowing — [cross − NSB, cross + NSA], clamped.
        // MANUAL §Mode 1 (~1321) and §Mode 2 (~1342).  Quality flags record
        // truncation so downstream consumers can reject pulses with biased
        // integrals.
        int wlo = cross - nsb;
        int whi = cross + nsa;
        if (wlo < 0)      { wlo = 0;       quality |= Q_NSB_TRUNCATED; }
        if (whi >= n)     { whi = n - 1;   quality |= Q_NSA_TRUNCATED; }
        if (i_peak >= n - 1) quality |= Q_PEAK_AT_BOUNDARY;

        // Mode 2 — Σ s' over the (clamped) window.
        float integral = 0.0f;
        for (int j = wlo; j <= whi; ++j) integral += s[j];

        // ── Emit pulse.
        DaqPeak &p = result.peaks[pulse_idx];
        p.pulse_id     = pulse_idx;
        p.vmin         = Vmin;
        p.vpeak        = Vpeak;
        p.va           = Va;
        p.coarse       = coarse;
        p.fine         = fine;
        p.time_units   = time_units;
        p.time_ns      = time_ns;
        p.cross_sample = cross;
        p.integral     = integral;
        p.window_lo    = wlo;
        p.window_hi    = whi;
        p.quality      = quality;
        ++pulse_idx;

        // ── Step 6 — End of pulse: descend below Vmin (= Vnoise).  Resume
        // the outer search at the first below-baseline sample.
        int j = i_peak + 1;
        while (j < n && s[j] > Vnoise) ++j;
        i = j;
    }

    result.npeaks = pulse_idx;
}

} // namespace fdec
