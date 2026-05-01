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
};

class WaveAnalyzer
{
public:
    explicit WaveAnalyzer(const WaveConfig &cfg = {}) : cfg(cfg) {}

    // Analyze one channel. Fills result in-place, no heap allocation.
    void Analyze(const uint16_t *samples, int nsamples, WaveResult &result) const;

    WaveConfig cfg;

private:
    // smooth into pre-allocated buffer (stack array from caller)
    void smooth(const uint16_t *raw, int n, float *buf) const;

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
