#pragma once
//=============================================================================
// Fadc250FwAnalyzer.h — firmware-faithful FADC250 Mode 1/2/3 emulation.
//
// Reproduces the on-board JLab 250 MHz Flash ADC pulse identification
// algorithm exactly as documented in the FADC250 User's Manual (Ed
// Jastrzembski, JLab) — Mode 3 (TDC) drives Modes 1 (Pulse Raw) and 2
// (Pulse Integral) by sharing the same threshold-crossing detector.
//
// Sister analyzer to WaveAnalyzer:
//   * WaveAnalyzer       — software-friendly local-maxima search with
//                          iterative-outlier-rejection pedestal estimate.
//                          Tunable, robust against pile-up, returns Peak[].
//   * Fadc250FwAnalyzer  — firmware bit-faithful: Vnoise from first 4 samples,
//                          mid-amplitude bracket interpolation for Tfine,
//                          NSB/NSA windowing.  Returns DaqPeak[] with
//                          coarse/fine/quality fields you can compare to
//                          firmware-reported pulse data.
//
// The two analyzers run side-by-side; the viewer toggles between them.
//
// PED handling: at PRad-II's data path the recorded waveforms have already
// passed firmware TET (channels below PED+TET are not written out), so the
// software emulator just receives the soft analyzer's pedestal mean as PED
// and proceeds.  No per-channel firmware register file in v1.
//
// Configuration source: evc::DaqConfig::Fadc250FwConfig in DaqConfig.h, populated from
// the optional "fadc250_firmware" block in daq_config.json.
//
// Manual citations are embedded inline in the .cpp — see the algorithm
// summary in docs/clas_fadc/FADC250_algorithms.md for the full discussion
// (including the two manual-typo resolutions: Va is the absolute mid value,
// and the leading-edge bracket search targets Va, not Vmin).
//=============================================================================

#include "Fadc250Data.h"
#include "DaqConfig.h"

namespace fdec
{

class Fadc250FwAnalyzer
{
public:
    explicit Fadc250FwAnalyzer(const evc::DaqConfig::Fadc250FwConfig &cfg = {})
        : cfg(cfg) {}

    // Analyze one channel's waveform.  Fills `result` in place.  No heap
    // allocation; one stack-resident buffer of MAX_SAMPLES floats holds the
    // pedestal-subtracted samples.
    //
    //   samples — raw 12-bit ADC values (PED has NOT been subtracted yet)
    //   n       — number of samples (PTW); must be ≤ MAX_SAMPLES
    //   PED     — per-channel pedestal value to subtract (pass the soft
    //             analyzer's WaveResult.ped.mean if no firmware PED register
    //             is available)
    //   result  — populated with vnoise + up to cfg.MAX_PULSES DaqPeak entries
    void Analyze(const uint16_t *samples, int n, float PED,
                 DaqWaveResult &result) const;

    evc::DaqConfig::Fadc250FwConfig cfg;
};

} // namespace fdec
