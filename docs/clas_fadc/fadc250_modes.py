"""
FADC250 firmware-faithful waveform analysis (Modes 1, 2, 3).

Reference: FADC250 User's Manual (Ed Jastrzembski, JLab).
Sampling clock: 250 MHz -> 4 ns per sample.
TDC fine-time resolution: 1/(250 MHz * 64) = 62.5 ps.

Algorithms implemented (per manual, ADC FPGA Functional Description):

  Pedestal subtraction (firmware-common pre-step on the trigger path):
      sample' = max(0, sample - PED)   per channel.
      The Mode-3 TDC, Mode-1 windowing, and Mode-2 sums all operate on
      these pedestal-subtracted samples.  TET comparisons are also
      against the pedestal-subtracted value.

  Mode 3 (TDC):
      1. Average of first 4 samples -> Vnoise (pedestal floor).
      2. Vmin = first sample > Vnoise (pulse start).
      3. Walk forward; Vp = peak (sample where the next sample drops:
         Vram < Vram_delay), and Vp must be > TET.
      4. Va = (Vp - Vmin) / 2 + Vmin            (mid value on leading edge)
      5. On the leading edge, find Vba (sample before crossing Va) and
         Vaa (sample after).
      6. coarse = index of Vba (in 4 ns clocks from window start).
         fine   = round( (Va - Vba) / (Vaa - Vba) * 64 ),  6-bit, 0..63.
         T_total = coarse * 64 + fine     (units of 62.5 ps).
      7. After end of pulse (sample falls back below Vmin), repeat
         search for next pulse, up to 4 pulses per window.

  Mode 1 (Pulse mode):
      For each found pulse, report the raw pedestal-subtracted samples
      in the window [Tcross - NSB, Tcross + NSA], where Tcross is the
      first sample exceeding TET on the leading edge.  Also reports
      Tcross itself as 'sample number from threshold' (the manual's
      "first sample number for pulse").

  Mode 2 (Integral mode):
      Same pulse identification + window as Mode 1, but the reported
      quantity is Sum_i samples[i] over the [Tcross-NSB, Tcross+NSA]
      window (pedestal subtraction already applied).  Tcross is also
      reported.

The manual states T1/T2 (the pulse times produced by Mode 3) are a
prerequisite for Modes 1 and 2; in firmware the same threshold-crossing
detector drives all three modes.  Here Mode 3 (TDC) is computed first
and its results are reused, exactly as the manual describes.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Sequence
import math


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

@dataclass
class FADC250Config:
    """FADC250 per-channel configuration registers (manual sec. ADC FPGA)."""
    PED: float = 0.0      # Pedestal subtract value (register 0x001E..0x002D)
    TET: float = 50.0     # Trigger Energy Threshold, 12-bit (register 0x000B..)
    NSB: int = 4          # Samples before threshold crossing  (>=2)
    NSA: int = 10         # Samples after threshold crossing
                          #   >=3 for modes 0/1, >=6 for mode 2
    MAX_PULSES: int = 4   # Manual: up to 4 identified pulses per channel

    # Sampling clock period.  At 250 MHz the firmware uses 4 ns/sample;
    # exposed here so the same code can be reused in simulation studies
    # at other rates without changing the algorithm.
    CLK_NS: float = 4.0

    def validate(self) -> None:
        if self.NSB < 2:
            raise ValueError("NSB must be >= 2 (manual minimum).")
        if self.NSA < 3:
            raise ValueError("NSA must be >= 3 (manual minimum, modes 0/1).")
        if self.MAX_PULSES < 1 or self.MAX_PULSES > 4:
            raise ValueError("MAX_PULSES must be in 1..4.")


# ----------------------------------------------------------------------
# Result containers
# ----------------------------------------------------------------------

@dataclass
class PulseTDC:
    """Mode 3 (TDC) result for a single identified pulse."""
    pulse_number: int          # 0..3
    Vmin: float                # ADC counts (pedestal-subtracted)
    Vpeak: float               # ADC counts (pedestal-subtracted)
    Va: float                  # mid value = (Vp - Vm)/2 + Vm
    coarse_clk: int            # coarse time, in 4-ns clocks (sample index of Vba)
    fine: int                  # fine time, 6-bit (0..63)
    T_units: int               # coarse*64 + fine, LSB = 62.5 ps
    T_ns: float                # T_units * 1/(250 MHz * 64)  in ns
    sample_cross: int          # first sample index where signal > TET (Tcross)
    quality: int               # 0 = good; nonzero = edge-case flag


@dataclass
class PulseMode1:
    """Mode 1 (Pulse Raw) result for a single pulse."""
    pulse_number: int
    sample_cross: int          # first sample over TET on this pulse
    samples: List[float]       # raw (pedestal-subtracted) samples in [cross-NSB, cross+NSA]
    window_start: int          # index of first sample in `samples`
    window_end: int            # index of last  sample in `samples` (inclusive)
    T_ns: float                # carried over from TDC for convenience


@dataclass
class PulseMode2:
    """Mode 2 (Integral) result for a single pulse."""
    pulse_number: int
    sample_cross: int
    integral: float            # Sum of samples over [cross-NSB, cross+NSA]
    n_samples: int             # how many samples actually went into the sum
    T_ns: float                # from TDC


@dataclass
class TriggerResult:
    """All-mode output for one trigger window."""
    pedestal_floor: float                    # Vnoise from first 4 samples
    pulses_tdc:    List[PulseTDC]    = field(default_factory=list)
    pulses_mode1:  List[PulseMode1]  = field(default_factory=list)
    pulses_mode2:  List[PulseMode2]  = field(default_factory=list)


# ----------------------------------------------------------------------
# Core
# ----------------------------------------------------------------------

class FADC250Analyzer:
    """Firmware-faithful analyzer.  One instance per channel."""

    # Quality factor bits (manual: 2-bit field on Pulse Integral / Pulse Time).
    Q_GOOD              = 0
    Q_PEAK_AT_BOUNDARY  = 1   # peak is at the last sample (pulse may extend past window)
    Q_NSB_TRUNCATED     = 2   # NSB ran into start of window
    Q_NSA_TRUNCATED     = 3   # NSA ran past end of window

    def __init__(self, config: FADC250Config):
        config.validate()
        self.cfg = config

    # ------------------------------------------------------------------
    # Step 0:  pedestal subtraction (trigger-path pre-step in firmware)
    # ------------------------------------------------------------------
    def _pedestal_subtract(self, raw: Sequence[float]) -> List[float]:
        """sample' = max(0, sample - PED).  Matches the firmware clamp."""
        ped = self.cfg.PED
        return [max(0.0, s - ped) for s in raw]

    # ------------------------------------------------------------------
    # Mode 3 (TDC) — runs first, drives modes 1 & 2
    # ------------------------------------------------------------------
    def _run_tdc(self, samples: Sequence[float]) -> (List[PulseTDC], float):
        """
        Implements the TDC Algorithm exactly as written in the manual
        ("TDC Algorithm for Mode 3", Data Processing section):

          1) Read four samples; Vnoise = average.
          2) Vmin = Vnoise   (manual step 1c).
          3) Walk forward; track the peak.  Vpeak = sample whose successor
             decreases (Vram < Vram_delay).  Vpeak must exceed TET.
          4) Va = (Vpeak - Vmin) / 2 + Vmin.
          5) Re-walk from pulse start: Vba is sample below Va, Vaa is the
             next sample (>= Va) on the leading edge.
             coarse = sample index of Vba (in 4-ns clocks from window start)
             fine   = round((Va - Vba)/(Vaa - Vba) * 64), clamped to [0,63].
             T_total = coarse * 64 + fine    (LSB = 62.5 ps).
          6) Walk past the trailing edge (descend below Vnoise) and search
             for the next pulse, up to MAX_PULSES.

        Returns (pulses, Vnoise).
        """
        N = len(samples)
        cfg = self.cfg
        pulses: List[PulseTDC] = []

        if N < 5:
            # Manual: "There must be at least 5 samples (background) before pulse."
            return pulses, 0.0

        # --- 1) Vnoise from first 4 samples ----------------------------------
        Vnoise = sum(samples[:4]) / 4.0

        # --- 2) Vmin = Vnoise (manual step 1c) -------------------------------
        # Note: this is *constant* across all pulses in the window per the
        # firmware's pseudocode.  The "first sample > Vnoise" line in the
        # requirements section just describes where the pulse begins.
        Vmin_global = Vnoise

        i = 4                       # next sample to inspect
        pulse_idx = 0

        while i < N - 1 and pulse_idx < cfg.MAX_PULSES:

            # --- find next sample above Vnoise (start of a candidate pulse) -
            while i < N and samples[i] <= Vnoise:
                i += 1
            if i >= N:
                break
            i_start = i

            # --- 3) Walk to peak: peak = sample whose successor decreases ---
            # Manual: "Read until Vram < Vram_delay.  Vpeak = Vram if Vram > TET."
            i_peak  = i_start
            Vpeak   = samples[i_start]
            i += 1
            while i < N:
                if samples[i] >= samples[i - 1]:
                    i_peak = i
                    Vpeak  = samples[i]
                    i += 1
                else:
                    # samples[i] < samples[i-1]  -> previous was the peak
                    break
            else:
                # Reached end of buffer without ever turning over.
                i_peak = N - 1
                Vpeak  = samples[N - 1]

            # Pulse only counts if peak > TET.
            if Vpeak <= cfg.TET:
                # Skip past this bump and keep searching.
                while i < N and samples[i] > Vnoise:
                    i += 1
                continue

            # First sample on the leading edge that exceeds TET (Tcross).
            sample_cross = i_start
            while sample_cross <= i_peak and samples[sample_cross] <= cfg.TET:
                sample_cross += 1
            if sample_cross > i_peak:
                sample_cross = i_peak

            # --- 4) Va = mid value -------------------------------------------
            Vmin = Vmin_global                              # per manual
            Va   = (Vpeak - Vmin) / 2.0 + Vmin

            # --- 5) Find Vba, Vaa on the leading edge ------------------------
            # Walk from i_start forward until samples[k] >= Va.
            k = i_start
            while k <= i_peak and samples[k] < Va:
                k += 1

            if k <= i_start or k > i_peak:
                # Va either below the very first pulse sample (very fast rise:
                # the leading edge skipped over Va between two samples) or
                # somehow above the peak.  Use the best pair we have.
                if k <= i_start:
                    # interpolate between the last background sample and i_start
                    if i_start - 1 >= 0:
                        Vba    = samples[i_start - 1]
                        Vaa    = samples[i_start]
                        coarse = i_start - 1
                    else:
                        Vba, Vaa = samples[i_start], samples[i_start]
                        coarse   = i_start
                else:
                    Vba, Vaa = samples[i_peak], samples[i_peak]
                    coarse   = i_peak
            else:
                Vba    = samples[k - 1]
                Vaa    = samples[k]
                coarse = k - 1                       # 4-ns clock index of Vba

            denom = Vaa - Vba
            if denom <= 0.0:
                fine = 0
            else:
                f = (Va - Vba) / denom * 64.0
                fine = int(round(f))
                if fine < 0:  fine = 0
                if fine > 63: fine = 63

            T_units = coarse * 64 + fine
            T_ns    = T_units * (cfg.CLK_NS / 64.0)   # 4 ns / 64 = 62.5 ps

            quality = self.Q_GOOD
            if i_peak >= N - 1:
                quality = self.Q_PEAK_AT_BOUNDARY

            pulses.append(PulseTDC(
                pulse_number = pulse_idx,
                Vmin         = Vmin,
                Vpeak        = Vpeak,
                Va           = Va,
                coarse_clk   = coarse,
                fine         = fine,
                T_units      = T_units,
                T_ns         = T_ns,
                sample_cross = sample_cross,
                quality      = quality,
            ))
            pulse_idx += 1

            # --- 6) Walk past trailing edge to look for the next pulse ------
            # Manual: "Read until Vram < Vmin.  End of first pulse."
            # Since Vmin == Vnoise, this is "wait until back at baseline".
            j = i_peak + 1
            while j < N and samples[j] > Vnoise:
                j += 1
            i = j      # next pulse search resumes here

        return pulses, Vnoise

    # ------------------------------------------------------------------
    # Mode 1 (Pulse mode)
    # ------------------------------------------------------------------
    def _run_mode1(self,
                   samples: Sequence[float],
                   tdc_pulses: List[PulseTDC]) -> List[PulseMode1]:
        N    = len(samples)
        nsb  = self.cfg.NSB
        nsa  = self.cfg.NSA
        out: List[PulseMode1] = []

        for p in tdc_pulses:
            cross = p.sample_cross
            ws = cross - nsb
            we = cross + nsa
            if ws < 0:  ws = 0          # NSB truncation at start of window
            if we >= N: we = N - 1      # NSA truncation at end of window

            out.append(PulseMode1(
                pulse_number = p.pulse_number,
                sample_cross = cross,
                samples      = list(samples[ws:we + 1]),
                window_start = ws,
                window_end   = we,
                T_ns         = p.T_ns,
            ))
        return out

    # ------------------------------------------------------------------
    # Mode 2 (Integral mode)
    # ------------------------------------------------------------------
    def _run_mode2(self,
                   samples: Sequence[float],
                   tdc_pulses: List[PulseTDC]) -> List[PulseMode2]:
        N    = len(samples)
        nsb  = self.cfg.NSB
        nsa  = self.cfg.NSA
        out: List[PulseMode2] = []

        for p in tdc_pulses:
            cross = p.sample_cross
            ws = max(0,     cross - nsb)
            we = min(N - 1, cross + nsa)
            window = samples[ws:we + 1]
            integral = float(sum(window))
            out.append(PulseMode2(
                pulse_number = p.pulse_number,
                sample_cross = cross,
                integral     = integral,
                n_samples    = len(window),
                T_ns         = p.T_ns,
            ))
        return out

    # ------------------------------------------------------------------
    # Public entry point: run all three modes in one shot
    # ------------------------------------------------------------------
    def analyze(self, raw_samples: Sequence[float]) -> TriggerResult:
        """
        Run TDC -> Mode 1 -> Mode 2 on a single-channel waveform.

        Parameters
        ----------
        raw_samples : sequence of float
            ADC samples across the trigger window.  In the firmware this
            is PTW samples wide (e.g. 100 samples = 400 ns at 250 MHz).

        Returns
        -------
        TriggerResult
            Aggregated TDC / Mode 1 / Mode 2 results, plus the noise floor.
        """
        s = self._pedestal_subtract(raw_samples)
        tdc, vnoise = self._run_tdc(s)
        m1  = self._run_mode1(s, tdc)
        m2  = self._run_mode2(s, tdc)
        return TriggerResult(
            pedestal_floor = vnoise,
            pulses_tdc     = tdc,
            pulses_mode1   = m1,
            pulses_mode2   = m2,
        )
