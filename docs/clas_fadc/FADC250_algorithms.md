# FADC250 Waveform Analysis — Modes 1, 2, 3

A faithful software reproduction of the JLab 250 MHz Flash ADC firmware
algorithms, as documented in *FADC250 User's Manual / FIRMWARE for FADC250
Ver2 ADC FPGA* (Ed Jastrzembski, JLab) — `docs/clas_fadc/FADC250UsersManual.pdf`.

This document is the **authoritative algorithm spec** for the C++ implementation
in `prad2dec/src/Fadc250FwAnalyzer.cpp`.  The Python reference in
`fadc250_modes.py` is a sketch that pre-dated this spec and may lag behind.
Manual line numbers below cite the plain-text extraction
`docs/clas_fadc/_manual.txt`.

## 1. Hardware context (manual §Mode 3 TDC overview ~1361, §Pulse Time data format ~2495)

| Parameter | Value |
|---|---|
| Sampling rate | 250 MSPS |
| Sample period | 4 ns |
| Trigger window (PTW) | up to 8 µs (manual register 0x0007) |
| Max pulses per channel per window | 1..4 (CONFIG 1 register 0x0003 bits 6-5) |
| Coarse time | 10-bit field (manual §Pulse Time, ~2502) |
| Fine time | 6-bit field (~2504) |
| TDC LSB | 1 / (250 MHz × 64) = **62.5 ps** |
| Per-channel pedestal register | manual register 0x001E..0x002D |
| Per-channel TET | 12-bit (registers 0x000B..0x001A) |

## 2. Configuration parameters (manual §registers ~777-815)

| Symbol | Manual register | Min | Notes |
|---|---|---|---|
| `PED` | 0x001E..0x002D | — | per-channel pedestal |
| `TET` | 0x000B..0x001A | — | trigger energy threshold, 12-bit, per channel |
| `NSB` | 0x0009 | 2 | samples before threshold crossing |
| `NSA` | 0x000A | 3 (modes 0/1), 6 (mode 2) | samples after threshold crossing |
| `MAX_PULSES` | CONFIG 1, bits 6-5 | 1 | up to 4 |

## 3. Pre-step — Pedestal subtraction (manual §Pedestal Subtraction ~1262)

Manual: *"programmable pedestal value in the Trigger Data Path. The result is
not allowed to go below zero. Each of the ADC has a separate pedestal value."*

```
s'[i] = max(0, s[i] − PED)
```

The TET comparison, Mode-3 TDC, and Mode-1/2 windowing all operate on these
pedestal-subtracted samples.

## 4. Mode 3 — TDC (must run first)

Mode 3 produces the timing reference (`Tcross`, the leading-edge threshold
crossing) and pulse identification that Modes 1 and 2 reuse.

### 4.1 Algorithm overview (manual ~1361)

> *"The TDC algorithm calculates time of the mid value (Va) of a pulse relative
> to the beginning of the look back window. Va is the value between the
> smallest and the peak value (Vp) of the pulse. The smallest value (Vm) is
> the beginning of the pulse. ... the fine time is the interpolating value
> of mid value away from next sample. The coarse value is 10 bits and the
> fine value is 6 bit."*

### 4.2 The seven steps (manual ~2371)

```
1) Search for Vaverage:
   a) Latch starting PTW_RAM ADR
   b) Read four samples. Vnoise = Average of 4 samples. sample_count += 4.
   c) Vmin = Vnoise. sample_count += 1.
   d) Read until Vram < Vram_delay. Vpeak = Vram if Vram > TET. sample_count++.
   e) Store PTW_RAM ADR for Vpeak.
   f) Vaverage = (Vpeak − Vmin) / 2          ← see [R1]

2) Search for sample before (Vba) and sample after (Vaa) Vaverage:
   a) Restore starting PTW_RAM ADR. Increment Pulse Timer per address advance.
   b) Read until Vram > Vmin. Vba = Vram     ← see [R2]
   c) Read one more for Vaa.
   d) Calculate Tfine.

3) Write Pulse Number and Pulse Timer to Processing RAM.
4) Increment Pulse Number.
5) Restore PTW_RAM ADR for Vpeak. Load sample count to Pulse Timer.
6) Read until Vram < Vmin. End of first pulse.
7) Go to step 1.
```

### 4.3 Manual ambiguities — resolved

The manual's literal step 1f and step 2b contradict the algorithm overview
(§4.1) and the firmware HDL block diagram (manual ~1704: `Linear_Interpolation`
→ `Divide_18By12` → `TDCSM`).  We resolve in favor of the overview — the
firmware HDL would not contain a divider block if step 2b literally just
identified the pulse-start sample.

| Tag | Manual literal | Resolved | Rationale |
|---|---|---|---|
| **[R1]** | `Va = (Vpeak − Vmin) / 2` (delta) | `Va = Vmin + (Vpeak − Vmin) / 2` (absolute mid) | Overview text + the `Vp / Va / Vm` figure (manual ~1372) place Va *between* Vm and Vp on the y-axis; the bracketing search in step 2 needs Va in the same coordinate system as the samples. The "/2" formula in step 1f reads as the *delta* between Vmin and Va, which still has to be added back to Vmin. |
| **[R2]** | `Read until Vram > Vmin. Vba = Vram` (Vba = first sample above noise = pulse start) | Walk `k` from `i_start` until `samples[k] ≥ Va`; `Vba = samples[k−1]`, `Vaa = samples[k]`, `coarse = k − 1` | Overview: *"fine time is the interpolating value of mid value away from next sample"* — the bracket is around **Va**, not Vmin. With literal "Vmin", every pulse would yield `Vba ≈ Vnoise` and `fine` would have no meaningful sub-sample resolution, defeating the entire 62.5 ps LSB. The HDL `Linear_Interpolation` block requires a non-trivial bracket pair. |

These two resolutions are cited inline in `prad2dec/src/Fadc250FwAnalyzer.cpp`
as `[R1]` / `[R2]` markers.

### 4.4 Implementation summary (after resolutions)

```
Vnoise = mean(s'[0..3])               # manual step 1b
Vmin   = Vnoise                       # manual step 1c (constant for the window)

i = 4
while pulses < MAX_PULSES:
    while s'[i] ≤ Vnoise: i++         # find pulse start
    i_start = i
    walk forward while s'[i] ≥ s'[i-1] (track i_peak / Vpeak)
    if Vpeak ≤ TET: skip past bump, continue        # manual step 1d (Vpeak > TET)

    cross = first leading-edge sample > TET
    Va    = Vmin + (Vpeak − Vmin) / 2               # [R1]

    walk k from i_start until s'[k] ≥ Va            # [R2]
    Vba    = s'[k − 1]
    Vaa    = s'[k]
    coarse = k − 1
    fine   = round((Va − Vba) / (Vaa − Vba) × 64)   # 0..63
    if fine == 64: fine = 0; coarse++               # carry into next clock
    T_units = coarse · 64 + fine
    T_ns    = T_units · CLK_NS / 64                 # 62.5 ps LSB at 250 MHz

    descend below Vnoise to find next pulse         # manual step 6
```

### 4.5 What Mode 3 emits per pulse

- `pulse_id` (0..MAX_PULSES − 1)
- `vmin`, `vpeak`, `va` — pedestal-subtracted ADC counts
- `coarse` (4-ns clocks; 10-bit field), `fine` (6-bit, 0..63)
- `time_units` (62.5 ps LSB), `time_ns` (convenience)
- `cross_sample` — leading-edge threshold-crossing index (manual's Pulse Raw
  Word 1 "first sample number for pulse" field, ~2468)
- `quality` — bitmask, see §7

## 5. Mode 1 — Pulse Raw (manual §Mode 1 ~1321)

For each pulse identified by Mode 3, return the raw pedestal-subtracted
samples in the window

```
[ cross_sample − NSB ,  cross_sample + NSA ]
```

clamped to `[0, N − 1]`.  Reported length is `NSB + NSA + 1` when the pulse
is fully contained.  Truncation cases set the corresponding quality bits
(see §7).

## 6. Mode 2 — Pulse Integral (manual §Mode 2 ~1342)

Same pulse identification + same window as Mode 1, but the reported quantity
is the scalar sum

```
integral = Σ s'[i]   for i in [cross − NSB, cross + NSA]
```

Pedestal subtraction is already baked into `s'`, so this is a pedestal-
subtracted charge integral.  By construction:

```
Mode2.integral  ==  sum(Mode1.samples)
```

## 7. Quality bitmask (`uint8_t`)

Manual reserves a 2-bit field (~2500); in software we have headroom and use a
bitmask so flags can compose:

| Constant | Bit | Meaning |
|---|---|---|
| `Q_GOOD` | 0 | no flags set |
| `Q_PEAK_AT_BOUNDARY` | 1 << 0 | `i_peak == N − 1` (pulse may extend past window) |
| `Q_NSB_TRUNCATED` | 1 << 1 | `cross − NSB < 0` |
| `Q_NSA_TRUNCATED` | 1 << 2 | `cross + NSA ≥ N` |
| `Q_VA_OUT_OF_RANGE` | 1 << 3 | Va not bracketed by leading-edge samples (very fast rise) |

## 8. Mode-execution order (manual ~1148, FIRMWARE block diagram)

```
                        ┌─────────────────────┐
   raw samples  ───────►│  pedestal subtract  │  s' = max(0, s − PED)
                        └─────────┬───────────┘
                                  │
                                  ▼
                        ┌─────────────────────┐
                        │  Mode 3 (TDC)       │  finds pulses, sets cross_sample
                        │  Vnoise → Vmin →    │  + coarse/fine via [R1]/[R2]
                        │  Vpeak → Va →       │
                        │  Vba/Vaa → fine     │
                        └─────────┬───────────┘
                                  │  pulses + cross_sample
                  ┌───────────────┴───────────────┐
                  ▼                               ▼
        ┌──────────────────┐           ┌──────────────────┐
        │  Mode 1: window  │           │  Mode 2: integral│
        │  raw samples in  │           │  Σ samples in    │
        │  [cross-NSB,     │           │  [cross-NSB,     │
        │   cross+NSA]     │           │   cross+NSA]     │
        └──────────────────┘           └──────────────────┘
```

Mode 3 always runs first because its `cross_sample` defines the integration
window for Modes 1 and 2.  The hardware does the same: a single threshold-
crossing detector drives all three processing options.

## 9. Subtle points

### 9.1 `Vmin = Vnoise`, not "first sample above Vnoise"

Manual §Requirements (~2366):
> *"There must be at least 5 samples (background) before pulse. ... The minimum
> value of the pulse is the first value that is greater than Vnoise."*

But the actual algorithm pseudocode (manual step 1c, ~2379):
> *"Vmin = Vnoise."*

The first quote describes *where a pulse begins*; the second is the assignment
the firmware actually performs.  These are **not** equivalent: setting `Vmin`
to the first rising-edge sample collapses Va onto Vpeak for fast pulses and
forces `fine = 0`, destroying sub-sample resolution.  Use `Vmin = Vnoise`.

### 9.2 End-of-pulse: walk to baseline, not to "old" Vmin

Manual step 6: *"Read until Vram < Vmin. End of first pulse."*  With
`Vmin = Vnoise` this is identical to descending below the baseline.  The
"first sample above noise" reading would re-trigger on the same pulse's tail.

### 9.3 Boundary truncation

If `cross + NSA ≥ N`, Mode 1 returns a shorter array and Mode 2 sums fewer
samples.  The manual's 2-bit quality field covers this case in firmware; in
software we set `Q_NSA_TRUNCATED` (and similarly `Q_NSB_TRUNCATED` /
`Q_PEAK_AT_BOUNDARY`).

### 9.4 Time offset between truth and reported `T_ns`

`Va` is the half-amplitude point on the leading edge, so the reported time is
naturally a fixed delay after the true pulse arrival time `t0`.  For a
`(1 − e^(−t/τ_rise))` rise, `Va` corresponds to roughly `t0 + τ_rise · ln 2`.
The offset is deterministic and cancels in any time-of-flight or coincidence
measurement.

### 9.5 Fine-time carry

When `Va == Vaa` exactly, `(Va − Vba) / (Vaa − Vba) = 1.0` and the literal
formula gives `fine = 64`, which doesn't fit in the firmware's 6-bit field.
We carry into coarse: `fine = 0`, `coarse++`.  This keeps `T_units` exact and
avoids the saturation case `fine = 63` (which would be off by one LSB).

## 10. Data-format mapping (manual §Data Formats ~933, ~2495)

The table below maps the C++ `DaqPeak` fields to the manual's 32-bit VME
readout words.  Useful when comparing software output to a real DAQ stream.

| Manual data type | Bit field | C++ field |
|---|---|---|
| Pulse Time (8) | (15..6) coarse | `coarse` |
| Pulse Time (8) | (5..0) fine | `fine` |
| Pulse Time (8) | (20..19) quality | `quality` (2 LSBs only — software has more flags) |
| Pulse Integral (7) | (18..0) | `integral` |
| Pulse Vmin/Vpeak (10) | (20..12) Vmin | `vmin` (firmware truncates to 9 bits; software keeps full float) |
| Pulse Vmin/Vpeak (10) | (11..0) Vpeak | `vpeak` |
| Pulse Raw (6) | (9..0) first-sample-number | `cross_sample` |
| Window Raw (4) | (11..0) PTW | `n_samples` |

## 11. References

- *FADC250 User's Manual*, Jefferson Lab Data Acquisition Group (R. Jones ed.;
  F.J. Barbosa, E. Jastrzembski, H. Dong, J. Wilson, C. Cuevas, D.J. Abbott
  authors).  PDF: `docs/clas_fadc/FADC250UsersManual.pdf`,
  text: `docs/clas_fadc/_manual.txt`.
- C++ implementation: `prad2dec/src/Fadc250FwAnalyzer.cpp` (manual citations
  inline as `MANUAL §...`).
- Hand-traced regression tests: `prad2dec/test/test_fadc_fw.cpp`
  (CTest target `prad2dec.fadc_fw`).
- Python sketch (pre-spec, may lag): `docs/clas_fadc/fadc250_modes.py`.
