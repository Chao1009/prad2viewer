# Software Waveform Analysis in `prad2dec`

`prad2dec` ships **two** offline analyzers that run on raw FADC250 samples
(`uint16_t[nsamples]`, 4 ns/sample at 250 MHz):

| Analyzer | Class | Purpose |
|---|---|---|
| Waveform | `fdec::WaveAnalyzer` (`WaveAnalyzer.{h,cpp}`) | Robust local-maxima peak finding for HyCal energy / time use. Tolerates noisy pedestals, finds multiple peaks per channel. |
| Firmware emulator | `fdec::Fadc250FwAnalyzer` (`Fadc250FwAnalyzer.{h,cpp}`) | Bit-faithful emulation of the JLab FADC250 firmware Mode 1/2/3 (Hall-D V3 + NSAT/NPED/MAXPED extensions). Used to compare offline reconstruction against on-board firmware output. |

Both are stack-allocated, zero-heap on the hot path, and run side-by-side
when `prad2ana_replay_rawdata` is invoked with the `-p` flag (see
[`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md)).

The remainder of this note walks through both algorithms on a single
example pulse, with parameter values matching the real run config in
[`database/daq_config.json`](../../../database/daq_config.json).

## Example waveform

100 samples × 4 ns. A short, bright pulse on top of a quiet ~146 ADC
baseline, followed by a long scintillation tail.

![overview](plots/fig1_overview.png)

| feature | value |
|---|---|
| length | 100 samples = 400 ns |
| baseline | ≈ 146 ADC (samples 0..29) |
| pulse onset | sample 30 (t = 120 ns) |
| peak | sample 32, 1393 ADC |
| rise time | 8 ns (cross → peak) |
| tail | slow exponential decay, still ~20 ADC above baseline at sample 99 |

This is representative of a HyCal PbWO₄ signal: a fast leading edge
(~10 ns) followed by a long PMT/scintillator tail.

## Waveform analyzer — `WaveAnalyzer`

Used for HyCal calibration / monitoring where we want a robust peak
height and a generous integral that follows the actual pulse shape rather
than a fixed firmware window.

### Parameters

All settings live in `fdec::WaveConfig` (see `WaveAnalyzer.h`); the values
shown are the defaults.

| field | default | unit | role |
|---|---:|---|---|
| `smooth_order`   |   2 | order | Triangular smoothing kernel order. `1` disables smoothing; `N` gives a `2N-1` tap kernel (effective half-width `N-1`). Larger values smear the rising edge but suppress per-sample noise. |
| `threshold`      | 5.0 | × pedestal RMS | Peak-acceptance threshold above the local baseline. |
| `min_threshold`  | 3.0 | ADC counts | Hard floor on the acceptance threshold — protects quiet channels where `5·rms` would underflow. |
| `min_peak_ratio` | 0.3 | fraction | When two peaks overlap, the secondary must be ≥ this fraction of the primary's height to survive. |
| `int_tail_ratio` | 0.1 | fraction | Integration stops when the pedsub waveform drops below `r × peak height`. Smaller values capture more of the tail. |
| `tail_break_n`   |   2 | samples | Integration only terminates after this many *consecutive* sub-threshold samples — a single noise dip in the tail no longer truncates the integral early. |
| `peak_pileup_gap`|   2 | samples | Two peaks whose integration windows are within this many samples of each other get the `Q_PEAK_PILED` flag set on both. |
| `ped_nsamples`   |  30 | samples | Window used for the pedestal estimate (start of the buffer). |
| `ped_flatness`   | 1.0 | ADC counts | Floor on the outlier-rejection band: samples are kept iff `|s − μ| < max(rms, ped_flatness)`. Prevents over-tight clipping on already-quiet baselines. |
| `ped_max_iter`   |   3 | iterations | Maximum outlier-rejection passes. Stops early if the mask doesn't change or fewer than 5 samples remain. |
| `overflow`       | 4095 | ADC counts | 12-bit overflow value — peaks at this height are tagged. |
| `clk_mhz`        | 250.0 | MHz | Sample rate for the time conversion `t_ns = pos · 1000 / clk_mhz`. |

### Pipeline

**1. Triangular smoothing.** `smooth_order = 2` (default) →
`buf[i] = (raw[i−1]·w + raw[i] + raw[i+1]·w) / (1 + 2w)` with
`w = 1 − 1/(smooth_order + 1)`. `smooth_order = 1` disables smoothing
(identity kernel); `smooth_order = N` uses a `2N − 1` tap kernel
(effective half-width `N − 1`).

**2. Iterative pedestal (median/MAD bootstrap).** First `ped_nsamples = 30`
samples of the *smoothed* trace.

- Seed `mean` with the **median** and `rms` with **MAD × 1.4826** of the
  window (robust against ≤50 % contamination — a previous-event tail or
  early ringing biases a simple-mean seed badly, loosens the σ-clip band,
  and the iteration can converge on a contaminated baseline).
- Iterate up to `ped_max_iter = 3` times: drop samples deviating more
  than `max(rms, ped_flatness)` from the mean, then re-compute mean/rms
  on the survivors.
- Track `nused` (surviving sample count) and a `Q_PED_*` quality bitmask
  (see [§ Pedestal quality](#pedestal-quality) below).
- Compute a least-squares `slope` (ADC/sample) on the survivors — catches
  baseline drift that the σ-clip alone hides.

For our trace: `mean = 145.61`, `rms = 0.45`, `nused = 28`, `slope ≈ 0`,
`quality = Q_PED_GOOD` after convergence.

**2b. Adaptive window.** If the leading window looks suspicious
(`Q_PED_NOT_CONVERGED`, `Q_PED_TOO_FEW_SAMPLES`, `Q_PED_OVERFLOW`, or
`nused < ped_nsamples / 2`) and the waveform is long enough that the
last `ped_nsamples` samples don't overlap the leading window
(`nsamples ≥ 2 · ped_nsamples`), the analyzer also estimates the
pedestal on the trailing samples and uses whichever has the lower RMS
(with `nused` as the tiebreaker).  When the trailing estimate wins,
`Q_PED_TRAILING_WINDOW` is set.

**3. Threshold.** `thr = max(threshold·rms, min_threshold)` =
`max(5·0.45, 3.0) = 3.0`. The hard floor `min_threshold` keeps the
threshold sane on quiet channels.

**4. Local-maxima search.** Walk smoothed buffer; a peak is accepted iff:

- it is a local max (the `trend()` flat-tolerance scales with `ped_rms`,
  so noise wiggles on quiet channels don't fragment plateaus into
  spurious mini-peaks)
- its height above the **local baseline** (linear interpolation between
  the surrounding minima) exceeds `thr`
- its height above the **pedestal mean** exceeds `thr` (the `thr` gate
  already implies `≥ 5·rms` in the default config; no separate
  `≥ 3·rms` check needed)

**5. Integration.** Walk outward from the peak, summing pedsub values.
A sample is treated as below-tail when `s′ < tail_cut =
int_tail_ratio · ped_height` (default 10 % of peak height) or
`s′ < ped_rms`.  Integration only terminates after **`tail_break_n`
consecutive** below-tail samples — a single noise dip is held in
`pending` and either committed (on recovery) or discarded (when the
run reaches the threshold).  This adapts to the pulse shape (wide
pulses get wide windows, narrow pulses get narrow ones) without being
sensitive to one-sample wiggles in the slow tail.

`peak.left` / `peak.right` are **INCLUSIVE** integration bounds: both
samples are in `peak.integral`, and the JS viewer's shading
(`waveform.js:185`, `for (j = p.l; j <= p.r; ++j)`) renders them
faithfully.

**6. Raw-position correction.** The recorded `pos` is the raw-sample
maximum near the smoothed peak (not the smoothed peak itself), so the
reported height equals the actual ADC at the peak rather than a smoothed
under-estimate.

**7. Sub-sample peak time.** A 3-point quadratic vertex through
`raw[pos−1], raw[pos], raw[pos+1]` lifts the time resolution from the
4 ns sample grid to ≪ 1 ns for clean peaks:
`δ = (h_-1 − h_+1) / (2·(h_-1 − 2·h_0 + h_+1))`, clamped to
`δ ∈ [−1, +1]` and applied only when the parabola is concave-down
(real maximum).  Reported as `peak.time = (pos + δ) · 1e3 / clk_mhz`
(ns).

**8. Pile-up flagging.** After each peak is accepted, its integration
window is compared to every previously-accepted peak's; if either
peak's bounds come within `cfg.peak_pileup_gap` samples of the other's,
both peaks get `Q_PEAK_PILED` set in `peak.quality`.  Diagnostic for
downstream cuts on isolated vs piled-up pulses.

For our trace:

| field | value |
|---|---|
| `peak.pos` | sample 32 |
| `peak.time` | 126.6 ns (raw + sub-sample δ = −0.36) |
| `peak.height` | 1247 ADC (raw − pedestal) |
| `peak.left, peak.right` | 29, 48 (INCLUSIVE — 20 samples in `peak.integral`) |
| `peak.integral` | 9600 (ADC·sample, pedsub) |
| `peak.quality` | `Q_PEAK_GOOD` (no pile-up) |

![waveform](plots/fig3_soft_analysis.png)

### Pedestal quality

The analyzer reports four scalars per channel that together describe how
trustworthy the pedestal estimate is — written to the `events` tree as
`hycal.ped_{mean,rms,nused,quality,slope}` (see
[`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md)):

| Field | Type | Use |
|---|---|---|
| `ped_mean`    | `float` | Pedestal mean after rejection |
| `ped_rms`     | `float` | RMS after rejection |
| `ped_nused`   | `uint8` | # samples that survived (compare to `ped_nsamples = 30`) |
| `ped_slope`   | `float` | LSQ drift (ADC/sample) on the survivors — non-zero suggests baseline tilt or pulse-tail contamination |
| `ped_quality` | `uint8` | `Q_PED_*` bitmask, see below |

Quality flags (defined in `prad2dec/include/Fadc250Data.h`):

| Bit | Flag | Set when |
|---|---|---|
| `0`     | `Q_PED_GOOD`             | clean estimate, no flags |
| `1<<0`  | `Q_PED_NOT_CONVERGED`    | `ped_max_iter` exhausted, kept-mask still moving |
| `1<<1`  | `Q_PED_FLOOR_ACTIVE`     | `rms < ped_flatness` — `ped_flatness` was the active band (typical for very quiet channels; informational) |
| `1<<2`  | `Q_PED_TOO_FEW_SAMPLES`  | < 5 samples survived rejection (rejection aborted) |
| `1<<3`  | `Q_PED_PULSE_IN_WINDOW`  | `findPeaks` returned a peak with `pos` inside the pedestal window we used |
| `1<<4`  | `Q_PED_OVERFLOW`         | a raw window sample hit `cfg.overflow` (4095) |
| `1<<5`  | `Q_PED_TRAILING_WINDOW`  | adaptive logic chose the trailing window over the leading one (informational, not a problem flag) |

A clean-event filter is `ped_quality == 0`.  For analyses that care
about the pedestal stability rather than the peak heights, cutting on
`Q_PED_NOT_CONVERGED | Q_PED_TOO_FEW_SAMPLES | Q_PED_PULSE_IN_WINDOW`
removes the events where the iterative cut couldn't settle on a clean
baseline.

### Parameter sensitivity

Two of the parameters above visibly change the analyzer's output on this
trace:

![params](plots/fig4_soft_parameters.png)

**Left — pedestal `ped_flatness` × `ped_max_iter`.** All 30 baseline
samples enter pass 1; samples deviating from the running mean by more
than `max(rms, ped_flatness) = 1.0` are dropped, the mean and rms are
recomputed, and the procedure repeats. After convergence the band has
collapsed onto the dominant cluster (146 / 147 ADC) and the outliers
(values 143 / 144 / 145 / 150) are off the kept set. `ped_flatness` sets
the noise floor below which the band stops shrinking — without it, a
particularly quiet channel could iterate down to a band of zero width and
reject everything except the modal value.

> Note: the demo runs the procedure on the *raw* samples for readability
> (integer values plot cleanly); the C++ runs it on the smoothed buffer,
> which is why the converged mean shown here (146.50) is slightly higher
> than the analyzer's reported value (145.61). The kept/rejected pattern
> is the same in both cases.

**Right — `int_tail_ratio`.** The integration walks outward from the
peak and stops when the pedsub waveform falls below `r × peak height`.
For this slow tail:

| `int_tail_ratio` | window | samples | integral |
|---:|:---:|---:|---:|
| 0.20 | [30, 41] | 12 | 8376 |
| 0.10 (default) | [30, 47] | 18 | 9477 |
| 0.05 | [30, 57] | 28 | 10332 |

The default of 0.10 is a good compromise: it captures the prompt peak
plus the first ~70 ns of the tail, missing only the slow scintillation
component (which is also where pile-up from the next event lives).
Smaller `r` recovers more tail energy but increases sensitivity to
baseline drift and downstream pulses.

**Smoothing — `smooth_order`.** On the bright pulse above, smoothing is
invisible at the scale of a 1247 ADC peak. It earns its keep on
small-signal channels where the per-sample fluctuation is comparable to
the pulse height. The figure below uses a different waveform — a small
~24 ADC bump on a baseline with ±3 ADC zig-zag — to show what the kernel
does:

![smoothing](plots/fig5_smoothing.png)

| `smooth_order` | spurious local maxima above +2 ADC | peak height (smoothed) |
|---:|---:|---:|
| 1 (raw) | 6 | 169 |
| 2 (default) | 3 | 166 |
| 4 | 1 | 162 |

`smooth_order = 1` keeps the raw zig-zag — six local maxima clear
+2 ADC, the peak finder needs the height-above-baseline +
peak-overlap rejections to find the real one.  `smooth_order = 2`
collapses the zig-zag without visibly attenuating the pulse.
`smooth_order = 4` removes essentially all baseline structure but
starts to clip the peak by ~7 ADC, so its use should be reserved for
very low-S/N channels where peak-finding robustness is worth more
than peak-height fidelity.

The remaining parameters affect the bright-pulse example only
marginally:

- `threshold = 5` × the post-clipping rms (0.45) gives 2.25 ADC, below
  the floor; `min_threshold = 3` wins. For a noisier channel with rms
  ≳ 1 ADC the `5·rms` rule kicks in.
- `min_peak_ratio` only matters when two peaks share an integration
  range — single-pulse waveforms never trigger it.

## Firmware emulator — `Fadc250FwAnalyzer`

The firmware analyzer reproduces the on-board pipeline so we can compare
offline analysis against firmware-reported values without re-running the
DAQ. The full algorithm spec lives in
[`docs/clas_fadc/FADC250_algorithms.md`](../../clas_fadc/FADC250_algorithms.md);
this section is a parameter-by-parameter walk-through.

### Parameters

Parameters live under the `fadc250_firmware` block in `daq_config.json`.
**`NSB` and `NSA` are in nanoseconds**, floored to whole 4 ns samples
inside the analyzer; everything else is unitless or in ADC counts.

| field | unit | role |
|---|---|---|
| `TET` | ADC counts | Trigger Energy Threshold above pedestal. A pulse is rejected if `Vpeak − Vmin ≤ TET`. |
| `NSB` | ns | Window before threshold crossing (Mode 2 integral). Floored to whole samples (`NSB / CLK_NS`). |
| `NSA` | ns | Window after threshold crossing. Same flooring as `NSB`. |
| `NPEAK` (= `MAX_PULSES`) | — | Max pulses kept per channel per readout window (1..4). |
| `NSAT` | samples | Consecutive-above-TET requirement after Tcross — rejects single-sample spikes. `NSAT=1` reproduces the legacy Mode 3 single-sample crossing. |
| `NPED` | samples | Number of leading samples summed for the `Vnoise` estimate. |
| `MAXPED` | ADC counts | Online outlier-rejection threshold: pedsub samples whose deviation exceeds `MAXPED` are dropped from the `Vnoise` sum. `0` disables. |
| `CLK_NS` | ns | Sample period (4 ns at 250 MHz). |

Run-config defaults (current `daq_config.json`):

```json
"TET": 10.0, "NSB": 8, "NSA": 128, "NPEAK": 1,
"NSAT": 4, "NPED": 3, "MAXPED": 1, "CLK_NS": 4.0
```

### Pipeline

Step-by-step, applied to the example waveform:

**1. Pedestal estimate (`Vnoise`).** Mean of the first `NPED = 3` samples
with `MAXPED = 1` outlier filter (drop any sample whose deviation from the
running mean exceeds 1 ADC). For our trace: `(146+147+144)/3 = 145.67`,
sample 1 (147) is filtered, refined mean = `145.0`.

**2. Pulse search.** `Vmin = Vnoise`. Walk the buffer starting at sample
`NPED`. The first pulse is detected as soon as a sample exceeds `Vnoise`
and walks monotonically up to a local maximum.

**3. Acceptance.** `Vpeak = 1393`. The pedestal-subtracted height is
`Vpeak − Vmin = 1247 ≫ TET = 10` → accepted.

**4. Tcross.** First leading-edge sample whose pedsub value exceeds `TET`:
sample 30, since `637 − 146 = 491 > 10`.

**5. NSAT gate.** `NSAT = 4` → samples 30, 31, 32, 33 must all be > TET.
They are (491, 1221, 1247, 1093) → accepted. With `NSAT = 1`, this gate is
a no-op.

**6. TDC — `Va`, bracket, fine time.**

```
Va  = Vmin + (Vpeak − Vmin) / 2
    = 146 + (1393 − 146) / 2
    = 769.5
```

Find the bracket on the rising edge: smallest `k` with `s[k] ≥ Va`. Here
`s[30] = 637 < 769.5`, `s[31] = 1367 ≥ 769.5` → `k = 31`. So
`Vba = s[30] = 637`, `Vaa = s[31] = 1367`. Fine time:

```
fine = round( (Va − Vba) / (Vaa − Vba) × 64 )
     = round( (769.5 − 637) / (1367 − 637) × 64 )
     = round( 0.1815 × 64 ) = 12

coarse     = k − 1 = 30
time_units = coarse·64 + fine = 1932    (LSB = 62.5 ps)
time_ns    = time_units × CLK_NS / 64 = 120.75 ns
```

Visually (`fig2`, left panel): the dot-dash `Va` line crosses the rising
edge between the diamond `Vba` and the square `Vaa` markers. The fine-time
arrow points from the `Vba` sample to the interpolated zero-crossing.

**7. Mode-2 integral.** Window `[cross − NSB_s, cross + NSA_s]` where
`NSB_s = NSB/4 = 2`, `NSA_s = NSA/4 = 32`. So `[28, 62]`, i.e. 35 samples
= 140 ns. The integrand is the pedestal-subtracted waveform `s′ = max(0,
s − Vnoise)`.

```
Σ s′[28..62] = 10589  (pedsub ADC·sample)
```

The shaded band in `fig2` (right panel) is exactly this sum.

![firmware](plots/fig2_firmware_analysis.png)

**8. Quality bitmask.** `0x00` = `Q_DAQ_GOOD`. Set bits would indicate:

| bit | flag | condition |
|---|---|---|
| `1 << 0` | `Q_DAQ_PEAK_AT_BOUNDARY` | peak landed on the last sample |
| `1 << 1` | `Q_DAQ_NSB_TRUNCATED`   | `cross − NSB_s < 0`, window clipped |
| `1 << 2` | `Q_DAQ_NSA_TRUNCATED`   | `cross + NSA_s ≥ N`, window clipped |
| `1 << 3` | `Q_DAQ_VA_OUT_OF_RANGE` | `Va` not bracketed on the rising edge (numerical edge case) |

## Side-by-side comparison

| field | waveform (`WaveAnalyzer`) | firmware (`Fadc250FwAnalyzer`) |
|---|---|---|
| pedestal | 145.61 ± 0.45 (28/30 samples, median/MAD seed + σ-clip) | 145.0 (3 samples, MAXPED filter) |
| time | 126.6 ns (raw peak + quadratic sub-sample δ) | 120.75 ns (TDC mid-amplitude interp.) |
| height | 1247 ADC (raw − ped) | 1247 ADC (`Vpeak − Vmin`) |
| integral window | [29, 48] (20 samples, tail-driven, INCLUSIVE) | [28, 62] (35 samples, fixed NSB/NSA) |
| integral | 9600 | 10589 |

The waveform analyzer's *time* is the peak vertex via 3-point quadratic
fit (sub-sample, ≪ 1 ns precision for clean peaks); the firmware's
*time* is mid-amplitude on the rising edge with 62.5 ps LSB — they
are intentionally different observables.

The firmware's wider window (140 ns vs the waveform analyzer's 80 ns)
picks up more of the slow scintillation tail. With `NSA = 128 ns` the
window stops at sample 62; the rest of the tail (samples 63..99) is
excluded.

## Reproducing the plots

Both algorithms are re-implemented in
[`scripts/plot_wave_analysis.py`](scripts/plot_wave_analysis.py)
(NumPy + Matplotlib only).

```bash
cd docs/technical_notes/waveform_analysis
python scripts/plot_wave_analysis.py
```

Regenerates `plots/fig1_overview.png`, `plots/fig2_firmware_analysis.png`,
`plots/fig3_soft_analysis.png`, `plots/fig4_soft_parameters.png`,
`plots/fig5_smoothing.png` and prints the numeric results above.

## See also

- [`docs/clas_fadc/FADC250_algorithms.md`](../../clas_fadc/FADC250_algorithms.md)
  — full firmware algorithm spec with manual cross-references
- [`prad2dec/include/WaveAnalyzer.h`](../../../prad2dec/include/WaveAnalyzer.h),
  [`WaveAnalyzer.cpp`](../../../prad2dec/src/WaveAnalyzer.cpp) — C++ source
- [`prad2dec/include/Fadc250FwAnalyzer.h`](../../../prad2dec/include/Fadc250FwAnalyzer.h),
  [`Fadc250FwAnalyzer.cpp`](../../../prad2dec/src/Fadc250FwAnalyzer.cpp) — C++ source
- [`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md) — branch
  layout for the replay tree (where both analyzer outputs land)
