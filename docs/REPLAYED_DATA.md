# Replayed Data Trees

Two ROOT trees are produced by the replay tools:

| Tool | Tree name | Contents |
|---|---|---|
| `prad2ana_replay_rawdata` | `events` | Raw FADC250 waveforms + GEM strip data + optional per-channel peak analysis |
| `prad2ana_replay_recon`   | `recon`  | HyCal clusters + GEM hits (lab frame) + HyCal↔GEM matches |

# `events` tree (raw)

Written by `prad2ana_replay_rawdata`.  Per-event scalars and per-channel
arrays sized by `hycal.nch` / `gem.nch`.  The `hycal.*` arrays cover **all**
FADC250 channels (HyCal + Veto + LMS) — distinguish via `hycal.module_type`:

| Value | Type | Source |
|---|---|---|
| 0 | Unknown | — |
| 1 | PbGlass | `hycal_modules.json` `t="PbGlass"` |
| 2 | PbWO4   | `hycal_modules.json` `t="PbWO4"`   |
| 3 | SCINT   | `hycal_modules.json` `t="SCINT"` (Veto V1..V4) |
| 4 | LMS     | `hycal_modules.json` `t="LMS"` (LMSPin, LMS1..3) |

`hycal.module_id` encoding (globally unique):
PbGlass = 1..1156, PbWO4 = 1001..2152, SCINT = 3001..3004,
LMS = 3100 (Pin) / 3101..3103.

## Branches

### Event header — always written

| Branch | Type | Meaning |
|---|---|---|
| `event_num`    | `int`   | Event number (TI / trigger bank) |
| `trigger_type` | `uint8` | Main trigger type |
| `trigger_bits` | `uint32` | FP trigger inputs (32-bit bitmask) |
| `timestamp`    | `int64` | 48-bit TI timestamp (250 MHz ticks) |
| `ssp_raw`      | `vector<uint32>` | Raw 0xE10C SSP trigger bank words |

### HyCal / Veto / LMS FADC250 — always written

| Branch | Type | Meaning |
|---|---|---|
| `hycal.nch`         | `int`            | Number of FADC250 channels in event |
| `hycal.module_id`   | `uint16[nch]`    | See ranges above |
| `hycal.module_type` | `uint8[nch]`     | Category enum (legend above) |
| `hycal.nsamples`    | `uint8[nch]`     | Samples per channel (≤ 200) |
| `hycal.samples`     | `uint16[nch][200]` | Raw 12-bit ADC samples |
| `hycal.gain_factor` | `float[nch]`     | Gain correction (1.0 for SCINT/LMS) |

### Soft-analyzer outputs — only with `-p`

Local-maxima search with median/MAD-bootstrapped iterative-outlier-rejection
pedestal (`fdec::WaveAnalyzer`).  Up to `MAX_PEAKS = 8` peaks per channel.
Without `-p`, the soft analyzer is skipped entirely (the pedestal estimate
that the firmware analyzer uses is also gated on `-p`).

| Branch | Type | Meaning |
|---|---|---|
| `hycal.ped_mean`       | `float[nch]`      | Pedestal mean (post-rejection) |
| `hycal.ped_rms`        | `float[nch]`      | Pedestal RMS  (post-rejection) |
| `hycal.ped_nused`      | `uint8[nch]`      | # samples surviving outlier rejection |
| `hycal.ped_quality`    | `uint8[nch]`      | `Q_PED_*` bitmask (legend below) |
| `hycal.ped_slope`      | `float[nch]`      | Linear drift across surviving samples (ADC/sample) |
| `hycal.npeaks`         | `uint8[nch]`      | Soft peaks found |
| `hycal.peak_height`    | `float[nch][8]`   | Peak height above pedestal |
| `hycal.peak_time`      | `float[nch][8]`   | Peak time (ns) — quadratic-vertex sub-sample interpolation around the raw peak |
| `hycal.peak_integral`  | `float[nch][8]`   | Peak integral over `[peak.left, peak.right]` (INCLUSIVE) |
| `hycal.peak_quality`   | `uint8[nch][8]`   | `Q_PEAK_*` bitmask: `1` = `Q_PEAK_PILED` (this peak's integration window touches/overlaps an adjacent peak's, within `cfg.peak_pileup_gap` samples) |

`hycal.ped_quality` bits (defined in `prad2dec/include/Fadc250Data.h`):

| Bit | Flag | Meaning |
|---|---|---|
| `0`     | `Q_PED_GOOD`             | clean pedestal — no flags set |
| `1<<0`  | `Q_PED_NOT_CONVERGED`    | `ped_max_iter` exhausted, kept-mask still moving |
| `1<<1`  | `Q_PED_FLOOR_ACTIVE`     | `rms < ped_flatness` — `ped_flatness` was the active band (typical for very quiet channels; informational) |
| `1<<2`  | `Q_PED_TOO_FEW_SAMPLES`  | < 5 samples survived rejection — estimate is unreliable |
| `1<<3`  | `Q_PED_PULSE_IN_WINDOW`  | a peak landed inside the pedestal window we used |
| `1<<4`  | `Q_PED_OVERFLOW`         | a raw sample in the window hit the 12-bit overflow (4095) |
| `1<<5`  | `Q_PED_TRAILING_WINDOW`  | adaptive logic preferred trailing samples over the leading window (informational, not a problem flag) |

A clean event filter is just `hycal.ped_quality == 0`; the
`PULSE_IN_WINDOW` and `TRAILING_WINDOW` flags are useful diagnostics for
events that fail it.

### Firmware (FADC250 Mode 1/2/3) peaks — only with `-p`

Bit-faithful emulation of the JLab FADC250 firmware (Hall-D V3
extensions: NSAT/NPED/MAXPED).  Configured via the
[`fadc250_waveform.firmware`](../database/daq_config.json) block in
`daq_config.json`.  See
[`docs/clas_fadc/FADC250_algorithms.md`](../docs/clas_fadc/FADC250_algorithms.md)
for the algorithm spec.

| Branch | Type | Meaning |
|---|---|---|
| `hycal.daq_npeaks`        | `uint8[nch]`    | Firmware pulses (≤ NPEAK) |
| `hycal.daq_peak_vp`       | `float[nch][8]` | Vpeak (pedestal-subtracted) |
| `hycal.daq_peak_integral` | `float[nch][8]` | Σ over [cross−NSB, cross+NSA] (Mode 2) |
| `hycal.daq_peak_time`     | `float[nch][8]` | Mid-amplitude time, ns (62.5 ps LSB) |
| `hycal.daq_peak_cross`    | `int[nch][8]`   | Tcross sample index |
| `hycal.daq_peak_pos`      | `int[nch][8]`   | Sample of Vpeak |
| `hycal.daq_peak_coarse`   | `int[nch][8]`   | 4-ns clock index of Vba (10-bit) |
| `hycal.daq_peak_fine`     | `int[nch][8]`   | Fine bits 0..63 (6-bit) |
| `hycal.daq_peak_quality`  | `uint8[nch][8]` | `Q_DAQ_*` bitmask: `1` = peak@boundary, `2` = NSB-trunc, `4` = NSA-trunc, `8` = Va out-of-range |

### GEM strips — always written

| Branch | Type | Meaning |
|---|---|---|
| `gem.nch`         | `int`             | Number of GEM strips fired |
| `gem.mpd_crate`   | `uint8[nch]`      | MPD crate ID |
| `gem.mpd_fiber`   | `uint8[nch]`      | MPD fiber ID |
| `gem.apv`         | `uint8[nch]`      | APV ADC channel |
| `gem.strip`       | `uint8[nch]`      | Strip number on the APV |
| `gem.ssp_samples` | `int16[nch][6]`   | 6 SSP time samples per strip |

# `recon` tree (reconstructed)

Written by `prad2ana_replay_recon`.  HyCal clustering, GEM hit
reconstruction, and per-cluster HyCal↔GEM matching.  All positions are in
the lab frame (target-centered, beam-aligned, mm).  Trigger filter applied
upstream — only physics events reach the tree.

### Event header

| Branch | Type | Meaning |
|---|---|---|
| `event_num`    | `int`            | Event number |
| `trigger_type` | `uint8`          | Main trigger type |
| `trigger_bits` | `uint32`         | FP trigger bits (32-bit bitmask) |
| `timestamp`    | `int64`          | 48-bit TI timestamp (250 MHz ticks) |
| `total_energy` | `float`          | Σ HyCal cluster energy (MeV) |
| `ssp_raw`      | `vector<uint32>` | Raw 0xE10C SSP trigger bank words |

### HyCal clusters

`n_clusters` ≤ 100.  Lab frame (target-centered, beam-aligned, mm).

| Branch | Type | Meaning |
|---|---|---|
| `n_clusters` | `int`                | Number of reconstructed clusters |
| `cl_x`       | `float[n_clusters]`  | Cluster x at HyCal face + shower depth |
| `cl_y`       | `float[n_clusters]`  | Cluster y |
| `cl_z`       | `float[n_clusters]`  | `hycal_z` + shower-depth correction |
| `cl_energy`  | `float[n_clusters]`  | Cluster energy (MeV) |
| `cl_nblocks` | `uint8[n_clusters]`  | Modules summed into the cluster |
| `cl_center`  | `uint16[n_clusters]` | Center module ID (same numbering as `hycal.module_id`) |
| `cl_flag`    | `uint32[n_clusters]` | HyCal cluster flags (bit field) |

### Per-cluster HyCal↔GEM match (all 4 GEMs)

For each HyCal cluster, the closest GEM hit on each of the 4 detectors
within the matching window is recorded (or `0`/`NaN` if none).  Use this
when you want a fixed-shape `[n_clusters][4]` view.

| Branch | Type | Meaning |
|---|---|---|
| `matchFlag` | `uint32[n_clusters]`    | Per-cluster match flags (which GEMs matched) |
| `matchGEMx` | `float[n_clusters][4]`  | Matched GEM x (det 0..3) |
| `matchGEMy` | `float[n_clusters][4]`  | Matched GEM y |
| `matchGEMz` | `float[n_clusters][4]`  | Matched GEM z |

### Quick-access matched pairs (clusters with ≥2 GEMs matched)

`match_num` ≤ 100. Convenient `[match_num][2]` view for analyses that only
care about clusters confirmed on at least two GEM planes.

| Branch | Type | Meaning |
|---|---|---|
| `match_num` | `int`                  | Number of clusters with ≥2 GEMs matched |
| `mHit_E`    | `float[match_num]`     | HyCal cluster energy (MeV) |
| `mHit_x`    | `float[match_num]`     | HyCal cluster x |
| `mHit_y`    | `float[match_num]`     | HyCal cluster y |
| `mHit_z`    | `float[match_num]`     | HyCal cluster z |
| `mHit_gx`   | `float[match_num][2]`  | First 2 matched GEM x |
| `mHit_gy`   | `float[match_num][2]`  | First 2 matched GEM y |
| `mHit_gz`   | `float[match_num][2]`  | First 2 matched GEM z |
| `mHit_gid`  | `float[match_num][2]`  | det_id (0..3) of those 2 GEM hits |

### GEM reconstructed hits

`n_gem_hits` ≤ 400.  All hits across all 4 detectors, lab frame.

| Branch | Type | Meaning |
|---|---|---|
| `n_gem_hits`   | `int`                  | Total GEM hits across all detectors |
| `det_id`       | `uint8[n_gem_hits]`    | GEM detector ID (0..3) |
| `gem_x`        | `float[n_gem_hits]`    | Hit x (lab) |
| `gem_y`        | `float[n_gem_hits]`    | Hit y |
| `gem_z`        | `float[n_gem_hits]`    | Hit z (per-detector plane) |
| `gem_x_charge` | `float[n_gem_hits]`    | Total ADC of the X cluster |
| `gem_y_charge` | `float[n_gem_hits]`    | Total ADC of the Y cluster |
| `gem_x_peak`   | `float[n_gem_hits]`    | Max-strip ADC, X plane |
| `gem_y_peak`   | `float[n_gem_hits]`    | Max-strip ADC, Y plane |
| `gem_x_size`   | `uint8[n_gem_hits]`    | Strips in X cluster |
| `gem_y_size`   | `uint8[n_gem_hits]`    | Strips in Y cluster |
| `gem_x_mTbin`  | `uint8[n_gem_hits]`    | Time-sample bin of max-ADC strip, X |
| `gem_y_mTbin`  | `uint8[n_gem_hits]`    | Time-sample bin of max-ADC strip, Y |

### Veto + LMS (peak summaries)

Lightweight tag of the best soft peak per Veto / LMS channel — full
waveforms live in the `events` tree, not here.

| Branch | Type | Meaning |
|---|---|---|
| `veto_nch`          | `int`              | Number of Veto channels with data |
| `veto_id`           | `uint8[veto_nch]`  | 0..3 (V1..V4) |
| `veto_npeaks`       | `int[veto_nch]`    | Soft peaks found |
| `veto_peak_time`    | `float[veto_nch][8]` | Peak time (ns) |
| `veto_peak_integral`| `float[veto_nch][8]` | Peak integral |
| `lms_nch`           | `int`              | Number of LMS channels with data |
| `lms_id`            | `uint8[lms_nch]`   | 0=Pin, 1..3 = LMS1..3 |
| `lms_npeaks`        | `int[lms_nch]`     | Soft peaks found |
| `lms_peak_time`     | `float[lms_nch][8]` | Peak time (ns) |
| `lms_peak_integral` | `float[lms_nch][8]` | Peak integral |

# Run example

```bash
ssh clasrun@clonfarm11
source ~/prad2_daq/prad2_env.csh
cd /data/replay_raw/
prad2ana_replay_rawdata /data/evio/data/prad_024154/prad_024154.evio.00000 -o ./ -p
```

Output: `/data/replay_raw/prad_024154.00000_raw.root`.

DAQ-emulation knobs (`TET` / `NSB` / `NSA` / `NPEAK` / `NSAT` / `NPED` /
`MAXPED`) are read from the `fadc250_waveform.firmware` block in
`daq_config.json`; the offline soft analyzer (`WaveAnalyzer`) is configured
from the sibling `fadc250_waveform.analyzer` block — override either to
match the actual run.
