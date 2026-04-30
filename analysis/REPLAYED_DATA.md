# Replayed Raw Data — `events` tree

ROOT tree written by `prad2ana_replay_rawdata`.  Per-event scalars and
per-channel arrays sized by `hycal.nch` / `gem.nch`.  The `hycal.*` arrays
cover **all** FADC250 channels (HyCal + Veto + LMS) — distinguish via
`hycal.module_type`:

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
| `hycal.nsamples`    | `int[nch]`       | Samples per channel (≤ 200) |
| `hycal.samples`     | `uint16[nch][200]` | Raw 12-bit ADC samples |
| `hycal.ped_mean`    | `float[nch]`     | Soft-analyzer pedestal mean |
| `hycal.ped_rms`     | `float[nch]`     | Soft-analyzer pedestal RMS |
| `hycal.gain_factor` | `float[nch]`     | Gain correction (1.0 for SCINT/LMS) |

### Soft-analyzer peaks — only with `-p`

Local-maxima search with iterative-outlier-rejection pedestal
(`fdec::WaveAnalyzer`).  Up to `MAX_PEAKS = 8` peaks per channel.

| Branch | Type | Meaning |
|---|---|---|
| `hycal.npeaks`         | `int[nch]`        | Soft peaks found |
| `hycal.peak_height`    | `float[nch][8]`   | Peak height above pedestal |
| `hycal.peak_time`      | `float[nch][8]`   | Peak time (ns) |
| `hycal.peak_integral`  | `float[nch][8]`   | Peak integral |

### Firmware (FADC250 Mode 1/2/3) peaks — only with `-p`

Bit-faithful emulation of the JLab FADC250 firmware (Hall-D V3
extensions: NSAT/NPED/MAXPED).  Configured via the
[`fadc250_firmware`](../database/daq_config.json) block in
`daq_config.json`.  See
[`docs/clas_fadc/FADC250_algorithms.md`](../docs/clas_fadc/FADC250_algorithms.md)
for the algorithm spec.

| Branch | Type | Meaning |
|---|---|---|
| `hycal.daq_npeaks`        | `int[nch]`      | Firmware pulses (≤ NPEAK) |
| `hycal.daq_peak_vp`       | `float[nch][8]` | Vpeak (pedestal-subtracted) |
| `hycal.daq_peak_integral` | `float[nch][8]` | Σ over [cross−NSB, cross+NSA] (Mode 2) |
| `hycal.daq_peak_time`     | `float[nch][8]` | Mid-amplitude time, ns (62.5 ps LSB) |
| `hycal.daq_peak_cross`    | `int[nch][8]`   | Tcross sample index |
| `hycal.daq_peak_pos`      | `int[nch][8]`   | Sample of Vpeak |
| `hycal.daq_peak_coarse`   | `int[nch][8]`   | 4-ns clock index of Vba (10-bit) |
| `hycal.daq_peak_fine`     | `int[nch][8]`   | Fine bits 0..63 (6-bit) |
| `hycal.daq_peak_quality`  | `uint8[nch][8]` | Bitmask: `1` = peak@boundary, `2` = NSB-trunc, `4` = NSA-trunc, `8` = Va out-of-range |

### GEM strips — always written

| Branch | Type | Meaning |
|---|---|---|
| `gem.nch`         | `int`             | Number of GEM strips fired |
| `gem.mpd_crate`   | `uint8[nch]`      | MPD crate ID |
| `gem.mpd_fiber`   | `uint8[nch]`      | MPD fiber ID |
| `gem.apv`         | `uint8[nch]`      | APV ADC channel |
| `gem.strip`       | `uint8[nch]`      | Strip number on the APV |
| `gem.ssp_samples` | `int16[nch][6]`   | 6 SSP time samples per strip |

## Run example

```bash
ssh clasrun@clonfarm11
source ~/prad2_daq/prad2_env.csh
cd /data/replay_raw/
prad2ana_replay_rawdata /data/evio/data/prad_024154/prad_024154.evio.00000 -o ./ -p
```

Output: `/data/replay_raw/prad_024154.00000_raw.root`.

DAQ-emulation knobs (`TET` / `NSB` / `NSA` / `NPEAK` / `NSAT` / `NPED` /
`MAXPED`) are read from the `fadc250_firmware` block in `daq_config.json` —
override there to match the actual run.
