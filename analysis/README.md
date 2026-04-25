# Offline Analysis Tools

Replay and physics analysis for PRad2. **Requires ROOT 6.0+.**

```bash
cmake -B build -DBUILD_ANALYSIS=ON
cmake --build build -j$(nproc)
```

## Tools

### Replay

> All analysis executables are installed with a `prad2ana_` prefix to avoid
> name collisions in a shared bindir.  The tool names below refer to the
> built binary — e.g. `prad2ana_replay_rawdata`.

**replay_rawdata** — Single EVIO file to ROOT tree with per-channel waveform data.
```bash
prad2ana_replay_rawdata <input.evio> [-o output.root] [-n max_events] [-p]
```

**replay_rawdata_m** — Multi-file, multi-threaded version of `replay_rawdata`. Processes all EVIO segments in a directory.
```bash
prad2ana_replay_rawdata_m <evio_dir> [-f max_files] [-n max_events] [-p] [-j num_threads] [-D daq_config.json] [-o merged.root]
```

**replay_recon** — HyCal reconstruction replay with clustering and per-module energy histograms.
```bash
prad2ana_replay_recon <input.evio> [-o output.root] [-c config.json] [-D daq_config.json] [-n N]
```

**replay_recon_m** — Multi-file, multi-threaded version of `replay_recon`. Supports GEM pedestal and zero-suppression options.
```bash
prad2ana_replay_recon_m <evio_dir> [-f max_files] [-n max_events] [-p] [-j num_threads] [-D daq_config.json] [-g gem_pedestal.json] [-z zerosup_threshold] [-o merged.root]
```
- `-p`  read PRad-I data format (no GEM)

### Calibration

**epCalib** — Elastic e-p calibration. Fits the elastic peak per module from rawdata ROOT files (peak mode) and writes gain correction constants.
```bash
prad2ana_epCalib <input.root> [-o output_calib_file] [-D daq_config.json] [-n max_events]
```

### Physics Analysis

**analysis_example** — Example offline analysis reading reconstructed ROOT trees. Fills energy, hit-position, and Moller-event histograms with optional GEM matching.
```bash
prad2ana_analysis_example <input_recon.root> [-o output.root] [-n max_events]
```

**cosmic_test** — Cosmic-ray analysis tool for commissioning. Reads raw waveform data and produces per-channel signal distributions.
```bash
prad2ana_cosmic_test <input.root> [-o output.root] [-D daq_config.json] [-n max_events]
```

## ACLiC Scripts (`scripts/`)

ROOT macros that compile against `libprad2dec` / `libprad2det` via ACLiC. They share one prelude — `rootlogon.C` — that resolves include paths, finds the static libraries in your build directory, and sets `PRAD2_DATABASE_DIR`. Run it once per ROOT session:

```bash
cd build
root -l ../analysis/scripts/rootlogon.C
```

Then `.x` any of the macros below. All accept a final `daq_config` argument that defaults to `$PRAD2_DATABASE_DIR/daq_config.json`.

### gem_raw_dump.C

Smallest GEM example — opens an EVIO file, finds every GEM raw bank (tags from `daq_config.json`'s `bank_tags.ssp_raw`, currently `0xE10C` / `0x0DE9`), and prints the first N raw 32-bit words per bank per event. No decoding. Useful for verifying firmware layout or feeding raw words into a custom parser.

```bash
.x ../analysis/scripts/gem_raw_dump.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.00000", 5, 8)
# args: evio_path, max_events (0=all), n_words_show
```

### gem_clusters_to_root.C

Full GEM analysis pipeline → ROOT tree.  Per event runs `EvChannel.DecodeEvent` → `GemSystem.ProcessEvent` (pedestal + CM + ZS) → `GemSystem.Reconstruct(GemCluster)`, then writes a flat `TTree` with all 1D strip clusters and their constituent strip hits including the full 6 time samples per hit.

Tree layout (one entry per physics event):

| Group   | Branches | Sized by |
|---------|----------|----------|
| event   | `event_num`, `trigger_bits` | scalar |
| cluster | `cl_det`, `cl_plane`, `cl_position`, `cl_peak_charge`, `cl_total_charge`, `cl_max_tb`, `cl_cross_talk`, `cl_nhits`, `cl_first` | `ncl` |
| hit     | `hit_cl`, `hit_strip`, `hit_position`, `hit_charge`, `hit_max_tb`, `hit_cross_talk`, `hit_ts0`…`hit_ts5` | `nhits` |

`cl_first` + `cl_nhits` slice the hit arrays per cluster; `hit_cl` is the back-pointer for hit→cluster joins.

```bash
# simplest — pedestal/CM auto-discovered from database/config.json's
# runinfo (latest entry; same path the live monitor uses):
.x ../analysis/scripts/gem_clusters_to_root.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.00000", \
    "gem_clusters_023867.root")

# explicit override (paths relative to PRAD2_DATABASE_DIR or absolute):
.x ../analysis/scripts/gem_clusters_to_root.C+( \
    "/data/.../prad_023867.evio.00000", "out.root", \
    "gem_peds/peds_23867.txt", "gem_peds/cm_23867.txt")
```

Full signature: `gem_clusters_to_root(evio_path, out_path, gem_ped_file=nullptr, gem_cm_file=nullptr, max_events=0, run_num=-1, daq_config=nullptr, gem_map_file=nullptr)`. Pass `nullptr` for the ped/cm args to auto-discover via `database/config.json` → `runinfo` → `LoadRunConfig(run_num)`. The script reads the GEM crate-remap from `daq_config.json` so it always matches what the live monitor reconstructs.

### tagger_hycal_correlation.C

Two-phase study of T10R↔E49…E58 tagger pairs vs HyCal PbWO4 sums. Phase 1 caches per-event TDC tuples and fits a Gaussian to each ΔT spectrum; Phase 2 applies an N-σ timing cut per pair and fills global W-sum height/integral histograms for events with at least one matched pair plus a W channel above threshold. Outputs a 12-panel summary canvas.

```bash
.x ../analysis/scripts/tagger_hycal_correlation.C+( \
    "/data/stage6/prad_023686/prad_023686.evio.00000", \
    "tagger_wsum_corr.root", 500000)
# args: evio_path, out_path, max_events
```

### lms_alpha_normalize.C

Scans a run's EVIO files (`prad_{run}.evio.00000`–`99999`), selects only **LMS** (bit 24) and **Alpha** (bit 25) trigger events via `trigger_bits`, and normalizes HyCal LMS signals using the Alpha source as a gain reference.

Uses the project's decoder (`EvChannel`, `WaveAnalyzer`, `DaqConfig`) — no manual FADC parsing needed. Channel identity (HyCal module vs LMS ref) is read from `daq_map.json`.

| Trigger | What fires | Purpose |
|---------|-----------|---------|
| LMS     | All HyCal modules + LMS 1/2/3 | Monitor module response via LED/laser pulser |
| Alpha   | LMS 1/2/3 only | Provide stable gain reference (Am-241 source) |

Normalization per module *i* (averaged across 3 references *j*):
```
norm_i = integral_i × mean( alpha_ref_j / lms_ref_j )
```
LMS events before the first Alpha event in the run are skipped. Each LMS event uses the most recent Alpha reading.

```bash
# one-liner, after rootlogon:
root -l rootlogon.C 'lms_alpha_normalize.C+("/path/to/data", 1234)'

# with explicit config overrides:
root -l rootlogon.C 'lms_alpha_normalize.C+("/path/to/data", 1234, "daq_config.json", "daq_map.json")'
```

**Output:** `lms_alpha_run{N}.root` (per-module normalized LMS, reference time-series TGraphs) and a 6-panel summary PNG.

## Adding a Tool

Create `tools/my_tool.cpp`, then add to `CMakeLists.txt`:
```cmake
add_analysis_tool(my_tool tools/my_tool.cpp)
```

Shared sources (`Replay.cpp`, `PhysicsTools.cpp`, `MatchingTools.cpp`) and dependencies linked automatically.

## Contributors
Yuan Li, Weizhi Xiong — Shandong University
Chao Peng - Argonne National Laboratory
