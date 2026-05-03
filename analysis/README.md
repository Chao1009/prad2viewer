# Offline Analysis Tools

Replay and physics analysis for PRad2. **Requires ROOT 6.0+.**

```bash
cmake -B build -DBUILD_ANALYSIS=ON
cmake --build build -j$(nproc)
cmake --install build --prefix /your/install/prefix    # optional
```

Builds three things:

| Artifact | Where |
|----------|-------|
| `prad2ana_*` executables (replay, calibration, …) | `<build>/bin/`, installed to `<prefix>/bin/` |
| **`libprad2ana.a`** — static library exposing `analysis::*` (Replay, PhysicsTools, MatchingTools) | `<build>/analysis/`, installed to `<prefix>/lib/` |
| Headers (`PhysicsTools.h`, `MatchingTools.h`, `Replay.h`, `ConfigSetup.h`) | installed to `<prefix>/include/prad2ana/` |

The library is what makes ACLiC scripts work in install mode — without it, `analysis::*` symbols would be unresolved at link time.

## Tools

### Replay

> All analysis executables are installed with a `prad2ana_` prefix to avoid
> name collisions in a shared bindir.  The tool names below refer to the
> built binary — e.g. `prad2ana_replay_rawdata`.

Both replay tools are multi-threaded — point them at one or more EVIO
files (or a directory of segments) and tune parallelism with `-j`.  Output
goes to `<output_dir>/<input_stem>_raw.root` (rawdata) or `_recon.root`
(recon), one ROOT file per input EVIO segment.

**replay_rawdata** — EVIO → per-channel waveform tree (`events`).
```bash
prad2ana_replay_rawdata <evio_or_dir> [more...] -o <output_dir> \
    [-f max_files] [-n max_events] [-j num_threads] \
    [-c daq_config.json] [-d hycal_map.json] [-p]
```
- `-p` adds peak branches: `hycal.peak_*` (soft-analyzer) **and**
  `hycal.daq_peak_*` (firmware Mode 1/2/3 emulation).  Without `-p`, only
  raw samples + module info are written; the soft analyzer is skipped
  entirely on the per-channel hot path.
- See [REPLAYED_DATA.md](../docs/REPLAYED_DATA.md) for the full branch list.

**replay_recon** — EVIO → reconstructed tree (`recon`): HyCal clustering,
GEM hits, HyCal↔GEM straight-line matching.
```bash
prad2ana_replay_recon <evio_or_dir> [more...] -o <output_dir> \
    [-f max_files] [-n max_events] [-j num_threads] \
    [-c daq_config.json] [-d hycal_map.json] \
    [-g gem_pedestal.json] [-z zerosup_threshold] [-p]
```
- `-p` here selects **PRad-I data format** (no GEM, ADC1881M Fastbus
  pedestals).  Different semantics from the `-p` on `replay_rawdata` —
  this is a legacy compatibility flag.
- `-g` overrides the GEM pedestal file from `runinfo`; `-z` overrides
  the zero-suppression sigma threshold.
- See [REPLAYED_DATA.md](../docs/REPLAYED_DATA.md) for the full branch list.

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

ROOT macros that compile against `libprad2dec` / `libprad2det` / `libprad2ana` via ACLiC. They share one prelude — `rootlogon.C` — that auto-detects whether you're in a build tree or an install tree and configures include paths + linker line accordingly.

**Build-tree mode** (preferred, picks up the freshest libs):
```bash
cd build
root -l ../analysis/scripts/rootlogon.C        # CMakeCache.txt in cwd → build mode
# or from anywhere:
PRAD2_BUILD_DIR=/path/to/build root -l rootlogon.C
```

**Install-tree mode** (after `cmake --install`):
```bash
source <prefix>/bin/prad2_setup.sh             # exports PRAD2_DATABASE_DIR
root -l <prefix>/share/prad2evviewer/analysis/scripts/rootlogon.C
```

Each path probe is logged as `[+] tag : path` (resolved) or `[-] tag : path` (skipped/missing) so a failed setup is one glance to debug. Set `PRAD2_ROOTLOGON_QUIET=1` to suppress the per-probe lines and keep just the section headers.

| Env var | What it overrides |
|---------|-------------------|
| `PRAD2_BUILD_DIR`        | build dir if not the cwd |
| `PRAD2_DATABASE_DIR`     | install-mode prefix anchor (set by `prad2_setup.sh`) |
| `PRAD2_EVIO_LIB`         | explicit `libevio.a` path (skips all evio probes) |
| `PRAD2_CODA_ROOT`        | non-default Hall-B CODA install root |
| `PRAD2_ROOTLOGON_QUIET`  | suppress per-probe `[+]/[-]` lines |

Then `.x` any of the macros below.

### gem_raw_dump.C

Smallest GEM example — opens an EVIO file, finds every GEM raw bank (tags from `daq_config.json`'s `bank_tags.ssp_raw`, currently `0xE10C` / `0x0DE9`), and prints the first N raw 32-bit words per bank per event. No decoding. Useful for verifying firmware layout or feeding raw words into a custom parser.

```bash
.x ../analysis/scripts/gem_raw_dump.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.00000", 5, 8)
# args: evio_path, max_events (0=all), n_words_show
```

### gem_hycal_matching.C

Full **HyCal + GEM** reconstruction pipeline with straight-line cluster matching → ROOT tree of matched HC↔GEM pairs and the constituent X/Y GEM strip waveforms.

> **Python counterpart**: `analysis/pyscripts/gem_hycal_matching.py` — same pipeline via `prad2py`, no ROOT, flat TSV/CSV out, identical best-match rule. See [Python counterparts](#python-counterparts-pyscripts) below.

**Trigger filter**: only events with `trigger_bits == 0x100` (production physics trigger) are reconstructed and written. Everything else (LMS / Alpha / cosmic / etc.) is skipped — the summary reports raw physics count vs. kept count.

**Multi-file** — chosen by the input path:
- `prad_023881.evio.*` → **glob mode**: enumerate every sibling `prad_023881.evio.<digits>` in the enclosing directory, process them all into one output tree (suffix-sorted). Warns to stderr about any gap in the suffix sequence — including a missing `.00000` start.
- A directory (e.g. `/data/prad_023881/`) → same enumeration, run number sniffed from the directory name.
- `prad_023881.evio.00000` → **single-file mode**: just that one split.

The summary reports `EVIO files opened : M / N`.

Per event (after the trigger cut):
- `EvChannel.DecodeEvent` → FADC + SSP buffers
- HyCal: `WaveAnalyzer` → `mod->energize` → `HyCalCluster.FormClusters() / ReconstructHits()`
- GEM: `GemSystem.ProcessEvent` (pedestal + CM + ZS) → `Reconstruct(GemCluster)` → 2D X×Y matched hits per detector
- Lab-frame transform via `RotateDetData` / `TransformDetData` (uses `runinfo` geometry for HyCal + each GEM)
- For each HyCal cluster, draw a line from `(0,0,0)` target through the lab-frame centroid (z = `hycal_z` + shower depth); intersect each GEM plane and find the closest GEM hit within `N · σ_total` of the projection
- **Best-match rule** (HyCal cluster as baseline): per (HC cluster, GEM detector) pair, keep at most ONE row — the GEM hit with the smallest 2D residual that's still inside the window. A given GEM hit can win against multiple HC clusters (no GEM-side exclusivity). The Python counterpart uses the same rule.
- For each match, look up the X & Y constituent `StripCluster` on the corresponding plane and copy every strip's full 6-sample waveform

Matching geometry (driven by `reconstruction_config.json:matching`):
```
σ_hc(face) = sqrt((A/sqrt(E_GeV))² + (B/E_GeV)² + C²)   [mm at HyCal face]
             A,B,C = matching.hycal_pos_res
σ_hc@gem   = σ_hc(face) · (z_gem / z_hc)
σ_gem      = matching.gem_pos_res[det_id]               [mm, per detector]
σ_total    = sqrt(σ_hc@gem² + σ_gem²)
match if  |residual| < N · σ_total                      [N defaults to 3, configurable]
```
The actual residual + `σ_total` are stored per match so downstream cuts can be tuned without re-running. The C++ formula lives on `HyCalSystem::PositionResolution(E)` (set via `SetPositionResolutionParams(A, B, C)`); the loader helper is `script_helpers.h::load_matching_config(...)`. Python counterpart: `_common.load_matching_config()` + `_common.hycal_pos_resolution(...)`.

Tree layout (`match`, one entry per physics event):

| Group | Branches | Sized by |
|-------|----------|----------|
| event   | `event_num`, `trigger_bits` | scalar |
| HyCal   | `hc_x`, `hc_y`, `hc_z`, `hc_energy`, `hc_center`, `hc_nblocks`, `hc_flag`, `hc_sigma` | `ncl` |
| match   | `m_hc_idx`, `m_det`, `m_gem_x`, `m_gem_y`, `m_gem_z`, `m_gem_x_charge`, `m_gem_y_charge`, `m_gem_x_size`, `m_gem_y_size`, `m_proj_x`, `m_proj_y`, `m_residual`, `m_sigma_total` | `nmatch` |
| match X cl | `m_xcl_position`, `m_xcl_total`, `m_xcl_peak`, `m_xcl_max_tb`, `m_xcl_first`, `m_xcl_nstrips` | `nmatch` |
| match Y cl | `m_ycl_position`, `m_ycl_total`, `m_ycl_peak`, `m_ycl_max_tb`, `m_ycl_first`, `m_ycl_nstrips` | `nmatch` |
| strips  | `s_match_idx`, `s_plane` (0=X, 1=Y), `s_strip`, `s_position`, `s_charge`, `s_max_tb`, `s_cross_talk`, `s_ts0`…`s_ts5` | `nstrips` |

`m_xcl_first` + `m_xcl_nstrips` slice the strip arrays per matched X cluster (and `m_ycl_*` for Y); `s_match_idx` is the back-pointer hit→match.

Pedestals, common-mode files, and HyCal calibration are auto-discovered from `database/reconstruction_config.json` → `runinfo` (the same path `app_state_init.cpp` uses on the live monitor), so the analysis tree's reconstruction matches what the monitor sees. Pass an explicit path to override.

```bash
# full-run replay — glob discovers every split, warns on gaps:
.x ../analysis/scripts/gem_hycal_matching.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.*", \
    "match_023867.root")

# debug a single split file:
.x ../analysis/scripts/gem_hycal_matching.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.00000", \
    "match_023867_seg0.root")

# tighter matching cut (2σ instead of 3σ default) — 5-arg overload:
.x ../analysis/scripts/gem_hycal_matching.C+( \
    "/data/.../prad_023867.evio.*", "out.root", \
    0L, -1, 2.0f)

# explicit overrides (paths relative to PRAD2_DATABASE_DIR or absolute):
.x ../analysis/scripts/gem_hycal_matching.C+( \
    "/data/.../prad_023867.evio.*", "out.root", \
    "gem_peds/peds_23867.txt", "gem_peds/cm_23867.txt", \
    "calibration/calibration_factor_0.json")
```

Convenience overloads (sidestep a cling default-arg-marshalling SEGV):
- `gem_hycal_matching(evio, out)`
- `gem_hycal_matching(evio, out, max_events)`
- `gem_hycal_matching(evio, out, max_events, run_num)`
- `gem_hycal_matching(evio, out, max_events, run_num, match_nsigma)`

Full 11-arg version (for explicit overrides — pass `""` to auto-discover any path):
`gem_hycal_matching(evio_path, out_path, gem_ped_file, gem_cm_file, hc_calib_file, max_events, run_num, match_nsigma, daq_config, gem_map_file, hc_map_file)`.

### plot_hits_at_hycal.C

Side-by-side 2D occupancy maps of **GEM hits projected to the HyCal surface** (left) and **HyCal cluster centroids on the HyCal surface** (right). Both plots share the same lab/target-centered, beam-aligned frame at z = `hycal_z`, so structure overlays directly between the two.

> **Python counterpart**: `analysis/pyscripts/plot_hits_at_hycal.py` — same pipeline via `prad2py`, no ROOT. Dumps a flat per-hit TSV/CSV (one row per HyCal cluster + per GEM hit projected to HyCal); plot externally with pandas/matplotlib. See [Python counterparts](#python-counterparts-pyscripts) below.

**Trigger filter**: only events with `trigger_bits == 0x100` (production physics trigger) contribute to the histograms.

**Multi-file** — chosen by the input path:
- `prad_023881.evio.*` → glob: enumerate all `prad_023881.evio.<digits>` siblings, fold into the same two histograms; gap warnings to stderr.
- A directory → same, run number sniffed from the directory name.
- `prad_023881.evio.00000` → single specific split.

Per event (after the trigger cut):
- HyCal: `WaveAnalyzer` → `mod->energize` → `HyCalCluster.FormClusters / ReconstructHits`. HC hits are built with `z = 0` (no shower-depth applied) so the lab transform places them at exactly z = `hycal_z`.
- GEM: `GemSystem.ProcessEvent` → `Reconstruct(GemCluster)` → per-detector hit list.
- Both go through the prad2det/prad2ana transforms: `RotateDetData` (per-detector tilt) → `TransformDetData` (per-detector position offset).
- GEM hits are then projected via `analysis::GetProjection(hits, hycal_z)` — straight line from target through each (x, y, z) intersected at z = `hycal_z`.
- All four GEM detectors fill a single combined left histogram.

Outputs the canvas in whatever format the extension implies (`.pdf`, `.png`, `.svg`, …) and a sibling `.root` file with both `TH2F`s and the canvas saved for re-plotting.

```bash
# full-run scan — glob discovers every split, warns on gaps:
.x ../analysis/scripts/plot_hits_at_hycal.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.*", \
    "hits_at_hycal.pdf")

# single split file (debugging):
.x ../analysis/scripts/plot_hits_at_hycal.C+( \
    "/data/stage6/prad_023867/prad_023867.evio.00000", \
    "hits_at_hycal_seg0.pdf")

# subset of events across the full run:
.x ../analysis/scripts/plot_hits_at_hycal.C+( \
    "/data/.../prad_023867.evio.*", "hits.png", 50000L)
```

Convenience overloads:
- `plot_hits_at_hycal(evio, out)`
- `plot_hits_at_hycal(evio, out, max_events)`
- `plot_hits_at_hycal(evio, out, max_events, run_num)`

Full 10-arg version (pass `""` to auto-discover any path):
`plot_hits_at_hycal(evio_path, out_path, max_events, run_num, gem_ped_file, gem_cm_file, hc_calib_file, daq_config, gem_map_file, hc_map_file)`.

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

Uses the project's decoder (`EvChannel`, `WaveAnalyzer`, `DaqConfig`) — no manual FADC parsing needed. Channel identity (HyCal module vs LMS ref) is read from `hycal_map.json`.

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
root -l rootlogon.C 'lms_alpha_normalize.C+("/path/to/data", 1234, "daq_config.json", "hycal_map.json")'
```

**Output:** `lms_alpha_run{N}.root` (per-module normalized LMS, reference time-series TGraphs) and a 6-panel summary PNG.

## Python counterparts (`pyscripts/`)

Two Python scripts mirror the ROOT analysis macros via the `prad2py` pybind11 module — same EVIO decode + HyCal/GEM reconstruction, but **no ROOT** and **flat TSV/CSV** output instead of ROOT trees / canvases. Run them anywhere `prad2py` is importable (`cmake -DBUILD_PYTHON=ON` and put the install dir on `PYTHONPATH`).

Both scripts share the same per-event pipeline (`_common.py` factors out runinfo loading, lab-frame transforms, file discovery, the trigger gate `0x100`, and the multi-file glob/dir/single mode).

### gem_hycal_matching.py

Same pipeline as `gem_hycal_matching.C`, same best-match rule (closest GEM hit per HC cluster × GEM detector pair within `N · σ_total`). Output is one row per matched tuple:

| Column | Notes |
|--------|-------|
| `event_num`, `trigger_bits` | event-level |
| `hc_idx`, `hc_x/y/z`, `hc_energy`, `hc_center`, `hc_nblocks`, `hc_sigma` | HyCal cluster (lab frame, z includes shower depth) |
| `det_id` (0..3) | which GEM won this row |
| `gem_x/y/z` | best-matched GEM hit, lab/target-centered mm |
| `gem_x_local`, `gem_y_local` | same hit in the GEM detector frame (no rotation/translation) |
| `gem_x_charge`, `gem_y_charge` | total ADC of the X / Y constituent cluster |
| `gem_x_peak`, `gem_y_peak` | max-strip ADC of the X / Y cluster |
| `gem_x_max_tb`, `gem_y_max_tb` | time sample (int) of the max-ADC strip — multiply by `ts_period` (default 25 ns) for ns |
| `gem_x_size`, `gem_y_size` | strip count of the X / Y cluster |
| `proj_x`, `proj_y`, `residual`, `sigma_total` | matching geometry |

```bash
# full run, glob discovers every split:
python analysis/pyscripts/gem_hycal_matching.py \
    /data/stage6/prad_023867/prad_023867.evio.* match_023867.tsv

# CSV, tighter cut, capped:
python analysis/pyscripts/gem_hycal_matching.py input.evio.* out.csv \
    --csv --match-nsigma 2.0 --max-events 50000
```

### plot_hits_at_hycal.py

Same pipeline as `plot_hits_at_hycal.C`. Dumps **all** hits (not just matched) — one row per HyCal cluster centroid and per GEM hit projected to z = `hycal_z`. Plot externally with pandas/matplotlib.

| Column | Notes |
|--------|-------|
| `event_num`, `trigger_bits` | event-level |
| `kind` | `"hycal"` or `"gem"` |
| `det_id` | 0..3 for GEM, -1 for HyCal |
| `x`, `y`, `z` | lab/target-centered mm at z = `hycal_z` |
| `energy` | MeV (HyCal); empty for GEM |

```bash
python analysis/pyscripts/plot_hits_at_hycal.py \
    /data/stage6/prad_023867/prad_023867.evio.* hits_023867.tsv

# minimal matplotlib:
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv("hits_023867.tsv", sep="\t")
for kind, sub in df.groupby("kind"):
    plt.hist2d(sub.x, sub.y, bins=260, range=[[-650, 650], [-650, 650]])
    plt.title(kind); plt.show()
```

Each script accepts the same path / overrides as its C++ counterpart (`--max-events`, `--run-num`, `--gem-ped-file`, etc.). Run with `--help` for the full list.

### plot_match_summary.py

Reads the per-match TSV/CSV from `gem_hycal_matching.py` (pandas) and emits four PNG plots (matplotlib) — no EVIO replay, just a fast post-processing pass on a table you already have.

| Output | What it shows |
|--------|---------------|
| `{prefix}_local_hits.png`  | 2×2 grid of 2D heatmaps — `(gem_x_local, gem_y_local)` per detector. Reveals acceptance / dead regions on each GEM in its own frame. |
| `{prefix}_lab_scatter.png` | Single scatter — `(gem_x, gem_y)` lab-frame, color-coded by `det_id`. All four GEMs overlaid. |
| `{prefix}_peak_adc.png`    | 2×2 histograms — `gem_x_peak` and `gem_y_peak` overlaid per detector. |
| `{prefix}_timing.png`      | 2×2 histograms — timing of the max-ADC strip per X / Y cluster, in ns (`gem_*_max_tb · ts_period`). |

Requires the post-2026-04 `gem_hycal_matching.py` columns (`gem_*_local`, `gem_*_peak`, `gem_*_max_tb`) — re-run the matcher if your TSV is older.

```bash
# default — saves four PNGs next to the input, then pops a GUI window:
python analysis/pyscripts/plot_match_summary.py match_023867.tsv

# headless / CI: explicit out-dir, finer binning, no GUI:
python analysis/pyscripts/plot_match_summary.py match_023867.tsv \
    --out-dir plots/ --bins 200 --no-show

# CSV input + non-default time-sample period (ns/sample):
python analysis/pyscripts/plot_match_summary.py match.csv --csv --ts-period 25.0
```

Flags: `--csv` (force CSV input), `--out-dir`, `--prefix` (default = input stem), `--bins` (default 120), `--ts-period` (default 25.0 ns), `--no-show`.

## Adding a Tool

Create `tools/my_tool.cpp`, then add to `CMakeLists.txt`:
```cmake
add_analysis_tool(my_tool tools/my_tool.cpp)
```

The helper takes care of the rest:
- compiles your `.cpp` into `prad2ana_my_tool` (binary prefix matches the install convention)
- links `libprad2ana.a` (transitively pulls in `prad2dec`, `prad2det`, ROOT)
- defines `DATABASE_DIR=...` so install-relative paths resolve
- routes the binary to `<build>/bin/`

If you add a *shared* source (something callable from multiple tools and from ACLiC scripts), put the implementation in `src/` and the declaration in `include/`, then list the new `.cpp` inside the `add_library(prad2ana STATIC ...)` call near the top of `CMakeLists.txt`. ACLiC scripts that link `libprad2ana.a` will pick up the new symbols automatically.

## Contributors
Yuan Li, Weizhi Xiong — Shandong University
Chao Peng - Argonne National Laboratory
