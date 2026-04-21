# Scripts

Python utilities for HyCal / tagger visualisation and monitoring.

All GUI tools in this directory share `hycal_geoview.py` — a reusable HyCal
module map widget (rectangles in physical coordinates, optional colour bar,
optional zoom/pan, hover tooltips) and a shared theme palette. New GUIs
subclass `HyCalMapWidget` and override the paint hooks; see
`hycal_geoview.py` for the extension points.

Theme: the module exposes `set_theme("dark" | "light")` and an
`apply_theme_palette(window)` helper. Scripts typically wire a `--theme`
argument in `main()` before constructing any widget.

Installed wrappers are generated when `-DBUILD_PYTHON=ON` is used at build
time — each GUI below can then be launched from `$PATH` by its bare name
(the `python scripts/...` examples shown work equally well from the source
checkout).

## HyCal GUIs

### hycal_scaler_map.py

PyQt6 live colour-coded HyCal FADC scaler map. Polls
`B_DET_HYCAL_FADC_<name>` EPICS channels every few seconds. Requires
`pyepics` for real EPICS; `--sim` works without it.

```bash
hycal_scaler_map          # real EPICS (default)
hycal_scaler_map --sim    # simulation (random values)
```

### hycal_pedestal_monitor.py

PyQt6 GUI for measuring and monitoring FADC250 pedestals on all 7 HyCal
crates. Reads pedestal means from `.cnf` files, parses per-channel RMS
from `faV3peds` stdout, and flags irregular channels.

```bash
hycal_pedestal_monitor          # view existing data
hycal_pedestal_monitor --sim    # test with simulated data
```

**Shift pedestal check (operator procedure).** Pedestals must be measured
before the first DAQ run of each shift while DAQ is idle.

1. Make sure the DAQ is **stopped**.
2. Launch `hycal_pedestal_monitor`.
3. Click **Measure Pedestals** and confirm. The tool SSHs to
   `adchycal1`..`adchycal7` and runs `faV3peds` (takes a few minutes).
4. Inspect the two maps (left: current mean, right: difference from
   configured) and the report panel for flagged channels:
   - `DEAD` — avg < 1, rms < 0.1
   - `OUT OF RANGE` — mean outside 50..300
   - `HIGH RMS` — sigma > 1.5
   - `DRIFT` — shifted > 3 counts from configured
5. Click **Save Report** to save for the shift log.
6. If new issues appear, notify the run coordinator before starting data
   taking. Thresholds are defined at the top of the script.

### hycal_gain_monitor.py

PyQt6 viewer for per-module LMS-based gain factors. Loads the text
`prad_{run:06d}_LMS.dat` files produced by the offline gain analysis,
plots the HyCal geo map, LMS reference-channel stability, a run-to-run
drift view, and a table of irregular `(module, run)` outliers.

```bash
hycal_gain_monitor
```

### hycal_map_builder.py

Generic HyCal geometry viewer that colour-maps user data loaded from JSON
or plain-text files. Supports a day/night theme toggle, PbGlass alpha
slider, zoom/pan, and palette cycling (click the colour bar).

- **JSON**: `{"<module_name>": {"<field>": <value>, ...}, ...}` — the last
  entry of a history list is used; nested dicts are flattened with dot
  notation; non-numeric fields (timestamps) are ignored.
- **Text**: whitespace/comma/tab-delimited rows `<module> <v1> <v2> ...`
  with optional header.

```bash
hycal_map_builder                     # empty map
hycal_map_builder mydata.json         # auto-load
hycal_map_builder mydata.txt --field rms
```

### waveform_viewer.py

Event-by-event browser for FADC250 waveforms. Opens an EVIO file in
**random-access** mode (native event-pointer table via
`evc::EvChannel::OpenRandomAccess`), does a single Scan-only pass to
index physics sub-events, and then lets you step through events with
Prev / Next / Jump — all fast in both directions.

The layout is a horizontal split: HyCal geo picker + raw-waveform plot on
the left, four stacked histograms on the right (peak height, integral,
time, peaks/event). For each viewed event the selected module's waveform
is drawn and its four histograms accumulate; clicking a module in the
geo picker (only lit modules are clickable) switches the selection.

A **Process next 10k** button runs a background pass that fills the
selected module's histograms without displaying waveforms, for fast
accumulation.

```bash
waveform_viewer                        # File → Open…
waveform_viewer run.evio.00000
```

File → Save writes the current per-module histograms to JSON for later
inspection.

## Tagger / TDC

### tagger_viewer.py

PyQt6 viewer for V1190 TDC hits from the tagger crate (ROC `0x008E`,
bank `0xE107`). Shows a per-slot bar chart of hits/channel, a
single-channel TDC-value histogram, plus event-wise correlation tabs
(Δt = A − B and a 2-D tdc(A) vs tdc(B) heatmap). The bar chart
auto-sizes its x-axis to 16 / 32 / 64 / 128 channels based on the
highest channel actually hit. Human-readable counter names come from
`database/tagger_map.json`.

**Offline — evio file** (decoded in-process by `prad2py`):

```bash
tagger_viewer /data/stage6/prad_023667/prad_023667.evio.00000 \
              -n 200000          # limit number of physics events (optional)
              --roc 0x8E         # restrict to the tagger ROC (optional)
```

**Online — live ET stream** from a running `prad2_server`. The server
only decodes tagger TDC hits when at least one client is subscribed, so
regular monitoring is unaffected:

```bash
# DAQ machine (one-time)
prad2_server --online --port 5051

# Viewer (anywhere with PyQt6 + QtWebSockets installed)
tagger_viewer --live ws://clondaq6:5051
```

On startup with `--live`, a fast subscribe/ack round-trip is done
*before* the main window opens — if the server is unreachable or the
protocol doesn't match, `tagger_viewer` exits with a clear error rather
than showing an empty window. Pass `--no-smoke-test` to skip it.

The File menu also has *Connect to prad2_server…* (Ctrl+L) and
*Disconnect*. Pause / Clear buttons sit next to the Bins spinner.
Memory is capped at 10 M hits (rolling — oldest half is dropped).

Wire protocol — WebSocket JSON messages `tagger_subscribe` /
`tagger_subscribed` / `tagger_unsubscribe`.  Binary frame format
(useful for writing a different client):

```
magic        "TGR1"   (4 bytes)
flags        u32      (bit 0 = some frames have been dropped)
n_hits       u32
first_seq    u32
last_seq     u32
dropped      u32      (cumulative since server start)
records      n_hits × 16-byte packed BinHit
               u32 event_num, u32 trigger_bits, u16 roc_tag,
               u8 slot, u8 channel_edge (bit 7 = edge, bits 6:0 = channel),
               u32 tdc
```

### extract_tagger_map.py

One-shot helper that parses a tagger cabling source and emits
`database/tagger_map.json` — the counter-name lookup used by
`tagger_viewer`. See the script header for the input format.

## DAQ dev tools (`scripts/daq_tool/`)

Dev-only GUIs that don't ship with the installation (run them from the
source checkout). The CMake install step excludes this directory.

### trigger_mask_editor.py

PyQt6 visual editor for FAV3 trigger masks. Displays a HyCal geo view
(with LMS1-3, LMSP, V1-V4 below) and lets you click or drag modules to
toggle channels off/on. Generates trigger mask `.cnf` files — only slots
with disabled channels are written. Unmapped DAQ channels (slot
positions with no module) are always masked off.

```bash
python scripts/daq_tool/trigger_mask_editor.py                     # start fresh
python scripts/daq_tool/trigger_mask_editor.py -i existing.cnf     # load existing mask
python scripts/daq_tool/trigger_mask_editor.py -o output.cnf       # set default save path
```

### fadc_gain_config.py

Generates a text-based `adchycal_gain.cnf` for the FADC250 DAQ
(`FAV3_ALLCH_GAIN` entries, one per 16-channel slot, grouped by
crate/slot). Gains come from a calibration JSON (`-c path.json`,
`{"name", "factor"}` entries) or uniform values per module type
(`--pbwo4-gain` / `--pbglass-gain`).

```bash
python scripts/daq_tool/fadc_gain_config.py
python scripts/daq_tool/fadc_gain_config.py -c database/calibration/adc_to_mev_factors_cosmic.json
python scripts/daq_tool/fadc_gain_config.py --pbwo4-gain 0.15 --pbglass-gain 0.12
```

## Operator shell scripts (`scripts/shell/`)

Not installed — run directly from the source checkout on DAQ machines.

- `run_gain_monitor.sh` — wrapper that parallelises `prad2ana_gain_monitor`
  across sub-files of a run and merges the outputs via `hadd`. Requires
  the installed analysis binaries on `$PATH` (source `prad2_setup.sh`
  first).

  ```bash
  scripts/shell/run_gain_monitor.sh <run_number> <num_cpus> [subfile_min] [subfile_max]
  ```

- `prad2_server_tmux_template` — tmux session template for running
  `prad2_server` alongside `prad2_client` on a split window.

## Using prad2py directly

`prad2py` exposes the decoder through a `dec` submodule and the
reconstruction through `det` — useful for custom offline analysis that
goes beyond what the viewers above do. Build it once with
`-DBUILD_PYTHON=ON`, then:

```python
import prad2py
from prad2py import dec                 # evio reader + event types

cfg = dec.load_daq_config()             # installed daq_config.json
ch  = dec.EvChannel(); ch.set_config(cfg)
ch.open("/data/.../prad_023671.evio.00000")

while ch.read() == dec.Status.success:
    if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
        continue
    for i in range(ch.get_n_events()):
        ch.select_event(i)                  # picks sub-event + clears cache
        info = ch.info()                    # cheapest: TI/trigger metadata
        fadc_evt = ch.fadc()                # FADC250 waveforms, cached
        tdc_evt  = ch.tdc()                 # V1190 tagger hits, cached
        # gem_evt = ch.gem();  vtp_evt = ch.vtp()
```

Random-access mode (also used internally by `waveform_viewer` and
`prad2_server`'s file mode):

```python
ch = dec.EvChannel(); ch.set_config(cfg)
ch.open_random_access("run.evio.00000")
n = ch.get_random_access_event_count()
# jump to any evio event in O(1)
for i in (0, n // 2, n - 1):
    ch.read_event_by_index(i)
    ch.scan(); ch.select_event(0)
    print(int(ch.info().event_number))
```

Full reconstruction helpers live in `det`:

```python
from prad2py import dec, det

# --- decoder -----------------------------------------------------------
cfg = dec.load_daq_config()
ch  = dec.EvChannel(); ch.set_config(cfg); ch.open("run.evio.00000")

# --- GEM reconstruction -----------------------------------------------
gsys = det.GemSystem()
gsys.init("database/gem_map.json")
gsys.load_pedestals("database/gem_ped.json")    # optional
gcl  = det.GemCluster()

# --- HyCal reconstruction ---------------------------------------------
hsys = det.HyCalSystem()
hsys.init("database/hycal_modules.json", "database/daq_map.json")
hsys.load_calibration("database/hycal_calib.json")
hcl  = det.HyCalCluster(hsys)

while ch.read() == dec.Status.success:
    if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
        continue
    for i in range(ch.get_n_events()):
        ch.select_event(i)

        # GEM 2-D hits
        gsys.clear()
        gsys.process_event(ch.gem())
        gsys.reconstruct(gcl)
        for h in gsys.get_all_hits():
            print("GEM", h.det_id, h.x, h.y, h.x_charge, h.y_charge)

        # HyCal clusters — feed per-module energies yourself (e.g. from
        # ch.fadc() + your calibration), then cluster:
        # hcl.clear()
        # for module_idx, energy_mev in my_hycal_hits(ch.fadc()):
        #     hcl.add_hit(module_idx, energy_mev)
        # hcl.form_clusters()
        # for c in hcl.reconstruct_hits():
        #     print("HyCal", c.center_id, c.x, c.y, c.energy)
```

The helper `prad2py.load_tdc_hits(path, ...)` remains for the common
"one-shot flat table of hits" workflow.

## GEM tools

GEM-specific scripts and the `gem_dump` C++ binary live in the top-level
[`gem/`](../gem/README.md) directory — see that README for details and
GEM detector reference notes.

## Tagger ↔ HyCal coincidence (ROOT macro)

See `analysis/scripts/tagger_hycal_correlation.C` — a self-contained
ROOT/ACLiC macro that builds ΔT histograms for (T10R, E49…E53) pairs,
Gaussian-fits each coincidence peak, applies a ±Nσ timing cut, and plots
the W1156 peak height/integral for the selected events.

```bash
cd build
root -l ../analysis/scripts/rootlogon.C
.x ../analysis/scripts/tagger_hycal_correlation.C+( \
     "/data/stage6/prad_023671/prad_023671.evio.00000", \
     "tagger_w1156_corr.root", 500000)
```
