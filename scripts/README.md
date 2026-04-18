# Scripts

Python utilities for detector visualization and monitoring.

## hycal_scaler_map.py

PyQt6 live colour-coded HyCal FADC scaler map. Polls `B_DET_HYCAL_FADC_<name>` EPICS channels every 10 s. Requires `pyepics` for real EPICS; `--sim` works without it.

```bash
python scripts/hycal_scaler_map.py          # real EPICS (default)
python scripts/hycal_scaler_map.py --sim    # simulation (random values)
```

## hycal_pedestal_monitor.py

PyQt6 GUI for measuring and monitoring FADC250 pedestals on all 7 HyCal crates. Reads pedestal means from `.cnf` files, parses per-channel RMS from `faV3peds` stdout, and flags irregular channels. No numpy/matplotlib -- only PyQt6.

```bash
python scripts/hycal_pedestal_monitor.py          # view existing data
python scripts/hycal_pedestal_monitor.py --sim     # test with simulated data
```

### Shift pedestal check (operator procedure)

Pedestals must be measured **before the first DAQ run of each shift** while DAQ is idle.

1. Make sure the DAQ is **stopped**.
2. Launch: `python scripts/hycal_pedestal_monitor.py`
3. Click **Measure Pedestals** and confirm. The tool SSHs to `adchycal1`--`adchycal7` and runs `faV3peds` (takes a few minutes).
4. Inspect the two maps (left: current mean, right: difference from configured) and the report panel for flagged channels:
   - `DEAD` -- avg < 1, rms < 0.1
   - `OUT OF RANGE` -- mean outside 50--300
   - `HIGH RMS` -- sigma > 1.5
   - `DRIFT` -- shifted > 3 counts from configured
5. Click **Save Report** to save for the shift log.
6. If new issues appear, notify the run coordinator before starting data taking.

Thresholds are defined at the top of the script and can be adjusted.

## trigger_mask_editor.py

PyQt6 visual editor for FAV3 trigger masks. Displays a HyCal geo view (with LMS1-3, LMSP, V1-V4 below) and lets you click or drag modules to toggle channels off/on. Generates trigger mask `.cnf` files -- only slots with disabled channels are written. Unmapped DAQ channels (slot positions with no module) are always masked off.

```bash
python scripts/trigger_mask_editor.py                     # start fresh
python scripts/trigger_mask_editor.py -i existing.cnf     # load existing mask
python scripts/trigger_mask_editor.py -o output.cnf       # set default save path
```

## gem_layout.py

Visualize GEM strip layout from `gem_map.json`.

```bash
python scripts/gem_layout.py [gem_map.json]
```

## gem_cluster_view.py

Visualize GEM clustering from `gem_dump -m evdump` JSON output.

```bash
python scripts/gem_cluster_view.py <event.json> [gem_map.json] [--det N] [-o file.png]
```

## tagger_viewer.py

PyQt6 viewer for V1190 TDC hits from the tagger crate (ROC `0x008E`, bank
`0xE107`). Shows a per-slot bar chart of hits/channel, a single-channel
TDC-value histogram, plus event-wise correlation tabs (Δt = A − B and a
2-D tdc(A) vs tdc(B) heatmap). No matplotlib / pyqtgraph; plots are
drawn with QPainter (numpy is the only scientific dep). The bar chart
auto-sizes its x-axis to 16 / 32 / 64 / 128 channels based on the
highest channel actually hit. Human-readable counter names come from
`database/tagger_map.json`.

Two data sources:

**Offline — evio file** (decoded in-process by `prad2py`):

```bash
cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build
# optional — the viewer auto-adds build/python/ to sys.path
export PYTHONPATH="$PWD/build/python:$PYTHONPATH"
python scripts/tagger_viewer.py /data/stage6/prad_023667/prad_023667.evio.00000 \
       -n 200000          # limit number of physics events (optional)
       --roc 0x8E         # restrict to the tagger ROC (optional)
```

**Online — live ET stream** from a running `prad2_server`. The server
only decodes tagger TDC hits when at least one client is subscribed, so
regular monitoring is unaffected:

```bash
# DAQ machine (one-time)
./build/bin/prad2_server --online --port 5051

# Viewer (anywhere with PyQt6 + QtWebSockets installed)
python scripts/tagger_viewer.py --live ws://clondaq6:5051
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
(useful for anyone writing a different client):

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

## Using prad2py directly (Phase 1 bindings)

`prad2py` exposes the decoder through a `dec` submodule — useful for custom
offline analysis that goes beyond what the tagger_viewer does. Build it with
`-DBUILD_PYTHON=ON` once, then:

```python
import prad2py
from prad2py import dec                         # evio reader + event types

cfg = dec.load_daq_config()                     # installed daq_config.json
ch  = dec.EvChannel()
ch.set_config(cfg)
st = ch.open("/data/.../prad_023671.evio.00000")
assert st == dec.Status.success

while ch.read() == dec.Status.success:
    if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
        continue
    for i in range(ch.get_n_events()):
        ch.select_event(i)                              # picks sub-event + clears cache

        # Cheapest accessor — TI/trigger metadata only (no detector decode):
        info = ch.info()
        # …do something with info.event_number / .trigger_bits / .timestamp …

        # Decode detector products on demand — each is cached for repeat calls
        # on the same sub-event, so asking for fadc() then gem() costs nothing
        # extra if you only need one of them:
        fadc_evt = ch.fadc()                            # FADC250/ADC1881M waveforms
        tdc_evt  = ch.tdc()                             # V1190 tagger hits
        # gem_evt = ch.gem();  vtp_evt = ch.vtp()       # when needed

        for roc_idx in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(roc_idx)
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    samples = slot.channel(c).samples   # numpy uint16 array
                    …
        for j in range(tdc_evt.n_hits):
            h = tdc_evt.hit(j)                          # TdcHit
            …
```

The helper `prad2py.load_tdc_hits(path, ...)` is still available for the
common "one-shot flat table of hits" workflow and lives on top of the
`dec` submodule. Phase 2 will add `prad2py.det` (HyCal / GEM).

### Tagger ↔ HyCal coincidence

See `analysis/scripts/tagger_hycal_correlation.C` — a self-contained
ROOT/ACLiC macro that builds ΔT histograms for (T10R, E49…E53) pairs,
Gaussian-fits each coincidence peak, applies a ±Nσ timing cut, and
plots the W1156 peak height/integral for the selected events.

```bash
cd build
root -l ../analysis/scripts/rootlogon.C
.x ../analysis/scripts/tagger_hycal_correlation.C+( \
     "/data/stage6/prad_023671/prad_023671.evio.00000", \
     "tagger_w1156_corr.root", 500000)
```
