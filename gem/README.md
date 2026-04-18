# GEM

Tools, scripts, and reference notes for the PRad-II GEM tracker.

All GEM-specific code lives here:

| File | Purpose |
|---|---|
| `gem_dump.cpp` | C++ CLI: raw/hits/clusters/evdump/summary/ped modes |
| `gem_event_viewer.py` | PyQt6 event-by-event GUI (evio file → live reconstruction) |
| `gem_cluster_view.py` | Static plotter for `gem_dump -m evdump` JSON output |
| `gem_layout.py` | Visualize strip geometry from `gem_map.json` |
| `gem_strip_map.py` | Thin wrapper over `prad2py.det.map_strip` (library) |
| `gem_view.py` | Matplotlib rendering + GemSystem adapters (library) |
| `check_strip_map.py` | Dev: cross-validate pipeline vs PRadAnalyzer and mpd_gem_view_ssp |
| `CMakeLists.txt` | Builds `gem_dump`, installs binary + Python scripts |

## Detector facts (as of 2026-04-18)

### Geometry
- **4 identical GEMs** — GEM0 + GEM1 (upstream plane), GEM2 + GEM3 (downstream).
  Each plane = two half-sensors overlapping at the beam line; one rotated 180°
  in XY so the beam holes align.
- **Active area:** 12 X-APVs × 128 ch @ 0.4 mm pitch (X), 24 Y-APVs × 128 ch @ 0.4 mm (Y).
- **Beam hole:** 52 × 52 mm rectangle on the beam side, centered in Y (exact
  location derived from the "match" APV strip positions at runtime — see
  `gem_map.json`'s `hole` block).

### Readout
- **2 crates** — `0x31` (`gemroc1`) and `0x34` (`gemroc2`), 72 APVs each
  → **144 APVs total** split evenly across the four detectors (36 APVs per GEM).
- **52 MPDs / event** (observed in run 001137). Each MPD hosts up to 16 APVs
  via individual fibers.
- **APV25 front-end:** 128 analog channels, multiplexed to a single ADC.
  SSP firmware samples 6 time bins per trigger (13-bit signed ADC).
- **Bank tag:** `0x0DE9` (Hall-A standard MPD raw format; confirmed in run
  001137). Nested under the ROC bank (`0x31` / `0x34`) inside the physics
  event bank (`0x82`).

### Data modes
| Mode | `ApvData.nstrips` | Typical use | Offline pipeline |
|---|---|---|---|
| **Full readout** | `== 128` (all channels sent) | Pedestal calibration runs, debug modes with CM passthrough | `processApv` runs full offline chain (pedestal subtract → sorting common-mode → `noise × zero_sup_threshold` ZS) |
| **Online ZS** | `< 128` (firmware dropped some) | Production runs | `processApv` short-circuits: every strip present in the bank is a surviving hit (firmware already did pedestal + CM + ZS); no ped file needed |

Auto-detected per APV. We use the **`nstrips` count**, not
`has_online_cm` — the MPD can emit its type-`0xD` CM debug headers
(which set `has_online_cm = true`) while still sending all 128 strips
raw. `nstrips < APV_STRIP_SIZE` is the only signal that reliably says
"firmware dropped some channels, so it also pedestal-subtracted."
See `prad2det/src/GemSystem.cpp`.

### Strip mapping (6-step pipeline)
Shared between online reconstruction (`gem::MapStrip` in
`prad2det/src/GemSystem.cpp`) and these scripts (`gem_strip_map.py`
wraps the same C++ entry via `prad2py.det`). Per-APV config in
`gem_map.json`:

| Field | Default | Meaning |
|---|---|---|
| `pos` | — | APV position on its plane (0 .. n_apvs−1) |
| `orient` | 0 | `1` flips the APV (`strip → 127 − strip`) |
| `pin_rotate` | 0 | Rotated connector pins — `16` for pos 11 near the beam hole |
| `shared_pos` | −1 | "Share this other APV's plane slot" (−1 = use `pos`) |
| `hybrid_board` | `true` | MPD hybrid-board pin conversion (`false` for SRS) |
| `match` | `""` | `"+Y"` / `"-Y"` marks split APVs above/below the beam hole |

Full geometry details + why certain parameters can't be collapsed live
in `memory/project_gem_geometry.md`.

## Typical workflows

### Open an event interactively

```bash
gem_event_viewer /volatile/hallb/prad/<run>/<file>.evio.00000
```

GUI pre-scans the file (progress bar), then lets you step through events
with Prev/Next/Goto. Threshold sliders re-run reconstruction on cached SSP
data — no disk I/O per slider change. Advanced tuning dock exposes every
`GemSystem` / `GemCluster` knob.

### Dump + visualize single events

```bash
# Pick events of interest (default: first 1 event with 2-D hits)
gem_dump -m evdump run.evio.00000 -n 10 -f clusters=2:3 -o /tmp/evt.json

# Render them to PNG
gem_cluster_view /tmp/evt_1.json
```

`-f` is a boolean filter — `clusters=2:3` = "≥2 clusters in ≥3 detectors".
See `gem_dump --help` for the full grammar.

### Pedestal calibration (only needed for full-readout test data)

```bash
# 1. Compute per-strip pedestals from a full-readout run.
gem_dump -m ped /volatile/.../gem0gem1_001137.evio.00000 \
         -o /volatile/.../gem0gem1_001137/gem_ped.json

# 2. Run any downstream analysis with those pedestals loaded.
gem_dump -m clusters /volatile/.../gem0gem1_001137.evio.00000 \
         -P /volatile/.../gem0gem1_001137/gem_ped.json -n 20

# gem_event_viewer takes the same file via --gem-ped:
gem_event_viewer /volatile/.../gem0gem1_001137.evio.00000 \
                 --gem-ped /volatile/.../gem0gem1_001137/gem_ped.json
```

Production runs with online ZS don't need step 1 — `processApv` skips the
offline pedestal/CM chain automatically.

**No auto-discovery.** Pedestals are per-run calibration products and a
wrong file is worse than none.  `gem_dump` and `gem_event_viewer` both
require the pedestal path to be passed explicitly (`-P` / `--gem-ped`).
If you run them against full-readout data without one, they'll warn
loudly on stderr / a modal dialog — not silently reconstruct empty
events.

### Summary diagnostics

```bash
gem_dump -m summary <file>.evio.00000 -n 100
```

Prints per-event MPD / APV / strip / hit / cluster / 2D-hit counts.
Useful first step after opening a new file to confirm the data makes
sense (non-zero strips → firmware is sending data; non-zero hits →
pedestals are working / online-ZS is active; non-zero 2D hits → XY
matching is succeeding).

### Visualize strip layout

```bash
gem_layout                    # uses database/gem_map.json
gem_layout path/to/alt_map.json
```

Draws every strip of one detector (all 4 are identical), overlays APV
boundaries and the beam hole, writes `gem_layout.png`.

### Dev sanity check

```bash
python gem/check_strip_map.py
```

For every APV in `gem_map.json`, maps all 128 channels through both
mpd_gem_view_ssp's and PRadAnalyzer's reference implementations and our
own, asserts they match. Run this after any change to the strip-mapping
pipeline.

## References

- `prad2det/src/GemSystem.cpp` — pedestal/CM/ZS + strip mapping implementation.
- `prad2det/src/GemCluster.cpp` — strip clustering + XY matching.
- `prad2dec/src/SspDecoder.cpp` — SSP/MPD/APV bitfield decoder.
- `database/gem_map.json` — detector geometry + per-APV mapping (source of truth).
- `docs/rols/banktags.md` — bank-tag reference including MPD `0x0DE9`.
- `memory/project_gem_geometry.md` — deep-dive on the 6-step mapping pipeline.
