# HyCal Calibration Tools

PRad-II, Jefferson Lab Hall B

---

## Gain Equalizer — Operator Manual

### Quick Start

On a counting house machine, logged in as **clasrun**:

```bash
cd /home/clasrun/prad2_daq/prad2evviewer
python calibration/hycal_gain_equalizer.py --expert
```

The default settings (server, HV address, password, target ADC, min counts, etc.)
should work out of the box. **Discuss with the Run Coordinator (RC) before changing
the settings.**

### Running a Gain Scan

1. Select **Path**: `snake-all-pwo-r`.
2. Set **Start** module. The plan is to finish the bottom part first, then
   return to the top. Ask the RC if you do not know which module to start from.
3. Click **Start**.

### During the Scan

Pay attention to:

- **Red marker** (expected beam position) vs. **hot region** on the scaler map —
  they should be at the same location.
- **Histogram** in the right panel — it should show a clear spectrum building up.
- **Event log** — watch for `WARN` and `ERROR` messages.

### On ERROR (Scan Stops)

An `ERROR` stops the scan automatically.

1. Note the module name and error message from the log.
2. The **Start** combo updates to show the failed module — verify it is correct.
3. Click **Start** to retry from that module.
4. If the same error repeats, **contact the RC**.

### After Each Row

1. Find the screenshot files in
   `/home/clasrun/prad2_daq/prad2evviewer/calibration/logs/`:
   ```
   GE_20260406_143025_W100_success.png
   GE_20260406_143530_W101_failure.png
   ```
   Files are named `GE_{timestamp}_{module}_{status}.png` and sort by time.
2. Post a log entry to **PRADLOG** with all screenshots from the completed row.
3. Failed modules need a redo or manual gain equalization — discuss with the RC
   if you are not sure how to proceed.

### Event Log Files

Full event logs are saved to `calibration/logs/gain_eq_YYYYMMDD_HHMMSS.log`.
Upload these with the PRADLOG entry.

---

## Snake Scan — Operator Manual

### Quick Start

```bash
cd /home/clasrun/prad2_daq/prad2evviewer
python calibration/hycal_snake_scan.py --expert
```

### Running a Scan

1. Select **Path** profile (`(autogen)` or a predefined path from `paths.json`).
2. For autogen: set **LG layers** (0 = PbWO4 only, 1--2 to include PbGlass).
3. Talk to RC if you're not sure about the scan path.
4. Set **Start** module and **Count** (0 = scan all from start to end).
5. Set **Dwell time** and **Pos. threshold**.
6. Click **Start Scan**.

### During the Scan

Verify that:

1. The **red marker** position matches the **scaler hot spot** — the current
   module should have the highest scaler reading.
2. **WARN/ERROR** messages in the event log are addressed promptly.

### Resume After Interruption

1. Find the last completed module (the **Done** colour on the map, or check
   the event log).
2. Select the next module as **Start**.
3. Click **Start Scan**.

---

## Tools

| Script | Purpose |
|--------|---------|
| `hycal_gain_equalizer.py` | Automatic gain equalization (expert/simulation/observer) |
| `hycal_snake_scan.py` | Snake scan with dwell time (expert/simulation/observer) |
| `scan_path_editor.py` | Manual path builder GUI |
| `gain_scanner.py` | Gain scan engine, spectrum analyzer, server/HV clients |
| `scan_geoview.py` | HyCal map widget with scaler overlay |
| `scan_epics.py` | EPICS PV utilities (motor, scaler) |
| `scan_engine.py` | Scan engine, path builder |
| `scan_utils.py` | Shared types, constants, coordinate transforms, theme |
| `paths.json` | Predefined scan path profiles |

### Command-Line Modes

```bash
python calibration/hycal_gain_equalizer.py             # simulation (read-only)
python calibration/hycal_gain_equalizer.py --expert    # expert (full control)
python calibration/hycal_gain_equalizer.py --observer  # observer (read-only)

python calibration/hycal_snake_scan.py                 # simulation
python calibration/hycal_snake_scan.py --expert        # expert
python calibration/hycal_snake_scan.py --observer      # observer

python calibration/scan_path_editor.py                 # path editor
```

## Coordinate System

| Beam at HyCal (0,0) | ptrans_x = **-126.75** | ptrans_y = **10.11** |
|---|---|---|

```
ptrans_x = -126.75 + module_x       (x same direction)
ptrans_y =   10.11 - module_y       (y inverted)
```

Travel limits: ptrans_x **-582.65** to **329.15** mm,
ptrans_y **-672.50** to **692.72** mm.

## Safety

- **Expert mode** writes to four motor PVs only:
  `ptrans_{x,y}.VAL` and `ptrans_{x,y}.SPMG`.
- **Gain equalizer** additionally sends HTTP commands to `prad2hvd` (HV set)
  and `prad2_server` (histogram clear). HV limits are enforced server-side
  by `prad2hvd`.
- **Observer mode** uses real EPICS reads but blocks all writes.
- **Simulation mode** uses a fake motor (no real EPICS writes) and blocks all
  HV/server writes.
