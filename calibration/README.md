# HyCal Calibration Tools

PRad-II, Jefferson Lab Hall B

## Operation Manual

### Quick Start

On **clonpc19**:

```bash
cd ~/prad2_daq/prad2evviewer
pip install -r calibration/requirements.txt   # first time only
./calibration/hycal_snake_scan.py --real
```

### Running a Scan

1. Select **Path** profile (`(autogen)` or a predefined path from `paths.json`).
2. For autogen: set **LG layers** (0 = PbWO4 only, 1--2 to include PbGlass).
3. Set **Start** module and **Count** (0 = scan all from start to end).
4. Set **Dwell time** and **Pos. threshold**.
5. Click **Start Scan**.

### Resume After Interruption

1. Find the last completed module (blue on map or in event log).
2. Select the next module as **Start**.
3. Click **Start Scan**.

### Shift Checklist

During the scan, verify that:

1. The **red crosshair** position matches the **occupancy plot** from the
   online event monitor -- the current module should have the highest occupancy.
2. The **FADC scalers** are consistent with the beam on the expected module.

If the occupancy does not match, **Stop** the scan, set **Start** to the
last good module, and try again. Log the observation and action on
**PRADLOG**. If the error persists, **contact the run coordinator**.

### Logging

Events are logged to `calibration/logs/snake_scan_YYYYMMDD_HHMMSS.log`.
**Upload log files when the scan or shift ends.**

### Troubleshooting

If a problem occurs, **try again first**. If it persists, contact the
run coordinator.

| Symptom | Likely cause |
|---------|--------------|
| PVs not connecting | IOC down, network / firewall |
| Motors don't move | Interlocks, motor enable, SPMG not Go |
| Move blocked | Target outside travel limits (see log) |
| Position errors | Motor speed, backlash, encoder |
| Move timeout (>300 s) | Motor stall, limit switch, IOC |

---

## Tools

| Script | Purpose |
|--------|---------|
| `hycal_snake_scan.py` | Snake scan GUI (simulation or real EPICS) |
| `scan_path_editor.py` | Manual path builder GUI |
| `scan_utils.py` | Shared types, constants, coordinate transforms |
| `paths.json` | Predefined scan path profiles |

### Snake Scan

```bash
python calibration/hycal_snake_scan.py              # simulation
python calibration/hycal_snake_scan.py --real        # real EPICS
python calibration/hycal_snake_scan.py --paths calibration/paths.json
```

Auto-generated path order: **Center (PbWO4)** -> **Bottom** -> nearest
**Side** -> **Top** -> other **Side**, minimising y-axis travel.

### Path Editor

```bash
python calibration/scan_path_editor.py
python calibration/scan_path_editor.py --paths calibration/paths.json
```

Click modules on the map to build a path. Intermediate modules on the
connecting line are auto-inserted. Profiles are saved to `paths.json`.

## Coordinate System

| Beam at HyCal (0,0) | ptrans_x = **-126.75** | ptrans_y = **10.11** |
|---|---|---|

```
ptrans_x = -126.75 + module_x       (x same direction)
ptrans_y =   10.11 - module_y       (y inverted)
```

Travel limits: ptrans_x **-582.65** to **329.15** mm,
ptrans_y **-672.50** to **692.72** mm.

**Safety:** Only four EPICS PVs are written: `ptrans_{x,y}.VAL` and
`ptrans_{x,y}.SPMG`. All others are read-only.
