#!/usr/bin/env python3
"""
_common.py — shared helpers for analysis/pyscripts/.

Both gem_hycal_matching.py and plot_hits_at_hycal.py do the same setup
(load DAQ config, runinfo, HyCal+GEM systems, discover EVIO splits) and
the same per-event boilerplate (waveform → cluster, GEM ProcessEvent +
Reconstruct, lab-frame transform).  This module factors that out so the
two scripts only differ in their per-event accumulation + output.

Mirrors the C++ helpers in:
  analysis/scripts/script_helpers.h     (path / file discovery)
  analysis/include/ConfigSetup.h        (RotateDetData + TransformDetData)

Requires:
  prad2py    (built from python/, exposes dec.* + det.*)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

try:
    from prad2py import dec, det
except ImportError as e:
    raise SystemExit(
        f"[ERROR] cannot import prad2py: {e}\n"
        "        Build the python bindings (cmake -DBUILD_PYTHON=ON) and "
        "ensure the install directory is on PYTHONPATH."
    )


# ============================================================================
# Path / run-number helpers (mirror script_helpers.h)
# ============================================================================

_RUN_PAT = re.compile(r"(?:prad|run)_0*(\d+)", re.IGNORECASE)


def extract_run_number(path: str) -> int:
    """Sniff the run number out of 'prad_NNNNNN.evio.*'-style names. -1 if none."""
    if not path:
        return -1
    m = _RUN_PAT.search(path)
    if not m:
        return -1
    try:
        return int(m.group(1))
    except ValueError:
        return -1


def resolve_db_path(p: str) -> str:
    """Resolve a possibly-relative database path against PRAD2_DATABASE_DIR."""
    if not p:
        return p
    if os.path.isabs(p):
        return p
    db = os.environ.get("PRAD2_DATABASE_DIR")
    if db is None:
        return p
    return os.path.join(db, p)


def discover_runinfo_path() -> Optional[str]:
    """Read database/reconstruction_config.json (under PRAD2_DATABASE_DIR or
    ./database) and return the resolved runinfo path; None if missing/malformed."""
    db = os.environ.get("PRAD2_DATABASE_DIR", "database")
    cfg_path = Path(db) / "reconstruction_config.json"
    try:
        with open(cfg_path, "r", encoding="utf-8") as f:
            j = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
    ri = j.get("runinfo")
    if not isinstance(ri, str):
        return None
    return resolve_db_path(ri)


def load_matching_config() -> tuple[tuple[float, float, float], list[float]]:
    """Read the 'matching' section from
    database/reconstruction_config.json and return ((A, B, C), gem_pos_res).
    Missing keys / file fall back to (2.6, 0, 0) and [0.1]*4 — the legacy
    inline values, so existing analysis behavior is preserved."""
    A, B, C = 2.6, 0.0, 0.0
    gem = [0.1, 0.1, 0.1, 0.1]
    db = os.environ.get("PRAD2_DATABASE_DIR", "database")
    try:
        with open(Path(db) / "reconstruction_config.json", "r", encoding="utf-8") as f:
            j = json.load(f)
    except (OSError, json.JSONDecodeError):
        return (A, B, C), gem
    m = j.get("matching")
    if not isinstance(m, dict):
        return (A, B, C), gem
    h = m.get("hycal_pos_res")
    if isinstance(h, list) and len(h) >= 3:
        A, B, C = float(h[0]), float(h[1]), float(h[2])
    g = m.get("gem_pos_res")
    if isinstance(g, list) and g:
        gem = [float(v) for v in g]
    return (A, B, C), gem


def hycal_pos_resolution(A: float, B: float, C: float, energy_mev: float) -> float:
    """sigma(E) at the HyCal face (mm), mirroring HyCalSystem::PositionResolution."""
    import math
    E_GeV = energy_mev / 1000.0 if energy_mev > 0 else 1e-6
    a = A / math.sqrt(E_GeV)
    b = B / E_GeV
    return math.sqrt(a * a + b * b + C * C)


def discover_split_files(any_path: str) -> list[str]:
    """Three modes by input shape (mirrors discover_split_files in
    script_helpers.h):
      * '*' in path  → glob mode: enumerate every sibling
        prad_<run>.evio.<digits>, sort by suffix, warn (stderr) on
        gaps from .00000 to highest.
      * directory    → same enumeration, sniff run from dir name.
      * anything else → return [any_path] unchanged (single-file mode)."""
    if not any_path:
        return []
    p = Path(any_path)
    wants_glob = "*" in any_path
    is_dir = p.is_dir()

    if not wants_glob and not is_dir:
        return [any_path]

    if is_dir:
        directory = p
        run = extract_run_number(p.name)
    else:
        directory = p.parent if str(p.parent) else Path(".")
        run = extract_run_number(p.name)
        if run < 0:
            run = extract_run_number(directory.name)

    if run < 0 or not directory.is_dir():
        sys.stderr.write(
            f"[WARN] discover_split_files: cannot resolve run/dir from "
            f"{any_path!r} — passing through as a single file.\n"
        )
        return [any_path]

    pat = re.compile(rf"^prad_0*{run}\.evio\.(\d+)$", re.IGNORECASE)
    matched: list[tuple[int, str]] = []
    for entry in directory.iterdir():
        m = pat.match(entry.name)
        if m:
            try:
                matched.append((int(m.group(1)), str(entry)))
            except ValueError:
                pass
    matched.sort()

    if matched:
        last = matched[-1][0]
        seen = {idx for idx, _ in matched}
        missing = [i for i in range(0, last + 1) if i not in seen]
        if missing:
            miss_str = " ".join(f".{i:05d}" for i in missing)
            sys.stderr.write(
                f"[WARN] split-file gaps in run {run} (found "
                f"{len(matched)} file(s), max suffix .{last:05d}): "
                f"missing {miss_str}\n"
            )

    if not matched:
        sys.stderr.write(
            f"[WARN] discover_split_files: no files matched 'prad_{run}.evio.*' "
            f"in {directory}\n"
        )
        return [any_path]

    return [path for _, path in matched]


# ============================================================================
# Runinfo (geometry + calibration paths)
# ============================================================================

@dataclass
class RunGeometry:
    """Subset of RunConfig the analysis scripts actually use.  Mirrors the
    fields read by ConfigSetup.h's RotateDetData / TransformDetData."""
    run_number:           int   = 0
    beam_energy:          float = 0.0
    hycal_x:              float = 0.0
    hycal_y:              float = 0.0
    hycal_z:              float = 0.0
    hycal_tilt_x:         float = 0.0
    hycal_tilt_y:         float = 0.0
    hycal_tilt_z:         float = 0.0
    gem_x:                list[float] = field(default_factory=lambda: [0.0]*4)
    gem_y:                list[float] = field(default_factory=lambda: [0.0]*4)
    gem_z:                list[float] = field(default_factory=lambda: [0.0]*4)
    gem_tilt_x:           list[float] = field(default_factory=lambda: [0.0]*4)
    gem_tilt_y:           list[float] = field(default_factory=lambda: [0.0]*4)
    gem_tilt_z:           list[float] = field(default_factory=lambda: [0.0]*4)
    energy_calib_file:    str   = ""
    gem_pedestal_file:    str   = ""
    gem_common_mode_file: str   = ""


def load_run_geometry(runinfo_path: str, eff_run: int) -> RunGeometry:
    """Parse a runinfo JSON and pick the entry with the largest run_number ≤
    eff_run (or the largest entry overall if eff_run ≤ 0)."""
    with open(runinfo_path, "r", encoding="utf-8") as f:
        info = json.load(f)
    configs = info.get("configurations", [])
    if not configs:
        raise RuntimeError(f"runinfo {runinfo_path}: no 'configurations' entries")

    if eff_run > 0:
        ordered = sorted(configs, key=lambda c: c.get("run_number", 0))
        cfg = None
        for c in ordered:
            if c.get("run_number", 0) <= eff_run:
                cfg = c
        if cfg is None:
            cfg = ordered[0]
            sys.stderr.write(
                f"[WARN] no runinfo entry for run {eff_run}; using oldest "
                f"(run {cfg.get('run_number', 0)}).\n"
            )
    else:
        cfg = max(configs, key=lambda c: c.get("run_number", 0))
        sys.stderr.write(
            f"[WARN] no run number specified; using largest runinfo entry "
            f"(run {cfg.get('run_number', 0)}).\n"
        )

    g = RunGeometry()
    g.run_number  = int(cfg.get("run_number", 0))
    g.beam_energy = float(cfg.get("beam_energy", 0.0))

    hc = cfg.get("hycal", {})
    pos  = hc.get("position", [0.0, 0.0, 0.0])
    tilt = hc.get("tilting",  [0.0, 0.0, 0.0])
    g.hycal_x, g.hycal_y, g.hycal_z = (float(v) for v in pos)
    g.hycal_tilt_x, g.hycal_tilt_y, g.hycal_tilt_z = (float(v) for v in tilt)

    gem = cfg.get("gem", {})
    g.gem_pedestal_file    = gem.get("pedestal_file", "")
    g.gem_common_mode_file = gem.get("common_mode_file", "")
    for d in gem.get("detectors", []):
        i = int(d.get("id", -1))
        if 0 <= i < 4:
            pos  = d.get("position", [0.0, 0.0, 0.0])
            tilt = d.get("tilting",  [0.0, 0.0, 0.0])
            g.gem_x[i], g.gem_y[i], g.gem_z[i] = (float(v) for v in pos)
            g.gem_tilt_x[i], g.gem_tilt_y[i], g.gem_tilt_z[i] = \
                (float(v) for v in tilt)

    cal = cfg.get("calibration", {})
    g.energy_calib_file = cal.get("file", "")
    return g


# ============================================================================
# Coordinate transforms (port of analysis::RotateDetData / TransformDetData)
# ============================================================================
# Sign convention matches ConfigSetup.h verbatim — the in-plane translation
# is `x -= ox; y -= oy`, the out-of-plane is `z += oz`.  This matches the
# C++ scripts' lab-frame outputs bit-for-bit.

_DEG2RAD = math.pi / 180.0


def rotate_xyz(x: float, y: float, z: float,
               rx_deg: float, ry_deg: float, rz_deg: float
               ) -> tuple[float, float, float]:
    """Apply Rz then Ry then Rx (extrinsic), each in degrees.  Each axis
    short-circuits when the angle is zero — same as RotateDetData()."""
    if rz_deg != 0.0:
        c = math.cos(rz_deg * _DEG2RAD); s = math.sin(rz_deg * _DEG2RAD)
        x, y = x * c - y * s, x * s + y * c
    if ry_deg != 0.0:
        c = math.cos(ry_deg * _DEG2RAD); s = math.sin(ry_deg * _DEG2RAD)
        x, z =  x * c + z * s, -x * s + z * c
    if rx_deg != 0.0:
        c = math.cos(rx_deg * _DEG2RAD); s = math.sin(rx_deg * _DEG2RAD)
        y, z = y * c - z * s, y * s + z * c
    return x, y, z


def transform_hycal(x: float, y: float, z: float, geo: RunGeometry
                    ) -> tuple[float, float, float]:
    """Detector-frame → lab-frame for HyCal: rotate around (0,0,0) then
    apply (-hycal_x, -hycal_y, +hycal_z)."""
    x, y, z = rotate_xyz(x, y, z,
                         geo.hycal_tilt_x, geo.hycal_tilt_y, geo.hycal_tilt_z)
    return (x - geo.hycal_x, y - geo.hycal_y, z + geo.hycal_z)


def transform_gem(x: float, y: float, z: float, det_id: int, geo: RunGeometry
                  ) -> tuple[float, float, float]:
    """Detector-frame → lab-frame for GEM detector det_id (0..3)."""
    if not (0 <= det_id < 4):
        return x, y, z
    x, y, z = rotate_xyz(x, y, z,
                         geo.gem_tilt_x[det_id],
                         geo.gem_tilt_y[det_id],
                         geo.gem_tilt_z[det_id])
    return (x - geo.gem_x[det_id],
            y - geo.gem_y[det_id],
            z + geo.gem_z[det_id])


def project_to_z(x: float, y: float, z: float, target_z: float
                 ) -> tuple[float, float, float]:
    """Straight-line target→hit projection to z = target_z.  Mirrors
    analysis::GetProjection (single-vertex assumption).  z must be > 0."""
    if z == 0.0:
        return x, y, target_z
    s = target_z / z
    return x * s, y * s, target_z


# ============================================================================
# Pipeline state + setup
# ============================================================================

@dataclass
class Pipeline:
    """Initialized HyCal + GEM systems plus geometry & EVIO file list.  One
    call to setup_pipeline() produces a ready-to-loop bundle."""
    cfg:           "dec.DaqConfig"       = None
    crate_map:     dict[int, int]        = field(default_factory=dict)
    geo:           RunGeometry           = field(default_factory=RunGeometry)
    hycal:         "det.HyCalSystem"     = None
    hc_clusterer:  "det.HyCalCluster"    = None
    wave_ana:      "dec.WaveAnalyzer"    = None
    gem_sys:       "det.GemSystem"       = None
    gem_clusterer: "det.GemCluster"      = None
    evio_files:    list[str]             = field(default_factory=list)


def _print(msg: str) -> None:
    print(msg, flush=True)


def setup_pipeline(*,
                   evio_path: str,
                   max_events: int = 0,
                   run_num: int = -1,
                   gem_ped_file: str = "",
                   gem_cm_file: str = "",
                   hc_calib_file: str = "",
                   daq_config: str = "",
                   gem_map_file: str = "",
                   hc_map_file: str = "",
                   ) -> Pipeline:
    """Mirror the C++ scripts' setup block.  All "_file" args accept "" to
    fall back to runinfo / database defaults.  Run number defaults to a
    sniff from the EVIO basename if -1."""
    p = Pipeline()

    # ---- DAQ config ------------------------------------------------------
    p.cfg = dec.load_daq_config(daq_config)  # "" → installed default
    _print(f"[setup] DAQ config : {daq_config or '(default)'}")
    p.crate_map = {int(roc.tag): int(roc.crate) for roc in p.cfg.roc_tags}

    # ---- runinfo ---------------------------------------------------------
    ri_path = discover_runinfo_path()
    if not ri_path:
        raise SystemExit("[ERROR] no runinfo pointer in database/reconstruction_config.json"
                         " — cannot resolve calibration / geometry.")
    eff_run = run_num
    if eff_run <= 0:
        sniff = extract_run_number(evio_path or "")
        if sniff > 0:
            eff_run = sniff
            _print(f"[setup] Run number : {eff_run} (extracted from filename)")
    else:
        _print(f"[setup] Run number : {eff_run} (caller-provided)")
    p.geo = load_run_geometry(ri_path, eff_run)
    _print(f"[setup] RunInfo    : {ri_path}  beam={p.geo.beam_energy:.0f} MeV  "
           f"hycal_z={p.geo.hycal_z:.1f} mm")

    # ---- HyCal -----------------------------------------------------------
    hc_map = hc_map_file or resolve_db_path("hycal_map.json")
    p.hycal = det.HyCalSystem()
    p.hycal.init(hc_map)

    hc_calib = hc_calib_file or resolve_db_path(p.geo.energy_calib_file)
    if hc_calib:
        n = p.hycal.load_calibration(hc_calib)
        _print(f"[setup] HC calib   : {hc_calib} ({n} modules)")
    else:
        _print("[WARN] no HyCal calibration file — energies will be wrong.")

    p.hc_clusterer = det.HyCalCluster(p.hycal)
    p.hc_clusterer.set_config(det.HyCalClusterConfig())

    p.wave_ana = dec.WaveAnalyzer(p.cfg.wave_cfg)

    # ---- GEM -------------------------------------------------------------
    gem_map = gem_map_file or resolve_db_path("gem_map.json")
    p.gem_sys = det.GemSystem()
    p.gem_sys.init(gem_map)
    _print(f"[setup] GEM map    : {gem_map}  "
           f"({p.gem_sys.get_n_detectors()} detectors)")

    ped_path = gem_ped_file or resolve_db_path(p.geo.gem_pedestal_file)
    cm_path  = gem_cm_file  or resolve_db_path(p.geo.gem_common_mode_file)
    if ped_path:
        p.gem_sys.load_pedestals(ped_path)
        _print(f"[setup] GEM peds   : {ped_path}")
    else:
        _print("[WARN] no GEM pedestal file — full-readout data reconstructs empty.")
    if cm_path:
        p.gem_sys.load_common_mode_range(cm_path)
        _print(f"[setup] GEM CM     : {cm_path}")

    p.gem_clusterer = det.GemCluster()

    # ---- EVIO discovery --------------------------------------------------
    p.evio_files = discover_split_files(evio_path or "")
    _print(f"[setup] EVIO       : {len(p.evio_files)} split file(s) for "
           f"input {evio_path or '(null)'}")
    for f in p.evio_files:
        _print(f"           {f}")
    return p


# ============================================================================
# Per-event reconstruction helpers
# ============================================================================
# Both scripts call these inside their own event loop to keep the pipeline
# idempotent — the user passes the already-decoded EventData / SspEventData
# in (so trigger filtering happens before we pay the reco cost).

# HyCal trigger window (ns) for picking the largest peak — matches the live
# monitor's "100..200" cut in app_state_init.cpp / the C++ scripts.
HC_TIME_LO = 100.0
HC_TIME_HI = 200.0


def reconstruct_hycal(p: Pipeline, fadc_evt) -> list:
    """Run HyCal waveform → energy → cluster on one decoded EventData.
    Returns a list of det.ClusterHit (detector frame)."""
    p.hc_clusterer.clear()
    for ri in range(fadc_evt.nrocs):
        roc = fadc_evt.roc(ri)
        if not roc.present:
            continue
        crate = p.crate_map.get(roc.tag)
        if crate is None:
            continue
        for s in roc.present_slots():
            slot = roc.slot(s)
            for c in slot.present_channels():
                mod = p.hycal.module_by_daq(crate, s, c)
                if mod is None or not mod.is_hycal():
                    continue
                cd = slot.channel(c)
                if cd.nsamples <= 0:
                    continue
                _, _, peaks = p.wave_ana.analyze(cd.samples)
                # Largest peak inside the trigger window.
                best = None
                best_h = -1.0
                for pk in peaks:
                    if HC_TIME_LO < pk.time < HC_TIME_HI and pk.height > best_h:
                        best = pk
                        best_h = pk.height
                if best is None:
                    continue
                p.hc_clusterer.add_hit(mod.index, mod.energize(best.integral))
    p.hc_clusterer.form_clusters()
    return p.hc_clusterer.reconstruct_hits()


def reconstruct_gem(p: Pipeline, ssp_evt) -> None:
    """Run GEM ProcessEvent + Reconstruct.  After this the per-detector
    hit lists are accessible via p.gem_sys.get_hits(d)."""
    p.gem_sys.clear()
    p.gem_sys.process_event(ssp_evt)
    p.gem_sys.reconstruct(p.gem_clusterer)


# ============================================================================
# Argparse helpers
# ============================================================================

def add_common_args(ap: argparse.ArgumentParser) -> None:
    """Register the args every analysis script accepts.  Mirrors the C++
    signatures of gem_hycal_matching / plot_hits_at_hycal."""
    ap.add_argument("evio_path",
                    help="EVIO input.  glob ('prad_NNNNNN.evio.*'), directory, "
                         "or single split (prad_NNNNNN.evio.00000).")
    ap.add_argument("out_path",
                    help="Output table path (.tsv / .csv).")
    ap.add_argument("--max-events", type=int, default=0,
                    help="Stop after N raw physics events (0 = all).")
    ap.add_argument("--run-num", type=int, default=-1,
                    help="Run number override (default -1 = sniff from filename).")
    ap.add_argument("--gem-ped-file",  default="",
                    help='GEM pedestal file (default "" = via runinfo).')
    ap.add_argument("--gem-cm-file",   default="",
                    help='GEM common-mode file (default "" = via runinfo).')
    ap.add_argument("--hc-calib-file", default="",
                    help='HyCal calibration file (default "" = via runinfo).')
    ap.add_argument("--daq-config",    default="",
                    help='DAQ config (default "" = installed default).')
    ap.add_argument("--gem-map-file",  default="",
                    help='GEM map (default "" = database/gem_map.json).')
    ap.add_argument("--hc-map-file",   default="",
                    help='HyCal map (default "" = database/hycal_map.json).')
    ap.add_argument("--csv", action="store_true",
                    help="Emit CSV instead of TSV.")
    ap.add_argument("--no-header", action="store_true",
                    help="Skip the column-name header row.")


# Trigger gate — only events with this exact bitmask are reconstructed.
# Same as the C++ scripts' `trigger_bits == 0x100` filter.
PHYSICS_TRIGGER_BITS = 0x100


def open_table_writer(out_path: str, csv_mode: bool):
    """Open `out_path` for writing.  Returns (file_handle, write_row callable)
    where write_row(seq) writes one row.  Caller closes the file."""
    import csv as _csv
    fh = open(out_path, "w", encoding="utf-8", newline="")
    if csv_mode:
        w = _csv.writer(fh, lineterminator="\n")
        return fh, w.writerow
    sep = "\t"
    def write_row(row):
        fh.write(sep.join("" if v is None else str(v) for v in row))
        fh.write("\n")
    return fh, write_row
