#!/usr/bin/env python3
"""
HyCal Snake Scan -- Module Scanner
===================================
Tkinter GUI that drives the HyCal transporter in a snake pattern so the
beam centres on each scanned module, dwells for a configurable time,
then advances to the next module.

Module positions are loaded from the HyCal module database JSON file,
which contains PbWO4 (inner), PbGlass (outer), and LMS modules.
The scan always includes all PbWO4 modules; the number of surrounding
PbGlass layers to include is configurable in the GUI (0--6).

Usage
-----
    python hycal_snake_scan.py                          # simulation
    python hycal_snake_scan.py --expert                  # expert operator
    python hycal_snake_scan.py --observer                # read-only monitor
    python hycal_snake_scan.py --database /path/to.json  # custom database

Coordinate system
-----------------
    ptrans_x, ptrans_y = (-126.75, 10.11)  -->  beam at HyCal centre (0,0)
    ptrans_x = BEAM_CENTER_X + module_x
    ptrans_y = BEAM_CENTER_Y - module_y

Writable PVs (the ONLY PVs this tool writes to):
    ptrans_x.VAL / ptrans_y.VAL    -- absolute set-point
    ptrans_x.SPMG / ptrans_y.SPMG  -- motor mode  Stop(0) Pause(1) Move(2) Go(3)

Position readback:
    ptrans_x.RBV / ptrans_y.RBV    -- true motor position (read-back value)
    hallb_ptrans_x_encoder / hallb_ptrans_y1_encoder -- raw encoder values

Note: bare ptrans_x / ptrans_y are aliases for ptrans_x.VAL / ptrans_y.VAL
(the set-point), NOT the actual position.

All other PVs are read-only for monitoring.

Observer mode (--observer)
--------------------------
    Read-only monitoring: shows motor position and PV values on the
    canvas, but all scan controls are disabled and PV writes are blocked.

Requirements
------------
    Python 3.8+
    pyepics  (only for --expert / --observer mode)
"""

from __future__ import annotations

import argparse
import math
import os
import random
import threading
import time
import tkinter as tk
from tkinter import ttk, messagebox
from datetime import datetime
from enum import IntEnum
from typing import Dict, List, Optional, Tuple

from scan_utils import (
    C, Module, load_modules, module_to_ptrans, ptrans_to_module,
    ptrans_in_limits, filter_scan_modules,
    BEAM_CENTER_X, BEAM_CENTER_Y, DEFAULT_DB_PATH,
    PTRANS_X_MIN, PTRANS_X_MAX, PTRANS_Y_MIN, PTRANS_Y_MAX,
)


# ============================================================================
#  SCAN-SPECIFIC CONSTANTS
# ============================================================================

DEFAULT_DWELL = 120.0    # seconds
DEFAULT_POS_THRESHOLD = 0.5   # mm
MOVE_TIMEOUT = 300.0     # seconds per single move
MAX_LG_LAYERS = 2
DEFAULT_VELO_X = 50.0    # mm/s  x-axis (fast)
DEFAULT_VELO_Y = 5.0     # mm/s  y-axis (slow)

SPMG_LABELS = {0: "Stop", 1: "Pause", 2: "Move", 3: "Go"}


class SPMG(IntEnum):
    STOP = 0
    PAUSE = 1
    MOVE = 2
    GO = 3


# -- EPICS PV names ----------------------------------------------------------

class PV:
    """All EPICS PV names used by this tool."""
    X_VAL  = "ptrans_x.VAL";   Y_VAL  = "ptrans_y.VAL"
    X_SPMG = "ptrans_x.SPMG";  Y_SPMG = "ptrans_y.SPMG"
    X_ENCODER = "hallb_ptrans_x_encoder"
    Y_ENCODER = "hallb_ptrans_y1_encoder"
    X_RBV  = "ptrans_x.RBV";   Y_RBV  = "ptrans_y.RBV"
    X_MOVN = "ptrans_x.MOVN";  Y_MOVN = "ptrans_y.MOVN"
    X_VELO = "ptrans_x.VELO";  Y_VELO = "ptrans_y.VELO"
    X_ACCL = "ptrans_x.ACCL";  Y_ACCL = "ptrans_y.ACCL"
    X_TDIR = "ptrans_x.TDIR";  Y_TDIR = "ptrans_y.TDIR"
    X_MSTA = "ptrans_x.MSTA";  Y_MSTA = "ptrans_y.MSTA"
    X_ATHM = "ptrans_x.ATHM";  Y_ATHM = "ptrans_y.ATHM"
    X_PREC = "ptrans_x.PREC";  Y_PREC = "ptrans_y.PREC"
    X_BVEL = "ptrans_x.BVEL";  Y_BVEL = "ptrans_y.BVEL"
    X_BACC = "ptrans_x.BACC";  Y_BACC = "ptrans_y.BACC"
    X_VBAS = "ptrans_x.VBAS";  Y_VBAS = "ptrans_y.VBAS"
    X_BDST = "ptrans_x.BDST";  Y_BDST = "ptrans_y.BDST"
    X_FRAC = "ptrans_x.FRAC";  Y_FRAC = "ptrans_y.FRAC"


def _snake_sector(modules: List[Module], going_right: bool,
                  top_to_bottom: bool) -> Tuple[List[Module], bool]:
    """Snake-scan a sector row by row.

    Rows are grouped by y (0.5 mm tolerance), ordered top-to-bottom or
    bottom-to-top.  Within each row, modules are sorted by x, alternating
    direction.  Returns (path, next_going_right).
    """
    if not modules:
        return [], going_right
    by_y = sorted(modules, key=lambda m: -m.y)
    rows: List[List[Module]] = []
    cur_row = [by_y[0]]
    for m in by_y[1:]:
        if abs(m.y - cur_row[0].y) < 0.5:
            cur_row.append(m)
        else:
            rows.append(cur_row)
            cur_row = [m]
    rows.append(cur_row)          # rows sorted top-to-bottom
    if not top_to_bottom:
        rows.reverse()
    path: List[Module] = []
    for row in rows:
        row.sort(key=lambda m: m.x, reverse=not going_right)
        path.extend(row)
        going_right = not going_right
    return path, going_right


def build_scan_path(scan_modules: List[Module]) -> Tuple[List[Module], int]:
    """Build scan path: PbWO4 centre first, then LG sectors.

    Order: Center (top→bottom) → LG sectors visited by nearest corner,
    with scan direction determined by the entry corner.
    This minimises y-axis travel (the slower motor).

    Returns (path, n_unoptimized).
    """
    if not scan_modules:
        return [], 0

    center = [m for m in scan_modules if m.sector == "Center"]
    lg_sectors: Dict[str, List[Module]] = {}
    for m in scan_modules:
        if m.sector in ("Top", "Right", "Bottom", "Left"):
            lg_sectors.setdefault(m.sector, []).append(m)

    path: List[Module] = []
    going = True                  # start going right

    # 1. PbWO4 centre — top to bottom
    if center:
        seg, going = _snake_sector(center, going_right=going,
                                   top_to_bottom=True)
        path.extend(seg)

    # 2. LG sectors — visit by nearest corner each time
    remaining = dict(lg_sectors)   # name -> module list
    while remaining:
        last = path[-1] if path else None
        # Find nearest corner across all remaining sectors
        best_name = None
        best_dist = float('inf')
        best_corner_right = True   # x-direction at entry
        best_corner_down = True    # y-direction at entry
        for name, mods in remaining.items():
            xs = [m.x for m in mods]
            ys = [m.y for m in mods]
            corners = [
                (min(xs), max(ys), True, True),    # top-left  → right, down
                (max(xs), max(ys), False, True),   # top-right → left,  down
                (min(xs), min(ys), True, False),   # bot-left  → right, up
                (max(xs), min(ys), False, False),  # bot-right → left,  up
            ]
            for cx, cy, go_r, go_d in corners:
                if last:
                    d = abs(last.x - cx) + abs(last.y - cy)
                else:
                    d = 0
                if d < best_dist:
                    best_dist = d
                    best_name = name
                    best_corner_right = go_r
                    best_corner_down = go_d

        mods = remaining.pop(best_name)
        seg, going = _snake_sector(mods,
                                   going_right=best_corner_right,
                                   top_to_bottom=best_corner_down)
        path.extend(seg)

    n_unopt = len(scan_modules) - len(path)
    return path, n_unopt


def estimate_scan_time(path: List[Module], start: int, count: int,
                       dwell: float, vx: float = DEFAULT_VELO_X,
                       vy: float = DEFAULT_VELO_Y) -> float:
    """Estimate total scan time from path geometry and motor velocities.

    Both axes move simultaneously so move_time = max(|dx|/vx, |dy|/vy).
    Total = sum(move_times) + n_modules * dwell_time.
    """
    end = min(start + count, len(path)) if count > 0 else len(path)
    if end <= start:
        return 0.0
    total_move = 0.0
    for i in range(start, end - 1):
        dx = abs(path[i + 1].x - path[i].x)
        dy = abs(path[i + 1].y - path[i].y)
        tx = dx / vx if vx > 0 else 0.0
        ty = dy / vy if vy > 0 else 0.0
        total_move += max(tx, ty)
    n_modules = end - start
    return total_move + n_modules * dwell


# ============================================================================
#  EPICS INTERFACES
# ============================================================================

# Map from short key to PV name
_PV_MAP: List[Tuple[str, str]] = [
    # writable
    ("x_val",     PV.X_VAL),     ("y_val",     PV.Y_VAL),
    ("x_spmg",    PV.X_SPMG),   ("y_spmg",    PV.Y_SPMG),
    # read-only
    ("x_encoder", PV.X_ENCODER), ("y_encoder", PV.Y_ENCODER),
    ("x_rbv",     PV.X_RBV),    ("y_rbv",     PV.Y_RBV),
    ("x_movn",    PV.X_MOVN),   ("y_movn",    PV.Y_MOVN),
    ("x_velo",    PV.X_VELO),   ("y_velo",    PV.Y_VELO),
    ("x_accl",    PV.X_ACCL),   ("y_accl",    PV.Y_ACCL),
    ("x_tdir",    PV.X_TDIR),   ("y_tdir",    PV.Y_TDIR),
    ("x_msta",    PV.X_MSTA),   ("y_msta",    PV.Y_MSTA),
    ("x_athm",    PV.X_ATHM),   ("y_athm",    PV.Y_ATHM),
    ("x_prec",    PV.X_PREC),   ("y_prec",    PV.Y_PREC),
    ("x_bvel",    PV.X_BVEL),   ("y_bvel",    PV.Y_BVEL),
    ("x_bacc",    PV.X_BACC),   ("y_bacc",    PV.Y_BACC),
    ("x_vbas",    PV.X_VBAS),   ("y_vbas",    PV.Y_VBAS),
    ("x_bdst",    PV.X_BDST),   ("y_bdst",    PV.Y_BDST),
    ("x_frac",    PV.X_FRAC),   ("y_frac",    PV.Y_FRAC),
]


class RealEPICS:
    """Channel-access interface using pyepics."""

    def __init__(self, writable: bool = False):
        import epics as _epics          # type: ignore
        self._epics = _epics
        self._pvs: Dict[str, object] = {}
        self._writable = writable

    def connect(self) -> Tuple[int, int]:
        for key, pvname in _PV_MAP:
            self._pvs[key] = self._epics.PV(pvname, connection_timeout=5.0)
        time.sleep(2.0)
        n = sum(1 for p in self._pvs.values() if p.connected)
        self._all_connected = (n == len(self._pvs))
        return n, len(self._pvs)

    def disconnected_pvs(self) -> List[str]:
        """Return PV names that failed to connect."""
        return [pvname for key, pvname in _PV_MAP
                if key in self._pvs and not self._pvs[key].connected]

    def get(self, key: str, default=None):
        pv = self._pvs.get(key)
        if pv and pv.connected:
            v = pv.get()
            return v if v is not None else default
        return default

    def put(self, key: str, value) -> bool:
        if not self._writable or not self._all_connected:
            return False
        pv = self._pvs.get(key)
        if pv and pv.connected:
            pv.put(value)
            return True
        return False


class SimulatedEPICS:
    """In-process motor simulation -- no EPICS needed."""

    def __init__(self):
        self._lock = threading.Lock()
        self._x = BEAM_CENTER_X
        self._y = BEAM_CENTER_Y
        self._tx = self._x
        self._ty = self._y
        self._x_spmg = int(SPMG.GO)
        self._y_spmg = int(SPMG.GO)
        self._x_movn = 0
        self._y_movn = 0
        self._x_speed = 0.5       # mm/s  x-axis
        self._y_speed = 0.2       # mm/s  y-axis
        self._moving = False
        self._thread: Optional[threading.Thread] = None

    def connect(self) -> Tuple[int, int]:
        return (0, 0)              # always "OK" in simulation

    def disconnected_pvs(self) -> List[str]:
        return []

    # -- read ----------------------------------------------------------------

    def get(self, key: str, default=None):
        with self._lock:
            return {
                "x_encoder": round(self._x + random.gauss(0, 0.002), 4),
                "y_encoder": round(self._y + random.gauss(0, 0.002), 4),
                "x_rbv": round(self._x, 3),
                "y_rbv": round(self._y, 3),
                "x_val": round(self._tx, 3),
                "y_val": round(self._ty, 3),
                "x_movn": self._x_movn,
                "y_movn": self._y_movn,
                "x_spmg": self._x_spmg,
                "y_spmg": self._y_spmg,
                "x_velo": self._x_speed,
                "y_velo": self._y_speed,
                "x_accl": 0.2,  "y_accl": 1.0,
                "x_tdir": 1 if self._tx >= self._x else 0,
                "y_tdir": 1 if self._ty >= self._y else 0,
                "x_msta": 0x10B, "y_msta": 0x10B,
                "x_athm": int(abs(self._x - BEAM_CENTER_X) < 1.0),
                "y_athm": int(abs(self._y - BEAM_CENTER_Y) < 1.0),
                "x_prec": 3, "y_prec": 3,
                "x_bvel": 1.0, "y_bvel": 1.0,
                "x_bacc": 1.0, "y_bacc": 1.0,
                "x_vbas": 0.5, "y_vbas": 0.5,
                "x_bdst": 0.0, "y_bdst": 0.0,
                "x_frac": 1.0, "y_frac": 1.0,
            }.get(key, default)

    # -- write (ONLY x_val, y_val, x_spmg, y_spmg) --------------------------

    def put(self, key: str, value) -> bool:
        with self._lock:
            if key == "x_val":
                self._tx = float(value)
            elif key == "y_val":
                self._ty = float(value)
            elif key == "x_spmg":
                self._x_spmg = int(value)
            elif key == "y_spmg":
                self._y_spmg = int(value)
            else:
                return False
            self._evaluate_motion()
        return True

    def _evaluate_motion(self):
        """Start / stop simulated motion based on SPMG state."""
        if self._x_spmg == SPMG.STOP or self._y_spmg == SPMG.STOP:
            self._moving = False
            self._x_movn = 0
            self._y_movn = 0
        elif self._x_spmg == SPMG.PAUSE or self._y_spmg == SPMG.PAUSE:
            self._moving = False
            self._x_movn = 0
            self._y_movn = 0
        elif self._x_spmg == SPMG.GO and self._y_spmg == SPMG.GO:
            if not self._moving:
                self._moving = True
                t = threading.Thread(target=self._run_move, daemon=True)
                t.start()

    def _run_move(self):
        dt = 0.02
        while True:
            with self._lock:
                if not self._moving:
                    self._x_movn = 0
                    self._y_movn = 0
                    return
                dx = self._tx - self._x
                dy = self._ty - self._y
                # Move each axis independently at its own speed
                if abs(dx) > 0.001:
                    step_x = min(self._x_speed * dt, abs(dx))
                    self._x += math.copysign(step_x, dx)
                if abs(dy) > 0.001:
                    step_y = min(self._y_speed * dt, abs(dy))
                    self._y += math.copysign(step_y, dy)
                self._x_movn = 1 if abs(self._tx - self._x) > 0.001 else 0
                self._y_movn = 1 if abs(self._ty - self._y) > 0.001 else 0
                if self._x_movn == 0 and self._y_movn == 0:
                    self._x, self._y = self._tx, self._ty
                    self._moving = False
                    return
            time.sleep(dt)


class ObserverEPICS:
    """Read-only wrapper around RealEPICS -- all puts are silently blocked."""

    def __init__(self):
        self._real = RealEPICS()

    def connect(self) -> Tuple[int, int]:
        return self._real.connect()

    def disconnected_pvs(self) -> List[str]:
        return self._real.disconnected_pvs()

    def get(self, key: str, default=None):
        return self._real.get(key, default)

    def put(self, key: str, value) -> bool:
        return False


# -- helpers shared by both interfaces --------------------------------------

def epics_move_to(ep, x: float, y: float) -> bool:
    """Command a move, clamped to transporter limits.

    Returns False (and does not move) if the target is outside limits.
    """
    if not ptrans_in_limits(x, y):
        return False
    ep.put("x_val", x)
    ep.put("y_val", y)
    ep.put("x_spmg", int(SPMG.GO))
    ep.put("y_spmg", int(SPMG.GO))
    return True

def epics_stop(ep):
    ep.put("x_spmg", int(SPMG.STOP))
    ep.put("y_spmg", int(SPMG.STOP))

def epics_pause(ep):
    ep.put("x_spmg", int(SPMG.PAUSE))
    ep.put("y_spmg", int(SPMG.PAUSE))

def epics_resume(ep):
    ep.put("x_spmg", int(SPMG.GO))
    ep.put("y_spmg", int(SPMG.GO))

def epics_is_moving(ep) -> bool:
    return bool(ep.get("x_movn", 0)) or bool(ep.get("y_movn", 0))

def epics_read_rbv(ep) -> Tuple[float, float]:
    """Read true motor position from RBV (read-back value) PVs."""
    return (ep.get("x_rbv", 0.0), ep.get("y_rbv", 0.0))


# ============================================================================
#  SCAN ENGINE  (runs in a background thread)
# ============================================================================

class ScanState:
    IDLE      = "IDLE"
    MOVING    = "MOVING"
    DWELLING  = "DWELLING"
    PAUSED    = "PAUSED"
    ERROR     = "ERROR"
    COMPLETED = "COMPLETED"


class ScanEngine:
    """Drives the snake scan in a background thread."""

    def __init__(self, epics, modules: List[Module], log_fn):
        self.ep = epics
        self.all_modules = modules
        self.path, n_unopt = build_scan_path(modules)
        self.log = log_fn                # log_fn(msg, level="info")
        if n_unopt:
            self.log(f"WARNING: {n_unopt} modules appended with unoptimized "
                     f"path (no structured route found)", level="warn")

        # -- tunables (set before start) --
        self.dwell_time: float = DEFAULT_DWELL
        self.pos_threshold: float = DEFAULT_POS_THRESHOLD

        # -- runtime state (read by GUI) --
        self.state: str = ScanState.IDLE
        self.current_idx: int = 0        # index into self.path
        self.dwell_remaining: float = 0.0
        self.completed: set = set()      # indices that finished successfully
        self.error_modules: set = set()  # indices that had position errors

        # -- thread control --
        self._thread: Optional[threading.Thread] = None
        self._stop = threading.Event()
        self._skip = threading.Event()
        self._paused = False
        self._ack_error = threading.Event()   # set by GUI to acknowledge error

    # -- public API (called from GUI / main thread) --------------------------

    @property
    def current_module(self) -> Optional[Module]:
        if 0 <= self.current_idx < len(self.path):
            return self.path[self.current_idx]
        return None

    @property
    def progress_text(self) -> str:
        done = len(self.completed)
        total = getattr(self, '_end_idx', len(self.path)) - self.current_idx \
                + len(self.completed)
        return f"{done}/{total}"

    @property
    def eta_seconds(self) -> float:
        if self.state == ScanState.IDLE:
            return 0.0
        end = getattr(self, '_end_idx', len(self.path))
        next_idx = self.current_idx + 1
        if next_idx >= end:
            return self.dwell_remaining
        remaining = end - next_idx
        return estimate_scan_time(
            self.path, next_idx, remaining, self.dwell_time,
        ) + self.dwell_remaining

    def start(self, start_idx: int = 0, count: int = 0):
        if self._thread and self._thread.is_alive():
            return
        self._stop.clear()
        self._skip.clear()
        self._paused = False
        self.current_idx = start_idx
        self._start_idx = start_idx
        self._end_idx = min(start_idx + count, len(self.path)) \
                        if count > 0 else len(self.path)
        self.completed.clear()
        self.error_modules.clear()
        self.state = ScanState.MOVING
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def resume_scan(self):
        """Resume from pause (during move or dwell)."""
        if self.state == ScanState.PAUSED:
            self._paused = False
            epics_resume(self.ep)
            self.state = ScanState.MOVING
            self.log("Scan resumed")

    def pause_scan(self):
        if self.state in (ScanState.MOVING, ScanState.DWELLING):
            self._paused = True
            epics_pause(self.ep)
            self.state = ScanState.PAUSED
            self.log("Scan paused", level="warn")

    def stop_scan(self):
        self._stop.set()
        self._skip.set()          # unblock any waits
        self._ack_error.set()     # unblock error wait
        self._paused = False
        epics_stop(self.ep)
        self.log("Scan stopped", level="warn")

    def skip_module(self):
        self._skip.set()

    def acknowledge_error(self):
        """User acknowledges position error -- scan continues."""
        self._ack_error.set()

    # -- background thread ---------------------------------------------------

    def _run(self):
        n_to_scan = self._end_idx - self.current_idx
        self.log(f"Scan started from {self.path[self.current_idx].name}, "
                 f"dwell {self.dwell_time:.0f}s, {n_to_scan} modules")
        try:
            for i in range(self.current_idx, self._end_idx):
                if self._stop.is_set():
                    break
                self.current_idx = i
                mod = self.path[i]
                px, py = module_to_ptrans(mod.x, mod.y)

                # -- move --
                self.state = ScanState.MOVING
                self.log(f"[{i+1}/{len(self.path)}] Moving to {mod.name} "
                         f"  ptrans({px:.3f}, {py:.3f})")
                if not epics_move_to(self.ep, px, py):
                    self.log(f"SKIPPED {mod.name}: ptrans({px:.3f}, {py:.3f}) "
                             f"outside travel limits", level="warn")
                    continue

                if not self._wait_move_done(px, py):
                    break

                # -- check position error --
                rbv_x, rbv_y = epics_read_rbv(self.ep)
                err = math.sqrt((rbv_x - px)**2 + (rbv_y - py)**2)
                if err > self.pos_threshold:
                    self.error_modules.add(i)
                    self.state = ScanState.ERROR
                    self.log(f"POSITION ERROR at {mod.name}: "
                             f"error={err:.3f} mm  (threshold {self.pos_threshold})",
                             level="error")
                    self._ack_error.clear()
                    self._ack_error.wait()        # block until user acknowledges
                    if self._stop.is_set():
                        break

                # -- dwell --
                self.state = ScanState.DWELLING
                self.dwell_remaining = self.dwell_time
                self.log(f"Dwelling at {mod.name} for {self.dwell_time:.0f}s")
                result = self._wait_dwell()
                if result == "stop":
                    break

                self.completed.add(i)
                self.dwell_remaining = 0.0

        finally:
            if self._stop.is_set():
                self.state = ScanState.IDLE
                self.log("Scan stopped by user")
            elif self.current_idx >= self._end_idx - 1:
                self.state = ScanState.COMPLETED
                self.log("Scan COMPLETE!", level="warn")
            else:
                self.state = ScanState.IDLE

    def _wait_move_done(self, target_x: float, target_y: float) -> bool:
        """Wait until MOVN=0 on both axes.  Returns False if stopped."""
        t0 = time.time()
        while not self._stop.is_set():
            # handle pause
            while self._paused and not self._stop.is_set():
                self.state = ScanState.PAUSED
                time.sleep(0.1)
            if self._stop.is_set():
                return False
            self.state = ScanState.MOVING
            if not epics_is_moving(self.ep):
                return True
            if time.time() - t0 > MOVE_TIMEOUT:
                self.log(f"MOVE TIMEOUT after {MOVE_TIMEOUT:.0f}s", level="error")
                return False
            time.sleep(0.1)
        return False

    def _wait_dwell(self) -> str:
        """Wait for dwell_time seconds.  Returns 'done', 'skip', or 'stop'.

        On resume after pause, the dwell timer restarts from the full
        dwell_time (beam trip recovery — partial data is discarded).
        """
        end = time.time() + self.dwell_time
        while time.time() < end:
            if self._stop.is_set():
                return "stop"
            if self._skip.is_set():
                self._skip.clear()
                self.log("Module skipped by user")
                return "skip"
            if self._paused:
                self.state = ScanState.PAUSED
                while self._paused and not self._stop.is_set():
                    time.sleep(0.1)
                if self._stop.is_set():
                    return "stop"
                # Restart full dwell after resume
                end = time.time() + self.dwell_time
                self.log("Dwell restarted after resume")
            self.state = ScanState.DWELLING
            self.dwell_remaining = max(0.0, end - time.time())
            time.sleep(0.1)
        return "done"


# ============================================================================
#  GUI
# ============================================================================

class SnakeScanGUI:

    CANVAS_SIZE = 680       # pixels
    CANVAS_PAD  = 8
    MOD_SHRINK  = 0.90      # render modules at 90% size for visual gaps

    AUTOGEN = "(autogen)"

    def __init__(self, root: tk.Tk, epics, simulation: bool,
                 all_modules: List[Module],
                 profiles: Optional[Dict[str, List[str]]] = None,
                 observer: bool = False):
        self.root = root
        self.ep = epics
        self.simulation = simulation
        self.observer = observer
        self.all_modules = all_modules
        self._profiles = profiles or {}
        self._active_profile = self.AUTOGEN
        self._lg_layers = 0

        # Precompute PbWO4 bounding box and PbGlass module size
        pwo4 = [m for m in all_modules if m.mod_type == "PbWO4"]
        self._pwo4_min_x = min(m.x for m in pwo4)
        self._pwo4_max_x = max(m.x for m in pwo4)
        self._pwo4_min_y = min(m.y for m in pwo4)
        self._pwo4_max_y = max(m.y for m in pwo4)
        glass = [m for m in all_modules if m.mod_type == "PbGlass"]
        self._lg_sx = glass[0].sx if glass else 38.15
        self._lg_sy = glass[0].sy if glass else 38.15

        self._log_lines: List[str] = []
        self._log_text = None   # set later by _build_ui

        # Open log file (one per session, timestamped)
        log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
        os.makedirs(log_dir, exist_ok=True)
        log_name = datetime.now().strftime("snake_scan_%Y%m%d_%H%M%S.log")
        self._log_file = open(os.path.join(log_dir, log_name), "w")

        # Split into scan targets vs display-only
        self.scan_modules = self._filter_scan_modules(0)

        self.engine = ScanEngine(epics, self.scan_modules, self._log)

        # module name -> path index (scan modules only)
        self._scan_name_to_idx: Dict[str, int] = {
            m.name: i for i, m in enumerate(self.engine.path)
        }
        self._scan_names: set = {m.name for m in self.scan_modules}

        self._selected_start_idx = 0

        # module name -> Module for tooltip lookup
        self._mod_by_name: Dict[str, Module] = {m.name: m for m in all_modules}

        # canvas item IDs:  module name -> rectangle id
        self._cell_ids: Dict[str, int] = {}
        self._display_greyed = False   # track greyed-out state for transitions

        # canvas coordinate mapping (computed in _build_canvas)
        self._scale = 1.0
        self._ox = 0.0
        self._oy = 0.0
        self._x_min = 0.0
        self._y_max = 0.0

        self._build_ui()
        if self.observer:
            self._disable_controls()
        # Expert mode: if any PVs are disconnected, puts are already blocked
        # in RealEPICS — disable the GUI controls and show which PVs failed.
        if not self.simulation and not self.observer:
            disconnected = self.ep.disconnected_pvs()
            if disconnected:
                self._disable_controls()
                msg = "The following PVs are not connected:\n\n"
                msg += "\n".join(f"  {pv}" for pv in disconnected)
                msg += "\n\nAll controls are disabled until every PV connects."
                messagebox.showerror("PV Connection Error", msg)
        self._poll()

    def _filter_scan_modules(self, lg_layers: int) -> List[Module]:
        return filter_scan_modules(self.all_modules, lg_layers,
                                   self._lg_sx, self._lg_sy)

    def _display_color(self, mod_type: str) -> str:
        """Static colour for display-only (non-scanned) modules."""
        if mod_type == "PbGlass":
            return C.MOD_GLASS
        elif mod_type == "PbWO4":
            return C.MOD_PWO4_BG
        return C.MOD_LMS

    # -----------------------------------------------------------------------
    #  UI construction
    # -----------------------------------------------------------------------

    def _build_ui(self):
        if self.observer:
            mode_suffix = "  [OBSERVER]"
        elif self.simulation:
            mode_suffix = "  [SIMULATION]"
        else:
            mode_suffix = "  [EXPERT OPERATOR]"
        self.root.title("HyCal Snake Scan" + mode_suffix)
        self.root.configure(bg=C.BG)
        self.root.resizable(True, True)

        style = ttk.Style()
        style.theme_use("clam")
        style.configure(".", background=C.BG, foreground=C.TEXT,
                         fieldbackground=C.PANEL, bordercolor=C.BORDER)
        style.configure("TLabel", background=C.BG, foreground=C.TEXT)
        style.configure("TLabelframe", background=C.BG, foreground=C.ACCENT)
        style.configure("TLabelframe.Label", background=C.BG,
                         foreground=C.ACCENT, font=("Consolas", 9, "bold"))
        style.configure("TButton", background=C.PANEL, foreground=C.TEXT,
                         padding=4)
        style.map("TButton",
                  background=[("active", C.BORDER)],
                  foreground=[("disabled", "#484f58")])
        style.configure("Accent.TButton", background="#1f6feb",
                         foreground="white")
        style.map("Accent.TButton",
                  background=[("active", "#388bfd")])
        style.configure("Danger.TButton", background="#da3633",
                         foreground="white")
        style.map("Danger.TButton",
                  background=[("active", "#f85149")])
        style.configure("Warn.TButton", background="#9e6a03",
                         foreground="white")
        style.map("Warn.TButton",
                  background=[("active", "#d29922")])
        style.configure("Green.TButton", background="#238636",
                         foreground="white")
        style.map("Green.TButton",
                  background=[("active", "#3fb950")])
        style.configure("TCombobox", fieldbackground=C.PANEL,
                         background=C.BORDER, foreground=C.TEXT,
                         selectbackground=C.BORDER,
                         selectforeground=C.TEXT)
        style.map("TCombobox",
                  fieldbackground=[("readonly", C.PANEL),
                                   ("disabled", C.BG)],
                  foreground=[("disabled", "#484f58")])

        # -- top status bar --------------------------------------------------
        top = tk.Frame(self.root, bg="#0d1520", height=32)
        top.pack(fill="x")
        tk.Label(top, text="  HYCAL SNAKE SCAN  ",
                 bg="#0d1520", fg=C.GREEN,
                 font=("Consolas", 13, "bold")).pack(side="left", padx=8)
        if self.observer:
            mode_text = "OBSERVER (read-only)"
            mode_fg = C.ORANGE
        elif self.simulation:
            mode_text = "SIMULATION"
            mode_fg = C.YELLOW
        else:
            mode_text = "EXPERT OPERATOR"
            mode_fg = C.GREEN
        tk.Label(top, text=mode_text, bg="#0d1520", fg=mode_fg,
                 font=("Consolas", 9, "bold")).pack(side="left", padx=4)

        self._lbl_state = tk.Label(top, text="IDLE", bg="#0d1520",
                                    fg=C.DIM, font=("Consolas", 10, "bold"))
        self._lbl_state.pack(side="right", padx=12)

        # -- main area -------------------------------------------------------
        main = tk.Frame(self.root, bg=C.BG)
        main.pack(fill="both", expand=True, padx=6, pady=4)

        # column 0: canvas
        left = tk.Frame(main, bg=C.BG)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 4))

        # column 1: controls + status
        right = tk.Frame(main, bg=C.BG)
        right.grid(row=0, column=1, sticky="nsew")

        main.columnconfigure(0, weight=0)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        self._build_canvas(left)
        self._build_controls(right)

        # -- bottom: log -----------------------------------------------------
        log_frame = ttk.LabelFrame(self.root, text=" Event Log ")
        log_frame.pack(fill="both", expand=False, padx=6, pady=(0, 6))

        self._log_text = tk.Text(log_frame, height=8, bg="#0d1117",
                                  fg=C.DIM, font=("Consolas", 9),
                                  wrap="word", state="disabled",
                                  borderwidth=0, highlightthickness=0)
        self._log_text.pack(fill="both", expand=True, padx=2, pady=2)
        self._log_text.tag_configure("info", foreground=C.TEXT)
        self._log_text.tag_configure("warn", foreground=C.YELLOW)
        self._log_text.tag_configure("error", foreground=C.RED)

    # -- observer: disable all controls ---------------------------------------

    def _disable_controls(self):
        """Disable all scan and direct-control widgets for observer mode."""
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_skip, self._btn_ack):
            w.configure(state="disabled")
        # Direct control buttons live inside the Direct Control frame;
        # walk the widget tree and disable all buttons there.
        for child in self.root.winfo_children():
            self._disable_buttons_recursive(child)
        # Disable input widgets
        for w in (self._start_combo, self._profile_combo):
            w.configure(state="disabled")
        for w in (self._lg_layers_spin, self._count_entry):
            w.configure(state="disabled")

    def _disable_buttons_recursive(self, widget):
        """Recursively disable all ttk.Button widgets."""
        if isinstance(widget, ttk.Button):
            widget.configure(state="disabled")
        for child in widget.winfo_children():
            self._disable_buttons_recursive(child)

    # -- canvas (module map) -------------------------------------------------

    def _build_canvas(self, parent):
        self._canvas_frame = ttk.LabelFrame(parent, text="")
        self._update_canvas_label()
        self._canvas_frame.pack(fill="both", expand=True)

        sz = self.CANVAS_SIZE
        self._canvas = tk.Canvas(self._canvas_frame, width=sz, height=sz,
                                  bg="#0a0e14", highlightthickness=0)
        self._canvas.pack(padx=4, pady=4)
        self._canvas.bind("<Button-1>", self._on_canvas_click)
        self._selected_mod_name = None

        self._compute_canvas_mapping()
        self._draw_modules()
        self._draw_motor_marker()

        # legend
        leg = tk.Frame(self._canvas_frame, bg=C.BG)
        leg.pack(fill="x", padx=4, pady=(0, 4))
        legend_items = [
            ("Todo", C.MOD_TODO), ("Skipped", C.MOD_SKIPPED),
            ("Moving", C.MOD_CURRENT), ("Dwell", C.MOD_DWELL),
            ("Done", C.MOD_DONE), ("Error", C.MOD_ERROR),
            ("Start", C.MOD_SELECTED), ("PbGlass", C.MOD_GLASS),
        ]
        for label, colour in legend_items:
            tk.Canvas(leg, width=10, height=10, bg=colour,
                      highlightthickness=0).pack(side="left", padx=(6, 1))
            tk.Label(leg, text=label, bg=C.BG, fg=C.DIM,
                     font=("Consolas", 8)).pack(side="left")

    def _update_canvas_label(self):
        n_pwo4 = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        if n_lg:
            base = f"Module Map ({n_pwo4} PbWO4 + {n_lg} PbGlass)"
        else:
            base = f"Module Map ({n_pwo4} PbWO4)"
        sel_name = getattr(self, '_selected_mod_name', None)
        if sel_name:
            mod = self._mod_by_name.get(sel_name)
            if mod:
                px, py = module_to_ptrans(mod.x, mod.y)
                base += (f" | {mod.name} ({mod.mod_type})"
                         f"  HyCal({mod.x:.1f}, {mod.y:.1f})"
                         f"  ptrans({px:.1f}, {py:.1f})")
        self._canvas_frame.configure(text=f" {base} ")

    def _compute_canvas_mapping(self):
        """Compute scale and offset to map HyCal mm -> canvas pixels."""
        # Exclude LMS from bounding box (they are off to the side and not drawn)
        drawn = [m for m in self.all_modules if m.mod_type != "LMS"]
        if not drawn:
            return
        x_min = min(m.x - m.sx / 2 for m in drawn)
        x_max = max(m.x + m.sx / 2 for m in drawn)
        y_min = min(m.y - m.sy / 2 for m in drawn)
        y_max = max(m.y + m.sy / 2 for m in drawn)

        usable = self.CANVAS_SIZE - 2 * self.CANVAS_PAD
        self._scale = min(usable / (x_max - x_min),
                          usable / (y_max - y_min))

        draw_w = (x_max - x_min) * self._scale
        draw_h = (y_max - y_min) * self._scale
        self._ox = self.CANVAS_PAD + (usable - draw_w) / 2
        self._oy = self.CANVAS_PAD + (usable - draw_h) / 2
        self._x_min = x_min
        self._y_max = y_max

    def _mod_to_canvas(self, m: Module) -> Tuple[float, float, float, float]:
        """Module -> canvas rectangle (x0, y0, x1, y1)."""
        cx = self._ox + (m.x - self._x_min) * self._scale
        cy = self._oy + (self._y_max - m.y) * self._scale
        hw = m.sx * self._scale * self.MOD_SHRINK / 2
        hh = m.sy * self._scale * self.MOD_SHRINK / 2
        return (cx - hw, cy - hh, cx + hw, cy + hh)

    def _mod_to_canvas_center(self, m: Module) -> Tuple[float, float]:
        """Module -> canvas centre point (cx, cy)."""
        cx = self._ox + (m.x - self._x_min) * self._scale
        cy = self._oy + (self._y_max - m.y) * self._scale
        return (cx, cy)

    def _draw_modules(self):
        """Draw all modules on the canvas."""
        self._canvas.delete("all")
        self._cell_ids.clear()

        # Draw display-only modules first (background layer)
        for m in self.all_modules:
            if m.name in self._scan_names:
                continue
            if m.mod_type == "LMS":
                continue  # skip LMS on map
            x0, y0, x1, y1 = self._mod_to_canvas(m)
            color = self._display_color(m.mod_type)
            rid = self._canvas.create_rectangle(
                x0, y0, x1, y1, fill=color, outline="", width=0,
                tags=(f"mod_{m.name}", "display"))
            self._cell_ids[m.name] = rid

        # Draw scan modules on top
        for m in self.scan_modules:
            x0, y0, x1, y1 = self._mod_to_canvas(m)
            rid = self._canvas.create_rectangle(
                x0, y0, x1, y1, fill=C.MOD_TODO, outline="", width=0,
                tags=(f"mod_{m.name}", "scan"))
            self._cell_ids[m.name] = rid

        # Draw transporter travel-limit boundary (red dashed)
        # Convert ptrans limits to HyCal coordinates, then to canvas
        lim_hx_min = PTRANS_X_MIN - BEAM_CENTER_X  # left-most beam x
        lim_hx_max = PTRANS_X_MAX - BEAM_CENTER_X  # right-most beam x
        lim_hy_min = BEAM_CENTER_Y - PTRANS_Y_MAX  # bottom-most beam y
        lim_hy_max = BEAM_CENTER_Y - PTRANS_Y_MIN  # top-most beam y
        bx0 = self._ox + (lim_hx_min - self._x_min) * self._scale
        by0 = self._oy + (self._y_max - lim_hy_max) * self._scale
        bx1 = self._ox + (lim_hx_max - self._x_min) * self._scale
        by1 = self._oy + (self._y_max - lim_hy_min) * self._scale
        self._canvas.create_rectangle(
            bx0, by0, bx1, by1,
            outline=C.RED, width=1, dash=(4, 4), tags=("limit_box",))

    def _draw_path_preview(self):
        """Draw the projected snake path line from start to end on the canvas."""
        self._clear_path_preview()
        path = self.engine.path
        start = self._selected_start_idx
        if start >= len(path):
            return
        count = self._count_var.get() if hasattr(self, '_count_var') else 0
        end = min(start + count, len(path)) if count > 0 else len(path)
        # Build coordinate list for the active segment
        coords = []
        for i in range(start, end):
            coords.extend(self._mod_to_canvas_center(path[i]))
        if len(coords) >= 4:
            self._canvas.create_line(
                *coords, fill=C.PATH_LINE, width=1.5,
                smooth=False, tags=("path_preview",))
            # Line on top of everything so it's visible over modules
            self._canvas.tag_raise("path_preview")

    def _clear_path_preview(self):
        self._canvas.delete("path_preview")

    def _on_canvas_click(self, event):
        items = self._canvas.find_closest(event.x, event.y)
        if not items:
            self._clear_selection()
            return
        tags = self._canvas.gettags(items[0])
        name = None
        for tag in tags:
            if tag.startswith("mod_"):
                name = tag[4:]
                break
        if name is None:
            self._clear_selection()
            return

        # Toggle: clicking the same module deselects
        if self._selected_mod_name == name:
            self._clear_selection()
            return

        mod = self._mod_by_name.get(name)
        if mod is None:
            self._clear_selection()
            return

        self._selected_mod_name = name

        # Set start module if it's a scan module
        if name in self._scan_name_to_idx:
            idx = self._scan_name_to_idx[name]
            self._selected_start_idx = idx
            self._start_var.set(name)
            self._log(f"Selected start module: {name}")
            self._draw_path_preview()

        # Highlight border
        self._canvas.delete("mod_highlight")
        x0, y0, x1, y1 = self._mod_to_canvas(mod)
        self._canvas.create_rectangle(
            x0, y0, x1, y1, outline=C.ACCENT, width=2,
            tags=("mod_highlight",))
        self._canvas.tag_raise("mod_highlight")

        # Update frame label with selected module info
        self._update_canvas_label()

    def _clear_selection(self):
        self._canvas.delete("mod_highlight")
        self._selected_mod_name = None
        self._update_canvas_label()

    def _update_canvas(self):
        if self.observer:
            self._update_canvas_observer()
            return

        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING,
                                ScanState.PAUSED, ScanState.ERROR)

        # Grey out / restore display-only modules on state transitions
        if running and not self._display_greyed:
            self._display_greyed = True
            for m in self.all_modules:
                if m.name not in self._scan_names and m.name in self._cell_ids:
                    self._canvas.itemconfigure(
                        self._cell_ids[m.name], fill=C.MOD_EXCLUDED)
        elif not running and self._display_greyed:
            self._display_greyed = False
            for m in self.all_modules:
                if m.name not in self._scan_names and m.name in self._cell_ids:
                    self._canvas.itemconfigure(
                        self._cell_ids[m.name],
                        fill=self._display_color(m.mod_type))

        # Update scan module colours
        idle = eng.state in (ScanState.IDLE, ScanState.COMPLETED)
        count = self._count_var.get() if hasattr(self, '_count_var') else 0
        start_idx = self._selected_start_idx
        end_idx = min(start_idx + count, len(eng.path)) \
                  if count > 0 else len(eng.path)
        # During scan, use the engine's range instead
        if not idle:
            start_idx = eng.current_idx  # already advancing
            end_idx = getattr(eng, '_end_idx', len(eng.path))
        for i, mod in enumerate(eng.path):
            rid = self._cell_ids.get(mod.name)
            if rid is None:
                continue
            if i == eng.current_idx and eng.state == ScanState.DWELLING:
                colour = C.MOD_DWELL
            elif i == eng.current_idx and eng.state in (ScanState.MOVING,
                                                         ScanState.PAUSED):
                colour = C.MOD_CURRENT
            elif i in eng.error_modules:
                colour = C.MOD_ERROR
            elif i in eng.completed:
                colour = C.MOD_DONE
            elif idle and i == self._selected_start_idx:
                colour = C.MOD_SELECTED
            elif i < start_idx or i >= end_idx:
                colour = C.MOD_SKIPPED
            else:
                colour = C.MOD_TODO
            self._canvas.itemconfigure(rid, fill=colour)

        # Path preview line
        self._canvas.delete("path_preview")
        if idle:
            self._draw_path_preview()
        elif running:
            # Dashed future path from next module to end
            coords = []
            for i in range(eng.current_idx + 1, end_idx):
                coords.extend(self._mod_to_canvas_center(eng.path[i]))
            if len(coords) >= 4:
                self._canvas.create_line(
                    *coords, fill=C.PATH_LINE, width=1,
                    dash=(3, 3), tags=("path_preview",))

        self._draw_motor_marker()

    def _update_canvas_observer(self):
        """Observer mode: only draw the motor position marker, no module colouring."""
        self._draw_motor_marker()

    def _draw_motor_marker(self):
        """Draw the red crosshair at the current motor position."""
        self._canvas.delete("motor_pos")
        rx = self.ep.get("x_rbv", BEAM_CENTER_X)
        ry = self.ep.get("y_rbv", BEAM_CENTER_Y)
        hx, hy = ptrans_to_module(rx, ry)
        cx = self._ox + (hx - self._x_min) * self._scale
        cy = self._oy + (self._y_max - hy) * self._scale
        r = 5
        self._canvas.create_line(
            cx - r, cy, cx + r, cy,
            fill=C.RED, width=1.5, tags=("motor_pos",))
        self._canvas.create_line(
            cx, cy - r, cx, cy + r,
            fill=C.RED, width=1.5, tags=("motor_pos",))
        self._canvas.tag_raise("motor_pos")

    # -- controls panel ------------------------------------------------------

    def _build_controls(self, parent):
        # === Scan Control ===
        sc = ttk.LabelFrame(parent, text=" Scan Control ")
        sc.pack(fill="x", pady=(0, 4))

        # Path profile selector
        r_pp = tk.Frame(sc, bg=C.BG)
        r_pp.pack(fill="x", padx=6, pady=2)
        tk.Label(r_pp, text="Path:", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        profile_names = [self.AUTOGEN] + sorted(self._profiles.keys())
        self._profile_var = tk.StringVar(value=self.AUTOGEN)
        self._profile_combo = ttk.Combobox(
            r_pp, textvariable=self._profile_var,
            values=profile_names, width=18, state="readonly",
            font=("Consolas", 9))
        self._profile_combo.pack(side="right")
        self._profile_combo.bind("<<ComboboxSelected>>",
                                  self._on_path_profile_changed)

        # LG layers (own row, only effective for autogen)
        r_lg = tk.Frame(sc, bg=C.BG)
        r_lg.pack(fill="x", padx=6, pady=2)
        tk.Label(r_lg, text="LG layers (0-2):", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._lg_layers_var = tk.IntVar(value=self._lg_layers)
        self._lg_layers_spin = tk.Spinbox(
            r_lg, from_=0, to=MAX_LG_LAYERS,
            textvariable=self._lg_layers_var,
            width=4, bg=C.PANEL, fg=C.TEXT, font=("Consolas", 9),
            buttonbackground=C.BORDER, insertbackground=C.TEXT,
            command=self._on_lg_layers_changed)
        self._lg_layers_spin.pack(side="right")

        # Start module + count (same row)
        r_sc = tk.Frame(sc, bg=C.BG)
        r_sc.pack(fill="x", padx=6, pady=2)
        tk.Label(r_sc, text="Start:", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        names = [m.name for m in self.engine.path]
        self._start_var = tk.StringVar(value=names[0] if names else "")
        self._start_combo = ttk.Combobox(r_sc, textvariable=self._start_var,
                                          values=names, width=8,
                                          font=("Consolas", 9))
        self._start_combo.pack(side="left", padx=(2, 8))
        tk.Label(r_sc, text="Count:", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._count_var = tk.IntVar(value=0)
        self._count_entry = tk.Spinbox(
            r_sc, from_=0, to=len(names),
            textvariable=self._count_var,
            width=6, bg=C.PANEL, fg=C.TEXT, font=("Consolas", 9),
            buttonbackground=C.BORDER, insertbackground=C.TEXT,
            command=self._draw_path_preview)
        self._count_entry.pack(side="right")
        # Dark theme for dropdown listbox
        self.root.option_add("*TCombobox*Listbox.background", C.PANEL)
        self.root.option_add("*TCombobox*Listbox.foreground", C.TEXT)
        self.root.option_add("*TCombobox*Listbox.selectBackground", C.ACCENT)
        self.root.option_add("*TCombobox*Listbox.selectForeground", "white")
        self._start_combo.bind("<<ComboboxSelected>>", self._on_start_selected)

        # dwell
        r = tk.Frame(sc, bg=C.BG)
        r.pack(fill="x", padx=6, pady=2)
        tk.Label(r, text="Dwell time (s):", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._dwell_var = tk.DoubleVar(value=DEFAULT_DWELL)
        tk.Spinbox(r, from_=1, to=9999, textvariable=self._dwell_var,
                   width=8, bg=C.PANEL, fg=C.TEXT, font=("Consolas", 9),
                   buttonbackground=C.BORDER,
                   insertbackground=C.TEXT).pack(side="right")

        # threshold
        r2 = tk.Frame(sc, bg=C.BG)
        r2.pack(fill="x", padx=6, pady=2)
        tk.Label(r2, text="Pos. threshold (mm):", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._thresh_var = tk.DoubleVar(value=DEFAULT_POS_THRESHOLD)
        tk.Spinbox(r2, from_=0.01, to=10.0, increment=0.1,
                   textvariable=self._thresh_var,
                   width=8, bg=C.PANEL, fg=C.TEXT, font=("Consolas", 9),
                   buttonbackground=C.BORDER,
                   insertbackground=C.TEXT).pack(side="right")

        # buttons
        bf = tk.Frame(sc, bg=C.BG)
        bf.pack(fill="x", padx=6, pady=6)

        self._btn_start = ttk.Button(bf, text="Start Scan",
                                      style="Green.TButton",
                                      command=self._cmd_start)
        self._btn_start.pack(side="left", expand=True, fill="x", padx=2)

        self._btn_pause = ttk.Button(bf, text="Pause",
                                      style="Warn.TButton",
                                      command=self._cmd_pause)
        self._btn_pause.pack(side="left", expand=True, fill="x", padx=2)

        self._btn_stop = ttk.Button(bf, text="Stop",
                                     style="Danger.TButton",
                                     command=self._cmd_stop)
        self._btn_stop.pack(side="left", expand=True, fill="x", padx=2)

        bf2 = tk.Frame(sc, bg=C.BG)
        bf2.pack(fill="x", padx=6, pady=(0, 6))

        self._btn_skip = ttk.Button(bf2, text="Skip Module",
                                     command=self._cmd_skip)
        self._btn_skip.pack(side="left", expand=True, fill="x", padx=2)

        self._btn_ack = ttk.Button(bf2, text="Ack Error",
                                    style="Warn.TButton",
                                    command=self._cmd_ack_error)
        self._btn_ack.pack(side="left", expand=True, fill="x", padx=2)

        # progress label + bar
        r_prog = tk.Frame(sc, bg=C.BG)
        r_prog.pack(fill="x", padx=6)
        self._lbl_progress = tk.Label(r_prog, text="Progress: --/--",
                                       bg=C.BG, fg=C.TEXT,
                                       font=("Consolas", 9))
        self._lbl_progress.pack(side="left")
        s = ttk.Style()
        s.configure("scan.Horizontal.TProgressbar",
                     troughcolor=C.PANEL, background=C.ACCENT,
                     bordercolor=C.BORDER, lightcolor=C.ACCENT,
                     darkcolor=C.ACCENT)
        self._progress_bar = ttk.Progressbar(
            r_prog, length=120, mode="determinate",
            style="scan.Horizontal.TProgressbar")
        self._progress_bar.pack(side="right", padx=(4, 0))
        self._lbl_current = tk.Label(sc, text="Current: --",
                                      bg=C.BG, fg=C.TEXT,
                                      font=("Consolas", 9))
        self._lbl_current.pack(padx=6, anchor="w")
        self._lbl_eta = tk.Label(sc, text="ETA: --",
                                  bg=C.BG, fg=C.DIM,
                                  font=("Consolas", 9))
        self._lbl_eta.pack(padx=6, anchor="w")
        self._lbl_dwell_cd = tk.Label(sc, text="",
                                       bg=C.BG, fg=C.GREEN,
                                       font=("Consolas", 9))
        self._lbl_dwell_cd.pack(padx=6, anchor="w", pady=(0, 4))

        # === Direct Control ===
        dc = ttk.LabelFrame(parent, text=" Direct Control ")
        dc.pack(fill="x", pady=(0, 4))

        dcb = tk.Frame(dc, bg=C.BG)
        dcb.pack(fill="x", padx=6, pady=6)

        ttk.Button(dcb, text="Move to Starting Point",
                   command=self._cmd_move_to_module
                   ).pack(fill="x", pady=1)
        ttk.Button(dcb, text="Reset to Beam Center",
                   style="Accent.TButton",
                   command=self._cmd_reset_center
                   ).pack(fill="x", pady=1)

        # === Motor Status  (scrollable) ===
        ms = ttk.LabelFrame(parent, text=" Motor Status ")
        ms.pack(fill="both", expand=True, pady=(0, 4))

        # use a canvas + frame for scrolling if needed
        inner = tk.Frame(ms, bg=C.BG)
        inner.pack(fill="both", expand=True, padx=4, pady=4)

        self._status_labels: Dict[str, tk.Label] = {}
        self._build_motor_block(inner, "X Motor", [
            ("Encoder",  "x_encoder"), ("RBV",  "x_rbv"),
            ("VAL",      "x_val"),     ("MOVN", "x_movn"),
            ("SPMG",     "x_spmg"),    ("VELO", "x_velo"),
            ("ACCL",     "x_accl"),    ("TDIR", "x_tdir"),
            ("MSTA",     "x_msta"),    ("ATHM", "x_athm"),
        ], row=0)

        self._build_motor_block(inner, "Y Motor", [
            ("Encoder",  "y_encoder"), ("RBV",  "y_rbv"),
            ("VAL",      "y_val"),     ("MOVN", "y_movn"),
            ("SPMG",     "y_spmg"),    ("VELO", "y_velo"),
            ("ACCL",     "y_accl"),    ("TDIR", "y_tdir"),
            ("MSTA",     "y_msta"),    ("ATHM", "y_athm"),
        ], row=1)

        # === Position Error ===
        pe = ttk.LabelFrame(parent, text=" Position Check ")
        pe.pack(fill="x", pady=(0, 4))
        pef = tk.Frame(pe, bg=C.BG)
        pef.pack(fill="x", padx=6, pady=4)

        self._lbl_expected = tk.Label(pef, text="Expected: --",
                                       bg=C.BG, fg=C.TEXT,
                                       font=("Consolas", 9))
        self._lbl_expected.pack(anchor="w")
        self._lbl_actual = tk.Label(pef, text="Actual:   --",
                                     bg=C.BG, fg=C.TEXT,
                                     font=("Consolas", 9))
        self._lbl_actual.pack(anchor="w")
        self._lbl_error = tk.Label(pef, text="Diff:     --",
                                    bg=C.BG, fg=C.TEXT,
                                    font=("Consolas", 9, "bold"))
        self._lbl_error.pack(anchor="w")

    def _build_motor_block(self, parent, title: str,
                           fields: List[Tuple[str, str]], row: int):
        frm = tk.Frame(parent, bg=C.BG)
        frm.grid(row=row, column=0, sticky="nsew", padx=2, pady=2)
        parent.rowconfigure(row, weight=1)
        parent.columnconfigure(0, weight=1)

        tk.Label(frm, text=title, bg=C.BG, fg=C.ACCENT,
                 font=("Consolas", 9, "bold")).grid(
            row=0, column=0, columnspan=4, sticky="w", pady=(0, 2))

        half = (len(fields) + 1) // 2
        for i, (label, key) in enumerate(fields):
            c = 0 if i < half else 2
            r = (i % half) + 1
            tk.Label(frm, text=f"{label}:", bg=C.BG, fg=C.DIM,
                     font=("Consolas", 8), anchor="e"
                     ).grid(row=r, column=c, sticky="e", padx=(4, 1))
            lbl = tk.Label(frm, text="--", bg=C.BG, fg=C.TEXT,
                           font=("Consolas", 9), anchor="w", width=12)
            lbl.grid(row=r, column=c + 1, sticky="w", padx=(0, 6))
            self._status_labels[key] = lbl

    # -----------------------------------------------------------------------
    #  Commands
    # -----------------------------------------------------------------------

    def _on_start_selected(self, _event=None):
        name = self._start_var.get()
        for i, m in enumerate(self.engine.path):
            if m.name == name:
                self._selected_start_idx = i
                self._draw_path_preview()
                break

    def _on_path_profile_changed(self, _event=None):
        """Switch between autogen and predefined path profiles."""
        name = self._profile_var.get()
        if name == self._active_profile:
            return
        self._active_profile = name

        if name == self.AUTOGEN:
            # Re-enable LG layers and rebuild with autogen
            self._lg_layers_spin.configure(state="normal")
            self._on_lg_layers_changed(force=True)
            return

        # Predefined path — disable LG layers (not applicable)
        self._lg_layers_spin.configure(state="disabled")
        mod_names = self._profiles.get(name, [])
        mod_by_name = {m.name: m for m in self.all_modules}
        # Resolve names to Module objects, skip unknown
        path_mods = [mod_by_name[n] for n in mod_names if n in mod_by_name]
        if not path_mods:
            self._log(f"Profile '{name}' has no valid modules", level="error")
            return

        self.scan_modules = path_mods
        self._scan_names = {m.name for m in path_mods}
        self.engine = ScanEngine(self.ep, path_mods, self._log)
        # Predefined path: engine.path IS the path (no autogen)
        self.engine.path = path_mods
        self._scan_name_to_idx = {
            m.name: i for i, m in enumerate(self.engine.path)
        }
        self._selected_start_idx = 0

        names = [m.name for m in self.engine.path]
        self._start_combo["values"] = names
        if names:
            self._start_var.set(names[0])
        self._count_entry.configure(to=len(names))
        self._count_var.set(0)

        self._update_canvas_label()
        self._display_greyed = False
        self._draw_modules()
        self._draw_motor_marker()
        self._log(f"Path profile: {name} ({len(path_mods)} modules)")

    def _on_lg_layers_changed(self, force: bool = False):
        """Rebuild scan engine when the user changes LG layers."""
        if self._active_profile != self.AUTOGEN:
            return
        new_layers = self._lg_layers_var.get()
        if new_layers == self._lg_layers and not force:
            return
        self._lg_layers = new_layers
        self.scan_modules = self._filter_scan_modules(new_layers)
        self._scan_names = {m.name for m in self.scan_modules}
        self.engine = ScanEngine(self.ep, self.scan_modules, self._log)
        self._scan_name_to_idx = {
            m.name: i for i, m in enumerate(self.engine.path)
        }
        self._selected_start_idx = 0

        # update start module dropdown and count max
        names = [m.name for m in self.engine.path]
        self._start_combo["values"] = names
        if names:
            self._start_var.set(names[0])
        self._count_entry.configure(to=len(names))
        self._count_var.set(0)

        # update canvas
        self._update_canvas_label()
        self._display_greyed = False
        self._draw_modules()
        self._draw_motor_marker()

        n_pwo4 = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        self._log(f"LG layers: {new_layers} "
                  f"({n_pwo4} PbWO4 + {n_lg} PbGlass = {len(self.scan_modules)})")

    def _cmd_start(self):
        self._on_start_selected()
        # Check boundary limits for the scan segment
        path = self.engine.path
        start = self._selected_start_idx
        count = self._count_var.get()
        end = min(start + count, len(path)) if count > 0 else len(path)
        oob = []
        for i in range(start, end):
            px, py = module_to_ptrans(path[i].x, path[i].y)
            if not ptrans_in_limits(px, py):
                oob.append(path[i].name)
        if oob:
            names = ", ".join(oob[:5])
            if len(oob) > 5:
                names += f", ... ({len(oob)} total)"
            self._log(f"BLOCKED: {len(oob)} modules outside travel limits: "
                      f"{names}", level="error")
            messagebox.showerror("Out of Bounds",
                                 f"{len(oob)} modules outside transporter "
                                 f"travel limits:\n{names}\n\n"
                                 f"Scan not started.")
            return
        self.engine.dwell_time = self._dwell_var.get()
        self.engine.pos_threshold = self._thresh_var.get()
        self.engine.start(self._selected_start_idx, count=count)

    def _cmd_pause(self):
        eng = self.engine
        if eng.state == ScanState.PAUSED:
            eng.resume_scan()
            self._btn_pause.configure(text="Pause")
        elif eng.state in (ScanState.MOVING, ScanState.DWELLING):
            eng.pause_scan()
            self._btn_pause.configure(text="Resume")

    def _cmd_stop(self):
        if self.engine.state != ScanState.IDLE:
            self.engine.stop_scan()
            self._btn_pause.configure(text="Pause")
        else:
            # No scan running — still send stop to halt any manual move,
            # then restore Go so future moves aren't blocked
            epics_stop(self.ep)
            epics_resume(self.ep)
            self._log("Motors stopped")

    def _cmd_skip(self):
        self.engine.skip_module()

    def _cmd_ack_error(self):
        self.engine.acknowledge_error()

    def _cmd_move_to_module(self):
        self._on_start_selected()
        mod = self.engine.path[self._selected_start_idx]
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if not epics_move_to(self.ep, px, py):
            self._log(f"BLOCKED: ptrans({px:.3f}, {py:.3f}) outside travel limits",
                      level="error")

    def _cmd_reset_center(self):
        self._log("Resetting to beam centre "
                  f"ptrans({BEAM_CENTER_X}, {BEAM_CENTER_Y})")
        epics_move_to(self.ep, BEAM_CENTER_X, BEAM_CENTER_Y)

    # -----------------------------------------------------------------------
    #  Logging
    # -----------------------------------------------------------------------

    def _log(self, msg: str, level: str = "info"):
        ts = datetime.now().strftime("%H:%M:%S")
        tag = level.upper().ljust(5)
        line = f"[{ts}] {tag} {msg}"
        self._log_lines.append(line)
        if self._log_file:
            self._log_file.write(line + "\n")
            self._log_file.flush()
        # schedule text widget update on the main thread
        self.root.after_idle(self._append_log, line, level)

    def _append_log(self, line: str, level: str):
        if self._log_text is None:
            return
        self._log_text.configure(state="normal")
        self._log_text.insert("end", line + "\n", level)
        self._log_text.see("end")
        self._log_text.configure(state="disabled")

    # -----------------------------------------------------------------------
    #  Periodic refresh (5 Hz)
    # -----------------------------------------------------------------------

    def _poll(self):
        self._update_status()
        self._update_canvas()
        self._update_scan_info()
        self._update_buttons()
        self.root.after(200, self._poll)

    def _update_status(self):
        """Refresh motor status labels from EPICS readback."""
        for key, lbl in self._status_labels.items():
            val = self.ep.get(key, "--")
            if val == "--" or val is None:
                lbl.configure(text="--", fg=C.DIM)
                continue

            # format value
            if key.endswith("_msta"):
                txt = f"0x{int(val):X}"
            elif key.endswith("_spmg"):
                txt = f"{SPMG_LABELS.get(int(val), '?')}({int(val)})"
            elif key.endswith("_movn") or key.endswith("_athm") or key.endswith("_tdir"):
                txt = str(int(val))
            elif isinstance(val, float):
                txt = f"{val:.3f}"
            else:
                txt = str(val)

            # colour moving indicators
            fg = C.TEXT
            if key.endswith("_movn") and int(val) == 1:
                fg = C.YELLOW
            elif key.endswith("_spmg"):
                spmg_val = int(val)
                if spmg_val == SPMG.STOP:
                    fg = C.RED
                elif spmg_val == SPMG.PAUSE:
                    fg = C.ORANGE
                else:  # MOVE or GO
                    fg = C.GREEN

            lbl.configure(text=txt, fg=fg)

        # position check
        eng = self.engine
        rx = self.ep.get("x_rbv", 0.0)
        ry = self.ep.get("y_rbv", 0.0)
        self._lbl_actual.configure(
            text=f"Actual:   ({rx:.3f}, {ry:.3f})")

        # Determine expected position: current scan target or selected module
        mod = eng.current_module if eng.state != ScanState.IDLE else None
        if mod is None and eng.path and 0 <= self._selected_start_idx < len(eng.path):
            mod = eng.path[self._selected_start_idx]

        scanning = eng.state in (ScanState.MOVING, ScanState.DWELLING,
                                  ScanState.PAUSED, ScanState.ERROR)
        if mod:
            px, py = module_to_ptrans(mod.x, mod.y)
            err = math.sqrt((rx - px)**2 + (ry - py)**2)
            self._lbl_expected.configure(
                text=f"Expected: {mod.name} ({px:.3f}, {py:.3f})")
            if scanning:
                err_fg = C.RED if err > eng.pos_threshold else C.GREEN
                self._lbl_error.configure(
                    text=f"Diff:     {err:.3f} mm", fg=err_fg)
            else:
                self._lbl_error.configure(
                    text=f"Diff:     {err:.3f} mm (not scanning)", fg=C.DIM)
        else:
            self._lbl_expected.configure(text="Expected: --")
            self._lbl_error.configure(text="Diff:     -- (not scanning)",
                                      fg=C.DIM)

    def _update_scan_info(self):
        eng = self.engine

        # state badge
        state_colours = {
            ScanState.IDLE:      C.DIM,
            ScanState.MOVING:    C.YELLOW,
            ScanState.DWELLING:  C.GREEN,
            ScanState.PAUSED:    C.ORANGE,
            ScanState.ERROR:     C.RED,
            ScanState.COMPLETED: C.ACCENT,
        }
        self._lbl_state.configure(
            text=eng.state,
            fg=state_colours.get(eng.state, C.DIM))

        # progress
        done = len(eng.completed)
        start = getattr(eng, '_start_idx', 0)
        end = getattr(eng, '_end_idx', len(eng.path))
        total = end - start
        self._lbl_progress.configure(text=f"Progress: {done}/{total}")
        self._progress_bar["maximum"] = max(total, 1)
        self._progress_bar["value"] = done

        mod = eng.current_module
        if mod:
            self._lbl_current.configure(text=f"Current:  {mod.name}")
        else:
            self._lbl_current.configure(text="Current:  --")

        # ETA
        if eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED):
            eta = eng.eta_seconds
            h, rem = divmod(int(eta), 3600)
            m, s = divmod(rem, 60)
            self._lbl_eta.configure(text=f"ETA:      {h}h {m:02d}m {s:02d}s")
        elif eng.state == ScanState.IDLE and eng.path:
            # Pre-scan estimate from path geometry and motor velocities
            vx = self.ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
            vy = self.ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
            eta = estimate_scan_time(
                eng.path, self._selected_start_idx,
                self._count_var.get(), self._dwell_var.get(), vx, vy)
            if eta > 0:
                h, rem = divmod(int(eta), 3600)
                m, s = divmod(rem, 60)
                self._lbl_eta.configure(
                    text=f"ETA:      ~{h}h {m:02d}m {s:02d}s")
            else:
                self._lbl_eta.configure(text="ETA:      --")
        else:
            self._lbl_eta.configure(text="ETA:      --")

        # dwell countdown
        if eng.state == ScanState.DWELLING:
            self._lbl_dwell_cd.configure(
                text=f"Dwell:    {eng.dwell_remaining:.1f}s remaining")
        else:
            self._lbl_dwell_cd.configure(text="")

    def _update_buttons(self):
        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING,
                                 ScanState.PAUSED, ScanState.ERROR)
        # Only reconfigure when state actually changes to avoid
        # clearing the ttk "active" (hover) flag every poll cycle.
        new_states = {
            "_btn_start":    "disabled" if running else "normal",
            "_btn_pause":    "normal" if running else "disabled",
            "_btn_stop":     "normal",
            "_btn_skip":     "normal" if eng.state == ScanState.DWELLING
                             else "disabled",
            "_btn_ack":      "normal" if eng.state == ScanState.ERROR
                             else "disabled",
            "_start_combo":  "readonly" if not running else "disabled",
            "_count_entry":  "normal" if not running else "disabled",
            "_profile_combo": "readonly" if not running else "disabled",
            "_lg_layers_spin": "normal" if (not running and
                                self._active_profile == self.AUTOGEN)
                               else "disabled",
        }
        for attr, desired in new_states.items():
            widget = getattr(self, attr)
            if str(widget.cget("state")) != desired:
                widget.configure(state=desired)


# ============================================================================
#  MAIN
# ============================================================================

PATHS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "paths.json")


def _load_path_profiles(path: str) -> Dict[str, List[str]]:
    if os.path.exists(path):
        import json
        with open(path) as f:
            return json.load(f)
    return {}


def main():
    parser = argparse.ArgumentParser(
        description="HyCal Snake Scan -- module scanner")
    parser.add_argument("--expert", action="store_true",
                        help="Expert operator mode: full control via real EPICS")
    parser.add_argument("--observer", action="store_true",
                        help="Observer mode: read-only, no PV writes")
    parser.add_argument("--database", default=DEFAULT_DB_PATH,
                        help="Path to hycal_modules.json")
    parser.add_argument("--paths", default=PATHS_FILE,
                        help="Path to paths.json for predefined scan paths")
    args = parser.parse_args()

    # Load modules from database
    all_modules = load_modules(args.database)
    by_type: Dict[str, int] = {}
    for m in all_modules:
        by_type[m.mod_type] = by_type.get(m.mod_type, 0) + 1
    print(f"Loaded {len(all_modules)} modules from {args.database}")
    for t, n in sorted(by_type.items()):
        print(f"  {t}: {n}")

    # Load predefined paths
    profiles = _load_path_profiles(args.paths)
    if profiles:
        print(f"Loaded {len(profiles)} path profiles from {args.paths}")

    observer = args.observer
    simulation = not args.expert and not observer

    if observer:
        ep = ObserverEPICS()
    elif simulation:
        ep = SimulatedEPICS()
    else:
        ep = RealEPICS(writable=True)

    n_ok, n_total = ep.connect()
    if not simulation:
        print(f"EPICS: connected {n_ok}/{n_total} PVs")
        if n_ok < n_total:
            disconnected = ep.disconnected_pvs()
            for pv in disconnected:
                print(f"  NOT connected: {pv}")
        if n_ok < n_total * 0.5:
            print("WARNING: many PVs not connected -- check IOC / network")

    root = tk.Tk()
    SnakeScanGUI(root, ep, simulation, all_modules, profiles,
                 observer=observer)
    root.mainloop()


if __name__ == "__main__":
    main()
