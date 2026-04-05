#!/usr/bin/env python3
"""
HyCal Snake Scan -- Module Scanner (PyQt6)
==========================================
PyQt6 GUI that drives the HyCal transporter in a snake pattern so the
beam centres on each scanned module, dwells for a configurable time,
then advances to the next module.

Usage
-----
    python hycal_snake_scan.py                          # simulation
    python hycal_snake_scan.py --expert                  # expert operator
    python hycal_snake_scan.py --observer                # read-only monitor

Coordinate system
-----------------
    ptrans_x, ptrans_y = (-126.75, 10.11)  -->  beam at HyCal centre (0,0)
    ptrans_x = BEAM_CENTER_X + module_x
    ptrans_y = BEAM_CENTER_Y - module_y

Writable PVs (the ONLY PVs this tool writes to):
    ptrans_x.VAL / ptrans_y.VAL    -- absolute set-point
    ptrans_x.SPMG / ptrans_y.SPMG  -- motor mode  Stop(0) Pause(1) Move(2) Go(3)

Requirements
------------
    Python 3.8+, PyQt6
    pyepics  (only for --expert / --observer mode)
"""

from __future__ import annotations

import argparse
import html as html_mod
import json
import math
import os
import random
import sys
import threading
import time
from datetime import datetime
from enum import IntEnum
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QGroupBox, QPushButton, QLabel, QComboBox, QSpinBox,
    QDoubleSpinBox, QTextEdit, QProgressBar, QMessageBox, QSplitter,
    QSizePolicy, QFrame, QDialog,
)
from PyQt6.QtCore import Qt, QRectF, QPointF, QTimer, pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QPen, QFont

from scan_utils import (
    C, Module, load_modules, module_to_ptrans, ptrans_to_module,
    ptrans_in_limits, filter_scan_modules, DARK_QSS,
    BEAM_CENTER_X, BEAM_CENTER_Y, DEFAULT_DB_PATH,
    PTRANS_X_MIN, PTRANS_X_MAX, PTRANS_Y_MIN, PTRANS_Y_MAX,
)


# ============================================================================
#  SCAN-SPECIFIC CONSTANTS
# ============================================================================

DEFAULT_DWELL = 120.0         # seconds
DEFAULT_POS_THRESHOLD = 0.5   # mm
DEFAULT_BEAM_THRESHOLD = 0.3  # nA
MOVE_TIMEOUT = 300.0          # seconds per single move
MAX_LG_LAYERS = 2
DEFAULT_VELO_X = 50.0         # mm/s
DEFAULT_VELO_Y = 5.0          # mm/s

SPMG_LABELS = {0: "Stop", 1: "Pause", 2: "Move", 3: "Go"}


class SPMG(IntEnum):
    STOP = 0; PAUSE = 1; MOVE = 2; GO = 3


class PV:
    X_VAL  = "ptrans_x.VAL";   Y_VAL  = "ptrans_y.VAL"
    X_SPMG = "ptrans_x.SPMG";  Y_SPMG = "ptrans_y.SPMG"
    X_ENCODER = "hallb_ptrans_x_encoder"
    Y_ENCODER = "hallb_ptrans_y1_encoder"
    BEAM_CUR = "hallb_IPM2C21A_CUR"
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


# -- path building -----------------------------------------------------------

def _snake_sector(modules, going_right, top_to_bottom):
    if not modules:
        return [], going_right
    by_y = sorted(modules, key=lambda m: -m.y)
    rows = []
    cur_row = [by_y[0]]
    for m in by_y[1:]:
        if abs(m.y - cur_row[0].y) < 0.5:
            cur_row.append(m)
        else:
            rows.append(cur_row)
            cur_row = [m]
    rows.append(cur_row)
    if not top_to_bottom:
        rows.reverse()
    path = []
    for row in rows:
        row.sort(key=lambda m: m.x, reverse=not going_right)
        path.extend(row)
        going_right = not going_right
    return path, going_right


def build_scan_path(scan_modules):
    if not scan_modules:
        return [], 0
    center = [m for m in scan_modules if m.sector == "Center"]
    lg_sectors = {}
    for m in scan_modules:
        if m.sector in ("Top", "Right", "Bottom", "Left"):
            lg_sectors.setdefault(m.sector, []).append(m)
    path = []
    going = True
    if center:
        seg, going = _snake_sector(center, going, True)
        path.extend(seg)
    remaining = dict(lg_sectors)
    while remaining:
        last = path[-1] if path else None
        best_name = None
        best_dist = float('inf')
        best_cr = True
        best_cd = True
        for name, mods in remaining.items():
            xs = [m.x for m in mods]
            ys = [m.y for m in mods]
            for cx, cy, gr, gd in [
                (min(xs), max(ys), True, True), (max(xs), max(ys), False, True),
                (min(xs), min(ys), True, False), (max(xs), min(ys), False, False),
            ]:
                d = (abs(last.x - cx) + abs(last.y - cy)) if last else 0
                if d < best_dist:
                    best_dist, best_name, best_cr, best_cd = d, name, gr, gd
        seg, going = _snake_sector(remaining.pop(best_name), best_cr, best_cd)
        path.extend(seg)
    return path, len(scan_modules) - len(path)


def estimate_scan_time(path, start, count, dwell, vx=DEFAULT_VELO_X, vy=DEFAULT_VELO_Y):
    end = min(start + count, len(path)) if count > 0 else len(path)
    if end <= start:
        return 0.0
    total_move = 0.0
    for i in range(start, end - 1):
        dx = abs(path[i + 1].x - path[i].x)
        dy = abs(path[i + 1].y - path[i].y)
        total_move += max(dx / vx if vx > 0 else 0, dy / vy if vy > 0 else 0)
    return total_move + (end - start) * dwell


# ============================================================================
#  EPICS INTERFACES
# ============================================================================

_PV_MAP = [
    ("x_val", PV.X_VAL), ("y_val", PV.Y_VAL),
    ("x_spmg", PV.X_SPMG), ("y_spmg", PV.Y_SPMG),
    ("x_encoder", PV.X_ENCODER), ("y_encoder", PV.Y_ENCODER),
    ("beam_cur", PV.BEAM_CUR),
    ("x_rbv", PV.X_RBV), ("y_rbv", PV.Y_RBV),
    ("x_movn", PV.X_MOVN), ("y_movn", PV.Y_MOVN),
    ("x_velo", PV.X_VELO), ("y_velo", PV.Y_VELO),
    ("x_accl", PV.X_ACCL), ("y_accl", PV.Y_ACCL),
    ("x_tdir", PV.X_TDIR), ("y_tdir", PV.Y_TDIR),
    ("x_msta", PV.X_MSTA), ("y_msta", PV.Y_MSTA),
    ("x_athm", PV.X_ATHM), ("y_athm", PV.Y_ATHM),
    ("x_prec", PV.X_PREC), ("y_prec", PV.Y_PREC),
    ("x_bvel", PV.X_BVEL), ("y_bvel", PV.Y_BVEL),
    ("x_bacc", PV.X_BACC), ("y_bacc", PV.Y_BACC),
    ("x_vbas", PV.X_VBAS), ("y_vbas", PV.Y_VBAS),
    ("x_bdst", PV.X_BDST), ("y_bdst", PV.Y_BDST),
    ("x_frac", PV.X_FRAC), ("y_frac", PV.Y_FRAC),
]


class RealEPICS:
    def __init__(self, writable=False):
        import epics as _epics
        self._epics = _epics
        self._pvs = {}
        self._writable = writable

    def connect(self):
        for key, pvname in _PV_MAP:
            self._pvs[key] = self._epics.PV(pvname, connection_timeout=5.0)
        time.sleep(2.0)
        n = sum(1 for p in self._pvs.values() if p.connected)
        self._all_connected = (n == len(self._pvs))
        return n, len(self._pvs)

    def disconnected_pvs(self):
        return [pv for k, pv in _PV_MAP if k in self._pvs and not self._pvs[k].connected]

    def get(self, key, default=None):
        pv = self._pvs.get(key)
        if pv and pv.connected:
            v = pv.get()
            return v if v is not None else default
        return default

    def put(self, key, value):
        if not self._writable or not self._all_connected:
            return False
        pv = self._pvs.get(key)
        if pv and pv.connected:
            pv.put(value)
            return True
        return False

    def stop(self):
        for key in ("x_spmg", "y_spmg"):
            pv = self._pvs.get(key)
            if pv and pv.connected:
                pv.put(int(SPMG.STOP))


class SimulatedEPICS:
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
        self._x_speed = 0.5
        self._y_speed = 0.2
        self._moving = False

    def connect(self):
        return (0, 0)

    def disconnected_pvs(self):
        return []

    def get(self, key, default=None):
        with self._lock:
            return {
                "x_encoder": round(self._x + random.gauss(0, 0.002), 4),
                "y_encoder": round(self._y + random.gauss(0, 0.002), 4),
                "x_rbv": round(self._x, 3), "y_rbv": round(self._y, 3),
                "x_val": round(self._tx, 3), "y_val": round(self._ty, 3),
                "x_movn": self._x_movn, "y_movn": self._y_movn,
                "x_spmg": self._x_spmg, "y_spmg": self._y_spmg,
                "x_velo": self._x_speed, "y_velo": self._y_speed,
                "x_accl": 0.2, "y_accl": 1.0,
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
                "beam_cur": 50.0,
            }.get(key, default)

    def put(self, key, value):
        with self._lock:
            if key == "x_val":      self._tx = float(value)
            elif key == "y_val":    self._ty = float(value)
            elif key == "x_spmg":   self._x_spmg = int(value)
            elif key == "y_spmg":   self._y_spmg = int(value)
            else: return False
            self._evaluate_motion()
        return True

    def stop(self):
        self.put("x_spmg", int(SPMG.STOP))
        self.put("y_spmg", int(SPMG.STOP))

    def _evaluate_motion(self):
        if self._x_spmg == SPMG.STOP or self._y_spmg == SPMG.STOP:
            self._moving = False; self._x_movn = 0; self._y_movn = 0
        elif self._x_spmg == SPMG.PAUSE or self._y_spmg == SPMG.PAUSE:
            self._moving = False; self._x_movn = 0; self._y_movn = 0
        elif self._x_spmg == SPMG.GO and self._y_spmg == SPMG.GO:
            if not self._moving:
                self._moving = True
                threading.Thread(target=self._run_move, daemon=True).start()

    def _run_move(self):
        dt = 0.02
        while True:
            with self._lock:
                if not self._moving:
                    self._x_movn = 0; self._y_movn = 0; return
                dx, dy = self._tx - self._x, self._ty - self._y
                if abs(dx) > 0.001:
                    self._x += math.copysign(min(self._x_speed * dt, abs(dx)), dx)
                if abs(dy) > 0.001:
                    self._y += math.copysign(min(self._y_speed * dt, abs(dy)), dy)
                self._x_movn = 1 if abs(self._tx - self._x) > 0.001 else 0
                self._y_movn = 1 if abs(self._ty - self._y) > 0.001 else 0
                if self._x_movn == 0 and self._y_movn == 0:
                    self._x, self._y = self._tx, self._ty
                    self._moving = False; return
            time.sleep(dt)


class ObserverEPICS:
    def __init__(self):
        self._real = RealEPICS()
    def connect(self):           return self._real.connect()
    def disconnected_pvs(self):  return self._real.disconnected_pvs()
    def get(self, key, default=None): return self._real.get(key, default)
    def put(self, key, value):   return False
    def stop(self):              pass


def epics_move_to(ep, x, y):
    if not ptrans_in_limits(x, y): return False
    ep.put("x_val", x); ep.put("y_val", y)
    ep.put("x_spmg", int(SPMG.GO)); ep.put("y_spmg", int(SPMG.GO))
    return True

def epics_stop(ep):     ep.stop()
def epics_pause(ep):    ep.put("x_spmg", int(SPMG.PAUSE)); ep.put("y_spmg", int(SPMG.PAUSE))
def epics_resume(ep):   ep.put("x_spmg", int(SPMG.GO));    ep.put("y_spmg", int(SPMG.GO))
def epics_is_moving(ep): return bool(ep.get("x_movn", 0)) or bool(ep.get("y_movn", 0))
def epics_read_rbv(ep): return (ep.get("x_rbv", 0.0), ep.get("y_rbv", 0.0))


# ============================================================================
#  SCAN ENGINE
# ============================================================================

class ScanState:
    IDLE = "IDLE"; MOVING = "MOVING"; DWELLING = "DWELLING"
    PAUSED = "PAUSED"; ERROR = "ERROR"; COMPLETED = "COMPLETED"


class ScanEngine:
    def __init__(self, epics, modules, log_fn):
        self.ep = epics
        self.all_modules = modules
        self.path, n_unopt = build_scan_path(modules)
        self.log = log_fn
        if n_unopt:
            self.log(f"WARNING: {n_unopt} modules with unoptimized path", level="warn")

        self.dwell_time = DEFAULT_DWELL
        self.pos_threshold = DEFAULT_POS_THRESHOLD
        self.beam_threshold = DEFAULT_BEAM_THRESHOLD

        self.state = ScanState.IDLE
        self.current_idx = 0
        self.dwell_remaining = 0.0
        self.completed = set()
        self.error_modules = set()
        self.beam_tripped = False

        self._thread = None
        self._stop = threading.Event()
        self._skip = threading.Event()
        self._paused = False
        self._ack_error = threading.Event()

    @property
    def current_module(self):
        if 0 <= self.current_idx < len(self.path):
            return self.path[self.current_idx]
        return None

    @property
    def eta_seconds(self):
        if self.state == ScanState.IDLE: return 0.0
        end = getattr(self, '_end_idx', len(self.path))
        nxt = self.current_idx + 1
        if nxt >= end: return self.dwell_remaining
        return estimate_scan_time(self.path, nxt, end - nxt, self.dwell_time) + self.dwell_remaining

    def start(self, start_idx=0, count=0):
        if self._thread and self._thread.is_alive(): return
        self._stop.clear(); self._skip.clear(); self._paused = False
        self.current_idx = start_idx
        self._start_idx = start_idx
        self._end_idx = min(start_idx + count, len(self.path)) if count > 0 else len(self.path)
        self.completed.clear(); self.error_modules.clear()
        self.state = ScanState.MOVING
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def resume_scan(self):
        if self.state == ScanState.PAUSED:
            self._paused = False; epics_resume(self.ep)
            self.state = ScanState.MOVING; self.log("Scan resumed")

    def pause_scan(self):
        if self.state in (ScanState.MOVING, ScanState.DWELLING):
            self._paused = True; epics_pause(self.ep)
            self.state = ScanState.PAUSED; self.log("Scan paused", level="warn")

    def stop_scan(self):
        self._stop.set(); self._skip.set(); self._ack_error.set()
        self._paused = False; epics_stop(self.ep)
        self.log("Scan stopped", level="warn")

    def skip_module(self):       self._skip.set()
    def acknowledge_error(self): self._ack_error.set()

    def _run(self):
        self.log(f"Scan started from {self.path[self.current_idx].name}, "
                 f"dwell {self.dwell_time:.0f}s, {self._end_idx - self.current_idx} modules")
        try:
            for i in range(self.current_idx, self._end_idx):
                if self._stop.is_set(): break
                self.current_idx = i
                mod = self.path[i]
                px, py = module_to_ptrans(mod.x, mod.y)
                self.state = ScanState.MOVING
                self.log(f"[{i+1}/{len(self.path)}] Moving to {mod.name}  ptrans({px:.3f}, {py:.3f})")
                if not epics_move_to(self.ep, px, py):
                    self.log(f"SKIPPED {mod.name}: outside limits", level="warn"); continue
                if not self._wait_move_done(): break
                rbv_x, rbv_y = epics_read_rbv(self.ep)
                err = math.sqrt((rbv_x - px)**2 + (rbv_y - py)**2)
                if err > self.pos_threshold:
                    self.error_modules.add(i); self.state = ScanState.ERROR
                    self.log(f"POSITION ERROR at {mod.name}: {err:.3f} mm", level="error")
                    self._ack_error.clear(); self._ack_error.wait()
                    if self._stop.is_set(): break
                self.state = ScanState.DWELLING
                self.dwell_remaining = self.dwell_time
                self.log(f"Dwelling at {mod.name} for {self.dwell_time:.0f}s")
                if self._wait_dwell() == "stop": break
                self.completed.add(i); self.dwell_remaining = 0.0
        finally:
            if self._stop.is_set():
                self.state = ScanState.IDLE; self.log("Scan stopped by user")
            elif self.current_idx >= self._end_idx - 1:
                self.state = ScanState.COMPLETED; self.log("Scan COMPLETE!", level="warn")
            else:
                self.state = ScanState.IDLE

    def _wait_move_done(self):
        t0 = time.time()
        while not self._stop.is_set():
            while self._paused and not self._stop.is_set():
                self.state = ScanState.PAUSED; time.sleep(0.1)
            if self._stop.is_set(): return False
            self.state = ScanState.MOVING
            if not epics_is_moving(self.ep): return True
            if time.time() - t0 > MOVE_TIMEOUT:
                self.log(f"MOVE TIMEOUT after {MOVE_TIMEOUT:.0f}s", level="error"); return False
            time.sleep(0.1)
        return False

    def _wait_dwell(self):
        end = time.time() + self.dwell_time
        while time.time() < end:
            if self._stop.is_set(): return "stop"
            if self._skip.is_set():
                self._skip.clear(); self.log("Module skipped"); return "skip"
            if self._paused:
                self.state = ScanState.PAUSED
                while self._paused and not self._stop.is_set(): time.sleep(0.1)
                if self._stop.is_set(): return "stop"
                end = time.time() + self.dwell_time
                self.log("Dwell restarted after resume")
            # -- beam trip pause (dwell only) --
            if self.beam_threshold > 0:
                bc = self.ep.get("beam_cur", None)
                if bc is not None and bc < self.beam_threshold:
                    self.beam_tripped = True
                    self.log(f"BEAM TRIP: {bc:.3f} nA < {self.beam_threshold} nA -- dwell paused", level="warn")
                    while not self._stop.is_set() and not self._skip.is_set():
                        bc2 = self.ep.get("beam_cur", 0.0)
                        if bc2 is not None and bc2 >= self.beam_threshold: break
                        time.sleep(0.5)
                    self.beam_tripped = False
                    if self._stop.is_set(): return "stop"
                    if self._skip.is_set():
                        self._skip.clear(); self.log("Module skipped"); return "skip"
                    bc2 = self.ep.get("beam_cur", 0.0)
                    self.log(f"BEAM RECOVERED: {bc2:.3f} nA -- restarting full dwell", level="warn")
                    end = time.time() + self.dwell_time
            self.state = ScanState.DWELLING
            self.dwell_remaining = max(0.0, end - time.time())
            time.sleep(0.1)
        return "done"


# ============================================================================
#  MAP WIDGET
# ============================================================================

class HyCalScanMapWidget(QWidget):
    moduleClicked = pyqtSignal(str)
    PAD = 8
    SHRINK = 0.90

    def __init__(self, all_modules, parent=None):
        super().__init__(parent)
        self._drawn = [m for m in all_modules if m.mod_type != "LMS"]
        self._mod_by_name = {m.name: m for m in all_modules}
        self._colors = {}
        self._path_line = []
        self._dash_line = []
        self._marker_hx = self._marker_hy = None
        self._highlight = None
        self._hover_name = None
        self._lim_hx_min = PTRANS_X_MIN - BEAM_CENTER_X
        self._lim_hx_max = PTRANS_X_MAX - BEAM_CENTER_X
        self._lim_hy_min = BEAM_CENTER_Y - PTRANS_Y_MAX
        self._lim_hy_max = BEAM_CENTER_Y - PTRANS_Y_MIN
        self._scale = 1.0
        self._ox = self._oy = 0.0
        self._x_min = self._y_max = 0.0
        self._rects = {}
        self.setMouseTracking(True)
        self.setMinimumSize(400, 400)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def setModuleColors(self, c):  self._colors = c; self.update()
    def setPathPreview(self, p):   self._path_line = p; self.update()
    def setDashPreview(self, p):   self._dash_line = p; self.update()
    def setHighlight(self, n):     self._highlight = n; self.update()

    def setMarkerPosition(self, hx, hy):
        self._marker_hx = hx; self._marker_hy = hy; self.update()

    def modCenter(self, m):
        cx = self._ox + (m.x - self._x_min) * self._scale
        cy = self._oy + (self._y_max - m.y) * self._scale
        return QPointF(cx, cy)

    def _updateLayout(self):
        if not self._drawn: return
        x_min = min(m.x - m.sx / 2 for m in self._drawn)
        x_max = max(m.x + m.sx / 2 for m in self._drawn)
        y_min = min(m.y - m.sy / 2 for m in self._drawn)
        y_max = max(m.y + m.sy / 2 for m in self._drawn)
        w, h = self.width(), self.height()
        uw, uh = w - 2 * self.PAD, h - 2 * self.PAD
        self._scale = min(uw / max(x_max - x_min, 1), uh / max(y_max - y_min, 1))
        dw = (x_max - x_min) * self._scale
        dh = (y_max - y_min) * self._scale
        self._ox = self.PAD + (uw - dw) / 2
        self._oy = self.PAD + (uh - dh) / 2
        self._x_min = x_min; self._y_max = y_max
        self._rects.clear()
        for m in self._drawn:
            cx = self._ox + (m.x - self._x_min) * self._scale
            cy = self._oy + (self._y_max - m.y) * self._scale
            hw = m.sx * self._scale * self.SHRINK / 2
            hh = m.sy * self._scale * self.SHRINK / 2
            self._rects[m.name] = QRectF(cx - hw, cy - hh, 2 * hw, 2 * hh)

    def paintEvent(self, event):
        self._updateLayout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        p.fillRect(self.rect(), QColor("#0a0e14"))
        for m in self._drawn:
            r = self._rects.get(m.name)
            if r: p.fillRect(r, QColor(self._colors.get(m.name, C.MOD_EXCLUDED)))
        # limit box
        bx0 = self._ox + (self._lim_hx_min - self._x_min) * self._scale
        by0 = self._oy + (self._y_max - self._lim_hy_max) * self._scale
        bx1 = self._ox + (self._lim_hx_max - self._x_min) * self._scale
        by1 = self._oy + (self._y_max - self._lim_hy_min) * self._scale
        p.setPen(QPen(QColor(C.RED), 1, Qt.PenStyle.DashLine))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(QRectF(bx0, by0, bx1 - bx0, by1 - by0))
        # path lines
        for pts, style, w in [(self._path_line, Qt.PenStyle.SolidLine, 1.5),
                               (self._dash_line, Qt.PenStyle.DashLine, 1.0)]:
            if len(pts) >= 2:
                p.setPen(QPen(QColor(C.PATH_LINE), w, style))
                for i in range(len(pts) - 1):
                    p.drawLine(pts[i], pts[i + 1])
        # highlight
        if self._highlight and self._highlight in self._rects:
            p.setPen(QPen(QColor(C.ACCENT), 2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._highlight])
        # motor crosshair
        if self._marker_hx is not None:
            cx = self._ox + (self._marker_hx - self._x_min) * self._scale
            cy = self._oy + (self._y_max - self._marker_hy) * self._scale
            p.setPen(QPen(QColor(C.RED), 1.5))
            p.drawLine(QPointF(cx - 5, cy), QPointF(cx + 5, cy))
            p.drawLine(QPointF(cx, cy - 5), QPointF(cx, cy + 5))
        p.end()

    def mousePressEvent(self, e):
        if e.button() != Qt.MouseButton.LeftButton: return
        pos = e.position()
        for name in reversed(list(self._rects)):
            if self._rects[name].contains(pos):
                self.moduleClicked.emit(name); return

    def mouseMoveEvent(self, e):
        pos = e.position()
        found = None
        for name in reversed(list(self._rects)):
            if self._rects[name].contains(pos): found = name; break
        if found != self._hover_name:
            self._hover_name = found
            if found:
                m = self._mod_by_name.get(found)
                if m:
                    px, py = module_to_ptrans(m.x, m.y)
                    self.setToolTip(f"{m.name} ({m.mod_type})\nHyCal({m.x:.1f}, {m.y:.1f})\nptrans({px:.1f}, {py:.1f})")
            else:
                self.setToolTip("")

    def resizeEvent(self, e):
        self._updateLayout(); super().resizeEvent(e)

    def sizeHint(self):
        from PyQt6.QtCore import QSize
        return QSize(680, 680)


# ============================================================================
#  MODULE INFO DIALOG
# ============================================================================

class ModuleInfoDialog(QDialog):
    """Pop-up showing module details with a Move To button."""

    def __init__(self, mod: Module, ep, log_fn, parent=None):
        super().__init__(parent)
        self._mod = mod
        self._ep = ep
        self._log = log_fn
        self.setWindowTitle(f"Module {mod.name}")
        self.setStyleSheet(DARK_QSS)
        self.setFixedWidth(320)

        lo = QVBoxLayout(self)

        # -- info grid --
        px, py = module_to_ptrans(mod.x, mod.y)
        in_limits = ptrans_in_limits(px, py)

        grid = QGridLayout()
        grid.setSpacing(4)
        rows = [
            ("Name",     mod.name),
            ("Type",     mod.mod_type),
            ("Sector",   mod.sector or "--"),
            ("Row/Col",  f"{mod.row} / {mod.col}" if mod.row else "--"),
            ("Size",     f"{mod.sx:.2f} x {mod.sy:.2f} mm"),
            ("HyCal",    f"({mod.x:.2f}, {mod.y:.2f}) mm"),
            ("Ptrans",   f"({px:.2f}, {py:.2f}) mm"),
            ("In limits", "Yes" if in_limits else "No"),
        ]
        for r, (label, value) in enumerate(rows):
            lk = QLabel(f"{label}:")
            lk.setStyleSheet(f"color: {C.DIM}; font: 9pt 'Consolas';")
            lk.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(lk, r, 0)
            lv = QLabel(str(value))
            if label == "In limits" and not in_limits:
                lv.setStyleSheet(f"color: {C.RED}; font: bold 9pt 'Consolas';")
            grid.addWidget(lv, r, 1)
        lo.addLayout(grid)

        lo.addSpacing(8)

        # -- buttons --
        btn_row = QHBoxLayout()
        btn_move = QPushButton(f"Move To {mod.name}")
        btn_move.setProperty("cssClass", "accent")
        btn_move.setEnabled(in_limits)
        btn_move.clicked.connect(self._doMove)
        btn_row.addWidget(btn_move)
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_row.addWidget(btn_close)
        lo.addLayout(btn_row)

    def _doMove(self):
        mod = self._mod
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if not epics_move_to(self._ep, px, py):
            self._log(f"BLOCKED: ptrans({px:.3f}, {py:.3f}) outside limits", level="error")
        self.close()


# ============================================================================
#  GUI
# ============================================================================

class SnakeScanWindow(QMainWindow):
    _logSignal = pyqtSignal(str, str)
    AUTOGEN = "(autogen)"
    NONE = "(none)"

    def __init__(self, epics, simulation, all_modules, profiles=None, observer=False):
        super().__init__()
        self.ep = epics
        self.simulation = simulation
        self.observer = observer
        self.all_modules = all_modules
        self._profiles = profiles or {}
        self._active_profile = self.AUTOGEN
        self._lg_layers = 0

        glass = [m for m in all_modules if m.mod_type == "PbGlass"]
        self._lg_sx = glass[0].sx if glass else 38.15
        self._lg_sy = glass[0].sy if glass else 38.15

        self._mod_by_name = {m.name: m for m in all_modules}
        self._log_lines = []

        log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
        os.makedirs(log_dir, exist_ok=True)
        self._log_file = open(os.path.join(log_dir,
            datetime.now().strftime("snake_scan_%Y%m%d_%H%M%S.log")), "w")

        self.scan_modules = filter_scan_modules(all_modules, 0, self._lg_sx, self._lg_sy)
        self.engine = ScanEngine(epics, self.scan_modules, self._log)
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(self.engine.path)}
        self._scan_names = {m.name for m in self.scan_modules}
        self._selected_start_idx = 0
        self._selected_mod_name = None
        self._status_labels = {}

        self._logSignal.connect(self._appendLog)
        self._buildUI()

        if self.observer:
            self._disableControls()
        if not self.simulation and not self.observer:
            disc = self.ep.disconnected_pvs()
            if disc:
                self._disableControls()
                QMessageBox.critical(self, "PV Connection Error",
                    "Not connected:\n" + "\n".join(f"  {p}" for p in disc))

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._poll)
        self._timer.start(200)

    # -----------------------------------------------------------------------
    #  Layout
    #
    #  [TOP BAR: title | mode | ===BEAM=== | state]
    #  [LEFT column              | RIGHT column (full height)       ]
    #  [  map + legend           |   Scan Control                   ]
    #  [  event log              |   Direct Control                 ]
    #  [                         |   Motor Status                   ]
    #  [                         |   Position Check                 ]
    # -----------------------------------------------------------------------

    def _buildUI(self):
        if self.observer:       suffix = "  [OBSERVER]"
        elif self.simulation:   suffix = "  [SIMULATION]"
        else:                   suffix = "  [EXPERT OPERATOR]"
        self.setWindowTitle("HyCal Snake Scan" + suffix)
        self.setStyleSheet(DARK_QSS)
        self.resize(1200, 820)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── top bar ──────────────────────────────────────────────────────
        top = QWidget()
        top.setFixedHeight(48)
        top.setStyleSheet("background: #0d1520;")
        tl = QHBoxLayout(top)
        tl.setContentsMargins(12, 0, 12, 0)

        lbl = QLabel("HYCAL SNAKE SCAN")
        lbl.setStyleSheet(f"color: {C.GREEN}; font: bold 13pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl)

        if self.observer:       mt, mf = "OBSERVER", C.ORANGE
        elif self.simulation:   mt, mf = "SIMULATION", C.YELLOW
        else:                   mt, mf = "EXPERT", C.GREEN
        lbl_mode = QLabel(mt)
        lbl_mode.setStyleSheet(f"color: {mf}; font: bold 9pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl_mode)

        tl.addSpacing(16)

        # ── prominent beam current ──
        beam_frame = QFrame()
        beam_frame.setStyleSheet(
            "QFrame { background: #161b22; border: 1px solid #30363d; border-radius: 4px; }")
        beam_frame.setFixedHeight(36)
        bf_layout = QHBoxLayout(beam_frame)
        bf_layout.setContentsMargins(10, 0, 10, 0)
        bf_layout.setSpacing(6)
        beam_icon = QLabel("BEAM")
        beam_icon.setStyleSheet("color: #8b949e; font: bold 8pt 'Consolas'; background: transparent; border: none;")
        bf_layout.addWidget(beam_icon)
        self._lbl_beam_val = QLabel("-- nA")
        self._lbl_beam_val.setStyleSheet(
            f"color: {C.GREEN}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
        self._lbl_beam_val.setMinimumWidth(140)
        bf_layout.addWidget(self._lbl_beam_val)
        self._lbl_beam_status = QLabel("")
        self._lbl_beam_status.setStyleSheet(
            "color: transparent; font: bold 9pt 'Consolas'; background: transparent; border: none;")
        bf_layout.addWidget(self._lbl_beam_status)
        tl.addWidget(beam_frame)

        tl.addStretch()

        self._lbl_state = QLabel("IDLE")
        self._lbl_state.setStyleSheet(
            f"color: {C.DIM}; font: bold 11pt 'Consolas'; background: transparent;")
        tl.addWidget(self._lbl_state)

        root.addWidget(top)

        # ── main area ────────────────────────────────────────────────────
        body = QHBoxLayout()
        body.setContentsMargins(6, 4, 6, 6)
        body.setSpacing(6)

        # LEFT column: map (top) + log (bottom), stacked in a splitter
        left_splitter = QSplitter(Qt.Orientation.Vertical)

        # map container
        map_container = QWidget()
        mc_layout = QVBoxLayout(map_container)
        mc_layout.setContentsMargins(0, 0, 0, 0)
        mc_layout.setSpacing(2)

        self._canvas_label = QLabel()
        self._canvas_label.setStyleSheet(f"color: {C.ACCENT}; font: bold 9pt 'Consolas';")
        mc_layout.addWidget(self._canvas_label)

        self._map = HyCalScanMapWidget(self.all_modules)
        self._map.moduleClicked.connect(self._onCanvasClick)
        mc_layout.addWidget(self._map, stretch=1)

        # legend
        leg = QHBoxLayout()
        leg.setSpacing(4); leg.setContentsMargins(0, 0, 0, 0)
        for label, colour in [("Todo", C.MOD_TODO), ("Skipped", C.MOD_SKIPPED),
                               ("Moving", C.MOD_CURRENT), ("Dwell", C.MOD_DWELL),
                               ("Done", C.MOD_DONE), ("Error", C.MOD_ERROR),
                               ("Start", C.MOD_SELECTED), ("PbGlass", C.MOD_GLASS)]:
            sw = QLabel(); sw.setFixedSize(10, 10)
            sw.setStyleSheet(f"background: {colour}; border: none;")
            leg.addWidget(sw)
            ll = QLabel(label); ll.setStyleSheet(f"color: {C.DIM}; font: 8pt 'Consolas';")
            leg.addWidget(ll)
        leg.addStretch()
        mc_layout.addLayout(leg)

        left_splitter.addWidget(map_container)

        # event log
        log_group = QGroupBox("Event Log")
        log_layout = QVBoxLayout(log_group)
        log_layout.setContentsMargins(4, 4, 4, 4)
        self._log_text = QTextEdit()
        self._log_text.setReadOnly(True)
        log_layout.addWidget(self._log_text)
        left_splitter.addWidget(log_group)

        left_splitter.setStretchFactor(0, 3)  # map gets more space
        left_splitter.setStretchFactor(1, 1)  # log gets less
        body.addWidget(left_splitter, stretch=1)

        # RIGHT column: full-height controls
        right = QVBoxLayout()
        right.setSpacing(4)
        right.setContentsMargins(0, 0, 0, 0)

        self._buildScanControl(right)
        self._buildDirectControl(right)
        self._buildMotorStatus(right)
        self._buildPositionCheck(right)
        right.addStretch()

        right_widget = QWidget()
        right_widget.setLayout(right)
        right_widget.setFixedWidth(340)
        body.addWidget(right_widget)

        root.addLayout(body, stretch=1)
        self._updateCanvasLabel()

    # -- Scan Control --------------------------------------------------------

    def _buildScanControl(self, parent):
        sc = QGroupBox("Scan Control")
        lo = QVBoxLayout(sc)

        r = QHBoxLayout(); r.addWidget(QLabel("Path:"))
        self._profile_combo = QComboBox()
        self._profile_combo.addItems([self.NONE, self.AUTOGEN] + sorted(self._profiles.keys()))
        self._profile_combo.setCurrentText(self.AUTOGEN)
        self._profile_combo.activated.connect(self._onPathProfileChanged)
        r.addWidget(self._profile_combo, stretch=1); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("LG layers (0-2):"))
        self._lg_spin = QSpinBox(); self._lg_spin.setRange(0, MAX_LG_LAYERS)
        self._lg_spin.valueChanged.connect(self._onLgLayersChanged)
        r.addWidget(self._lg_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Start:"))
        names = [m.name for m in self.engine.path]
        self._start_combo = QComboBox(); self._start_combo.setEditable(True)
        self._start_combo.addItems(names); self._start_combo.setMinimumWidth(80)
        self._start_combo.activated.connect(self._onStartSelected)
        r.addWidget(self._start_combo)
        r.addWidget(QLabel("Count:"))
        self._count_spin = QSpinBox(); self._count_spin.setRange(0, len(names))
        self._count_spin.setSpecialValueText("all")
        self._count_spin.valueChanged.connect(lambda _: self._drawPathPreview())
        r.addWidget(self._count_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Dwell (s):"))
        self._dwell_spin = QDoubleSpinBox(); self._dwell_spin.setRange(1, 9999)
        self._dwell_spin.setValue(DEFAULT_DWELL); self._dwell_spin.setDecimals(0)
        r.addWidget(self._dwell_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Pos. threshold (mm):"))
        self._thresh_spin = QDoubleSpinBox(); self._thresh_spin.setRange(0.01, 10.0)
        self._thresh_spin.setValue(DEFAULT_POS_THRESHOLD); self._thresh_spin.setSingleStep(0.1)
        self._thresh_spin.setDecimals(2)
        r.addWidget(self._thresh_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Beam threshold (nA):"))
        self._beam_thresh_spin = QDoubleSpinBox(); self._beam_thresh_spin.setRange(0.0, 1000.0)
        self._beam_thresh_spin.setValue(DEFAULT_BEAM_THRESHOLD)
        self._beam_thresh_spin.setSingleStep(0.1); self._beam_thresh_spin.setDecimals(2)
        self._beam_thresh_spin.setSpecialValueText("off")
        r.addWidget(self._beam_thresh_spin); lo.addLayout(r)

        bf = QHBoxLayout()
        self._btn_start = QPushButton("Start Scan"); self._btn_start.setProperty("cssClass", "green")
        self._btn_start.clicked.connect(self._cmdStart); bf.addWidget(self._btn_start)
        self._btn_pause = QPushButton("Pause"); self._btn_pause.setProperty("cssClass", "warn")
        self._btn_pause.clicked.connect(self._cmdPause); bf.addWidget(self._btn_pause)
        self._btn_stop = QPushButton("Stop"); self._btn_stop.setProperty("cssClass", "danger")
        self._btn_stop.clicked.connect(self._cmdStop); bf.addWidget(self._btn_stop)
        lo.addLayout(bf)

        bf2 = QHBoxLayout()
        self._btn_skip = QPushButton("Skip Module"); self._btn_skip.clicked.connect(self._cmdSkip)
        bf2.addWidget(self._btn_skip)
        self._btn_ack = QPushButton("Ack Error"); self._btn_ack.setProperty("cssClass", "warn")
        self._btn_ack.clicked.connect(self._cmdAckError); bf2.addWidget(self._btn_ack)
        lo.addLayout(bf2)

        r = QHBoxLayout()
        self._lbl_progress = QLabel("Progress: --/--"); r.addWidget(self._lbl_progress)
        self._progress_bar = QProgressBar(); self._progress_bar.setMaximumWidth(140)
        r.addWidget(self._progress_bar); lo.addLayout(r)

        self._lbl_current = QLabel("Current:  --"); lo.addWidget(self._lbl_current)
        self._lbl_eta = QLabel("ETA:      --")
        self._lbl_eta.setStyleSheet(f"color: {C.DIM};"); lo.addWidget(self._lbl_eta)
        self._lbl_dwell_cd = QLabel("")
        self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};"); lo.addWidget(self._lbl_dwell_cd)

        parent.addWidget(sc)

    def _buildDirectControl(self, parent):
        dc = QGroupBox("Direct Control"); lo = QVBoxLayout(dc)
        self._btn_move = QPushButton("Move to Starting Point")
        self._btn_move.clicked.connect(self._cmdMoveToModule); lo.addWidget(self._btn_move)
        self._btn_reset = QPushButton("Reset to Beam Center")
        self._btn_reset.setProperty("cssClass", "accent")
        self._btn_reset.clicked.connect(self._cmdResetCenter); lo.addWidget(self._btn_reset)
        parent.addWidget(dc)

    def _buildMotorStatus(self, parent):
        ms = QGroupBox("Motor Status"); lo = QVBoxLayout(ms)
        for title, fields in [("X Motor", [
                ("Encoder", "x_encoder"), ("RBV", "x_rbv"), ("VAL", "x_val"),
                ("MOVN", "x_movn"), ("SPMG", "x_spmg"), ("VELO", "x_velo"),
                ("ACCL", "x_accl"), ("TDIR", "x_tdir"), ("MSTA", "x_msta"), ("ATHM", "x_athm")]),
            ("Y Motor", [
                ("Encoder", "y_encoder"), ("RBV", "y_rbv"), ("VAL", "y_val"),
                ("MOVN", "y_movn"), ("SPMG", "y_spmg"), ("VELO", "y_velo"),
                ("ACCL", "y_accl"), ("TDIR", "y_tdir"), ("MSTA", "y_msta"), ("ATHM", "y_athm")])]:
            tl = QLabel(title)
            tl.setStyleSheet(f"color: {C.ACCENT}; font: bold 9pt 'Consolas';")
            lo.addWidget(tl)
            g = QGridLayout(); g.setSpacing(2)
            half = (len(fields) + 1) // 2
            for i, (label, key) in enumerate(fields):
                c = 0 if i < half else 2; r = i % half
                ln = QLabel(f"{label}:")
                ln.setStyleSheet(f"color: {C.DIM}; font: 8pt 'Consolas';")
                ln.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                g.addWidget(ln, r, c)
                lv = QLabel("--"); lv.setMinimumWidth(90)
                g.addWidget(lv, r, c + 1)
                self._status_labels[key] = lv
            lo.addLayout(g)
        parent.addWidget(ms)

    def _buildPositionCheck(self, parent):
        pe = QGroupBox("Position Check"); lo = QVBoxLayout(pe)
        self._lbl_expected = QLabel("Expected: --"); lo.addWidget(self._lbl_expected)
        self._lbl_actual = QLabel("Actual:   --"); lo.addWidget(self._lbl_actual)
        self._lbl_error = QLabel("Diff:     --")
        self._lbl_error.setStyleSheet("font: bold 9pt 'Consolas';"); lo.addWidget(self._lbl_error)
        parent.addWidget(pe)

    def _disableControls(self):
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_skip, self._btn_ack, self._btn_move, self._btn_reset):
            w.setEnabled(False)
        for w in (self._start_combo, self._profile_combo, self._lg_spin,
                  self._count_spin, self._dwell_spin, self._thresh_spin,
                  self._beam_thresh_spin):
            w.setEnabled(False)

    # -- canvas helpers ------------------------------------------------------

    def _updateCanvasLabel(self):
        n_pwo4 = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        base = f"Module Map ({n_pwo4} PbWO4 + {n_lg} PbGlass)" if n_lg else f"Module Map ({n_pwo4} PbWO4)"
        if self._selected_mod_name:
            mod = self._mod_by_name.get(self._selected_mod_name)
            if mod:
                px, py = module_to_ptrans(mod.x, mod.y)
                base += f" | {mod.name} ({mod.mod_type})  HyCal({mod.x:.1f}, {mod.y:.1f})  ptrans({px:.1f}, {py:.1f})"
        self._canvas_label.setText(f" {base} ")

    def _updateCanvas(self):
        if self.observer:
            rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
            self._map.setMarkerPosition(*ptrans_to_module(rx, ry)); return

        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        colors = {}
        for m in self.all_modules:
            if m.name in self._scan_names or m.mod_type == "LMS": continue
            colors[m.name] = C.MOD_EXCLUDED if running else (
                C.MOD_GLASS if m.mod_type == "PbGlass" else C.MOD_PWO4_BG if m.mod_type == "PbWO4" else C.MOD_LMS)
        idle = eng.state in (ScanState.IDLE, ScanState.COMPLETED)
        count = self._count_spin.value()
        si = self._selected_start_idx
        ei = min(si + count, len(eng.path)) if count > 0 else len(eng.path)
        if not idle:
            si = eng.current_idx; ei = getattr(eng, '_end_idx', len(eng.path))
        for i, mod in enumerate(eng.path):
            if i == eng.current_idx and eng.state == ScanState.DWELLING:
                colors[mod.name] = C.MOD_DWELL
            elif i == eng.current_idx and eng.state in (ScanState.MOVING, ScanState.PAUSED):
                colors[mod.name] = C.MOD_CURRENT
            elif i in eng.error_modules: colors[mod.name] = C.MOD_ERROR
            elif i in eng.completed: colors[mod.name] = C.MOD_DONE
            elif idle and i == self._selected_start_idx: colors[mod.name] = C.MOD_SELECTED
            elif i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
            else: colors[mod.name] = C.MOD_TODO
        self._map.setModuleColors(colors)
        if idle:
            self._drawPathPreview(); self._map.setDashPreview([])
        elif running:
            self._map.setPathPreview([])
            self._map.setDashPreview([self._map.modCenter(eng.path[i]) for i in range(eng.current_idx + 1, ei)])
        rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
        self._map.setMarkerPosition(*ptrans_to_module(rx, ry))

    def _drawPathPreview(self):
        path = self.engine.path
        s = self._selected_start_idx
        if s >= len(path): self._map.setPathPreview([]); return
        c = self._count_spin.value()
        e = min(s + c, len(path)) if c > 0 else len(path)
        self._map.setPathPreview([self._map.modCenter(path[i]) for i in range(s, e)])

    def _onCanvasClick(self, name):
        if self._selected_mod_name == name:
            self._selected_mod_name = None; self._map.setHighlight(None)
            self._updateCanvasLabel(); return
        self._selected_mod_name = name
        if name in self._scan_name_to_idx:
            self._selected_start_idx = self._scan_name_to_idx[name]
            idx = self._start_combo.findText(name)
            if idx >= 0: self._start_combo.setCurrentIndex(idx)
            self._drawPathPreview()
        self._map.setHighlight(name); self._updateCanvasLabel()

        # When idle, pop up module info dialog
        idle = self.engine.state in (ScanState.IDLE, ScanState.COMPLETED)
        if idle and not self.observer:
            mod = self._mod_by_name.get(name)
            if mod:
                dlg = ModuleInfoDialog(mod, self.ep, self._log, parent=self)
                dlg.exec()

    # -- commands ------------------------------------------------------------

    def _onStartSelected(self, _):
        name = self._start_combo.currentText()
        for i, m in enumerate(self.engine.path):
            if m.name == name:
                self._selected_start_idx = i; self._drawPathPreview(); break

    def _onPathProfileChanged(self, _):
        name = self._profile_combo.currentText()
        if name == self._active_profile: return
        self._active_profile = name
        if name == self.AUTOGEN:
            self._lg_spin.setEnabled(True); self._onLgLayersChanged(force=True); return
        self._lg_spin.setEnabled(False)
        if name == self.NONE:
            self.scan_modules = []; self._scan_names = set()
            self.engine = ScanEngine(self.ep, [], self._log)
            self._scan_name_to_idx = {}; self._selected_start_idx = 0
            self._start_combo.clear(); self._count_spin.setMaximum(0); self._count_spin.setValue(0)
            self._map.setPathPreview([]); self._map.setDashPreview([])
            self._map.setHighlight(None); self._selected_mod_name = None
            self._updateCanvasLabel()
            self._log("Path: none (direct control only)")
            return
        mod_by_name = {m.name: m for m in self.all_modules}
        path_mods = [mod_by_name[n] for n in self._profiles.get(name, []) if n in mod_by_name]
        if not path_mods:
            self._log(f"Profile '{name}' empty", level="error"); return
        self.scan_modules = path_mods; self._scan_names = {m.name for m in path_mods}
        self.engine = ScanEngine(self.ep, path_mods, self._log)
        self.engine.path = path_mods
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(self.engine.path)}
        self._selected_start_idx = 0
        ns = [m.name for m in self.engine.path]
        self._start_combo.clear(); self._start_combo.addItems(ns)
        self._count_spin.setMaximum(len(ns)); self._count_spin.setValue(0)
        self._updateCanvasLabel()
        self._log(f"Path profile: {name} ({len(path_mods)} modules)")

    def _onLgLayersChanged(self, value=0, force=False):
        if self._active_profile != self.AUTOGEN: return
        nl = self._lg_spin.value()
        if nl == self._lg_layers and not force: return
        self._lg_layers = nl
        self.scan_modules = filter_scan_modules(self.all_modules, nl, self._lg_sx, self._lg_sy)
        self._scan_names = {m.name for m in self.scan_modules}
        self.engine = ScanEngine(self.ep, self.scan_modules, self._log)
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(self.engine.path)}
        self._selected_start_idx = 0
        ns = [m.name for m in self.engine.path]
        self._start_combo.clear(); self._start_combo.addItems(ns)
        self._count_spin.setMaximum(len(ns)); self._count_spin.setValue(0)
        self._updateCanvasLabel()
        np_ = sum(1 for m in self.scan_modules if m.mod_type == "PbWO4")
        ng = sum(1 for m in self.scan_modules if m.mod_type == "PbGlass")
        self._log(f"LG layers: {nl} ({np_} PbWO4 + {ng} PbGlass = {len(self.scan_modules)})")

    def _cmdStart(self):
        self._onStartSelected(0)
        path = self.engine.path; s = self._selected_start_idx; c = self._count_spin.value()
        e = min(s + c, len(path)) if c > 0 else len(path)
        oob = [path[i].name for i in range(s, e) if not ptrans_in_limits(*module_to_ptrans(path[i].x, path[i].y))]
        if oob:
            ns = ", ".join(oob[:5]) + (f" ... ({len(oob)} total)" if len(oob) > 5 else "")
            self._log(f"BLOCKED: {len(oob)} modules outside limits: {ns}", level="error")
            QMessageBox.critical(self, "Out of Bounds", f"{len(oob)} outside limits:\n{ns}"); return
        self.engine.dwell_time = self._dwell_spin.value()
        self.engine.pos_threshold = self._thresh_spin.value()
        self.engine.beam_threshold = self._beam_thresh_spin.value()
        self.engine.start(self._selected_start_idx, count=c)

    def _cmdPause(self):
        eng = self.engine
        if eng.state == ScanState.PAUSED:
            eng.resume_scan(); self._btn_pause.setText("Pause")
        elif eng.state in (ScanState.MOVING, ScanState.DWELLING):
            eng.pause_scan(); self._btn_pause.setText("Resume")

    def _cmdStop(self):
        if self.engine.state != ScanState.IDLE:
            self.engine.stop_scan(); self._btn_pause.setText("Pause")
        else:
            epics_stop(self.ep); self._log("Motors stopped")

    def _cmdSkip(self):      self.engine.skip_module()
    def _cmdAckError(self):  self.engine.acknowledge_error()

    def _cmdMoveToModule(self):
        self._onStartSelected(0)
        if not self.engine.path: return
        mod = self.engine.path[self._selected_start_idx]
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if not epics_move_to(self.ep, px, py):
            self._log(f"BLOCKED: outside limits", level="error")

    def _cmdResetCenter(self):
        self._log(f"Resetting to beam centre ptrans({BEAM_CENTER_X}, {BEAM_CENTER_Y})")
        epics_move_to(self.ep, BEAM_CENTER_X, BEAM_CENTER_Y)

    # -- logging -------------------------------------------------------------

    def _log(self, msg, level="info"):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {level.upper().ljust(5)} {msg}"
        self._log_lines.append(line)
        if self._log_file:
            self._log_file.write(line + "\n"); self._log_file.flush()
        self._logSignal.emit(line, level)

    def _appendLog(self, line, level):
        colors = {"info": C.TEXT, "warn": C.YELLOW, "error": C.RED}
        c = colors.get(level, C.DIM)
        self._log_text.append(
            f'<span style="color:{c}; font-family:Consolas; font-size:9pt;">'
            f'{html_mod.escape(line)}</span>')

    # -- polling (5 Hz) ------------------------------------------------------

    def _poll(self):
        self._updateStatus()
        self._updateCanvas()
        self._updateScanInfo()
        self._updateButtons()
        self._updateBeamDisplay()

    def _updateBeamDisplay(self):
        bc = self.ep.get("beam_cur", None)
        if bc is None:
            self._lbl_beam_val.setText("-- nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.DIM}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("")
            return
        thresh = self._beam_thresh_spin.value()
        tripped = self.engine.beam_tripped
        if tripped:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.RED}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("TRIP")
            self._lbl_beam_status.setStyleSheet(f"color: {C.RED}; font: bold 10pt 'Consolas'; background: transparent; border: none;")
        elif thresh > 0 and bc < thresh:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.YELLOW}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("LOW")
            self._lbl_beam_status.setStyleSheet(f"color: {C.YELLOW}; font: bold 10pt 'Consolas'; background: transparent; border: none;")
        else:
            self._lbl_beam_val.setText(f"{bc:.2f} nA")
            self._lbl_beam_val.setStyleSheet(f"color: {C.GREEN}; font: bold 14pt 'Consolas'; background: transparent; border: none;")
            self._lbl_beam_status.setText("")

    def _updateStatus(self):
        for key, lbl in self._status_labels.items():
            val = self.ep.get(key, "--")
            if val == "--" or val is None:
                lbl.setText("--"); lbl.setStyleSheet(f"color: {C.DIM};"); continue
            if key.endswith("_msta"):        txt = f"0x{int(val):X}"
            elif key.endswith("_spmg"):      txt = f"{SPMG_LABELS.get(int(val), '?')}({int(val)})"
            elif key.endswith(("_movn", "_athm", "_tdir")): txt = str(int(val))
            elif isinstance(val, float):     txt = f"{val:.3f}"
            else:                            txt = str(val)
            fg = C.TEXT
            if key.endswith("_movn") and int(val) == 1: fg = C.YELLOW
            elif key.endswith("_spmg"):
                sv = int(val)
                fg = C.RED if sv == SPMG.STOP else C.ORANGE if sv == SPMG.PAUSE else C.GREEN
            lbl.setText(txt); lbl.setStyleSheet(f"color: {fg};")
        # position check
        eng = self.engine
        rx, ry = self.ep.get("x_rbv", 0.0), self.ep.get("y_rbv", 0.0)
        self._lbl_actual.setText(f"Actual:   ({rx:.3f}, {ry:.3f})")
        mod = eng.current_module if eng.state != ScanState.IDLE else None
        if mod is None and eng.path and 0 <= self._selected_start_idx < len(eng.path):
            mod = eng.path[self._selected_start_idx]
        scanning = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        if mod:
            px, py = module_to_ptrans(mod.x, mod.y)
            err = math.sqrt((rx - px)**2 + (ry - py)**2)
            self._lbl_expected.setText(f"Expected: ({px:.3f}, {py:.3f}) {mod.name}")
            if scanning:
                ef = C.RED if err > eng.pos_threshold else C.GREEN
                self._lbl_error.setText(f"Diff:     {err:.3f} mm")
                self._lbl_error.setStyleSheet(f"color: {ef}; font: bold 9pt 'Consolas';")
            else:
                self._lbl_error.setText(f"Diff:     {err:.3f} mm (idle)")
                self._lbl_error.setStyleSheet(f"color: {C.DIM}; font: bold 9pt 'Consolas';")
        else:
            self._lbl_expected.setText("Expected: --")
            self._lbl_error.setText("Diff:     --")
            self._lbl_error.setStyleSheet(f"color: {C.DIM}; font: bold 9pt 'Consolas';")

    def _updateScanInfo(self):
        eng = self.engine
        sc = {ScanState.IDLE: C.DIM, ScanState.MOVING: C.YELLOW, ScanState.DWELLING: C.GREEN,
              ScanState.PAUSED: C.ORANGE, ScanState.ERROR: C.RED, ScanState.COMPLETED: C.ACCENT}
        self._lbl_state.setText(eng.state)
        self._lbl_state.setStyleSheet(f"color: {sc.get(eng.state, C.DIM)}; font: bold 11pt 'Consolas'; background: transparent;")
        done = len(eng.completed)
        s = getattr(eng, '_start_idx', 0); e = getattr(eng, '_end_idx', len(eng.path))
        total = e - s
        self._lbl_progress.setText(f"Progress: {done}/{total}")
        self._progress_bar.setMaximum(max(total, 1)); self._progress_bar.setValue(done)
        mod = eng.current_module
        self._lbl_current.setText(f"Current:  {mod.name}" if mod else "Current:  --")
        if eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED):
            eta = eng.eta_seconds; h, rem = divmod(int(eta), 3600); m, s = divmod(rem, 60)
            self._lbl_eta.setText(f"ETA:      {h}h {m:02d}m {s:02d}s")
        elif eng.state == ScanState.IDLE and eng.path:
            vx = self.ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
            vy = self.ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
            eta = estimate_scan_time(eng.path, self._selected_start_idx,
                                     self._count_spin.value(), self._dwell_spin.value(), vx, vy)
            if eta > 0:
                h, rem = divmod(int(eta), 3600); m, s = divmod(rem, 60)
                self._lbl_eta.setText(f"ETA:      ~{h}h {m:02d}m {s:02d}s")
            else: self._lbl_eta.setText("ETA:      --")
        else: self._lbl_eta.setText("ETA:      --")
        if eng.beam_tripped:
            self._lbl_dwell_cd.setText("Dwell:    BEAM TRIP -- waiting")
            self._lbl_dwell_cd.setStyleSheet(f"color: {C.RED}; font: bold 9pt 'Consolas';")
        elif eng.state == ScanState.DWELLING:
            self._lbl_dwell_cd.setText(f"Dwell:    {eng.dwell_remaining:.1f}s remaining")
            self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};")
        else:
            self._lbl_dwell_cd.setText(""); self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};")

    def _updateButtons(self):
        if self.observer: return
        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        has_path = len(eng.path) > 0
        self._btn_start.setEnabled(not running and has_path)
        self._btn_pause.setEnabled(running)
        self._btn_stop.setEnabled(True)
        self._btn_skip.setEnabled(eng.state == ScanState.DWELLING)
        self._btn_ack.setEnabled(eng.state == ScanState.ERROR)
        self._start_combo.setEnabled(not running and has_path)
        self._count_spin.setEnabled(not running and has_path)
        self._profile_combo.setEnabled(not running)
        self._lg_spin.setEnabled(not running and self._active_profile == self.AUTOGEN)
        self._btn_move.setEnabled(not running and has_path)
        self._btn_reset.setEnabled(not running)

    def closeEvent(self, e):
        if self._log_file: self._log_file.close()
        super().closeEvent(e)


# ============================================================================
#  MAIN
# ============================================================================

PATHS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paths.json")

def main():
    parser = argparse.ArgumentParser(description="HyCal Snake Scan")
    parser.add_argument("--expert", action="store_true")
    parser.add_argument("--observer", action="store_true")
    parser.add_argument("--database", default=DEFAULT_DB_PATH)
    parser.add_argument("--paths", default=PATHS_FILE)
    args = parser.parse_args()

    all_modules = load_modules(args.database)
    by_type = {}
    for m in all_modules:
        by_type[m.mod_type] = by_type.get(m.mod_type, 0) + 1
    print(f"Loaded {len(all_modules)} modules from {args.database}")
    for t, n in sorted(by_type.items()):
        print(f"  {t}: {n}")

    profiles = {}
    if os.path.exists(args.paths):
        with open(args.paths) as f:
            profiles = json.load(f)
        print(f"Loaded {len(profiles)} path profiles")

    observer = args.observer
    simulation = not args.expert and not observer

    if observer:        ep = ObserverEPICS()
    elif simulation:    ep = SimulatedEPICS()
    else:               ep = RealEPICS(writable=True)

    n_ok, n_total = ep.connect()
    if not simulation:
        print(f"EPICS: {n_ok}/{n_total} PVs connected")
        for pv in ep.disconnected_pvs():
            print(f"  NOT connected: {pv}")

    app = QApplication(sys.argv)
    win = SnakeScanWindow(ep, simulation, all_modules, profiles, observer=observer)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
