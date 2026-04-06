"""
Scan engine for HyCal snake-scan calibration.

Builds a snake-pattern scan path over HyCal modules and drives the
transporter through it, handling dwell timing, beam-trip recovery,
pause/resume, and position-error acknowledgement.
"""

from __future__ import annotations
import math
import threading
import time
from typing import List, Optional, Set

from scan_utils import Module, module_to_ptrans, ptrans_in_limits, filter_scan_modules
from scan_epics import (
    SPMG, epics_move_to, epics_stop, epics_pause, epics_resume,
    epics_is_moving, epics_read_rbv,
)


# ============================================================================
#  CONSTANTS
# ============================================================================

DEFAULT_DWELL = 120.0         # seconds
DEFAULT_POS_THRESHOLD = 0.5   # mm
DEFAULT_BEAM_THRESHOLD = 0.3  # nA
MOVE_TIMEOUT = 900.0          # seconds per single move (15 min)
DEFAULT_VELO_X = 50.0         # mm/s
DEFAULT_VELO_Y = 5.0          # mm/s
MAX_LG_LAYERS = 2


# ============================================================================
#  PATH BUILDING
# ============================================================================

def build_scan_path(scan_modules):
    from collections import defaultdict
    sectors = defaultdict(list)
    for m in scan_modules:
        sectors[m.sector].append(m)
    center_mods = sectors.pop("Center", [])
    lg_sectors = {k: v for k, v in sectors.items() if k != "LMS"}

    def _snake_sector(mods, col_right, row_down):
        rows = defaultdict(list)
        for m in mods:
            rows[round(m.y, 1)].append(m)
        sorted_rows = sorted(rows.keys(), reverse=row_down)
        path = []
        going = col_right
        for ry in sorted_rows:
            cols = sorted(rows[ry], key=lambda m: m.x, reverse=not going)
            path.extend(cols)
            going = not going
        return path, going

    path = []
    going = True
    if center_mods:
        seg, going = _snake_sector(center_mods, going, True)
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
#  SCAN STATE
# ============================================================================

class ScanState:
    IDLE = "IDLE"; MOVING = "MOVING"; DWELLING = "DWELLING"
    PAUSED = "PAUSED"; ERROR = "ERROR"; COMPLETED = "COMPLETED"


# ============================================================================
#  SCAN ENGINE
# ============================================================================

class ScanEngine:
    def __init__(self, epics, modules, log_fn):
        self.ep = epics
        self.all_modules = modules
        # ``modules`` is already the ordered path (caller-prepared)
        self.path = list(modules)
        self.log = log_fn

        self.dwell_time = DEFAULT_DWELL
        self.pos_threshold = DEFAULT_POS_THRESHOLD
        self.beam_threshold = DEFAULT_BEAM_THRESHOLD

        self.state = ScanState.IDLE
        self.current_idx = 0
        self.dwell_remaining = 0.0
        self.completed: Set[int] = set()
        self.error_modules: Set[int] = set()
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
