#!/usr/bin/env python3
"""
Scintillator–HyCal Coincidence Monitor
=======================================
Connects to a running prad2_server (HTTP REST API, port 5051 by default),
iterates through every event in the loaded file, and accumulates per-module
coincidence statistics between the four upstream scintillators (V1–V4) and
every HyCal module.

Two event-selection modes are provided:

  AND mode — rate(V_i, M) = N(V_i fired AND M is best) / N(M is best)
             Only events where at least one scintillator fired are counted.
             Gives P(V_i fired | M was the highest-ADC module).

  OR  mode — rate(V_i, M) = N(V_i fired AND M is best) / N(V_i fired)
             All HyCal cluster events counted regardless of scintillator.
             Gives P(M is best | V_i fired).

A channel "fired" when it has at least one FADC peak whose height (above
pedestal) exceeds the user-specified threshold.

The bottom half of the window shows individual waveforms: the currently
selected scintillator (V1–V4) and the HyCal module last clicked on the map,
both fetched from the server for the event number entered in the Event Browser.

Usage
-----
    python scripts/coincidence_monitor.py [--url http://HOST:PORT]
                                          [--theme dark|light]
"""
from __future__ import annotations

import argparse
import json as json_mod
import math
import sys
import threading
import time
import multiprocessing
import os
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

import urllib.request
import urllib.error

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QDoubleSpinBox, QLineEdit, QProgressBar,
    QSplitter, QSizePolicy, QButtonGroup, QRadioButton, QGroupBox,
    QSpinBox, QFrame, QMessageBox, QCheckBox, QFileDialog, QScrollArea,
    QListWidget, QAbstractItemView,
)
from PyQt6.QtCore import Qt, QThread, QTimer, pyqtSignal, QRectF, QPointF
from PyQt6.QtGui import (
    QColor, QFont, QPen, QPainter, QPolygonF,
)

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget, PALETTES, PALETTE_NAMES,
    cmap_qcolor,
    apply_theme_palette, set_theme, available_themes, THEME, themed,
)

try:
    import prad2py as _prad2py          # type: ignore
    _HAVE_PRAD2PY = True
except ImportError:
    _prad2py = None                     # type: ignore
    _HAVE_PRAD2PY = False


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON  = DB_DIR / "hycal_modules.json"
DAQ_MAP_JSON  = DB_DIR / "hycal_daq_map.json"
DAQ_CFG_JSON  = DB_DIR / "daq_config.json"

DEFAULT_URL = "http://localhost:5051"

# Scintillator names (keys resolved at runtime from daq_map + daq_config)
SCINTILLATORS: Dict[str, str] = {
    "V1": "V1",
    "V2": "V2",
    "V3": "V3",
    "V4": "V4",
}

# Default thresholds (FADC peak height above pedestal, ADC counts)
DEFAULT_SCINT_THR = 500.0
DEFAULT_HYCAL_THR = 500.0

# Default HyCal signal time window (ns) — same FADC window as scintillators
DEFAULT_HYCAL_TMIN = 160.0
DEFAULT_HYCAL_TMAX = 200.0

# Default max number of local ADC maxima allowed in HyCal (1 = single cluster only)
DEFAULT_MAX_LOCAL_MAXIMA = 1

# Event-selection modes
MODE_AND = "AND"   # rate = N(Vi AND M_best) / N(M_best)  — P(Vi fired | M is best)
MODE_OR  = "OR"    # rate = N(Vi AND M_best) / N(Vi fired) — P(M is best | Vi fired)

# Skip calibration/background events by trigger bit OR by module occupancy.
LMS_TRIGGER_BIT   = 1 << 24   # bit 24 — LMS light source
ALPHA_TRIGGER_BIT = 1 << 25   # bit 25 — alpha source
SKIP_TRIGGER_MASK = LMS_TRIGGER_BIT | ALPHA_TRIGGER_BIT
LMS_MAX_MODULES   = 1000      # more than this many modules above threshold → LMS

# Cluster-size cut: events with fewer than this many HyCal modules above threshold
# are discarded as isolated electronic noise / discharge artefacts.
DEFAULT_MIN_CLUSTER_MODS = 2

# Top-level display modes
DISPLAY_COINC   = "coinc"    # coincidence rate / occupancy statistics
DISPLAY_INSTANT = "instant"  # live event-by-event ADC display

# Map display modes (used within DISPLAY_COINC)
VIEW_COINC   = "coincidence"  # colour = coincidence rate with selected scintillator
VIEW_OCC     = "occupancy"    # colour = number of events module fired above threshold
VIEW_INSTANT = "instant"      # colour = per-module ADC signal for current event

# Veto scintillator motor PVs  (prad:vetoN.VAL = setpoint, .RBV = read-back)
VETO_PV_BASES: Dict[str, str] = {
    "V1": "prad:veto1",
    "V2": "prad:veto2",
    "V3": "prad:veto3",
    "V4": "prad:veto4",
}
VETO_POLL_MS = 500   # 2 Hz

N_WORKERS  = 16    # parallel HTTP workers for event scanning
BATCH_SIZE = 64    # events per processing batch (small keeps in-flight JSON bounded)

CLK_MHZ = 250.0    # FADC clock (for x-axis in ns)


# ===========================================================================
#  Veto motor position reader
# ===========================================================================

class VetoMotorController:
    """Reads and writes position PVs for the four veto scintillator motors.

    PVs per motor:
        prad:vetoN.VAL  — setpoint (write to command a move)
        prad:vetoN.RBV  — actual read-back position
        prad:vetoN.MOVN — 1 while moving, 0 when at rest

    Degrades gracefully when pyepics is not installed.
    """

    def __init__(self) -> None:
        self._pvs: Dict[str, Any] = {}
        self._epics_ok = False
        try:
            import epics as _epics
            for vname, base in VETO_PV_BASES.items():
                self._pvs[f"{vname}_val"]  = _epics.PV(f"{base}.VAL")
                self._pvs[f"{vname}_rbv"]  = _epics.PV(f"{base}.RBV")
                self._pvs[f"{vname}_movn"] = _epics.PV(f"{base}.MOVN")
            self._epics_ok = True
        except ImportError:
            pass

    def get(self, key: str) -> Optional[float]:
        pv = self._pvs.get(key)
        if pv is None or not pv.connected:
            return None
        v = pv.get()
        return float(v) if v is not None else None

    def put(self, key: str, value: float) -> bool:
        """Write value to a VAL PV. Returns True on success."""
        pv = self._pvs.get(key)
        if pv is None or not pv.connected:
            return False
        pv.put(value)
        return True

    @property
    def available(self) -> bool:
        return self._epics_ok


# ===========================================================================
#  HTTP helpers
# ===========================================================================

def _load_crate_to_roc(path: Path) -> Dict[int, int]:
    """Return {crate_index: roc_tag} from daq_config.json roc_tags list.
    Only 'roc' type entries are included (not ti_slave, tdc, gem, etc.).
    """
    mapping: Dict[int, int] = {}
    try:
        with open(path) as f:
            cfg = json_mod.load(f)
        for entry in cfg.get("roc_tags", []):
            if entry.get("type") != "roc":
                continue
            crate = entry.get("crate")
            tag_raw = entry.get("tag", "")
            try:
                tag = int(tag_raw, 16) if isinstance(tag_raw, str) else int(tag_raw)
                mapping[int(crate)] = tag
            except (ValueError, TypeError):
                pass
    except Exception:
        pass
    return mapping


def _load_daq_map(path: Path,
                  crate_to_roc: Optional[Dict[int, int]] = None) -> Dict[str, str]:
    """Return {module_name: "roc_tag_slot_channel"} from hycal_daq_map.json.

    The event JSON produced by the C++ server uses the actual ROC tag (not the
    sequential crate index) as the first component of the channel key.  Pass
    crate_to_roc so the keys produced here match what the server emits.
    """
    with open(path) as f:
        entries = json_mod.load(f)
    result: Dict[str, str] = {}
    for e in entries:
        crate = e["crate"]
        roc = crate_to_roc.get(crate, crate) if crate_to_roc else crate
        result[e["name"]] = f"{roc}_{e['slot']}_{e['channel']}"
    return result


def _build_w_module_layers(modules_path: Path) -> Dict[str, int]:
    """Return {W_module_name: layer} where layer 1 is immediately around the beam hole.

    The central hole occupies rows 17-18, cols 17-18.  A module's layer is its
    Chebyshev distance to the nearest hole cell:
        layer = max(max(0, 17-row, row-18), max(0, 17-col, col-18))
    """
    with open(modules_path) as f:
        mods = json_mod.load(f)
    result: Dict[str, int] = {}
    for m in mods:
        name = m['n']
        if not name.startswith('W'):
            continue
        dr = max(0, 17 - m['row'], m['row'] - 18)
        dc = max(0, 17 - m['col'], m['col'] - 18)
        result[name] = max(dr, dc)
    return result


def _http_get(url: str, timeout: float = 5.0) -> Optional[dict]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json_mod.loads(resp.read())
    except Exception:
        return None


def _build_neighbor_map(modules) -> Dict[str, frozenset]:
    """Return name→frozenset-of-neighbor-names for HyCal modules.

    Two modules are neighbors when they share an edge (touching sides, not just
    a corner).  Works for mixed W/PbGl geometries because it checks actual
    center-distance against the sum of half-widths with a 1 mm tolerance.
    """
    phys = [(m.name, m.x, m.y, m.sx, m.sy)
            for m in modules if m.mod_type != "LMS"]
    neighbors: Dict[str, List[str]] = {name: [] for name, *_ in phys}
    tol = 1.0   # mm
    for i, (n1, x1, y1, sx1, sy1) in enumerate(phys):
        for n2, x2, y2, sx2, sy2 in phys[i + 1:]:
            dx = abs(x1 - x2)
            dy = abs(y1 - y2)
            share_v = (abs(dx - (sx1 + sx2) / 2) < tol
                       and dy < min(sy1, sy2) / 2 + tol)
            share_h = (abs(dy - (sy1 + sy2) / 2) < tol
                       and dx < min(sx1, sx2) / 2 + tol)
            if share_v or share_h:
                neighbors[n1].append(n2)
                neighbors[n2].append(n1)
    return {name: frozenset(nbs) for name, nbs in neighbors.items()}


def _count_local_maxima(module_adc: Dict[str, float],
                        neighbor_map: Dict[str, frozenset]) -> int:
    """Count HyCal modules that are strictly higher than all their neighbors.

    Neighbors absent from module_adc are treated as ADC = 0 (below threshold).
    """
    count = 0
    for name, adc in module_adc.items():
        if all(adc > module_adc.get(nb, 0.0) for nb in neighbor_map.get(name, ())):
            count += 1
    return count


def _build_tuple_role_map(scint_keys: Dict[str, str],
                          module_keys: Dict[str, str]) -> Dict[tuple, tuple]:
    """Build (roc_tag, slot, chan) → (kind, name) lookup.

    kind=0 → scintillator, kind=1 → HyCal module.
    Avoids f-string formatting in the per-channel hot loop.
    """
    role: Dict[tuple, tuple] = {}
    for sname, key in scint_keys.items():
        roc, s, c = key.split("_")
        role[(int(roc), int(s), int(c))] = (0, sname)
    for mname, key in module_keys.items():
        roc, s, c = key.split("_")
        role[(int(roc), int(s), int(c))] = (1, mname)
    return role


# Worker-process global for the update queue.  ``forkserver``/``spawn``
# refuse to pickle a Queue in apply_async args ("Queue objects should
# only be shared between processes through inheritance"), so we install
# it once per worker via Pool(initializer=...).
_WORKER_UPDATE_Q = None


def _init_worker(q):
    global _WORKER_UPDATE_Q
    _WORKER_UPDATE_Q = q


def _process_files(args: dict) -> None:
    """Process a list of EVIO files in a subprocess.

    Periodically pushes incremental count deltas to ``args["update_queue"]``
    so the parent QThread can update the UI while subprocesses are still
    running.  Returns None — all results flow through the queue.
    """
    paths           = args["paths"]
    tuple_to_role   = args["tuple_to_role"]
    scint_thr       = args["scint_thr"]
    hycal_thr       = args["hycal_thr"]
    scint_t_min     = args["scint_t_min"]
    scint_t_max     = args["scint_t_max"]
    hycal_t_min     = args["hycal_t_min"]
    hycal_t_max     = args["hycal_t_max"]
    neighbor_map    = args["neighbor_map"]
    max_lm          = args["max_lm"]
    mode            = args["mode"]
    min_mods        = args["min_mods"]
    scint_names     = args["scint_names"]
    mod_names       = args["mod_names"]
    skip_mask       = args["skip_mask"]
    lms_max         = args["lms_max"]
    worker_id       = args["worker_id"]
    n_files_in_chunk = len(paths)
    update_q        = _WORKER_UPDATE_Q
    push_interval   = args.get("push_interval", 0.5)

    # Delta accumulators — reset to zero after each push to the parent.
    d_module_hits    = {m: 0 for m in mod_names}
    d_scint_hits     = {s: 0 for s in scint_names}
    d_scint_hits_any = {s: 0 for s in scint_names}
    d_coincidences   = {s: {m: 0 for m in mod_names} for s in scint_names}
    d_processed      = 0

    # Per-file progress state — included in every push so the parent UI can
    # render one progress bar per worker.
    cur_file_idx     = -1
    cur_basename     = ""
    cur_records_done = 0
    cur_records_tot  = 0
    finished_chunk   = False

    # Split the channel list into a tiny scint table (visited first so we can
    # short-circuit non-coincidence events) and a HyCal table.  The ROC tag
    # is bucket-indexed up front to save dict.get() overhead per channel.
    scint_locs: List[tuple] = []   # [(roc_tag, slot, chan, name), ...]
    hycal_locs_by_roc: Dict[int, List[tuple]] = {}
    for (rt, sn, cn), (kind, name) in tuple_to_role.items():
        if kind == 0:
            scint_locs.append((rt, sn, cn, name))
        else:
            hycal_locs_by_roc.setdefault(rt, []).append((sn, cn, name))
    # Sort each ROC's HyCal channel list by (slot, chan) so consecutive
    # iterations stay on the same slot — lets us cache slot.channel() lookups.
    for _lst in hycal_locs_by_roc.values():
        _lst.sort()

    # In AND mode, or when the scint time cut is active, we never count an
    # event whose scintillators don't fire — so skipping HyCal in that case
    # is safe and saves all per-channel work for those events.
    require_scint = (mode == MODE_AND or scint_t_min > -math.inf)

    # Pre-compute sample-index bounds for the time-window peak search.
    # FADC250 ticks at 250 MHz → 4 ns per sample.  N_PED is the number of
    # leading samples used to estimate pedestal — chosen large enough to be
    # statistically stable but small enough never to overlap a physics
    # signal (those arrive after ~100 ns at the earliest).
    N_PED = 30
    _ns_per_sample = 4.0
    _NO_UPPER = 1 << 30

    if scint_t_min > -math.inf:
        s_lo = max(N_PED, int(scint_t_min / _ns_per_sample))
    else:
        s_lo = N_PED
    if scint_t_max < math.inf:
        s_hi = int(scint_t_max / _ns_per_sample) + 1
    else:
        s_hi = _NO_UPPER

    if hycal_t_min > -math.inf:
        h_lo = max(N_PED, int(hycal_t_min / _ns_per_sample))
    else:
        h_lo = N_PED
    if hycal_t_max < math.inf:
        h_hi = int(hycal_t_max / _ns_per_sample) + 1
    else:
        h_hi = _NO_UPPER

    def _flush():
        nonlocal d_processed
        if update_q is None:
            return
        update_q.put({
            "module_hits":     dict(d_module_hits),
            "scint_hits":      dict(d_scint_hits),
            "scint_hits_any":  dict(d_scint_hits_any),
            "coincidences":    {s: dict(c) for s, c in d_coincidences.items()},
            "processed":       d_processed,
            "worker_id":       worker_id,
            "n_files":         n_files_in_chunk,
            "file_idx":        cur_file_idx,
            "file_basename":   cur_basename,
            "records_done":    cur_records_done,
            "records_total":   cur_records_tot,
            "finished":        finished_chunk,
        })
        for m in mod_names:
            d_module_hits[m] = 0
        for s in scint_names:
            d_scint_hits[s] = 0
            d_scint_hits_any[s] = 0
            for m in mod_names:
                d_coincidences[s][m] = 0
        d_processed = 0

    try:
        import prad2py as _p2       # type: ignore
    except ImportError as e:
        sys.stderr.write(f"[worker {worker_id} pid={os.getpid()}] "
                         f"prad2py import failed: {e}\n")
        sys.stderr.flush()
        finished_chunk = True
        _flush()
        return None

    dec      = _p2.dec
    cfg      = dec.load_daq_config()
    ch       = dec.EvChannel()
    ch.set_config(cfg)
    # WaveAnalyzer not needed in the parallel worker — the per-event hot
    # loop uses a numpy-only peak-in-window heuristic instead.

    # If prad2py was built with the slot-batched fast path, use it.  The
    # method returns a single float32 array per slot (one C call instead
    # of 16 separate `.samples` numpy allocations), which is the only
    # remaining ~10× speedup at this layer.
    has_batch = hasattr(dec.SlotData, "peak_in_window")

    last_push = time.monotonic()

    for fidx, path in enumerate(paths):
        cur_file_idx     = fidx
        cur_basename     = os.path.basename(path)
        cur_records_done = 0
        cur_records_tot  = 0

        st = ch.open_auto(path)
        if st != dec.Status.success:
            sys.stderr.write(f"[worker {worker_id} pid={os.getpid()}] "
                             f"open_auto({path}) -> {st}\n")
            sys.stderr.flush()
            ch.close()
            _flush()   # let parent know this file was skipped
            continue

        if ch.is_random_access():
            cur_records_tot = ch.get_random_access_event_count()

        # Push initial progress so the bar appears immediately.
        _flush()
        last_push = time.monotonic()

        while True:
            if ch.read() != dec.Status.success:
                break
            cur_records_done += 1
            if not ch.scan():
                continue
            if ch.get_event_type() != dec.EventType.Physics:
                continue

            for si in range(ch.get_n_events()):
                ch.select_event(si)
                info = ch.info()
                d_processed += 1

                if int(info.trigger_bits) & skip_mask:
                    continue

                fadc_evt = ch.fadc()
                sf: Dict[str, bool]  = {s: False for s in scint_names}
                ma: Dict[str, float] = {}

                # Build a roc_tag -> roc-object map once per event so the
                # scint and HyCal passes can both look up ROCs without
                # iterating fadc_evt twice.
                roc_by_tag: Dict[int, object] = {}
                for r in range(fadc_evt.nrocs):
                    rr = fadc_evt.roc(r)
                    roc_by_tag[int(rr.tag)] = rr

                if has_batch:
                    # ---------- Fast C++ batched path ----------
                    # One C call per slot returns a (MAX_CHANNELS,) float32
                    # array of (max-in-window − pedestal) values.  Python
                    # only does dict lookups; no per-channel numpy alloc.
                    scint_any = False
                    for rt, sn, cn, name in scint_locs:
                        rr = roc_by_tag.get(rt)
                        if rr is None:
                            continue
                        heights = rr.slot(sn).peak_in_window(s_lo, s_hi, N_PED)
                        if heights[cn] > scint_thr:
                            sf[name]  = True
                            scint_any = True

                    if require_scint and not scint_any:
                        continue

                    for rt, hy_locs in hycal_locs_by_roc.items():
                        rr = roc_by_tag.get(rt)
                        if rr is None:
                            continue
                        last_slot_num = -1
                        heights = None
                        for sn, cn, name in hy_locs:
                            if sn != last_slot_num:
                                heights = rr.slot(sn).peak_in_window(
                                    h_lo, h_hi, N_PED)
                                last_slot_num = sn
                            h = heights[cn]
                            if h > hycal_thr:
                                ma[name] = float(h)
                else:
                    # ---------- Pure-Python fallback path ------
                    # Used when prad2py hasn't been rebuilt with
                    # SlotData.peak_in_window.  ~10× slower per event,
                    # but functionally equivalent.
                    scint_any = False
                    for rt, sn, cn, name in scint_locs:
                        rr = roc_by_tag.get(rt)
                        if rr is None:
                            continue
                        chan = rr.slot(sn).channel(cn)
                        if chan.nsamples < N_PED + 1:
                            continue
                        samples = chan.samples
                        if samples.max() - samples.min() < scint_thr:
                            continue
                        ped    = samples[:N_PED].mean()
                        window = samples[s_lo:s_hi]
                        if window.size == 0:
                            continue
                        if window.max() - ped > scint_thr:
                            sf[name]  = True
                            scint_any = True

                    if require_scint and not scint_any:
                        continue

                    for rt, hy_locs in hycal_locs_by_roc.items():
                        rr = roc_by_tag.get(rt)
                        if rr is None:
                            continue
                        last_slot_num = -1
                        slot = None
                        for sn, cn, name in hy_locs:
                            if sn != last_slot_num:
                                slot = rr.slot(sn)
                                last_slot_num = sn
                            chan = slot.channel(cn)
                            if chan.nsamples < N_PED + 1:
                                continue
                            samples = chan.samples
                            if samples.max() - samples.min() < hycal_thr:
                                continue
                            ped    = samples[:N_PED].mean()
                            window = samples[h_lo:h_hi]
                            if window.size == 0:
                                continue
                            height = float(window.max() - ped)
                            if height > hycal_thr:
                                ma[name] = height

                n_ma = len(ma)
                if n_ma > lms_max or n_ma < min_mods:
                    continue

                if (max_lm > 0 and n_ma > 1
                        and _count_local_maxima(ma, neighbor_map) > max_lm):
                    continue

                best = max(ma, key=ma.__getitem__)

                for sname, fired in sf.items():
                    if fired:
                        d_scint_hits_any[sname] += 1

                any_scint = sf.get("V1") or sf.get("V2") or sf.get("V3") or sf.get("V4")
                if not any_scint and (mode == MODE_AND or scint_t_min > -math.inf):
                    continue

                d_module_hits[best] += 1
                for sname, fired in sf.items():
                    if fired:
                        d_scint_hits[sname] += 1
                        d_coincidences[sname][best] += 1

            # Time-based partial push (between records, not per event,
            # to keep time.monotonic() out of the hottest inner loop).
            now = time.monotonic()
            if now - last_push >= push_interval:
                _flush()
                last_push = now

        # Final push for this file so the parent sees the 100% mark.
        _flush()
        ch.close()

    finished_chunk = True
    _flush()
    return None


def _channel_fired(ch_data: dict, threshold: float,
                   t_min: float = -math.inf, t_max: float = math.inf) -> bool:
    """True if any peak has height > threshold and time within [t_min, t_max] ns."""
    for pk in ch_data.get("pk", []):
        if pk.get("h", 0.0) > threshold and t_min <= pk.get("t", 0.0) <= t_max:
            return True
    return False


# ===========================================================================
#  Statistics container
# ===========================================================================

class Stats:
    """Thread-safe coincidence accumulator."""

    def __init__(self, module_names):
        self._lock = threading.Lock()
        self.module_hits:    Dict[str, int] = {m: 0 for m in module_names}
        self.scint_hits:     Dict[str, int] = {s: 0 for s in SCINTILLATORS}
        self.scint_hits_any: Dict[str, int] = {s: 0 for s in SCINTILLATORS}
        self.coincidences: Dict[str, Dict[str, int]] = {
            s: {m: 0 for m in module_names} for s in SCINTILLATORS
        }
        self.processed = 0

    def update(self, scint_fired: Dict[str, bool],
               module_fired: Dict[str, bool]) -> None:
        """Accumulate one event.

        For every module Mj that fired:
            module_hits[Mj]          += 1
            coincidences[Vi][Mj]     += 1  for each Vi that also fired

        This gives  rate(Vi, Mj) = coincidences[Vi][Mj] / module_hits[Mj]
                                 = N(Vi fired AND Mj fired) / N(Mj fired)
        """
        with self._lock:
            for sname, sf in scint_fired.items():
                if sf:
                    self.scint_hits[sname] += 1
            for mname, mf in module_fired.items():
                if mf:
                    self.module_hits[mname] += 1
                    for sname, sf in scint_fired.items():
                        if sf:
                            self.coincidences[sname][mname] += 1
            self.processed += 1

    def update_scint_any(self, scint_fired: Dict[str, bool]) -> None:
        """Count per-scintillator fires for events that passed the cluster cut,
        regardless of AND/OR mode.  Call this before the AND filter."""
        with self._lock:
            for sname, sf in scint_fired.items():
                if sf:
                    self.scint_hits_any[sname] += 1

    def snapshot(self, mode: str = MODE_AND) -> dict:
        with self._lock:
            rates = {}
            for sname in SCINTILLATORS:
                rates[sname] = {}
                s_denom = self.scint_hits[sname]   # N(V_i fired AND HyCal cluster)
                for mname, m_denom in self.module_hits.items():
                    ncoinc = self.coincidences[sname][mname]
                    if mode == MODE_AND:
                        # P(V_i fired | M is best): how often V_i coincides with M
                        rates[sname][mname] = (ncoinc / m_denom
                                               if m_denom > 0 else math.nan)
                    else:
                        # P(M is best | V_i fired): spatial dist. of HyCal given V_i
                        rates[sname][mname] = (ncoinc / s_denom
                                               if s_denom > 0 else math.nan)
            return {
                "rates": rates,
                "mode": mode,
                "module_hits": dict(self.module_hits),
                "scint_hits": dict(self.scint_hits),
                "scint_hits_any": dict(self.scint_hits_any),
                "processed": self.processed,
            }


# ===========================================================================
#  Waveform collector
# ===========================================================================

class WaveformCollector:
    """Thread-safe accumulator of per-event waveform records for coincidences.

    For each qualifying event, stores ADC samples from the fired scintillator
    and the HyCal module with the highest ADC.  Samples are read from the
    inline 's' field if present (ET mode); otherwise a separate
    /api/waveform/{ev}/{key} request is made (file mode).

    Call save(path) to write a compressed .npz file when done.
    """

    def __init__(self, server_url: str, max_records: int) -> None:
        self._url  = server_url
        self._max  = max_records
        self._buf: List[dict] = []
        self._lock = threading.Lock()

    @property
    def count(self) -> int:
        with self._lock:
            return len(self._buf)

    @property
    def full(self) -> bool:
        with self._lock:
            return len(self._buf) >= self._max

    def _get_samples(self, ev_num: int, ch_data: dict, key: str) -> List[int]:
        sp = ch_data.get("s", [])
        if sp:
            return list(sp)
        wf = _http_get(f"{self._url}/api/waveform/{ev_num}/{key}")
        return wf.get("s", []) if wf and "error" not in wf else []

    def collect(self, data: dict,
                scint_name: str, scint_key: str,
                hycal_name: str, hycal_key: str) -> bool:
        """Try to add one record.  Returns True if a record was stored."""
        with self._lock:
            if len(self._buf) >= self._max:
                return False
        channels  = data.get("channels", {})
        ev_num    = data.get("event_number", data.get("event", 0))
        scint_ch  = channels.get(scint_key, {})
        hycal_ch  = channels.get(hycal_key, {})
        scint_sp  = self._get_samples(ev_num, scint_ch, scint_key)
        hycal_sp  = self._get_samples(ev_num, hycal_ch, hycal_key)
        scint_int = max((pk.get("i", 0.0) for pk in scint_ch.get("pk", [])), default=0.0)
        hycal_int = max((pk.get("i", 0.0) for pk in hycal_ch.get("pk", [])), default=0.0)
        with self._lock:
            if len(self._buf) >= self._max:
                return False
            self._buf.append({
                "event":         ev_num,
                "scint_name":    scint_name,
                "scint_key":     scint_key,
                "scint_int":     scint_int,
                "scint_pm":      float(scint_ch.get("pm", 0) or 0),
                "scint_samples": scint_sp,
                "hycal_name":    hycal_name,
                "hycal_key":     hycal_key,
                "hycal_int":     hycal_int,
                "hycal_pm":      float(hycal_ch.get("pm", 0) or 0),
                "hycal_samples": hycal_sp,
            })
            return True

    def save(self, path: Path) -> int:
        """Write accumulated records to a compressed .npz file.
        Returns the number of records written."""
        with self._lock:
            buf = list(self._buf)
        if not buf:
            return 0
        path.parent.mkdir(parents=True, exist_ok=True)
        n = len(buf)
        sp_len = max(
            max((len(e["scint_samples"]) for e in buf), default=0),
            max((len(e["hycal_samples"]) for e in buf), default=0),
            1,
        )
        scint_sp = np.zeros((n, sp_len), dtype=np.int16)
        hycal_sp = np.zeros((n, sp_len), dtype=np.int16)
        for i, e in enumerate(buf):
            ss, hs = e["scint_samples"], e["hycal_samples"]
            scint_sp[i, :len(ss)] = ss
            hycal_sp[i, :len(hs)] = hs
        np.savez_compressed(
            str(path),
            event_numbers   = np.array([e["event"]      for e in buf], dtype=np.int64),
            scint_names     = np.array([e["scint_name"]  for e in buf]),
            scint_keys      = np.array([e["scint_key"]   for e in buf]),
            scint_integrals = np.array([e["scint_int"]   for e in buf], dtype=np.float32),
            scint_ped_means = np.array([e["scint_pm"]    for e in buf], dtype=np.float32),
            scint_samples   = scint_sp,
            hycal_names     = np.array([e["hycal_name"]  for e in buf]),
            hycal_keys      = np.array([e["hycal_key"]   for e in buf]),
            hycal_integrals = np.array([e["hycal_int"]   for e in buf], dtype=np.float32),
            hycal_ped_means = np.array([e["hycal_pm"]    for e in buf], dtype=np.float32),
            hycal_samples   = hycal_sp,
        )
        return n


# ===========================================================================
#  Coincidence scan worker
# ===========================================================================

def _fetch_event(server_url: str, ev: int) -> Optional[dict]:
    return _http_get(f"{server_url}/api/event/{ev}")


def _collect_wfm(collector: "WaveformCollector", data: dict,
                 scint_fired: Dict[str, bool],
                 scint_keys: Dict[str, str],
                 hycal_name: str, hycal_key: str) -> None:
    """Pick the highest-ADC fired scintillator and hand off to the collector."""
    channels = data.get("channels", {})
    best_sname, best_skey, best_sint = None, None, -1.0
    for sname, skey in scint_keys.items():
        if not scint_fired.get(sname):
            continue
        integral = max(
            (pk.get("i", 0.0) for pk in channels.get(skey, {}).get("pk", [])),
            default=0.0,
        )
        if integral > best_sint:
            best_sint, best_sname, best_skey = integral, sname, skey
    if best_sname is None:
        return
    collector.collect(data,
                      scint_name=best_sname, scint_key=best_skey,
                      hycal_name=hycal_name, hycal_key=hycal_key)


class ProcessWorker(QThread):
    progress         = pyqtSignal(int, int)
    stats_update     = pyqtSignal(dict)
    finished         = pyqtSignal(str)
    waveforms_saved  = pyqtSignal(int, str)   # (count, file_path)

    def __init__(self, server_url: str, n_events: int,
                 module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 scint_thr: float, hycal_thr: float,
                 mode: str = MODE_AND,
                 min_mods: int = DEFAULT_MIN_CLUSTER_MODS,
                 scint_t_min: float = -math.inf,
                 scint_t_max: float = math.inf,
                 hycal_t_min: float = -math.inf,
                 hycal_t_max: float = math.inf,
                 neighbor_map: Optional[Dict[str, frozenset]] = None,
                 max_local_maxima: int = 0,
                 wfm_collector: Optional["WaveformCollector"] = None,
                 wfm_save_path: Optional[Path] = None,
                 parent=None):
        super().__init__(parent)
        self._url        = server_url
        self._n          = n_events
        self._mod_keys   = module_keys
        self._scint_keys = scint_keys
        self._scint_thr  = scint_thr
        self._hycal_thr  = hycal_thr
        self._mode       = mode
        self._min_mods   = min_mods
        self._scint_t_min  = scint_t_min
        self._scint_t_max  = scint_t_max
        self._hycal_t_min  = hycal_t_min
        self._hycal_t_max  = hycal_t_max
        self._neighbor_map = neighbor_map or {}
        self._max_lm       = max_local_maxima
        self._wfm_coll   = wfm_collector
        self._wfm_path   = wfm_save_path
        self._stop_evt   = threading.Event()

    def stop(self):
        self._stop_evt.set()

    def run(self):
        stats     = Stats(list(self._mod_keys.keys()))
        last_emit = time.monotonic()

        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            batch_start = 1
            while batch_start <= self._n and not self._stop_evt.is_set():
                batch_end = min(batch_start + BATCH_SIZE - 1, self._n)
                futures = {
                    pool.submit(_fetch_event, self._url, ev): ev
                    for ev in range(batch_start, batch_end + 1)
                }

                for fut in as_completed(futures):
                    if self._stop_evt.is_set():
                        break
                    data = fut.result()
                    if not data or "error" in data:
                        with stats._lock:
                            stats.processed += 1
                        continue

                    # Reject LMS / alpha events by trigger bit (fast path)
                    if data.get("trigger_bits", 0) & SKIP_TRIGGER_MASK:
                        with stats._lock:
                            stats.processed += 1
                        continue

                    channels = data.get("channels", {})
                    scint_fired = {
                        sname: (skey in channels
                                and _channel_fired(channels[skey], self._scint_thr,
                                               self._scint_t_min, self._scint_t_max))
                        for sname, skey in self._scint_keys.items()
                    }

                    # Compute per-module ADC (max peak height within time cut).
                    module_adc: Dict[str, float] = {}
                    for mname, mkey in self._mod_keys.items():
                        if mkey in channels:
                            peaks = channels[mkey].get("pk", [])
                            adc = max(
                                (float(pk.get("h", 0.0)) for pk in peaks
                                 if self._hycal_t_min <= pk.get("t", 0.0)
                                                      <= self._hycal_t_max),
                                default=0.0)
                            if adc > self._hycal_thr:
                                module_adc[mname] = adc

                    # Reject LMS events by module occupancy (>1000 modules fired)
                    if len(module_adc) > LMS_MAX_MODULES:
                        with stats._lock:
                            stats.processed += 1
                        continue

                    # Reject isolated single-module noise: require a cluster.
                    if len(module_adc) < self._min_mods:
                        with stats._lock:
                            stats.processed += 1
                        continue

                    # Reject multi-cluster events (>1 local ADC maximum).
                    if (self._max_lm > 0
                            and _count_local_maxima(module_adc, self._neighbor_map)
                            > self._max_lm):
                        with stats._lock:
                            stats.processed += 1
                        continue

                    # Only the module with the highest ADC gets event credit.
                    # The cluster cut above ensures it is a genuine cluster, not
                    # an isolated discharge.
                    best = max(module_adc, key=module_adc.get)
                    module_fired = {best: True}

                    # Track individual scintillator fires before the AND filter so
                    # the stats panel always shows per-scintillator rates.
                    stats.update_scint_any(scint_fired)

                    # Skip if no scint fired in-window (AND mode always; any
                    # mode when a time cut is active).
                    if not any(scint_fired.values()) and (
                            self._mode == MODE_AND
                            or self._scint_t_min > -math.inf):
                        with stats._lock:
                            stats.processed += 1
                        continue

                    stats.update(scint_fired, module_fired)

                    if self._wfm_coll and not self._wfm_coll.full:
                        _collect_wfm(self._wfm_coll, data, scint_fired,
                                     self._scint_keys, best, self._mod_keys[best])

                batch_start = batch_end + 1
                now = time.monotonic()
                if now - last_emit > 0.5:
                    snap = stats.snapshot(self._mode)
                    self.progress.emit(snap["processed"], self._n)
                    self.stats_update.emit(snap)
                    last_emit = now

        snap = stats.snapshot(self._mode)
        self.progress.emit(snap["processed"], self._n)
        self.stats_update.emit(snap)
        if self._wfm_coll and self._wfm_path and self._wfm_coll.count > 0:
            n = self._wfm_coll.save(self._wfm_path)
            self.waveforms_saved.emit(n, str(self._wfm_path))
        self.finished.emit("" if not self._stop_evt.is_set() else "stopped")


class ProcessWorkerET(QThread):
    """Accumulates coincidence statistics from live ET events.

    Polls /api/ring at ~5 Hz, fetches each new sequence number from the
    ring buffer, and processes it exactly once.  Runs until stopped.
    """
    progress         = pyqtSignal(int, int)   # (processed, -1)  — -1 signals ET mode
    stats_update     = pyqtSignal(dict)
    finished         = pyqtSignal(str)
    waveforms_saved  = pyqtSignal(int, str)   # (count, file_path)

    POLL_INTERVAL = 0.05  # seconds between /api/ring polls (20 Hz)

    def __init__(self, server_url: str,
                 module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 scint_thr: float, hycal_thr: float,
                 mode: str = MODE_AND,
                 min_mods: int = DEFAULT_MIN_CLUSTER_MODS,
                 max_rate_hz: float = 0.0,
                 scint_t_min: float = -math.inf,
                 scint_t_max: float = math.inf,
                 hycal_t_min: float = -math.inf,
                 hycal_t_max: float = math.inf,
                 neighbor_map: Optional[Dict[str, frozenset]] = None,
                 max_local_maxima: int = 0,
                 wfm_collector: Optional["WaveformCollector"] = None,
                 wfm_save_path: Optional[Path] = None,
                 parent=None):
        super().__init__(parent)
        self._url          = server_url
        self._mod_keys     = module_keys
        self._scint_keys   = scint_keys
        self._scint_thr    = scint_thr
        self._hycal_thr    = hycal_thr
        self._mode         = mode
        self._min_mods     = min_mods
        self._max_rate     = max_rate_hz
        self._scint_t_min  = scint_t_min
        self._scint_t_max  = scint_t_max
        self._hycal_t_min  = hycal_t_min
        self._hycal_t_max  = hycal_t_max
        self._neighbor_map = neighbor_map or {}
        self._max_lm       = max_local_maxima
        self._wfm_coll    = wfm_collector
        self._wfm_path    = wfm_save_path
        self._stop_evt    = threading.Event()

    def stop(self):
        self._stop_evt.set()

    def _process_event(self, data: dict, stats: "Stats") -> None:
        """Process one decoded event dict and accumulate into stats."""
        if data.get("trigger_bits", 0) & SKIP_TRIGGER_MASK:
            with stats._lock:
                stats.processed += 1
            return

        channels = data.get("channels", {})
        scint_fired = {
            sname: (skey in channels
                    and _channel_fired(channels[skey], self._scint_thr,
                                               self._scint_t_min, self._scint_t_max))
            for sname, skey in self._scint_keys.items()
        }

        module_adc: Dict[str, float] = {}
        for mname, mkey in self._mod_keys.items():
            if mkey in channels:
                peaks = channels[mkey].get("pk", [])
                adc = max(
                    (float(pk.get("h", 0.0)) for pk in peaks
                     if self._hycal_t_min <= pk.get("t", 0.0) <= self._hycal_t_max),
                    default=0.0)
                if adc > self._hycal_thr:
                    module_adc[mname] = adc

        if len(module_adc) > LMS_MAX_MODULES:
            with stats._lock:
                stats.processed += 1
            return

        if len(module_adc) < self._min_mods:
            with stats._lock:
                stats.processed += 1
            return

        if (self._max_lm > 0
                and _count_local_maxima(module_adc, self._neighbor_map)
                > self._max_lm):
            with stats._lock:
                stats.processed += 1
            return

        best = max(module_adc, key=module_adc.get)
        module_fired = {best: True}

        stats.update_scint_any(scint_fired)

        if not any(scint_fired.values()) and (
                self._mode == MODE_AND or self._scint_t_min > -math.inf):
            with stats._lock:
                stats.processed += 1
            return

        stats.update(scint_fired, module_fired)

        if self._wfm_coll and not self._wfm_coll.full:
            _collect_wfm(self._wfm_coll, data, scint_fired,
                         self._scint_keys, best, self._mod_keys[best])

    def run(self):
        stats     = Stats(list(self._mod_keys.keys()))
        last_seq  = 0    # max sequence number seen; avoids unbounded seen_seqs set
        last_emit = time.monotonic()

        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            while not self._stop_evt.is_set():
                batch_t0  = time.monotonic()
                ring_data = _http_get(f"{self._url}/api/ring", timeout=2.0)
                if ring_data is None:
                    time.sleep(self.POLL_INTERVAL)
                    continue

                # Only process sequence numbers we haven't seen yet.
                # Sequences are monotonically increasing so a single integer
                # suffices (no need for a growing seen_seqs set).
                new_seqs = sorted(
                    s for s in ring_data.get("ring", []) if s > last_seq)

                if new_seqs:
                    futures = {
                        pool.submit(_fetch_event, self._url, seq): seq
                        for seq in new_seqs
                        if not self._stop_evt.is_set()
                    }
                    for fut in as_completed(futures):
                        if self._stop_evt.is_set():
                            break
                        data = fut.result()
                        if not data or "error" in data:
                            with stats._lock:
                                stats.processed += 1
                        else:
                            self._process_event(data, stats)

                    last_seq = max(new_seqs)

                now = time.monotonic()
                if now - last_emit > 0.5:
                    snap = stats.snapshot(self._mode)
                    self.progress.emit(snap["processed"], -1)
                    self.stats_update.emit(snap)
                    last_emit = now

                # Rate limiting: sleep at least POLL_INTERVAL; sleep longer
                # if the user set a max rate and the batch finished too fast.
                batch_elapsed = time.monotonic() - batch_t0
                if self._max_rate > 0 and new_seqs:
                    target = len(new_seqs) / self._max_rate
                    sleep_t = max(self.POLL_INTERVAL, target - batch_elapsed)
                else:
                    sleep_t = self.POLL_INTERVAL
                time.sleep(sleep_t)

        snap = stats.snapshot(self._mode)
        self.progress.emit(snap["processed"], -1)
        self.stats_update.emit(snap)
        if self._wfm_coll and self._wfm_path and self._wfm_coll.count > 0:
            n = self._wfm_coll.save(self._wfm_path)
            self.waveforms_saved.emit(n, str(self._wfm_path))
        self.finished.emit("" if not self._stop_evt.is_set() else "stopped")


# ===========================================================================
#  Local EVIO worker  (no server — uses prad2py directly)
# ===========================================================================

class ProcessWorkerLocal(QThread):
    """Reads an EVIO file directly via prad2py and accumulates coincidence
    statistics.  Requires prad2py to be installed (_HAVE_PRAD2PY == True).

    Optimised for throughput: uses reverse key-lookup dicts so only channels
    that belong to a known scintillator or HyCal module are analysed; raw
    waveform samples are never copied to Python lists in the stats path.
    Local plain-Python counters replace the lock-based Stats object in the
    hot loop; a snapshot dict is built on demand for periodic UI updates.
    """
    progress          = pyqtSignal(int, int)   # (processed, total); total=0 → unknown
    stats_update      = pyqtSignal(dict)
    finished          = pyqtSignal(str)
    waveforms_saved   = pyqtSignal(int, str)
    coincidence_event = pyqtSignal(dict)       # step-through: per-event waveform data
    worker_progress   = pyqtSignal(dict)       # parallel mode: per-worker file progress
    workers_setup     = pyqtSignal(int)        # parallel mode: number of workers about to run

    def __init__(self, evio_paths: List[str],
                 module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 scint_thr: float, hycal_thr: float,
                 mode: str = MODE_AND,
                 min_mods: int = DEFAULT_MIN_CLUSTER_MODS,
                 scint_t_min: float = -math.inf,
                 scint_t_max: float = math.inf,
                 hycal_t_min: float = -math.inf,
                 hycal_t_max: float = math.inf,
                 neighbor_map: Optional[Dict[str, frozenset]] = None,
                 max_local_maxima: int = 0,
                 wfm_collector: Optional[WaveformCollector] = None,
                 wfm_save_path: Optional[Path] = None,
                 step_through: bool = False,
                 n_workers: int = 1,
                 parent=None):
        super().__init__(parent)
        self._paths        = list(evio_paths)
        self._mod_keys     = module_keys
        self._scint_keys   = scint_keys
        self._scint_thr    = scint_thr
        self._hycal_thr    = hycal_thr
        self._mode         = mode
        self._min_mods     = min_mods
        self._scint_t_min  = scint_t_min
        self._scint_t_max  = scint_t_max
        self._hycal_t_min  = hycal_t_min
        self._hycal_t_max  = hycal_t_max
        self._neighbor_map = neighbor_map or {}
        self._max_lm       = max_local_maxima
        self._wfm_coll     = wfm_collector
        self._wfm_path     = wfm_save_path
        self._step_through = step_through
        self._n_workers    = max(1, int(n_workers))
        self._stop_evt     = threading.Event()
        self._continue_evt = threading.Event()

        # Reverse lookups: channel-key → name, built once so the inner loop
        # does O(1) dict lookups instead of iterating over all module keys.
        self._scint_key_to_name: Dict[str, str] = {v: k for k, v in scint_keys.items()}
        self._mod_key_to_name:   Dict[str, str] = {v: k for k, v in module_keys.items()}
        self._all_keys = frozenset(scint_keys.values()) | frozenset(module_keys.values())
        # Tuple-key lookup: (roc_tag, slot, chan) → (kind, name)
        # Avoids f-string allocation per channel in the hot loop.
        self._tuple_to_role = _build_tuple_role_map(scint_keys, module_keys)

    def stop(self):
        self._stop_evt.set()
        self._continue_evt.set()   # unblock any pending step-through wait

    def step_continue(self):
        self._continue_evt.set()

    @staticmethod
    def _fetch_wfm_channels(fadc_evt, analyzer, needed_keys: frozenset) -> Dict[str, dict]:
        """Build a minimal channel dict (with raw samples) for a small set of keys.
        Called only for coincidence events when waveform saving is enabled."""
        channels: Dict[str, dict] = {}
        for r in range(fadc_evt.nrocs):
            roc     = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = f"{roc_tag}_{s}_{c}"
                    if key not in needed_keys:
                        continue
                    samples = slot.channel(c).samples
                    if samples.size < 4:
                        continue
                    ped_mean, ped_rms, peaks = analyzer.analyze(samples)
                    channels[key] = {
                        "pm": float(ped_mean),
                        "pr": float(ped_rms),
                        "s":  list(samples),
                        "pk": [{"i": float(p.integral), "h": float(p.height),
                                 "t": float(p.time),    "p": int(p.pos),
                                 "l": int(p.left),      "r": int(p.right),
                                 "o": int(p.overflow)}
                               for p in peaks],
                    }
        return channels

    def run(self):
        if not _HAVE_PRAD2PY:
            self.finished.emit("error: prad2py not available")
            return

        mod_names      = list(self._mod_keys.keys())
        processed      = 0
        scint_hits     = {s: 0 for s in SCINTILLATORS}
        scint_hits_any = {s: 0 for s in SCINTILLATORS}
        module_hits    = {m: 0 for m in mod_names}
        coincidences   = {s: {m: 0 for m in mod_names} for s in SCINTILLATORS}

        scint_thr     = self._scint_thr
        hycal_thr     = self._hycal_thr
        scint_t_min   = self._scint_t_min
        scint_t_max   = self._scint_t_max
        hycal_t_min   = self._hycal_t_min
        hycal_t_max   = self._hycal_t_max
        neighbor_map  = self._neighbor_map
        max_lm        = self._max_lm
        mode          = self._mode
        min_mods      = self._min_mods
        wfm_coll      = self._wfm_coll
        skip_mask     = SKIP_TRIGGER_MASK
        scint_names   = list(SCINTILLATORS)
        tuple_to_role = self._tuple_to_role

        def _make_snap() -> dict:
            rates: Dict[str, Dict[str, float]] = {}
            for sname in scint_names:
                rates[sname] = {}
                s_denom = scint_hits[sname]
                for mname in mod_names:
                    m_denom = module_hits[mname]
                    ncoinc  = coincidences[sname][mname]
                    if mode == MODE_AND:
                        rates[sname][mname] = ncoinc / m_denom if m_denom > 0 else math.nan
                    else:
                        rates[sname][mname] = ncoinc / s_denom if s_denom > 0 else math.nan
            return {
                "rates":          rates,
                "mode":           mode,
                "module_hits":    dict(module_hits),
                "scint_hits":     dict(scint_hits),
                "scint_hits_any": dict(scint_hits_any),
                "processed":      processed,
            }

        # ------------------------------------------------------------------
        # Parallel path: one subprocess per file.  Skips when waveform
        # collection or step-through is enabled (those need per-event UI).
        # ------------------------------------------------------------------
        if (len(self._paths) > 1
                and not self._wfm_coll
                and not self._step_through):
            n_files = len(self._paths)
            # Partition files into n_workers near-equal chunks.  With 10 files
            # and 3 workers → sizes [4, 3, 3]; the first ``extra`` chunks each
            # take one extra file.
            n_workers = max(1, min(self._n_workers, n_files))
            base, extra = divmod(n_files, n_workers)
            chunks: List[List[str]] = []
            start = 0
            for i in range(n_workers):
                size = base + (1 if i < extra else 0)
                if size > 0:
                    chunks.append(list(self._paths[start:start + size]))
                start += size

            # "fork" inherits Qt state and is fragile from a Qt thread.
            # "forkserver" forks from a clean intermediate process, which
            # is much more reliable for GUI apps.
            try:
                ctx = multiprocessing.get_context("forkserver")
            except ValueError:
                ctx = multiprocessing.get_context("fork")
            update_q = ctx.Queue()
            pool     = ctx.Pool(processes=len(chunks),
                                initializer=_init_worker,
                                initargs=(update_q,))

            # Stash worker errors so they can be surfaced as a finished("error: …")
            # message — silent crashes are how a "processing ends immediately"
            # bug usually looks.
            worker_errors: List[str] = []
            def _on_err(exc):
                import traceback
                msg = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
                sys.stderr.write(f"[coinc worker error]\n{msg}\n")
                sys.stderr.flush()
                worker_errors.append(repr(exc))

            file_args = [
                {
                    "paths":         chunk,
                    "tuple_to_role": tuple_to_role,
                    "scint_thr":     scint_thr,
                    "hycal_thr":     hycal_thr,
                    "scint_t_min":   scint_t_min,
                    "scint_t_max":   scint_t_max,
                    "hycal_t_min":   hycal_t_min,
                    "hycal_t_max":   hycal_t_max,
                    "neighbor_map":  neighbor_map,
                    "max_lm":        max_lm,
                    "mode":          mode,
                    "min_mods":      min_mods,
                    "scint_names":   scint_names,
                    "mod_names":     mod_names,
                    "skip_mask":     skip_mask,
                    "lms_max":       LMS_MAX_MODULES,
                    "push_interval": 0.5,
                    "worker_id":     i,
                }
                for i, chunk in enumerate(chunks)
            ]
            # Tell the UI how many workers to draw progress rows for.
            self.workers_setup.emit(len(chunks))
            async_results = [pool.apply_async(_process_files, (a,),
                                              error_callback=_on_err)
                             for a in file_args]
            pool.close()

            def _merge(msg):
                nonlocal processed
                for m, v in msg["module_hits"].items():
                    if v: module_hits[m] += v
                for s, v in msg["scint_hits"].items():
                    if v: scint_hits[s] += v
                for s, v in msg["scint_hits_any"].items():
                    if v: scint_hits_any[s] += v
                for s, ccol in msg["coincidences"].items():
                    row = coincidences[s]
                    for m, v in ccol.items():
                        if v: row[m] += v
                processed += msg["processed"]

            # Coalesce per-worker progress so we emit at most one signal per
            # worker per UI cycle even when the worker pushes more often.
            latest_progress: Dict[int, dict] = {}

            last_emit = time.monotonic()
            try:
                while not self._stop_evt.is_set():
                    drained = False
                    while True:
                        try:
                            msg = update_q.get(timeout=0.1)
                        except Exception:
                            break
                        _merge(msg)
                        if "worker_id" in msg:
                            latest_progress[msg["worker_id"]] = msg
                        drained = True

                    n_chunks_total = len(async_results)
                    n_done = sum(1 for ar in async_results if ar.ready())

                    now = time.monotonic()
                    if drained and now - last_emit > 0.5:
                        self.progress.emit(n_done, n_chunks_total)
                        self.stats_update.emit(_make_snap())
                        for prog in latest_progress.values():
                            self.worker_progress.emit(prog)
                        latest_progress.clear()
                        last_emit = now

                    if n_done >= n_chunks_total:
                        break
            finally:
                if self._stop_evt.is_set():
                    pool.terminate()

                # Wait for worker processes to fully exit.  Each worker's queue
                # feeder thread flushes pending messages on process exit, so
                # this is when the *final* per-worker _flush() lands in the
                # OS pipe buffer.  Without this wait, the parent would see
                # ar.ready()==True and proceed before the worker's last delta
                # message arrived, silently losing the tail of each chunk.
                pool.join()

                # Drain everything that's now sitting in the queue, including
                # the final flushes from every worker.
                while True:
                    try:
                        msg = update_q.get_nowait()
                    except Exception:
                        break
                    _merge(msg)
                    if "worker_id" in msg:
                        latest_progress[msg["worker_id"]] = msg
                for prog in latest_progress.values():
                    self.worker_progress.emit(prog)
                latest_progress.clear()

                # Surface any exceptions raised inside the workers.
                for ar in async_results:
                    if ar.ready():
                        try:
                            ar.get(timeout=0)
                        except Exception as e:
                            if not worker_errors:
                                _on_err(e)

            snap = _make_snap()
            n_chunks_total = len(async_results)
            self.progress.emit(n_chunks_total, n_chunks_total)
            self.stats_update.emit(snap)
            if worker_errors and not self._stop_evt.is_set():
                self.finished.emit(f"error: worker failed — {worker_errors[0]}")
            else:
                self.finished.emit("" if not self._stop_evt.is_set() else "stopped")
            return

        # ------------------------------------------------------------------
        # Sequential path: single file, or wfm/step-through requested.
        # ------------------------------------------------------------------
        dec      = _prad2py.dec
        cfg      = dec.load_daq_config()
        ch       = dec.EvChannel()
        ch.set_config(cfg)
        analyzer = dec.WaveAnalyzer(cfg.wave_cfg)

        # Pre-scan to get total EVIO record count for the progress bar.
        total_records = 0
        for path in self._paths:
            if ch.open_auto(path) == dec.Status.success:
                if ch.is_random_access():
                    total_records += ch.get_random_access_event_count()
                else:
                    total_records = 0
                    ch.close()
                    break
                ch.close()
            else:
                ch.close()

        last_emit = time.monotonic()
        n_records = 0

        for path in self._paths:
            if self._stop_evt.is_set():
                break

            if ch.open_auto(path) != dec.Status.success:
                continue

            while not self._stop_evt.is_set():
                if ch.read() != dec.Status.success:
                    break
                n_records += 1
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for si in range(ch.get_n_events()):
                    ch.select_event(si)
                    info = ch.info()
                    processed += 1

                    if int(info.trigger_bits) & skip_mask:
                        continue

                    fadc_evt = ch.fadc()
                    sf: Dict[str, bool]  = {}
                    ma: Dict[str, float] = {}

                    for r in range(fadc_evt.nrocs):
                        roc     = fadc_evt.roc(r)
                        roc_tag = int(roc.tag)
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                role = tuple_to_role.get((roc_tag, s, c))
                                if role is None:
                                    continue
                                samples = slot.channel(c).samples
                                if samples.size < 4:
                                    continue
                                # Cheap pre-filter — peak height above
                                # pedestal is bounded by max - min, so if
                                # that is below threshold we can skip the
                                # full waveform analysis entirely.
                                kind, name = role
                                thr_quick = scint_thr if kind == 0 else hycal_thr
                                if samples.max() - samples.min() < thr_quick:
                                    continue
                                _, _, peaks = analyzer.analyze(samples)
                                if kind == 0:
                                    sf[name] = any(
                                        p.height > scint_thr
                                        and scint_t_min <= p.time <= scint_t_max
                                        for p in peaks
                                    )
                                else:
                                    height = max(
                                        (p.height for p in peaks
                                         if hycal_t_min <= p.time <= hycal_t_max),
                                        default=0.0)
                                    if height > hycal_thr:
                                        ma[name] = height

                    for sname in scint_names:
                        if sname not in sf:
                            sf[sname] = False

                    n_ma = len(ma)
                    if n_ma > LMS_MAX_MODULES or n_ma < min_mods:
                        continue

                    # Single-module clusters are trivially 1 local max — skip the scan.
                    if (max_lm > 0 and n_ma > 1
                            and _count_local_maxima(ma, neighbor_map) > max_lm):
                        continue

                    best = max(ma, key=ma.__getitem__)

                    for sname, fired in sf.items():
                        if fired:
                            scint_hits_any[sname] += 1

                    any_scint = sf["V1"] or sf["V2"] or sf["V3"] or sf["V4"]
                    if not any_scint and (mode == MODE_AND
                                          or scint_t_min > -math.inf):
                        continue

                    module_hits[best] += 1
                    for sname, fired in sf.items():
                        if fired:
                            scint_hits[sname] += 1
                            coincidences[sname][best] += 1

                    if wfm_coll and not wfm_coll.full and any_scint:
                        best_key   = self._mod_keys[best]
                        wfm_needed = frozenset(self._scint_keys.values()) | {best_key}
                        wfm_ch     = self._fetch_wfm_channels(fadc_evt, analyzer, wfm_needed)
                        data_stub  = {
                            "trigger_bits": int(info.trigger_bits),
                            "event_number": int(info.event_number),
                            "event":        int(info.event_number),
                            "channels":     wfm_ch,
                        }
                        _collect_wfm(wfm_coll, data_stub, sf,
                                     self._scint_keys, best, best_key)

                    if self._step_through and any_scint:
                        best_key   = self._mod_keys[best]
                        wfm_needed = frozenset(self._scint_keys.values()) | {best_key}
                        wfm_ch     = self._fetch_wfm_channels(fadc_evt, analyzer, wfm_needed)
                        snap = _make_snap()
                        self.stats_update.emit(snap)
                        self.progress.emit(n_records, total_records)
                        self._continue_evt.clear()
                        self.coincidence_event.emit({
                            "event_number": int(info.event_number),
                            "best":         best,
                            "best_adc":     ma[best],
                            "scint_fired":  dict(sf),
                            "wfm_channels": wfm_ch,
                        })
                        while not self._stop_evt.is_set():
                            if self._continue_evt.wait(timeout=0.05):
                                break
                        if self._stop_evt.is_set():
                            break

                now = time.monotonic()
                if now - last_emit > 0.5:
                    snap = _make_snap()
                    self.progress.emit(n_records, total_records)
                    self.stats_update.emit(snap)
                    last_emit = now

            ch.close()

        snap = _make_snap()
        self.progress.emit(total_records if total_records else n_records,
                           total_records)
        self.stats_update.emit(snap)
        if self._wfm_coll and self._wfm_path and self._wfm_coll.count > 0:
            n = self._wfm_coll.save(self._wfm_path)
            self.waveforms_saved.emit(n, str(self._wfm_path))
        self.finished.emit("" if not self._stop_evt.is_set() else "stopped")


class InstantDisplayWorker(QThread):
    """Live event-by-event display from the ET ring buffer.

    Polls /api/ring at ~10 Hz, fetches each new latest event, computes
    per-module ADC values (max peak integral above pedestal), and emits
    event_ready with the raw channels dict so the main thread can render
    both the map and the waveform panels.
    """
    event_ready = pyqtSignal(dict)   # {adc_vals, channels, seq, max_mod}
    finished    = pyqtSignal(str)

    POLL_INTERVAL = 0.1   # 10 Hz

    def __init__(self, server_url: str,
                 module_keys: Dict[str, str],
                 hycal_thr: float = 0.0,
                 min_mods: int = DEFAULT_MIN_CLUSTER_MODS,
                 parent=None):
        super().__init__(parent)
        self._url       = server_url
        self._mod_keys  = module_keys
        self._hycal_thr = hycal_thr
        self._min_mods  = min_mods
        self._stop_evt  = threading.Event()

    def stop(self):
        self._stop_evt.set()

    def run(self):
        last_seq = 0

        while not self._stop_evt.is_set():
            ring_data = _http_get(f"{self._url}/api/ring", timeout=2.0)
            if ring_data is None:
                time.sleep(self.POLL_INTERVAL)
                continue

            latest = ring_data.get("latest", 0)
            if latest == 0 or latest <= last_seq:
                time.sleep(self.POLL_INTERVAL)
                continue

            last_seq = latest
            data = _http_get(f"{self._url}/api/event/{latest}")
            if not data or "error" in data:
                time.sleep(self.POLL_INTERVAL)
                continue

            # Reject LMS / alpha events by trigger bit.
            if data.get("trigger_bits", 0) & SKIP_TRIGGER_MASK:
                time.sleep(self.POLL_INTERVAL)
                continue

            channels = data.get("channels", {})

            # Per-module ADC = max peak height above pedestal, zeroed if below threshold.
            adc_vals: Dict[str, float] = {}
            for mname, mkey in self._mod_keys.items():
                ch = channels.get(mkey, {})
                peaks = ch.get("pk", [])
                val = max((float(pk.get("h", 0.0)) for pk in peaks), default=0.0)
                adc_vals[mname] = val if val > self._hycal_thr else 0.0

            # Skip events with too few modules (noise) or too many (LMS occupancy).
            n_above = sum(1 for v in adc_vals.values() if v > 0.0)
            if n_above < max(self._min_mods, 1) or n_above > LMS_MAX_MODULES:
                time.sleep(self.POLL_INTERVAL)
                continue

            max_mod = max(adc_vals, key=adc_vals.get) if adc_vals else ""

            self.event_ready.emit({
                "adc_vals":  adc_vals,
                "channels":  channels,
                "seq":       latest,
                "max_mod":   max_mod,
            })

            time.sleep(self.POLL_INTERVAL)

        self.finished.emit("" if not self._stop_evt.is_set() else "stopped")


# ===========================================================================
#  Waveform fetcher
# ===========================================================================

class WaveformFetcher(QThread):
    """Fetches scintillator and HyCal module waveforms for one event.

    File mode: calls /api/waveform/{ev}/{key} for each channel.
    ET mode:   calls /api/event/latest and extracts channel data from
               the full event JSON (ring events include samples).
    """
    waveform_ready = pyqtSignal(str, dict)   # label ("scint"|"module"), data

    def __init__(self, url: str, event_n: int,
                 scint_keys: Dict[str, str], module_key: str,
                 et_mode: bool = False,
                 parent=None):
        super().__init__(parent)
        self._url         = url
        self._ev          = event_n
        self._scint_keys  = scint_keys   # name → channel key, all 4 scints
        self._module_key  = module_key
        self._et_mode     = et_mode

    def run(self):
        if self._et_mode:
            event_data = _http_get(f"{self._url}/api/event/latest")
            if not event_data or "error" in event_data:
                err = {"error": "no ET event available"}
                self.waveform_ready.emit("scint",  {n: err for n in self._scint_keys})
                self.waveform_ready.emit("module", err)
                return
            channels = event_data.get("channels", {})
            scint_wfms = {}
            for name, key in self._scint_keys.items():
                if not key:
                    scint_wfms[name] = {"error": "no channel"}
                else:
                    ch = channels.get(key)
                    scint_wfms[name] = ch if ch else {"error": f"channel {key} not in event"}
            self.waveform_ready.emit("scint", scint_wfms)
            if not self._module_key:
                self.waveform_ready.emit("module", {"error": "no channel"})
            else:
                ch = channels.get(self._module_key)
                self.waveform_ready.emit(
                    "module",
                    ch if ch else {"error": f"channel {self._module_key} not in event"})
        else:
            scint_wfms = {}
            for name, key in self._scint_keys.items():
                if not key:
                    scint_wfms[name] = {"error": "no channel"}
                else:
                    data = _http_get(f"{self._url}/api/waveform/{self._ev}/{key}")
                    scint_wfms[name] = data or {"error": "no response"}
            self.waveform_ready.emit("scint", scint_wfms)
            if not self._module_key:
                self.waveform_ready.emit("module", {"error": "no channel"})
            else:
                data = _http_get(f"{self._url}/api/waveform/{self._ev}/{self._module_key}")
                self.waveform_ready.emit("module", data or {"error": "no response"})


# ===========================================================================
#  Waveform plot widget
# ===========================================================================

class WavePanel(QWidget):
    """Simple FADC waveform display driven by the server's waveform JSON.

    Expected input (from /api/waveform/<n>/<key>):
        {"s": [int, ...], "pm": float, "pr": float,
         "pk": [{"p": int, "h": float, "i": float,
                 "l": int, "r": int, "t": float, "o": int}, ...]}
    """

    PAD_L, PAD_R, PAD_T, PAD_B = 52, 14, 28, 32

    _PEAK_COLORS = (
        "#00b4d8", "#ff6b6b", "#51cf66", "#ffd43b",
        "#cc5de8", "#ff922b", "#20c997", "#f06595",
    )

    def __init__(self, label: str = "", parent=None):
        super().__init__(parent)
        self._label       = label
        self._title       = label
        self._samples: List[int] = []
        self._peaks:   List[dict] = []
        self._ped_mean    = 0.0
        self._ped_rms     = 0.0
        self._threshold   = 0.0
        self._t_min       = -math.inf
        self._t_max       = math.inf
        self._fired       = False
        self._placeholder = "No data — select an event and click Fetch"

        self.setMinimumHeight(120)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def set_data(self, wave_json: dict, threshold: float,
                 title: Optional[str] = None,
                 t_min: float = -math.inf,
                 t_max: float = math.inf) -> None:
        """Load waveform from the server JSON response."""
        if "error" in wave_json:
            self.clear(title or self._label,
                       wave_json.get("error", "channel not found"))
            return

        self._title     = title or self._label
        self._samples   = list(wave_json.get("s", []))
        self._peaks     = list(wave_json.get("pk", []))
        self._ped_mean  = float(wave_json.get("pm", 0))
        self._ped_rms   = float(wave_json.get("pr", 0))
        self._threshold = threshold
        self._t_min     = t_min
        self._t_max     = t_max
        self._fired     = any(
            pk.get("h", 0) > threshold
            and t_min <= pk.get("t", 0.0) <= t_max
            for pk in self._peaks
        )
        self.update()

    def clear(self, title: Optional[str] = None,
              placeholder: Optional[str] = None) -> None:
        self._title     = title or self._label
        self._samples   = []
        self._peaks     = []
        self._fired     = False
        self._t_min     = -math.inf
        self._t_max     = math.inf
        if placeholder is not None:
            self._placeholder = placeholder
        self.update()

    # ------------------------------------------------------------------
    #  Painting
    # ------------------------------------------------------------------

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = QRectF(self.PAD_L, self.PAD_T,
                   max(1.0, self.width() - self.PAD_L - self.PAD_R),
                   max(1.0, self.height() - self.PAD_T - self.PAD_B))
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        self._draw_title(p, r)

        n = len(self._samples)
        if n < 2:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, self._placeholder)
            return

        ymin, ymax = self._y_range()

        def sx(i: float) -> float:
            return r.left() + i / (n - 1) * r.width()

        def sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        y_ped = sy(self._ped_mean) if self._ped_mean != 0 else None
        self._draw_pedestal(p, r, y_ped)
        self._draw_peak_fills(p, sx, sy, y_ped, n)
        self._draw_waveform(p, sx, sy, n)
        self._draw_peak_markers(p, sx, sy, n)
        self._draw_axes(p, r, ymin, ymax, n, sx)
        self._draw_info(p, r)

    def _y_range(self):
        s = self._samples
        ymin, ymax = float(min(s)), float(max(s))
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad = (ymax - ymin) * 0.06
        return ymin - pad, ymax + pad

    def _draw_title(self, p: QPainter, r: QRectF):
        f = QFont("Monospace", 10)
        f.setBold(True)
        p.setFont(f)
        p.setPen(QColor(THEME.TEXT))
        p.drawText(int(r.left()), int(r.top() - 8), self._title)
        if self._samples:
            fm = p.fontMetrics()
            tw = fm.horizontalAdvance(self._title)
            fired_txt = "  [FIRED]" if self._fired else "  [—]"
            p.setPen(QColor(THEME.SUCCESS) if self._fired else QColor(THEME.TEXT_DIM))
            p.drawText(int(r.left() + tw), int(r.top() - 8), fired_txt)

    def _draw_pedestal(self, p: QPainter, r: QRectF, y_ped: Optional[float]):
        if y_ped is None:
            return
        p.setPen(QPen(QColor(THEME.TEXT_DIM), 1, Qt.PenStyle.DashLine))
        p.drawLine(int(r.left()), int(y_ped), int(r.right()), int(y_ped))

    def _draw_peak_fills(self, p, sx, sy, y_ped, n):
        if y_ped is None:
            return
        for i, pk in enumerate(self._peaks):
            base  = QColor(self._PEAK_COLORS[i % len(self._PEAK_COLORS)])
            # dim peaks that fail height threshold or fall outside the time window
            if (pk.get("h", 0) <= self._threshold
                    or not (self._t_min <= pk.get("t", 0.0) <= self._t_max)):
                base.setAlphaF(0.4)
            fill = QColor(base)
            fill.setAlphaF(fill.alphaF() * 0.25)
            lft = max(0, int(pk.get("l", pk["p"])))
            rgt = min(n - 1, int(pk.get("r", pk["p"])))
            poly = QPolygonF()
            for k in range(lft, rgt + 1):
                poly.append(QPointF(sx(k), sy(self._samples[k])))
            poly.append(QPointF(sx(rgt), y_ped))
            poly.append(QPointF(sx(lft), y_ped))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(fill)
            p.drawPolygon(poly)

    def _draw_waveform(self, p, sx, sy, n):
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        p.setBrush(Qt.BrushStyle.NoBrush)
        s = self._samples
        for i in range(n - 1):
            p.drawLine(int(sx(i)),     int(sy(s[i])),
                       int(sx(i + 1)), int(sy(s[i + 1])))

    def _draw_peak_markers(self, p, sx, sy, n):
        s = self._samples
        for i, pk in enumerate(self._peaks):
            pos = int(pk.get("p", 0))
            if pos < 0 or pos >= n:
                continue
            col = QColor(self._PEAK_COLORS[i % len(self._PEAK_COLORS)])
            if (pk.get("h", 0) <= self._threshold
                    or not (self._t_min <= pk.get("t", 0.0) <= self._t_max)):
                col.setAlphaF(0.4)
            p.setPen(QPen(col, 1.2))
            p.setBrush(col)
            cx, cy = sx(pos), sy(float(s[pos]))
            diamond = QPolygonF([
                QPointF(cx,     cy - 4),
                QPointF(cx + 4, cy),
                QPointF(cx,     cy + 4),
                QPointF(cx - 4, cy),
            ])
            p.drawPolygon(diamond)

            # height label above the diamond
            ht = pk.get("h", 0)
            p.setFont(QFont("Monospace", 8))
            p.setPen(col)
            p.drawText(int(cx - 16), int(cy - 7), f"{ht:.0f}")

    def _draw_axes(self, p, r, ymin, ymax, n, sx):
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4), f"{val:.0f}")
        tick_every = max(1, n // 8)
        for i in range(0, n, tick_every):
            x = r.left() + i / max(1, n - 1) * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            ns_val = i * 1e3 / CLK_MHZ
            p.drawText(int(x - 18), int(r.bottom() + 14), f"{ns_val:g}")
        p.setFont(QFont("Monospace", 9))
        p.drawText(int(r.left() + r.width() / 2 - 10),
                   int(r.bottom() + 26), "ns")

    def _draw_info(self, p, r):
        above = sum(1 for pk in self._peaks
                    if pk.get("h", 0) > self._threshold
                    and self._t_min <= pk.get("t", 0.0) <= self._t_max)
        info = (f"ped={self._ped_mean:.1f}  rms={self._ped_rms:.2f}"
                f"  peaks={len(self._peaks)} ({above} above thr)")
        p.setFont(QFont("Monospace", 9))
        fm = p.fontMetrics()
        tw = fm.horizontalAdvance(info)
        th = fm.height()
        box = QRectF(r.right() - tw - 8, r.top() + 4, tw + 6, th + 2)
        bg = QColor(THEME.BG)
        bg.setAlphaF(0.75)
        p.fillRect(box, bg)
        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(box, Qt.AlignmentFlag.AlignCenter, info)


# ===========================================================================
#  Multi-channel waveform panel (all scintillators overlaid)
# ===========================================================================

class MultiWavePanel(QWidget):
    """FADC waveform display: multiple channels overlaid with a colour legend.

    Primary use: show all scintillators (V1-V4) on one shared axis so the
    operator can compare timing and amplitude across all channels at once.

    API mirrors WavePanel for backward compat:
      set_multi_data(channels, threshold, scint_order, title, t_min, t_max)
      set_data(wave_json, threshold, title, t_min, t_max)   # single-channel compat
      clear(title, placeholder)
    """

    PAD_L, PAD_R, PAD_T, PAD_B = 52, 14, 28, 32

    _CHAN_COLORS = ["#4a9eff", "#ff6b6b", "#51cf66", "#ffd43b",
                    "#cc5de8", "#ff922b", "#20c997", "#f06595"]

    def __init__(self, label: str = "", parent=None):
        super().__init__(parent)
        self._label       = label
        self._title       = label
        self._channels: List[dict] = []
        self._threshold   = 0.0
        self._t_min       = -math.inf
        self._t_max       = math.inf
        self._placeholder = "No data — select an event and click Fetch"
        self.setMinimumHeight(120)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def set_multi_data(self, channels: dict, threshold: float,
                       scint_order: list,
                       title: Optional[str] = None,
                       t_min: float = -math.inf,
                       t_max: float = math.inf) -> None:
        """Overlay multiple waveforms. *channels* maps name → wave_json dict."""
        self._title     = title or self._label
        self._threshold = threshold
        self._t_min     = t_min
        self._t_max     = t_max
        self._channels  = []
        for i, name in enumerate(scint_order):
            wj    = channels.get(name, {})
            color = self._CHAN_COLORS[i % len(self._CHAN_COLORS)]
            if "error" in wj or "s" not in wj:
                self._channels.append({
                    "name": name, "color": color,
                    "samples": [], "peaks": [],
                    "ped_mean": 0.0, "ped_rms": 0.0, "fired": False,
                })
            else:
                samples = list(wj.get("s", []))
                peaks   = list(wj.get("pk", []))
                fired   = any(
                    pk.get("h", 0) > threshold
                    and t_min <= pk.get("t", 0.0) <= t_max
                    for pk in peaks
                )
                self._channels.append({
                    "name": name, "color": color,
                    "samples": samples, "peaks": peaks,
                    "ped_mean": float(wj.get("pm", 0)),
                    "ped_rms":  float(wj.get("pr", 0)),
                    "fired":    fired,
                })
        self.update()

    def set_data(self, wave_json: dict, threshold: float,
                 title: Optional[str] = None,
                 t_min: float = -math.inf,
                 t_max: float = math.inf) -> None:
        """Backward-compat: display a single channel (no multi-overlay)."""
        self.set_multi_data({self._label: wave_json}, threshold,
                            [self._label], title, t_min, t_max)

    def clear(self, title: Optional[str] = None,
              placeholder: Optional[str] = None) -> None:
        self._title    = title or self._label
        self._channels = []
        if placeholder is not None:
            self._placeholder = placeholder
        self.update()

    # ------------------------------------------------------------------
    #  Painting
    # ------------------------------------------------------------------

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = QRectF(self.PAD_L, self.PAD_T,
                   max(1.0, self.width() - self.PAD_L - self.PAD_R),
                   max(1.0, self.height() - self.PAD_T - self.PAD_B))
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        self._draw_title(p, r)

        all_samples = [ch["samples"] for ch in self._channels if ch["samples"]]
        if not all_samples:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, self._placeholder)
            self._draw_legend(p, r)
            return

        n = max(len(s) for s in all_samples)
        ymin, ymax = self._y_range(all_samples)

        def sx(i: float) -> float:
            return r.left() + i / max(1, n - 1) * r.width()

        def sy(v: float) -> float:
            return r.bottom() - (v - ymin) / max(1e-6, ymax - ymin) * r.height()

        for ch in self._channels:
            if ch["samples"]:
                self._draw_peak_fills(p, sx, sy, ch)
        for ch in self._channels:
            if ch["samples"]:
                self._draw_waveform(p, sx, sy, ch)
        for ch in self._channels:
            if ch["samples"]:
                self._draw_peak_markers(p, sx, sy, ch)
        self._draw_axes(p, r, ymin, ymax, n, sx)
        self._draw_legend(p, r)

    def _y_range(self, all_samples):
        flat = [v for s in all_samples for v in s]
        ymin, ymax = float(min(flat)), float(max(flat))
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad = (ymax - ymin) * 0.06
        return ymin - pad, ymax + pad

    def _draw_title(self, p: QPainter, r: QRectF):
        f = QFont("Monospace", 10)
        f.setBold(True)
        p.setFont(f)
        p.setPen(QColor(THEME.TEXT))
        p.drawText(int(r.left()), int(r.top() - 8), self._title)

    def _draw_peak_fills(self, p, sx, sy, ch):
        ped = ch["ped_mean"]
        if not ped:
            return
        y_ped    = sy(ped)
        samples  = ch["samples"]
        cn       = len(samples)
        for pk in ch["peaks"]:
            base = QColor(ch["color"])
            if (pk.get("h", 0) <= self._threshold
                    or not (self._t_min <= pk.get("t", 0.0) <= self._t_max)):
                base.setAlphaF(0.3)
            fill = QColor(base)
            fill.setAlphaF(fill.alphaF() * 0.25)
            lft  = max(0, int(pk.get("l", pk["p"])))
            rgt  = min(cn - 1, int(pk.get("r", pk["p"])))
            poly = QPolygonF()
            for k in range(lft, rgt + 1):
                poly.append(QPointF(sx(k), sy(samples[k])))
            poly.append(QPointF(sx(rgt), y_ped))
            poly.append(QPointF(sx(lft), y_ped))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(fill)
            p.drawPolygon(poly)

    def _draw_waveform(self, p, sx, sy, ch):
        samples = ch["samples"]
        cn      = len(samples)
        p.setPen(QPen(QColor(ch["color"]), 1.4))
        p.setBrush(Qt.BrushStyle.NoBrush)
        for i in range(cn - 1):
            p.drawLine(int(sx(i)),     int(sy(samples[i])),
                       int(sx(i + 1)), int(sy(samples[i + 1])))

    def _draw_peak_markers(self, p, sx, sy, ch):
        samples = ch["samples"]
        cn      = len(samples)
        for pk in ch["peaks"]:
            pos = int(pk.get("p", 0))
            if pos < 0 or pos >= cn:
                continue
            col = QColor(ch["color"])
            if (pk.get("h", 0) <= self._threshold
                    or not (self._t_min <= pk.get("t", 0.0) <= self._t_max)):
                col.setAlphaF(0.4)
            p.setPen(QPen(col, 1.2))
            p.setBrush(col)
            cx, cy = sx(pos), sy(float(samples[pos]))
            diamond = QPolygonF([
                QPointF(cx,     cy - 4),
                QPointF(cx + 4, cy),
                QPointF(cx,     cy + 4),
                QPointF(cx - 4, cy),
            ])
            p.drawPolygon(diamond)
            p.setFont(QFont("Monospace", 8))
            p.setPen(col)
            p.drawText(int(cx - 16), int(cy - 7), f"{pk.get('h', 0):.0f}")

    def _draw_axes(self, p, r, ymin, ymax, n, sx):
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4), f"{val:.0f}")
        tick_every = max(1, n // 8)
        for i in range(0, n, tick_every):
            x = r.left() + i / max(1, n - 1) * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            ns_val = i * 1e3 / CLK_MHZ
            p.drawText(int(x - 18), int(r.bottom() + 14), f"{ns_val:g}")
        p.setFont(QFont("Monospace", 9))
        p.drawText(int(r.left() + r.width() / 2 - 10),
                   int(r.bottom() + 26), "ns")

    def _draw_legend(self, p: QPainter, r: QRectF):
        if not self._channels:
            return
        f = QFont("Monospace", 9)
        p.setFont(f)
        fm = p.fontMetrics()
        lh = fm.height() + 3

        entries = []
        for ch in self._channels:
            peaks_in = [pk for pk in ch.get("peaks", [])
                        if pk.get("h", 0) > self._threshold
                        and self._t_min <= pk.get("t", 0.0) <= self._t_max]
            t_str     = f" t={peaks_in[0].get('t', 0):.0f}ns" if peaks_in else ""
            fired_txt = "[FIRED]" if ch["fired"] else "[—]"
            label     = f"{ch['name']} {fired_txt}{t_str}"
            entries.append((label, ch["color"], ch["fired"]))

        max_w   = max(fm.horizontalAdvance(e[0]) for e in entries) + 20
        total_h = lh * len(entries) + 8

        lx = int(r.right() - max_w - 6)
        ly = int(r.top() + 4)

        bg = QColor(THEME.BG)
        bg.setAlphaF(0.82)
        p.fillRect(QRectF(lx - 2, ly, max_w + 4, total_h), bg)

        y = ly + 4
        for label, color, fired in entries:
            col = QColor(color)
            p.fillRect(QRectF(lx, y + 1, 12, lh - 2), col)
            p.setPen(QColor(THEME.SUCCESS) if fired else QColor(THEME.TEXT_DIM))
            p.setFont(f)
            p.drawText(int(lx + 16), int(y + fm.ascent()), label)
            y += lh


# ===========================================================================
#  Coincidence map widget
# ===========================================================================

class CoincidenceMapWidget(HyCalMapWidget):
    """HyCal map coloured by coincidence rate, with informative tooltip."""

    def __init__(self, parent=None):
        super().__init__(parent, enable_zoom_pan=True, show_colorbar=True)
        self._snapshot: dict = {}
        self._active_scint = "V1"
        self._view_mode = VIEW_COINC
        self._instant_max_mod = ""
        self.set_palette("rainbow")

    def set_snapshot(self, snap: dict, scint: str,
                     view_mode: str = VIEW_COINC) -> None:
        self._snapshot    = snap
        self._active_scint = scint
        self._view_mode   = view_mode

        if view_mode == VIEW_OCC:
            hits = snap.get("module_hits", {})
            self.set_values(hits if hits else {})
            max_hits = max(hits.values()) if hits else 1
            self.set_range(0, max_hits)
        else:
            rates = snap.get("rates", {}).get(scint, {})
            valid = {m: v for m, v in rates.items() if not math.isnan(v)}
            self.set_values(valid if valid else {})
            self.set_range(0.0, 1.0)
        self.update()

    def set_instant_event(self, adc_vals: dict, max_mod: str = "") -> None:
        """Display per-module ADC values for a single live event."""
        self._view_mode = VIEW_INSTANT
        self._instant_max_mod = max_mod
        self.set_values(adc_vals)
        max_val = max(adc_vals.values(), default=0.0)
        self.set_range(0.0, max(max_val, 1.0))
        self.update()

    # Dark purple used for modules that have been seen (module_hits > 0) but
    # have zero coincidence rate or zero occupancy count.  Distinct from the
    # NO_DATA_COLOR (dark grey) which means the module was never seen at all.
    _SEEN_ZERO = QColor(55, 0, 75)

    def _paint_modules(self, p: QPainter):
        if self._view_mode not in (VIEW_COINC, VIEW_OCC):
            super()._paint_modules(p)
            return
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        seen_zero = self._SEEN_ZERO
        vmin, vmax = self._vmin, self._vmax
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                p.fillRect(rect, no_data)
            elif v == 0.0:
                p.fillRect(rect, seen_zero)
            else:
                t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                p.fillRect(rect, cmap_qcolor(max(t, 0.0), stops))

    def _tooltip_text(self, name: str) -> str:
        if self._view_mode == VIEW_INSTANT:
            adc = self._values.get(name, 0.0)
            suffix = "  ★ max" if name == self._instant_max_mod else ""
            return f"{name}\nADC: {adc:.0f}{suffix}"

        snap  = self._snapshot
        scint = self._active_scint
        nhits = snap.get("module_hits", {}).get(name, 0)

        if self._view_mode == VIEW_OCC:
            processed = snap.get("processed", 0)
            frac = nhits / processed if processed > 0 else float("nan")
            if processed == 0:
                return f"{name}\nNo data"
            return (f"{name}\n"
                    f"Module hits: {nhits:,}\n"
                    f"Occupancy: {frac:.4f}")

        rate = snap.get("rates", {}).get(scint, {}).get(name, math.nan)
        if math.isnan(rate):
            return f"{name}\nNo data"
        mode = snap.get("mode", MODE_AND)
        if mode == MODE_AND:
            denom = nhits
            denom_label = "Module hits"
        else:
            denom = snap.get("scint_hits", {}).get(scint, 0)
            denom_label = f"{scint} hits"
        ncoinc = round(rate * denom)
        return (f"{name}\n"
                f"Coincidence rate: {rate:.4f}\n"
                f"Coincidences: {ncoinc:,}\n"
                f"{denom_label}: {denom:,}")


# ===========================================================================
#  Main window
# ===========================================================================

class MainWindow(QMainWindow):

    def __init__(self, server_url: str):
        super().__init__()
        self.setWindowTitle("Scintillator–HyCal Coincidence Monitor")
        self.resize(1400, 900)

        self._server_url      = server_url
        self._n_events        = 0
        self._local_evio_paths: List[str] = []   # files queued for local analysis
        self._stats_worker:   Optional[ProcessWorker]        = None
        self._instant_worker: Optional[InstantDisplayWorker] = None
        self._fetcher:  Optional[WaveformFetcher] = None
        self._snapshot: dict = {}
        self._selected_module: str = ""   # module last clicked on the map

        crate_to_roc = _load_crate_to_roc(DAQ_CFG_JSON)
        daq = _load_daq_map(DAQ_MAP_JSON, crate_to_roc)
        self._modules = load_modules(MODULES_JSON)
        self._w_layers: Dict[str, int] = _build_w_module_layers(MODULES_JSON)
        self._neighbor_map = _build_neighbor_map(self._modules)
        physics_names = {m.name for m in self._modules if m.mod_type != "LMS"}
        self._mod_keys: Dict[str, str] = {
            name: key for name, key in daq.items()
            if name in physics_names and name not in SCINTILLATORS
        }
        # Scintillator channel keys resolved from daq_map (V1-V4 are listed there)
        self._scint_keys: Dict[str, str] = {
            sname: daq[sname] for sname in SCINTILLATORS if sname in daq
        }

        self._veto_ctrl    = VetoMotorController()
        self._et_mode      = False         # True when server is in online/ET mode
        self._map_view     = VIEW_COINC
        self._display_mode = DISPLAY_COINC # top-level mode switch
        self._display_paused = False       # True while instant display is frozen

        self._build_ui()
        self._map.set_modules(self._modules)
        self._map.set_palette("rainbow")
        self._map.set_range(0.0, 1.0)

        self._veto_timer = QTimer(self)
        self._veto_timer.timeout.connect(self._poll_veto_positions)
        self._veto_timer.start(VETO_POLL_MS)

    # ------------------------------------------------------------------
    #  UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        h_split = QSplitter(Qt.Orientation.Horizontal)

        # ---- left control panel ----------------------------------------
        left = QWidget()
        left.setStyleSheet(themed(f"QWidget{{background:{THEME.PANEL};}}"))
        lv = QVBoxLayout(left)
        lv.setContentsMargins(10, 10, 10, 10)
        lv.setSpacing(8)

        # Server
        srv_box = QGroupBox("Server")
        srv_box.setStyleSheet(self._groupbox_style())
        sv = QVBoxLayout(srv_box)
        sv.setSpacing(4)
        self._url_edit = QLineEdit(self._server_url)
        self._url_edit.setStyleSheet(self._input_style())
        sv.addWidget(self._url_edit)
        self._connect_btn = QPushButton("Connect")
        self._connect_btn.setStyleSheet(self._btn_style())
        self._connect_btn.clicked.connect(self._on_connect)
        sv.addWidget(self._connect_btn)
        self._conn_label = QLabel("Not connected")
        self._conn_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        sv.addWidget(self._conn_label)
        self._mode_btn = QPushButton("Switch to ET Mode")
        self._mode_btn.setEnabled(False)
        self._mode_btn.setStyleSheet(self._btn_style())
        self._mode_btn.clicked.connect(self._on_mode_switch)
        sv.addWidget(self._mode_btn)
        lv.addWidget(srv_box)

        # Local EVIO files (prad2py direct — no server needed)
        local_box = QGroupBox("Local EVIO Files (no server)")
        local_box.setStyleSheet(self._groupbox_style())
        loc = QVBoxLayout(local_box)
        loc.setSpacing(4)

        # Button row: Add / Remove selected / Clear all
        loc_btn_row = QHBoxLayout()
        loc_btn_row.setSpacing(4)
        self._browse_btn = QPushButton("Add files…")
        self._browse_btn.setStyleSheet(self._btn_style())
        if not _HAVE_PRAD2PY:
            self._browse_btn.setEnabled(False)
            self._browse_btn.setToolTip("prad2py not found — rebuild with -DBUILD_PYTHON=ON")
        self._browse_btn.clicked.connect(self._on_add_evio_files)
        loc_btn_row.addWidget(self._browse_btn)
        self._remove_file_btn = QPushButton("Remove")
        self._remove_file_btn.setStyleSheet(self._btn_style())
        self._remove_file_btn.setEnabled(False)
        self._remove_file_btn.clicked.connect(self._on_remove_evio_file)
        loc_btn_row.addWidget(self._remove_file_btn)
        self._clear_files_btn = QPushButton("Clear")
        self._clear_files_btn.setStyleSheet(self._btn_style())
        self._clear_files_btn.setEnabled(False)
        self._clear_files_btn.clicked.connect(self._on_clear_evio_files)
        loc_btn_row.addWidget(self._clear_files_btn)
        loc.addLayout(loc_btn_row)

        # List showing selected files (basename; full path in tooltip)
        self._file_list = QListWidget()
        self._file_list.setMaximumHeight(80)
        self._file_list.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self._file_list.setStyleSheet(
            themed(f"QListWidget{{background:{THEME.PANEL};"
                   f"color:{THEME.TEXT};border:1px solid {THEME.BORDER};"
                   f"font-size:10px;}}"
                   f"QListWidget::item:selected{{background:{THEME.ACCENT};}}"))
        self._file_list.itemSelectionChanged.connect(self._on_file_selection_changed)
        loc.addWidget(self._file_list)

        self._local_file_label = QLabel("No files selected")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._local_file_label.setWordWrap(True)
        loc.addWidget(self._local_file_label)
        lv.addWidget(local_box)

        # Display Mode
        disp_box = QGroupBox("Display Mode")
        disp_box.setStyleSheet(self._groupbox_style())
        dpv = QVBoxLayout(disp_box)
        dpv.setSpacing(2)
        self._disp_group = QButtonGroup(self)
        self._rb_disp_coinc   = QRadioButton("Coincidence Stats")
        self._rb_disp_instant = QRadioButton("Instant Event Display")
        for rb in (self._rb_disp_coinc, self._rb_disp_instant):
            rb.setStyleSheet(f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._disp_group.addButton(rb)
            dpv.addWidget(rb)
        self._rb_disp_coinc.setChecked(True)
        self._disp_group.buttonClicked.connect(self._on_display_mode_changed)
        lv.addWidget(disp_box)

        # Scintillator selector
        sci_box = QGroupBox("Scintillator")
        sci_box.setStyleSheet(self._groupbox_style())
        scv = QVBoxLayout(sci_box)
        scv.setSpacing(2)
        self._scint_group = QButtonGroup(self)
        for name in SCINTILLATORS:
            rb = QRadioButton(name)
            rb.setStyleSheet(f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._scint_group.addButton(rb)
            scv.addWidget(rb)
            if name == "V1":
                rb.setChecked(True)
        self._scint_group.buttonClicked.connect(self._on_scint_changed)
        lv.addWidget(sci_box)

        # Map view mode (Coincidence Stats only)
        self._mapview_box = QGroupBox("Map View")
        self._mapview_box.setStyleSheet(self._groupbox_style())
        mvv = QVBoxLayout(self._mapview_box)
        mvv.setSpacing(2)
        self._mapview_group = QButtonGroup(self)
        self._rb_coinc = QRadioButton("Coincidence Rate")
        self._rb_occ   = QRadioButton("Occupancy (hit count)")
        for rb in (self._rb_coinc, self._rb_occ):
            rb.setStyleSheet(
                f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._mapview_group.addButton(rb)
            mvv.addWidget(rb)
        self._rb_coinc.setChecked(True)
        self._mapview_group.buttonClicked.connect(self._on_mapview_changed)
        lv.addWidget(self._mapview_box)

        # Thresholds
        thr_box = QGroupBox("Thresholds (ADC peak height)")
        thr_box.setStyleSheet(self._groupbox_style())
        tv = QVBoxLayout(thr_box)
        tv.setSpacing(4)
        lbl_s = QLabel("Scintillator:")
        lbl_s.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_s)
        self._scint_thr_spin = QDoubleSpinBox()
        self._scint_thr_spin.setRange(0, 100000)
        self._scint_thr_spin.setValue(DEFAULT_SCINT_THR)
        self._scint_thr_spin.setDecimals(0)
        self._scint_thr_spin.setSingleStep(50)
        self._scint_thr_spin.setStyleSheet(self._input_style())
        tv.addWidget(self._scint_thr_spin)

        lbl_tcut = QLabel("Scint time cut (ns):")
        lbl_tcut.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_tcut)

        tcut_row = QHBoxLayout()
        tcut_row.setSpacing(4)
        self._scint_tcut_min = QDoubleSpinBox()
        self._scint_tcut_min.setRange(0, 10000)
        self._scint_tcut_min.setValue(160.0)
        self._scint_tcut_min.setDecimals(0)
        self._scint_tcut_min.setSingleStep(4)
        self._scint_tcut_min.setSuffix(" ns")
        self._scint_tcut_min.setEnabled(True)
        self._scint_tcut_min.setStyleSheet(self._input_style())
        tcut_row.addWidget(self._scint_tcut_min)
        lbl_to = QLabel("to")
        lbl_to.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tcut_row.addWidget(lbl_to)
        self._scint_tcut_max = QDoubleSpinBox()
        self._scint_tcut_max.setRange(0, 10000)
        self._scint_tcut_max.setValue(200.0)
        self._scint_tcut_max.setDecimals(0)
        self._scint_tcut_max.setSingleStep(4)
        self._scint_tcut_max.setSuffix(" ns")
        self._scint_tcut_max.setEnabled(True)
        self._scint_tcut_max.setStyleSheet(self._input_style())
        tcut_row.addWidget(self._scint_tcut_max)
        tv.addLayout(tcut_row)

        lbl_h = QLabel("HyCal module:")
        lbl_h.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_h)
        self._hycal_thr_spin = QDoubleSpinBox()
        self._hycal_thr_spin.setRange(0, 100000)
        self._hycal_thr_spin.setValue(DEFAULT_HYCAL_THR)
        self._hycal_thr_spin.setDecimals(0)
        self._hycal_thr_spin.setSingleStep(10)
        self._hycal_thr_spin.setStyleSheet(self._input_style())
        tv.addWidget(self._hycal_thr_spin)

        lbl_htcut = QLabel("HyCal time cut (ns):")
        lbl_htcut.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_htcut)
        htcut_row = QHBoxLayout()
        htcut_row.setSpacing(4)
        self._hycal_tcut_min = QDoubleSpinBox()
        self._hycal_tcut_min.setRange(0, 10000)
        self._hycal_tcut_min.setValue(DEFAULT_HYCAL_TMIN)
        self._hycal_tcut_min.setDecimals(0)
        self._hycal_tcut_min.setSingleStep(4)
        self._hycal_tcut_min.setSuffix(" ns")
        self._hycal_tcut_min.setStyleSheet(self._input_style())
        htcut_row.addWidget(self._hycal_tcut_min)
        lbl_hto = QLabel("to")
        lbl_hto.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        htcut_row.addWidget(lbl_hto)
        self._hycal_tcut_max = QDoubleSpinBox()
        self._hycal_tcut_max.setRange(0, 10000)
        self._hycal_tcut_max.setValue(DEFAULT_HYCAL_TMAX)
        self._hycal_tcut_max.setDecimals(0)
        self._hycal_tcut_max.setSingleStep(4)
        self._hycal_tcut_max.setSuffix(" ns")
        self._hycal_tcut_max.setStyleSheet(self._input_style())
        htcut_row.addWidget(self._hycal_tcut_max)
        tv.addLayout(htcut_row)

        lbl_maxlm = QLabel("Max HyCal local maxima (1 = single cluster):")
        lbl_maxlm.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_maxlm)
        self._max_lm_spin = QSpinBox()
        self._max_lm_spin.setRange(1, 20)
        self._max_lm_spin.setValue(DEFAULT_MAX_LOCAL_MAXIMA)
        self._max_lm_spin.setStyleSheet(self._input_style())
        tv.addWidget(self._max_lm_spin)

        n_cpu_max = max(1, os.cpu_count() or 1)
        lbl_cpu = QLabel(f"Parallel CPUs for local files (1–{n_cpu_max}):")
        lbl_cpu.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_cpu)
        self._cpu_spin = QSpinBox()
        self._cpu_spin.setRange(1, n_cpu_max)
        self._cpu_spin.setValue(min(4, n_cpu_max))
        self._cpu_spin.setStyleSheet(self._input_style())
        tv.addWidget(self._cpu_spin)

        lbl_excl = QLabel("Inner W exclusion layers (0 = none, 16 = all W):")
        lbl_excl.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_excl)
        self._excl_spin = QSpinBox()
        self._excl_spin.setRange(0, 16)
        self._excl_spin.setValue(0)
        self._excl_spin.setStyleSheet(self._input_style())
        tv.addWidget(self._excl_spin)
        lbl_minmods = QLabel("Min HyCal modules fired (cluster cut, 1 = off):")
        lbl_minmods.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        tv.addWidget(lbl_minmods)
        self._min_mods_spin = QSpinBox()
        self._min_mods_spin.setRange(1, 50)
        self._min_mods_spin.setValue(DEFAULT_MIN_CLUSTER_MODS)
        self._min_mods_spin.setStyleSheet(self._input_style())
        tv.addWidget(self._min_mods_spin)
        lv.addWidget(thr_box)

        # Event selection mode (Coincidence Stats only)
        self._selmode_box = QGroupBox("Event Selection Mode")
        self._selmode_box.setStyleSheet(self._groupbox_style())
        mv = QVBoxLayout(self._selmode_box)
        mv.setSpacing(4)
        self._mode_grp = QButtonGroup(self)
        self._rb_and = QRadioButton("AND — selected veto + HyCal fire")
        self._rb_or  = QRadioButton("OR  — all HyCal events (no veto gate)")
        for rb in (self._rb_and, self._rb_or):
            rb.setStyleSheet(
                f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._mode_grp.addButton(rb)
            mv.addWidget(rb)
        self._rb_and.setChecked(True)
        lv.addWidget(self._selmode_box)

        # Processing
        proc_box = QGroupBox("Processing")
        proc_box.setStyleSheet(self._groupbox_style())
        pv = QVBoxLayout(proc_box)
        pv.setSpacing(6)

        self._start_btn = QPushButton("Start")
        self._start_btn.setEnabled(False)
        self._start_btn.setStyleSheet(self._btn_style(accent=True))
        self._start_btn.clicked.connect(self._on_start_stop)
        pv.addWidget(self._start_btn)

        pause_row = QHBoxLayout()
        pause_row.setSpacing(4)
        self._pause_btn = QPushButton("Pause")
        self._pause_btn.setEnabled(False)
        self._pause_btn.setStyleSheet(self._btn_style())
        self._pause_btn.clicked.connect(self._on_pause)
        pause_row.addWidget(self._pause_btn)
        self._resume_btn = QPushButton("Resume")
        self._resume_btn.setEnabled(False)
        self._resume_btn.setStyleSheet(self._btn_style())
        self._resume_btn.clicked.connect(self._on_resume)
        pause_row.addWidget(self._resume_btn)
        pv.addLayout(pause_row)

        self._screenshot_btn = QPushButton("Save screenshot…")
        self._screenshot_btn.setStyleSheet(self._btn_style())
        self._screenshot_btn.setToolTip(
            "Save a PNG snapshot of the entire window.")
        self._screenshot_btn.clicked.connect(self._on_screenshot)
        pv.addWidget(self._screenshot_btn)

        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setValue(0)
        self._progress.setStyleSheet(
            f"QProgressBar{{background:{THEME.PANEL};border:1px solid "
            f"{THEME.BORDER};border-radius:4px;height:14px;}}"
            f"QProgressBar::chunk{{background:{THEME.ACCENT_STRONG};"
            f"border-radius:3px;}}")
        pv.addWidget(self._progress)

        # Per-worker progress rows (parallel local-EVIO mode).  Hidden until a
        # parallel job kicks off; populated with one row per worker process.
        self._worker_prog_container = QWidget()
        wpc_layout = QVBoxLayout(self._worker_prog_container)
        wpc_layout.setSpacing(2)
        wpc_layout.setContentsMargins(0, 0, 0, 0)
        self._worker_prog_layout: QVBoxLayout = wpc_layout
        self._worker_prog_rows: List[tuple] = []   # [(label, bar), ...]
        self._worker_prog_container.setVisible(False)
        pv.addWidget(self._worker_prog_container)

        rate_row = QHBoxLayout()
        rate_row.setSpacing(6)
        lbl_rate = QLabel("Max ET rate:")
        lbl_rate.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        rate_row.addWidget(lbl_rate)
        self._max_rate_spin = QSpinBox()
        self._max_rate_spin.setRange(0, 10000)
        self._max_rate_spin.setValue(0)
        self._max_rate_spin.setSuffix(" ev/s")
        self._max_rate_spin.setSpecialValueText("unlimited")
        self._max_rate_spin.setSingleStep(50)
        self._max_rate_spin.setStyleSheet(self._input_style())
        self._max_rate_spin.setToolTip(
            "Maximum events per second consumed from the ET ring buffer.\n"
            "0 = process as fast as possible (unlimited).\n"
            "Only effective in ET mode.")
        rate_row.addWidget(self._max_rate_spin, 1)
        pv.addLayout(rate_row)

        # Waveform saving
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet(f"color:{THEME.BORDER};")
        pv.addWidget(sep)

        self._save_wfm_chk = QCheckBox("Save coincidence waveforms")
        self._save_wfm_chk.setStyleSheet(
            f"QCheckBox{{color:{THEME.TEXT};font-size:12px;}}")
        pv.addWidget(self._save_wfm_chk)

        wfm_row = QHBoxLayout()
        wfm_row.setSpacing(6)
        lbl_wfm = QLabel("Max records:")
        lbl_wfm.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        wfm_row.addWidget(lbl_wfm)
        self._wfm_max_spin = QSpinBox()
        self._wfm_max_spin.setRange(1, 100000)
        self._wfm_max_spin.setValue(200)
        self._wfm_max_spin.setSingleStep(100)
        self._wfm_max_spin.setStyleSheet(self._input_style())
        self._wfm_max_spin.setToolTip(
            "Maximum number of coincidence waveform pairs to save.")
        wfm_row.addWidget(self._wfm_max_spin, 1)
        pv.addLayout(wfm_row)

        self._wfm_label = QLabel("")
        self._wfm_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._wfm_label.setWordWrap(True)
        pv.addWidget(self._wfm_label)

        self._status_label = QLabel("Ready")
        self._status_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        self._status_label.setWordWrap(True)
        pv.addWidget(self._status_label)
        lv.addWidget(proc_box)

        # Event browser
        ev_box = QGroupBox("Event Browser")
        ev_box.setStyleSheet(self._groupbox_style())
        ev = QVBoxLayout(ev_box)
        ev.setSpacing(4)

        ev_row = QHBoxLayout()
        ev_row.setSpacing(4)
        lbl_ev = QLabel("Event #")
        lbl_ev.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        ev_row.addWidget(lbl_ev)
        self._ev_spin = QSpinBox()
        self._ev_spin.setRange(1, 1)
        self._ev_spin.setValue(1)
        self._ev_spin.setStyleSheet(self._input_style())
        ev_row.addWidget(self._ev_spin, 1)
        ev.addLayout(ev_row)

        self._fetch_btn = QPushButton("Fetch waveforms")
        self._fetch_btn.setEnabled(False)
        self._fetch_btn.setStyleSheet(self._btn_style())
        self._fetch_btn.clicked.connect(self._on_fetch_waveforms)
        ev.addWidget(self._fetch_btn)

        self._sel_mod_label = QLabel("Module: (click map)")
        self._sel_mod_label.setStyleSheet(
            f"color:{THEME.TEXT_DIM};font-size:11px;")
        self._sel_mod_label.setWordWrap(True)
        ev.addWidget(self._sel_mod_label)
        lv.addWidget(ev_box)

        # Statistics
        stats_box = QGroupBox("Statistics")
        stats_box.setStyleSheet(self._groupbox_style())
        stv = QVBoxLayout(stats_box)
        self._stats_label = QLabel("—")
        self._stats_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        self._stats_label.setWordWrap(True)
        stv.addWidget(self._stats_label)
        lv.addWidget(stats_box)

        # Veto motor positions — two sub-rows per motor:
        #   row A: name  |  "RBV:"  |  readback value  |  moving indicator
        #   row B: ""    |  "Set:"  |  spinbox          |  Move button
        veto_box = QGroupBox("Veto Motor Positions")
        veto_box.setStyleSheet(self._groupbox_style())
        vg = QGridLayout(veto_box)
        vg.setHorizontalSpacing(6)
        vg.setVerticalSpacing(2)
        vg.setContentsMargins(6, 8, 6, 8)
        self._veto_widgets: Dict[str, dict] = {}
        for i, vname in enumerate(VETO_PV_BASES):
            ra = i * 3       # row A: readback
            rb = i * 3 + 1  # row B: setpoint
            rc = i * 3 + 2  # row C: thin separator

            # -- row A: readback -----------------------------------------
            lbl_name = QLabel(vname)
            lbl_name.setStyleSheet(
                f"color:{THEME.ACCENT};font-size:12px;font-weight:bold;")
            vg.addWidget(lbl_name, ra, 0)

            lbl_rbv_key = QLabel("RBV:")
            lbl_rbv_key.setStyleSheet(
                f"color:{THEME.TEXT_DIM};font-size:11px;")
            vg.addWidget(lbl_rbv_key, ra, 1)

            lbl_rbv = QLabel("—")
            lbl_rbv.setStyleSheet(
                f"color:{THEME.TEXT};font-size:12px;font-weight:bold;")
            lbl_rbv.setMinimumWidth(70)
            vg.addWidget(lbl_rbv, ra, 2)

            lbl_movn = QLabel("● moving")
            lbl_movn.setStyleSheet(
                f"color:{THEME.TEXT_DIM};font-size:10px;")
            vg.addWidget(lbl_movn, ra, 3)

            # -- row B: setpoint + move ----------------------------------
            lbl_set_key = QLabel("Set:")
            lbl_set_key.setStyleSheet(
                f"color:{THEME.TEXT_DIM};font-size:11px;")
            vg.addWidget(lbl_set_key, rb, 1)

            spin = QDoubleSpinBox()
            spin.setRange(-9999, 9999)
            spin.setDecimals(2)
            spin.setSingleStep(0.1)
            spin.setValue(0.0)
            spin.setStyleSheet(self._input_style())
            vg.addWidget(spin, rb, 2)

            btn = QPushButton("Move")
            btn.setStyleSheet(self._btn_style(accent=True))
            btn.clicked.connect(lambda _, v=vname: self._move_veto(v))
            vg.addWidget(btn, rb, 3)

            # -- row C: separator ----------------------------------------
            sep = QFrame()
            sep.setFrameShape(QFrame.Shape.HLine)
            sep.setStyleSheet(f"color:{THEME.BORDER};")
            sep.setFixedHeight(6)
            vg.addWidget(sep, rc, 0, 1, 4)

            self._veto_widgets[vname] = {
                "rbv": lbl_rbv, "movn": lbl_movn,
                "spin": spin, "btn": btn,
            }

        if not self._veto_ctrl.available:
            warn = QLabel("pyepics not available")
            warn.setStyleSheet(f"color:{THEME.WARN};font-size:10px;")
            vg.addWidget(warn, len(VETO_PV_BASES) * 3, 0, 1, 4)
        lv.addWidget(veto_box)

        lv.addStretch(1)

        # Wrap left panel in a scroll area so the window can resize freely
        # vertically even when the panel content is taller than the window.
        left_scroll = QScrollArea()
        left_scroll.setWidget(left)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        left_scroll.setMinimumWidth(340)
        left_scroll.setMaximumWidth(500)
        left_scroll.setStyleSheet(
            themed(f"QScrollArea{{background:{THEME.PANEL};"
                   f"border:none;}}"))

        # ---- right side: map on top, waveforms on bottom ---------------
        v_split = QSplitter(Qt.Orientation.Vertical)

        # Thin toolbar above the HyCal map for step-through controls.
        map_container = QWidget()
        map_vl = QVBoxLayout(map_container)
        map_vl.setContentsMargins(0, 0, 0, 0)
        map_vl.setSpacing(2)

        step_bar = QWidget()
        step_bar.setStyleSheet(
            themed(f"QWidget{{background:{THEME.PANEL};"
                   f"border-bottom:1px solid {THEME.BORDER};}}"))
        step_hl = QHBoxLayout(step_bar)
        step_hl.setContentsMargins(8, 4, 8, 4)
        step_hl.setSpacing(8)

        self._step_chk = QCheckBox("Step through coincidence events")
        self._step_chk.setStyleSheet(
            f"QCheckBox{{color:{THEME.TEXT};font-size:12px;}}"
            f"QCheckBox::indicator{{width:14px;height:14px;"
            f"background:#ffffff;border:2px solid #888888;border-radius:3px;}}"
            f"QCheckBox::indicator:checked{{background:#4a9eff;"
            f"border:2px solid #4a9eff;}}"
        )
        self._step_chk.setToolTip(
            "Pause on each coincidence event found in the local EVIO file\n"
            "and show its waveforms in the panels below.\n"
            "Click 'Continue →' to advance to the next coincidence.\n"
            "Only available with local EVIO files.")
        step_hl.addWidget(self._step_chk)

        self._step_continue_btn = QPushButton("Continue →")
        self._step_continue_btn.setEnabled(False)
        self._step_continue_btn.setStyleSheet(self._btn_style(accent=True))
        self._step_continue_btn.setToolTip("Advance to the next coincidence event.")
        self._step_continue_btn.clicked.connect(self._on_step_continue)
        step_hl.addWidget(self._step_continue_btn)
        step_hl.addStretch(1)

        map_vl.addWidget(step_bar)

        self._map = CoincidenceMapWidget()
        self._map.moduleClicked.connect(self._on_module_clicked)
        map_vl.addWidget(self._map)

        v_split.addWidget(map_container)

        wave_row = QWidget()
        wave_row.setMinimumHeight(280)
        wrl = QHBoxLayout(wave_row)
        wrl.setContentsMargins(0, 0, 0, 0)
        wrl.setSpacing(4)
        self._wave_scint  = MultiWavePanel("Scintillators")
        self._wave_module = WavePanel("HyCal Module")
        wrl.addWidget(self._wave_scint)
        wrl.addWidget(self._wave_module)
        v_split.addWidget(wave_row)

        v_split.setStretchFactor(0, 1)
        v_split.setStretchFactor(1, 1)

        h_split.addWidget(left_scroll)
        h_split.addWidget(v_split)
        h_split.setStretchFactor(0, 2)
        h_split.setStretchFactor(1, 3)

        self.setCentralWidget(h_split)

    # ------------------------------------------------------------------
    #  Style helpers
    # ------------------------------------------------------------------

    def _groupbox_style(self) -> str:
        return themed(
            f"QGroupBox{{color:{THEME.TEXT_DIM};font-size:11px;font-weight:bold;"
            f"border:1px solid {THEME.BORDER};border-radius:6px;margin-top:6px;"
            f"padding-top:4px;}}"
            f"QGroupBox::title{{subcontrol-origin:margin;left:8px;padding:0 4px;}}")

    def _input_style(self) -> str:
        return themed(
            f"QLineEdit,QDoubleSpinBox,QSpinBox{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT};border:1px solid {THEME.BORDER};"
            f"border-radius:4px;padding:3px 6px;font-size:12px;}}"
            f"QLineEdit:focus,QDoubleSpinBox:focus,QSpinBox:focus{{"
            f"border-color:{THEME.ACCENT_BORDER};}}")

    def _btn_style(self, accent: bool = False) -> str:
        bg  = THEME.ACCENT_STRONG if accent else THEME.BUTTON
        hov = THEME.ACCENT        if accent else THEME.BUTTON_HOVER
        fg  = "#ffffff"           if accent else THEME.TEXT
        return themed(
            f"QPushButton{{background:{bg};color:{fg};border:1px solid "
            f"{THEME.BORDER};border-radius:6px;padding:5px 10px;font-size:12px;}}"
            f"QPushButton:hover{{background:{hov};}}"
            f"QPushButton:disabled{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT_MUTED};}}")

    # ------------------------------------------------------------------
    #  Slots — connection & scan
    # ------------------------------------------------------------------

    def _on_connect(self):
        url = self._url_edit.text().rstrip("/")
        self._server_url = url
        self._conn_label.setText("Connecting…")
        self._connect_btn.setEnabled(False)
        QApplication.processEvents()

        cfg = _http_get(f"{url}/api/config")
        if cfg is None:
            self._conn_label.setText("Connection failed")
            self._conn_label.setStyleSheet(f"color:{THEME.DANGER};font-size:11px;")
            self._connect_btn.setEnabled(True)
            return

        # Server connection takes priority — clear any local file selection
        # so Start uses the server, not the previously queued EVIO files.
        self._local_evio_paths.clear()
        self._file_list.clear()
        self._local_file_label.setText("No files selected")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._remove_file_btn.setEnabled(False)
        self._clear_files_btn.setEnabled(False)

        self._apply_config(cfg)
        self._connect_btn.setEnabled(True)

    def _on_add_evio_files(self):
        """Open a multi-select file dialog and add EVIO files to the queue."""
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add EVIO Files", "",
            "EVIO files (*.evio *.evio.*);;All files (*)")
        if not paths:
            return

        dec     = _prad2py.dec
        cfg_dec = dec.load_daq_config()
        ch      = dec.EvChannel()
        ch.set_config(cfg_dec)

        added = 0
        for path in paths:
            if path in self._local_evio_paths:
                continue   # skip duplicates
            if ch.open_auto(path) != dec.Status.success:
                ch.close()
                continue   # skip unreadable files silently
            n_records = (ch.get_random_access_event_count()
                         if ch.is_random_access() else 0)
            ch.close()

            self._local_evio_paths.append(path)
            fname = Path(path).name
            n_str = f"{n_records:,} rec" if n_records else "seq"
            item_text = f"{fname}  ({n_str})"
            from PyQt6.QtWidgets import QListWidgetItem
            item = QListWidgetItem(item_text)
            item.setToolTip(path)
            self._file_list.addItem(item)
            added += 1

        if not self._local_evio_paths:
            return

        n = len(self._local_evio_paths)
        self._local_file_label.setText(f"{n} file(s) queued")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.SUCCESS};font-size:10px;")
        self._clear_files_btn.setEnabled(True)

        # Local file mode takes priority — clear any server state.
        self._n_events = 0
        self._et_mode  = False
        self._conn_label.setText(f"Using local files (no server)")
        self._conn_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        self._start_btn.setEnabled(True)
        self._fetch_btn.setEnabled(False)
        self._mode_btn.setEnabled(False)

    def _on_remove_evio_file(self):
        """Remove the currently selected file from the queue."""
        row = self._file_list.currentRow()
        if row < 0 or row >= len(self._local_evio_paths):
            return
        self._file_list.takeItem(row)
        self._local_evio_paths.pop(row)
        self._on_file_selection_changed()
        n = len(self._local_evio_paths)
        if n == 0:
            self._local_file_label.setText("No files selected")
            self._local_file_label.setStyleSheet(
                f"color:{THEME.TEXT_DIM};font-size:10px;")
            self._clear_files_btn.setEnabled(False)
            self._start_btn.setEnabled(
                self._n_events > 0 or self._et_mode)
        else:
            self._local_file_label.setText(f"{n} file(s) queued")

    def _on_clear_evio_files(self):
        """Remove all files from the queue."""
        self._local_evio_paths.clear()
        self._file_list.clear()
        self._local_file_label.setText("No files selected")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._remove_file_btn.setEnabled(False)
        self._clear_files_btn.setEnabled(False)
        self._start_btn.setEnabled(self._n_events > 0 or self._et_mode)
        self._conn_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")

    def _on_file_selection_changed(self):
        self._remove_file_btn.setEnabled(
            self._file_list.currentRow() >= 0)

    def _apply_config(self, cfg: dict) -> None:
        """Update UI state from a /api/config response dict."""
        n            = cfg.get("event_count", 0)
        mode         = cfg.get("mode", "unknown")
        et_connected = cfg.get("et_connected", False)
        self._n_events = n
        self._et_mode  = (mode == "online")

        if n > 0:
            self._conn_label.setText(f"Connected  ({n:,} events, mode={mode})")
            self._conn_label.setStyleSheet(f"color:{THEME.SUCCESS};font-size:11px;")
            self._start_btn.setEnabled(True)
            self._fetch_btn.setEnabled(True)
            self._ev_spin.setRange(1, n)
        elif self._et_mode and et_connected:
            self._conn_label.setText("Connected — ET online (live accumulation)")
            self._conn_label.setStyleSheet(f"color:{THEME.SUCCESS};font-size:11px;")
            self._start_btn.setEnabled(True)
            self._fetch_btn.setEnabled(True)   # fetches latest ring event
        elif self._et_mode:
            self._conn_label.setText("Connected — ET online (waiting for DAQ)")
            self._conn_label.setStyleSheet(f"color:{THEME.WARN};font-size:11px;")
            self._start_btn.setEnabled(True)
            self._fetch_btn.setEnabled(False)
        else:
            self._conn_label.setText(f"Connected — no file loaded (mode={mode})")
            self._conn_label.setStyleSheet(f"color:{THEME.WARN};font-size:11px;")
            self._start_btn.setEnabled(False)
            self._fetch_btn.setEnabled(False)

        self._mode_btn.setEnabled(True)
        if self._et_mode:
            self._mode_btn.setText("Switch to File Mode")
        else:
            self._mode_btn.setText("Switch to ET Mode")

    def _any_worker_running(self) -> bool:
        return (
            (self._stats_worker   is not None and self._stats_worker.isRunning()) or
            (self._instant_worker is not None and self._instant_worker.isRunning())
        )

    def _on_mode_switch(self):
        """Toggle the server between online (ET) and file mode."""
        if self._any_worker_running():
            self._stop_worker()
            if self._stats_worker:   self._stats_worker.wait(2000)
            if self._instant_worker: self._instant_worker.wait(2000)

        endpoint = "/api/mode/file" if self._et_mode else "/api/mode/online"
        self._mode_btn.setEnabled(False)
        self._mode_btn.setText("Switching…")
        QApplication.processEvents()

        result = _http_get(f"{self._server_url}{endpoint}")
        if result is None:
            self._conn_label.setText("Mode switch failed")
            self._conn_label.setStyleSheet(f"color:{THEME.DANGER};font-size:11px;")
            self._mode_btn.setEnabled(True)
            self._mode_btn.setText(
                "Switch to File Mode" if self._et_mode else "Switch to ET Mode")
            return

        cfg = _http_get(f"{self._server_url}/api/config")
        if cfg:
            self._apply_config(cfg)

    def _active_scint(self) -> str:
        btn = self._scint_group.checkedButton()
        return btn.text() if btn else "V1"


    def _on_display_mode_changed(self):
        instant = self._rb_disp_instant.isChecked()
        self._display_mode = DISPLAY_INSTANT if instant else DISPLAY_COINC

        # Pause/Resume only active in Instant mode while the instant worker runs.
        worker_live = self._instant_worker is not None and self._instant_worker.isRunning()
        if instant and worker_live:
            if not self._display_paused:
                self._pause_btn.setEnabled(True)
                self._resume_btn.setEnabled(False)
            else:
                self._pause_btn.setEnabled(False)
                self._resume_btn.setEnabled(True)
        else:
            self._display_paused = False
            self._pause_btn.setEnabled(False)
            self._resume_btn.setEnabled(False)

        # Map view (Coinc Rate / Occupancy) is only meaningful in stats mode.
        # AND/OR mode stays visible because it governs the background stats worker
        # that accumulates in both display modes.
        self._mapview_box.setVisible(not instant)

        # Switch the map view without stopping workers or resetting stats.
        if instant:
            self._wave_module.clear("Max ADC Module")
            # If no live instant data yet, blank the map until the next event.
            if not (self._instant_worker and self._instant_worker.isRunning()):
                self._map.set_values({})
                self._map.update()
        else:
            self._wave_module.clear("HyCal Module")
            # Restore the accumulated stats immediately.
            if self._snapshot:
                self._map.set_snapshot(self._snapshot, self._active_scint(),
                                       self._map_view)
            else:
                self._map.set_values({})
                self._map.update()

    def _on_scint_changed(self):
        if self._snapshot and self._display_mode == DISPLAY_COINC:
            self._map.set_snapshot(self._snapshot, self._active_scint(), self._map_view)
            self._update_stats_label(self._snapshot)
        # refresh scintillator waveform panel title & re-fetch if event is available
        if self._n_events > 0 or self._et_mode:
            self._wave_scint.clear("Scintillators")
            self._on_fetch_waveforms()

    def _on_mapview_changed(self):
        self._map_view = VIEW_OCC if self._rb_occ.isChecked() else VIEW_COINC
        if self._snapshot:
            self._map.set_snapshot(self._snapshot, self._active_scint(), self._map_view)
            self._update_stats_label(self._snapshot)

    def _on_start_stop(self):
        if self._any_worker_running():
            self._stop_worker()
        else:
            self._start_worker()

    def _on_pause(self):
        self._display_paused = True
        self._pause_btn.setEnabled(False)
        self._resume_btn.setEnabled(True)

    def _on_resume(self):
        self._display_paused = False
        self._pause_btn.setEnabled(True)
        self._resume_btn.setEnabled(False)

    def _on_screenshot(self):
        import datetime
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        default_path = str(Path.cwd() / f"coinc_monitor_{ts}.png")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save screenshot", default_path,
            "PNG image (*.png);;All files (*)")
        if not path:
            return
        if not path.lower().endswith(".png"):
            path += ".png"
        pixmap = self.grab()
        if pixmap.save(path, "PNG"):
            self._status_label.setText(f"Screenshot saved → {path}")
        else:
            QMessageBox.warning(self, "Screenshot failed",
                                f"Could not write PNG to:\n{path}")

    def _start_worker(self):
        self._snapshot = {}
        self._display_paused = False
        self._map.set_values({})
        self._map.update()

        self._connect_btn.setEnabled(False)
        self._mode_btn.setEnabled(False)
        self._browse_btn.setEnabled(False)
        self._start_btn.setText("Stop")
        self._start_btn.setStyleSheet(self._btn_style(accent=False))

        # Apply inner-layer exclusion: drop W modules whose layer ≤ excl_layers.
        excl = self._excl_spin.value()
        if excl > 0:
            active_mod_keys = {
                name: key for name, key in self._mod_keys.items()
                if not name.startswith('W') or self._w_layers.get(name, 0) > excl
            }
        else:
            active_mod_keys = self._mod_keys

        sel_mode = MODE_AND if self._rb_and.isChecked() else MODE_OR

        min_mods = self._min_mods_spin.value()

        scint_t_min = self._scint_tcut_min.value()
        scint_t_max = self._scint_tcut_max.value()

        hycal_t_min      = self._hycal_tcut_min.value()
        hycal_t_max      = self._hycal_tcut_max.value()
        max_local_maxima = self._max_lm_spin.value()
        n_workers        = self._cpu_spin.value()

        # --- Waveform collector (optional) ---
        wfm_coll = wfm_path = None
        if self._save_wfm_chk.isChecked():
            import datetime
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            wfm_path = Path(f"coinc_wfm_{ts}.npz")
            wfm_coll = WaveformCollector(
                server_url=self._server_url,
                max_records=self._wfm_max_spin.value(),
            )
            self._wfm_label.setText(f"Will save to: {wfm_path}")

        # --- Stats worker (always) ---
        if self._local_evio_paths:
            self._stats_worker = ProcessWorkerLocal(
                evio_paths=self._local_evio_paths,
                module_keys=active_mod_keys,
                scint_keys=self._scint_keys,
                scint_thr=self._scint_thr_spin.value(),
                hycal_thr=self._hycal_thr_spin.value(),
                mode=sel_mode,
                min_mods=min_mods,
                scint_t_min=scint_t_min,
                scint_t_max=scint_t_max,
                hycal_t_min=hycal_t_min,
                hycal_t_max=hycal_t_max,
                neighbor_map=self._neighbor_map,
                max_local_maxima=max_local_maxima,
                wfm_collector=wfm_coll,
                wfm_save_path=wfm_path,
                step_through=self._step_chk.isChecked(),
                n_workers=n_workers,
                parent=self,
            )
            self._progress.setRange(0, 100)
        elif self._et_mode:
            self._stats_worker = ProcessWorkerET(
                server_url=self._server_url,
                module_keys=active_mod_keys,
                scint_keys=self._scint_keys,
                scint_thr=self._scint_thr_spin.value(),
                hycal_thr=self._hycal_thr_spin.value(),
                mode=sel_mode,
                min_mods=min_mods,
                max_rate_hz=float(self._max_rate_spin.value()),
                scint_t_min=scint_t_min,
                scint_t_max=scint_t_max,
                hycal_t_min=hycal_t_min,
                hycal_t_max=hycal_t_max,
                neighbor_map=self._neighbor_map,
                max_local_maxima=max_local_maxima,
                wfm_collector=wfm_coll,
                wfm_save_path=wfm_path,
                parent=self,
            )
            self._progress.setRange(0, 0)
        else:
            self._stats_worker = ProcessWorker(
                server_url=self._server_url,
                n_events=self._n_events,
                module_keys=active_mod_keys,
                scint_keys=self._scint_keys,
                scint_thr=self._scint_thr_spin.value(),
                hycal_thr=self._hycal_thr_spin.value(),
                mode=sel_mode,
                min_mods=min_mods,
                scint_t_min=scint_t_min,
                scint_t_max=scint_t_max,
                hycal_t_min=hycal_t_min,
                hycal_t_max=hycal_t_max,
                neighbor_map=self._neighbor_map,
                max_local_maxima=max_local_maxima,
                wfm_collector=wfm_coll,
                wfm_save_path=wfm_path,
                parent=self,
            )
            self._progress.setRange(0, 100)

        self._stats_worker.progress.connect(self._on_progress)
        self._stats_worker.stats_update.connect(self._on_stats_update)
        self._stats_worker.finished.connect(self._on_finished)
        self._stats_worker.waveforms_saved.connect(self._on_waveforms_saved)
        if isinstance(self._stats_worker, ProcessWorkerLocal):
            self._stats_worker.coincidence_event.connect(self._on_coincidence_event)
            self._stats_worker.workers_setup.connect(self._setup_worker_prog_rows)
            self._stats_worker.worker_progress.connect(self._on_worker_progress)
        self._stats_worker.start()

        # Instant display worker: ET mode only (not available for local files).
        if self._et_mode and not self._local_evio_paths:
            self._instant_worker = InstantDisplayWorker(
                server_url=self._server_url,
                module_keys=active_mod_keys,
                hycal_thr=self._hycal_thr_spin.value(),
                min_mods=min_mods,
                parent=self,
            )
            self._instant_worker.event_ready.connect(self._on_instant_event)
            self._instant_worker.finished.connect(self._on_finished)
            self._instant_worker.start()
        else:
            self._instant_worker = None

        # Pause/Resume only meaningful in Instant Event Display mode (ET).
        if self._et_mode and not self._local_evio_paths and self._display_mode == DISPLAY_INSTANT:
            self._pause_btn.setEnabled(True)
            self._resume_btn.setEnabled(False)
        else:
            self._pause_btn.setEnabled(False)
            self._resume_btn.setEnabled(False)

        n_local = len(self._local_evio_paths)
        src = (f"{n_local} local file(s)" if self._local_evio_paths
               else ("ET live" if self._et_mode else "file"))
        self._status_label.setText(
            f"{'Accumulating' if self._et_mode else 'Processing'}… "
            f"({src}, mode: {sel_mode})")
        self._scint_thr_spin.setEnabled(False)
        self._scint_tcut_min.setEnabled(False)
        self._scint_tcut_max.setEnabled(False)
        self._hycal_thr_spin.setEnabled(False)
        self._hycal_tcut_min.setEnabled(False)
        self._hycal_tcut_max.setEnabled(False)
        self._max_lm_spin.setEnabled(False)
        self._cpu_spin.setEnabled(False)
        self._excl_spin.setEnabled(False)
        self._min_mods_spin.setEnabled(False)
        self._max_rate_spin.setEnabled(False)
        self._rb_and.setEnabled(False)
        self._rb_or.setEnabled(False)
        self._save_wfm_chk.setEnabled(False)
        self._wfm_max_spin.setEnabled(False)
        self._step_chk.setEnabled(False)

    def _stop_worker(self):
        if self._stats_worker:   self._stats_worker.stop()
        if self._instant_worker: self._instant_worker.stop()
        self._status_label.setText("Stopping…")
        self._start_btn.setEnabled(False)

    # ------------------------------------------------------------------
    #  Per-worker progress (parallel local-EVIO mode)
    # ------------------------------------------------------------------

    def _setup_worker_prog_rows(self, n_workers: int):
        self._clear_worker_prog_rows()
        if n_workers <= 0:
            return
        for i in range(n_workers):
            row = QWidget()
            rh  = QHBoxLayout(row)
            rh.setSpacing(4)
            rh.setContentsMargins(0, 0, 0, 0)

            lbl = QLabel(f"W{i}: waiting…")
            lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:10px;")
            lbl.setMinimumWidth(180)
            lbl.setMaximumWidth(220)
            lbl.setWordWrap(False)
            lbl.setToolTip("")
            rh.addWidget(lbl, 0)

            bar = QProgressBar()
            bar.setRange(0, 0)   # indeterminate until first update
            bar.setValue(0)
            bar.setTextVisible(True)
            bar.setFormat("%p%")
            bar.setFixedHeight(12)
            bar.setStyleSheet(
                f"QProgressBar{{background:{THEME.PANEL};border:1px solid "
                f"{THEME.BORDER};border-radius:3px;font-size:9px;}}"
                f"QProgressBar::chunk{{background:{THEME.ACCENT};"
                f"border-radius:2px;}}")
            rh.addWidget(bar, 1)

            self._worker_prog_layout.addWidget(row)
            self._worker_prog_rows.append((lbl, bar))
        self._worker_prog_container.setVisible(True)

    def _clear_worker_prog_rows(self):
        while self._worker_prog_layout.count():
            item = self._worker_prog_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self._worker_prog_rows = []
        self._worker_prog_container.setVisible(False)

    def _on_worker_progress(self, info: dict):
        wid = int(info.get("worker_id", -1))
        if wid < 0 or wid >= len(self._worker_prog_rows):
            return
        lbl, bar = self._worker_prog_rows[wid]
        n_files   = int(info.get("n_files", 0))
        fidx      = int(info.get("file_idx", -1))
        basename  = str(info.get("file_basename", ""))
        rec_done  = int(info.get("records_done", 0))
        rec_tot   = int(info.get("records_total", 0))
        finished  = bool(info.get("finished", False))

        if finished:
            lbl.setText(f"W{wid}: done")
            bar.setRange(0, 100)
            bar.setValue(100)
            return

        if fidx >= 0 and basename:
            lbl.setText(f"W{wid} [{fidx + 1}/{n_files}] {basename}")
            lbl.setToolTip(basename)
        else:
            lbl.setText(f"W{wid}: starting…")

        if rec_tot > 0:
            bar.setRange(0, rec_tot)
            bar.setValue(min(rec_done, rec_tot))
        else:
            bar.setRange(0, 0)   # indeterminate
            bar.setValue(0)

    def _on_progress(self, processed: int, total: int):
        if self._display_mode == DISPLAY_INSTANT:
            return   # status label owned by _on_instant_event in this mode
        if total < 0:   # ET mode — no total
            self._status_label.setText(f"Accumulated {processed:,} events")
        elif total == 0:   # local EVIO sequential (unknown total)
            self._status_label.setText(f"Processed {processed:,} events…")
        else:
            pct = min(100, int(100 * processed / total))
            self._progress.setValue(pct)
            self._status_label.setText(f"Processed {processed:,} / {total:,}")

    def _on_stats_update(self, snap: dict):
        self._snapshot = snap
        # Only repaint the map when the user is viewing coinc/occ stats.
        if self._display_mode == DISPLAY_COINC:
            self._map.set_snapshot(snap, self._active_scint(), self._map_view)
        self._update_stats_label(snap)

    def _update_stats_label(self, snap: dict):
        scint = self._active_scint()
        mhits = snap.get("module_hits", {})
        processed = snap.get("processed", 0)
        if self._map_view == VIEW_OCC:
            total_hits = sum(mhits.values())
            max_hits   = max(mhits.values()) if mhits else 0
            nonzero    = sum(1 for v in mhits.values() if v > 0)
            mean_hits  = total_hits / nonzero if nonzero else 0.0
            self._stats_label.setText(
                f"Events processed: {processed:,}\n"
                f"Module hits (total): {total_hits:,}\n"
                f"Non-zero modules: {nonzero}\n"
                f"Mean hits: {mean_hits:.1f}\n"
                f"Max hits:  {max_hits:,}"
            )
        else:
            shits     = snap.get("scint_hits",     {}).get(scint, 0)
            scint_hits_any = snap.get("scint_hits_any", {})
            rates = snap.get("rates", {}).get(scint, {})
            valid = [v for v in rates.values() if not math.isnan(v) and v > 0]
            mean_rate = sum(valid) / len(valid) if valid else 0.0
            max_rate  = max(valid)              if valid else 0.0
            # Show all 4 scintillators' individual (pre-AND) hit counts at once
            # so the user can immediately spot which scintillator is not firing.
            indiv_parts = "  ".join(
                f"{s}={scint_hits_any.get(s, 0):,}" for s in SCINTILLATORS
            )
            and_label = (f"{scint} AND HyCal hits: {shits:,}"
                         if self._rb_and.isChecked()
                         else f"{scint} hits: {shits:,}")
            self._stats_label.setText(
                f"Indiv. hits: {indiv_parts}\n"
                f"{and_label}\n"
                f"Module hits (total): {sum(mhits.values()):,}\n"
                f"Non-zero modules: {len(valid)}\n"
                f"Mean rate: {mean_rate:.4f}\n"
                f"Max rate:  {max_rate:.4f}"
            )

    def _on_coincidence_event(self, data: dict):
        """Step-through mode: show waveforms for a single coincidence event and wait."""
        sf         = data["scint_fired"]
        scint_name = self._active_scint()

        # Only pause when the currently-selected scintillator itself fired in-window.
        if not sf.get(scint_name, False):
            if isinstance(self._stats_worker, ProcessWorkerLocal):
                self._stats_worker.step_continue()
            return

        ev_num  = data["event_number"]
        best    = data["best"]
        best_adc = data.get("best_adc", 0.0)
        wfm_ch  = data["wfm_channels"]
        t_min = self._scint_tcut_min.value()
        t_max = self._scint_tcut_max.value()
        scint_wfms = {
            sname: wfm_ch.get(skey, {"error": "no waveform"})
            for sname, skey in self._scint_keys.items()
        }
        self._wave_scint.set_multi_data(
            scint_wfms, self._scint_thr_spin.value(),
            list(SCINTILLATORS),
            f"Scintillators  (event {ev_num})",
            t_min=t_min, t_max=t_max)

        best_key   = self._mod_keys.get(best, "")
        module_wfm = wfm_ch.get(best_key, {"error": "no waveform"})
        self._wave_module.set_data(
            module_wfm, self._hycal_thr_spin.value(),
            f"{best}  ADC={best_adc:.0f}  (event {ev_num})")

        fired_names = [s for s, f in sf.items() if f]
        self._status_label.setText(
            f"[step] Event {ev_num}  |  best: {best}  |  fired: "
            f"{', '.join(fired_names) if fired_names else 'none'}\n"
            f"Click 'Continue →' for next coincidence.")
        self._step_continue_btn.setEnabled(True)

    def _on_step_continue(self):
        self._step_continue_btn.setEnabled(False)
        if isinstance(self._stats_worker, ProcessWorkerLocal):
            self._stats_worker.step_continue()

    def _on_instant_event(self, data: dict):
        """Handle a new live event from InstantDisplayWorker."""
        if self._display_mode != DISPLAY_INSTANT:
            return   # accumulating in background; don't touch the map or panels
        if self._display_paused:
            return   # paused: freeze display while stats keep accumulating

        adc_vals = data.get("adc_vals", {})
        channels = data.get("channels", {})
        seq      = data.get("seq", 0)
        max_mod  = data.get("max_mod", "")

        self._map.set_instant_event(adc_vals, max_mod)

        scint_wfms = {}
        for sname, skey in self._scint_keys.items():
            ch = channels.get(skey, {}) if skey else {}
            scint_wfms[sname] = ch if ch else {"error": "no signal"}
        self._wave_scint.set_multi_data(
            scint_wfms, self._scint_thr_spin.value(),
            list(SCINTILLATORS),
            f"Scintillators  (seq {seq})")

        if max_mod:
            max_key = self._mod_keys.get(max_mod, "")
            max_ch  = channels.get(max_key, {}) if max_key else {}
            max_adc = adc_vals.get(max_mod, 0.0)
            self._wave_module.set_data(
                max_ch if max_ch else {"error": "no waveform"},
                self._hycal_thr_spin.value(),
                f"{max_mod}  ADC={max_adc:.0f}  (seq {seq})")

        self._status_label.setText(
            f"Seq {seq}  |  max: {max_mod}  {adc_vals.get(max_mod, 0):.0f}")

    def _on_finished(self, msg: str):
        # Both workers emit finished; wait until neither is running.
        if self._any_worker_running():
            return
        self._display_paused = False
        self._pause_btn.setEnabled(False)
        self._resume_btn.setEnabled(False)
        self._progress.setRange(0, 100)
        self._progress.setValue(100 if msg == "" else self._progress.value())
        self._clear_worker_prog_rows()
        self._start_btn.setText("Start")
        self._start_btn.setStyleSheet(self._btn_style(accent=True))
        self._start_btn.setEnabled(True)
        self._connect_btn.setEnabled(True)
        self._mode_btn.setEnabled(True)
        if _HAVE_PRAD2PY:
            self._browse_btn.setEnabled(True)
            self._remove_file_btn.setEnabled(self._file_list.currentRow() >= 0)
            self._clear_files_btn.setEnabled(len(self._local_evio_paths) > 0)
        self._scint_thr_spin.setEnabled(True)
        self._scint_tcut_min.setEnabled(True)
        self._scint_tcut_max.setEnabled(True)
        self._hycal_thr_spin.setEnabled(True)
        self._hycal_tcut_min.setEnabled(True)
        self._hycal_tcut_max.setEnabled(True)
        self._max_lm_spin.setEnabled(True)
        self._cpu_spin.setEnabled(True)
        self._excl_spin.setEnabled(True)
        self._min_mods_spin.setEnabled(True)
        self._max_rate_spin.setEnabled(True)
        self._rb_and.setEnabled(True)
        self._rb_or.setEnabled(True)
        self._save_wfm_chk.setEnabled(True)
        self._wfm_max_spin.setEnabled(True)
        self._step_chk.setEnabled(True)
        self._step_continue_btn.setEnabled(False)
        if msg == "stopped":
            self._status_label.setText("Stopped.")
        else:
            n = self._snapshot.get("processed", 0)
            self._status_label.setText(f"Done — {n:,} events processed.")

    def _on_waveforms_saved(self, count: int, path: str):
        self._wfm_label.setText(f"Saved {count} waveform pairs → {path}")

    # ------------------------------------------------------------------
    #  Slots — waveform browser
    # ------------------------------------------------------------------

    def _on_module_clicked(self, name: str):
        if not name:
            return
        if name in SCINTILLATORS:
            # Scintillator clicked on map → select its radio button and refresh left panel
            for btn in self._scint_group.buttons():
                if btn.text() == name:
                    btn.setChecked(True)
                    break
            self._on_scint_changed()
        else:
            # HyCal module clicked → update right waveform panel
            self._selected_module = name
            self._sel_mod_label.setText(f"Module: {name}")
            self._on_fetch_waveforms()

    def _on_fetch_waveforms(self):
        if self._n_events == 0 and not self._et_mode:
            return

        # Stop any in-flight fetcher
        if self._fetcher and self._fetcher.isRunning():
            self._fetcher.wait(500)

        ev         = self._ev_spin.value()
        module_key = self._mod_keys.get(self._selected_module, "")

        ev_label = "latest ET event" if self._et_mode else f"event {ev}"
        self._wave_scint.clear("Scintillators", f"Fetching {ev_label}…")
        self._wave_module.clear(
            self._selected_module or "HyCal Module",
            f"Fetching {ev_label}…" if self._selected_module
            else "Click a module on the map")

        self._fetcher = WaveformFetcher(
            self._server_url, ev, self._scint_keys, module_key,
            et_mode=self._et_mode, parent=self)

        self._fetcher.waveform_ready.connect(self._on_waveform_ready)
        self._fetcher.start()

    def _on_waveform_ready(self, label: str, data: dict):
        ev = self._ev_spin.value()
        if label == "scint":
            t_min = self._scint_tcut_min.value()
            t_max = self._scint_tcut_max.value()
            title = f"Scintillators  (event {ev})"
            self._wave_scint.set_multi_data(
                data, self._scint_thr_spin.value(),
                list(SCINTILLATORS), title,
                t_min=t_min, t_max=t_max)
        else:
            mod_name = self._selected_module or "HyCal Module"
            title = f"{mod_name}  (event {ev})"
            self._wave_module.set_data(
                data, self._hycal_thr_spin.value(), title)

    # ------------------------------------------------------------------
    #  Veto position polling
    # ------------------------------------------------------------------

    def _poll_veto_positions(self):
        for vname in VETO_PV_BASES:
            rbv  = self._veto_ctrl.get(f"{vname}_rbv")
            movn = self._veto_ctrl.get(f"{vname}_movn")
            w = self._veto_widgets[vname]

            if rbv is not None:
                w["rbv"].setText(f"{rbv:.2f}")
                w["rbv"].setStyleSheet(f"color:{THEME.TEXT};font-size:11px;")
            else:
                w["rbv"].setText("—")
                w["rbv"].setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")

            moving = movn is not None and int(movn) == 1
            w["movn"].setText("● moving" if moving else "")
            w["movn"].setStyleSheet(
                f"color:{THEME.SUCCESS};font-size:10px;" if moving
                else f"color:transparent;font-size:10px;")

    def _move_veto(self, vname: str):
        """Write the spinbox setpoint to the motor VAL PV."""
        target = self._veto_widgets[vname]["spin"].value()
        ok = self._veto_ctrl.put(f"{vname}_val", target)
        if not ok:
            self._status_label.setText(
                f"{vname} move failed — PV not connected")

    # ------------------------------------------------------------------
    #  Cleanup
    # ------------------------------------------------------------------

    def closeEvent(self, event):
        self._veto_timer.stop()
        if self._any_worker_running():
            self._stop_worker()
            if self._stats_worker:   self._stats_worker.wait(3000)
            if self._instant_worker: self._instant_worker.wait(3000)
        if self._fetcher and self._fetcher.isRunning():
            self._fetcher.wait(1000)
        event.accept()


# ===========================================================================
#  Entry point
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Scintillator–HyCal coincidence rate monitor")
    parser.add_argument("--url", default=DEFAULT_URL,
                        help=f"prad2_server base URL (default: {DEFAULT_URL})")
    parser.add_argument("--theme", choices=available_themes(), default="dark")
    args = parser.parse_args()

    app = QApplication(sys.argv)
    app.setApplicationName("Coincidence Monitor")

    set_theme(args.theme)

    win = MainWindow(server_url=args.url)
    apply_theme_palette(win)
    win.setStyleSheet(themed(f"QMainWindow{{background:{THEME.BG};}}"))
    win.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
