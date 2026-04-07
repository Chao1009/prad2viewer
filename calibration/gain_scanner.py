"""
Automatic Gain Equalizer for HyCal
===================================
Moves the beam to each crystal module, collects peak height histogram data
from prad2_server, finds the right edge of the Bremsstrahlung spectrum, and
adjusts HV via prad2hvd until the edge aligns with a target ADC value.

Classes:
    SpectrumAnalyzer  — pure analysis (edge finding, voltage step)
    ServerClient      — HTTP client for prad2_server
    HVClient          — WebSocket client for prad2hvd
    GainScanEngine    — orchestrates the full scan loop
"""

from __future__ import annotations

import json
import math
import os
import threading
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from scan_utils import Module, module_to_ptrans, ptrans_in_limits
from scan_epics import epics_move_to, epics_is_moving, epics_stop, SPMG


# ============================================================================
#  Spectrum Analysis
# ============================================================================

class SpectrumAnalyzer:
    """Analyse peak height histograms to find the spectrum right edge."""

    def __init__(self, target_adc: float = 3200.0,
                 bin_step: float = 10.0, bin_min: float = 0.0,
                 smooth_window: int = 5,
                 edge_fraction: float = 0.05,
                 use_log_cumul: bool = True,
                 pedestal_adc: float = 200.0):
        self.target_adc = target_adc
        self.bin_step = bin_step
        self.bin_min = bin_min
        self.smooth_window = smooth_window  # moving-average window
        self.edge_fraction = edge_fraction  # cumulative fraction threshold
        self.use_log_cumul = use_log_cumul  # use log(1+count) for cumulative
        self.pedestal_adc = pedestal_adc    # exclude bins below this ADC

    def find_right_edge(self, bins: List[int]) -> Optional[int]:
        """Find the right-edge bin of the Bremsstrahlung continuum spectrum.

        Algorithm:

        1. **Exclude pedestal**: bins below ``pedestal_adc`` (default 200).
        2. **Smooth** with moving average (window=5).
        3. **Log-scale cumulative from right** (default): accumulate
           ``log(1 + count)`` from the rightmost bin leftward.  The log
           transform compresses high-count body bins so sparse tail bins
           carry more relative weight.  When the cumulative reaches
           ``edge_fraction`` (default 2%) of the log-total, we've found
           the edge.  If ``use_log_cumul`` is False, raw counts are used.
        4. **Confirm falling edge** via smoothed slope.

        Returns the bin index, or None if no edge found.
        """
        n = len(bins)
        if n < self.smooth_window + 3:
            return None

        # step 1: exclude pedestal
        ped_bin = int((self.pedestal_adc - self.bin_min) / self.bin_step) \
            if self.bin_step > 0 else 0
        ped_bin = max(0, min(ped_bin, n))

        # step 2: smooth
        hw = self.smooth_window // 2
        smooth = [0.0] * n
        for i in range(n):
            lo = max(0, i - hw)
            hi = min(n, i + hw + 1)
            smooth[i] = sum(bins[lo:hi]) / (hi - lo)

        # step 3: cumulative from right
        if self.use_log_cumul:
            total = sum(math.log1p(bins[i]) for i in range(ped_bin, n))
            weight_fn = math.log1p
        else:
            total = sum(bins[ped_bin:])
            weight_fn = float
        if total <= 0:
            return None

        threshold = total * self.edge_fraction
        cumul = 0.0
        candidate = None
        for i in range(n - 1, ped_bin - 1, -1):
            cumul += weight_fn(bins[i])
            if cumul >= threshold:
                candidate = i
                break
        if candidate is None:
            return None

        # step 3: confirm falling edge — walk right from candidate to find
        # where the smoothed value starts dropping consistently
        # (candidate might be slightly inside the body; refine outward)
        peak_smooth = max(smooth)
        if peak_smooth <= 0:
            return candidate
        drop_threshold = peak_smooth * 0.05  # 5% of peak = noise floor
        for i in range(candidate, min(n - 1, candidate + 20)):
            if smooth[i] < drop_threshold:
                return max(candidate, i - 1)
            # check if we're on falling slope: this bin < left neighbor
            if i > 0 and smooth[i] < smooth[i - 1] * 0.7:
                return i
        return candidate

    def edge_to_adc(self, bin_index: int) -> float:
        """Convert bin index to ADC value (centre of bin)."""
        return self.bin_min + (bin_index + 0.5) * self.bin_step

    def compute_voltage_step(self, edge_adc: float,
                             current_vset: float = 0.0) -> float:
        """Compute HV adjustment based on the ADC difference from target.

        | ADC difference | Step  |
        |----------------|-------|
        | > 1000         | 50 V  |
        | 500 – 1000     | 30 V  |
        | 200 – 500      | 20 V  |
        | 100 – 200      |  5 V  |
        | < 100          | diff/20 V |

        Sign: positive = edge too low → increase voltage.
        """
        diff = self.target_adc - edge_adc  # positive = need more gain
        sign = 1 if diff > 0 else -1
        ad = abs(diff)
        if ad > 1000:    dv = 50.0
        elif ad > 500:   dv = 30.0
        elif ad > 200:   dv = 20.0
        elif ad > 100:   dv = 5.0
        else:            dv = ad / 20.0
        return round(sign * dv, 1)


# ============================================================================
#  prad2_server HTTP Client
# ============================================================================

class ServerClient:
    """HTTP client for prad2_server histogram and occupancy APIs."""

    def __init__(self, url: str = "http://clondaq6:5051", log_fn=None,
                 read_only: bool = False):
        self.url = url.rstrip("/")
        self._key_map: Dict[str, str] = {}
        self._log = log_fn or (lambda msg, **kw: None)
        self._read_only = read_only

    def _get(self, path: str) -> Any:
        import urllib.request
        with urllib.request.urlopen(f"{self.url}{path}", timeout=10) as r:
            return json.loads(r.read())

    def _post(self, path: str) -> Any:
        import urllib.request
        req = urllib.request.Request(f"{self.url}{path}", method="POST", data=b"")
        with urllib.request.urlopen(req, timeout=10) as r:
            return json.loads(r.read())

    def get_config(self) -> dict:
        self._log("Server: GET /api/config")
        return self._get("/api/config")

    def clear_histograms(self):
        if self._read_only:
            self._log("Server: POST /api/hist/clear [BLOCKED — read-only]", level="warn")
            return {}
        self._log("Server: POST /api/hist/clear")
        return self._post("/api/hist/clear")

    def get_occupancy(self) -> dict:
        return self._get("/api/occupancy")

    def get_height_histogram(self, key: str, quiet: bool = False) -> dict:
        return self._get(f"/api/heighthist/{key}")

    def build_key_map(self) -> Dict[str, str]:
        """Build module_name -> histogram_key mapping from server config."""
        cfg = self.get_config()
        daq = cfg.get("daq", [])
        crate_roc = cfg.get("crate_roc", {})
        key_map: Dict[str, str] = {}
        for entry in daq:
            name = entry.get("name", "")
            crate = str(entry.get("crate", ""))
            slot = entry.get("slot", 0)
            ch = entry.get("channel", 0)
            roc = crate_roc.get(crate, crate)
            key_map[name] = f"{roc}_{slot}_{ch}"
        self._key_map = key_map
        return key_map

    def get_module_counts(self, module_name: str) -> int:
        """Get hit count for a single module from occupancy."""
        key = self._key_map.get(module_name, "")
        if not key:
            return 0
        occ = self.get_occupancy()
        return occ.get("occ", {}).get(key, 0)


# ============================================================================
#  prad2hvd HTTP Client
# ============================================================================

class HVClient:
    """HTTP client for prad2hvd voltage control.

    Stateless HTTP — no WebSocket, no background threads, no broken pipes.

    API:
        GET  /api/voltage?name=W1124           → read voltage (no auth)
        POST /api/voltage  {name, value}        → set voltage (X-Auth header)
        POST /api/auth     {password}           → test auth (returns granted level)
    """

    def __init__(self, url: str = "http://clonpc19:8765", log_fn=None,
                 read_only: bool = False):
        self.url = url.rstrip("/")
        self._log = log_fn or (lambda msg, **kw: None)
        self._read_only = read_only
        self._password = ""

    def connect(self, password: str = ""):
        """Verify connectivity and authenticate."""
        import urllib.request, urllib.error
        self._password = password
        self._log(f"HV: connecting to {self.url}")
        # test connectivity — 404 is OK (means server is up, module not found)
        try:
            self._http_get("/api/voltage?name=_test_")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                pass  # server is reachable, module just doesn't exist
            elif self._read_only:
                self._log(f"HV: not reachable (HTTP {e.code}) — running read-only", level="warn")
                return
            else:
                raise
        except Exception as e:
            if self._read_only:
                self._log(f"HV: not reachable ({e}) — running read-only", level="warn")
                return
            raise
        # test auth if password given
        if password:
            self._log("HV: authenticating (Expert)")
            resp = self._http_post("/api/auth", {"password": password})
            granted = resp.get("granted", 0)
            if granted < 2:
                raise RuntimeError(f"HV authentication failed (granted={granted})")
            self._log(f"HV: authenticated OK (level {granted})")
        self._log("HV: connected")

    def _http_get(self, path: str) -> Any:
        import urllib.request
        with urllib.request.urlopen(f"{self.url}{path}", timeout=10) as r:
            return json.loads(r.read())

    def _http_post(self, path: str, data: dict) -> Any:
        import urllib.request
        body = json.dumps(data).encode()
        req = urllib.request.Request(
            f"{self.url}{path}", data=body, method="POST",
            headers={"Content-Type": "application/json"})
        if self._password:
            req.add_header("X-Auth", self._password)
        with urllib.request.urlopen(req, timeout=10) as r:
            return json.loads(r.read())

    def get_voltage(self, name: str) -> Optional[dict]:
        """Read voltage info for a module by name."""
        import urllib.error
        try:
            resp = self._http_get(f"/api/voltage?name={name}")
            self._log(f"HV: GET {name} → vset={resp.get('vset'):.1f} "
                      f"vmon={resp.get('vmon'):.1f} limit={resp.get('limit'):.0f}")
            return resp
        except urllib.error.HTTPError as e:
            if e.code == 404:
                self._log(f"HV: GET {name} → not found", level="warn")
            else:
                self._log(f"HV: GET {name} → HTTP {e.code}", level="error")
            return None
        except Exception as e:
            self._log(f"HV: GET {name} → failed: {e}", level="error")
            return None

    def set_voltage(self, name: str, value: float,
                    old_value: float = 0.0) -> bool:
        """Set voltage by module name. Returns True on success."""
        dv = value - old_value
        direction = "increase" if dv > 0 else "decrease" if dv < 0 else "no change"
        if self._read_only:
            self._log(f"HV: SET {name} {direction} {old_value:.1f} → {value:.1f} V "
                      f"(ΔV={dv:+.1f}) [BLOCKED — read-only]", level="warn")
            return True
        self._log(f"HV: SET {name} {direction} {old_value:.1f} → {value:.1f} V (ΔV={dv:+.1f})")
        import urllib.error
        try:
            resp = self._http_post("/api/voltage",
                                    {"name": name, "value": round(value, 2)})
            self._log(f"HV: SET {name} → {resp.get('status', 'ok')}")
            return True
        except urllib.error.HTTPError as e:
            if e.code == 403:
                self._log(f"HV: SET {name} → forbidden (bad password?)", level="error")
            else:
                self._log(f"HV: SET {name} → HTTP {e.code}", level="error")
            return False
        except Exception as e:
            self._log(f"HV: SET {name} → failed: {e}", level="error")
            return False

    def close(self):
        pass  # stateless — nothing to close


# ============================================================================
#  Gain Scan Engine
# ============================================================================

class GainScanState:
    IDLE       = "IDLE"
    MOVING     = "MOVING"
    COLLECTING = "COLLECTING"
    ANALYZING  = "ANALYZING"
    ADJUSTING  = "ADJUSTING"
    CONVERGED  = "CONVERGED"
    FAILED     = "FAILED"
    COMPLETED  = "COMPLETED"


class GainScanEngine:
    """Orchestrates the automatic gain equalization scan.

    Each module is processed as a self-contained step: fresh server/HV
    connections are created, the module is equalized (or fails), then
    connections are closed.  This avoids stale WebSocket state between
    modules.
    """

    # defaults (configurable before start)
    target_adc: float = 3200.0
    min_counts: int = 10000
    max_iterations: int = 8
    convergence_tol: float = 50.0    # ADC units
    hv_settle_time: float = 10.0     # seconds (legacy, unused)
    pos_threshold: float = 0.5       # mm
    beam_threshold: float = 0.3      # nA
    collect_poll_sec: float = 2.0    # occupancy poll interval
    move_timeout: float = 900.0      # seconds (15 min — long y-axis traversals)
    edge_adc_min: float = 500.0      # reject edges below this
    edge_adc_max: float = 3900.0     # reject edges above this
    use_log_y: bool = True           # log y scale for report snapshots

    # Linear-response model
    RESPONSE_CHANGE_TOL: float = 0.20    # max relative change in response_constant before fallback

    # VMon settle wait
    VMON_INITIAL_WAIT: float = 3.0       # seconds — mandatory wait before checking
    VMON_POLL_INTERVAL: float = 1.5      # seconds between HV reads
    VMON_TIMEOUT: float = 40.0           # total timeout
    VMON_RTOL: float = 0.10              # |Δvmon - Δvset| / |Δvset| relative tolerance
    VMON_ATOL: float = 2.0               # absolute tolerance (V)

    def __init__(self, motor_ep,
                 server_url: str, hv_url: str, hv_password: str,
                 read_only: bool,
                 modules: List[Module], log_fn,
                 key_map: Dict[str, str]):
        self.ep = motor_ep
        self._server_url = server_url
        self._hv_url = hv_url
        self._hv_password = hv_password
        self._read_only = read_only
        # ``modules`` is already the ordered path (caller-prepared)
        self.path = list(modules)
        self.log = log_fn
        self.key_map = key_map
        self.analyzer = SpectrumAnalyzer(target_adc=self.target_adc)
        # current step's connections (created per module, visible to UI)
        self.server: Optional[ServerClient] = None
        self.hv: Optional[HVClient] = None

        self.state = GainScanState.IDLE
        self.current_idx = 0
        self.current_iteration = 0
        self.last_edge_adc: Optional[float] = None
        self.last_edge_bin: Optional[int] = None
        self.last_dv: Optional[float] = None
        self.last_bins: List[int] = []
        self.last_vset: Optional[float] = None
        self.last_vmon: Optional[float] = None
        self.module_counts = 0
        self.collect_rate: float = 0.0  # Hz, updated during collection

        # linear response model state (per module)
        self._linear_response: bool = False
        self._response_constant: Optional[float] = None  # ADC per V slope
        self._prev_vmon: Optional[float] = None
        self._prev_edge: Optional[float] = None
        self._prev_vset: Optional[float] = None

        # per-iteration history for the current module (for screenshot report)
        self.iteration_history: List[Dict] = []
        # directory for saving reports
        self.report_dir: str = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "logs")

        self.converged: set = set()
        self.failed: set = set()

        self._thread: Optional[threading.Thread] = None
        self._stop = threading.Event()
        self._skip = threading.Event()
        self._redo = threading.Event()
        self._paused = False
        self._has_run = False  # True after first start, for resume detection

    @property
    def current_module(self) -> Optional[Module]:
        if 0 <= self.current_idx < len(self.path):
            return self.path[self.current_idx]
        return None

    def start(self, start_idx: int = 0, count: int = 0):
        """Start or resume the gain scan.

        If resuming after FAILED or STOP: retries from ``current_idx``
        (the failed/stopped module), keeping previous converged/failed results.
        Otherwise: fresh start from ``start_idx``, clearing all results.
        """
        if self._thread and self._thread.is_alive():
            return
        self._stop.clear()
        self._skip.clear()
        self._redo.clear()
        self._paused = False

        resuming = self._has_run and self.state in (GainScanState.FAILED, GainScanState.IDLE)
        if resuming and 0 <= self.current_idx < len(self.path):
            # resume: retry from current_idx, keep _end_idx and results
            self._start_idx = self.current_idx
            mod = self.path[self.current_idx]
            self.log(f"Resuming from {mod.name} [{self._start_idx + 1}/{len(self.path)}]")
        else:
            # fresh start
            self.current_idx = start_idx
            self._start_idx = start_idx
            self._end_idx = min(start_idx + count, len(self.path)) if count > 0 else len(self.path)
            self.converged.clear()
            self.failed.clear()

        self._has_run = True
        self.state = GainScanState.MOVING
        self.analyzer.target_adc = self.target_adc
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def pause(self):
        if self.state not in (GainScanState.IDLE, GainScanState.COMPLETED):
            self._paused = True
            self.log("Gain scan paused", level="warn")

    def resume(self):
        if self._paused:
            self._paused = False
            self.log("Gain scan resumed")

    def stop(self):
        self._stop.set()
        self._skip.set()
        self._paused = False
        epics_stop(self.ep)
        self.log("Gain scan stop requested", level="warn")

    def skip_module(self):
        self._skip.set()

    def redo_module(self):
        self._redo.set()

    # -- main loop ----------------------------------------------------------

    def _run(self):
        n = self._end_idx - self._start_idx
        self.log(f"Gain scan started: {n} modules, target ADC {self.target_adc:.0f}")

        # both clients are stateless HTTP — create once, reuse for all modules
        self.server = ServerClient(self._server_url,
                                   log_fn=self.log, read_only=self._read_only)
        self.hv = HVClient(self._hv_url,
                           log_fn=self.log, read_only=self._read_only)
        try:
            self.hv.connect(password=self._hv_password)
        except Exception as e:
            self.log(f"HV connection failed: {e}", level="error")
            self.state = GainScanState.FAILED
            return

        try:
            for i in range(self._start_idx, self._end_idx):
                if self._stop.is_set():
                    break

                # -- reset per-module state --
                self.current_idx = i
                self.current_iteration = 0
                self.last_edge_adc = None
                self.last_edge_bin = None
                self.last_dv = None
                self.last_bins = []
                self.last_vset = None
                self.last_vmon = None
                self.module_counts = 0
                self.collect_rate = 0.0
                self.iteration_history = []
                # reset linear response state
                self._linear_response = False
                self._response_constant = None
                self._prev_vmon = None
                self._prev_edge = None
                self._prev_vset = None
                mod = self.path[i]

                self.log(f"── [{i+1}/{len(self.path)}] {mod.name} ──")

                self._process_module(i, mod)
                self._skip.clear()

                # if this module failed, pause the scan for user intervention
                if i in self.failed:
                    self.log(f"{mod.name}: STOPPED — click Resume to retry this module, or Reset to clear and pick a new starting point",
                             level="error")
                    return

        finally:
            if self._stop.is_set():
                self.state = GainScanState.IDLE
                mod = self.path[self.current_idx] if 0 <= self.current_idx < len(self.path) else None
                if mod:
                    self.log(f"Stopped at {mod.name} — click Resume to retry, "
                             f"or Reset to clear and pick a new starting point",
                             level="warn")
                else:
                    self.log("Gain scan stopped", level="warn")
            elif self.state != GainScanState.FAILED:
                self.state = GainScanState.COMPLETED
                self.log(f"Gain scan COMPLETE: {len(self.converged)} converged, "
                         f"{len(self.failed)} failed", level="warn")

    def _process_module(self, i: int, mod: Module):
        """Process one module: move → (collect → analyze → adjust) × N."""
        px, py = module_to_ptrans(mod.x, mod.y)

        # -- move --
        self.state = GainScanState.MOVING
        # estimate move time from current position and motor velocities
        rx = self.ep.get("x_rbv", 0.0) or 0.0
        ry = self.ep.get("y_rbv", 0.0) or 0.0
        vx = self.ep.get("x_velo", 50.0) or 50.0
        vy = self.ep.get("y_velo", 5.0) or 5.0
        eta = max(abs(px - rx) / max(vx, 0.1), abs(py - ry) / max(vy, 0.1))
        timeout = eta + 300.0  # ETA + 5 min margin
        self.log(f"Moving to {mod.name}  ptrans({px:.3f}, {py:.3f})  "
                 f"ETA={eta:.0f}s  timeout={timeout:.0f}s")
        if not epics_move_to(self.ep, px, py):
            self.log(f"SKIPPED {mod.name}: outside limits", level="warn")
            return
        if not self._wait_move_done(timeout):
            return

        # check DAQ key
        key = self.key_map.get(mod.name)
        if not key:
            self.log(f"SKIPPED {mod.name}: no DAQ mapping", level="warn")
            return

        # -- iterate: collect → analyze → adjust --
        iteration = 0
        while iteration < self.max_iterations:
            if self._stop.is_set():
                return
            if self._skip.is_set():
                self._skip.clear()
                self.log(f"Module {mod.name} skipped")
                return
            if self._redo.is_set():
                self._redo.clear()
                iteration = 0
                self.iteration_history = []
                self.last_edge_adc = None
                self.last_edge_bin = None
                self.last_dv = None
                self.last_bins = []
                self.log(f"Module {mod.name} redo — restarting iterations")
                continue
            iteration += 1
            self.current_iteration = iteration

            # collect
            self.state = GainScanState.COLLECTING
            self.module_counts = 0
            try:
                self.server.clear_histograms()
            except Exception as e:
                self.log(f"Server error (clear): {e}", level="error")
                self._mark_failed(i, mod); return
            if not self._wait_for_counts(mod, key):
                if self._stop.is_set():
                    return
                self._mark_failed(i, mod); return

            # analyze
            self.state = GainScanState.ANALYZING
            try:
                hist = self.server.get_height_histogram(key)
            except Exception as e:
                self.log(f"Server error (hist): {e}", level="error")
                self._mark_failed(i, mod); return
            bins = hist.get("bins", [])
            self.last_bins = bins
            edge = self.analyzer.find_right_edge(bins)
            self.last_edge_bin = edge
            if edge is None:
                self.log(f"{mod.name} iter {iteration+1}: no edge found", level="error")
                self._mark_failed(i, mod); return
            edge_adc = self.analyzer.edge_to_adc(edge)
            self.last_edge_adc = edge_adc

            # read HV (single call per iteration — for display + voltage step)
            info = self.hv.get_voltage(mod.name)
            if info is None:
                if self._read_only:
                    info = {"vset": 1000.0, "vmon": 0.0, "limit": 2000.0}
                else:
                    self.log(f"{mod.name}: HV read failed", level="error")
                    self._mark_failed(i, mod); return
            self.last_vset = info.get("vset")
            self.last_vmon = info.get("vmon")

            # record iteration snapshot
            self.iteration_history.append({
                "iteration": iteration + 1,
                "bins": list(bins),
                "edge_bin": edge,
                "edge_adc": edge_adc,
                "vset": self.last_vset,
                "vmon": self.last_vmon,
                "time": datetime.now().strftime("%H:%M:%S"),
            })

            if edge_adc < self.edge_adc_min or edge_adc > self.edge_adc_max:
                self.log(f"{mod.name} iter {iteration+1}: edge {edge_adc:.0f} "
                         f"out of range [{self.edge_adc_min:.0f}, {self.edge_adc_max:.0f}]",
                         level="error")
                self._mark_failed(i, mod); return

            # check convergence
            if abs(edge_adc - self.target_adc) <= self.convergence_tol:
                self.converged.add(i)
                self.state = GainScanState.CONVERGED
                self.log(f"{mod.name}: CONVERGED at {edge_adc:.0f} (iter {iteration+1})")
                self._save_module_report(mod, "success")
                return

            # adjust HV
            self.state = GainScanState.ADJUSTING
            current_v = info.get("vset", 0)
            current_vmon = info.get("vmon", current_v)
            limit_v = info.get("limit", 99999)

            # -- update linear response model --
            mode_tag = "lookup"
            if (self._prev_vmon is not None and self._prev_edge is not None):
                dvmon_meas = current_vmon - self._prev_vmon
                if abs(dvmon_meas) > 1e-6:
                    new_response = (edge_adc - self._prev_edge) / dvmon_meas
                    if new_response <= 0:
                        # non-physical
                        self._linear_response = False
                        mode_tag = (f"fallback: response={new_response:.2f} "
                                    f"ADC/V (non-physical)")
                    elif self._linear_response and self._response_constant:
                        rel = abs(new_response - self._response_constant) \
                            / abs(self._response_constant)
                        if rel > self.RESPONSE_CHANGE_TOL:
                            self._linear_response = False
                            mode_tag = (f"fallback: response "
                                        f"{self._response_constant:.2f}→"
                                        f"{new_response:.2f} ADC/V "
                                        f"(Δ{rel*100:.0f}%)")
                        else:
                            self._response_constant = new_response
                            self._linear_response = True
                    else:
                        # first slope measurement (or coming back from fallback)
                        self._response_constant = new_response
                        self._linear_response = True

            # -- compute dV --
            if self._linear_response and self._response_constant \
                    and self._response_constant > 0:
                dv = (self.target_adc - edge_adc) / self._response_constant
                # round to 0.1 V for HV crate granularity
                dv = round(dv, 1)
                if mode_tag == "lookup":
                    mode_tag = f"linear: response={self._response_constant:.2f} ADC/V"
            else:
                dv = self.analyzer.compute_voltage_step(edge_adc, current_v)

            self.last_dv = dv
            new_v = current_v + dv

            if new_v > limit_v:
                self.log(f"{mod.name}: would exceed limit ({new_v:.1f} > {limit_v:.1f})",
                         level="error")
                self._mark_failed(i, mod); return

            self.log(f"{mod.name} iter {iteration+1}: edge={edge_adc:.0f} "
                     f"ΔV={dv:+.1f} ({current_v:.1f}→{new_v:.1f})  [{mode_tag}]")

            # save prev state for next iteration's slope calc
            self._prev_vmon = current_vmon
            self._prev_edge = edge_adc
            self._prev_vset = current_v

            if not self.hv.set_voltage(mod.name, new_v, old_value=current_v):
                self.log(f"{mod.name}: HV set failed", level="error")
                self._mark_failed(i, mod); return

            # wait for VMon to catch up to the new VSet
            if not self._wait_vmon_settle(mod, dv, current_vmon):
                self._mark_failed(i, mod); return
            if self._stop.is_set():
                return

        # max iterations exhausted
        if i not in self.converged and i not in self.failed:
            self.log(f"{mod.name}: max iterations reached "
                     f"(last edge={self.last_edge_adc:.0f})", level="warn")
            self._mark_failed(i, mod)

    def _mark_failed(self, idx: int, mod: Module):
        self.failed.add(idx)
        self.state = GainScanState.FAILED
        self._save_module_report(mod, "failure")

    def _save_module_report(self, mod: Module, status: str):
        """Save a vertically concatenated histogram screenshot for a module.

        Each iteration is drawn as a small histogram panel, stacked top to
        bottom in time order.  Filename: GE_{time}_{name}_{status}.png
        """
        history = self.iteration_history
        if not history:
            return
        try:
            from PyQt6.QtGui import QImage, QPainter, QColor, QPen, QFont
            from PyQt6.QtCore import Qt, QRectF
        except ImportError:
            self.log("Cannot save report: PyQt6 not available in thread", level="warn")
            return

        PANEL_W, PANEL_H = 600, 160
        PAD_L, PAD_R, PAD_T, PAD_B = 50, 12, 32, 20
        img_h = PANEL_H * len(history)
        img = QImage(PANEL_W, img_h, QImage.Format.Format_RGB32)
        img.fill(QColor("#0d1117"))

        p = QPainter(img)
        target_adc = self.target_adc
        bin_step = self.analyzer.bin_step
        bin_min = self.analyzer.bin_min

        for panel_idx, snap in enumerate(history):
            y0 = panel_idx * PANEL_H
            bins = snap["bins"]
            n = len(bins)
            if n == 0:
                continue
            vmax = max(bins) if bins else 1
            if vmax == 0: vmax = 1

            pw = PANEL_W - PAD_L - PAD_R
            ph = PANEL_H - PAD_T - PAD_B

            # header
            p.setPen(QColor("#58a6ff"))
            p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
            total_counts = sum(bins)
            mode = "LOG" if self.analyzer.use_log_cumul else "LIN"
            frac_pct = self.analyzer.edge_fraction * 100
            title = (f"{mod.name}  iter {snap['iteration']}  "
                     f"N={total_counts}  "
                     f"({snap['edge_adc']:.0f}:{frac_pct:.0f}%-{mode})")
            p.drawText(QRectF(PAD_L, y0 + 4, pw, PAD_T - 4),
                       Qt.AlignmentFlag.AlignLeft, title)
            p.setPen(QColor("#8b949e"))
            p.setFont(QFont("Consolas", 9))
            info = f"{snap['time']}"
            if snap.get("vmon") is not None:
                info += f"  VMon={snap['vmon']:.1f}"
            if snap.get("vset") is not None:
                info += f"  VSet={snap['vset']:.1f}"
            p.drawText(QRectF(PAD_L, y0 + 4, pw, PAD_T - 4),
                       Qt.AlignmentFlag.AlignRight, info)

            # axes
            ax, ay = PAD_L, y0 + PAD_T
            p.setPen(QPen(QColor("#30363d"), 1))
            p.drawLine(ax, ay, ax, ay + ph)
            p.drawLine(ax, ay + ph, ax + pw, ay + ph)

            # bars (log or linear y)
            import math as _rmath
            use_log = self.use_log_y
            log_vmax = _rmath.log10(max(vmax, 1)) if use_log else 0
            bar_w = pw / n
            p.setPen(Qt.PenStyle.NoPen)
            for bi, v in enumerate(bins):
                if v <= 0: continue
                if use_log:
                    frac = _rmath.log10(v) / log_vmax if log_vmax > 0 else 0
                else:
                    frac = v / vmax
                bh = frac * ph
                bx = ax + bi * bar_w
                by = ay + ph - bh
                p.fillRect(QRectF(bx, by, max(bar_w - 0.3, 0.3), bh),
                           QColor("#58a6ff"))

            # target line (red dashed)
            target_bin = int((target_adc - bin_min) / bin_step) if bin_step > 0 else -1
            if 0 <= target_bin < n:
                tx = ax + (target_bin + 0.5) * bar_w
                p.setPen(QPen(QColor("#f85149"), 1.5, Qt.PenStyle.DashLine))
                p.drawLine(int(tx), ay, int(tx), ay + ph)

            # edge line (green)
            edge_bin = snap.get("edge_bin")
            if edge_bin is not None and 0 <= edge_bin < n:
                ex = ax + (edge_bin + 0.5) * bar_w
                p.setPen(QPen(QColor("#3fb950"), 2))
                p.drawLine(int(ex), ay, int(ex), ay + ph)

            # y-axis grid
            if use_log:
                decade = 10
                while decade < vmax:
                    gf = _rmath.log10(decade) / log_vmax
                    gy = ay + ph - gf * ph
                    p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
                    p.drawLine(ax + 1, int(gy), ax + pw, int(gy))
                    p.setPen(QColor("#8b949e"))
                    p.setFont(QFont("Consolas", 7))
                    p.drawText(QRectF(0, gy - 5, PAD_L - 4, 10),
                               Qt.AlignmentFlag.AlignRight, f"{decade}")
                    decade *= 10
            else:
                for gi in range(1, 5):
                    gf = gi / 5
                    gy = ay + ph - gf * ph
                    gval = int(vmax * gf)
                    p.setPen(QPen(QColor("#21262d"), 1, Qt.PenStyle.DotLine))
                    p.drawLine(ax + 1, int(gy), ax + pw, int(gy))
                    p.setPen(QColor("#8b949e"))
                    p.setFont(QFont("Consolas", 7))
                    p.drawText(QRectF(0, gy - 5, PAD_L - 4, 10),
                               Qt.AlignmentFlag.AlignRight, f"{gval}")

            # separator
            p.setPen(QPen(QColor("#30363d"), 1))
            p.drawLine(0, y0 + PANEL_H - 1, PANEL_W, y0 + PANEL_H - 1)

        p.end()

        os.makedirs(self.report_dir, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"GE_{ts}_{mod.name}_{status}.png"
        path = os.path.join(self.report_dir, fname)
        img.save(path)
        self.log(f"Report saved: {fname}")

    def _wait_move_done(self, timeout: float = None) -> bool:
        if timeout is None:
            timeout = self.move_timeout
        t0 = time.time()
        while not self._stop.is_set():
            self._check_paused()
            if self._stop.is_set():
                return False
            if not epics_is_moving(self.ep):
                return True
            if time.time() - t0 > timeout:
                self.log(f"MOVE TIMEOUT after {timeout:.0f}s", level="error")
                return False
            time.sleep(0.1)
        return False

    LOW_RATE_THRESHOLD = 10.0  # Hz — warn if collection rate drops below this
    LOW_RATE_RESET_POLLS = 10  # reset data after this many consecutive low-rate polls

    def _wait_for_counts(self, mod: Module, key: str) -> bool:
        """Poll occupancy until the target module has min_counts hits."""
        retries = 0
        prev_counts = 0
        prev_time = time.time()
        low_rate_streak = 0
        self.collect_rate = 0.0
        while not self._stop.is_set() and not self._skip.is_set() and not self._redo.is_set():
            self._check_paused()
            if self._stop.is_set():
                return False
            # beam trip check
            if self.beam_threshold > 0:
                bc = self.ep.get("beam_cur", None)
                if bc is not None and bc < self.beam_threshold:
                    self.log(f"BEAM TRIP: {bc:.3f} nA — waiting", level="warn")
                    while not self._stop.is_set():
                        bc2 = self.ep.get("beam_cur", 0.0)
                        if bc2 is not None and bc2 >= self.beam_threshold:
                            break
                        time.sleep(0.5)
                    if self._stop.is_set():
                        return False
                    self.log("BEAM RECOVERED — restarting collection", level="warn")
                    try:
                        self.server.clear_histograms()
                    except Exception:
                        pass
                    prev_counts = 0; prev_time = time.time(); low_rate_streak = 0
            try:
                hist = self.server.get_height_histogram(key, quiet=True)
                counts = sum(hist.get("bins", []))
                self.module_counts = counts
                # compute rate
                now = time.time()
                dt = now - prev_time
                if dt > 0.5:
                    self.collect_rate = (counts - prev_counts) / dt
                    if self.collect_rate < self.LOW_RATE_THRESHOLD and counts > 0:
                        low_rate_streak += 1
                        if low_rate_streak == 1:
                            self.log(f"{mod.name}: low rate {self.collect_rate:.1f} Hz "
                                     f"(< {self.LOW_RATE_THRESHOLD:.0f} Hz)", level="warn")
                        if low_rate_streak >= self.LOW_RATE_RESET_POLLS:
                            self.log(f"{mod.name}: low rate persisted for "
                                     f"{low_rate_streak} polls — resetting data",
                                     level="warn")
                            try:
                                self.server.clear_histograms()
                            except Exception:
                                pass
                            prev_counts = 0; low_rate_streak = 0
                    else:
                        low_rate_streak = 0
                    prev_counts = counts; prev_time = now
                if counts >= self.min_counts:
                    return True
                retries = 0
            except Exception as e:
                retries += 1
                if retries > 3:
                    self.log(f"Server unreachable: {e}", level="error")
                    return False
            time.sleep(self.collect_poll_sec)
        return False

    def _check_paused(self):
        while self._paused and not self._stop.is_set():
            time.sleep(0.1)

    def _wait_paused(self, seconds: float):
        """Sleep for *seconds*, respecting pause, stop, redo, and skip."""
        end = time.time() + seconds
        while time.time() < end:
            if self._stop.is_set() or self._redo.is_set() or self._skip.is_set():
                return
            self._check_paused()
            time.sleep(min(0.2, end - time.time()))

    def _wait_vmon_settle(self, mod: Module, expected_dvset: float,
                          prev_vmon: float) -> bool:
        """Wait until VMon catches up to the new VSet.

        Always waits at least VMON_INITIAL_WAIT seconds (HV ramp delay),
        then polls every VMON_POLL_INTERVAL until the agreement check passes
        or VMON_TIMEOUT is reached.

        Agreement: |delta_vmon - expected_dvset| <= max(VMON_RTOL*|dvset|, VMON_ATOL)

        Returns False on timeout (caller should mark module failed).
        Returns True on success or interrupt (stop/skip/redo).
        """
        # mandatory initial wait — let HV ramp begin
        self._wait_paused(self.VMON_INITIAL_WAIT)
        if self._stop.is_set() or self._redo.is_set() or self._skip.is_set():
            return True

        tol = max(self.VMON_RTOL * abs(expected_dvset), self.VMON_ATOL)
        deadline = time.time() + self.VMON_TIMEOUT - self.VMON_INITIAL_WAIT
        dvmon = 0.0
        vmon = prev_vmon
        while time.time() < deadline:
            if self._stop.is_set() or self._redo.is_set() or self._skip.is_set():
                return True
            self._check_paused()
            info = self.hv.get_voltage(mod.name)
            if info is not None:
                vmon = info.get("vmon", 0)
                dvmon = vmon - prev_vmon
                if abs(dvmon - expected_dvset) <= tol:
                    self.last_vmon = vmon
                    self.log(f"{mod.name}: VMon settled at {vmon:.1f} V "
                             f"(Δ={dvmon:+.1f}, target Δ={expected_dvset:+.1f})")
                    return True
            time.sleep(self.VMON_POLL_INTERVAL)
        self.log(f"{mod.name}: VMon did not settle within "
                 f"{self.VMON_TIMEOUT:.0f}s "
                 f"(target Δ={expected_dvset:+.1f}, got Δ={dvmon:+.1f})",
                 level="error")
        return False
