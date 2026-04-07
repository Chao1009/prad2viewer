"""
Shared utilities for the offline HyCal gain-equalization workflow.

Used by:
    gain_eq_analyze.py      — pass 1: ROOT → fit results, append to JSON
    gain_eq_distribution.py — pass 2: visualize the per-channel mean distribution
    gain_eq_propose.py      — pass 3: history → batch VSet JSON for prad2hvd

Unified history JSON format (gain_eq_history_v1)
================================================
One file per equalization campaign. Per channel, an array of iteration
entries — one entry appended each time gain_eq_analyze runs:

    {
        "format": "gain_eq_history_v1",
        "channels": {
            "W3": [
                {
                    "peak_height_mean":   2540.3,
                    "peak_height_sigma":  85.2,
                    "peak_integral_mean": 12500.0,
                    "peak_integral_sigma": 450.0,
                    "count":   12543,
                    "VMon":    1523.9,
                    "VSet":    1525.0,
                    "status":  "OK",
                    "timestamp": "2026-04-06 14:00:00",
                    "input":   "prad_023527.root"
                },
                { ... iter 2 ... }
            ],
            "G1": [ ... ],
            ...
        }
    }
"""

from __future__ import annotations

import fnmatch
import json
import os
from typing import Any, Dict, List, Optional, Tuple


# ============================================================================
#  Path defaults — relative to analysis/scripts/
# ============================================================================

DEFAULT_DAQ_MAP        = "../../database/daq_map.json"
DEFAULT_HYCAL_MODULES  = "../../database/hycal_modules.json"
DEFAULT_VOLTAGE_LIMITS = "../../../prad2hvmon/database/voltage_limits.json"


def resolve_db(path: str) -> str:
    """Resolve a database path relative to this script's directory."""
    if os.path.isabs(path) or os.path.exists(path):
        return path
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, path))


# ============================================================================
#  Database loaders
# ============================================================================

def load_daq_map(path: str = DEFAULT_DAQ_MAP) -> Dict[Tuple[int, int, int], str]:
    """Load daq_map.json → {(crate, slot, channel): module_name}.

    crate is the logical crate index (0..6 for HyCal) as stored by Replay.cpp.
    """
    with open(resolve_db(path)) as f:
        entries = json.load(f)
    out: Dict[Tuple[int, int, int], str] = {}
    for e in entries:
        key = (int(e["crate"]), int(e["slot"]), int(e["channel"]))
        out[key] = e["name"]
    return out


def load_hycal_modules(path: str = DEFAULT_HYCAL_MODULES) -> Dict[str, dict]:
    """Load hycal_modules.json → {module_name: entry_dict}."""
    with open(resolve_db(path)) as f:
        entries = json.load(f)
    return {e["n"]: e for e in entries}


def load_voltage_limits(path: str = DEFAULT_VOLTAGE_LIMITS) -> List[Tuple[str, float]]:
    """Load voltage_limits.json → [(pattern, max_voltage)] in declaration order.

    Returns [] if the file isn't reachable (so propose script can warn instead
    of crashing on machines without prad2hvmon checked out).
    """
    try:
        with open(resolve_db(path)) as f:
            data = json.load(f)
    except FileNotFoundError:
        return []
    return [(item["pattern"], float(item["voltage"])) for item in data["limits"]]


def voltage_limit_for(name: str,
                      limits: List[Tuple[str, float]]) -> float:
    """Pattern match (first wins) → max voltage. 0.0 = unknown."""
    for pattern, vmax in limits:
        if fnmatch.fnmatchcase(name, pattern):
            return vmax
    return 0.0


def filter_modules_by_type(names: List[str],
                           modules: Dict[str, dict],
                           wanted: List[str]) -> List[str]:
    """Return only modules whose 't' field matches one of the wanted prefixes."""
    out = []
    for n in names:
        m = modules.get(n)
        if not m:
            continue
        t = m.get("t", "")
        if any(t == w or t.startswith(w) for w in wanted):
            out.append(n)
    return out


def natural_module_sort_key(name: str):
    """Sort G/W/L/P modules by type prefix then by numeric suffix."""
    if not name:
        return ("", 0)
    prefix = name[0]
    try:
        num = int(name[1:])
    except ValueError:
        num = 0
    return (prefix, num)


# ============================================================================
#  History JSON I/O (the single source of truth across iterations)
# ============================================================================

HISTORY_FORMAT = "gain_eq_history_v1"


def load_history(path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Load history JSON → {channel_name: [iter1, iter2, ...]}.

    Returns an empty dict if the file doesn't exist.
    """
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        data = json.load(f)
    if not isinstance(data, dict):
        return {}
    if data.get("format") != HISTORY_FORMAT:
        raise RuntimeError(
            f"{path}: unexpected format '{data.get('format')}', "
            f"expected '{HISTORY_FORMAT}'")
    channels = data.get("channels", {})
    if not isinstance(channels, dict):
        return {}
    return channels


def save_history(path: str,
                 channels: Dict[str, List[Dict[str, Any]]]) -> None:
    """Save the channels dict back to a history JSON file."""
    with open(path, "w") as f:
        json.dump({
            "format": HISTORY_FORMAT,
            "channels": channels,
        }, f, indent=2)


def latest_iteration(channels: Dict[str, List[Dict[str, Any]]],
                     name: str) -> Optional[Dict[str, Any]]:
    """Return the most recent iteration entry for a channel, or None."""
    entries = channels.get(name, [])
    return entries[-1] if entries else None


def n_iterations(channels: Dict[str, List[Dict[str, Any]]]) -> int:
    """Return the maximum iteration count across all channels."""
    return max((len(v) for v in channels.values()), default=0)


# ============================================================================
#  Batch settings JSON (prad2hvmon_settings_v1) for prad2hvd /api/load_settings
# ============================================================================

def build_batch_settings(timestamp: str,
                         channels: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Wrap a list of channel entries in the prad2hvmon_settings_v1 envelope.

    Each channel dict should have:
        crate (str), slot (int), channel (int), name (str), params (dict)
    """
    return {
        "format": "prad2hvmon_settings_v1",
        "timestamp": timestamp,
        "channels": channels,
    }


# ============================================================================
#  Minimal HV client (HTTP only — no pyepics, no threads)
# ============================================================================

class HVClient:
    """Minimal HTTP client for prad2hvd voltage read/write.

    Standalone — does NOT depend on calibration/gain_scanner.py (which pulls
    in pyepics).

    Endpoints:
        GET  /api/voltage?name=G1   → no auth; returns vset/vmon/crate/slot/channel
        POST /api/voltage           → expert auth; body {name, value}
        POST /api/auth              → returns {granted: <level>}
        POST /api/load_settings     → expert auth; body = batch settings JSON
    """

    def __init__(self, url: str = "http://clonpc19:8765",
                 password: str = "", verbose: bool = False):
        self.url = url.rstrip("/")
        self.password = password
        self.verbose = verbose

    def _get(self, path: str) -> Any:
        import urllib.request
        with urllib.request.urlopen(f"{self.url}{path}", timeout=10) as r:
            return json.loads(r.read())

    def _post(self, path: str, data: dict) -> Any:
        import urllib.request
        body = json.dumps(data).encode()
        req = urllib.request.Request(
            f"{self.url}{path}", data=body, method="POST",
            headers={"Content-Type": "application/json"})
        if self.password:
            req.add_header("X-Auth", self.password)
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read())

    def authenticate(self) -> int:
        if not self.password:
            return 0
        resp = self._post("/api/auth", {"password": self.password})
        granted = int(resp.get("granted", 0))
        if self.verbose:
            print(f"HV: authenticated, level={granted}")
        return granted

    def get_voltage(self, name: str) -> Optional[dict]:
        """Read VSet/VMon and crate addressing for one module by name."""
        import urllib.error
        try:
            return self._get(f"/api/voltage?name={name}")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            raise

    def set_voltage(self, name: str, value: float) -> bool:
        try:
            self._post("/api/voltage", {"name": name, "value": round(value, 2)})
            return True
        except Exception as e:
            print(f"HV: SET {name} failed: {e}")
            return False

    def load_settings(self, settings: Dict[str, Any]) -> dict:
        return self._post("/api/load_settings", settings)
