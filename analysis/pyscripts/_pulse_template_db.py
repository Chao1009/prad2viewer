"""
_pulse_template_db.py — load a per-channel pulse-shape JSON (the output of
fit_pulse_template.py) and classify each channel into one of four states
that the pile-up deconvolver can act on:

    state              meaning
    ----------------   --------------------------------------------------
    "good"             channel passed every gate; use its own template
    "fallback_global"  channel didn't pass, but config asked us to fall
                       back to a global-median template
    "no_template"      channel name absent from JSON or no template available
    "bad"              channel present but failed quality gate (and no
                       fallback requested)

The class is constructed once per run; lookups are O(1).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Template:
    """Per-channel pulse template parameters (medians from fit_pulse_template.py)."""
    name:        str       # channel name as it appears in the JSON
    tau_r_ns:    float
    tau_f_ns:    float
    chi2_med:    float
    n_used:      int
    is_global:   bool      # True for the synthesized global-median template


@dataclass(frozen=True)
class DeconvConfig:
    """Knobs that drive template selection and deconv behaviour.  Source
    of truth is `fadc250_waveform.analyzer.nnls_deconv` in
    `database/daq_config.json` — see `load_nnls_deconv_config()`.
    The `*_override` fields override the template JSON's own
    min_pulses / chi2_max gates when non-None."""
    enabled:                     bool  = False
    template_file:               str   = ""
    fallback_to_global_template: bool  = False
    apply_to_all_peaks:          bool  = False
    tau_r_range_ns:  Tuple[float, float] = (0.5, 10.0)
    tau_f_range_ns:  Tuple[float, float] = (2.0, 100.0)
    cond_number_max: float = 1.0e6
    pre_samples:     int = 8
    post_samples:    int = 40
    min_pulses_override: Optional[int]   = None
    chi2_max_override:   Optional[float] = None


def load_nnls_deconv_config(daq_config_path: Optional[str]) -> DeconvConfig:
    """Read `fadc250_waveform.analyzer.nnls_deconv` out of a daq_config
    JSON.  Empty path or missing block → defaults (deconv disabled)."""
    if not daq_config_path:
        return DeconvConfig()
    p = Path(daq_config_path)
    if not p.is_file():
        raise FileNotFoundError(f"daq_config not found: {p}")
    with open(p, "r", encoding="utf-8") as f:
        j = json.load(f)
    sec = (j.get("fadc250_waveform", {})
            .get("analyzer", {})
            .get("nnls_deconv", {}))
    if not sec:
        # Block absent → deconv simply disabled.
        return DeconvConfig()

    known = {"enabled", "template_file",
             "fallback_to_global_template", "apply_to_all_peaks",
             "tau_r_range_ns", "tau_f_range_ns", "cond_number_max",
             "pre_samples", "post_samples",
             "min_pulses_override", "chi2_max_override"}
    payload = {k: v for k, v in sec.items() if not k.startswith("_")}
    extra = set(payload) - known
    if extra:
        print(f"[WARN] nnls_deconv: unknown keys ignored: {sorted(extra)}",
              flush=True)
    return DeconvConfig(
        enabled                     = bool(payload.get("enabled", False)),
        template_file               = str(payload.get("template_file", "")),
        fallback_to_global_template = bool(payload.get("fallback_to_global_template", False)),
        apply_to_all_peaks          = bool(payload.get("apply_to_all_peaks", False)),
        tau_r_range_ns  = tuple(payload.get("tau_r_range_ns",  (0.5, 10.0))),
        tau_f_range_ns  = tuple(payload.get("tau_f_range_ns",  (2.0, 100.0))),
        cond_number_max = float(payload.get("cond_number_max", 1.0e6)),
        pre_samples     = int(payload.get("pre_samples",  8)),
        post_samples    = int(payload.get("post_samples", 40)),
        min_pulses_override = payload.get("min_pulses_override"),
        chi2_max_override   = payload.get("chi2_max_override"),
    )


# ---------------------------------------------------------------------------
# TemplateDB
# ---------------------------------------------------------------------------

class TemplateDB:
    """Loaded view of a fit_pulse_template.py JSON, with per-channel
    classification and an on-demand global-median template for fallback."""

    def __init__(self, json_path: str, cfg: DeconvConfig):
        with open(json_path, "r", encoding="utf-8") as f:
            j = json.load(f)
        self.cfg  = cfg
        self.meta = j.get("_meta", {})

        # Effective gates: CLI/config overrides win over the fit JSON's
        # own defaults (so users can re-classify without re-running the fit).
        self.min_pulses = int(cfg.min_pulses_override
                              if cfg.min_pulses_override is not None
                              else self.meta.get("min_pulses", 50))
        self.chi2_max = float(cfg.chi2_max_override
                              if cfg.chi2_max_override is not None
                              else self.meta.get("chi2_max", 3.0))

        self._templates: Dict[str, Template] = {}
        self._is_good:   Dict[str, bool]      = {}

        good_taur, good_tauf = [], []
        for name, rec in j.items():
            if name.startswith("_"):
                continue
            try:
                tr = float(rec["tau_r_ns"]["median"])
                tf = float(rec["tau_f_ns"]["median"])
                chi2 = float(rec["chi2_per_dof"]["median"])
                n_used = int(rec["n_pulses_used"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (np.isfinite(tr) and np.isfinite(tf) and np.isfinite(chi2)):
                continue

            tmpl = Template(name=name, tau_r_ns=tr, tau_f_ns=tf,
                            chi2_med=chi2, n_used=n_used, is_global=False)
            self._templates[name] = tmpl
            good = self._passes_gates(tmpl)
            self._is_good[name] = good
            if good:
                good_taur.append(tr)
                good_tauf.append(tf)

        # Global-median template, computed once.  Only meaningful when at
        # least a handful of channels passed gates.
        self._global: Optional[Template] = None
        if len(good_taur) >= 5:
            self._global = Template(
                name="__global__",
                tau_r_ns=float(np.median(good_taur)),
                tau_f_ns=float(np.median(good_tauf)),
                chi2_med=float("nan"),
                n_used=len(good_taur),
                is_global=True,
            )

    # ------------------------------------------------------------------

    def _passes_gates(self, t: Template) -> bool:
        if t.n_used < self.min_pulses:        return False
        if not np.isfinite(t.chi2_med):       return False
        if t.chi2_med >= self.chi2_max:       return False
        lo_r, hi_r = self.cfg.tau_r_range_ns
        lo_f, hi_f = self.cfg.tau_f_range_ns
        if not (lo_r <= t.tau_r_ns <= hi_r):  return False
        if not (lo_f <= t.tau_f_ns <= hi_f):  return False
        return True

    # ------------------------------------------------------------------

    def lookup(self, channel_name: str) -> Tuple[Optional[Template], str]:
        """Return (template, state) for `channel_name`.

        state ∈ {"good", "fallback_global", "no_template", "bad"}.
        When state is "no_template" or "bad", the caller should skip
        deconv unless cfg.fallback_to_global_template AND a global
        template exists — in which case state will already be
        "fallback_global" instead.
        """
        t = self._templates.get(channel_name)
        if t is not None and self._is_good.get(channel_name, False):
            return t, "good"

        if self.cfg.fallback_to_global_template and self._global is not None:
            return self._global, "fallback_global"

        if t is None:
            return None, "no_template"
        return None, "bad"

    # ------------------------------------------------------------------

    @property
    def n_channels(self) -> int:        return len(self._templates)

    @property
    def n_good(self) -> int:            return sum(self._is_good.values())

    @property
    def has_global(self) -> bool:       return self._global is not None

    @property
    def global_template(self) -> Optional[Template]:
        return self._global

    def summary_str(self) -> str:
        return (f"TemplateDB: {self.n_channels} channels loaded, "
                f"{self.n_good} pass gates "
                f"(min_pulses={self.min_pulses}, chi2_max={self.chi2_max}, "
                f"τ_r∈{self.cfg.tau_r_range_ns}, τ_f∈{self.cfg.tau_f_range_ns}); "
                f"global template: "
                f"{'τ_r=%.2f τ_f=%.2f' % (self._global.tau_r_ns, self._global.tau_f_ns)
                  if self._global else 'unavailable'}")
