#!/usr/bin/env python3
"""
fit_pulse_template.py — per-channel pulse-shape characterisation for the
PRad-II FADC250 readout.

For every channel that carries waveform samples in the input EVIO file(s)
we collect "clean" isolated pulses (single peak, no pile-up flag, no
overflow, peak well inside the buffer) and fit each one with a two-tau
pulse model:

    p(t) = A · (1 − exp(−(t−t0)/τ_r)) · exp(−(t−t0)/τ_f)    for t ≥ t0
           0                                                 otherwise

Per-pulse fits give a population of (τ_r, τ_f, A, t0, χ²/dof) values per
channel.  We summarise them with median + MAD (robust against the noisy
tail of the population) and emit one entry per channel into a JSON file.

The output is intended as the input to the future pile-up deconvolution
step — see docs/technical_notes/waveform_analysis/wave_analysis.md.

Plotting
--------
With `--plot-dir <dir>` set, we:
  * keep a small cache of raw pulses per channel during fitting,
  * auto-pick the N best-fitting and N worst-fitting channels (criterion:
    median χ²/dof; "good" = passed the good_fit gate, "bad" = enough
    pulses but failed the gate),
  * write one diagnostic PNG per picked channel (raw pulses + median-fit
    overlay) under <plot-dir>/{good,bad,explicit}/<name>.png,
  * write a summary PNG with the global τ_r / τ_f / χ² distributions.

Use `--plot-channels W001,W002,...` to force-plot specific channels in
addition to the auto-picked ones.

Usage
-----
    python fit_pulse_template.py <evio_path> [<evio_path> ...] -o out.json
                                 [--max-events N]
                                 [--min-pulses N]
                                 [--chi2-max X]
                                 [--height-min ADC]
                                 [--channels W001,W002,…]
                                 [--plot-dir plots/]
                                 [--n-plot-good N] [--n-plot-bad N]
                                 [--plot-channels W001,W002,…]
                                 [--no-summary-plot]

Each `<evio_path>` accepts the same shapes as the other analysis scripts —
glob, directory, or single split file.  Multiple paths are concatenated.
Run `--help` for the full list.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import _common as C
from _common import dec  # prad2py.dec re-export


# Number of raw pulses we keep per channel for the diagnostic plots.
# 20 × 50 samples × 8 bytes × ~2000 channels ≈ 16 MB — negligible.
PULSE_CACHE = 20


# ---------------------------------------------------------------------------
# Pulse model
# ---------------------------------------------------------------------------
#
# Fitting itself runs in C++ via dec.WaveAnalyzer.fit_pulse_shape() —
# three-parameter Levenberg-Marquardt on the unit-amplitude two-tau model
# (T(t; t0, τ_r, τ_f) / T_max), with the per-pulse peak-height
# normalisation that decouples shape from amplitude.  See
# prad2dec/src/WaveAnalyzer.cpp::FitPulseShape for the implementation.
#
# We keep a small Python `two_tau_unit` here only for the diagnostic
# plot overlay (drawn once per channel from the median fit params) — the
# fitting hot path no longer touches it.

def two_tau_unit(t: np.ndarray, t0: float, tau_r: float, tau_f: float) -> np.ndarray:
    """Two-tau pulse normalized to unit peak height — for plot overlays.
    The C++ fitter has its own copy."""
    out = np.zeros_like(t, dtype=np.float64)
    mask = t > t0
    if not mask.any():
        return out
    dt = t[mask] - t0
    raw = (1.0 - np.exp(-dt / tau_r)) * np.exp(-dt / tau_f)
    u = tau_r / (tau_r + tau_f)
    t_max = (1.0 - u) * (u ** (tau_r / tau_f))
    out[mask] = raw / t_max
    return out


def two_tau_p_unit(t: np.ndarray, t0: float, tau_r: float, tau_f: float,
                   p: float) -> np.ndarray:
    """[1 - exp(-(t-t0)/τ_r)]^p · exp(-(t-t0)/τ_f), peak normalised to 1.
    For plot overlays only; C++ FitPulseShapeTwoTauP is the real fitter."""
    out = np.zeros_like(t, dtype=np.float64)
    mask = t > t0
    if not mask.any():
        return out
    dt = t[mask] - t0
    u = 1.0 - np.exp(-dt / tau_r)
    u = np.clip(u, 0.0, None)
    raw = (u ** p) * np.exp(-dt / tau_f)
    denom = tau_r + p * tau_f
    u_pk = p * tau_f / denom
    t_max = (u_pk ** p) * ((tau_r / denom) ** (tau_r / tau_f))
    out[mask] = raw / t_max
    return out


# ---------------------------------------------------------------------------
# Per-channel accumulator
# ---------------------------------------------------------------------------

@dataclass
class ChannelStats:
    name: str
    channel_id: str            # roc_<tag>_<slot>_<channel>
    n_attempted: int = 0
    n_used: int = 0
    tau_r:    List[float] = field(default_factory=list)
    tau_f:    List[float] = field(default_factory=list)
    t0:       List[float] = field(default_factory=list)
    p_list:   List[float] = field(default_factory=list)   # populated only when --model=two_tau_p
    chi2:     List[float] = field(default_factory=list)
    peak_amp: List[float] = field(default_factory=list)   # original ADC peak height per pulse
    ped_mean_sum: float = 0.0
    ped_rms_sum:  float = 0.0
    ped_n:        int   = 0
    # First few pedsub pulses, kept for the diagnostic plot only.  Stored
    # *un-normalised* — plotting code divides by the matching peak_amp on
    # the fly so users see the raw shapes alongside the normalised stack.
    sample_pulses: List[np.ndarray] = field(default_factory=list)
    sample_peak_amps: List[float]   = field(default_factory=list)


def _median_mad(values: List[float]) -> Tuple[float, float]:
    if not values:
        return (float("nan"), float("nan"))
    arr = np.asarray(values, dtype=np.float64)
    med = float(np.median(arr))
    mad = float(np.median(np.abs(arr - med)))
    return med, mad


def finalize_channel(s: ChannelStats, min_pulses: int, chi2_max: float
                     ) -> Dict:
    tr_med, tr_mad = _median_mad(s.tau_r)
    tf_med, tf_mad = _median_mad(s.tau_f)
    t0_med, t0_mad = _median_mad(s.t0)
    pk_med, pk_mad = _median_mad(s.peak_amp)
    chi2_med = float(np.median(s.chi2)) if s.chi2 else float("nan")
    chi2_mean = float(np.mean(s.chi2)) if s.chi2 else float("nan")
    good = bool(s.n_used >= min_pulses and chi2_med < chi2_max)
    ped_mean = s.ped_mean_sum / s.ped_n if s.ped_n else float("nan")
    ped_rms  = s.ped_rms_sum  / s.ped_n if s.ped_n else float("nan")
    rec = {
        "channel_id": s.channel_id,
        "n_pulses_attempted": s.n_attempted,
        "n_pulses_used":      s.n_used,
        "tau_r_ns": {"median": tr_med, "mad": tr_mad},
        "tau_f_ns": {"median": tf_med, "mad": tf_mad},
        "t0_ns":    {"median": t0_med, "mad": t0_mad},
        "peak_amp_adc": {"median": pk_med, "mad": pk_mad},
        "chi2_per_dof": {"median": chi2_med, "mean": chi2_mean},
        "ped": {"mean": float(ped_mean), "rms": float(ped_rms)},
        "good_fit": good,
    }
    if s.p_list:                     # only present for --model=two_tau_p
        p_med, p_mad = _median_mad(s.p_list)
        rec["p"] = {"median": p_med, "mad": p_mad}
    return rec


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _import_pyplot():
    """Headless matplotlib import.  None if matplotlib is unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None


def plot_channel(plt, st: ChannelStats, clk_ns: float, pre: int, post: int,
                 out_png: Path, label: str = "") -> None:
    """Two-panel diagnostic — top: raw pulses + raw-scale median fit;
    bottom: normalised stack + unit-amplitude median-shape fit (the view
    the LM actually sees)."""

    fig, (ax_raw, ax_norm) = plt.subplots(
        2, 1, figsize=(8, 8),
        sharex=True,
        gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08))

    has_fit = bool(st.tau_r)

    if has_fit:
        med_t0    = float(np.median(st.t0))
        med_tau_r = float(np.median(st.tau_r))
        med_tau_f = float(np.median(st.tau_f))
        med_pk    = float(np.median(st.peak_amp))
        t_dense   = np.linspace(0, (pre + post) * clk_ns, 400)
        if st.p_list:
            med_p = float(np.median(st.p_list))
            unit_curve = two_tau_p_unit(t_dense, med_t0, med_tau_r,
                                        med_tau_f, med_p)
        else:
            med_p = None
            unit_curve = two_tau_unit(t_dense, med_t0, med_tau_r, med_tau_f)

        # Top — raw scale.
        for s in st.sample_pulses:
            t = np.arange(len(s)) * clk_ns
            ax_raw.plot(t, s, color="C0", alpha=0.25, lw=0.7)
        ax_raw.plot(t_dense, unit_curve * med_pk,
                    color="C3", lw=2.0, label="median shape × median peak")

        # Bottom — normalised stack.
        for s, pk in zip(st.sample_pulses, st.sample_peak_amps):
            if pk <= 0:
                continue
            t = np.arange(len(s)) * clk_ns
            ax_norm.plot(t, s / pk, color="C0", alpha=0.25, lw=0.7)
        ax_norm.plot(t_dense, unit_curve, color="C3", lw=2.0,
                     label="unit-amplitude fit")
        ax_norm.axhline(1.0, color="0.6", lw=0.6, ls="--")

        chi2_med = float(np.median(st.chi2)) if st.chi2 else float("nan")
        p_str = f"  p={med_p:.2f}" if med_p is not None else ""
        title = (f"{label + ' ' if label else ''}{st.name}  "
                 f"({st.channel_id})\n"
                 f"n_used={st.n_used}  "
                 f"τ_r={med_tau_r:.2f} ns  τ_f={med_tau_f:.2f} ns{p_str}  "
                 f"peak={med_pk:.0f} ADC  χ²/dof med={chi2_med:.2f}")
    else:
        for s in st.sample_pulses:
            t = np.arange(len(s)) * clk_ns
            ax_raw.plot(t, s, color="C0", alpha=0.25, lw=0.7)
        title = f"{st.name}: no successful fits"

    ax_raw.set_ylabel("pedsub ADC")
    ax_raw.grid(True, alpha=0.3)
    if has_fit:
        ax_raw.legend(loc="best", fontsize=8)
    ax_raw.set_title(title)

    ax_norm.set_xlabel("time within window (ns)")
    ax_norm.set_ylabel("pedsub / peak  (unit)")
    ax_norm.grid(True, alpha=0.3)
    if has_fit:
        ax_norm.legend(loc="best", fontsize=8)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_summary(plt, summaries: List[Dict], min_pulses: int, chi2_max: float,
                 out_png: Path) -> None:
    """Global τ_r / τ_f / χ² histograms across all fitted channels."""
    have = [s for s in summaries if s["n_pulses_used"] >= min_pulses]
    if not have:
        print("[plot] no channels with enough pulses for summary", flush=True)
        return
    tau_r = np.array([s["tau_r_ns"]["median"] for s in have])
    tau_f = np.array([s["tau_f_ns"]["median"] for s in have])
    chi2  = np.array([s["chi2_per_dof"]["median"] for s in have])
    good  = np.array([s["good_fit"] for s in have])

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    ax.hist(tau_r, bins=60, color="C0", alpha=0.85)
    ax.set_xlabel("τ_r (ns) — channel median")
    ax.set_ylabel("channels")
    ax.set_title(f"rise time   (N={len(have)})")
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.hist(tau_f, bins=60, color="C2", alpha=0.85)
    ax.set_xlabel("τ_f (ns) — channel median")
    ax.set_ylabel("channels")
    ax.set_title("fall time")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    bins = np.logspace(np.log10(max(chi2.min(), 1e-2)),
                       np.log10(max(chi2.max(), 10.0)), 60)
    ax.hist(chi2[good],  bins=bins, color="C2", alpha=0.7, label="good_fit")
    ax.hist(chi2[~good], bins=bins, color="C3", alpha=0.7, label="bad")
    ax.axvline(chi2_max, color="k", ls="--", lw=1, label=f"χ²_max={chi2_max}")
    ax.set_xscale("log")
    ax.set_xlabel("χ²/dof — channel median")
    ax.set_ylabel("channels")
    ax.set_title(f"fit quality   (good={int(good.sum())}, "
                 f"bad={int((~good).sum())})")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    sc = ax.scatter(tau_r, tau_f, c=np.log10(np.clip(chi2, 1e-2, None)),
                    s=8, cmap="viridis")
    fig.colorbar(sc, ax=ax, label="log₁₀(χ²/dof)")
    ax.set_xlabel("τ_r (ns)")
    ax.set_ylabel("τ_f (ns)")
    ax.set_title("τ_r vs τ_f, coloured by χ²/dof")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def select_plot_targets(stats: Dict[str, ChannelStats], min_pulses: int,
                        chi2_max: float, n_good: int, n_bad: int
                        ) -> Tuple[List[str], List[str]]:
    """Auto-pick the N best- and N worst-fitting channels by median χ²/dof.
    "Good" = enough pulses AND median χ² < chi2_max (then sort ascending).
    "Bad" = enough pulses but median χ² >= chi2_max (then sort descending)."""
    good, bad = [], []
    for name, s in stats.items():
        if s.n_used < min_pulses or not s.chi2:
            continue
        med = float(np.median(s.chi2))
        (good if med < chi2_max else bad).append((med, name))
    good.sort()
    bad.sort(reverse=True)
    return [n for _, n in good[:n_good]], [n for _, n in bad[:n_bad]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Per-channel pulse-shape fit on FADC250 waveforms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("evio_paths", nargs="+",
                    help="One or more EVIO inputs (each may be a glob, "
                         "directory, or single split file).  Splits found "
                         "for every input are concatenated.")
    ap.add_argument("-o", "--out", required=True,
                    help="Output JSON path.")
    ap.add_argument("--max-events", type=int, default=0,
                    help="Stop after N raw physics sub-events across all "
                         "files (0 = process everything).")
    ap.add_argument("--max-pulses-per-channel", type=int, default=500,
                    help="Stop fitting a channel once it has this many "
                         "accepted pulses.")
    ap.add_argument("--min-pulses", type=int, default=50,
                    help="Channel needs at least this many converged fits "
                         "to be marked good_fit=true (also gates summary "
                         "plot inclusion).")
    ap.add_argument("--chi2-max", type=float, default=3.0,
                    help="Median chi2/dof threshold for good_fit=true.")
    ap.add_argument("--height-min", type=float, default=20.0,
                    help="Minimum peak height (ADC, pedsub) for fitting.")
    ap.add_argument("--height-rms-mult", type=float, default=10.0,
                    help="Also require peak height >= this × ped.rms.")
    ap.add_argument("--model", choices=["two_tau", "two_tau_p"],
                    default="two_tau",
                    help="Pulse-shape model.  'two_tau': "
                         "(1 - exp(-(t-t0)/τ_r)) · exp(-(t-t0)/τ_f).  "
                         "'two_tau_p': adds a rise-edge exponent so the "
                         "rise *shape* is free, not just the timescale: "
                         "[1 - exp(-(t-t0)/τ_r)]^p · exp(-(t-t0)/τ_f).  "
                         "Use two_tau_p when the standard model's onset "
                         "is too sharp on the data (high-amp PMT pulses "
                         "in particular).")
    ap.add_argument("--model-err-floor", type=float, default=0.01,
                    help="Floor on per-sample sigma in the normalised fit, as "
                         "a fraction of unit peak.  Caps how much a small "
                         "model-vs-data mismatch (sub-sample alignment, slow "
                         "second-tail component) can inflate chi2/dof on "
                         "high-amplitude pulses.  0 disables the floor.")
    ap.add_argument("--pre-samples",  type=int, default=8,
                    help="Samples kept before the peak for the fit window.")
    ap.add_argument("--post-samples", type=int, default=40,
                    help="Samples kept after the peak for the fit window.")
    ap.add_argument("--channels", default="",
                    help="Comma-separated channel name filter "
                         "(e.g. W001,W002).  Empty = all channels.")
    ap.add_argument("--plot-dir", default="",
                    help="Directory for diagnostic PNGs.  Empty = no plots.")
    ap.add_argument("--n-plot-good", type=int, default=5,
                    help="Auto-pick this many best-fitting channels for "
                         "diagnostic plots (sorted ascending χ²).")
    ap.add_argument("--n-plot-bad", type=int, default=5,
                    help="Auto-pick this many worst-fitting channels for "
                         "diagnostic plots (sorted descending χ²).")
    ap.add_argument("--plot-channels", default="",
                    help="Comma-separated extra channels to plot regardless "
                         "of fit quality (e.g. W001,W002).")
    ap.add_argument("--no-summary-plot", action="store_true",
                    help="Skip the global τ_r/τ_f/χ² summary PNG.")
    ap.add_argument("--progress-every", type=int, default=2000,
                    help="Print a progress line every N physics sub-events "
                         "(file index, events, fits, rate, elapsed).")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config (default: installed default).")
    ap.add_argument("--hc-map-file", default="",
                    help="HyCal modules map (default: database lookup).")
    args = ap.parse_args()

    chan_filter   = {s.strip() for s in args.channels.split(",") if s.strip()}
    extra_plot    = {s.strip() for s in args.plot_channels.split(",") if s.strip()}
    plotting      = bool(args.plot_dir)
    plot_dir      = Path(args.plot_dir) if plotting else None

    # Discover EVIO files for every input arg, concatenate, dedupe order-
    # preserving so users see the same iteration order as the CLI list.
    all_files: List[str] = []
    seen = set()
    for inp in args.evio_paths:
        for fp in C.discover_split_files(inp):
            if fp not in seen:
                seen.add(fp)
                all_files.append(fp)
    if not all_files:
        raise SystemExit(f"[ERROR] no EVIO splits found for inputs: "
                         f"{args.evio_paths}")
    print(f"[setup] {len(all_files)} EVIO split(s) across "
          f"{len(args.evio_paths)} input arg(s)", flush=True)

    # setup_pipeline still does its own discovery, but we override the
    # resulting evio_files with the concatenated multi-input list below.
    p = C.setup_pipeline(
        evio_path=args.evio_paths[0],
        max_events=args.max_events,
        daq_config=args.daq_config,
        hc_map_file=args.hc_map_file,
    )
    p.evio_files = all_files

    clk_mhz = float(p.cfg.wave_cfg.clk_mhz)
    clk_ns  = (1000.0 / clk_mhz) if clk_mhz > 0 else 4.0
    pre, post = args.pre_samples, args.post_samples

    stats: Dict[str, ChannelStats] = {}
    name_cache: Dict[Tuple[int, int, int], Tuple[str, str]] = {}

    def _key_for(roc_tag: int, crate: Optional[int], s: int, c: int) -> Tuple[str, str]:
        cached = name_cache.get((roc_tag, s, c))
        if cached is not None:
            return cached
        chan_id = f"{roc_tag}_{s}_{c}"
        name = chan_id
        if crate is not None:
            mod = p.hycal.module_by_daq(crate, s, c)
            if mod is not None:
                name = mod.name
        name_cache[(roc_tag, s, c)] = (name, chan_id)
        return name, chan_id

    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    n_phys = n_pulses_attempted = n_pulses_used = 0
    n_files_open = 0
    t0_wall = time.monotonic()
    progress_every = max(1, int(args.progress_every))
    next_progress = progress_every

    def _emit_progress(file_idx: int, file_total: int, fname: str) -> None:
        elapsed = time.monotonic() - t0_wall
        rate = n_phys / elapsed if elapsed > 0 else 0.0
        print(f"  [progress] file {file_idx}/{file_total} {Path(fname).name}  "
              f"phys={n_phys}  fit_attempted={n_pulses_attempted}  "
              f"fit_used={n_pulses_used}  "
              f"rate={rate:.0f} ev/s  elapsed={elapsed:.1f}s",
              flush=True)

    try:
        for fpath in p.evio_files:
            if ch.open_auto(fpath) != dec.Status.success:
                print(f"[WARN] skip (cannot open): {fpath}", flush=True)
                continue
            n_files_open += 1
            print(f"[file {n_files_open}/{len(p.evio_files)}] {fpath}", flush=True)
            done = False

            while ch.read() == dec.Status.success:
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for i in range(ch.get_n_events()):
                    decoded = ch.decode_event(i, with_ssp=False)
                    if not decoded["ok"]:
                        continue
                    n_phys += 1
                    if args.max_events and n_phys >= args.max_events:
                        done = True

                    fadc_evt = decoded["event"]
                    for ri in range(fadc_evt.nrocs):
                        roc = fadc_evt.roc(ri)
                        if not roc.present:
                            continue
                        crate = p.crate_map.get(roc.tag)
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                cd = slot.channel(c)
                                if cd.nsamples <= 0:
                                    continue
                                name, chan_id = _key_for(roc.tag, crate, s, c)
                                if chan_filter and name not in chan_filter:
                                    continue

                                samples = np.asarray(cd.samples, dtype=np.uint16)
                                ped, rms, peaks = p.wave_ana.analyze(samples)

                                if len(peaks) != 1:
                                    continue
                                pk = peaks[0]
                                if pk.quality != 0:
                                    continue
                                if pk.height < args.height_min:
                                    continue
                                if pk.height < args.height_rms_mult * rms:
                                    continue
                                if pk.overflow:
                                    continue

                                lo, hi = pk.pos - pre, pk.pos + post + 1
                                if lo < 0 or hi > samples.shape[0]:
                                    continue

                                # Hand the raw uint16 slice straight to the
                                # C++ fitter — no float64 conversion, no
                                # pedsub in Python.
                                slice_u16 = samples[lo:hi]
                                rel_peak  = pk.pos - lo

                                st = stats.get(name)
                                if st is None:
                                    st = ChannelStats(name=name, channel_id=chan_id)
                                    stats[name] = st
                                if st.n_used >= args.max_pulses_per_channel:
                                    continue

                                st.n_attempted += 1
                                n_pulses_attempted += 1
                                st.ped_mean_sum += float(ped)
                                st.ped_rms_sum  += float(rms)
                                st.ped_n        += 1

                                if args.model == "two_tau_p":
                                    fit = dec.WaveAnalyzer.fit_pulse_shape_two_tau_p(
                                        slice_u16, rel_peak,
                                        float(ped), float(rms), clk_ns,
                                        args.model_err_floor)
                                    if not fit.ok:
                                        continue
                                    st.tau_r.append(fit.tau_r_ns)
                                    st.tau_f.append(fit.tau_f_ns)
                                    st.t0.append(fit.t0_ns)
                                    st.p_list.append(fit.p)
                                    st.peak_amp.append(fit.peak_amp)
                                    st.chi2.append(fit.chi2_per_dof)
                                else:  # two_tau
                                    fit = dec.WaveAnalyzer.fit_pulse_shape(
                                        slice_u16, rel_peak,
                                        float(ped), float(rms), clk_ns,
                                        args.model_err_floor)
                                    if not fit.ok:
                                        continue
                                    st.tau_r.append(fit.tau_r_ns)
                                    st.tau_f.append(fit.tau_f_ns)
                                    st.t0.append(fit.t0_ns)
                                    st.peak_amp.append(fit.peak_amp)
                                    st.chi2.append(fit.chi2_per_dof)
                                st.n_used += 1
                                n_pulses_used += 1

                                # Cache raw pulses + their peak amps for the
                                # diagnostic plots (normalised stack vs
                                # median fit).  Pedsub on the fly here so
                                # the cache lives in absolute ADC like the
                                # fit's input.  Only when plotting is on.
                                if plotting and len(st.sample_pulses) < PULSE_CACHE:
                                    st.sample_pulses.append(
                                        slice_u16.astype(np.float64) - ped)
                                    st.sample_peak_amps.append(fit.peak_amp)

                if n_phys >= next_progress:
                    _emit_progress(n_files_open, len(p.evio_files), fpath)
                    # Bump past every threshold this CODA read crossed, so
                    # we always have a fresh next-target rather than firing
                    # repeatedly on the same plateau.
                    while next_progress <= n_phys:
                        next_progress += progress_every

                if done:
                    break

            ch.close()
            if done:
                break
    except KeyboardInterrupt:
        print("\n[interrupted — writing partial results]", flush=True)

    elapsed = time.monotonic() - t0_wall
    print(f"[done] {n_phys} phys events  /  {n_pulses_attempted} fits attempted"
          f"  /  {n_pulses_used} converged  /  {len(stats)} channels"
          f"  /  {elapsed:.1f}s", flush=True)

    # Aggregate + write JSON.
    out: Dict = {
        "_meta": {
            "inputs": list(args.evio_paths),
            "n_evio_splits": len(p.evio_files),
            "n_phys_events": n_phys,
            "n_pulses_attempted": n_pulses_attempted,
            "n_pulses_used": n_pulses_used,
            "n_channels": len(stats),
            "model": args.model,
            "model_formula": (
                "((1 - exp(-(t-t0)/tau_r)) * exp(-(t-t0)/tau_f)) / T_max(tau_r, tau_f),  t > t0"
                if args.model == "two_tau" else
                "([1 - exp(-(t-t0)/tau_r)]^p * exp(-(t-t0)/tau_f)) / T_max(tau_r, tau_f, p),  t > t0   "
                "[two-tau with rise-edge exponent]"
            ),
            "fit_params": (["t0_ns", "tau_r_ns", "tau_f_ns"]
                           if args.model == "two_tau"
                           else ["t0_ns", "tau_r_ns", "tau_f_ns", "p"]),
            "amplitude_handling": "each pulse divided by its own pedsub peak height before the fit; peak_amp_adc records the un-normalised peak amplitude per channel for downstream diagnostics",
            "sigma_per_sample": "ped.rms / pulse_peak_amp (i.e. relative noise on the normalised pulse) — keeps chi2/dof comparable across channels with very different signal amplitudes",
            "clk_ns": clk_ns,
            "pre_samples": pre,
            "post_samples": post,
            "min_pulses": args.min_pulses,
            "chi2_max": args.chi2_max,
            "height_min": args.height_min,
            "height_rms_mult": args.height_rms_mult,
            "model_err_floor": args.model_err_floor,
        }
    }
    summaries: List[Dict] = []
    for name in sorted(stats):
        rec = finalize_channel(stats[name], args.min_pulses, args.chi2_max)
        out[name] = rec
        summaries.append({"name": name, **rec})

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=False)
    print(f"[write] {out_path}", flush=True)

    # ---- plotting ----------------------------------------------------------
    if not plotting:
        return
    plt = _import_pyplot()
    if plt is None:
        print("[plot] matplotlib unavailable — skipping", flush=True)
        return

    # Auto-pick best/worst, then add the user's explicit list.
    good_names, bad_names = select_plot_targets(
        stats, args.min_pulses, args.chi2_max,
        args.n_plot_good, args.n_plot_bad)

    print(f"[plot] dir={plot_dir}  "
          f"good={len(good_names)}  bad={len(bad_names)}  "
          f"explicit={len(extra_plot)}", flush=True)

    for name in good_names:
        st = stats.get(name)
        if st is None or not st.sample_pulses:
            continue
        plot_channel(plt, st, clk_ns, pre, post,
                     plot_dir / "good" / f"{name}.png", label="[good]")
    for name in bad_names:
        st = stats.get(name)
        if st is None or not st.sample_pulses:
            continue
        plot_channel(plt, st, clk_ns, pre, post,
                     plot_dir / "bad" / f"{name}.png", label="[bad]")
    for name in sorted(extra_plot):
        st = stats.get(name)
        if st is None or not st.sample_pulses:
            print(f"[plot] no pulses for explicit channel {name}", flush=True)
            continue
        plot_channel(plt, st, clk_ns, pre, post,
                     plot_dir / "explicit" / f"{name}.png")

    if not args.no_summary_plot:
        plot_summary(plt, summaries, args.min_pulses, args.chi2_max,
                     plot_dir / "summary.png")
        print(f"[plot] wrote {plot_dir / 'summary.png'}", flush=True)


if __name__ == "__main__":
    main()
