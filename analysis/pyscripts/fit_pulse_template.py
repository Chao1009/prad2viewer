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

try:
    from scipy.optimize import curve_fit
except ImportError:
    raise SystemExit("[ERROR] scipy is required (pip install scipy)")

import _common as C
from _common import dec  # prad2py.dec re-export


# Number of raw pulses we keep per channel for the diagnostic plots.
# 20 × 50 samples × 8 bytes × ~2000 channels ≈ 16 MB — negligible.
PULSE_CACHE = 20


# ---------------------------------------------------------------------------
# Pulse model
# ---------------------------------------------------------------------------

def two_tau(t: np.ndarray, A: float, t0: float, tau_r: float, tau_f: float) -> np.ndarray:
    """Two-time-constant pulse (rise * fall).  Zero for t < t0."""
    out = np.zeros_like(t, dtype=np.float64)
    mask = t > t0
    if not mask.any():
        return out
    dt = t[mask] - t0
    out[mask] = A * (1.0 - np.exp(-dt / tau_r)) * np.exp(-dt / tau_f)
    return out


def fit_pulse(samples_pedsub: np.ndarray, peak_idx: int, clk_ns: float
              ) -> Optional[Tuple[Dict[str, float], Dict[str, float], float]]:
    """Fit one pedsub pulse.  Time axis is in ns.  Returns (params, perr,
    chi2_per_dof) on success, None on convergence failure or unphysical fit.
    """
    n = samples_pedsub.shape[0]
    if n < 8:
        return None
    t = np.arange(n, dtype=np.float64) * clk_ns

    h_peak = float(samples_pedsub[peak_idx])
    if h_peak <= 0.0:
        return None
    t_peak_ns = peak_idx * clk_ns

    p0 = (2.5 * h_peak, t_peak_ns - 2.0 * clk_ns, 1.0 * clk_ns, 5.0 * clk_ns)
    lower = (0.1 * h_peak, -2.0 * clk_ns, 0.2 * clk_ns, 1.0 * clk_ns)
    upper = (50.0 * h_peak, t[-1],         10.0 * clk_ns, 80.0 * clk_ns)

    tail = samples_pedsub[max(1, int(0.8 * n)) :]
    sigma = max(float(tail.std()), 1.0)

    try:
        popt, pcov = curve_fit(
            two_tau, t, samples_pedsub.astype(np.float64),
            p0=p0, bounds=(lower, upper),
            sigma=np.full(n, sigma), absolute_sigma=True,
            maxfev=2000,
        )
    except (RuntimeError, ValueError):
        return None

    if not (np.isfinite(popt).all() and np.isfinite(pcov).all()):
        return None
    perr = np.sqrt(np.clip(np.diag(pcov), 0.0, None))

    resid = samples_pedsub - two_tau(t, *popt)
    dof = max(1, n - 4)
    chi2_per_dof = float(((resid / sigma) ** 2).sum() / dof)

    A, t0, tau_r, tau_f = popt
    params = {"A": float(A), "t0_ns": float(t0),
              "tau_r_ns": float(tau_r), "tau_f_ns": float(tau_f)}
    errors = {"A": float(perr[0]), "t0_ns": float(perr[1]),
              "tau_r_ns": float(perr[2]), "tau_f_ns": float(perr[3])}
    return params, errors, chi2_per_dof


# ---------------------------------------------------------------------------
# Per-channel accumulator
# ---------------------------------------------------------------------------

@dataclass
class ChannelStats:
    name: str
    channel_id: str            # roc_<tag>_<slot>_<channel>
    n_attempted: int = 0
    n_used: int = 0
    tau_r: List[float] = field(default_factory=list)
    tau_f: List[float] = field(default_factory=list)
    A:     List[float] = field(default_factory=list)
    t0:    List[float] = field(default_factory=list)
    chi2:  List[float] = field(default_factory=list)
    ped_mean_sum: float = 0.0
    ped_rms_sum:  float = 0.0
    ped_n:        int   = 0
    # First few pedsub pulses, kept for the diagnostic plot only.
    sample_pulses: List[np.ndarray] = field(default_factory=list)


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
    a_med,  a_mad  = _median_mad(s.A)
    t0_med, t0_mad = _median_mad(s.t0)
    chi2_med = float(np.median(s.chi2)) if s.chi2 else float("nan")
    chi2_mean = float(np.mean(s.chi2)) if s.chi2 else float("nan")
    good = bool(s.n_used >= min_pulses and chi2_med < chi2_max)
    ped_mean = s.ped_mean_sum / s.ped_n if s.ped_n else float("nan")
    ped_rms  = s.ped_rms_sum  / s.ped_n if s.ped_n else float("nan")
    return {
        "channel_id": s.channel_id,
        "n_pulses_attempted": s.n_attempted,
        "n_pulses_used":      s.n_used,
        "tau_r_ns": {"median": tr_med, "mad": tr_mad},
        "tau_f_ns": {"median": tf_med, "mad": tf_mad},
        "A":        {"median": a_med,  "mad": a_mad},
        "t0_ns":    {"median": t0_med, "mad": t0_mad},
        "chi2_per_dof": {"median": chi2_med, "mean": chi2_mean},
        "ped": {"mean": float(ped_mean), "rms": float(ped_rms)},
        "good_fit": good,
    }


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
    """Per-channel diagnostic: raw pulses + median-parameter fit overlay."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for s in st.sample_pulses:
        t = np.arange(len(s)) * clk_ns
        ax.plot(t, s, color="C0", alpha=0.25, lw=0.7)
    if st.A:
        med_p = (float(np.median(st.A)), float(np.median(st.t0)),
                 float(np.median(st.tau_r)), float(np.median(st.tau_f)))
        t_dense = np.linspace(0, (pre + post) * clk_ns, 400)
        ax.plot(t_dense, two_tau(t_dense, *med_p),
                color="C3", lw=2.0, label="median fit")
        chi2_med = float(np.median(st.chi2)) if st.chi2 else float("nan")
        title = (f"{label + ' ' if label else ''}{st.name}  "
                 f"({st.channel_id})\n"
                 f"n_used={st.n_used}  "
                 f"τ_r={med_p[2]:.2f} ns  τ_f={med_p[3]:.2f} ns  "
                 f"χ²/dof med={chi2_med:.2f}")
    else:
        title = f"{st.name}: no successful fits"
    ax.set_xlabel("time within window (ns)")
    ax.set_ylabel("pedsub ADC")
    ax.set_title(title)
    if st.A:
        ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def plot_summary(plt, summaries: List[Dict], min_pulses: int, chi2_max: float,
                 out_png: Path) -> None:
    """Global τ_r / τ_f / χ² histograms across all fitted channels.
    Channels without enough pulses are excluded from the histograms."""
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

                                slc_pedsub = samples[lo:hi].astype(np.float64) - ped
                                rel_peak = pk.pos - lo

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

                                fit = fit_pulse(slc_pedsub, rel_peak, clk_ns)
                                if fit is None:
                                    continue
                                params, _perr, chi2 = fit
                                st.tau_r.append(params["tau_r_ns"])
                                st.tau_f.append(params["tau_f_ns"])
                                st.A.append(params["A"])
                                st.t0.append(params["t0_ns"])
                                st.chi2.append(chi2)
                                st.n_used += 1
                                n_pulses_used += 1

                                # Cache raw pulses for later plotting.
                                # Only when plotting is on — otherwise skip
                                # to keep memory tight on big runs.
                                if plotting and len(st.sample_pulses) < PULSE_CACHE:
                                    st.sample_pulses.append(slc_pedsub.copy())

                if done:
                    break
                if n_phys and n_phys % 5000 == 0:
                    print(f"  phys={n_phys}  fit_attempted={n_pulses_attempted}"
                          f"  fit_used={n_pulses_used}", flush=True)

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
            "model": "two_tau",
            "model_formula": "A * (1 - exp(-(t-t0)/tau_r)) * exp(-(t-t0)/tau_f), t > t0",
            "clk_ns": clk_ns,
            "pre_samples": pre,
            "post_samples": post,
            "min_pulses": args.min_pulses,
            "chi2_max": args.chi2_max,
            "height_min": args.height_min,
            "height_rms_mult": args.height_rms_mult,
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
