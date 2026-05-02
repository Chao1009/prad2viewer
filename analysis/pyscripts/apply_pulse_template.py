#!/usr/bin/env python3
"""
apply_pulse_template.py — drive the C++ NNLS pile-up deconvolution
(`fdec::WaveAnalyzer::Deconvolve`) from Python and visualise the result.

The deconvolver itself lives in `prad2dec/src/WaveAnalyzer.cpp`; this
script just loads templates, walks an EVIO file, calls the binding for
every channel that needs deconv, accumulates per-channel stats, and (when
`--plot-pileup-dir` is given) writes one before/after PNG per piled event.

All knobs come from the `fadc250_waveform.analyzer.nnls_deconv` block in
`daq_config.json` — the same JSON the C++ analyzer reads.

Usage
-----
    python apply_pulse_template.py <evio_path> [<evio_path> …] \
        [--daq-config database/daq_config.json] \
        [--template templates.json] \
        [--out-summary deconv_summary.json] \
        [--plot-pileup-dir plots/pileup] \
        [--max-events N] \
        [--max-plots N] \
        [--enable | --no-enable] \
        [--fallback-global] [--all-peaks]

`--template` overrides the `template_file` field in the config; the path
is otherwise resolved against PRAD2_DATABASE_DIR.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import _common as C
from _common import dec
from _pulse_template_db import (DeconvConfig, TemplateDB, Template,
                                load_nnls_deconv_config)


# ---------------------------------------------------------------------------
# Per-channel accumulator
# ---------------------------------------------------------------------------

@dataclass
class ChanSummary:
    name: str
    channel_id: str
    template_state: str = "no_template"

    n_with_peaks:             int = 0
    n_with_piled:             int = 0
    n_deconv_applied:         int = 0
    n_deconv_singular:        int = 0
    n_deconv_no_template:     int = 0
    n_deconv_bad_template:    int = 0
    n_deconv_fallback_global: int = 0

    height_ratio: List[float] = field(default_factory=list)
    chi2_per_dof: List[float] = field(default_factory=list)

    def finalize(self) -> Dict:
        def _stats(xs: List[float]) -> Dict[str, float]:
            if not xs:
                return {"median": float("nan"), "mean": float("nan"), "n": 0}
            arr = np.asarray(xs, dtype=np.float64)
            return {"median": float(np.median(arr)),
                    "mean":   float(np.mean(arr)),
                    "n":      int(arr.size)}
        return {
            "channel_id":     self.channel_id,
            "template_state": self.template_state,
            "n_with_peaks":   self.n_with_peaks,
            "n_with_piled":   self.n_with_piled,
            "deconv": {
                "applied":              self.n_deconv_applied,
                "singular":             self.n_deconv_singular,
                "skipped_no_template":  self.n_deconv_no_template,
                "skipped_bad_template": self.n_deconv_bad_template,
                "fallback_global":      self.n_deconv_fallback_global,
            },
            "height_ratio_dec_over_wa": _stats(self.height_ratio),
            "chi2_per_dof_global":      _stats(self.chi2_per_dof),
        }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _import_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None


# Continuous-time two-tau template (matches the C++ side).  Used only for
# drawing the dense overlay on the diagnostic PNGs.
def _two_tau_dense(t: np.ndarray, t0: float, tau_r: float, tau_f: float
                   ) -> np.ndarray:
    out = np.zeros_like(t, dtype=np.float64)
    mask = t > t0
    if mask.any():
        dt = t[mask] - t0
        out[mask] = (1.0 - np.exp(-dt / tau_r)) * np.exp(-dt / tau_f)
    return out


def _template_offset_ns(tau_r: float, tau_f: float) -> float:
    return tau_r * np.log((tau_r + tau_f) / tau_r)


def plot_pileup_event(plt, samples_pedsub: np.ndarray,
                      ped_rms: float,
                      peaks: List, deconv,
                      tmpl: Template, clk_ns: float,
                      title: str, out_png: Path) -> None:
    """Two-panel diagnostic: top = waveform + deconv overlay, bottom =
    residual.  `deconv` may be None (e.g. singular case) — then we draw
    just the WaveAnalyzer view."""
    n = samples_pedsub.shape[0]
    t = np.arange(n) * clk_ns

    fig, axes = plt.subplots(2, 1, figsize=(10, 6),
                              sharex=True,
                              gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    ax  = axes[0]
    axR = axes[1]

    # Waveform
    ax.plot(t, samples_pedsub, color="k", lw=1.0, label="pedsub samples")
    ax.axhline(0, color="0.6", lw=0.7, ls="--")
    ax.fill_between(t, -ped_rms, ped_rms, color="0.85", zorder=0,
                     label=f"±ped.rms ({ped_rms:.1f})")

    # WaveAnalyzer-reported peaks (down-triangles, red if pile-up flag set).
    has_label_wa = False
    for pk in peaks:
        color = "C3" if (pk.quality & dec.Q_PEAK_PILED) else "C0"
        ax.axvline(pk.time, color=color, ls=":", lw=0.8, alpha=0.5)
        ax.plot([pk.time], [pk.height], marker="v", ms=9, color=color,
                mec="k", mew=0.5,
                label=("WaveAnalyzer height" if not has_label_wa else None))
        has_label_wa = True

    # Deconvolved overlay
    if deconv is not None and deconv.state in (dec.Q_DECONV_APPLIED,
                                               dec.Q_DECONV_FALLBACK_GLOBAL):
        offset = _template_offset_ns(tmpl.tau_r_ns, tmpl.tau_f_ns)
        t_dense = np.linspace(0, t[-1], 800)
        amp = list(deconv.amplitude)
        total = np.zeros_like(t_dense)
        for k, pk in enumerate(peaks):
            t0_k = float(pk.time) - offset
            comp = amp[k] * _two_tau_dense(t_dense, t0_k,
                                           tmpl.tau_r_ns, tmpl.tau_f_ns)
            total += comp
            ax.plot(t_dense, comp, color="C2", lw=0.7, alpha=0.45,
                    label=("template component" if k == 0 else None))
        ax.plot(t_dense, total, color="C2", lw=1.6,
                label=f"deconv sum  χ²/dof={deconv.chi2_per_dof:.2f}")
        # Up-triangles for deconv heights at the peak times.
        for k, pk in enumerate(peaks):
            ax.plot([pk.time], [deconv.height[k]], marker="^", ms=9,
                    color="C2", mec="k", mew=0.5,
                    label=("deconv height" if k == 0 else None))

        # Residual panel: samples - deconv sum (sampled at sample times).
        # The "deconv sum at sample times" is built by summing the C++
        # template columns analytically (same formula).
        fit = np.zeros_like(samples_pedsub)
        for k, pk in enumerate(peaks):
            t0_k = float(pk.time) - offset
            fit += amp[k] * _two_tau_dense(t.astype(np.float64), t0_k,
                                           tmpl.tau_r_ns, tmpl.tau_f_ns)
        resid = samples_pedsub - fit
        axR.plot(t, resid, color="C2", lw=0.9)
    else:
        axR.text(0.5, 0.5, "no deconv (state="
                            f"{deconv.state if deconv is not None else '—'})",
                  transform=axR.transAxes, ha="center", va="center",
                  color="C3")

    axR.axhline(0, color="0.6", lw=0.7, ls="--")
    axR.fill_between(t, -ped_rms, ped_rms, color="0.85", zorder=0)

    ax.set_ylabel("pedsub ADC")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    axR.set_xlabel("time (ns)")
    axR.set_ylabel("residual")
    axR.grid(True, alpha=0.3)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Map the C++ Q_DECONV_* state byte to a stat-field name on ChanSummary.
def _stat_field_for(state_byte: int) -> str:
    if state_byte == dec.Q_DECONV_APPLIED:         return "n_deconv_applied"
    if state_byte == dec.Q_DECONV_FALLBACK_GLOBAL: return "n_deconv_fallback_global"
    if state_byte == dec.Q_DECONV_SINGULAR:        return "n_deconv_singular"
    if state_byte == dec.Q_DECONV_BAD_TEMPLATE:    return "n_deconv_bad_template"
    if state_byte == dec.Q_DECONV_NO_TEMPLATE:     return "n_deconv_no_template"
    return "n_deconv_no_template"  # Q_DECONV_NOT_RUN — group with no_template


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Drive C++ NNLS pile-up deconvolution and dump diagnostics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("evio_paths", nargs="+",
                    help="EVIO inputs (glob, dir, or split file).")
    ap.add_argument("--template", default="",
                    help="Pulse-template JSON (output of fit_pulse_template.py). "
                         "Overrides template_file in nnls_deconv config.")
    ap.add_argument("--out-summary", default="",
                    help="Per-channel summary JSON.  Empty = none.")
    ap.add_argument("--plot-pileup-dir", default="",
                    help="Directory for before/after PNGs of piled events. "
                         "Empty = no plots.")
    ap.add_argument("--max-events", type=int, default=0,
                    help="Stop after N raw physics sub-events (0 = all).")
    ap.add_argument("--max-plots", type=int, default=200,
                    help="Cap total pile-up PNGs across the whole run.")
    ap.add_argument("--max-plots-per-channel", type=int, default=4,
                    help="Per-channel plot cap so one bright channel doesn't "
                         "eat the entire --max-plots budget.")
    ap.add_argument("--progress-every", type=int, default=2000,
                    help="Print progress line every N physics sub-events.")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config (default: installed default).")
    ap.add_argument("--hc-map-file", default="",
                    help="HyCal modules map (default: database lookup).")
    ap.add_argument("--enable",    dest="enable_override", action="store_true",
                    default=None,
                    help="Force enabled=true regardless of config.")
    ap.add_argument("--no-enable", dest="enable_override", action="store_false",
                    help="Force enabled=false regardless of config.")
    ap.add_argument("--fallback-global", action="store_true",
                    help="Force fallback_to_global_template=true.")
    ap.add_argument("--all-peaks", action="store_true",
                    help="Force apply_to_all_peaks=true.")
    args = ap.parse_args()

    # ---- config ----------------------------------------------------------
    daq_cfg_path = args.daq_config or C.resolve_db_path("daq_config.json")
    cfg = load_nnls_deconv_config(daq_cfg_path)
    if args.enable_override is not None:
        cfg = replace(cfg, enabled=args.enable_override)
    if args.fallback_global:
        cfg = replace(cfg, fallback_to_global_template=True)
    if args.all_peaks:
        cfg = replace(cfg, apply_to_all_peaks=True)

    if not cfg.enabled:
        print("[setup] nnls_deconv disabled — set "
              "fadc250_waveform.analyzer.nnls_deconv.enabled=true "
              "(or pass --enable). Exiting without doing any work.",
              flush=True)
        return

    template_path = args.template or C.resolve_db_path(cfg.template_file)
    if not template_path:
        raise SystemExit("[ERROR] no template JSON given — pass --template "
                         "or set fadc250_waveform.analyzer.nnls_deconv."
                         "template_file in daq_config.json.")
    db = TemplateDB(template_path, cfg)
    print(f"[setup] {db.summary_str()}", flush=True)

    # ---- EVIO discovery (multi-input, dedup) ----------------------------
    all_files: List[str] = []
    seen = set()
    for inp in args.evio_paths:
        for fp in C.discover_split_files(inp):
            if fp not in seen:
                seen.add(fp)
                all_files.append(fp)
    if not all_files:
        raise SystemExit(f"[ERROR] no EVIO splits found for {args.evio_paths}")
    print(f"[setup] {len(all_files)} EVIO split(s)", flush=True)

    # ---- pipeline (loads daq_cfg into the C++ analyzer's wave_cfg) -------
    p = C.setup_pipeline(
        evio_path=args.evio_paths[0],
        max_events=args.max_events,
        daq_config=daq_cfg_path,
        hc_map_file=args.hc_map_file,
    )
    p.evio_files = all_files

    # NB: this script never binds a PulseTemplateStore to the analyzer,
    # so the auto-deconv path (inside Analyze) stays off — we deliberately
    # use the explicit deconvolve() API below to compute the "after" values
    # alongside the unchanged "before" peaks.  Production code paths
    # (Replay / viewer servers) get auto-deconv via the C++ store binding;
    # see analysis/src/Replay.cpp and src/app_state_init.cpp.

    clk_mhz = float(p.cfg.wave_cfg.clk_mhz)
    clk_ns  = (1000.0 / clk_mhz) if clk_mhz > 0 else 4.0

    # ---- name-cache (ROC tag, slot, channel) → (display name, raw key) --
    name_cache: Dict[Tuple[int, int, int], Tuple[str, str]] = {}

    def _key_for(roc_tag: int, crate: Optional[int], s: int, c: int
                 ) -> Tuple[str, str]:
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

    chans: Dict[str, ChanSummary] = {}
    plt = _import_pyplot() if args.plot_pileup_dir else None
    plot_dir = Path(args.plot_pileup_dir) if args.plot_pileup_dir else None
    n_plots = 0
    plots_per_channel: Dict[str, int] = {}

    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    n_phys = 0
    n_files_open = 0
    t0_wall = time.monotonic()
    progress_every = max(1, int(args.progress_every))
    next_progress = progress_every

    def _emit_progress(file_idx: int, file_total: int, fname: str) -> None:
        elapsed = time.monotonic() - t0_wall
        rate = n_phys / elapsed if elapsed > 0 else 0.0
        print(f"  [progress] file {file_idx}/{file_total} {Path(fname).name}  "
              f"phys={n_phys}  channels={len(chans)}  "
              f"plots={n_plots}/{args.max_plots}  "
              f"rate={rate:.0f} ev/s  elapsed={elapsed:.1f}s",
              flush=True)

    try:
        for fpath in p.evio_files:
            if ch.open_auto(fpath) != dec.Status.success:
                print(f"[WARN] skip (cannot open): {fpath}", flush=True)
                continue
            n_files_open += 1
            print(f"[file {n_files_open}/{len(p.evio_files)}] {fpath}",
                  flush=True)
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
                            for c_ in slot.present_channels():
                                cd = slot.channel(c_)
                                if cd.nsamples <= 0:
                                    continue
                                name, chan_id = _key_for(roc.tag, crate, s, c_)

                                samples = np.asarray(cd.samples,
                                                     dtype=np.uint16)
                                wres = p.wave_ana.analyze_result(samples)
                                peaks = wres.peaks
                                if not peaks:
                                    continue

                                summary = chans.get(name)
                                if summary is None:
                                    _, init_state = db.lookup(name)
                                    summary = ChanSummary(name=name,
                                                          channel_id=chan_id,
                                                          template_state=init_state)
                                    chans[name] = summary
                                summary.n_with_peaks += 1

                                has_pile = any(pk.quality & dec.Q_PEAK_PILED
                                               for pk in peaks)
                                if has_pile:
                                    summary.n_with_piled += 1

                                if not (cfg.apply_to_all_peaks or has_pile):
                                    continue

                                tmpl_obj, state = db.lookup(name)
                                if state == "no_template":
                                    summary.n_deconv_no_template += 1
                                    continue
                                if state == "bad":
                                    summary.n_deconv_bad_template += 1
                                    continue
                                # state ∈ {"good", "fallback_global"}

                                cpp_tmpl = dec.PulseTemplate(
                                    tmpl_obj.tau_r_ns,
                                    tmpl_obj.tau_f_ns,
                                    bool(tmpl_obj.is_global))
                                out = p.wave_ana.deconvolve(samples, wres,
                                                            cpp_tmpl)

                                fld = _stat_field_for(int(out.state))
                                setattr(summary, fld,
                                        getattr(summary, fld) + 1)

                                if out.state in (dec.Q_DECONV_APPLIED,
                                                 dec.Q_DECONV_FALLBACK_GLOBAL):
                                    summary.chi2_per_dof.append(
                                        float(out.chi2_per_dof))
                                    h_dec = list(out.height)
                                    for k, pk in enumerate(peaks):
                                        h_wa = float(pk.height)
                                        if h_wa > 0:
                                            summary.height_ratio.append(
                                                float(h_dec[k]) / h_wa)

                                # Plot if pile-up + budget left
                                if (plt is not None and has_pile
                                        and n_plots < args.max_plots):
                                    pcc = plots_per_channel.get(name, 0)
                                    if pcc < args.max_plots_per_channel:
                                        title = (
                                            f"{name}  ev={fadc_evt.info.event_number}  "
                                            f"state={int(out.state)}  "
                                            f"τ_r={tmpl_obj.tau_r_ns:.2f}ns  "
                                            f"τ_f={tmpl_obj.tau_f_ns:.2f}ns  "
                                            f"K={wres.npeaks}  "
                                            f"src={state}")
                                        out_png = (plot_dir / name /
                                                   f"ev{fadc_evt.info.event_number}.png")
                                        pedsub = (samples.astype(np.float64)
                                                  - wres.ped.mean)
                                        plot_pileup_event(
                                            plt, pedsub, float(wres.ped.rms),
                                            peaks, out, tmpl_obj, clk_ns,
                                            title, out_png)
                                        n_plots += 1
                                        plots_per_channel[name] = pcc + 1

                if n_phys >= next_progress:
                    _emit_progress(n_files_open, len(p.evio_files), fpath)
                    while next_progress <= n_phys:
                        next_progress += progress_every
                if done:
                    break

            ch.close()
            if done:
                break
    except KeyboardInterrupt:
        print("\n[interrupted — writing partial summary]", flush=True)

    elapsed = time.monotonic() - t0_wall
    print(f"[done] {n_phys} phys events  /  {len(chans)} channels touched  "
          f"/  {n_plots} pile-up PNGs  /  {elapsed:.1f}s", flush=True)

    # ---- aggregation + summary I/O --------------------------------------
    totals = {
        "n_phys_events": n_phys,
        "n_evio_splits": len(p.evio_files),
        "n_channels":    len(chans),
        "deconv": {"applied": 0, "singular": 0,
                   "skipped_no_template": 0, "skipped_bad_template": 0,
                   "fallback_global": 0},
        "n_with_piled": 0,
    }
    state_counts = {"good": 0, "fallback_global": 0,
                    "no_template": 0, "bad": 0}
    for su in chans.values():
        totals["n_with_piled"] += su.n_with_piled
        totals["deconv"]["applied"]              += su.n_deconv_applied
        totals["deconv"]["singular"]             += su.n_deconv_singular
        totals["deconv"]["skipped_no_template"]  += su.n_deconv_no_template
        totals["deconv"]["skipped_bad_template"] += su.n_deconv_bad_template
        totals["deconv"]["fallback_global"]      += su.n_deconv_fallback_global
        state_counts[su.template_state] = state_counts.get(su.template_state, 0) + 1

    summary_doc: Dict = {
        "_meta": {
            "inputs":          list(args.evio_paths),
            "template_json":   template_path,
            "daq_config_json": daq_cfg_path,
            "config":          {k: list(v) if isinstance(v, tuple) else v
                                for k, v in cfg.__dict__.items()},
            "template_db_summary": db.summary_str(),
            "clk_ns":          clk_ns,
        },
        "totals":               totals,
        "channel_state_counts": state_counts,
        "per_channel": {name: chans[name].finalize()
                        for name in sorted(chans)},
    }

    d = totals["deconv"]
    print("[summary]")
    print(f"  channels touched: {totals['n_channels']}")
    print(f"  channel template state: " +
          "  ".join(f"{k}={v}" for k, v in state_counts.items()))
    print(f"  events with piled peaks: {totals['n_with_piled']}")
    print(f"  deconv applied={d['applied']}  singular={d['singular']}  "
          f"fallback_global={d['fallback_global']}  "
          f"skipped(no_template)={d['skipped_no_template']}  "
          f"skipped(bad_template)={d['skipped_bad_template']}", flush=True)

    if args.out_summary:
        out_path = Path(args.out_summary)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(summary_doc, f, indent=2, sort_keys=False)
        print(f"[write] {out_path}", flush=True)


if __name__ == "__main__":
    main()
