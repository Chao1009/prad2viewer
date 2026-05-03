#!/usr/bin/env python3
"""
deconv_pileup_demo.py — pile-up deconvolution demo plots.

Walks an EVIO file, finds channel-events that the soft WaveAnalyzer
flagged with `Q_PEAK_PILED`, runs the explicit C++ deconvolution
(`fdec::WaveAnalyzer::Deconvolve`) on each one, and writes a small
gallery of before/after PNGs.

The C++ deconv core is the same code production replay uses when
`daq_config.json:nnls_deconv.enabled = true`.  Production currently
ships with `enabled=false` (conservative default while the per-type
templates are validated); this script bypasses the master switch by
calling `WaveAnalyzer.deconvolve()` directly, so it always runs given
a usable template.

Usage
-----
    python deconv_pileup_demo.py <evio_path> \\
        [--daq-config database/daq_config.json] \\
        [--template templates.json] \\
        [--out-dir plots/deconv] \\
        [--max-events N] [--max-plots K]

`--template` overrides `nnls_deconv.template_file` from the daq config.
Both paths resolve against `PRAD2_DATABASE_DIR` if not absolute.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np

try:
    from prad2py import dec
except ImportError as exc:
    raise SystemExit(
        f"[ERROR] cannot import prad2py: {exc}\n"
        "        Build the python bindings (cmake -DBUILD_PYTHON=ON) and "
        "ensure the install directory is on PYTHONPATH."
    )


# --------------------------------------------------------------------------
# Path helpers
# --------------------------------------------------------------------------

def resolve_db_path(p: str) -> str:
    """Resolve a possibly-relative path against PRAD2_DATABASE_DIR."""
    if not p or os.path.isabs(p):
        return p
    db = os.environ.get("PRAD2_DATABASE_DIR")
    return os.path.join(db, p) if db else p


# --------------------------------------------------------------------------
# Plot one channel-event (before / after overlay)
# --------------------------------------------------------------------------

def plot_one(samples: np.ndarray,
             ped: float,
             wres,
             dec_out,
             tmpl,
             title: str,
             out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = samples.size
    t_ns = np.arange(n) * 4.0     # 250 MHz FADC

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_ns, samples.astype(np.float64) - ped,
            color="0.3", lw=1.2, label="raw (pedsub)")
    ax.axhline(0, color="0.85", lw=0.6)

    # WaveAnalyzer peaks (before deconv): vertical lines at the peak time,
    # height marker at peak.height (raw, not deconvolved).
    wa_peaks = list(wres.peaks)
    for pk in wa_peaks:
        ax.plot([pk.time], [pk.height], "v", color="C3", ms=8, mfc="white",
                label="_nolegend_")
        ax.axvline(pk.time, color="C3", ls=":", lw=0.7, alpha=0.7)

    # Deconvolution result: per-peak amplitude + integral, plus the model
    # reconstruction (Σ a_k T(t - τ_k)).
    npk = dec_out.npeaks
    amps = np.array(dec_out.amplitude[:npk], dtype=np.float64)
    heights = np.array(dec_out.height[:npk],  dtype=np.float64)
    times = [wa_peaks[k].time for k in range(npk)]
    for tk, hk in zip(times, heights):
        ax.plot([tk], [hk], "o", color="C0", ms=6, label="_nolegend_")

    # Synthesise the deconvolved model on a dense grid for visual
    # comparison — same continuous-time form as the C++ template.
    t_dense = np.linspace(0, t_ns[-1], 4 * n)
    tau_r = float(tmpl.tau_r_ns)
    tau_f = float(tmpl.tau_f_ns)
    # Analytic peak offset of the unit-amplitude template.
    t_peak_offset = tau_r * np.log((tau_r + tau_f) / tau_r)
    t_max = (1.0 - np.exp(-t_peak_offset / tau_r)) * np.exp(-t_peak_offset / tau_f)
    model = np.zeros_like(t_dense)
    for ak, tk in zip(amps, times):
        # The C++ deconv anchors each pulse so its analytic peak lands on
        # tk (= WaveAnalyzer peak time); shift the template the same way.
        t0 = tk - t_peak_offset
        m = t_dense > t0
        if m.any():
            dt = t_dense[m] - t0
            model[m] += (ak * t_max
                         * (1.0 - np.exp(-dt / tau_r))
                         * np.exp(-dt / tau_f))
    ax.plot(t_dense, model, color="C0", lw=1.4, alpha=0.85,
            label=f"deconv model  τ_r={tau_r:.1f}  τ_f={tau_f:.1f} ns")

    # Legend bookkeeping — proxy artists for the per-peak markers.
    from matplotlib.lines import Line2D
    extra = [
        Line2D([], [], marker="v", color="C3", ls="None",
               mfc="white", ms=8, label=f"WA peaks (n={len(wa_peaks)})"),
        Line2D([], [], marker="o", color="C0", ls="None", ms=6,
               label=f"deconv heights (n={npk})"),
    ]
    ax.legend(handles=ax.lines[:1] + extra + ax.lines[-1:], loc="upper right",
              fontsize=8)

    ax.set_xlabel("time (ns)")
    ax.set_ylabel("ADC − ped")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Find piled-up channel-events in an EVIO file and "
                    "produce before/after deconvolution PNGs.")
    ap.add_argument("evio_path", help="EVIO file (or first split, e.g. .00000)")
    ap.add_argument("--daq-config", default="",
                    help='DAQ config (default "" → installed default).')
    ap.add_argument("--template", default="",
                    help="Override pulse_templates.json path.")
    ap.add_argument("--out-dir", type=Path, default=Path("plots/deconv"),
                    help="Directory for output PNGs (created if missing).")
    ap.add_argument("--max-events", type=int, default=2000,
                    help="Stop after this many physics events (0 = all).")
    ap.add_argument("--max-plots", type=int, default=20,
                    help="Stop after writing this many PNGs.")
    ap.add_argument("--min-peaks", type=int, default=2,
                    help="Only consider channels with at least this many peaks.")
    args = ap.parse_args()

    # ---- DAQ config + WaveAnalyzer ---------------------------------------
    daq_cfg_path = args.daq_config or resolve_db_path("daq_config.json")
    cfg = dec.load_daq_config(daq_cfg_path) if daq_cfg_path else dec.load_daq_config()
    print(f"[setup] DAQ config : {daq_cfg_path or '(default)'}", flush=True)

    wcfg = dec.WaveConfig(cfg.wave_cfg)
    # Honour command-line template override if given.
    tmpl_rel = args.template or cfg.wave_cfg.nnls_deconv.template_file
    tmpl_path = resolve_db_path(tmpl_rel)
    if not tmpl_path or not Path(tmpl_path).is_file():
        print(f"[ERROR] template file not found: {tmpl_path or '(none)'}",
              file=sys.stderr)
        return 1

    store = dec.PulseTemplateStore()
    if not store.load_from_file(tmpl_path, wcfg):
        print(f"[ERROR] PulseTemplateStore.load_from_file({tmpl_path}) failed",
              file=sys.stderr)
        return 1
    print(f"[setup] templates  : {tmpl_path}  "
          f"({store.n_types_loaded} types, "
          f"{store.n_channels_known} channels typed)",
          flush=True)

    wave_ana = dec.WaveAnalyzer(wcfg)

    # ---- Crate map (roc tag → crate index) for naming output files ------
    roc_to_crate: dict[int, int] = {}
    for entry in cfg.roc_tags:
        if entry.type in ("roc", "gem") and entry.crate >= 0:
            roc_to_crate[entry.tag] = entry.crate

    # ---- EVIO loop -------------------------------------------------------
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ch = dec.EvChannel()
    ch.set_config(cfg)
    if ch.open_auto(args.evio_path) != dec.Status.success:
        print(f"[ERROR] cannot open {args.evio_path}", file=sys.stderr)
        return 1

    n_events = n_piled = n_plotted = 0
    while n_plotted < args.max_plots:
        if args.max_events and n_events >= args.max_events:
            break
        if ch.read() != dec.Status.success:
            break
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for ei in range(ch.get_n_events()):
            ch.select_event(ei)
            n_events += 1
            fadc = ch.fadc()    # FADC composite event (HyCal + Veto + LMS)
            for ri in range(fadc.nrocs):
                roc = fadc.roc(ri)
                if not roc.present:
                    continue
                crate = roc_to_crate.get(roc.tag, roc.tag)
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        cd = slot.channel(c)
                        if cd.nsamples <= 0:
                            continue
                        samples = np.array(cd.samples[:cd.nsamples],
                                           dtype=np.uint16)
                        wres = wave_ana.analyze_result(samples)
                        peaks = list(wres.peaks)
                        if len(peaks) < args.min_peaks:
                            continue
                        # Pile-up = at least one peak with Q_PEAK_PILED set.
                        if not any(pk.quality & dec.Q_PEAK_PILED
                                   for pk in peaks):
                            continue
                        n_piled += 1
                        tmpl = store.lookup(roc.tag, s, c)
                        if tmpl is None:
                            continue
                        dec_out = wave_ana.deconvolve(samples, wres, tmpl)
                        if dec_out.state != dec.Q_DECONV_OK:
                            continue
                        # Plot.
                        title = (f"crate{crate} slot{s} ch{c}  "
                                 f"event #{n_events}  "
                                 f"({len(peaks)} peaks, deconv n={dec_out.npeaks})")
                        out_path = (args.out_dir
                                    / f"deconv_c{crate}_s{s:02d}_ch{c:02d}"
                                      f"_ev{n_events:06d}.png")
                        plot_one(samples, wres.ped.mean,
                                 wres, dec_out, tmpl, title, out_path)
                        n_plotted += 1
                        print(f"[plot {n_plotted:3d}/{args.max_plots}] "
                              f"{out_path.name}", flush=True)
                        if n_plotted >= args.max_plots:
                            break
                    if n_plotted >= args.max_plots:
                        break
                if n_plotted >= args.max_plots:
                    break
            if n_plotted >= args.max_plots:
                break

    print(f"[done] events={n_events}  piled-channels-seen={n_piled}  "
          f"plots={n_plotted}  out={args.out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
