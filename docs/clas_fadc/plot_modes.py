"""
Plot the FADC250 mode 1/2/3 algorithms on representative waveforms,
in the same style as the figures in the FADC250 User's Manual.

For each mode we annotate:
  - TET (threshold)
  - NSB and NSA windows around the threshold crossing
  - Vp (peak)
  - T (algorithm-determined time)
And, where relevant:
  - Vmin / Va (mid-amplitude)        for Mode 3
  - integration window + sum         for Mode 2

Saves four PNGs under ./plots/.
"""

import math
import os
import random

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from fadc250_modes import FADC250Config, FADC250Analyzer


# ---------------------------------------------------------------------------
# Synthetic waveform helpers
# ---------------------------------------------------------------------------
def make_pulse(t_ns_truth, amplitude, n=100, clk_ns=4.0,
               tau_rise=2.0, tau_fall=15.0, pedestal=100.0, noise=0.5,
               seed=None):
    """PMT-like pulse: pedestal + Σ A·(1−e^(−Δt/τr))·e^(−Δt/τf)."""
    if seed is not None:
        random.seed(seed)
    out = []
    for i in range(n):
        t = i * clk_ns
        v = pedestal
        for t0, A in zip(t_ns_truth, amplitude):
            dt = t - t0
            if dt > 0:
                v += A * (1.0 - math.exp(-dt / tau_rise)) * math.exp(-dt / tau_fall)
        if noise > 0:
            v += random.gauss(0.0, noise)
        out.append(v)
    return out


# ---------------------------------------------------------------------------
# Plot styling
# ---------------------------------------------------------------------------
PED_LINE_KW   = dict(color="0.55", lw=1.0, ls="--")
TET_LINE_KW   = dict(color="tab:red", lw=1.4, ls="--")
PEAK_KW       = dict(marker="o", ms=8, mfc="white", mec="tab:red", mew=2.0,
                     ls="none", zorder=5)
TIME_KW       = dict(color="tab:green", lw=1.8)
SAMPLE_KW     = dict(color="tab:blue", lw=1.2, marker=".", ms=4)
RAW_FILL_KW   = dict(color="tab:orange", alpha=0.18)
INT_FILL_KW   = dict(color="tab:purple",  alpha=0.22)
NSB_BRACE_KW  = dict(color="tab:orange", lw=1.5)
NSA_BRACE_KW  = dict(color="tab:cyan",   lw=1.5)


def setup_axes(ax, title, t_ns, samples_pedsub, ped_in_plot, tet_pedsub):
    """Common axis setup. Plots samples in pedestal-subtracted space."""
    ax.plot(t_ns, samples_pedsub, **SAMPLE_KW, label="ADC samples (ped-subtracted)")
    ax.axhline(0.0, **PED_LINE_KW, label="Vnoise (baseline)")
    ax.axhline(tet_pedsub, **TET_LINE_KW, label=f"TET = {tet_pedsub:.0f}")
    ax.set_xlabel("time  [ns]")
    ax.set_ylabel("ADC counts (pedestal subtracted)")
    ax.set_title(title)
    ax.grid(alpha=0.25)


def annotate_NSB_NSA(ax, t_cross_ns, nsb_ns, nsa_ns, y, ymin, ymax):
    """Draw NSB / NSA brackets below the waveform with arrows + labels."""
    # NSB bracket (Tcross - NSB·dt  ...  Tcross)
    y_brace = ymin + 0.05 * (ymax - ymin)
    ax.annotate(
        "", xy=(t_cross_ns - nsb_ns, y_brace),
        xytext=(t_cross_ns, y_brace),
        arrowprops=dict(arrowstyle="<->", color=NSB_BRACE_KW["color"], lw=1.5),
    )
    ax.text(t_cross_ns - nsb_ns / 2.0, y_brace + 0.04 * (ymax - ymin),
            f"NSB = {int(nsb_ns / 4)} samples",
            ha="center", va="bottom", color=NSB_BRACE_KW["color"], fontsize=9)

    # NSA bracket (Tcross  ...  Tcross + NSA·dt)
    ax.annotate(
        "", xy=(t_cross_ns, y_brace),
        xytext=(t_cross_ns + nsa_ns, y_brace),
        arrowprops=dict(arrowstyle="<->", color=NSA_BRACE_KW["color"], lw=1.5),
    )
    ax.text(t_cross_ns + nsa_ns / 2.0, y_brace + 0.04 * (ymax - ymin),
            f"NSA = {int(nsa_ns / 4)} samples",
            ha="center", va="bottom", color=NSA_BRACE_KW["color"], fontsize=9)


# ---------------------------------------------------------------------------
# Build a single waveform + run the analyzer
# ---------------------------------------------------------------------------
def analyze_waveform(t_truth_list, amp_list, *, cfg, n=100, clk=4.0,
                    tau_rise=2.0, tau_fall=15.0, pedestal=100.0, seed=1):
    raw = make_pulse(t_truth_list, amp_list, n=n, clk_ns=clk,
                     tau_rise=tau_rise, tau_fall=tau_fall,
                     pedestal=pedestal, seed=seed)
    ana = FADC250Analyzer(cfg)
    res = ana.analyze(raw)
    pedsub = [max(0.0, s - cfg.PED) for s in raw]
    t_axis = [i * clk for i in range(n)]
    return raw, pedsub, t_axis, res


# ---------------------------------------------------------------------------
# Per-mode plotters
# ---------------------------------------------------------------------------
def plot_mode3(savepath, cfg):
    """Mode 3 / TDC: show Vmin, Vp, Va, Vba/Vaa, fine-time interpolation."""
    raw, pedsub, t_axis, res = analyze_waveform(
        [40.0], [800], cfg=cfg, seed=1
    )
    p = res.pulses_tdc[0]
    m1 = res.pulses_mode1[0]

    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=110)
    ymin, ymax = -50, max(pedsub) * 1.20
    setup_axes(ax, "Mode 3 (TDC) — fine-time interpolation",
               t_axis, pedsub, cfg.PED, cfg.TET)

    # Mark Vmin (= Vnoise) and Va (mid-amplitude)
    ax.axhline(p.Vmin, color="0.35", lw=1.0, ls=":", label=f"Vmin = {p.Vmin:.1f}")
    ax.axhline(p.Va,   color="tab:olive", lw=1.2, ls="--",
               label=f"Va = (Vp − Vmin)/2 = {p.Va:.1f}")

    # Peak marker
    t_peak_ns = (t_axis[1] - t_axis[0])  # ensure float
    # Find peak sample index in pedsub
    i_peak = pedsub.index(p.Vpeak) if p.Vpeak in pedsub else int(round(p.T_ns / 4.0))
    ax.plot([t_axis[i_peak]], [p.Vpeak], **PEAK_KW)
    ax.annotate(f"Vp = {p.Vpeak:.1f}",
                xy=(t_axis[i_peak], p.Vpeak),
                xytext=(t_axis[i_peak] + 8, p.Vpeak),
                color="tab:red", fontsize=10,
                arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.0))

    # Vba and Vaa — bracketing samples on leading edge for fine-time
    coarse = p.coarse_clk
    Vba    = pedsub[coarse]
    Vaa    = pedsub[coarse + 1]
    ax.plot([t_axis[coarse]],     [Vba], "s", ms=8, mfc="white", mec="tab:purple",
            mew=2.0, label=f"Vba = {Vba:.1f}")
    ax.plot([t_axis[coarse + 1]], [Vaa], "s", ms=8, mfc="tab:purple", mec="tab:purple",
            mew=2.0, label=f"Vaa = {Vaa:.1f}")
    # Linear interpolation line between Vba and Vaa
    ax.plot([t_axis[coarse], t_axis[coarse + 1]], [Vba, Vaa],
            color="tab:purple", lw=1.2, ls="-", alpha=0.7)

    # T (algorithm-determined time) — vertical line
    ax.axvline(p.T_ns, **TIME_KW,
               label=f"T = coarse·64 + fine = {p.T_units}  →  {p.T_ns:.3f} ns")
    ax.plot([p.T_ns], [p.Va], "X", ms=11, color="tab:green", mec="black", mew=0.8,
            zorder=6)

    # Tcross marker
    t_cross_ns = p.sample_cross * 4.0
    ax.axvline(t_cross_ns, color="tab:red", lw=0.9, ls=":", alpha=0.6)
    ax.text(t_cross_ns + 1, ymax * 0.92, "Tcross\n(first > TET)",
            color="tab:red", fontsize=8.5, va="top")

    # Annotate coarse / fine separately
    box = (
        f"coarse = {p.coarse_clk}  (×4 ns)\n"
        f"fine    = {p.fine}/64\n"
        f"T_units = {p.T_units}   (LSB = 62.5 ps)\n"
        f"T          = {p.T_ns:.3f} ns"
    )
    ax.text(0.98, 0.98, box, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6"))

    ax.set_xlim(0, 100)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="upper right", fontsize=8.5, ncol=1,
              bbox_to_anchor=(1.0, 0.78))
    fig.tight_layout()
    fig.savefig(savepath)
    plt.close(fig)
    print(f"[mode 3] saved {savepath}")


def plot_mode1(savepath, cfg):
    """Mode 1: pulse-raw window [Tcross-NSB, Tcross+NSA] is highlighted."""
    raw, pedsub, t_axis, res = analyze_waveform(
        [40.0], [800], cfg=cfg, seed=1
    )
    p  = res.pulses_tdc[0]
    m1 = res.pulses_mode1[0]

    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=110)
    ymin, ymax = -50, max(pedsub) * 1.20
    setup_axes(ax, "Mode 1 (Pulse Raw) — samples in [Tcross − NSB, Tcross + NSA]",
               t_axis, pedsub, cfg.PED, cfg.TET)

    # Highlight raw window
    t_lo = m1.window_start * 4.0
    t_hi = m1.window_end   * 4.0
    ax.axvspan(t_lo - 2.0, t_hi + 2.0, **RAW_FILL_KW,
               label=f"Mode 1 window  ({len(m1.samples)} samples reported)")

    # Peak
    i_peak = pedsub.index(p.Vpeak) if p.Vpeak in pedsub else int(round(p.T_ns / 4.0))
    ax.plot([t_axis[i_peak]], [p.Vpeak], **PEAK_KW)
    ax.annotate(f"Vp = {p.Vpeak:.1f}",
                xy=(t_axis[i_peak], p.Vpeak),
                xytext=(t_axis[i_peak] + 8, p.Vpeak),
                color="tab:red", fontsize=10,
                arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.0))

    # T marker
    ax.axvline(p.T_ns, **TIME_KW, label=f"T = {p.T_ns:.3f} ns")

    # Tcross + NSB/NSA brackets
    t_cross_ns = p.sample_cross * 4.0
    ax.axvline(t_cross_ns, color="tab:red", lw=0.9, ls=":", alpha=0.6)
    ax.text(t_cross_ns - 1, ymax * 0.92, "Tcross",
            color="tab:red", fontsize=9, ha="right", va="top")
    annotate_NSB_NSA(ax, t_cross_ns, cfg.NSB * 4, cfg.NSA * 4, p.Vpeak, ymin, ymax)

    # Mark each "reported" sample with a filled marker
    for k, v in enumerate(m1.samples):
        ax.plot([(m1.window_start + k) * 4.0], [v], "o",
                ms=5, mfc="tab:orange", mec="tab:orange")

    ax.set_xlim(0, 100)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(savepath)
    plt.close(fig)
    print(f"[mode 1] saved {savepath}")


def plot_mode2(savepath, cfg):
    """Mode 2: integration window highlighted, integral value in box."""
    raw, pedsub, t_axis, res = analyze_waveform(
        [40.0], [800], cfg=cfg, seed=1
    )
    p  = res.pulses_tdc[0]
    m1 = res.pulses_mode1[0]
    m2 = res.pulses_mode2[0]

    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=110)
    ymin, ymax = -50, max(pedsub) * 1.20
    setup_axes(ax, "Mode 2 (Pulse Integral) — Σ samples in [Tcross − NSB, Tcross + NSA]",
               t_axis, pedsub, cfg.PED, cfg.TET)

    # Fill area under the curve in the integration window
    win_idx = list(range(m1.window_start, m1.window_end + 1))
    win_t   = [t_axis[i]  for i in win_idx]
    win_y   = [pedsub[i] for i in win_idx]
    ax.fill_between(win_t, 0.0, win_y, **INT_FILL_KW,
                    label=f"Σ = {m2.integral:.1f}  ({m2.n_samples} samples)")

    # Peak
    i_peak = pedsub.index(p.Vpeak) if p.Vpeak in pedsub else int(round(p.T_ns / 4.0))
    ax.plot([t_axis[i_peak]], [p.Vpeak], **PEAK_KW)
    ax.annotate(f"Vp = {p.Vpeak:.1f}",
                xy=(t_axis[i_peak], p.Vpeak),
                xytext=(t_axis[i_peak] + 8, p.Vpeak),
                color="tab:red", fontsize=10,
                arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.0))

    # T marker
    ax.axvline(p.T_ns, **TIME_KW, label=f"T = {p.T_ns:.3f} ns")

    # NSB / NSA brackets around Tcross
    t_cross_ns = p.sample_cross * 4.0
    ax.axvline(t_cross_ns, color="tab:red", lw=0.9, ls=":", alpha=0.6)
    ax.text(t_cross_ns - 1, ymax * 0.92, "Tcross",
            color="tab:red", fontsize=9, ha="right", va="top")
    annotate_NSB_NSA(ax, t_cross_ns, cfg.NSB * 4, cfg.NSA * 4, p.Vpeak, ymin, ymax)

    # Stats box
    box = (
        f"NSB = {cfg.NSB} samples = {cfg.NSB*4} ns\n"
        f"NSA = {cfg.NSA} samples = {cfg.NSA*4} ns\n"
        f"TET = {cfg.TET:.0f}\n"
        f"Σ   = {m2.integral:.1f}\n"
        f"T   = {p.T_ns:.3f} ns"
    )
    ax.text(0.98, 0.98, box, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6"))

    ax.set_xlim(0, 100)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="upper right", fontsize=8.5,
              bbox_to_anchor=(1.0, 0.84))
    fig.tight_layout()
    fig.savefig(savepath)
    plt.close(fig)
    print(f"[mode 2] saved {savepath}")


def plot_two_pulses(savepath, cfg):
    """Multi-pulse: shows the manual's '4-pulse-per-window' capability."""
    raw, pedsub, t_axis, res = analyze_waveform(
        [30.0, 230.0], [600, 400], cfg=cfg, seed=2
    )

    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=110)
    ymax = max(max(pedsub) * 1.20, 100)
    ymin = -50
    setup_axes(ax, "Multi-pulse window — Mode 1/2/3 applied per pulse",
               t_axis, pedsub, cfg.PED, cfg.TET)

    palette = ["tab:orange", "tab:cyan", "tab:olive", "tab:brown"]
    for idx, (p, m1, m2) in enumerate(zip(res.pulses_tdc,
                                          res.pulses_mode1,
                                          res.pulses_mode2)):
        col = palette[idx % len(palette)]

        # M2 integration fill
        win_idx = list(range(m1.window_start, m1.window_end + 1))
        win_t   = [t_axis[i]  for i in win_idx]
        win_y   = [pedsub[i] for i in win_idx]
        ax.fill_between(win_t, 0.0, win_y, color=col, alpha=0.22)

        # Peak marker
        # locate peak by index nearest to T
        i_peak = max(win_idx, key=lambda j: pedsub[j])
        ax.plot([t_axis[i_peak]], [pedsub[i_peak]], "o", ms=8,
                mfc="white", mec=col, mew=2.0)
        ax.annotate(f"Vp{idx} = {pedsub[i_peak]:.0f}",
                    xy=(t_axis[i_peak], pedsub[i_peak]),
                    xytext=(t_axis[i_peak] + 6, pedsub[i_peak] + 30),
                    color=col, fontsize=9,
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.0))

        # T marker
        ax.axvline(p.T_ns, color=col, lw=1.6)
        # Place label above the peak so it doesn't collide with the TET line
        ax.text(p.T_ns + 1, pedsub[i_peak] * 1.05 + 30,
                f"T{idx}={p.T_ns:.2f} ns\nΣ={m2.integral:.0f}",
                color=col, fontsize=8.5, va="bottom",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=col, lw=0.8))

    ax.set_xlim(0, 100 * 4)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(savepath)
    plt.close(fig)
    print(f"[multi]  saved {savepath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    os.makedirs("plots", exist_ok=True)
    cfg = FADC250Config(PED=100.0, TET=50.0, NSB=4, NSA=10)

    plot_mode3      ("plots/mode3_tdc.png",     cfg)
    plot_mode1      ("plots/mode1_pulse_raw.png", cfg)
    plot_mode2      ("plots/mode2_integral.png",  cfg)
    plot_two_pulses ("plots/multi_pulse.png",     cfg)

    print("\nAll plots written to ./plots/")
