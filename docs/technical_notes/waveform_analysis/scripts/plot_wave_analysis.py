#!/usr/bin/env python3
"""
plot_wave_analysis.py — visualisations for the prad2dec waveform analyzers.

Calls the actual `fdec::WaveAnalyzer` and `fdec::Fadc250FwAnalyzer`
through the prad2py Python bindings, so every figure and printed
number tracks the C++ implementation exactly — no separate Python
re-implementation that can drift from the production code.

Build the project with the python bindings enabled
(`-DBUILD_PYTHON=ON` for CMake) so `import prad2py` resolves.

Outputs seven PNGs into ../plots/ (relative to this script):
  plots/fig1_overview.png           — full waveform with key markers
  plots/fig2_firmware_analysis.png  — Vnoise / TET / Tcross / Va bracket / NSB / NSA
  plots/fig3_soft_analysis.png      — pedestal / smoothing / peak / integration
  plots/fig4_soft_parameters.png    — pedestal-iteration + int_tail_ratio sensitivity
  plots/fig5_smoothing.png          — smoothing on a low-S/N pulse
  plots/fig6_robustness.png         — median+MAD vs. simple-mean pedestal
                                      seed on a synthetic contaminated trace
  plots/fig7_crowded.png            — three closely-spaced pulses, pile-up
                                      flagging (Q_PEAK_PILED)

Run:
  cd docs/technical_notes/waveform_analysis
  python scripts/plot_wave_analysis.py
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import prad2py.dec as dec     # actual C++ analyzers via pybind11

# ---------------------------------------------------------------------------
# Example waveform (100 samples × 4 ns/sample = 400 ns)
# ---------------------------------------------------------------------------
RAW = np.array([
    146,147,144,147,146,147,143,147,146,146,
    144,145,145,146,144,147,146,146,145,147,
    145,146,143,147,145,147,144,146,147,150,
    637,1367,1393,1239,1135,911,767,685,572,514,
    474,429,393,361,332,310,293,286,268,258,
    246,238,228,224,214,215,211,209,204,202,
    198,198,191,192,187,186,180,181,178,178,
    175,177,173,175,171,173,170,172,167,167,
    166,169,165,167,164,169,166,167,165,171,
    164,168,170,168,165,167,168,169,166,167,
], dtype=np.float64)

CLK_NS = 4.0
N = len(RAW)
TIME = np.arange(N) * CLK_NS

# Smaller signal where the per-sample fluctuation is comparable to the
# pulse height — chosen to show what the smoothing kernel does.  Pulse
# peaks at sample 39 (~24 ADC above baseline), comparable to the ±3 ADC
# zig-zag noise on the baseline.
SMALL = np.array([
    147,145,146,144,145,144,145,144,146,145,
    145,144,145,145,145,142,149,143,146,143,
    145,143,144,143,147,144,145,144,145,144,
    145,144,146,143,145,144,146,146,159,169,
    167,161,161,157,152,150,149,147,149,146,
    149,146,148,145,146,146,146,146,146,146,
    146,146,145,143,145,145,145,145,145,145,
    145,145,145,144,148,144,146,143,146,144,
    144,145,145,144,146,145,145,143,145,143,
    145,144,145,144,147,144,145,143,145,143,
], dtype=np.float64)
N_SMALL = len(SMALL)
TIME_SMALL = np.arange(N_SMALL) * CLK_NS


def triangular_smooth(s, smooth_order):
    """Calls fdec::WaveAnalyzer::smooth via the binding so the displayed
    smoothed buffer is bit-identical to what the production peak finder
    sees."""
    cfg = dec.WaveConfig(); cfg.smooth_order = smooth_order
    return np.asarray(
        dec.WaveAnalyzer(cfg).smooth(np.ascontiguousarray(s, dtype=np.uint16)),
        dtype=np.float64,
    )

# ---------------------------------------------------------------------------
# Fadc250FwAnalyzer — firmware Mode 1/2/3 emulation (binding adapter)
# ---------------------------------------------------------------------------
def firmware_analyze(s, *, TET=10.0, NSB_ns=8, NSA_ns=128, NPED=3,
                     MAXPED=1, NSAT=4, MAX_PULSES=4, CLK_NS=4.0):
    """Run fdec::Fadc250FwAnalyzer through the binding on `s`.  The C++
    works in pedsub coordinates internally; this adapter computes the
    matching raw-ADC pedestal, passes it as `PED`, then converts the
    DaqPeak's pedsub fields back to raw so the existing plotting code
    (which renders on a raw-ADC y-axis) keeps working."""
    n = len(s)
    s_u16 = np.ascontiguousarray(s, dtype=np.uint16)

    # raw-ADC pedestal: MAXPED-filtered mean of raw[:NPED] — matches the
    # online firmware's vnoise computation in raw coordinates.
    nped = min(NPED, n)
    vnoise_initial = float(np.mean(s[:nped]))
    if MAXPED > 0:
        kept = [v for v in s[:nped] if abs(v - vnoise_initial) <= MAXPED]
        vnoise = float(np.mean(kept)) if kept else vnoise_initial
    else:
        vnoise = vnoise_initial

    cfg = dec.Fadc250FwConfig()
    cfg.TET, cfg.NSB, cfg.NSA = TET, NSB_ns, NSA_ns
    cfg.NPED, cfg.MAXPED, cfg.NSAT = NPED, MAXPED, NSAT
    cfg.MAX_PULSES, cfg.CLK_NS = MAX_PULSES, CLK_NS

    _, peaks = dec.Fadc250FwAnalyzer(cfg).analyze_full(s_u16, vnoise)
    if not peaks:
        return None
    pk = peaks[0]   # this demo only inspects the first pulse

    # DaqPeak.{vmin, vpeak, va} are pedsub; convert back to raw for plots.
    Vpeak = pk.vpeak + vnoise
    Va    = pk.va    + vnoise
    # Bracket samples (Vba, Vaa) — read directly from raw at coarse / coarse+1.
    k = pk.coarse + 1
    Vba = float(s[pk.coarse]) if 0 <= pk.coarse < n else vnoise
    Vaa = float(s[k]) if 0 <= k < n else vnoise

    return dict(
        vnoise=vnoise, vnoise_initial=vnoise_initial,
        sp=np.maximum(0.0, s.astype(float) - vnoise),
        nsb_samples=int(NSB_ns / CLK_NS),
        nsa_samples=int(NSA_ns / CLK_NS),
        cross=pk.cross_sample, peak=pk.peak_sample,
        Vpeak=Vpeak, Vmin=vnoise, Va=Va, Vba=Vba, Vaa=Vaa, k=k,
        coarse=pk.coarse, fine=pk.fine,
        time_units=pk.time_units, time_ns=pk.time_ns,
        window_lo=pk.window_lo, window_hi=pk.window_hi,
        integral=pk.integral, quality=pk.quality,
    )


# ---------------------------------------------------------------------------
# WaveAnalyzer — soft analyzer (smoothing + iterative pedestal + local maxima)
# ---------------------------------------------------------------------------
def soft_analyze(s, *, smooth_order=2, peak_nsigma=5.0, min_peak_height=10.0,
                 ped_nsamples=30, ped_flatness=1.0, ped_max_iter=3,
                 int_tail_ratio=0.1, tail_break_n=2, peak_pileup_gap=2,
                 clk_mhz=250.0):
    """Run fdec::WaveAnalyzer through the binding on `s` and pack the
    first peak's fields into a dict that matches the existing plotting
    code's expectations.  ``Pedestal`` and ``Peak`` quality bitmasks
    come straight from the C++ struct (no Python re-implementation)."""
    n = len(s)
    s_u16 = np.ascontiguousarray(s, dtype=np.uint16)

    cfg = dec.WaveConfig()
    cfg.smooth_order    = smooth_order
    cfg.peak_nsigma     = peak_nsigma
    cfg.min_peak_height = min_peak_height
    cfg.ped_nsamples    = ped_nsamples
    cfg.ped_flatness    = ped_flatness
    cfg.ped_max_iter    = ped_max_iter
    cfg.int_tail_ratio  = int_tail_ratio
    cfg.tail_break_n    = tail_break_n
    cfg.peak_pileup_gap = peak_pileup_gap
    cfg.clk_mhz         = clk_mhz

    ana = dec.WaveAnalyzer(cfg)
    ped, peaks = ana.analyze_full(s_u16)
    buf = np.asarray(ana.smooth(s_u16), dtype=np.float64)
    if not peaks:
        return None
    pk = peaks[0]

    ped_height  = float(buf[pk.pos]) - ped.mean
    tail_cut    = ped_height * int_tail_ratio
    thr         = max(threshold * ped.rms, min_threshold)
    t_subsample = pk.time * clk_mhz / 1000.0 - pk.pos

    return dict(
        ped_mean=ped.mean, ped_rms=ped.rms,
        ped_nused=ped.nused, ped_quality=ped.quality, ped_slope=ped.slope,
        threshold=thr, buf=buf,
        pos=pk.pos, height=float(pk.height),
        int_left=pk.left, int_right=pk.right,
        integral=float(pk.integral), time_ns=pk.time,
        tail_cut=tail_cut, t_subsample=t_subsample,
        peak_quality=pk.quality,
    )


# ---------------------------------------------------------------------------
# Run analyzers
# ---------------------------------------------------------------------------
fw = firmware_analyze(RAW)
sa = soft_analyze(RAW)

print("=" * 60)
print("Firmware analyzer  (Fadc250FwAnalyzer)")
print("=" * 60)
print(f"  Vnoise (NPED+MAXPED): {fw['vnoise']:.2f} ADC")
print(f"  TET window           : NSB=8 ns -> {fw['nsb_samples']} samples,"
      f" NSA=128 ns -> {fw['nsa_samples']} samples")
print(f"  Tcross               : sample {fw['cross']} (t={fw['cross']*CLK_NS:.1f} ns)")
print(f"  Vpeak (raw)          : {fw['Vpeak']:.0f} ADC at sample {fw['peak']}")
print(f"  Va (mid)             : {fw['Va']:.2f} ADC")
print(f"  bracket              : Vba={fw['Vba']:.0f}@s{fw['k']-1}, Vaa={fw['Vaa']:.0f}@s{fw['k']}")
print(f"  coarse / fine        : {fw['coarse']} / {fw['fine']}  (LSB = 62.5 ps)")
print(f"  time_ns              : {fw['time_ns']:.4f} ns")
print(f"  integration window   : [{fw['window_lo']}, {fw['window_hi']}]"
      f" = {(fw['window_hi']-fw['window_lo']+1)*CLK_NS:.0f} ns")
print(f"  Mode-2 integral      : {fw['integral']:.0f} (pedsub ADC·sample)")
print(f"  quality bitmask      : 0x{fw['quality']:02x}")

print()
print("=" * 60)
print("Soft analyzer  (WaveAnalyzer)")
print("=" * 60)
print(f"  pedestal             : mean={sa['ped_mean']:.2f}, rms={sa['ped_rms']:.2f},"
      f" nused={sa['ped_nused']}/30 (median/MAD bootstrap)")
print(f"  threshold            : {sa['threshold']:.2f} (= max(5*rms, 3))")
print(f"  raw peak             : sample {sa['pos']} (h={sa['height']:.0f})")
print(f"  integration          : [{sa['int_left']}, {sa['int_right']}]"
      f" inclusive (tail-cut={sa['tail_cut']:.0f}, N_break=2)")
print(f"  integral             : {sa['integral']:.0f}")
print(f"  time_ns              : {sa['time_ns']:.3f} ns"
      f" (sub-sample δ = {sa['t_subsample']:+.3f})")


# ---------------------------------------------------------------------------
# Plot 1 — overview
# ---------------------------------------------------------------------------
out_dir = Path(__file__).parent.parent / 'plots'
out_dir.mkdir(exist_ok=True)
fig, ax = plt.subplots(figsize=(10, 4.2))
ax.plot(TIME, RAW, color='#1f77b4', lw=1.0, label='raw samples')
ax.plot(TIME, RAW, 'o', color='#1f77b4', ms=2.5)
ax.axhline(fw['vnoise'], color='#888', lw=0.8, ls='--',
           label=f"Vnoise = {fw['vnoise']:.1f}")
ax.axhline(fw['vnoise'] + 10, color='#d62728', lw=0.8, ls=':',
           label='Vnoise + TET')
ax.axvspan(0, fw['cross'] * CLK_NS, color='#bbb', alpha=0.18,
           label='pedestal region')
ax.axvspan(fw['cross'] * CLK_NS, (fw['peak'] + 1) * CLK_NS,
           color='#d62728', alpha=0.12, label='leading edge')
ax.axvline(fw['cross'] * CLK_NS, color='#d62728', lw=0.8, ls='--')
ax.axvline(fw['peak'] * CLK_NS, color='#2ca02c', lw=0.8, ls='--')
ax.annotate(f"peak\n{int(fw['Vpeak'])} ADC",
            xy=(fw['peak'] * CLK_NS, fw['Vpeak']),
            xytext=(fw['peak'] * CLK_NS + 30, fw['Vpeak'] - 50),
            fontsize=9, ha='left',
            arrowprops=dict(arrowstyle='->', lw=0.7))
ax.set_xlabel('time (ns)')
ax.set_ylabel('ADC counts')
ax.set_title('Example waveform — 100 samples × 4 ns')
ax.legend(loc='upper right', fontsize=8, ncol=2, framealpha=0.9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(out_dir / 'fig1_overview.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 2 — firmware analysis (zoomed near pulse)
# ---------------------------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.5),
                                gridspec_kw={'width_ratios': [1.1, 1]})

# Panel A: rising edge with TDC bracket
zoom_lo, zoom_hi = (fw['cross'] - 4) * CLK_NS, (fw['peak'] + 4) * CLK_NS
xs = np.arange(max(0, fw['cross'] - 4), min(N, fw['peak'] + 5))
axA.plot(xs * CLK_NS, RAW[xs], 'o-', color='#1f77b4', ms=4, lw=1.0,
         label='raw samples')

# TET line + Vnoise
axA.axhline(fw['vnoise'], color='#888', lw=0.8, ls='--', label='Vnoise')
axA.axhline(fw['vnoise'] + 10, color='#d62728', lw=0.8, ls=':',
            label='Vnoise + TET')

# Va horizontal
axA.axhline(fw['Va'], color='#9467bd', lw=1.0, ls='-.',
            label=f"Va = ½(Vpeak+Vmin) = {fw['Va']:.0f}")

# Vba / Vaa markers
axA.plot([(fw['k'] - 1) * CLK_NS], [fw['Vba']], 'D', color='#9467bd',
         ms=8, mfc='white', mew=1.5, label=f"Vba @ s{fw['k']-1}")
axA.plot([fw['k'] * CLK_NS],       [fw['Vaa']], 's', color='#9467bd',
         ms=8, mfc='white', mew=1.5, label=f"Vaa @ s{fw['k']}")

# Tcross vertical
axA.axvline(fw['cross'] * CLK_NS, color='#d62728', lw=0.8, ls='--',
            label=f"Tcross = s{fw['cross']}")
# Coarse-time tick (sample of Vba) and fine arrow at Va
axA.annotate(
    '', xy=(fw['time_ns'], fw['Va']),
    xytext=((fw['k'] - 1) * CLK_NS, fw['Va']),
    arrowprops=dict(arrowstyle='->', lw=1.2, color='#9467bd'))
axA.text(fw['time_ns'] + 1, fw['Va'] + 60,
         f"time = coarse·4 + fine·(4/64)\n"
         f"     = {fw['coarse']}·4 + {fw['fine']}·(1/16)\n"
         f"     = {fw['time_ns']:.3f} ns",
         fontsize=8, color='#222',
         bbox=dict(facecolor='white', edgecolor='#bbb', alpha=0.92, pad=4))

axA.set_xlim(zoom_lo, zoom_hi)
axA.set_xlabel('time (ns)')
axA.set_ylabel('ADC counts')
axA.set_title('Firmware TDC: leading-edge bracket → coarse + fine')
axA.legend(loc='upper left', fontsize=7, framealpha=0.9)
axA.grid(alpha=0.3)

# Panel B: full window with NSB/NSA + integration
axB.plot(TIME, RAW, color='#1f77b4', lw=1.0)
axB.fill_between(TIME[fw['window_lo']:fw['window_hi'] + 1],
                 fw['vnoise'],
                 RAW[fw['window_lo']:fw['window_hi'] + 1],
                 color='#1f77b4', alpha=0.18,
                 label=f"Σ s′ = {fw['integral']:.0f}")
axB.axhline(fw['vnoise'], color='#888', lw=0.8, ls='--')
axB.axvline(fw['cross'] * CLK_NS, color='#d62728', lw=0.8, ls='--',
            label='Tcross')

# NSB / NSA brackets just below the pulse base
yb = fw['vnoise'] + 0.05 * (fw['Vpeak'] - fw['vnoise'])
nsb_x0 = fw['window_lo'] * CLK_NS
nsa_x1 = fw['window_hi'] * CLK_NS
tcross_ns = fw['cross'] * CLK_NS
axB.plot([nsb_x0, tcross_ns], [yb, yb], color='#ffa94d', lw=2)
axB.plot([tcross_ns, nsa_x1], [yb, yb], color='#22b8cf', lw=2)
axB.text((nsb_x0 + tcross_ns) / 2, yb + 25,
         f"NSB = {fw['nsb_samples']*CLK_NS:.0f} ns",
         color='#ffa94d', fontsize=9, ha='center')
axB.text((tcross_ns + nsa_x1) / 2, yb + 25,
         f"NSA = {fw['nsa_samples']*CLK_NS:.0f} ns",
         color='#22b8cf', fontsize=9, ha='center')
axB.set_xlim(0, N * CLK_NS)
axB.set_xlabel('time (ns)')
axB.set_ylabel('ADC counts')
axB.set_title(
    f"Mode-2 integration window  [cross−NSB, cross+NSA]"
    f"\n{fw['nsb_samples']*CLK_NS:.0f} + {fw['nsa_samples']*CLK_NS:.0f}"
    f" = {(fw['window_hi']-fw['window_lo']+1)*CLK_NS:.0f} ns")
axB.legend(loc='upper right', fontsize=8)
axB.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / 'fig2_firmware_analysis.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 3 — soft analyzer
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 4.2))
ax.plot(TIME, RAW, 'o', color='#1f77b4', ms=2.5, label='raw samples')
ax.plot(TIME, sa['buf'], color='#ff7f0e', lw=1.2,
        label=f"smoothed (smooth_order = 2)")
ax.axhline(sa['ped_mean'], color='#888', lw=0.8, ls='--',
           label=f"pedestal = {sa['ped_mean']:.2f} ± {sa['ped_rms']:.2f}")
ax.axhline(sa['ped_mean'] + sa['threshold'], color='#d62728', lw=0.8,
           ls=':', label=f"threshold = max(5·rms, 3) = {sa['threshold']:.2f}")
ax.axhline(sa['ped_mean'] + sa['tail_cut'], color='#2ca02c', lw=0.8,
           ls=':', label=f"tail cut-off = 10% × peak height")

# integration shading
xs = np.arange(sa['int_left'], sa['int_right'] + 1)
ax.fill_between(xs * CLK_NS, sa['ped_mean'], RAW[xs],
                color='#ff7f0e', alpha=0.18,
                label=f"Σ (raw − ped) = {sa['integral']:.0f}")

# peak marker
ax.plot([sa['pos'] * CLK_NS], [RAW[sa['pos']]], '*', color='#d62728',
        ms=12, label=f"peak (raw) @ sample {sa['pos']}")

ax.set_xlim(60, 280)
ax.set_xlabel('time (ns)')
ax.set_ylabel('ADC counts')
ax.set_title('Soft analyzer — smoothing, pedestal, peak, tail-cutoff integration')
ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(out_dir / 'fig3_soft_analysis.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 4 — soft analyzer parameter sensitivity
# ---------------------------------------------------------------------------
# Panel A — pedestal iterative outlier rejection: replays the algorithm on
# the first ped_nsamples samples, marking which got dropped.
# Panel B — int_tail_ratio sensitivity: integrate from the peak walking
# outward for r ∈ {0.05, 0.10, 0.20}; show resulting bounds and integrals.
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.5))

# --- Panel A ---------------------------------------------------------------
PED_FLATNESS = 1.0
nped = 30
peds = RAW[:nped]
xped = np.arange(nped) * CLK_NS

kept_mask = np.ones(nped, dtype=bool)
# Median + MAD·1.4826 robust seed (matches the C++ analyzer).
mean = float(np.median(peds))
rms  = float(np.median(np.abs(peds - mean))) * 1.4826
init_mean, init_rms = mean, rms
for _ in range(3):
    band = max(rms, PED_FLATNESS)
    new_mask = np.abs(peds - mean) < band
    new_kept = kept_mask & new_mask
    if new_kept.sum() == kept_mask.sum() or new_kept.sum() < 5:
        break
    kept_mask = new_kept
    sc = peds[kept_mask]
    mean = float(np.mean(sc)); rms = float(np.std(sc))

band = max(rms, PED_FLATNESS)

# median+MAD seed (no σ-clip iteration yet) for comparison
axA.axhline(init_mean, color='#bbb', lw=0.8, ls=':',
            label=f"median seed = {init_mean:.2f} (MAD·1.4826 = {init_rms:.2f})")
# converged band
axA.fill_between(xped, mean - band, mean + band, color='#1f77b4',
                 alpha=0.10,
                 label=f"final ±max(rms, ped_flatness) = ±{band:.2f}")
axA.axhline(mean, color='#1f77b4', lw=1.2,
            label=f"converged mean = {mean:.2f} ({rms:.2f} rms)")

axA.plot(xped[kept_mask], peds[kept_mask], 'o', color='#1f77b4', ms=5,
         label=f"kept ({int(kept_mask.sum())})")
if (~kept_mask).any():
    axA.plot(xped[~kept_mask], peds[~kept_mask], 'X', color='#d62728',
             ms=10, mew=1.2,
             label=f"rejected ({int((~kept_mask).sum())})")

axA.set_xlabel('time (ns)')
axA.set_ylabel('ADC counts')
axA.set_title(f"Pedestal — iterative outlier rejection\n"
              f"first ped_nsamples = {nped}, ped_max_iter = 3, "
              f"ped_flatness = {PED_FLATNESS}")
axA.legend(loc='lower right', fontsize=7.5, framealpha=0.92)
axA.grid(alpha=0.3)

# --- Panel B ---------------------------------------------------------------
ratios = [0.05, 0.10, 0.20]
colors = ['#2ca02c', '#ff7f0e', '#d62728']
peak_pos = sa['pos']
ped_mean = sa['ped_mean']
ped_height = float(sa['height'])

axB.plot(TIME, RAW, color='#1f77b4', lw=1.0)
axB.axhline(ped_mean, color='#888', lw=0.8, ls='--', alpha=0.8)

for r, col in zip(ratios, colors):
    cut = ped_height * r
    j = peak_pos
    while j + 1 < N and (RAW[j + 1] - ped_mean) >= cut:
        j += 1
    right = j
    j = peak_pos
    while j - 1 >= 0 and (RAW[j - 1] - ped_mean) >= cut:
        j -= 1
    left = j
    integ = float(np.sum(RAW[left:right + 1] - ped_mean))
    nsamp = right - left + 1
    # Cut-off line, only over the integration span
    axB.plot([left * CLK_NS, right * CLK_NS],
             [ped_mean + cut, ped_mean + cut],
             color=col, lw=1.0, ls=':', alpha=0.9)
    axB.axvline(right * CLK_NS, color=col, lw=0.6, ls='--', alpha=0.4)
    axB.text(right * CLK_NS + 1, ped_mean + cut,
             f" r={r:.2f}\n [{left},{right}]\n Σ={integ:.0f} ({nsamp}s)",
             color=col, fontsize=7.5, va='center')

axB.set_xlim(60, 320)
axB.set_xlabel('time (ns)')
axB.set_ylabel('ADC counts')
axB.set_title("int_tail_ratio sensitivity\n"
              "integration stops when (raw − ped) < r × peak height")
axB.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / 'fig4_soft_parameters.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 5 — smoothing on a small-signal zig-zag waveform
# ---------------------------------------------------------------------------
# Demonstrates the `smooth_order` parameter on a low-S/N signal where the
# per-sample fluctuation is comparable to the pulse height.  Local-maxima
# search on the raw trace would trip on the zig-zag; smoothing collapses
# the noise into a single clean peak.
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.4),
                                gridspec_kw={'width_ratios': [1.2, 1]})

variants = [
    (1, '#888',    1.0, '-',  'raw (smooth_order = 1, no smoothing)'),
    (2, '#ff7f0e', 1.6, '-',  'smooth_order = 2  (default)'),
    (4, '#1f77b4', 1.6, '--', 'smooth_order = 4'),
]

# count local maxima above the noise floor for each smoothing setting
def count_local_max(buf, threshold_above=2.0):
    """count points that are strict local maxima at least `threshold_above`
    ADC above the surrounding mean — a coarse stand-in for the analyzer's
    thresholded peak search."""
    n = len(buf); cnt = 0
    bg = np.median(buf)
    for i in range(1, n - 1):
        if buf[i] > buf[i - 1] and buf[i] >= buf[i + 1] \
           and buf[i] - bg > threshold_above:
            cnt += 1
    return cnt

# Panel A: full waveform overlay
axA.plot(TIME_SMALL, SMALL, 'o', color='#bbb', ms=2.5, label='raw samples')
maxima_counts = []
for res, col, lw, ls, lbl in variants:
    sm = triangular_smooth(SMALL, res)
    axA.plot(TIME_SMALL, sm, ls=ls, color=col, lw=lw, label=lbl)
    maxima_counts.append((res, count_local_max(sm)))

axA.axhline(np.median(SMALL), color='#444', lw=0.7, ls=':',
            label=f"baseline ≈ {np.median(SMALL):.0f} ADC")
axA.set_xlabel('time (ns)')
axA.set_ylabel('ADC counts')
axA.set_title("Smoothing on a low-S/N pulse\n"
              "(peak ≈ 24 ADC above ±3 ADC zig-zag baseline)")
axA.legend(loc='upper right', fontsize=8, framealpha=0.92)
axA.grid(alpha=0.3)

# Panel B: zoom near the pulse + spurious-maxima count
zoom_lo, zoom_hi = 100, 240
axB.plot(TIME_SMALL, SMALL, 'o', color='#bbb', ms=3.5, label='raw')
for res, col, lw, ls, lbl in variants:
    sm = triangular_smooth(SMALL, res)
    axB.plot(TIME_SMALL, sm, ls=ls, color=col, lw=lw, label=f'smooth_order = {res}')
axB.set_xlim(zoom_lo, zoom_hi)
axB.set_xlabel('time (ns)')
axB.set_ylabel('ADC counts')

# Maxima count as table inset
table_lines = ["local maxima > +2 ADC:"] + \
              [f"  smooth_order={r}:  {n}" for r, n in maxima_counts]
axB.text(0.97, 0.40, "\n".join(table_lines), transform=axB.transAxes,
         fontsize=8.5, ha='right', va='top', family='monospace',
         bbox=dict(facecolor='white', edgecolor='#bbb', alpha=0.95, pad=4))

axB.set_title("Zoom on pulse region")
axB.legend(loc='upper right', fontsize=8, framealpha=0.92)
axB.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / 'fig5_smoothing.png', dpi=130)
plt.close(fig)

print()
print("Smoothing demo (small-signal waveform):")
for r, n in maxima_counts:
    print(f"  smooth_order={r}: {n} local maxima > +2 ADC")


# ---------------------------------------------------------------------------
# Plot 6 — pedestal seed robustness (median+MAD vs simple-mean)
#
# Synthetic 30-sample window: a previous-event scintillation tail
# contaminates samples 0..9 (exponential decay from ~14 ADC over the
# baseline), then clean ±0.4 ADC noise around the true baseline of 146.
# Both seeds feed the same iterative σ-clip; the median+MAD seed lands
# on the true baseline, the simple-mean seed gets dragged high by the
# contamination and locks the σ-clip band around the wrong value.
# ---------------------------------------------------------------------------
# The OLD code path (simple-mean seed) is intentionally kept as a small
# Python σ-clip helper here — the production binding only exposes the
# new median+MAD seed, so we can't recover the old behaviour through it.
# This is the only re-implementation in this file; everything else
# (including the NEW pedestal below) goes through prad2py.
def _sigma_clip_with_seed(samples, mean, rms, max_iter=3, flatness=1.0):
    sc = samples.copy()
    for _ in range(max_iter):
        keep = np.abs(sc - mean) < max(rms, flatness)
        if keep.sum() == sc.size or keep.sum() < 5:
            break
        sc = sc[keep]
        mean = float(np.mean(sc)); rms = float(np.std(sc))
    return mean, rms, len(sc)


np.random.seed(7)
NPED_DEMO = 30
TRUE_BASE = 146.0
N_CONTAM  = 14   # < half the window so the median is still on the baseline
CONTAM_OFFSET = 2.5  # ADC above baseline — close enough that σ-clip locks
                     # onto the wrong group from a simple-mean seed
demo = np.full(NPED_DEMO, TRUE_BASE) + np.random.normal(0, 0.4, NPED_DEMO)
# Synthetic previous-event tail: a flat-ish 2.5-ADC bias on the first 14
# samples, with light additional noise so the contaminated samples don't
# all sit on a single line.
for i in range(N_CONTAM):
    demo[i] += CONTAM_OFFSET + np.random.normal(0, 0.15)

# OLD seed: simple mean / std → Python σ-clip
old_seed_mean = float(np.mean(demo)); old_seed_rms = float(np.std(demo))
old_mean, old_rms, old_nused = _sigma_clip_with_seed(
    demo, old_seed_mean, old_seed_rms)

# NEW: drop the demo (rounded to uint16) into the actual binding and
# read back the pedestal it computed.  ped_nsamples=len(demo) keeps the
# adaptive trailing-window logic from kicking in.
demo_u16 = np.ascontiguousarray(np.round(demo), dtype=np.uint16)
ncfg = dec.WaveConfig(); ncfg.ped_nsamples = NPED_DEMO
new_ped, _ = dec.WaveAnalyzer(ncfg).analyze_full(demo_u16)
new_seed_mean = float(np.median(demo))      # for legend annotation only
new_seed_rms  = float(np.median(np.abs(demo - new_seed_mean))) * 1.4826
new_mean, new_rms, new_nused = new_ped.mean, new_ped.rms, new_ped.nused

fig, ax = plt.subplots(figsize=(11, 4.6))
xs = np.arange(NPED_DEMO)
contam_mask = xs < N_CONTAM

ax.plot(xs[~contam_mask], demo[~contam_mask], 'o', color='#1f77b4', ms=6,
        label='clean baseline samples')
ax.plot(xs[contam_mask], demo[contam_mask], 'X', color='#d62728', ms=10,
        mew=1.5,
        label='contaminated samples (synthetic previous-event tail)')

ax.axhline(TRUE_BASE, color='#444', lw=1.2, ls=':',
           label=f'true baseline = {TRUE_BASE:.1f}')
ax.axhline(old_mean, color='#d62728', lw=1.4, ls='--',
           label=(f'OLD seed (mean ± σ): converged μ = {old_mean:.2f}, '
                  f'σ = {old_rms:.2f}, nused = {old_nused}/30  '
                  f'(bias = {old_mean - TRUE_BASE:+.2f})'))
ax.axhline(new_mean, color='#2ca02c', lw=1.8,
           label=(f'NEW seed (median + MAD): converged μ = {new_mean:.2f}, '
                  f'σ = {new_rms:.2f}, nused = {new_nused}/30  '
                  f'(bias = {new_mean - TRUE_BASE:+.2f})'))

# annotate the seed values for orientation
ax.annotate(f'mean seed = {old_seed_mean:.2f}',
            xy=(NPED_DEMO - 1, old_seed_mean),
            xytext=(NPED_DEMO - 1.5, old_seed_mean + 1.2),
            fontsize=8, color='#d62728', ha='right',
            arrowprops=dict(arrowstyle='->', color='#d62728', lw=0.8))
ax.annotate(f'median seed = {new_seed_mean:.2f}',
            xy=(NPED_DEMO - 1, new_seed_mean),
            xytext=(NPED_DEMO - 1.5, new_seed_mean - 1.2),
            fontsize=8, color='#2ca02c', ha='right',
            arrowprops=dict(arrowstyle='->', color='#2ca02c', lw=0.8))

ax.set_xlabel('sample index')
ax.set_ylabel('ADC counts')
ax.set_title('Pedestal seed robustness — '
             'a contaminated leading window biases the simple-mean seed,\n'
             'while the median + MAD·1.4826 seed lands on the true baseline')
ax.legend(loc='upper right', fontsize=8.5, framealpha=0.94)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(out_dir / 'fig6_robustness.png', dpi=130)
plt.close(fig)

print()
print("Pedestal seed robustness (synthetic contaminated window):")
print(f"  OLD (mean+σ seed):     converged μ={old_mean:.2f}, σ={old_rms:.2f}, "
      f"nused={old_nused}/30  bias={old_mean - TRUE_BASE:+.2f}")
print(f"  NEW (median+MAD seed): converged μ={new_mean:.2f}, σ={new_rms:.2f}, "
      f"nused={new_nused}/30  bias={new_mean - TRUE_BASE:+.2f}")


# ---------------------------------------------------------------------------
# Plot 7 — crowded window: 3 closely-spaced pulses, pile-up flagging
#
# Synthesises three PbWO₄-like pulses 15 samples (60 ns) and 15 samples
# apart in an 80-sample window, then runs the actual binding.  Each
# pulse rises in ~3 samples, decays with τ ≈ 12 samples — so adjacent
# integration windows touch and trigger Q_PEAK_PILED on both peaks of
# each pair.
# ---------------------------------------------------------------------------
def _pbwo_pulse(t0, height, n, rise_tau=2.5, decay_tau=12.0):
    out = np.zeros(n)
    rise_done = height * (1 - np.exp(-6 / rise_tau))
    for i in range(n):
        if i < t0:
            continue
        if i - t0 < 6:
            out[i] = height * (1 - np.exp(-(i - t0) / rise_tau))
        else:
            out[i] = rise_done * np.exp(-(i - t0 - 6) / decay_tau)
    return out


np.random.seed(13)
N_CROWD = 80
CROWD_PED = 146.0
crowd = (CROWD_PED + np.random.normal(0, 0.4, N_CROWD)
         + _pbwo_pulse(20, 800.0, N_CROWD)
         + _pbwo_pulse(35, 600.0, N_CROWD)
         + _pbwo_pulse(50, 350.0, N_CROWD))
crowd_u16   = np.round(crowd).astype(np.uint16)
TIME_CROWD  = np.arange(N_CROWD) * CLK_NS

crowd_cfg = dec.WaveConfig()
crowd_ana = dec.WaveAnalyzer(crowd_cfg)
crowd_ped, crowd_peaks = crowd_ana.analyze_full(crowd_u16)
crowd_buf = np.asarray(crowd_ana.smooth(crowd_u16), dtype=np.float64)

PEAK_COLORS = ['#d62728', '#2ca02c', '#9467bd', '#8c564b', '#e377c2']
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.plot(TIME_CROWD, crowd, 'o', color='#1f77b4', ms=2.8, label='raw samples')
ax.plot(TIME_CROWD, crowd_buf, color='#ff7f0e', lw=1.0, alpha=0.75,
        label='smoothed (smooth_order = 2)')
ax.axhline(crowd_ped.mean, color='#888', lw=0.8, ls='--',
           label=f'pedestal = {crowd_ped.mean:.1f} ± {crowd_ped.rms:.2f}')

n_piled = 0
for i, pk in enumerate(crowd_peaks):
    col = PEAK_COLORS[i % len(PEAK_COLORS)]
    is_piled = bool(pk.quality & dec.Q_PEAK_PILED)
    n_piled += int(is_piled)
    xs_band = np.arange(pk.left, pk.right + 1)
    ax.fill_between(xs_band * CLK_NS, crowd_ped.mean, crowd[xs_band],
                    color=col, alpha=0.18,
                    label=(f'peak {i}: t={pk.time:6.2f} ns,  h={pk.height:5.0f},'
                           f'  Σ={pk.integral:6.0f},  '
                           f'{"PILED" if is_piled else "isolated"}'))
    ax.plot(pk.pos * CLK_NS, crowd[pk.pos], '*', color=col, ms=14,
            mec='black', mew=0.6)
    ax.annotate(f'{i}', xy=(pk.pos * CLK_NS, crowd[pk.pos]),
                xytext=(pk.pos * CLK_NS, crowd[pk.pos] + 25),
                fontsize=10, ha='center', color=col, fontweight='bold')

ax.set_xlabel('time (ns)')
ax.set_ylabel('ADC counts')
ax.set_title(
    f'Crowded window — three closely-spaced PbWO₄ pulses\n'
    f'analyzer found {len(crowd_peaks)} peaks; '
    f'{n_piled} flagged Q_PEAK_PILED (within '
    f'peak_pileup_gap = {crowd_cfg.peak_pileup_gap} samples)')
ax.legend(loc='upper right', fontsize=8.5, framealpha=0.94)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(out_dir / 'fig7_crowded.png', dpi=130)
plt.close(fig)

print()
print(f"Crowded waveform — {len(crowd_peaks)} peaks found "
      f"({n_piled} piled-up):")
for i, pk in enumerate(crowd_peaks):
    flags = 'PILED' if pk.quality & dec.Q_PEAK_PILED else 'isolated'
    print(f"  peak {i}: pos={pk.pos:3d}, t={pk.time:7.2f} ns, "
          f"h={pk.height:6.1f}, Σ={pk.integral:7.1f}, "
          f"window=[{pk.left:2d}, {pk.right:2d}] ({flags})")
print()
print(f"Wrote 7 PNGs to {out_dir}")
