#!/usr/bin/env python3
"""
plot_wave_analysis.py — visualisations for the prad2dec waveform analyzers.

Re-implements `fdec::WaveAnalyzer` (soft) and `fdec::Fadc250FwAnalyzer`
(firmware Mode 1/2/3 emulation) in pure Python and runs both on the
example waveform shipped in this directory.

Outputs five PNGs into ../plots/ (relative to this script):
  plots/fig1_overview.png           — full waveform with key markers
  plots/fig2_firmware_analysis.png  — Vnoise / TET / Tcross / Va bracket / NSB / NSA
  plots/fig3_soft_analysis.png      — pedestal / smoothing / peak / integration
  plots/fig4_soft_parameters.png    — pedestal-iteration + int_tail_ratio sensitivity
  plots/fig5_smoothing.png          — smoothing on a low-S/N pulse

Run:
  cd docs/technical_notes/waveform_analysis
  python scripts/plot_wave_analysis.py
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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
    """Same triangular kernel as fdec::WaveAnalyzer::smooth (res = 1 → no-op)."""
    n = len(s)
    if smooth_order <= 1:
        return s.astype(np.float64).copy()
    buf = np.empty(n)
    for i in range(n):
        val = float(s[i]); wsum = 1.0
        for j in range(1, smooth_order):
            if j > i or i + j >= n:
                continue
            w = 1.0 - j / float(smooth_order + 1)
            val += w * (float(s[i - j]) + float(s[i + j]))
            wsum += 2.0 * w
        buf[i] = val / wsum
    return buf

# ---------------------------------------------------------------------------
# Fadc250FwAnalyzer — firmware Mode 1/2/3 emulation (ns convention)
# ---------------------------------------------------------------------------
def firmware_analyze(s, *, TET=10.0, NSB_ns=8, NSA_ns=128, NPED=3,
                     MAXPED=1, NSAT=4, MAX_PULSES=4, CLK_NS=4.0):
    """Returns dict with vnoise + first-pulse fields, mirroring DaqPeak."""
    # --- Vnoise: mean of first NPED samples, optionally with online outlier
    #     rejection if MAXPED>0 (sample dropped if pedsub > MAXPED)
    n = len(s)
    nped = min(NPED, n)
    vnoise_initial = float(np.mean(s[:nped]))
    if MAXPED > 0:
        kept = [v for v in s[:nped] if abs(v - vnoise_initial) <= MAXPED]
        vnoise = float(np.mean(kept)) if kept else vnoise_initial
    else:
        vnoise = vnoise_initial

    # --- pedestal-subtracted samples
    sp = np.maximum(0.0, s - vnoise)

    # --- NSB/NSA from ns to integer samples (floor)
    nsb_s = int(NSB_ns / CLK_NS)
    nsa_s = int(NSA_ns / CLK_NS)

    # --- Walk for the first pulse
    i = NPED
    Vmin = vnoise
    while i < n - 1:
        # find pulse start (first sample > Vnoise)
        while i < n and s[i] <= vnoise:
            i += 1
        if i >= n:
            break

        # walk to peak
        i_peak = i
        while i_peak + 1 < n and s[i_peak + 1] > s[i_peak]:
            i_peak += 1
        Vpeak = float(s[i_peak])

        if (Vpeak - vnoise) <= TET:    # below threshold → keep searching
            i = i_peak + 1
            continue

        # Tcross: first leading-edge sample whose pedsub value > TET
        cross = i
        while cross <= i_peak and (s[cross] - vnoise) <= TET:
            cross += 1
        if cross > i_peak:
            i = i_peak + 1
            continue

        # NSAT consecutive samples above TET
        ok = all((s[k] - vnoise) > TET for k in range(cross, min(cross + NSAT, n)))
        if not ok:
            i = i_peak + 1
            continue

        # Va — manual mid value
        Va = Vmin + (Vpeak - Vmin) / 2.0

        # bracket: smallest k with s[k] >= Va on leading edge
        k = cross
        while k <= i_peak and s[k] < Va:
            k += 1
        if k > i_peak:
            quality = 0x08    # Q_DAQ_VA_OUT_OF_RANGE
            coarse = i_peak
            fine = 0
            Vba = Vaa = float(s[i_peak])
        else:
            Vaa = float(s[k])
            Vba = float(s[k - 1]) if k > 0 else Vmin
            denom = (Vaa - Vba) if (Vaa > Vba) else 1.0
            fine = int(round((Va - Vba) / denom * 64.0))
            coarse = k - 1
            if fine >= 64:
                fine -= 64
                coarse += 1
            quality = 0x00

        time_units = coarse * 64 + fine
        time_ns = time_units * (CLK_NS / 64.0)

        wlo = cross - nsb_s
        whi = cross + nsa_s
        if wlo < 0:
            wlo = 0
            quality |= 0x02    # Q_DAQ_NSB_TRUNCATED
        if whi >= n:
            whi = n - 1
            quality |= 0x04    # Q_DAQ_NSA_TRUNCATED
        if i_peak >= n - 1:
            quality |= 0x01    # Q_DAQ_PEAK_AT_BOUNDARY

        integral = float(np.sum(sp[wlo:whi + 1]))

        return dict(
            vnoise=vnoise, vnoise_initial=vnoise_initial,
            sp=sp, nsb_samples=nsb_s, nsa_samples=nsa_s,
            cross=cross, peak=i_peak, Vpeak=Vpeak,
            Vmin=Vmin, Va=Va, Vba=Vba, Vaa=Vaa, k=k,
            coarse=coarse, fine=fine, time_units=time_units, time_ns=time_ns,
            window_lo=wlo, window_hi=whi, integral=integral, quality=quality,
        )
    return None


# ---------------------------------------------------------------------------
# WaveAnalyzer — soft analyzer (smoothing + iterative pedestal + local maxima)
# ---------------------------------------------------------------------------
def soft_analyze(s, *, smooth_order=2, threshold=5.0, min_threshold=3.0,
                 ped_nsamples=30, ped_flatness=1.0, ped_max_iter=3,
                 int_tail_ratio=0.1, tail_break_n=2, clk_mhz=250.0):
    n = len(s)

    # triangular smoothing
    if smooth_order <= 1:
        buf = s.astype(np.float64).copy()
    else:
        buf = np.empty(n)
        for i in range(n):
            val = s[i]; wsum = 1.0
            for j in range(1, smooth_order):
                if j > i or i + j >= n:
                    continue
                w = 1.0 - j / float(smooth_order + 1)
                val += w * (s[i - j] + s[i + j])
                wsum += 2.0 * w
            buf[i] = val / wsum

    # ── pedestal: median+MAD seed → iterative σ-clip
    nped = min(ped_nsamples, n)
    sc = buf[:nped].copy()
    mean = float(np.median(sc))
    rms  = float(np.median(np.abs(sc - mean))) * 1.4826
    for _ in range(ped_max_iter):
        keep = np.abs(sc - mean) < max(rms, ped_flatness)
        if keep.sum() == sc.size or keep.sum() < 5:
            break
        sc = sc[keep]
        mean = float(np.mean(sc)); rms = float(np.std(sc))
    nused = int(sc.size)

    thr = max(threshold * rms, min_threshold)

    # find first local maximum above threshold
    for i in range(1, n - 1):
        if buf[i] < buf[i - 1] or buf[i] < buf[i + 1]:
            continue
        # walk left/right
        left = i; right = i
        while left > 0 and buf[left] > buf[left - 1]:
            left -= 1
        while right < n - 1 and buf[right] >= buf[right + 1]:
            right += 1
        height_smooth = buf[i] - mean
        if height_smooth < thr:
            continue
        # (the redundant `height_smooth < 3 * rms` check has been removed;
        #  thr ≥ 5*rms in the default config already binds first.)

        # raw position correction
        raw_pos = i; raw_h = s[i] - mean
        for j in range(1, smooth_order + 1):
            if i - j >= 0 and (s[i - j] - mean) > raw_h:
                raw_pos = i - j; raw_h = s[i - j] - mean
            if i + j < n and (s[i + j] - mean) > raw_h:
                raw_pos = i + j; raw_h = s[i + j] - mean

        # ── integration: walk outward with N-consecutive tail termination,
        # INCLUSIVE int_left/int_right (only updated on commit).
        ped_height = buf[i] - mean
        tail_cut = ped_height * int_tail_ratio
        int_left = i; int_right = i
        integ = ped_height
        N_break = max(1, tail_break_n)

        # leftward
        below_run = 0; pending = 0.0
        for j in range(i - 1, left - 1, -1):
            v = buf[j] - mean
            if v < tail_cut or v < rms:
                below_run += 1; pending += v
                if below_run >= N_break: break
            else:
                integ += pending + v
                pending = 0.0; below_run = 0
                int_left = j
        # rightward
        below_run = 0; pending = 0.0
        for j in range(i + 1, right + 1):
            v = buf[j] - mean
            if v < tail_cut or v < rms:
                below_run += 1; pending += v
                if below_run >= N_break: break
            else:
                integ += pending + v
                pending = 0.0; below_run = 0
                int_right = j

        # ── sub-sample peak time: 3-point quadratic vertex around raw_pos
        t_subsample = 0.0
        if 0 < raw_pos < n - 1:
            h_minus = float(s[raw_pos - 1])
            h_zero  = float(s[raw_pos])
            h_plus  = float(s[raw_pos + 1])
            denom = h_minus - 2.0 * h_zero + h_plus
            if denom < -1e-3:
                delta = 0.5 * (h_minus - h_plus) / denom
                t_subsample = max(-1.0, min(1.0, delta))

        return dict(
            ped_mean=mean, ped_rms=rms, ped_nused=nused,
            threshold=thr, buf=buf,
            pos=raw_pos, height=float(raw_h), t_subsample=t_subsample,
            int_left=int_left, int_right=int_right, integral=float(integ),
            time_ns=(raw_pos + t_subsample) * 1e3 / clk_mhz,
            tail_cut=tail_cut,
        )
    return None


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
print()
print(f"Wrote 5 PNGs to {out_dir}")
