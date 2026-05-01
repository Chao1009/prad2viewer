#!/usr/bin/env python3
"""
plot_gem_clustering.py — visualisations for gem::GemCluster.

Illustrates the strip-level clustering pipeline + 2-D X/Y matching:

  figs/gem_fig1_layout.png         — 4 GEM detectors, X/Y planes, beam hole
  figs/gem_fig2_strip_clustering.png — group + split + charge-weighted position
  figs/gem_fig3_xy_matching.png    — Cartesian product with cuts vs ADC-sorted
  figs/gem_fig4_params.png         — split_thres effect + cross-talk distances

Run:
  python plot_gem_clustering.py
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from pathlib import Path

HERE = Path(__file__).parent
FIGS = HERE / 'figs'
FIGS.mkdir(exist_ok=True)
DB   = HERE.parent.parent / 'database' / 'gem_daq_map.json'


# ---------------------------------------------------------------------------
# Strip-mapping (mirrors gem::buildStripMap in prad2det) — used to derive
# the actual X/Y plane size and beam-hole position from gem_daq_map.json.
# ---------------------------------------------------------------------------
def map_strip(ch, pos, orient, pin_rotate=0, shared_pos=-1,
              hybrid_board=True, N=128, readout_center=32):
    readout_off = readout_center + pin_rotate
    eff_pos = shared_pos if shared_pos >= 0 else pos
    plane_shift = (eff_pos - pos) * N - pin_rotate
    strip = 32 * (ch % 4) + 8 * (ch // 4) - 31 * (ch // 16)
    if hybrid_board:
        strip = strip + 1 + strip % 4 - 5 * ((strip // 4) % 2)
    if readout_off > 0:
        if strip & 1:
            strip = readout_off - (strip + 1) // 2
        else:
            strip = readout_off + strip // 2
    strip &= (N - 1)
    if orient == 1:
        strip = (N - 1) - strip
    strip += plane_shift + pos * N
    return strip


_raw = json.loads(DB.read_text(encoding='utf-8'))
_apv_ch = _raw.get('apv_channels', 128)
_ro_c   = _raw.get('readout_center', 32)
_apvs   = [a for a in _raw['apvs'] if isinstance(a, dict) and 'crate' in a]

# All four GEMs are identical — derive geometry from det 0
_det0 = [a for a in _apvs if a.get('det') == 0]
def _strips_of(plane):
    out = []
    match = []
    for a in (x for x in _det0 if x.get('plane') == plane):
        pin_r = a.get('pin_rotate') or 0
        sp    = a.get('shared_pos') if a.get('shared_pos') is not None else -1
        for ch in range(_apv_ch):
            s = map_strip(ch, a['pos'], a.get('orient', 0), pin_r, sp,
                          a.get('hybrid_board', True), _apv_ch, _ro_c)
            out.append(s)
            if a.get('match'):
                match.append(s)
    return out, match

_xs, _xs_match = _strips_of('X')
_ys, _ys_match = _strips_of('Y')

PITCH  = _raw['layers'][0].get('x_pitch', 0.4)
HOLE_W = _raw['hole']['width']
HOLE_H = _raw['hole']['height']
DET_W  = (max(_xs) + 1) * PITCH                              # 563.2 mm
DET_H  = (max(_ys) + 1) * PITCH                              # 1228.8 mm
HOLE_CX = (min(_xs_match) + max(_xs_match) + 1) / 2 * PITCH  # 534.4 mm — off-centre
HOLE_CY = DET_H / 2                                          # 614.4 mm


# ---------------------------------------------------------------------------
# GEM cluster algorithm — Python implementation (matches GemCluster.cpp)
# ---------------------------------------------------------------------------
def group_and_split(hits, *, consecutive_thres=1, split_thres=14.0,
                    min_hits=1, max_hits=20):
    """Run group + recursive valley split over a list of (strip, charge) tuples.
    Returns list of clusters; each cluster is a dict with hits + position +
    peak/total charge."""
    if not hits:
        return []
    hits = sorted(hits, key=lambda h: h[0])

    # group consecutive
    groups = []
    cur = [hits[0]]
    for h in hits[1:]:
        if h[0] - cur[-1][0] > consecutive_thres:
            groups.append(cur)
            cur = [h]
        else:
            cur.append(h)
    groups.append(cur)

    out = []
    for g in groups:
        out.extend(_split_recursive(g, split_thres))

    # filter by size
    out = [c for c in out if min_hits <= len(c['hits']) <= max_hits]

    # charge-weighted position
    for c in out:
        tot = sum(h[1] for h in c['hits'])
        if tot > 0:
            c['position'] = sum(h[0] * h[1] for h in c['hits']) / tot * PITCH
        else:
            c['position'] = c['hits'][0][0] * PITCH
        c['total_charge'] = tot
        c['peak_charge']  = max(h[1] for h in c['hits'])
    return out


def _split_recursive(group, thres):
    """Recursively split a group at the first local minimum (valley) where
    the upturn exceeds `thres`."""
    if len(group) < 3:
        return [{'hits': group}]

    descending = False
    minimum_idx = 0
    extremum = False

    for i in range(len(group) - 1):
        c, cn = group[i][1], group[i + 1][1]
        if descending:
            if c < group[minimum_idx][1]:
                minimum_idx = i
            if cn - c > thres:        # upturn confirms valley
                extremum = True
                break
        else:
            if c - cn > thres:        # descending → potential valley ahead
                descending = True
                minimum_idx = i + 1

    if not extremum:
        return [{'hits': group}]

    # halve the charge of the valley strip (overlap)
    valley = list(group[minimum_idx])
    valley[1] /= 2.0

    left  = group[:minimum_idx]
    right = [tuple(valley)] + group[minimum_idx + 1:]
    return [{'hits': left}] + _split_recursive(right, thres)


# ---------------------------------------------------------------------------
# Plot 1 — GEM detector layout
# ---------------------------------------------------------------------------
# Coordinate convention: origin at the bottom-left of the detector, matching
# the strip-numbering convention (strip 0 ≡ x = 0 / y = 0).  The beam hole
# is *not* at the centre — it sits at (HOLE_CX, HOLE_CY) on the right side
# of the X plane because APV positions 10/11 carry the half-strip "match"
# pair that wraps around it.
fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6.0),
                                gridspec_kw={'width_ratios': [1, 1.4]})

ax = axA
# detector outline
ax.add_patch(Rectangle((0, 0), DET_W, DET_H,
                        facecolor='#e7f5ff', edgecolor='#1971c2', lw=1.5))

# X-strip indications: vertical lines (constant x), broken at the hole
n_xlines = 9
for x in np.linspace(20, DET_W - 20, n_xlines):
    if HOLE_CX - HOLE_W / 2 <= x <= HOLE_CX + HOLE_W / 2:
        ax.plot([x, x], [0, HOLE_CY - HOLE_H / 2], color='#74c0fc', lw=0.4, alpha=0.7)
        ax.plot([x, x], [HOLE_CY + HOLE_H / 2, DET_H], color='#74c0fc', lw=0.4, alpha=0.7)
    else:
        ax.plot([x, x], [0, DET_H], color='#74c0fc', lw=0.4, alpha=0.7)

# Y-strip indications: horizontal lines (constant y), broken at the hole
n_ylines = 17
for y in np.linspace(40, DET_H - 40, n_ylines):
    if HOLE_CY - HOLE_H / 2 <= y <= HOLE_CY + HOLE_H / 2:
        ax.plot([0, HOLE_CX - HOLE_W / 2], [y, y], color='#fcc2d7', lw=0.4, alpha=0.7)
        ax.plot([HOLE_CX + HOLE_W / 2, DET_W], [y, y], color='#fcc2d7', lw=0.4, alpha=0.7)
    else:
        ax.plot([0, DET_W], [y, y], color='#fcc2d7', lw=0.4, alpha=0.7)

# beam hole
ax.add_patch(Rectangle((HOLE_CX - HOLE_W / 2, HOLE_CY - HOLE_H / 2),
                        HOLE_W, HOLE_H,
                        facecolor='#fff3bf', edgecolor='#fab005', lw=1.5,
                        zorder=5))
ax.annotate('beam hole\n52 × 52 mm',
            xy=(HOLE_CX, HOLE_CY),
            xytext=(HOLE_CX - 220, HOLE_CY + 80),
            fontsize=8, color='#444', ha='center',
            arrowprops=dict(arrowstyle='->', lw=0.9, color='#666'))

# strip-direction labels
ax.annotate('', xy=(-22, 60), xytext=(-22, 200),
            arrowprops=dict(arrowstyle='<->', lw=1.2, color='#1971c2'))
ax.text(-30, 130, 'X strips\n(vertical,\nmeasure x)',
        ha='right', va='center', fontsize=8, color='#1971c2')
ax.annotate('', xy=(40, DET_H + 22), xytext=(180, DET_H + 22),
            arrowprops=dict(arrowstyle='<->', lw=1.2, color='#e64980'))
ax.text(110, DET_H + 35, 'Y strips (horizontal, measure y)',
        ha='center', va='bottom', fontsize=8, color='#e64980')

ax.set_xlim(-90, DET_W + 30)
ax.set_ylim(-30, DET_H + 70)
ax.set_aspect('equal')
ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
ax.set_title(f"Single GEM detector — beam-direction view\n"
             f"X: 12 APVs × 128 ch ({DET_W:.1f} mm)  ·  "
             f"Y: 24 APVs × 128 ch ({DET_H:.1f} mm)")
ax.grid(False)

# Panel B — beam-axis arrangement of 4 detectors.  Vertical extent is
# the full Y plane size (DET_H = 1228.8 mm); the beam crosses each
# detector at the hole y-centre (HOLE_CY = DET_H/2).
ax = axB
detectors = [
    ('GEM0', 5407 + 39.71 / 2, '#1971c2'),
    ('GEM1', 5407 - 39.71 / 2, '#74c0fc'),
    ('GEM2', 5807 + 39.71 / 2, '#e64980'),
    ('GEM3', 5807 - 39.71 / 2, '#faa2c1'),
]
for i, (name, z, col) in enumerate(detectors):
    ax.add_patch(Rectangle((z - 5, 0), 10, DET_H,
                            facecolor=col, edgecolor='black', lw=0.8, alpha=0.6))
    # alternate label heights so paired detectors don't collide
    yoff = DET_H + (15 if i % 2 == 0 else 90)
    ax.text(z, yoff, name, ha='center', va='bottom',
            fontsize=9, color=col, fontweight='bold')
    # mark the beam-hole opening in each detector
    ax.add_patch(Rectangle((z - 5, HOLE_CY - HOLE_H / 2), 10, HOLE_H,
                            facecolor='#fff3bf', edgecolor='#fab005', lw=0.8))

# beam line (passes through the hole y-centre)
ax.annotate('', xy=(6050, HOLE_CY), xytext=(5200, HOLE_CY),
            arrowprops=dict(arrowstyle='->', lw=1.5, color='#444'))
ax.text(5300, HOLE_CY + 40, 'beam', fontsize=9, color='#444')

ax.text(5407, -100, 'Layer 1\n(z ≈ 5407 mm)',
        ha='center', va='top', fontsize=9, color='#666')
ax.text(5807, -100, 'Layer 2\n(z ≈ 5807 mm)',
        ha='center', va='top', fontsize=9, color='#666')

ax.set_xlim(5200, 6050)
ax.set_ylim(-180, DET_H + 200)
ax.set_xlabel('z (mm, lab)')
ax.set_ylabel('y (mm, detector frame)')
ax.set_title("4-GEM stack — paired layers (Δz = 39.7 mm)\n"
             "yellow strip = beam hole (52 mm tall, on beam axis)")
ax.grid(alpha=0.3)
ax.set_aspect('auto')

fig.tight_layout()
fig.savefig(FIGS / 'gem_fig1_layout.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 2 — strip clustering: group + split + position
# ---------------------------------------------------------------------------
# Build a synthetic strip ADC trace with three groups:
#   Group A (strips 110-114): single 4-strip cluster, peak at 112
#   Group B (strips 130-139): multi-peak cluster (split at valley)
#   Group C (strips 158-160): small 3-strip cluster
def gauss(s, mu, sigma, amp):
    return amp * np.exp(-((s - mu) ** 2) / (2 * sigma ** 2))


strips = np.arange(100, 175)
adc = np.zeros_like(strips, dtype=float)
# Group A — single peak
for i, s in enumerate(strips):
    adc[i] += gauss(s, 112.2, 1.1, 800)
# Group B — multi-peak (two overlapping showers, valley at 134)
for i, s in enumerate(strips):
    adc[i] += gauss(s, 132.0, 1.2, 600) + gauss(s, 137.0, 1.0, 480)
# Group C — small
for i, s in enumerate(strips):
    adc[i] += gauss(s, 159.0, 0.8, 300)
# noise floor — random walk small fluctuation
rng = np.random.default_rng(0)
adc += rng.normal(0, 4, size=adc.size)
# Threshold of 30 ADC (realistic zero-suppression cut)
hits = [(int(s), float(a)) for s, a in zip(strips, adc) if a >= 30.0]

clusters = group_and_split(hits, consecutive_thres=1, split_thres=14.0)

# render
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.0),
                                gridspec_kw={'width_ratios': [1.6, 1]})

# Panel A — full strip ADC + cluster spans
ax = axA
ax.bar([h[0] for h in hits], [h[1] for h in hits], width=0.85,
       color='#74c0fc', edgecolor='#1971c2', lw=0.5, label='strip ADC ≥ 30')

# colour-shade each cluster span
clu_colors = ['#fa8072', '#ffa94d', '#5cb85c', '#9775fa', '#f06595']
for i, c in enumerate(clusters):
    s_lo = min(h[0] for h in c['hits'])
    s_hi = max(h[0] for h in c['hits'])
    ax.axvspan(s_lo - 0.5, s_hi + 0.5, color=clu_colors[i % len(clu_colors)],
               alpha=0.2)
    # mark charge-weighted position
    ax.axvline(c['position'] / PITCH, color=clu_colors[i % len(clu_colors)],
               lw=1.5, ls='--')
    ax.text(c['position'] / PITCH, max(adc) * 0.95,
            f"x = {c['position']:.2f} mm\nΣ = {c['total_charge']:.0f}",
            ha='center', va='top', fontsize=8,
            color=clu_colors[i % len(clu_colors)],
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=2))

ax.set_xlabel('strip number (pitch 0.4 mm)')
ax.set_ylabel('ADC counts')
ax.set_title(f"Strip clustering — {len(clusters)} clusters from "
             f"{len(hits)} above-threshold strips\n"
             f"shaded = consecutive group, dashed = charge-weighted position")
ax.set_xlim(105, 165)
ax.grid(alpha=0.3, axis='y')
ax.legend(loc='upper right', fontsize=8)

# Panel B — zoom on the multi-peak group, showing the split
ax = axB
g_strips = [h[0] for h in hits if 128 <= h[0] <= 140]
g_adc    = [h[1] for h in hits if 128 <= h[0] <= 140]
ax.bar(g_strips, g_adc, width=0.85, color='#74c0fc',
       edgecolor='#1971c2', lw=0.5)

# find the valley strip
valley = None
for i in range(1, len(g_adc) - 1):
    if g_adc[i] < g_adc[i - 1] and g_adc[i] < g_adc[i + 1]:
        if g_adc[i + 1] - g_adc[i] > 14.0:
            valley = g_strips[i]
            break

if valley is not None:
    ax.axvline(valley + 0.5, color='#d62728', lw=2.0, ls='--',
               label=f'split at strip {valley} (valley)')
    ax.annotate('', xy=(valley + 0.5, g_adc[g_strips.index(valley)] + 100),
                xytext=(valley + 0.5, g_adc[g_strips.index(valley)] + 250),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='#d62728'))
    # split_thres marker
    if g_strips.index(valley) + 1 < len(g_adc):
        rise = g_adc[g_strips.index(valley) + 1] - g_adc[g_strips.index(valley)]
        ax.text(valley + 1.6, g_adc[g_strips.index(valley)] + rise / 2,
                f"upturn = {rise:.0f}\n> split_thres = 14",
                fontsize=8, color='#d62728', va='center')

# overlay the two reconstructed positions for this group
g_clusters = [c for c in clusters
              if 128 <= min(h[0] for h in c['hits']) <= 140]
for i, c in enumerate(g_clusters):
    ax.axvline(c['position'] / PITCH, color='#2ca02c', lw=1.5, ls=':')
    ax.text(c['position'] / PITCH, max(g_adc) * 0.92,
            f"x = {c['position']:.2f} mm",
            ha='center', va='top', fontsize=8, color='#2ca02c',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=2))

ax.set_xlabel('strip number')
ax.set_ylabel('ADC counts')
ax.set_title("Recursive valley split\n(left + right become two clusters)")
ax.legend(loc='lower right', fontsize=8)
ax.grid(alpha=0.3, axis='y')

fig.tight_layout()
fig.savefig(FIGS / 'gem_fig2_strip_clustering.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 3 — X/Y matching (Cartesian + cuts vs ADC-sorted)
# ---------------------------------------------------------------------------
# Three X clusters and three Y clusters with peak charges + mean times
xs = [
    {'pos': -180.0, 'peak': 1200.0, 't_ns': 60.0},   # prompt big shower
    {'pos':   45.0, 'peak':  900.0, 't_ns': 75.0},   # prompt big shower
    {'pos':  120.0, 'peak':  150.0, 't_ns': 145.0},  # late, small (out-of-time noise)
]
ys = [
    {'pos':   80.0, 'peak': 1100.0, 't_ns': 65.0},
    {'pos': -160.0, 'peak':  850.0, 't_ns': 70.0},
    {'pos':   30.0, 'peak':  140.0, 't_ns': 150.0},
]

ADC_ASYM_CUT = 0.8
TIME_CUT     = 50.0

# Cartesian + cuts — apply both cuts
matches_cart = []
for ix, x in enumerate(xs):
    for iy, y in enumerate(ys):
        sum_p = x['peak'] + y['peak']
        asym  = abs(x['peak'] - y['peak']) / sum_p
        dt    = abs(x['t_ns'] - y['t_ns'])
        ok    = (asym <= ADC_ASYM_CUT) and (dt <= TIME_CUT)
        matches_cart.append({'x': x, 'y': y, 'asym': asym, 'dt': dt, 'ok': ok,
                             'ix': ix, 'iy': iy})

# ADC-sorted: pair by descending peak charge
xs_sorted = sorted(range(len(xs)), key=lambda i: -xs[i]['peak'])
ys_sorted = sorted(range(len(ys)), key=lambda i: -ys[i]['peak'])
matches_adc = list(zip(xs_sorted, ys_sorted))[:min(len(xs), len(ys))]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.5))

# Panel A: Cartesian + cuts
ax = axA
# layout: x-axis = X-cluster position, top-margin = peak; y-axis = Y-cluster position
for x in xs:
    ax.axvline(x['pos'], color='#1971c2', lw=0.6, alpha=0.4)
for y in ys:
    ax.axhline(y['pos'], color='#e64980', lw=0.6, alpha=0.4)
# label X clusters
for i, x in enumerate(xs):
    ax.text(x['pos'], -260, f"X{i}\nQ={x['peak']:.0f}\nt={x['t_ns']:.0f}",
            ha='center', va='top', fontsize=8, color='#1971c2')
# label Y clusters
for i, y in enumerate(ys):
    ax.text(-260, y['pos'], f"Y{i}\nQ={y['peak']:.0f}\nt={y['t_ns']:.0f}",
            ha='right', va='center', fontsize=8, color='#e64980')

# all 9 candidate hits
for m in matches_cart:
    px, py = m['x']['pos'], m['y']['pos']
    if m['ok']:
        ax.plot(px, py, 'o', color='#2ca02c', ms=12, mec='black', mew=0.8)
        ax.text(px + 6, py + 6,
                f"asym={m['asym']:.2f}\nΔt={m['dt']:.0f} ns",
                fontsize=7, color='#2ca02c')
    else:
        ax.plot(px, py, 'x', color='#d62728', ms=12, mew=2)
        # annotate why it failed
        reason = []
        if m['asym'] > ADC_ASYM_CUT:
            reason.append(f"asym {m['asym']:.2f}>{ADC_ASYM_CUT}")
        if m['dt'] > TIME_CUT:
            reason.append(f"Δt {m['dt']:.0f}>{TIME_CUT}")
        ax.text(px + 6, py + 6, '\n'.join(reason),
                fontsize=7, color='#d62728')

ax.set_xlim(-260, 200); ax.set_ylim(-260, 180)
ax.set_xlabel('X cluster position (mm)')
ax.set_ylabel('Y cluster position (mm)')
ax.set_title(f"Mode 1 — Cartesian + cuts\n"
             f"|asym| ≤ {ADC_ASYM_CUT}, |Δt| ≤ {TIME_CUT} ns "
             f"({sum(1 for m in matches_cart if m['ok'])} hits / "
             f"{len(matches_cart)} candidates)")
ax.grid(alpha=0.3)

leg = [Line2D([], [], marker='o', color='w', markerfacecolor='#2ca02c',
              markeredgecolor='black', markersize=10, label='accepted (2-D hit)'),
       Line2D([], [], marker='x', color='#d62728', linestyle='', markersize=10,
              mew=2, label='rejected by a cut')]
ax.legend(handles=leg, loc='upper right', fontsize=8)

# Panel B: ADC-sorted matching
ax = axB
ax.set_xlim(-260, 200); ax.set_ylim(-260, 180)

# show the same X/Y cluster grid faintly
for x in xs:
    ax.axvline(x['pos'], color='#1971c2', lw=0.6, alpha=0.3)
for y in ys:
    ax.axhline(y['pos'], color='#e64980', lw=0.6, alpha=0.3)
for i, x in enumerate(xs):
    rank = xs_sorted.index(i) + 1
    ax.text(x['pos'], -260, f"X{i}\nQ={x['peak']:.0f}\nrank={rank}",
            ha='center', va='top', fontsize=8, color='#1971c2')
for i, y in enumerate(ys):
    rank = ys_sorted.index(i) + 1
    ax.text(-260, y['pos'], f"Y{i}\nQ={y['peak']:.0f}\nrank={rank}",
            ha='right', va='center', fontsize=8, color='#e64980')

for ix, iy in matches_adc:
    px, py = xs[ix]['pos'], ys[iy]['pos']
    ax.plot(px, py, 's', color='#fab005', ms=14, mec='black', mew=0.8)
    ax.text(px + 6, py + 6, f"X{ix} ↔ Y{iy}\npaired by rank",
            fontsize=7, color='#a16207')

ax.set_xlabel('X cluster position (mm)')
ax.set_ylabel('Y cluster position (mm)')
ax.set_title(f"Mode 0 — ADC-sorted 1:1\n"
             "rank by descending peak ADC, pair index-wise "
             f"({len(matches_adc)} hits)")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(FIGS / 'gem_fig3_xy_matching.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 4 — parameter sensitivity: split_thres + cross-talk
# ---------------------------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 4.8))

# Panel A: split_thres effect — synthetic data with two valleys of
# different depths so each threshold yields a distinct cluster count.
ax = axA
B_strips = list(range(128, 139))
B_adc    = [40, 110, 600, 200, 60, 240, 70, 100, 80, 50, 35]
B_hits   = list(zip(B_strips, B_adc))

ax.bar(B_strips, B_adc, width=0.85, color='#dee2e6',
       edgecolor='#868e96', lw=0.5)

scenarios = [(200.0, '#2ca02c'), (50.0, '#1971c2'), (5.0, '#d62728')]
labels = []
for thr, col in scenarios:
    cls = group_and_split(B_hits, split_thres=thr)
    labels.append(f"split_thres = {thr:>4.0f} → {len(cls)} cluster"
                  f"{'s' if len(cls) != 1 else ''}")
    for c in cls:
        s_lo = min(h[0] for h in c['hits'])
        s_hi = max(h[0] for h in c['hits'])
        ax.plot([s_lo - 0.4, s_hi + 0.4],
                [max(B_adc) * (1.05 + scenarios.index((thr, col)) * 0.10)] * 2,
                color=col, lw=4, solid_capstyle='butt')
        ax.plot([(s_lo + s_hi) / 2.0],
                [max(B_adc) * (1.05 + scenarios.index((thr, col)) * 0.10)],
                'v', color=col, ms=6)

ax.set_xlabel('strip number')
ax.set_ylabel('ADC counts')
ax.set_title("split_thres on a multi-peak group\n"
             "lower threshold → more aggressive splitting")

leg_lines = [Line2D([], [], color=col, lw=4, label=lbl)
             for (_, col), lbl in zip(scenarios, labels)]
ax.legend(handles=leg_lines, loc='lower right', fontsize=8)
ax.grid(alpha=0.3, axis='y')
ax.set_ylim(0, max(B_adc) * 1.42)

# Panel B: cross-talk identification at characteristic distances
ax = axB
# Synthesize a primary cluster + ghost cluster at one of the characteristic
# distances (24.4 mm = 61 strips at 0.4 mm pitch)
charac = [6.4, 17.6, 24.4, 24.8, 25.2, 25.6, 26.0, 26.4, 26.8, 33.6, 44.8]
prim_pos  = 50.0          # mm
ghost_pos = prim_pos + 24.4
strips2   = np.arange(50, 250)
adc2      = (gauss(strips2, prim_pos / PITCH,  1.1, 1500)
             + gauss(strips2, ghost_pos / PITCH, 0.9, 90)        # cross-talk
             + rng.normal(0, 3, size=strips2.size))
adc2 = np.maximum(adc2, 0)

ax.bar(strips2 * PITCH, adc2, width=PITCH * 0.9, color='#dee2e6',
       edgecolor='#868e96', lw=0.4)
# above-threshold strips that pass to clustering — use 30 ADC threshold
mask = adc2 >= 30.0
ax.bar((strips2 * PITCH)[mask], adc2[mask], width=PITCH * 0.9,
       color='#74c0fc', edgecolor='#1971c2', lw=0.4)

# mark primary and ghost
ax.axvline(prim_pos,  color='#2ca02c', lw=1.5, ls='--',
           label=f'primary  Q≈{1500:.0f}')
ax.axvline(ghost_pos, color='#d62728', lw=1.5, ls='--',
           label=f'ghost (CT)  at Δ = 24.4 mm')

# show all characteristic distances as faint vertical lines; only the
# well-separated ones get a numeric label (the 24.4..26.8 cluster is
# annotated as a single range to avoid overlap)
labelled = {6.4: '6.4', 17.6: '17.6', 25.6: '24.4–26.8', 33.6: '33.6', 44.8: '44.8'}
for d in charac:
    ax.axvline(prim_pos + d, color='#fa5252', lw=0.5, ls=':', alpha=0.4)
for d, txt in labelled.items():
    ax.text(prim_pos + d, max(adc2) * 1.02, txt,
            rotation=0, fontsize=7.5, ha='center', va='bottom',
            color='#c92a2a', alpha=0.85)

ax.set_xlim(40, 110)
ax.set_xlabel('position (mm)')
ax.set_ylabel('ADC counts')
ax.set_title("Cross-talk at characteristic distances\n"
             "small clusters at Δ ∈ {6.4, 17.6, 24.4 … 44.8} mm flagged + dropped")
ax.legend(loc='upper right', fontsize=8)
ax.grid(alpha=0.3, axis='y')

fig.tight_layout()
fig.savefig(FIGS / 'gem_fig4_params.png', dpi=130)
plt.close(fig)

print(f"Strip clustering: {len(clusters)} clusters")
for c in clusters:
    print(f"  hits {len(c['hits']):2d}  pos {c['position']:7.2f} mm"
          f"  Σ={c['total_charge']:6.0f}  peak={c['peak_charge']:6.0f}")

print(f"\nWrote 4 PNGs to {FIGS}")
