#!/usr/bin/env python3
"""
plot_hycal_clustering.py — visualisations for fdec::HyCalCluster.

Loads the actual HyCal geometry from `database/hycal_modules.json` and
illustrates the four steps of the Island clustering algorithm:

  hycal_fig1_layout.png         — full HyCal sectors + module types
  hycal_fig2_single_cluster.png — DFS grouping, seed, log-weighted centroid
  hycal_fig3_split.png          — two overlapping showers, profile-based split
  hycal_fig4_params.png         — log_weight_thres effect + shower-depth curves

Run:
  python plot_hycal_clustering.py
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
from pathlib import Path

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
HERE   = Path(__file__).parent
DB     = HERE.parent.parent / 'database' / 'hycal_modules.json'
MODS   = json.loads(DB.read_text())
HYCAL  = [m for m in MODS if m['t'] in ('PbGlass', 'PbWO4')]


def find_mod(x, y, t=None):
    """Return module dict whose centre is closest to (x, y) (optionally
    filtered by type)."""
    pool = [m for m in HYCAL if t is None or m['t'] == t]
    return min(pool, key=lambda m: (m['x'] - x) ** 2 + (m['y'] - y) ** 2)


def neighbors(center, radius=1):
    """Return modules whose grid offset from `center` is within `radius`."""
    sx, sy = center['sx'], center['sy']
    out = []
    for m in HYCAL:
        if m['t'] != center['t']:
            continue
        dx = (m['x'] - center['x']) / sx
        dy = (m['y'] - center['y']) / sy
        if abs(dx) <= radius + 0.01 and abs(dy) <= radius + 0.01:
            out.append((m, dx, dy))
    return out


# ---------------------------------------------------------------------------
# Synthetic shower model: 2-D Gaussian, total ~E in MeV
# ---------------------------------------------------------------------------
def shower_energy(modules_with_dxy, x0, y0, E_tot, sigma_mm):
    """Distribute `E_tot` over the listed modules using a Gaussian profile
    centred at (x0, y0) with width `sigma_mm`."""
    raw = []
    for m, _, _ in modules_with_dxy:
        d2 = (m['x'] - x0) ** 2 + (m['y'] - y0) ** 2
        raw.append(np.exp(-d2 / (2 * sigma_mm ** 2)))
    raw = np.array(raw)
    return E_tot * raw / raw.sum()


def log_weighted_centroid(positions, energies, E_tot, thres=3.6):
    """Return (x, y) using the log-weight scheme in HyCalCluster::get_weight."""
    wsum = 0.0; wx = 0.0; wy = 0.0
    for (px, py), E in zip(positions, energies):
        if E <= 0 or E_tot <= 0:
            continue
        w = thres + np.log(E / E_tot)
        if w > 0:
            wsum += w
            wx   += w * px
            wy   += w * py
    if wsum <= 0:
        return None
    return wx / wsum, wy / wsum


def energy_weighted_centroid(positions, energies):
    e = np.array(energies)
    if e.sum() <= 0:
        return None
    px = np.array([p[0] for p in positions])
    py = np.array([p[1] for p in positions])
    return float((px * e).sum() / e.sum()), float((py * e).sum() / e.sum())


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
COL_GLASS = '#fff3bf'
COL_PWO4  = '#dbeafe'
COL_EDGE  = '#888'

def draw_module(ax, m, facecolor=None, edgecolor=None, lw=0.4, alpha=1.0,
                hatch=None):
    fc = facecolor if facecolor else (COL_PWO4 if m['t'] == 'PbWO4' else COL_GLASS)
    ec = edgecolor if edgecolor else COL_EDGE
    ax.add_patch(Rectangle(
        (m['x'] - m['sx'] / 2, m['y'] - m['sy'] / 2),
        m['sx'], m['sy'], facecolor=fc, edgecolor=ec, lw=lw, alpha=alpha,
        hatch=hatch))


# ---------------------------------------------------------------------------
# Plot 1 — HyCal layout overview
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 7.5))
for m in HYCAL:
    draw_module(ax, m)

# beam hole — central ~80×80 mm² gap (no PbWO4 modules in the inner block)
# inferred from the module layout, but we'll just leave the empty space.

# annotate sector names
sector_pos = {'Top': (0, 410), 'Bottom': (0, -410),
              'Left': (-410, 0), 'Right': (410, 0), 'Center': (0, -100)}
for name, (px, py) in sector_pos.items():
    ax.text(px, py, name, fontsize=11, ha='center', va='center',
            color='#333', fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))

# legend swatches
swatch_glass = Rectangle((0, 0), 1, 1, facecolor=COL_GLASS, edgecolor=COL_EDGE)
swatch_pwo4  = Rectangle((0, 0), 1, 1, facecolor=COL_PWO4,  edgecolor=COL_EDGE)
ax.legend([swatch_glass, swatch_pwo4],
          ['PbGlass (576 modules, 38.15 mm)',
           'PbWO₄ (1152 modules, 20.77 mm)'],
          loc='lower right', fontsize=9, framealpha=0.95)

ax.set_xlim(-600, 600); ax.set_ylim(-600, 600)
ax.set_aspect('equal')
ax.set_xlabel('x (mm, lab)'); ax.set_ylabel('y (mm, lab)')
ax.set_title('HyCal module layout — 5 sectors, two crystal types')
ax.grid(False)
fig.tight_layout()
fig.savefig(HERE / 'hycal_fig1_layout.png', dpi=130)
plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 2 — single cluster: shower in 5×5 PbWO4 region
# ---------------------------------------------------------------------------
# pick a centre well inside the PbWO4 region
center = find_mod(60, 100, t='PbWO4')
nbrs   = neighbors(center, radius=2)   # 5×5 region

# synthesise an EM shower with peak slightly off-centre to make centroid
# reconstruction visible
x0, y0 = center['x'] + 5.0, center['y'] - 4.0
energies = shower_energy(nbrs, x0, y0, E_tot=1100.0, sigma_mm=15.0)

# threshold (mimic min_module_energy = 1 MeV)
THR = 1.0
hit_mask = energies >= THR
seed_idx = int(np.argmax(energies))
seed     = nbrs[seed_idx][0]

# log-weighted centroid using all hits in the 3×3 around the seed
pos_3x3 = []; e_3x3 = []
for (m, dx, dy), E in zip(nbrs, energies):
    if not hit_mask[nbrs.index((m, dx, dy))] and E < THR:
        continue
    if abs(m['x'] - seed['x']) < seed['sx'] * 1.01 and \
       abs(m['y'] - seed['y']) < seed['sy'] * 1.01:
        pos_3x3.append((m['x'], m['y']))
        e_3x3.append(E)

E_total = float(energies[hit_mask].sum())
xc, yc = log_weighted_centroid(pos_3x3, e_3x3, E_total)
xe, ye = energy_weighted_centroid(pos_3x3, e_3x3)

fig, ax = plt.subplots(figsize=(7.5, 7.5))

# colour each module by deposited energy (log scale)
emax = energies.max()
for (m, _, _), E in zip(nbrs, energies):
    if E < THR:
        draw_module(ax, m, facecolor='#f8f9fa')
        continue
    intensity = (np.log10(E) - np.log10(THR)) / (np.log10(emax) - np.log10(THR))
    intensity = max(0.05, min(0.95, intensity))
    fc = plt.cm.YlOrRd(intensity)
    draw_module(ax, m, facecolor=fc, edgecolor='#666', lw=0.6)
    ax.text(m['x'], m['y'], f"{E:.0f}", ha='center', va='center',
            fontsize=8, color='#333')

# 3×3 box around seed (shower-position window)
ax.add_patch(Rectangle((seed['x'] - seed['sx'] * 1.5,
                        seed['y'] - seed['sy'] * 1.5),
                       seed['sx'] * 3, seed['sy'] * 3,
                       facecolor='none', edgecolor='#1f77b4', lw=2.0,
                       linestyle='--', label='3×3 position window'))

# seed star
ax.plot(seed['x'], seed['y'], '*', color='#d62728', ms=18, mec='white',
        mew=1.2, zorder=10)

# centroids
ax.plot(xc, yc, 'o', color='#1f77b4', ms=10, mec='white', mew=1.5,
        zorder=11)
ax.plot(xe, ye, 'D', color='#2ca02c', ms=8, mec='white', mew=1.0,
        zorder=11)
ax.plot(x0, y0, 'x', color='#000', ms=12, mew=2.0, zorder=12)

# x/y axis range: a bit wider than the 5×5 region
xlims = (center['x'] - seed['sx'] * 3, center['x'] + seed['sx'] * 3)
ylims = (center['y'] - seed['sy'] * 3, center['y'] + seed['sy'] * 3)
ax.set_xlim(xlims); ax.set_ylim(ylims)
ax.set_aspect('equal')

# legend
legend_elems = [
    Line2D([], [], marker='*', color='w', markerfacecolor='#d62728',
           markersize=15, label='seed (local max)'),
    Line2D([], [], marker='o', color='w', markerfacecolor='#1f77b4',
           markersize=10, label=f'log-weighted centroid (T={3.6})'),
    Line2D([], [], marker='D', color='w', markerfacecolor='#2ca02c',
           markersize=8, label='energy-weighted centroid'),
    Line2D([], [], marker='x', color='#000', linestyle='', markersize=10,
           label=f'true shower (E={1100:.0f} MeV)'),
    Line2D([], [], color='#1f77b4', linestyle='--', lw=2.0,
           label='3×3 position window'),
]
ax.legend(handles=legend_elems, loc='upper left', fontsize=8, framealpha=0.95)

ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
ax.set_title("Single cluster — energy in MeV, seed + position window\n"
             f"E_total = {E_total:.0f} MeV, σ_shower = 15 mm")
fig.tight_layout()
fig.savefig(HERE / 'hycal_fig2_single_cluster.png', dpi=130)
plt.close(fig)

print(f"Single cluster:")
print(f"  true position    : ({x0:.2f}, {y0:.2f}) mm")
print(f"  seed module      : {seed['n']} at ({seed['x']:.2f}, {seed['y']:.2f})")
print(f"  energy-weighted  : ({xe:.2f}, {ye:.2f}) mm  "
      f"|err| = {np.hypot(xe-x0, ye-y0):.2f} mm")
print(f"  log-weighted     : ({xc:.2f}, {yc:.2f}) mm  "
      f"|err| = {np.hypot(xc-x0, yc-y0):.2f} mm")


# ---------------------------------------------------------------------------
# Plot 3 — two showers, profile-based split
# ---------------------------------------------------------------------------
# Two showers separated by ~2 modules.  Use a 6×6 region of PbWO4.
center = find_mod(0, 0, t='PbWO4')
region = neighbors(center, radius=3)        # 7×7 region

x1, y1 = center['x'] - 25, center['y'] + 5
x2, y2 = center['x'] + 25, center['y'] - 8
E1, E2 = 1500.0, 900.0
SIG = 14.0

e1 = shower_energy(region, x1, y1, E1, SIG)
e2 = shower_energy(region, x2, y2, E2, SIG)
total = e1 + e2

# find local maxima — mod with energy > all PbWO4 grid neighbors
def local_max_indices(modules_with_dxy, energies, min_seed_E=10.0):
    out = []
    for i, (m, _, _) in enumerate(modules_with_dxy):
        if energies[i] < min_seed_E:
            continue
        is_max = True
        for j, (mn, _, _) in enumerate(modules_with_dxy):
            if i == j:
                continue
            dxg = abs(mn['x'] - m['x']) / m['sx']
            dyg = abs(mn['y'] - m['y']) / m['sy']
            if dxg <= 1.01 and dyg <= 1.01 and energies[j] > energies[i]:
                is_max = False
                break
        if is_max:
            out.append(i)
    return out

THR = 1.0
hit_mask = total >= THR
maxima = local_max_indices(region, total, min_seed_E=10.0)
seed1, seed2 = region[maxima[0]][0], region[maxima[1]][0]

# profile-based split fractions (the get_profile_frac scheme):
# weight by exp(-r²/2σ²) using SimpleProfile-equivalent σ ≈ 0.36 modules
SIG_GRID = 0.36

def profile_weight(m, mc):
    dx = (m['x'] - mc['x']) / mc['sx']
    dy = (m['y'] - mc['y']) / mc['sy']
    d2 = dx * dx + dy * dy
    return np.exp(-d2 / (2 * SIG_GRID ** 2))

E_seed1 = total[maxima[0]]; E_seed2 = total[maxima[1]]
fracs = np.zeros((len(region), 2))
for i, (m, _, _) in enumerate(region):
    if total[i] < THR:
        continue
    w1 = profile_weight(m, seed1) * E_seed1
    w2 = profile_weight(m, seed2) * E_seed2
    s = w1 + w2
    if s > 0:
        fracs[i, 0] = w1 / s
        fracs[i, 1] = w2 / s

# resulting cluster energies
E_c1 = float(np.sum(total * fracs[:, 0]))
E_c2 = float(np.sum(total * fracs[:, 1]))

# panels: (a) input total energy, (b) split assignment
fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6.5))

# Panel A
emax = total.max()
for (m, _, _), E in zip(region, total):
    if E < THR:
        draw_module(axA, m, facecolor='#f8f9fa')
        continue
    intensity = max(0.1, min(0.95, E / emax))
    fc = plt.cm.YlOrRd(intensity)
    draw_module(axA, m, facecolor=fc, edgecolor='#666', lw=0.5)
    axA.text(m['x'], m['y'], f"{E:.0f}", ha='center', va='center',
             fontsize=7.5, color='#222')

# mark the two maxima
for s, lbl, col in [(seed1, '★1', '#d62728'), (seed2, '★2', '#1f77b4')]:
    axA.plot(s['x'], s['y'], '*', color=col, ms=22, mec='white', mew=1.3,
             zorder=10)
    axA.text(s['x'], s['y'] + s['sy'] * 0.65, lbl, ha='center', va='center',
             fontsize=11, color=col, fontweight='bold')
axA.plot(x1, y1, 'x', color=col, ms=10, mew=1.8, zorder=11)

xlims = (center['x'] - center['sx'] * 4, center['x'] + center['sx'] * 4)
ylims = (center['y'] - center['sy'] * 4, center['y'] + center['sy'] * 4)
axA.set_xlim(xlims); axA.set_ylim(ylims)
axA.set_aspect('equal')
axA.set_xlabel('x (mm)'); axA.set_ylabel('y (mm)')
axA.set_title(f"Input — total energy in MeV\n"
              f"two local maxima found ({E_seed1:.0f}, {E_seed2:.0f} MeV)")

# Panel B — colour-coded by which cluster a module is mostly assigned to
# (alpha proportional to dominant fraction).
for i, (m, _, _) in enumerate(region):
    if total[i] < THR:
        draw_module(axB, m, facecolor='#f8f9fa')
        continue
    f1, f2 = fracs[i]
    if f1 > f2:
        col = '#d62728'; frac = f1
    else:
        col = '#1f77b4'; frac = f2
    # blend with white based on dominance (1.0 → solid, 0.5 → faint)
    alpha = 0.25 + 0.7 * (frac - 0.5) / 0.5
    alpha = max(0.10, min(0.95, alpha))
    draw_module(axB, m, facecolor=col, alpha=alpha, edgecolor='#666',
                lw=0.5)
    # show split fractions for the boundary modules
    if 0.10 < f1 < 0.90:
        axB.text(m['x'], m['y'],
                 f"{f1*100:.0f}/{f2*100:.0f}",
                 ha='center', va='center', fontsize=7, color='#000',
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='none',
                           pad=0.5))

for s, lbl, col in [(seed1, '★1', '#d62728'), (seed2, '★2', '#1f77b4')]:
    axB.plot(s['x'], s['y'], '*', color=col, ms=22, mec='white', mew=1.3,
             zorder=10)
axB.set_xlim(xlims); axB.set_ylim(ylims)
axB.set_aspect('equal')
axB.set_xlabel('x (mm)'); axB.set_ylabel('y (mm)')
axB.set_title(
    f"Profile split — fractions in % shown on boundary modules\n"
    f"E(★1) = {E_c1:.0f} MeV   E(★2) = {E_c2:.0f} MeV   (input {E1:.0f}, {E2:.0f})")

fig.tight_layout()
fig.savefig(HERE / 'hycal_fig3_split.png', dpi=130)
plt.close(fig)

print(f"\nTwo-shower split:")
print(f"  injected E1 = {E1:.0f}, recovered E(seed1) = {E_c1:.0f} MeV"
      f" ({E_c1/E1*100-100:+.1f}%)")
print(f"  injected E2 = {E2:.0f}, recovered E(seed2) = {E_c2:.0f} MeV"
      f" ({E_c2/E2*100-100:+.1f}%)")


# ---------------------------------------------------------------------------
# Plot 4 — parameter sensitivity: log_weight_thres + shower depth
# ---------------------------------------------------------------------------
# Reuse the single-cluster shower data
region_pos = [(m['x'], m['y']) for m, _, _ in nbrs]

# scan T over the typical range
t_values = np.linspace(2.0, 6.0, 9)
positions = []
for T in t_values:
    res = log_weighted_centroid(pos_3x3, e_3x3, E_total, thres=T)
    positions.append(res)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.6))

# Panel A: log_weight_thres effect
errs = [np.hypot(p[0] - x0, p[1] - y0) for p in positions]
axA.plot(t_values, errs, 'o-', color='#1f77b4', ms=6)
axA.axvline(3.6, color='#d62728', lw=0.8, ls='--',
            label='default T = 3.6')
axA.axhline(0, color='#888', lw=0.5)
axA.set_xlabel('log_weight_thres  T')
axA.set_ylabel('|reconstructed − true|  (mm)')
axA.set_title("Position error vs log-weight threshold\n"
              "(single 1.1 GeV PbWO₄ shower from fig 2)")
axA.grid(alpha=0.3)
axA.legend(fontsize=9)
# annotate energy-weighted equivalent (T → ∞ limit)
err_e = np.hypot(xe - x0, ye - y0)
axA.axhline(err_e, color='#2ca02c', lw=0.8, ls=':',
            label='energy-weighted error')
axA.text(t_values[-1], err_e + 0.08, 'energy-weighted',
         color='#2ca02c', fontsize=8, ha='right')

# Panel B: shower depth curves
def shower_depth(E_mev, kind):
    if E_mev <= 0:
        return 0.0
    if kind == 'PbWO4':
        return 8.6 * (np.log(E_mev / 1.1) - 0.5)
    return 26.7 * (np.log(E_mev / 2.84) - 0.5)

E_axis = np.logspace(np.log10(20), np.log10(3000), 200)
axB.plot(E_axis, [shower_depth(e, 'PbWO4')   for e in E_axis],
         color='#1f77b4', lw=1.6,
         label='PbWO₄  (X₀ = 8.6 mm,  Eᶜ = 1.1 MeV)')
axB.plot(E_axis, [shower_depth(e, 'PbGlass') for e in E_axis],
         color='#d62728', lw=1.6,
         label='PbGlass (X₀ = 26.7 mm, Eᶜ = 2.84 MeV)')
axB.set_xscale('log')
axB.set_xlabel('cluster energy  E (MeV)')
axB.set_ylabel('shower-max depth  t (mm)')
axB.set_title("Shower depth  t = X₀ · (ln(E/Eᶜ) − 0.5)\n"
              "added to z_HyCal-face for cl_z")
axB.legend(loc='lower right', fontsize=8.5)
axB.grid(alpha=0.3, which='both')

fig.tight_layout()
fig.savefig(HERE / 'hycal_fig4_params.png', dpi=130)
plt.close(fig)

print(f"\nWrote 4 PNGs to {HERE}")
