# HyCal Clustering in `prad2det`

`fdec::HyCalCluster` (in [`prad2det/include/HyCalCluster.h`](../../prad2det/include/HyCalCluster.h))
implements the **Island clustering** algorithm originally written for
PRad-I (`PRadIslandCluster` / `PRadHyCalReconstructor`), now ported to
operate on the index-based `HyCalSystem` geometry. It groups adjacent
above-threshold modules, splits multi-maximum islands using a transverse
shower profile, and reconstructs each cluster's position with a
log-weighted centroid.

The algorithm is invoked once per event. Inputs come from the
calorimeter's per-module calibrated energies; outputs are a list of
`ClusterHit` records (centre module, x/y position at the HyCal face,
energy, block count, flags).

## HyCal geometry

HyCal is a hybrid sandwich: an inner ~700 × 700 mm² PbWO₄ matrix wrapped
by an outer PbGlass ring. Modules sit on a pixel-perfect grid within
each of five sectors (Center, Top, Right, Bottom, Left); the sector
boundary is a known geometric step in module size and is handled by the
`HyCalSystem::qdist()` quantised-distance helper.

![layout](figs/hycal_fig1_layout.png)

| sector | type | module size | count |
|---|---|---:|---:|
| Center | PbWO₄ | 20.77 × 20.75 mm | 1152 |
| Top / Right / Bottom / Left | PbGlass | 38.15 × 38.15 mm | 576 (total) |

`HyCalSystem` builds an O(1) per-sector grid + a pre-computed
cross-sector neighbour list at `Init()` time, so the per-event hot path
never touches the JSON or recomputes geometry.

## Algorithm

### Step 1 — DFS grouping

Each above-threshold hit (`energy > min_module_energy`) is loaded with
`AddHit(module_index, energy)`. The grouping pass walks every hit and,
via depth-first search through `HyCalSystem::for_each_neighbor`, builds
**connected components** — a "group" is a maximal set of hits where each
member shares an edge (or, optionally, a corner) with another member.

`corner_conn = false` (the default) means only edge-sharing modules
join. Setting it to `true` adds diagonal neighbours, which can be
useful for very narrow showers but tends to merge accidentally close
clusters.

The neighbour iteration uses the row/col grid for same-sector
adjacency (O(8)) and falls through to a small precomputed list across
sector boundaries (typically 0–3 entries per module).

### Step 2 — Find local maxima

Within each group, `find_maxima()` walks the hits and keeps every
module whose energy is

- ≥ `min_center_energy` (rejects soft hits that can't seat a cluster), and
- strictly greater than every grid neighbour (corners always included
  for the maximality test, regardless of `corner_conn`).

Single-maximum groups become a single cluster; multi-maximum groups go
through Step 3.

### Step 3 — Split

When a group has more than one local maximum, hits that lie between
the maxima must be partitioned. `split_hits()` does this iteratively:

1. **Initial fractions.** For each maximum *i* and each hit *j* in the
   group, set `frac[j][i] = profile(d_ij) · E_max[i]` where `profile` is
   a transverse-shower-shape lookup (see `IClusterProfile`) and `d_ij`
   is the quantised module distance.
2. **Refine.** `eval_fraction()` runs `split_iter` passes (default 6).
   Each pass uses the current normalised fractions to reconstruct a
   provisional centroid for every maximum, then re-evaluates
   `frac[j][i] = profile(|hit − centroid_i|) · E_total_i` against that
   centroid. The 3×3 neighbourhood around each centroid is what feeds
   the position reconstruction inside the loop.
3. **Emit.** Hits whose normalised fraction `< least_split` (default
   0.01) are dropped from that cluster; surviving hits are added with
   their fraction-weighted energy. Modules shared by two clusters end
   up split between them (clusters get `kSplit` set in their flag).

### Step 4 — Position reconstruction

For every accepted cluster, `reconstruct_pos()` recomputes the centroid
once more, this time using the **log-weighted** scheme:

```
w_i = max(0, log_weight_thres + ln(E_i / E_total))
x = x_seed + (Σ w_i · dx_i) / Σ w_i  · seed.size_x
y = y_seed + (Σ w_i · dy_i) / Σ w_i  · seed.size_y
```

Only hits within the seed's 3×3 grid neighbourhood contribute (≤ 9
entries; the seed itself plus 8 neighbours, capped at
`POS_RECON_HITS = 15`). The threshold `log_weight_thres = 3.6` zeros
out modules that carry less than `e^(−3.6) ≈ 2.7 %` of the cluster
energy, which suppresses long shower tails and noise-driven hits while
still letting the dominant 3×3 modules drive the position.

The output `ClusterHit` carries the centre-module's PrimEx ID, the
reconstructed (x, y) at the HyCal face (mm, lab frame), total energy,
number of blocks, and a flag bitmask. Downstream code adds a
shower-depth z offset (see [Parameters → Shower depth](#shower-depth))
to get the 3-D cluster position.

## Parameters

All settings live in `fdec::ClusterConfig`. Defaults match
[`database/reconstruction_config.json`](../../database/reconstruction_config.json).

| field | default | unit | role |
|---|---:|---|---|
| `min_module_energy` | 1.0 | MeV | Hit threshold — modules below this are dropped before grouping. |
| `min_center_energy` | 10.0 | MeV | Seed threshold — a module needs at least this energy to be a local maximum. |
| `min_cluster_energy` | 50.0 | MeV | Cluster acceptance — clusters totaling less than this are dropped at `ReconstructHits()`. |
| `min_cluster_size` | 1 | modules | Cluster acceptance — minimum block count. |
| `corner_conn` | `false` | — | When `true`, diagonal neighbours join during DFS grouping. |
| `split_iter` | 6 | iterations | Passes through `eval_fraction()` for multi-maximum groups. |
| `least_split` | 0.01 | fraction | Modules whose normalised split fraction falls below this are dropped from the secondary cluster. |
| `log_weight_thres` | 3.6 | — | `T` in `w = max(0, T + ln(E_i / E_total))`. Higher values include more low-energy modules; lower values pull the centroid toward the seed. |

## Worked example — single cluster

A 1.1 GeV photon shower placed at (57, 89) mm, σ_shower ≈ 15 mm
(Molière radius for PbWO₄ is ~20 mm, so most of the energy lives in
the central 3×3).

![single](figs/hycal_fig2_single_cluster.png)

| quantity | value |
|---|---:|
| seed module (W428) | (51.92, 93.38) mm |
| 3×3 sum / cluster total | 1097 / 1097 MeV |
| true position | (56.92, 89.38) mm |
| energy-weighted | (56.04, 90.07) mm — 1.13 mm error |
| log-weighted (T = 3.6) | (56.88, 89.60) mm — 0.23 mm error |

Two things to notice. First, the seed module is *not* the true shower
position — that's expected: the seed always lands on a module centre,
and the actual shower is offset within the module. Position
reconstruction's job is to recover the offset from neighbours' energy
sharing. Second, the log-weighted centroid is ~5× closer to the truth
than the energy-weighted centroid here, because the linear weighting
gives non-trivial weight to modules carrying < 5 % of the energy and
their position uncertainty pulls the reconstructed point toward the
sampling cell centres.

## Worked example — two-shower split

Two showers (1500 MeV at (−25, 5) mm, 900 MeV at (25, −8) mm) drop
into the same connected island. Single DFS pass returns one group of
~25 modules with two local maxima.

![split](figs/hycal_fig3_split.png)

The left panel shows the input total energy (per module, MeV); the two
local maxima — at the seed candidates ★1 / ★2 — are clearly separated
by ~3 modules. The right panel shows each module's dominant cluster
(red = ★1, blue = ★2). Modules along the shared boundary (faint pink
strip) carry mixed contributions, with their split percentages printed
on top. Recovered cluster energies are 1479 / 918 MeV vs the injected
1500 / 900 MeV — a ~1–2 % bias from the simple analytical profile.

The `kSplit` flag is set on both resulting clusters so downstream code
can choose to trust them less or apply a leakage correction.

## Parameter sensitivity

![params](figs/hycal_fig4_params.png)

**Left — `log_weight_thres`.** Sweeping `T` over [2, 6] on the single
cluster from above shows a clear minimum near the default 3.6.
- `T < 3` cuts too many neighbours: the centroid collapses onto the
  seed cell centre and tracks the seed module's grid position rather
  than the true shower offset.
- `T > 4` admits too many low-energy modules: their weights are no
  longer dominated by the bright neighbours, and the centroid drifts
  toward the cell-centre average — eventually approaching the
  energy-weighted result (green dotted line).
- The default 3.6 lands at the bottom of the well by design (it was
  tuned against MC and beam-test data in the original PRad
  reconstruction).

**Right — shower depth.** The cluster's z-position uses
`shower_depth(center_id, energy)`:

```
t = X₀ · (ln(E / Eᶜ) − Cf)         Cf = 0.5 (photon)
PbWO₄  : X₀ = 8.6  mm,  Eᶜ = 1.1  MeV
PbGlass: X₀ = 26.7 mm,  Eᶜ = 2.84 MeV
```

PbGlass showers reach ~3× deeper into the calorimeter than PbWO₄ at the
same energy because the radiation length is longer. Replay code adds
this offset to `hycal_z` so `cl_z` reflects shower-max, not the front
face — important for matching to GEM tracks since the shower-max
depth is energy-dependent.

## Output — `ClusterHit`

`ReconstructHits()` returns one record per accepted cluster:

| field | type | meaning |
|---|---|---|
| `center_id` | `int` | PrimEx ID of the seed module (1..1156 for PbGlass G-modules, 1001..2152 for PbWO₄ W-modules + 1000) |
| `x`, `y` | `float` | Lab-frame mm at the HyCal face + shower-depth correction |
| `energy` | `float` | Total cluster energy (MeV) |
| `nblocks` | `int` | Modules contributing to this cluster (post-split) |
| `npos` | `int` | Modules used in the log-weighted position (≤ 9) |
| `flag` | `uint32_t` | Bitmask of layout + algorithm flags (see `LayoutFlag` in `HyCalSystem.h`) |

Useful flag bits (defined in `HyCalSystem.h`):

| bit | flag | meaning |
|---|---|---|
| 2 | `kTransition` | seed sits on the PbWO₄ ↔ PbGlass boundary |
| 3 | `kInnerBound` | seed touches the beam hole |
| 4 | `kOuterBound` | seed touches HyCal's outer edge |
| 5 | `kDeadModule` | seed is flagged dead in the geometry config |
| 6 | `kDeadNeighbor` | a neighbour is flagged dead — leakage correction may be needed |
| 7 | `kSplit` | cluster came out of `split_hits()` |
| 8 | `kLeakCorr` | leakage correction has been applied |

## Reproducing the plots

Both the geometry (loaded from `hycal_modules.json`) and the algorithm
illustrations live in
[`plot_hycal_clustering.py`](plot_hycal_clustering.py) (NumPy +
Matplotlib only). The script implements the seed-finding and
profile-based split in pure Python — close enough to the C++ for
illustration, not bit-exact.

```bash
cd docs/prad2det
python plot_hycal_clustering.py
```

Regenerates `figs/hycal_fig1_layout.png`, `figs/hycal_fig2_single_cluster.png`,
`figs/hycal_fig3_split.png`, `figs/hycal_fig4_params.png` and prints reconstructed positions
+ split fractions to stdout.

## See also

- [`prad2det/include/HyCalCluster.h`](../../prad2det/include/HyCalCluster.h),
  [`HyCalCluster.cpp`](../../prad2det/src/HyCalCluster.cpp) — algorithm source
- [`prad2det/include/HyCalSystem.h`](../../prad2det/include/HyCalSystem.h) —
  geometry, neighbour grids, sector helpers
- [`database/hycal_modules.json`](../../database/hycal_modules.json) —
  per-module geometry (used by the plot script)
- [`database/reconstruction_config.json`](../../database/reconstruction_config.json) —
  per-run cluster-config defaults
- [`docs/REPLAYED_DATA.md`](../REPLAYED_DATA.md) —
  branch layout for the recon tree (where `ClusterHit` lands as `cl_*`)
- PRad-I lineage: `PRadIslandCluster` / `PRadHyCalReconstructor` in
  PRadAnalyzer
