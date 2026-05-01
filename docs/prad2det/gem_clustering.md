# GEM Clustering in `prad2det`

`gem::GemCluster` (in [`prad2det/include/GemCluster.h`](../../prad2det/include/GemCluster.h))
takes per-plane strip charges from the SSP/MPD readout and produces 2-D
GEM hits. It runs in two stages: a **per-plane 1-D clustering**
pipeline (group consecutive strips → recursively split at charge
valleys → charge-weighted position → cross-talk filter), followed by
an **X/Y matching** pass that turns paired X- and Y-clusters into
3-vector `GEMHit`s.

The algorithm is ported from `mpd_gem_view_ssp::GEMCluster` (Bai /
Gnanvo / Peng) and is invoked once per detector per event from
`gem::GemSystem::Reconstruct()`.

## GEM detector layout

PRad-II uses **four GEM detectors** in **two paired layers** for
redundant 2-D tracking and ghost-hit rejection. Each detector reads
out two orthogonal strip planes through APV25 chips at 0.4 mm pitch:

![layout](gem_fig1_layout.png)

| field | value |
|---|---|
| detectors | GEM0, GEM1 (layer 1, z ≈ 5407 mm); GEM2, GEM3 (layer 2, z ≈ 5807 mm) |
| inter-detector spacing | 39.7 mm (within a layer) |
| X plane | 12 APVs × 128 ch, pitch 0.4 mm — 1 408 strips ≈ 563.2 mm. APV positions 10/11 share strips around the hole (`shared_pos: 10`, `pin_rotate: 16`) |
| Y plane | 24 APVs × 128 ch, pitch 0.4 mm — 3 072 strips ≈ 1 228.8 mm. Strips that cross the hole y-band are split into top + bottom segments |
| beam hole | 52 × 52 mm², centred vertically (y = 614.4 mm) and offset along x (x = 534.4 mm) — sits inside the APV pos 10/11 strip range |

Strip indices come from `GemSystem::ProcessEvent()` after pedestal
subtraction, common-mode correction, and zero-suppression. Each
`StripHit` carries the plane-wise strip number, the maximum
time-sample charge (`charge`), the time bin where that maximum lives
(`max_timebin`), the physical position in mm, and the full
time-sample ADC vector (used downstream for time-difference cuts).

## Algorithm

### Step 1 — Group consecutive strips

`groupHits()` sorts the input `StripHit`s by strip number, then walks
the sorted list and starts a new group whenever the gap to the next
strip exceeds `consecutive_thres` (default 1, i.e. only strictly
adjacent strips group). Each group is a candidate 1-D cluster but may
still need to be split if multiple showers share neighbouring strips.

### Step 2 — Recursive valley split

For each group `splitCluster()` looks for an internal local minimum
preceded by a sufficiently steep descent and followed by a
sufficiently steep ascent — both gated by `split_thres` (default 14
ADC counts):

1. Walk the group left-to-right. Set `descending = true` as soon as
   `charge[i] − charge[i+1] > split_thres`.
2. Once descending, track the running minimum.
3. The first ascent step where `charge[i+1] − charge[i] > split_thres`
   confirms the valley.
4. Halve the charge of the valley strip (it is shared between the
   left and right sub-clusters), emit the left sub-cluster, and
   recurse on the right.

Groups smaller than 3 strips never split. Groups with multiple
valleys split recursively, producing one cluster per "peak".

### Step 3 — Charge-weighted position

Each cluster's position comes from `reconstructCluster()`:

```
position    = Σ (strip_pos_i · charge_i) / Σ charge_i      [mm]
peak_charge = max_i charge_i
total_charge = Σ charge_i
max_timebin  = time bin of the seed (highest-charge) strip
```

This is a plain centroid, not log-weighted — at 0.4 mm pitch the
strip density is high enough that linear weighting tracks the shower
position well, and unlike HyCal there's no module-scale grid to
worry about.

### Step 4 — Cross-talk identification

APV25 capacitive coupling and electronic cross-talk produce **ghost
clusters** at characteristic distances from a true cluster (a known
function of the chip's pin-out and sampling pattern). `setCrossTalk()`
sorts clusters by ascending `peak_charge` and walks the list: any
cluster all of whose strips are flagged `cross_talk` (set upstream by
`GemSystem` based on inter-strip ADC ratios) and whose distance to a
larger cluster matches one of the `charac_dists` entries within
±`cross_talk_width` mm gets its own `cross_talk` flag set.

The default characteristic distances mirror the original
`mpd_gem_view_ssp` values:

```cpp
charac_dists = { 6.4, 17.6, 24.4, 24.8, 25.2, 25.6,
                 26.0, 26.4, 26.8, 33.6, 44.8 };   // mm
```

### Step 5 — Filter

`filterClusters()` drops:

- clusters with fewer than `min_cluster_hits` (default 1) strips,
- clusters with more than `max_cluster_hits` (default 20) strips
  (these are typically noise bursts, not real showers),
- clusters flagged as cross-talk.

### Step 6 — Cartesian X/Y matching

`CartesianReconstruct()` turns a list of accepted X-clusters and a
list of accepted Y-clusters into 2-D `GEMHit`s. Two modes are
supported:

**Mode 0 — ADC-sorted 1:1.** Sort both lists by descending
`peak_charge` and pair index-wise. Always produces
`min(N_x, N_y)` hits with no rejected combinations. Useful when the
upstream strip clustering is clean enough that the brightest X
genuinely matches the brightest Y.

**Mode 1 — full Cartesian + cuts (default).** Form every X×Y pair,
then drop any pair that fails either of:

- **ADC asymmetry** — `|Q_X_peak − Q_Y_peak| / (Q_X_peak + Q_Y_peak) ≤ match_adc_asymmetry` (default 0.8). Real GEM hits deposit similar charge on both planes; large asymmetries flag ghosts from accidental coincidence.
- **Time difference** — the seed-strip ADC-weighted mean times must satisfy `|⟨t⟩_X − ⟨t⟩_Y| ≤ match_time_diff` (default 50 ns). The seed mean time is `Σ(adc_i · t_i) / Σ adc_i` over the time samples, where `t_i = (i + 1) · ts_period` ns.

Any pair passing both cuts becomes a `GEMHit` with `det_id` and
per-plane charge / size / max-timebin recorded for downstream use.

## Parameters

All settings live in `gem::ClusterConfig` (in
[`prad2det/include/GemSystem.h`](../../prad2det/include/GemSystem.h)).
The system stores one `ClusterConfig` per detector via
`SetReconConfigs()`, so different detectors can have different cuts —
useful when the four GEMs have different APV gains or noise levels.

| field | default | unit | role |
|---|---:|---|---|
| `min_cluster_hits` | 1 | strips | Lower bound on cluster size — drops single-noisy-strip "clusters" if set ≥ 2. |
| `max_cluster_hits` | 20 | strips | Upper bound on cluster size — kills runaway clusters caused by noise bursts or beam halo. |
| `consecutive_thres` | 1 | strips | Maximum gap between adjacent strips that still keep them in the same group. `1` = strictly adjacent only. |
| `split_thres` | 14 | ADC counts | Charge-difference threshold that gates both descent detection and valley-confirmation in `splitCluster()`. Lower values split more aggressively. |
| `cross_talk_width` | 2 | mm | Tolerance on the characteristic-distance match in `setCrossTalk()`. |
| `charac_dists` | `{6.4, 17.6, 24.4..26.8, 33.6, 44.8}` | mm | APV25 cross-talk characteristic distances. |
| `match_mode` | 1 | — | `0` = ADC-sorted 1:1 matching, `1` = Cartesian product with cuts. |
| `match_adc_asymmetry` | 0.8 | fraction | Cap on `|Q_X − Q_Y|/(Q_X + Q_Y)` (mode 1). Set negative to disable. |
| `match_time_diff` | 50 | ns | Cap on `|⟨t⟩_X − ⟨t⟩_Y|` (mode 1). Set negative to disable. |
| `ts_period` | 25 | ns | Time-sample period (default = 1 / 40 MHz APV clock). |

## Worked example — strip clustering

A 75-strip window with three real showers and ~4 ADC noise:

![strip](gem_fig2_strip_clustering.png)

The left panel shows the full above-threshold strip distribution
(threshold 30 ADC). DFS-style grouping with `consecutive_thres = 1`
produces three groups (110-114, 130-139, 158-160). The middle group
has two local maxima — `splitCluster()` walks from the 132 peak,
detects a descent into the 134/135 valley, then sees a 198 ADC upturn
into the 137 peak (»`split_thres = 14`) and partitions the group at
strip 135 (right panel). The valley strip's charge is halved and
contributes to both sub-clusters, so neither side double-counts the
shared edge.

| cluster | strips | position (mm) | Σ ADC | peak ADC |
|---|:---:|---:|---:|---:|
| 1 | 110-114 | 44.87 | 2143 | 778 |
| 2 (left of split) | 130-134 | 52.81 | 1750 | 599 |
| 2 (right of split) | 135-139 | 54.82 | 1182 | 480 |
| 3 | 158-160 | 63.60 |  566 | 297 |

Note that the position (charge-weighted centroid) lands between
strips at sub-pitch resolution — for cluster 1, x = 44.87 mm sits
between strips 112 (44.8 mm) and 113 (45.2 mm), reflecting the actual
shower offset within the strip pitch.

## Worked example — X/Y matching

Three X-plane clusters and three Y-plane clusters per detector, with
two prompt big showers and one late, small "out-of-time" pair (e.g.
backsplash or a delayed accidental):

![matching](gem_fig3_xy_matching.png)

**Left panel — Mode 1 (default).** All 9 X×Y candidates are listed.
The big-prompt × big-prompt pairings (X0/X1 with Y0/Y1) easily pass
both cuts. The big × late pairings fail the time cut (`Δt > 50 ns`).
Note that the small × small (X2 ↔ Y2) pair *does* pass — both have
similar (small) charge and similar (late) timing, so it gets
reconstructed as a 2-D hit even though it is most likely noise. If
that's a problem in production, tighter `match_adc_asymmetry` plus a
minimum on `peak_charge` upstream filters it out.

**Right panel — Mode 0.** No physical cuts; just sort both lists by
peak ADC and pair X[rank] ↔ Y[rank]. Three hits, by construction, but
with no defence against ghost-pair formation when accidental
coincidence is significant. Mode 0 is mostly useful for debugging or
for very low-occupancy runs where the cuts of mode 1 would just throw
away good hits.

## Parameter sensitivity

![params](gem_fig4_params.png)

**Left — `split_thres`.** Same multi-peak group, three different
thresholds:

- `split_thres = 200`: too coarse — the two real peaks aren't
  resolved (1 cluster).
- `split_thres = 50`: catches the deep valley between the two main
  peaks (2 clusters), correct for this trace.
- `split_thres = 5`: also splits a small secondary fluctuation (3
  clusters), some of which are spurious.

The default of 14 ADC is calibrated for the typical noise RMS (~5
ADC) plus a margin to avoid splitting on noise; lower it on quieter
detectors, raise it on noisy ones.

**Right — cross-talk.** A bright primary cluster at 50 mm with a
small ghost cluster 24.4 mm away. The dotted lines show all 11
characteristic distances — the densely-clustered 24.4–26.8 mm group
covers the most common APV25 cross-talk pattern. `setCrossTalk()`
identifies the small cluster as a cross-talk match (it sits within
`cross_talk_width = 2 mm` of one of the characteristic distances and
has only `cross_talk`-flagged strips), and `filterClusters()` drops
it. The primary survives.

## Output — `GEMHit`

Each X/Y match produces one `GEMHit`:

| field | type | meaning |
|---|---|---|
| `x`, `y`, `z` | `float` | Hit position (mm). `z` is set by the application from per-detector geometry; `GemCluster` itself sets `z = 0`. |
| `det_id` | `int` | 0..3 (GEM0..GEM3). |
| `x_charge`, `y_charge` | `float` | Total ADC of the X / Y cluster. |
| `x_peak`, `y_peak` | `float` | Max-strip ADC of the X / Y cluster. |
| `x_max_timebin`, `y_max_timebin` | `short` | Time-sample bin of the max-ADC strip on each plane. |
| `x_size`, `y_size` | `int` | Number of strips in the X / Y cluster. |

These map directly onto the `gem_*` branches of the recon tree (see
[`docs/REPLAYED_DATA.md`](../REPLAYED_DATA.md)).

## Reproducing the plots

The detector geometry is read from
[`database/gem_daq_map.json`](../../database/gem_daq_map.json); the
strip-clustering and matching algorithms are re-implemented in pure
Python in [`plot_gem_clustering.py`](plot_gem_clustering.py) (NumPy +
Matplotlib only).

```bash
cd docs/prad2det
python plot_gem_clustering.py
```

Regenerates `gem_fig1_layout.png`, `gem_fig2_strip_clustering.png`,
`gem_fig3_xy_matching.png`, `gem_fig4_params.png` and prints the
reconstructed cluster table to stdout.

## See also

- [`prad2det/include/GemCluster.h`](../../prad2det/include/GemCluster.h),
  [`GemCluster.cpp`](../../prad2det/src/GemCluster.cpp) — algorithm source
- [`prad2det/include/GemSystem.h`](../../prad2det/include/GemSystem.h) —
  hierarchy, pedestal/CM/zero-suppression, strip mapping, per-detector
  `ClusterConfig` storage
- [`database/gem_daq_map.json`](../../database/gem_daq_map.json) —
  APV mapping + plane / pitch / hole geometry
- [`database/reconstruction_config.json`](../../database/reconstruction_config.json) —
  per-run cluster-config defaults
- [`docs/REPLAYED_DATA.md`](../REPLAYED_DATA.md) —
  branch layout for the recon tree (where `GEMHit`s land as `gem_*`)
- mpd_gem_view_ssp `GEMCluster` — original implementation lineage
