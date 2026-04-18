#!/usr/bin/env python3
"""
Visualize PRad-II GEM strip layout from gem_map.json.

Shows:
- X strips (vertical lines) in blue — shortened near beam hole for split APVs
- Y strips (horizontal lines) in red — split at beam hole boundary
- APV boundaries as dashed lines
- Beam hole as a yellow rectangle
- One subplot per detector (GEM0-GEM3)

Usage:
    python gem_layout.py [path/to/gem_map.json]

Defaults to database/gem_map.json if no argument given.
"""

import argparse
import json
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from gem_strip_map import map_apv_strips


def load_gem_map(path):
    with open(path, encoding="utf-8") as f:
        raw = json.load(f)

    layers = {l["id"]: l for l in raw["layers"]}
    apvs = [e for e in raw["apvs"] if "crate" in e]
    hole = raw.get("hole", None)

    return layers, apvs, hole, raw


def build_strip_layout(layers, apvs, hole, raw):
    """Build per-detector strip positions, accounting for beam hole.

    X strips near the hole are shortened (split APVs with 16 disconnected channels).
    Y strips are split at the hole Y boundary into top/bottom segments.
    """
    detectors = {}
    strips_per_apv = raw.get("apv_channels", 128)

    for det_id, layer in layers.items():
        x_pitch = layer["x_pitch"]
        y_pitch = layer["y_pitch"]
        y_size = layer["y_apvs"] * strips_per_apv * y_pitch

        detectors[det_id] = {
            "name": layer["name"],
            "x_size": 0,       # computed from actual strip positions below
            "y_size": y_size,
            "x_pitch": x_pitch,
            "y_pitch": y_pitch,
            "x_strips": [],
            "y_strips": [],
            "x_apv_edges": set(),
            "y_apv_edges": set(),
        }

    apv_ch = raw.get("apv_channels", 128)
    ro_center = raw.get("readout_center", 32)

    # first pass: compute strip positions and derive x_size + hole position
    all_x_strips = {det_id: set() for det_id in detectors}
    match_strips = {det_id: [] for det_id in detectors}  # strips from match APVs
    apv_data = []

    for apv in apvs:
        det_id = apv["det"]
        if det_id not in detectors:
            continue
        plane = apv["plane"]
        plane_strips = map_apv_strips(apv, apv_channels=apv_ch, readout_center=ro_center)
        apv_data.append((apv, det_id, plane, plane_strips))

        if plane == "X":
            all_x_strips[det_id].update(plane_strips)
            if apv.get("match", ""):
                match_strips[det_id].extend(plane_strips)

    # set x_size from actual strip range
    for det_id, strips in all_x_strips.items():
        if strips:
            det = detectors[det_id]
            x_max_strip = max(strips)
            det["x_size"] = (x_max_strip + 1) * det["x_pitch"]

    # derive hole position from match APV strips and detector geometry
    # x_center: center of the match strips range (where pos 10/11 overlap)
    # y_center: center of detector Y
    if hole:
        hw = hole["width"]
        hh = hole["height"]
    else:
        hw = hh = 0

    # compute per-detector (use first detector as reference since all identical)
    ref_det_id = min(detectors.keys())
    ref_det = detectors[ref_det_id]
    if hole and match_strips[ref_det_id]:
        ms = match_strips[ref_det_id]
        hx = (min(ms) + max(ms) + 1) / 2 * ref_det["x_pitch"]
        hy = ref_det["y_size"] / 2
        hole_x0, hole_x1 = hx - hw / 2, hx + hw / 2
        hole_y0, hole_y1 = hy - hh / 2, hy + hh / 2
        # update hole dict for display
        hole["x_center"] = hx
        hole["y_center"] = hy
    else:
        hole_x0 = hole_x1 = hole_y0 = hole_y1 = -1

    # second pass: build strip lines
    for apv, det_id, plane, plane_strips in apv_data:
        det = detectors[det_id]
        match = apv.get("match", "")

        if plane == "X":
            pitch = det["x_pitch"]
            strip_positions = sorted(set(plane_strips))
            x_min = min(strip_positions) * pitch
            x_max = (max(strip_positions) + 1) * pitch
            if match == "+Y" and hole:
                y0_edge, y1_edge = hole_y1, det["y_size"]
            elif match == "-Y" and hole:
                y0_edge, y1_edge = 0, hole_y0
            else:
                y0_edge, y1_edge = 0, det["y_size"]

            det["x_apv_edges"].add((x_min, y0_edge, y1_edge))
            det["x_apv_edges"].add((x_max, y0_edge, y1_edge))

            for s in plane_strips:
                strip_x = s * pitch
                det["x_strips"].append((strip_x, y0_edge, y1_edge))

        elif plane == "Y":
            pitch = det["y_pitch"]
            strip_positions = sorted(set(plane_strips))
            y_min = min(strip_positions) * pitch
            y_max = (max(strip_positions) + 1) * pitch
            det["y_apv_edges"].add(y_min)
            det["y_apv_edges"].add(y_max)

            for s in plane_strips:
                strip_y = s * pitch

                # Y strips split at hole boundary (strictly inside hole)
                if hole and hole_y0 < strip_y < hole_y1:
                    det["y_strips"].append((strip_y, 0, hole_x0))
                    det["y_strips"].append((strip_y, hole_x1, det["x_size"]))
                else:
                    det["y_strips"].append((strip_y, 0, det["x_size"]))

    return detectors


def plot_detector(ax, det, det_id, hole, show_every=8):
    """Plot one GEM detector's strip layout."""
    name = det["name"]
    x_size = det["x_size"]
    y_size = det["y_size"]

    # detector outline
    ax.add_patch(plt.Rectangle((0, 0), x_size, y_size,
                                fill=False, edgecolor="gray", linewidth=1.5))

    # beam hole
    if hole:
        hx = hole["x_center"]
        hy = hole["y_center"]
        hw = hole["width"]
        hh = hole["height"]
        ax.add_patch(plt.Rectangle((hx - hw/2, hy - hh/2), hw, hh,
                                    fill=True, facecolor="#44442200",
                                    edgecolor="#ffcc00", linewidth=2, linestyle="-",
                                    zorder=5))

    # X strips (vertical lines) — blue
    # group by Y extent to preserve show_every sampling within each group
    x_by_extent = {}
    for (x, y0, y1) in det["x_strips"]:
        key = (y0, y1)
        x_by_extent.setdefault(key, []).append((x, y0, y1))
    x_lines = []
    for key in sorted(x_by_extent):
        group = sorted(x_by_extent[key])
        for i, (x, y0, y1) in enumerate(group):
            if i % show_every == 0:
                x_lines.append([(x, y0), (x, y1)])
    if x_lines:
        ax.add_collection(LineCollection(x_lines, colors="steelblue",
                                          linewidths=0.3, alpha=0.6))

    # Y strips (horizontal lines) — red
    y_by_extent = {}
    for (y, x0, x1) in det["y_strips"]:
        key = (x0, x1)
        y_by_extent.setdefault(key, []).append((y, x0, x1))
    y_lines = []
    for key in sorted(y_by_extent):
        group = sorted(y_by_extent[key])
        for i, (y, x0, x1) in enumerate(group):
            if i % show_every == 0:
                y_lines.append([(x0, y), (x1, y)])
    if y_lines:
        ax.add_collection(LineCollection(y_lines, colors="indianred",
                                          linewidths=0.3, alpha=0.6))

    # APV boundaries — dashed lines
    for (x, y0, y1) in sorted(det["x_apv_edges"]):
        ax.plot([x, x], [y0, y1], color="steelblue", linewidth=0.8, alpha=0.5, linestyle="--")
    for y in sorted(det["y_apv_edges"]):
        ax.axhline(y, color="indianred", linewidth=0.8, alpha=0.5, linestyle="--")

    ax.set_xlim(-x_size * 0.05, x_size * 1.05)
    ax.set_ylim(-y_size * 0.05, y_size * 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")

    n_x = len(det["x_strips"])
    n_y = len([s for s in det["y_strips"]])
    ax.set_title(f"{name} — {n_x} X strips, {n_y} Y strip segments")

    handles = [
        mpatches.Patch(color="steelblue", alpha=0.6, label=f"X strips ({n_x})"),
        mpatches.Patch(color="indianred", alpha=0.6, label=f"Y strips"),
    ]
    if hole:
        handles.append(mpatches.Patch(facecolor="#ffcc0044", edgecolor="#ffcc00",
                                       label="Beam hole"))
    ax.legend(handles=handles, loc="upper right", fontsize=8)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize GEM strip layout from gem_map.json.")
    parser.add_argument("-G", "--gem-map", default=None,
                        help="GEM map JSON path (default: auto-search "
                             "common locations).")
    args = parser.parse_args()

    gem_map_path = args.gem_map
    if not gem_map_path:
        for candidate in [
            "database/gem_map.json",
            "../database/gem_map.json",
            "../../database/gem_map.json",
            "gem_map.json",
        ]:
            if os.path.exists(candidate):
                gem_map_path = candidate
                break
        else:
            print("Error: cannot find gem_map.json (pass -G <path>)")
            sys.exit(1)

    print(f"Loading: {gem_map_path}")
    layers, apvs, hole, raw = load_gem_map(gem_map_path)
    detectors = build_strip_layout(layers, apvs, hole, raw)

    if hole:
        print(f"Beam hole: {hole['width']}x{hole['height']} mm "
              f"at ({hole['x_center']}, {hole['y_center']})")

    print(f"Detectors: {len(detectors)}")
    for det_id, det in sorted(detectors.items()):
        print(f"  {det['name']}: {det['x_size']:.1f} x {det['y_size']:.1f} mm, "
              f"{len(det['x_strips'])} X strips, {len(det['y_strips'])} Y strip segments")

    # all 4 GEMs are identical — show only the first one
    det_id = min(detectors.keys())
    det = detectors[det_id]

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    fig.suptitle(f"PRad-II GEM Strip Layout ({det['name']})", fontsize=14)
    plot_detector(ax, det, det_id, hole)

    plt.tight_layout()
    plt.savefig("gem_layout.png", dpi=150, bbox_inches="tight")
    print("Saved: gem_layout.png")
    plt.show()


if __name__ == "__main__":
    main()
