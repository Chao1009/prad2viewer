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

import json
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection


def load_gem_map(path):
    with open(path, encoding="utf-8") as f:
        raw = json.load(f)

    layers = {l["id"]: l for l in raw["layers"]}
    apvs = [e for e in raw["apvs"] if "crate" in e]
    hole = raw.get("hole", None)

    return layers, apvs, hole


def build_strip_layout(layers, apvs, hole):
    """Build per-detector strip positions, accounting for beam hole.

    X strips near the hole are shortened (split APVs with 16 disconnected channels).
    Y strips are split at the hole Y boundary into top/bottom segments.
    """
    detectors = {}
    strips_per_apv = 128

    # hole geometry (same for all detectors in local frame)
    if hole:
        hx = hole["x_center"]
        hy = hole["y_center"]
        hw = hole["width"]
        hh = hole["height"]
        hole_x0, hole_x1 = hx - hw / 2, hx + hw / 2
        hole_y0, hole_y1 = hy - hh / 2, hy + hh / 2
    else:
        hole_x0 = hole_x1 = hole_y0 = hole_y1 = -1

    for det_id, layer in layers.items():
        x_pitch = layer["x_pitch"]
        y_pitch = layer["y_pitch"]
        x_size = layer["x_apvs"] * strips_per_apv * x_pitch
        y_size = layer["y_apvs"] * strips_per_apv * y_pitch

        detectors[det_id] = {
            "name": layer["name"],
            "x_size": x_size,
            "y_size": y_size,
            "x_pitch": x_pitch,
            "y_pitch": y_pitch,
            "x_strips": [],
            "y_strips": [],
            "x_apv_edges": set(),
            "y_apv_edges": set(),
        }

    for apv in apvs:
        det_id = apv["det"]
        if det_id not in detectors:
            continue
        det = detectors[det_id]
        plane = apv["plane"]
        pos = apv["pos"]
        orient = apv["orient"]
        readout_offset = apv.get("readout_offset", 32)
        strip_offset = apv.get("strip_offset", 0)
        hybrid_board = apv.get("hybrid_board", True)

        # compute plane strip positions using the MapStrip pipeline
        plane_strips = []
        for ch in range(strips_per_apv):
            # step 1: APV25 internal channel mapping
            s = 32 * (ch % 4) + 8 * (ch // 4) - 31 * (ch // 16)
            # step 2: hybrid board pin conversion
            if hybrid_board:
                s = s + 1 + s % 4 - 5 * ((s // 4) % 2)
            # step 3: readout strip mapping (skip if readout_offset == 0)
            if readout_offset > 0:
                if s & 1:
                    s = readout_offset - (s + 1) // 2
                else:
                    s = readout_offset + s // 2
            # step 4: 7-bit mask
            s &= 0x7f

            # step 5: orient flip
            if orient == 1:
                s = 127 - s

            # step 6: plane-wide strip number
            s += strip_offset + pos * strips_per_apv

            plane_strips.append(s)

        if plane == "X":
            pitch = det["x_pitch"]
            strip_positions = sorted(set(plane_strips))
            x_min = min(strip_positions) * pitch
            x_max = (max(strip_positions) + 1) * pitch
            det["x_apv_edges"].add(x_min)
            det["x_apv_edges"].add(x_max)

            for s in plane_strips:
                strip_x = s * pitch
                det["x_strips"].append((strip_x, 0, det["y_size"]))

        elif plane == "Y":
            pitch = det["y_pitch"]
            strip_positions = sorted(set(plane_strips))
            y_min = min(strip_positions) * pitch
            y_max = (max(strip_positions) + 1) * pitch
            det["y_apv_edges"].add(y_min)
            det["y_apv_edges"].add(y_max)

            for s in plane_strips:
                strip_y = s * pitch

                # Y strips are split at hole boundary
                if hole and hole_y0 <= strip_y <= hole_y1:
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
    x_lines = []
    for i, (x, y0, y1) in enumerate(sorted(det["x_strips"])):
        if i % show_every == 0:
            x_lines.append([(x, y0), (x, y1)])
    if x_lines:
        ax.add_collection(LineCollection(x_lines, colors="steelblue",
                                          linewidths=0.3, alpha=0.6))

    # Y strips (horizontal lines) — red
    y_lines = []
    sorted_y = sorted(det["y_strips"])
    for i, (y, x0, x1) in enumerate(sorted_y):
        if i % show_every == 0:
            y_lines.append([(x0, y), (x1, y)])
    if y_lines:
        ax.add_collection(LineCollection(y_lines, colors="indianred",
                                          linewidths=0.3, alpha=0.6))

    # APV boundaries — dashed lines
    for x in sorted(det["x_apv_edges"]):
        ax.axvline(x, color="steelblue", linewidth=0.8, alpha=0.5, linestyle="--")
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
    if len(sys.argv) > 1:
        gem_map_path = sys.argv[1]
    else:
        for candidate in [
            "database/gem_map.json",
            "../database/gem_map.json",
            "gem_map.json",
        ]:
            if os.path.exists(candidate):
                gem_map_path = candidate
                break
        else:
            print("Usage: python gem_layout.py [path/to/gem_map.json]")
            sys.exit(1)

    print(f"Loading: {gem_map_path}")
    layers, apvs, hole = load_gem_map(gem_map_path)
    detectors = build_strip_layout(layers, apvs, hole)

    if hole:
        print(f"Beam hole: {hole['width']}x{hole['height']} mm "
              f"at ({hole['x_center']}, {hole['y_center']})")

    print(f"Detectors: {len(detectors)}")
    for det_id, det in sorted(detectors.items()):
        print(f"  {det['name']}: {det['x_size']:.1f} x {det['y_size']:.1f} mm, "
              f"{len(det['x_strips'])} X strips, {len(det['y_strips'])} Y strip segments")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("PRad-II GEM Strip Layout", fontsize=14)

    for det_id in sorted(detectors.keys()):
        row = det_id // 2
        col = det_id % 2
        plot_detector(axes[row][col], detectors[det_id], det_id, hole)

    plt.tight_layout()
    plt.savefig("gem_layout.png", dpi=150, bbox_inches="tight")
    print("Saved: gem_layout.png")
    plt.show()


if __name__ == "__main__":
    main()
