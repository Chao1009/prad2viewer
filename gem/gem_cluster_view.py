#!/usr/bin/env python3
"""
Visualize GEM clustering results from ``gem_dump -m evdump`` JSON output.

Thin front-end: reads the per-event JSON file(s) produced by gem_dump and
hands the shaped data to :mod:`gem_view` for rendering.  For live EVIO-
file stepping (no JSON roundtrip), see ``gem_event_viewer.py``.

Usage:
    python gem_cluster_view.py <event.json> [gem_map.json] [--det N] [-o file.png]
"""

import argparse
import json
import os
import sys

import matplotlib.pyplot as plt

from gem_layout import build_strip_layout, load_gem_map
from gem_view import build_apv_map, draw_event, process_zs_hits


def load_event(path):
    raw = open(path, "rb").read()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        text = raw.decode("utf-16")
    elif raw[:3] == b"\xef\xbb\xbf":
        text = raw.decode("utf-8-sig")
    else:
        text = raw.decode("utf-8")
    return json.loads(text)


def print_event_summary(det_list, det_hits):
    for dd in det_list:
        did = dd["id"]
        hits = det_hits.get(did, {"x": [], "y": []})
        xcl = dd.get("x_clusters", [])
        ycl = dd.get("y_clusters", [])
        print(f"\n  {dd['name']}: {len(hits['x'])} X hits, {len(hits['y'])} Y hits, "
              f"{len(xcl)}+{len(ycl)} clusters, {len(dd.get('hits_2d',[]))} 2D hits")
        if xcl or ycl:
            print(f"  {'plane':>5} {'pos(mm)':>8} {'peak':>8} {'total':>8} "
                  f"{'size':>4} {'tbin':>4} {'xtalk':>5}  strips")
            print(f"  {'-'*5:>5} {'-'*8:>8} {'-'*8:>8} {'-'*8:>8} "
                  f"{'-'*4:>4} {'-'*4:>4} {'-'*5:>5}  {'-'*10}")
            for plane, cls in [("X", xcl), ("Y", ycl)]:
                for cl in cls:
                    strips = cl.get("hit_strips", [])
                    srange = f"{min(strips)}-{max(strips)}" if strips else ""
                    print(f"  {plane:>5} {cl['position']:>8.2f} {cl['peak_charge']:>8.1f} "
                          f"{cl['total_charge']:>8.1f} {cl['size']:>4} "
                          f"{cl['max_timebin']:>4} {'y' if cl.get('cross_talk') else '':>5}  "
                          f"{srange}")
        if dd.get("hits_2d"):
            print("  2D hits: " +
                  "  ".join(f"({h['x']:.1f}, {h['y']:.1f})" for h in dd["hits_2d"]))


def render_event_to_file(event_path, detectors, apv_map, hole, raw,
                         det_filter=-1, output=None):
    """Render one event JSON to a PNG file. Returns output path."""
    event = load_event(event_path)
    det_list = event.get("detectors", [])
    if det_filter >= 0:
        det_list = [d for d in det_list if d["id"] == det_filter]
    if not det_list:
        return None

    det_hits = process_zs_hits(event.get("zs_apvs", []), apv_map,
                               detectors, hole, raw)

    n = len(det_list)
    ref = detectors[min(detectors.keys())]
    cell_w = 6
    cell_h = cell_w * ref["y_size"] / ref["x_size"]
    fig = plt.figure(figsize=(cell_w * n, cell_h + 1.5), constrained_layout=True)

    ev_num = event.get("event_number", "?")
    draw_event(fig, detectors, det_list, det_hits, hole,
               title=f"GEM Cluster View -- Event #{ev_num}",
               det_filter=det_filter)

    if output is None:
        output = os.path.splitext(event_path)[0] + ".png"
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output


def main():
    import glob as globmod

    parser = argparse.ArgumentParser(
        description="Visualize GEM clustering from gem_dump -m evdump JSON. "
                    "Accepts a single file, directory, or glob pattern.")
    parser.add_argument("event_json", nargs="+",
                        help="Event JSON file(s), directory, or glob pattern. "
                             "Accepts multiple files (shell-expanded globs).")
    parser.add_argument("-G", "--gem-map", default=None,
                        help="GEM map JSON (default: auto-search common paths)")
    parser.add_argument("--det", type=int, default=-1,
                        help="Show only detector N (default: all)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output PNG (single file mode only)")
    args = parser.parse_args()

    gem_map_path = args.gem_map
    if not gem_map_path:
        # Env var (set by setup.sh / setup.csh) wins, then common source-tree
        # and CWD-relative fallbacks for dev use.
        env_db = os.environ.get("PRAD2_DATABASE_DIR")
        cands = []
        if env_db:
            cands.append(os.path.join(env_db, "gem_map.json"))
        here = os.path.dirname(os.path.abspath(__file__))
        cands.append(os.path.join(here, "..", "database", "gem_map.json"))
        cands += ["database/gem_map.json", "../database/gem_map.json",
                  "../../database/gem_map.json", "gem_map.json"]
        for c in cands:
            if os.path.exists(c):
                gem_map_path = c; break
    if not gem_map_path:
        print("Error: cannot find gem_map.json (pass -G <path>)"); sys.exit(1)

    # Accept three input styles, in any mix:
    #   - one or more explicit JSON files (shell-expanded globs land here)
    #   - a directory → pick up gem_event*.json inside
    #   - a quoted glob pattern → expand ourselves
    files = []
    for arg in args.event_json:
        if os.path.isdir(arg):
            files += sorted(globmod.glob(os.path.join(arg, "gem_event*.json")))
        elif "*" in arg or "?" in arg:
            files += sorted(globmod.glob(arg))
        else:
            files.append(arg)
    files = [f for f in files if f.lower().endswith(".json")]
    if not files:
        print("Error: no JSON files found"); sys.exit(1)

    print(f"GEM map    : {gem_map_path}")
    layers, gem_map_apvs, hole, raw = load_gem_map(gem_map_path)
    detectors = build_strip_layout(layers, gem_map_apvs, hole, raw)
    apv_map = build_apv_map(gem_map_apvs)

    print(f"Files      : {len(files)}")

    for i, fpath in enumerate(files):
        fname = os.path.basename(fpath)
        event = load_event(fpath)
        if not isinstance(event, dict) or "detectors" not in event:
            print(f"\n[{i+1}/{len(files)}] {fname} -- skipped (not an event file)")
            continue
        det_list = event.get("detectors", [])
        if args.det >= 0:
            det_list = [d for d in det_list if d["id"] == args.det]
        det_hits = process_zs_hits(event.get("zs_apvs", []), apv_map,
                                   detectors, hole, raw)

        print(f"\n[{i+1}/{len(files)}] {fname}")
        print_event_summary(det_list, det_hits)

        out = args.output if (args.output and len(files) == 1) else None
        result = render_event_to_file(fpath, detectors, apv_map, hole,
                                      raw, args.det, out)
        if result:
            print(f"  -> {result}")

    print(f"\nDone: {len(files)} file(s) processed.")


if __name__ == "__main__":
    main()
