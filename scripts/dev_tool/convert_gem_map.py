#!/usr/bin/env python3
"""
convert_gem_map.py — Convert the upstream APV mapping text format
(docs/gem/gem_map_prad2*.txt) into the gem_map.json structure we
consume in database/gem_map.json.

Usage
-----
    python scripts/dev_tool/convert_gem_map.py INPUT.txt OUTPUT.json \
        [--crate-map OLD=NEW,OLD=NEW,...] \
        [--template database/gem_map.json]

Example (HallB layout: crate 146 -> logical crate 1, crate 147 -> 2)
    python scripts/dev_tool/convert_gem_map.py \
        docs/gem/gem_map_prad2_hallB.txt \
        database/gem_map.json \
        --crate-map 146=1,147=2 \
        --template database/gem_map.json

Text format (per upstream mpd_gem_view_ssp)
-------------------------------------------
    Layer, LayerID, ChambersPerLayer, readout, XOffset, YOffset, GEMType,
           x_apvs, y_apvs, x_pitch, y_pitch, x_flip, y_flip
    APV,   CrateID, Layer, FiberID, GEMID, dim(0=X,1=Y), ADCId, I2C,
           Pos, Invert, other(normal/split), backplane, GEMPOS

Mapping to JSON
---------------
    crate    <- crate_map.get(CrateID, CrateID)
    mpd      <- FiberID
    adc      <- ADCId
    det      <- GEMID - 1            (chamber id, 0..3)
    plane    <- "X" / "Y"             (from dim)
    orient   <- Invert                (0 = non-inverted side, 1 = inverted)
    pos      <- Pos                   (APV position in plane, 0..11 for X)
    det_pos  <- GEMPOS

Split-APV fields are added automatically by pattern:
    plane=X, orient=0, pos=11  ->  pin_rotate=16, shared_pos=10, match="-Y"
    plane=X, orient=1, pos=10  ->  match="+Y"

These are the two APVs in each chamber that straddle the beam hole.

The script preserves header config (apv_channels, hole, thresholds, ...)
from --template if given; otherwise it writes built-in defaults.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Defaults — used only when --template is not provided.
# ---------------------------------------------------------------------------
DEFAULT_HEADER: "OrderedDict[str, Any]" = OrderedDict([
    ("apv_channels", 128),
    ("readout_center", 32),
    ("common_mode_threshold", 20.0),
    ("zero_suppression_threshold", 5.0),
    ("cross_talk_threshold", 8.0),
    ("position_resolution", 0.08),
    ("reject_first_timebin", True),
    ("reject_last_timebin", True),
    ("min_peak_adc", 30.0),
    ("min_sum_adc", 60.0),
    ("match_mode", 1),
    ("match_adc_asymmetry", 0.8),
    ("match_time_diff", 50.0),
    ("match_ts_period", 25.0),
    ("hole", OrderedDict([("width", 52.0), ("height", 52.0)])),
])


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

# Match `# X-dimension MPD slot 2` or `# Y-dimension MPD17` etc.
MPD_COMMENT_RE = re.compile(r"^\s*#\s*([XY])-dimension\s+(.*?)\s*$", re.IGNORECASE)


def parse_text(path: Path) -> tuple[list[dict], list[dict], list[tuple]]:
    """Parse a gem_map_prad2*.txt file.

    Returns (layers, apvs, mpd_comments) where:
        layers        — list of physical layer descriptors
        apvs          — list of APV records (in document order)
        mpd_comments  — list of ((line_no, key) -> comment) tuples; key is
                        the FiberID of the *next* APV after the comment.
                        The emitter walks APVs in document order and
                        prints whichever pending comment matches.
    """
    layers: list[dict] = []
    apvs: list[dict] = []
    mpd_comments: list[tuple] = []  # (next_apv_index, plane, slot_label)
    pending_comment: tuple[str, str] | None = None

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            m = MPD_COMMENT_RE.match(raw)
            if m:
                pending_comment = (m.group(1).upper(), m.group(2).strip())
                continue

            if line.startswith("#"):
                continue

            parts = [p.strip() for p in line.rstrip(",").split(",")]
            kind = parts[0]

            if kind == "Layer" and len(parts) >= 13:
                lid, chambers, readout, xoff, yoff, gemtype, \
                    xapvs, yapvs, xpitch, ypitch, xflip, yflip = parts[1:13]
                layers.append({
                    "layer_id": int(lid),
                    "chambers_per_layer": int(chambers),
                    "readout": readout,
                    "x_offset": float(xoff),
                    "y_offset": float(yoff),
                    "gem_type": gemtype,
                    "x_apvs": int(xapvs),
                    "y_apvs": int(yapvs),
                    "x_pitch": float(xpitch),
                    "y_pitch": float(ypitch),
                    "x_flip": int(xflip),
                    "y_flip": int(yflip),
                })
            elif kind == "APV" and len(parts) >= 13:
                vals = parts[1:14]
                rec = {
                    "_crate_orig": int(vals[0]),
                    "_layer_phys": int(vals[1]),
                    "mpd":   int(vals[2]),  # FiberID
                    "_gemid": int(vals[3]),
                    "_dim":   int(vals[4]),
                    "adc":   int(vals[5]),
                    "_i2c":  int(vals[6]),
                    "pos":   int(vals[7]),
                    "orient": int(vals[8]),
                    "_other": vals[9],
                    "_backplane": int(vals[10]),
                    "det_pos": int(vals[11]) if len(vals) > 11 else 0,
                }
                if pending_comment is not None:
                    mpd_comments.append((len(apvs),
                                         pending_comment[0],
                                         pending_comment[1]))
                    pending_comment = None
                apvs.append(rec)

    return layers, apvs, mpd_comments


# ---------------------------------------------------------------------------
# Transform
# ---------------------------------------------------------------------------

def remap_and_finalize(apvs: list[dict],
                       crate_map: dict[int, int]) -> list[dict]:
    """Apply crate remap, drop helper fields, and add split-APV fields."""
    out: list[dict] = []
    for r in apvs:
        crate = crate_map.get(r["_crate_orig"], r["_crate_orig"])
        plane = "X" if r["_dim"] == 0 else "Y"
        rec = OrderedDict([
            ("crate", crate),
            ("mpd",   r["mpd"]),
            ("adc",   r["adc"]),
            ("det",   r["_gemid"] - 1),
            ("plane", plane),
            ("orient", r["orient"]),
            ("pos",   r["pos"]),
            ("det_pos", r["det_pos"]),
        ])
        # Split APV adjacent to beam hole — non-inverted side, X plane, top.
        if plane == "X" and rec["orient"] == 0 and rec["pos"] == 11:
            rec["pin_rotate"] = 16
            rec["shared_pos"] = 10
            rec["match"] = "-Y"
        # Inverted side, X plane, top (paired with the above for matching).
        elif plane == "X" and rec["orient"] == 1 and rec["pos"] == 10:
            rec["match"] = "+Y"
        out.append(rec)
    return out


def build_layers_section(text_layers: list[dict],
                         apvs: list[dict]) -> list[OrderedDict]:
    """Build the per-chamber `layers` list that the C++ loader consumes.

    The text file describes physical layers (each holding 1 or more chambers);
    the JSON `layers` array historically holds one entry per chamber. We
    derive chamber count from the GEMIDs actually present in the APV table,
    and copy x_apvs / y_apvs / pitches from the matching physical layer.
    """
    # Group GEMIDs by physical layer (from APV records).
    layer_to_gems: dict[int, set[int]] = {}
    for r in apvs:
        layer_to_gems.setdefault(r.get("_layer_phys", 0), set()).add(r["_gemid"])

    out: list[OrderedDict] = []
    layer_lookup = {L["layer_id"]: L for L in text_layers}
    for layer_id in sorted(layer_to_gems.keys()):
        L = layer_lookup.get(layer_id, {
            "x_apvs": 12, "y_apvs": 24, "x_pitch": 0.4, "y_pitch": 0.4,
            "gem_type": "PRADGEM",
        })
        for gemid in sorted(layer_to_gems[layer_id]):
            det_id = gemid - 1
            out.append(OrderedDict([
                ("id", det_id),
                ("name", f"GEM{det_id}"),
                ("layer_id", layer_id),
                ("type", L.get("gem_type", "PRADGEM")),
                ("x_apvs", L.get("x_apvs", 12)),
                ("y_apvs", L.get("y_apvs", 24)),
                ("x_pitch", L.get("x_pitch", 0.4)),
                ("y_pitch", L.get("y_pitch", 0.4)),
            ]))
    return out


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def _compact_dict(d: dict) -> str:
    """Render a dict on a single line, preserving key order."""
    parts = []
    for k, v in d.items():
        parts.append(f'"{k}": {json.dumps(v, ensure_ascii=False)}')
    return "{" + ", ".join(parts) + "}"


def render_apvs_array(apvs: list[dict],
                      raw_apvs: list[dict],
                      mpd_comments: list[tuple],
                      indent: str = "        ") -> str:
    """Pretty-print the apvs array with chamber + MPD comments interleaved."""
    # Build mapping from APV index -> (plane, slot_label) from the raw stream.
    idx_comment: dict[int, tuple[str, str]] = {}
    for idx, plane, label in mpd_comments:
        idx_comment[idx] = (plane, label)

    # Per-chamber crate / fiber summary for the chamber header.
    chamber_crates: dict[int, dict[int, list[int]]] = {}
    for rec in apvs:
        det = rec["det"]
        chamber_crates.setdefault(det, {}).setdefault(rec["crate"], []).append(rec["mpd"])

    def chamber_header(det: int) -> str:
        crates = chamber_crates[det]
        bits = []
        for crate in sorted(crates):
            fibers = sorted(set(crates[crate]))
            if len(fibers) == 1:
                rng = f"fiber {fibers[0]}"
            else:
                rng = f"fibers {fibers[0]}-{fibers[-1]}"
            bits.append(f"crate {crate} ({rng})")
        return (f'{indent}{{"// ====== Chamber {det} (GEM{det}) — '
                f'{", ".join(bits)} ======": ""}},')

    lines: list[str] = []
    last_det: int | None = None
    last_mpd_key: tuple[int, int] | None = None  # (crate, mpd)

    for i, rec in enumerate(apvs):
        det = rec["det"]
        if det != last_det:
            if lines:
                lines.append("")
            lines.append(chamber_header(det))
            lines.append("")
            last_det = det
            last_mpd_key = None

        mpd_key = (rec["crate"], rec["mpd"])
        if mpd_key != last_mpd_key:
            if last_mpd_key is not None:
                lines.append("")
            comment = idx_comment.get(i)
            if comment is not None:
                plane, label = comment
                text = f"{plane}-dimension {label} (fiber {rec['mpd']})"
            else:
                text = f"{rec['plane']}-plane fiber {rec['mpd']}"
            lines.append(f'{indent}{{"// {text}": ""}},')
            last_mpd_key = mpd_key

        lines.append(f'{indent}{_compact_dict(rec)},')

    # Drop trailing comma on the last entry.
    if lines:
        last = lines[-1]
        if last.endswith(","):
            lines[-1] = last[:-1]

    return "\n".join(lines)


def render_layers_array(layers: list[dict], indent: str = "        ") -> str:
    rows = [f"{indent}{_compact_dict(L)}," for L in layers]
    if rows:
        rows[-1] = rows[-1][:-1]
    return "\n".join(rows)


def emit_json(header: dict,
              layers: list[dict],
              apvs: list[dict],
              raw_apvs: list[dict],
              mpd_comments: list[tuple],
              source: Path,
              crate_map: dict[int, int]) -> str:
    """Assemble the final JSON document as a string."""
    lines: list[str] = []
    lines.append("{")
    lines.append(f'    "// GEM APV mapping for PRad-II": "",')
    lines.append(f'    "// Generated by scripts/dev_tool/convert_gem_map.py from {source.as_posix()}": "",')
    if crate_map:
        lines.append(f'    "// Crate remap applied: '
                     f'{", ".join(f"{k}->{v}" for k, v in sorted(crate_map.items()))}": "",')
    lines.append(f'    "// Edit the source text + re-run the script; do not hand-edit '
                 f'unless you also update the source.": "",')
    lines.append("")

    # Header config
    for key, val in header.items():
        if key == "hole" and isinstance(val, dict):
            lines.append(f'    "hole": {{')
            inner = []
            for k, v in val.items():
                inner.append(f'        "{k}": {json.dumps(v)}')
            lines.append(",\n".join(inner))
            lines.append("    },")
        elif isinstance(val, str) and key.startswith("//"):
            lines.append(f'    {json.dumps(key)}: {json.dumps(val)},')
        else:
            lines.append(f'    {json.dumps(key)}: {json.dumps(val)},')

    lines.append("")
    lines.append('    "layers": [')
    lines.append(render_layers_array(layers, indent="        "))
    lines.append("    ],")
    lines.append("")
    lines.append('    "apvs": [')
    lines.append(render_apvs_array(apvs, raw_apvs, mpd_comments, indent="        "))
    lines.append("    ]")
    lines.append("}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_template_header(path: Path) -> "OrderedDict[str, Any]":
    """Pull header config (everything except `layers` and `apvs`) from an
    existing JSON file. Stripped of its `// comment` keys so we can re-emit
    fresh provenance comments."""
    with path.open("r", encoding="utf-8") as f:
        j = json.load(f, object_pairs_hook=OrderedDict)
    header: "OrderedDict[str, Any]" = OrderedDict()
    for k, v in j.items():
        if k in ("layers", "apvs"):
            continue
        if isinstance(k, str) and k.startswith("//"):
            continue
        header[k] = v
    return header


def parse_crate_map(spec: str) -> dict[int, int]:
    if not spec:
        return {}
    out: dict[int, int] = {}
    for pair in spec.split(","):
        pair = pair.strip()
        if not pair:
            continue
        old, new = pair.split("=")
        out[int(old.strip())] = int(new.strip())
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert gem_map_prad2*.txt -> gem_map.json")
    ap.add_argument("input_txt", type=Path,
                    help="Input text file (e.g. docs/gem/gem_map_prad2_hallB.txt)")
    ap.add_argument("output_json", type=Path,
                    help="Output JSON file (e.g. database/gem_map.json)")
    ap.add_argument("--crate-map", default="",
                    help="Comma-separated OLD=NEW pairs (e.g. 146=1,147=2)")
    ap.add_argument("--template", type=Path, default=None,
                    help="Existing JSON to source header config from")
    args = ap.parse_args()

    crate_map = parse_crate_map(args.crate_map)
    text_layers, raw_apvs, mpd_comments = parse_text(args.input_txt)
    apvs = remap_and_finalize(raw_apvs, crate_map)
    layers = build_layers_section(text_layers, raw_apvs)

    header = (load_template_header(args.template)
              if args.template else OrderedDict(DEFAULT_HEADER))

    out = emit_json(header, layers, apvs, raw_apvs, mpd_comments,
                    args.input_txt, crate_map)
    args.output_json.write_text(out, encoding="utf-8")

    # Summary
    crate_counts: dict[int, int] = {}
    fiber_counts: dict[int, set] = {}
    for rec in apvs:
        crate_counts[rec["crate"]] = crate_counts.get(rec["crate"], 0) + 1
        fiber_counts.setdefault(rec["crate"], set()).add(rec["mpd"])
    print(f"Wrote {args.output_json}")
    print(f"  detectors: {len(layers)}  apvs: {len(apvs)}")
    for c in sorted(crate_counts):
        print(f"  crate {c}: {crate_counts[c]} APVs across "
              f"{len(fiber_counts[c])} fibers (mpd ids: "
              f"{min(fiber_counts[c])}..{max(fiber_counts[c])})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
