#!/usr/bin/env python3
"""
merge_hycal_map.py — Join the legacy hycal_modules.json + hycal_daq_map.json
into the unified hycal_map.json schema.

Output schema (one entry per HyCal/LMS/Veto/booster element):

    [
      {"n": "G1",
       "t": "PbGlass",
       "geo": {"sx":..., "sy":..., "x":..., "y":..., "sec":..., "row":..., "col":...},
       "daq": {"crate":..., "slot":..., "channel":...}},
      {"n": "pradbst4",
       "t": "booster",
       "geo": {"sx":..., "sy":..., "x":..., "y":...},
       "bst": {"ip":..., "port":...}},
      ...
    ]

Notes
-----
* `geo`/`daq`/`bst` are each optional.  Detector modules carry geo+daq;
  prad2hvmon boosters carry geo+bst; partially-mapped modules drop daq.
* The DAQ side joins on the module name (the legacy schema used `name`
  rather than `n`).  Mismatches are warned but do not fail the merge.
* The output preserves the input ordering of hycal_modules.json so diffs
  remain readable; each top-level entry serializes onto a single line.

Usage
-----
    # standard PRad-II merge
    python scripts/dev_tool/merge_hycal_map.py \
        --modules database/hycal_modules.json \
        --daq     database/hycal_daq_map.json \
        --out     database/hycal_map.json

    # PRad-1 legacy daq side, modern modules
    python scripts/dev_tool/merge_hycal_map.py \
        --modules database/hycal_modules.json \
        --daq     database/prad1/prad_hycal_daq_map.json \
        --out     database/prad1/prad_hycal_map.json

    # prad2hvmon: no daq side, boosters in `bst` block
    python scripts/dev_tool/merge_hycal_map.py \
        --modules ../prad2hvmon/database/hycal_modules.json \
        --out     ../prad2hvmon/database/hycal_map.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any


# Field groupings — anything outside these stays at the top level.
GEO_KEYS = ("sx", "sy", "x", "y", "sec", "row", "col")
BST_KEYS = ("ip", "port")
DAQ_KEYS = ("crate", "slot", "channel")
TOP_KEYS = ("n", "t")


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def build_daq_index(daq_entries: list[dict]) -> dict[str, dict]:
    """name -> {crate, slot, channel} from a legacy daq map array."""
    out: dict[str, dict] = {}
    for entry in daq_entries:
        name = entry.get("name")
        if not name:
            print(f"WARN: daq entry without 'name': {entry}", file=sys.stderr)
            continue
        if name in out:
            print(f"WARN: duplicate daq entry for {name!r}", file=sys.stderr)
        out[name] = {k: entry[k] for k in DAQ_KEYS if k in entry}
    return out


def merge_module(mod: dict, daq_index: dict[str, dict]) -> "OrderedDict[str, Any]":
    """Reshape one legacy module record into the new schema."""
    name = mod.get("n")
    if not name:
        raise ValueError(f"module record missing 'n': {mod}")

    out: "OrderedDict[str, Any]" = OrderedDict()
    out["n"] = name
    if "t" in mod:
        out["t"] = mod["t"]

    geo = OrderedDict((k, mod[k]) for k in GEO_KEYS if k in mod)
    if geo:
        out["geo"] = geo

    bst = OrderedDict((k, mod[k]) for k in BST_KEYS if k in mod)
    if bst:
        out["bst"] = bst

    daq = daq_index.get(name)
    if daq:
        out["daq"] = OrderedDict((k, daq[k]) for k in DAQ_KEYS if k in daq)

    # Surface any unexpected keys so we don't silently lose them.
    leftover = [k for k in mod
                if k not in TOP_KEYS + GEO_KEYS + BST_KEYS]
    for k in leftover:
        out[k] = mod[k]
        print(f"INFO: preserving unrecognised key {k!r} on {name!r}",
              file=sys.stderr)

    return out


def write_one_line_per_record(records: list[dict], out_path: Path) -> None:
    """Emit `[\\n  rec,\\n  rec,\\n  ...\\n  rec\\n]\\n`.

    Matches the existing hycal_daq_map.json layout (compact rows in a
    pretty-printed array) so git diffs stay row-aligned.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("[\n")
        for i, rec in enumerate(records):
            line = json.dumps(rec, ensure_ascii=False, separators=(", ", ":"))
            tail = "," if i + 1 < len(records) else ""
            f.write(f"  {line}{tail}\n")
        f.write("]\n")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modules", required=True, type=Path,
                    help="Path to hycal_modules.json")
    ap.add_argument("--daq", type=Path, default=None,
                    help="Path to hycal_daq_map.json (omit for modules-only "
                         "merges, e.g. prad2hvmon)")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output hycal_map.json path")
    args = ap.parse_args(argv)

    modules = load_json(args.modules)
    if not isinstance(modules, list):
        print(f"ERROR: {args.modules} is not a JSON array", file=sys.stderr)
        return 2

    daq_index: dict[str, dict] = {}
    if args.daq is not None:
        daq_entries = load_json(args.daq)
        if not isinstance(daq_entries, list):
            print(f"ERROR: {args.daq} is not a JSON array", file=sys.stderr)
            return 2
        daq_index = build_daq_index(daq_entries)

    module_names = {m.get("n") for m in modules if m.get("n")}
    orphan_daq = sorted(n for n in daq_index if n not in module_names)
    if orphan_daq:
        print(f"WARN: {len(orphan_daq)} daq entries have no matching module "
              f"and will be dropped: {orphan_daq}", file=sys.stderr)

    records = [merge_module(m, daq_index) for m in modules]

    if args.daq is not None:
        modules_without_daq = [r["n"] for r in records if "daq" not in r]
        if modules_without_daq:
            print(f"WARN: {len(modules_without_daq)} modules have no daq "
                  f"entry: {modules_without_daq}", file=sys.stderr)

    write_one_line_per_record(records, args.out)
    print(f"Wrote {len(records)} records to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
