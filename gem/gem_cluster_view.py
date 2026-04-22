#!/usr/bin/env python3
"""
Back-compat wrapper: forwards to ``gem_event_viewer.py --json`` so the
rendering code lives in a single place.

Keeps the original CLI shape:
    gem_cluster_view.py <event.json> [<event.json>...]
    gem_cluster_view.py <dir_with_gem_event*.json>
    gem_cluster_view.py '*.json'
    gem_cluster_view.py -G gem_map.json --det N -o out.png <event.json>
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    viewer = script_dir / "gem_event_viewer.py"

    # gem_event_viewer.py --json accepts the same positional file/dir/glob list
    # that this wrapper used to accept directly, and the same -G / --det / -o
    # flags.  We route everything through --json.
    argv = sys.argv[1:]
    # Split flags (start with '-') from positional paths, insert --json in
    # front of the positional list so argparse can parse it.
    flags: list[str] = []
    paths: list[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("-"):
            flags.append(a)
            # Value follows for every flag we forward here.
            if "=" not in a and i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                flags.append(argv[i + 1])
                i += 1
        else:
            paths.append(a)
        i += 1

    cmd = [sys.executable, str(viewer), *flags]
    if paths:
        cmd += ["--json", *paths]
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    sys.exit(main())
