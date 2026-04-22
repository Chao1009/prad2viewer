#!/usr/bin/env python3
"""
Back-compat wrapper: forwards to ``gem_event_viewer.py --layout`` so the
rendering code lives in a single place.

Keeps the original CLI shape:
    gem_layout.py                      # uses default gem_map.json, writes gem_layout.png
    gem_layout.py -G path/to/gem_map.json
    gem_layout.py -o layout.png
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    viewer = script_dir / "gem_event_viewer.py"
    cmd = [sys.executable, str(viewer), "--layout", *sys.argv[1:]]
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    sys.exit(main())
