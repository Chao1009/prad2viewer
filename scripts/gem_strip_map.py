"""
Shared GEM strip mapping — thin wrapper over prad2py.det.map_strip /
map_apv_strips.

The 6-step pipeline that maps an APV25 channel index to a plane-wide strip
number lives in C++ (prad2det/src/GemSystem.cpp: gem::MapStrip) and is
shared between on-line reconstruction and these off-line scripts.  This
module exists so existing callers can continue to do

    from gem_strip_map import map_strip, map_apv_strips

without caring about the binding path.  If you're writing new code, feel
free to import from ``prad2py.det`` directly.

Build requirement: the ``prad2py`` pybind11 module must be built and
importable.  Configure with ``-DBUILD_PYTHON=ON`` and either install it
(``cmake --build build --target install``) or prepend ``build/python/`` to
``PYTHONPATH``.  When the module can't be found we raise a clear
ImportError — we used to fall back to a pure-Python copy here, but that
was prone to drifting from the C++ implementation.
"""

from __future__ import annotations

import os as _os
import sys as _sys


def _resolve_prad2py():
    """Import prad2py.det, auto-discovering build/python/ next to the repo."""
    try:
        from prad2py import det
        return det
    except ImportError:
        pass

    # Repo-local fallback: add ../build/python to sys.path and retry once.
    here = _os.path.dirname(_os.path.abspath(__file__))
    candidate = _os.path.abspath(_os.path.join(here, "..", "build", "python"))
    if _os.path.isdir(candidate) and candidate not in _sys.path:
        _sys.path.insert(0, candidate)
        try:
            from prad2py import det
            return det
        except ImportError:
            pass

    raise ImportError(
        "gem_strip_map requires the prad2py pybind11 module (prad2py.det).\n"
        "Build it with:\n"
        "    cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build\n"
        "then either install it or add build/python/ to PYTHONPATH."
    )


_det = _resolve_prad2py()


def map_strip(ch, plane_index, orient, pin_rotate=0, shared_pos=-1,
              hybrid_board=True, apv_channels=128, readout_center=32):
    """Map APV channel to plane-wide strip number.

    Thin wrapper over ``prad2py.det.map_strip``.  Returns ``(local_strip,
    plane_strip)`` to preserve the original API — `local_strip` here is
    just the final plane-wide strip with the plane offset undone (useful
    for some callers that wanted the per-plane local value).  If you only
    need the plane-wide number, ``prad2py.det.map_strip`` returns it
    directly.
    """
    plane = _det.map_strip(ch=ch,
                           plane_index=plane_index,
                           orient=orient,
                           pin_rotate=pin_rotate,
                           shared_pos=shared_pos,
                           hybrid_board=hybrid_board,
                           apv_channels=apv_channels,
                           readout_center=readout_center)

    # Reconstruct the local (pre-offset) strip number so legacy callers that
    # unpack `(local, plane_strip)` keep working.  Mirrors the plane_shift
    # math in gem::MapStrip (steps 4-6).
    eff_pos = shared_pos if shared_pos >= 0 else plane_index
    plane_shift = (eff_pos - plane_index) * apv_channels - pin_rotate
    local = plane - (plane_shift + plane_index * apv_channels)
    return local, plane


def map_apv_strips(apv, apv_channels=128, readout_center=32):
    """Map all channels of an APV entry (from gem_map.json) to plane strip
    numbers.

    Returns a list of length ``apv_channels`` — the plane-wide strip number
    for each APV channel.  Accepts the same dict shape the JSON loader
    produces (keys: ``pos``, ``orient``, optional ``pin_rotate``,
    ``shared_pos``, ``hybrid_board``).
    """
    return _det.map_apv_strips(
        plane_index=apv["pos"],
        orient=apv["orient"],
        pin_rotate=apv.get("pin_rotate", 0),
        shared_pos=apv.get("shared_pos", -1),
        hybrid_board=apv.get("hybrid_board", True),
        apv_channels=apv_channels,
        readout_center=readout_center,
    )
