"""
Shared types, constants, and helpers for HyCal calibration tools.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple


# ============================================================================
#  CONSTANTS
# ============================================================================

# Transporter coordinates when the beam hits HyCal centre (0, 0)
BEAM_CENTER_X: float = -126.75   # mm
BEAM_CENTER_Y: float = 10.11     # mm

# Transporter travel limits (symmetric about centre)
_LIMIT_RB_X = -582.65
_LIMIT_RB_Y = -672.50
PTRANS_X_MIN = _LIMIT_RB_X
PTRANS_X_MAX = 2 * BEAM_CENTER_X - _LIMIT_RB_X
PTRANS_Y_MIN = _LIMIT_RB_Y
PTRANS_Y_MAX = 2 * BEAM_CENTER_Y - _LIMIT_RB_Y

# Default database path (relative to this file → ../database/)
DEFAULT_DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "database", "hycal_map.json")


# ============================================================================
#  COLOUR PALETTE (dark control-room theme)
# ============================================================================

class C:
    BG       = "#0d1117"
    PANEL    = "#161b22"
    BORDER   = "#30363d"
    TEXT     = "#c9d1d9"
    DIM      = "#8b949e"
    ACCENT   = "#58a6ff"
    GREEN    = "#3fb950"
    YELLOW   = "#d29922"
    RED      = "#f85149"
    ORANGE   = "#db6d28"
    # canvas module states
    MOD_TODO      = "#21262d"
    MOD_CURRENT   = "#d29922"
    MOD_DWELL     = "#3fb950"
    MOD_DONE      = "#1f6feb"
    MOD_ERROR     = "#f85149"
    MOD_SELECTED  = "#db6d28"
    MOD_INPATH    = "#3fb950"
    # display-only / inactive
    MOD_GLASS     = "#162230"
    MOD_PWO4_BG   = "#1a2a1a"
    MOD_LMS       = "#2d1f3d"
    MOD_EXCLUDED  = "#111418"
    MOD_SKIPPED   = "#15181d"
    PATH_LINE     = "#ffffff"


# ============================================================================
#  MODULE DATA
# ============================================================================

@dataclass
class Module:
    name: str
    mod_type: str      # "PbWO4", "PbGlass", "LMS"
    x: float           # centre x in HyCal frame (mm)
    y: float           # centre y in HyCal frame (mm)
    sx: float          # module width  (mm)
    sy: float          # module height (mm)
    sector: str = ""   # "Center", "Top", "Right", "Bottom", "Left", "LMS"
    row: int = 0       # row within sector (1-indexed)
    col: int = 0       # col within sector (1-indexed)


def load_modules(json_path: str = DEFAULT_DB_PATH) -> List[Module]:
    """Load all modules from hycal_map.json (entries with a "geo" block)."""
    with open(json_path) as f:
        data = json.load(f)
    modules: List[Module] = []
    for entry in data:
        g = entry.get("geo")
        if not g:
            continue
        modules.append(Module(
            name=entry["n"], mod_type=entry["t"],
            x=g["x"], y=g["y"], sx=g["sx"], sy=g["sy"],
            sector=g.get("sec", ""),
            row=g.get("row", 0), col=g.get("col", 0),
        ))
    return modules


# ============================================================================
#  COORDINATE TRANSFORMS
# ============================================================================

def module_to_ptrans(mx: float, my: float) -> Tuple[float, float]:
    """HyCal-frame module centre --> transporter set-point."""
    return (BEAM_CENTER_X + mx, BEAM_CENTER_Y - my)


def ptrans_to_module(px: float, py: float) -> Tuple[float, float]:
    """Transporter position --> beam position on HyCal (HyCal-frame)."""
    return (px - BEAM_CENTER_X, BEAM_CENTER_Y - py)


def ptrans_in_limits(x: float, y: float) -> bool:
    """Check if a transporter position is within travel limits."""
    return (PTRANS_X_MIN <= x <= PTRANS_X_MAX and
            PTRANS_Y_MIN <= y <= PTRANS_Y_MAX)


# ============================================================================
#  LG LAYER FILTER
# ============================================================================

def filter_scan_modules(all_modules: List[Module], lg_layers: int,
                        lg_sx: float = 38.15, lg_sy: float = 38.15
                        ) -> List[Module]:
    """All PbWO4 + PbGlass within lg_layers of PbWO4 bounding box."""
    pwo4 = [m for m in all_modules if m.mod_type == "PbWO4"]
    if not pwo4:
        return list(all_modules)
    min_x = min(m.x for m in pwo4)
    max_x = max(m.x for m in pwo4)
    min_y = min(m.y for m in pwo4)
    max_y = max(m.y for m in pwo4)
    scan = list(pwo4)
    if lg_layers > 0:
        mx = lg_layers * lg_sx
        my = lg_layers * lg_sy
        for m in all_modules:
            if m.mod_type == "PbGlass" and \
               min_x - mx <= m.x <= max_x + mx and \
               min_y - my <= m.y <= max_y + my:
                scan.append(m)
    return scan


# ============================================================================
#  DARK QSS THEME (PyQt6 calibration tools)
# ============================================================================

DARK_QSS = """
QMainWindow, QWidget { background: #0d1117; color: #c9d1d9; }
QLabel { color: #c9d1d9; font: 13pt 'Consolas'; }
QGroupBox {
    color: #58a6ff; font: bold 13pt 'Consolas';
    border: 1px solid #30363d; border-radius: 4px;
    margin-top: 8px; padding-top: 14px;
}
QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; }
QPushButton {
    background: #161b22; color: #c9d1d9; border: 1px solid #30363d;
    border-radius: 3px; padding: 5px 12px; font: 13pt 'Consolas';
}
QPushButton:hover { background: #21262d; border-color: #58a6ff; }
QPushButton:pressed { background: #0d1117; }
QPushButton:disabled { color: #6e7681; border-color: #21262d; }
QPushButton[cssClass="green"]  { background: #238636; color: white; border-color: #2ea043; }
QPushButton[cssClass="green"]:hover  { background: #2ea043; }
QPushButton[cssClass="green"]:disabled  { background: #1a3a22; color: #5a7060; border-color: #1f4528; }
QPushButton[cssClass="danger"] { background: #da3633; color: white; border-color: #f85149; }
QPushButton[cssClass="danger"]:hover { background: #f85149; }
QPushButton[cssClass="danger"]:disabled { background: #3a1818; color: #785050; border-color: #4a1f1f; }
QPushButton[cssClass="warn"]   { background: #9e6a03; color: white; border-color: #d29922; }
QPushButton[cssClass="warn"]:hover   { background: #d29922; }
QPushButton[cssClass="warn"]:disabled   { background: #3a2a08; color: #786850; border-color: #4a3510; }
QPushButton[cssClass="accent"] { background: #1f6feb; color: white; border-color: #388bfd; }
QPushButton[cssClass="accent"]:hover { background: #388bfd; }
QPushButton[cssClass="accent"]:disabled { background: #16294a; color: #506578; border-color: #1d365e; }
QComboBox {
    background: #161b22; color: #c9d1d9; border: 1px solid #30363d;
    border-radius: 3px; padding: 3px 8px; font: 13pt 'Consolas'; min-height: 20px;
}
QComboBox::drop-down { border: none; width: 20px; }
QComboBox::down-arrow { image: none; border-left: 4px solid transparent;
    border-right: 4px solid transparent; border-top: 5px solid #8b949e; margin-right: 6px; }
QComboBox QAbstractItemView {
    background: #161b22; color: #c9d1d9; border: 1px solid #30363d;
    selection-background-color: #1f6feb; selection-color: white;
}
QSpinBox, QDoubleSpinBox {
    background: #161b22; color: #c9d1d9; border: 1px solid #30363d;
    border-radius: 3px; padding: 3px 4px; font: 13pt 'Consolas'; min-height: 20px;
}
QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    background: #21262d; border: none; width: 16px;
}
QTextEdit, QPlainTextEdit {
    background: #0d1117; color: #8b949e; border: 1px solid #30363d;
    border-radius: 3px; font: 13pt 'Consolas';
}
QProgressBar {
    background: #161b22; border: 1px solid #30363d; border-radius: 3px;
    text-align: center; color: #c9d1d9; font: 12pt 'Consolas'; max-height: 18px;
}
QProgressBar::chunk { background: #1f6feb; border-radius: 2px; }
QScrollArea { border: none; background: transparent; }
QScrollBar:vertical { background: #0d1117; width: 8px; border: none; }
QScrollBar::handle:vertical { background: #30363d; border-radius: 3px; min-height: 20px; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QToolTip { background: #161b22; color: #c9d1d9; border: 1px solid #30363d;
    font: 13pt 'Consolas'; padding: 4px; }
"""
