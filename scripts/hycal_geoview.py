"""Shared HyCal geo-view widget for PyQt6 scripts.

Provides a common ``Module`` dataclass, colour palettes, and an
extensible ``HyCalMapWidget`` base class.  Scripts in this directory
(hycal_scaler_map, hycal_pedestal_monitor, hycal_map_builder,
hycal_gain_monitor, trigger_mask_editor) subclass the widget to add
overlays, custom fills, or different mouse behaviour.

Typical usage:

    class MyMap(HyCalMapWidget):
        def _paint_modules(self, p):
            # optional custom fill; default uses value colormap
            ...

    w = MyMap(enable_zoom_pan=True)
    w.set_modules(load_modules(MODULES_JSON))
    w.set_values({name: value, ...})
    w.set_range(vmin, vmax)
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, Union

from PyQt6.QtWidgets import (
    QWidget, QPushButton, QSizePolicy, QToolTip,
    QLineEdit, QLabel, QHBoxLayout, QVBoxLayout, QApplication,
)
from PyQt6.QtCore import Qt, QRectF, QPointF, QSize, QTimer, QObject, QEvent, pyqtSignal
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QFontMetricsF, QLinearGradient,
    QPalette, QDoubleValidator,
)


# ===========================================================================
#  Shared theme
# ===========================================================================
#
# ``THEME`` is a class whose *class attributes* hold the currently active
# colours. Stylesheets elsewhere in ``scripts/`` read them through f-strings,
# so calling ``set_theme(name)`` before any widgets are constructed swaps
# every colour in one shot. Add new theme names to ``_THEMES`` and new fields
# to ``THEME`` together; every theme must define every field.
#
# Scripts typically wire this up via a ``--theme`` CLI flag:
#     parser.add_argument("--theme", choices=available_themes(),
#                         default="dark")
#     ...
#     set_theme(args.theme)           # before constructing any window
#     apply_theme_palette(window)


class THEME:
    """Active palette — class attrs are overwritten by :func:`set_theme`.

    Colour vocabulary follows the Apple-inspired design system in
    ``DESIGN-apple.md``: binary dark / light surfaces with a single Apple Blue
    accent reserved for interactive elements.
    """

    # --- surfaces ---
    BG            = "#000000"   # window background (Pure Black)
    CANVAS        = "#000000"   # chart / HyCal map canvas
    PANEL         = "#1d1d1f"   # input surfaces (text edits, tables, combos)
    BUTTON        = "#1d1d1f"   # raised controls (Primary Dark)
    BUTTON_HOVER  = "#28282a"   # button :hover background (Dark Surface 3)
    ALT_BASE      = "#242426"   # alternating table rows (Dark Surface 5)
    TOOLTIP       = "#2a2a2d"   # hover/info tooltip background (Dark Surface 4)

    # --- lines ---
    BORDER        = "#424245"   # subtle border — Apple rarely uses borders
    GRID          = "#1d1d1f"   # chart gridlines (very faint on dark)

    # --- text ---
    TEXT          = "#ffffff"
    TEXT_STRONG   = "#ffffff"
    TEXT_DIM      = "#86868b"   # Apple secondary grey
    TEXT_MUTED    = "#6e6e73"   # tertiary / disabled

    # --- semantic / state ---
    ACCENT        = "#2997ff"   # Bright Blue — links/highlights on dark
    ACCENT_STRONG = "#0071e3"   # Apple Blue — primary CTA
    ACCENT_BORDER = "#0071e3"   # focus ring
    SUCCESS       = "#30d158"   # iOS green (system green dark)
    WARN          = "#ff9f0a"   # iOS orange
    DANGER        = "#ff453a"   # iOS red
    HIGHLIGHT     = "#ff9f0a"   # orange emphasis / drift
    NO_DATA       = "#1d1d1f"   # map-fill when no value

    # --- misc ---
    SELECT_BORDER = "#ffffff"   # selected-module border (white on dark)


_THEMES: Dict[str, Dict[str, str]] = {
    "dark": {
        "BG":            "#000000",
        "BG_SUBTLE":     "#161b22",     # inset plot/panel tile
        "CANVAS":        "#000000",
        "PANEL":         "#1d1d1f",
        "BUTTON":        "#1d1d1f",
        "BUTTON_HOVER":  "#28282a",
        "ALT_BASE":      "#242426",
        "TOOLTIP":       "#2a2a2d",
        "BORDER":        "#424245",
        "GRID":          "#1d1d1f",
        "TEXT":          "#ffffff",
        "TEXT_STRONG":   "#ffffff",
        "TEXT_DIM":      "#86868b",
        "TEXT_MUTED":    "#6e6e73",
        "ACCENT":        "#2997ff",
        "ACCENT_STRONG": "#0071e3",
        "ACCENT_BORDER": "#0071e3",
        "SUCCESS":       "#30d158",
        "WARN":          "#ff9f0a",
        "DANGER":        "#ff453a",
        "HIGHLIGHT":     "#ff9f0a",
        "NO_DATA":       "#1d1d1f",
        "SELECT_BORDER": "#ffffff",
    },
    "light": {
        # Apple light: #ffffff / #f5f5f7 section alternation, #1d1d1f text,
        # #0066cc inline links, #0071e3 CTA blue.
        "BG":            "#ffffff",
        "BG_SUBTLE":     "#f5f5f7",     # inset plot/panel tile
        "CANVAS":        "#f5f5f7",
        "PANEL":         "#ffffff",
        "BUTTON":        "#fafafc",
        "BUTTON_HOVER":  "#ededf2",
        "ALT_BASE":      "#f5f5f7",
        "TOOLTIP":       "#ffffff",
        "BORDER":        "#d2d2d7",
        "GRID":          "#e5e5ea",
        "TEXT":          "#1d1d1f",
        "TEXT_STRONG":   "#000000",
        "TEXT_DIM":      "#6e6e73",
        "TEXT_MUTED":    "#86868b",
        "ACCENT":        "#0066cc",
        "ACCENT_STRONG": "#0071e3",
        "ACCENT_BORDER": "#0071e3",
        "SUCCESS":       "#248a3d",
        "WARN":          "#c93400",
        "DANGER":        "#d70015",
        "HIGHLIGHT":     "#bf5700",
        "NO_DATA":       "#e5e5ea",
        "SELECT_BORDER": "#000000",
    },
}


def available_themes() -> List[str]:
    """Theme names acceptable to :func:`set_theme`."""
    return list(_THEMES.keys())


def set_theme(name: str) -> None:
    """Activate one of :data:`_THEMES` by mutating :class:`THEME` in place.

    Call **before** any window is constructed; stylesheets in this project
    are plain f-strings that read ``THEME.*`` once, at widget creation time.
    Switching themes after the UI is built will not re-render existing
    stylesheets.
    """
    try:
        values = _THEMES[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown theme {name!r}. Available: {available_themes()}"
        ) from exc
    for key, value in values.items():
        setattr(THEME, key, value)


def apply_theme_palette(widget) -> None:
    """Install the active theme's :class:`QPalette` on ``widget``.

    Sets QPalette roles used by top-level windows throughout the
    ``scripts/`` GUIs. Idempotent.
    """
    pal = widget.palette()
    for role, colour in (
        (QPalette.ColorRole.Window,     THEME.BG),
        (QPalette.ColorRole.WindowText, THEME.TEXT),
        (QPalette.ColorRole.Base,       THEME.PANEL),
        (QPalette.ColorRole.Text,       THEME.TEXT),
        (QPalette.ColorRole.Button,     THEME.BUTTON),
        (QPalette.ColorRole.ButtonText, THEME.TEXT),
        (QPalette.ColorRole.Highlight,  THEME.ACCENT),
    ):
        pal.setColor(role, QColor(colour))
    widget.setPalette(pal)


# Back-compat alias for the previous public name.
apply_dark_palette = apply_theme_palette


# ---------------------------------------------------------------------------
#  themed(qss) — rewrite hard-coded dark hex codes to the active theme
# ---------------------------------------------------------------------------
#
# Scripts can keep their existing Qt stylesheet strings unchanged and simply
# wrap ``setStyleSheet(...)`` in ``themed(...)``. The helper rewrites the
# historical dark-theme hex codes to whatever the active :class:`THEME`
# resolves to. New stylesheets may also use this form so adding a new theme
# only requires adding an entry to ``_THEMES``.
#
#   label.setStyleSheet(themed("QLabel{background:#161b22;color:#c9d1d9;}"))
#
# Codes outside ``_QSS_MAP`` are left untouched — useful for one-off accents.

_QSS_MAP: Dict[str, str] = {
    # base surfaces
    "#0a0e14": "CANVAS",
    "#0d1117": "BG",
    "#161b22": "PANEL",
    "#21262d": "BUTTON",
    "#30363d": "BORDER",
    "#131820": "ALT_BASE",
    # text
    "#c9d1d9": "TEXT",
    "#e6edf3": "TEXT_STRONG",
    "#8b949e": "TEXT_DIM",
    "#555555": "TEXT_MUTED",
    "#555":    "TEXT_MUTED",
    # accents / state
    "#58a6ff": "ACCENT",
    "#1f6feb": "ACCENT_STRONG",
    "#388bfd": "ACCENT_BORDER",
    "#3fb950": "SUCCESS",
    "#d29922": "WARN",
    "#f85149": "DANGER",
    "#f97316": "HIGHLIGHT",
    "#ff2222": "DANGER",
    # widget fills
    "#1a1a2e": "NO_DATA",
}


def themed(qss: str) -> str:
    """Rewrite well-known dark-theme hex codes to the active :class:`THEME`.

    Intended for wrapping Qt stylesheets so a single pass of theme swapping
    (``set_theme('light')``) re-colours every widget without per-call
    f-strings. Unknown colour literals pass through unchanged.
    """
    out = qss
    for hex_code, key in _QSS_MAP.items():
        # plain lookup — THEME values are already strings.
        out = out.replace(hex_code, getattr(THEME, key))
    return out


# ===========================================================================
#  Module dataclass
# ===========================================================================

class Module:
    """A HyCal detector module with geometric size and position."""
    __slots__ = ("name", "mod_type", "x", "y", "sx", "sy")

    def __init__(self, name: str, mod_type: str,
                 x: float, y: float, sx: float, sy: float):
        self.name = name
        self.mod_type = mod_type
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy


def load_modules(path: Path) -> List[Module]:
    """Load modules from a JSON file (format: list of {n, t, x, y, sx, sy})."""
    with open(path) as f:
        data = json.load(f)
    return [Module(e["n"], e["t"], e["x"], e["y"], e["sx"], e["sy"])
            for e in data]


# ===========================================================================
#  Colour palettes
# ===========================================================================

# First two palettes (``rainbow`` and ``blue-yellow``) match the
# corresponding palettes in prad2hvmon's web monitor (resources/
# monitor_geo_view.js: ``rainbow``, ``darkblue``) so the desktop and web
# views use the same color language.  Order matters — palette cycle
# starts at index 0.
PALETTES: Dict[str, List[Tuple[float, Tuple[int, int, int]]]] = {
    "rainbow": [
        (0.00, (30,   58,  95)), (0.25, (59,  130, 246)),
        (0.50, (45,  212, 160)), (0.75, (234, 179,   8)),
        (1.00, (245, 101, 101)),
    ],
    "blue-yellow": [
        (0.00, (11,   22,  40)),
        (0.50, (59,  158, 255)),
        (1.00, (234, 179,   8)),
    ],
    "viridis": [
        (0.00, (68,   1,  84)), (0.25, (59,  82, 139)),
        (0.50, (33, 145, 140)), (0.75, (94, 201,  98)),
        (1.00, (253, 231,  37)),
    ],
    "inferno": [
        (0.00, (0,     0,   4)), (0.25, (120,  28, 109)),
        (0.50, (229,  89,  52)), (0.75, (253, 198,  39)),
        (1.00, (252, 255, 164)),
    ],
    "coolwarm": [
        (0.00, (59,   76, 192)), (0.25, (141, 176, 254)),
        (0.50, (221, 221, 221)), (0.75, (245, 148, 114)),
        (1.00, (180,   4,  38)),
    ],
    "hot": [
        (0.00, (11,   0,   0)), (0.33, (230,   0,   0)),
        (0.66, (255, 210,   0)), (1.00, (255, 255, 255)),
    ],
    "blue-orange": [
        (0.00, (10,   42, 110)), (0.25, (30,   90, 180)),
        (0.50, (80,   80,  80)), (0.75, (220, 120,  30)),
        (1.00, (249, 115,  22)),
    ],
    "greyscale": [
        (0.00, (20,   20,  20)), (1.00, (240, 240, 240)),
    ],
}
PALETTE_NAMES: List[str] = list(PALETTES.keys())


def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)


def cmap_qcolor(t: float, stops) -> QColor:
    """Map ``t`` in [0, 1] to a QColor along the given palette stops."""
    t = max(0.0, min(1.0, t))
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1:
            s = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return QColor(_lerp(c0[0], c1[0], s),
                          _lerp(c0[1], c1[1], s),
                          _lerp(c0[2], c1[2], s))
    _, c = stops[-1]
    return QColor(*c)


# ===========================================================================
#  Per-tab view state
# ===========================================================================
#
# ``MapViewState`` is a lightweight snapshot of everything that can vary
# "per view" in a HyCalMapWidget — colormap range, palette, log/linear
# mapping, and (optionally) the data dict itself.  Scripts that show
# multiple tabs/modes against a single (heavy) widget keep one
# ``MapViewState`` per tab and call ``widget.apply_view(state)`` on tab
# switch instead of reconstructing the widget.
#
# See ``HyCalMapWidget.apply_view`` / ``capture_view`` and the
# ``ColorRangeControl`` widget further down.


@dataclass
class MapViewState:
    """Per-tab/per-mode profile of a HyCalMapWidget's view settings.

    ``values=None`` means "leave the widget's current data dict alone" —
    useful when several tabs share the same data and only the *display*
    differs (range, palette, log scale).
    """
    values:      Optional[Dict[str, float]] = None
    vmin:        float = 0.0
    vmax:        float = 1.0
    palette_idx: int   = 0
    log_scale:   bool  = False
    label:       str   = ""


# ===========================================================================
#  HyCal map base widget
# ===========================================================================

class HyCalMapWidget(QWidget):
    """Extensible HyCal geometry view with value → colour mapping.

    Features
    --------
    * Automatic layout: modules laid out in physical coordinates, axis-correct
      (y flipped so positive y is up).
    * Optional colour bar at the bottom (click to cycle palette).
    * Optional zoom/pan (mouse wheel + drag, middle click to reset, overlay
      Reset button top-right).
    * Optional log-scale value mapping.
    * Hover tooltip and module click signal.

    Subclass hooks (override to customise)
    --------------------------------------
    * ``_paint_modules(p)``          — per-module fill loop (default uses
                                        ``set_values`` + current palette).
    * ``_paint_before_modules(p, w, h)`` — drawn after background, before modules.
    * ``_paint_overlays(p, w, h)``   — drawn after modules, before colour bar.
                                        Default paints the hover highlight.
    * ``_paint_after_colorbar(p, w, h)`` — drawn last (legends etc.).
    * ``_colorbar_center_text()``    — palette name line; default shows palette
                                        name and "[log]" flag.
    * ``_fmt_value(v)``              — vmin/vmax label format.
    * ``_tooltip_text(name)``        — tooltip when hovering a module.
    """

    moduleHovered = pyqtSignal(str)
    moduleClicked = pyqtSignal(str)   # "" means deselect
    paletteClicked = pyqtSignal()
    rangeEdited   = pyqtSignal(float, float)   # user-edited via inline editor

    _CLICK_THRESHOLD = 4

    # Colour roles are resolved from the module-level :class:`THEME` at paint
    # time. The property accessors below return the currently active colours
    # so subclasses that used to override ``BG_COLOR`` etc. continue to work
    # via plain attribute assignment.

    @property
    def BG_COLOR(self) -> QColor:
        return QColor(THEME.CANVAS)

    @property
    def NO_DATA_COLOR(self) -> QColor:
        return QColor(THEME.NO_DATA)

    @property
    def HOVER_COLOR(self) -> QColor:
        return QColor(THEME.ACCENT)

    @property
    def CB_BORDER(self) -> QColor:
        return QColor(THEME.ACCENT)

    @property
    def CB_TEXT(self) -> QColor:
        return QColor(THEME.TEXT_DIM)

    @property
    def EMPTY_TEXT(self) -> QColor:
        return QColor(THEME.TEXT_MUTED)

    def __init__(self, parent=None, *,
                 shrink: float = 0.92,
                 margin: int = 12,
                 margin_top: int = 10,
                 margin_bottom: int = 50,
                 include_lms: bool = False,
                 show_colorbar: bool = True,
                 enable_zoom_pan: bool = False,
                 enable_inline_range_edit: bool = True,
                 min_size: Tuple[int, int] = (400, 400)):
        super().__init__(parent)
        self._shrink = shrink
        self._margin = margin
        self._margin_top = margin_top
        self._margin_bottom = margin_bottom
        self._include_lms = include_lms
        self._show_colorbar = show_colorbar
        self._enable_zoom_pan = enable_zoom_pan
        self._enable_inline_range_edit = enable_inline_range_edit

        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(*min_size)

        self._modules: List[Module] = []
        self._values: Dict[str, float] = {}
        self._vmin = 0.0
        self._vmax = 1.0
        self._log_scale = False
        self._palette_idx = 0
        self._hovered: Optional[str] = None
        self._rects: Dict[str, QRectF] = {}
        self._rect_names_rev: List[str] = []
        self._geo_bounds: Tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0)
        self._cb_rect: Optional[QRectF] = None
        self._layout_dirty = True

        # Inline range editor (set up lazily on first colorbar paint).
        # ``_cb_min_hit`` / ``_cb_max_hit`` are screen-space rects covering the
        # vmin / vmax label + pencil glyph; ``_inline_editor`` is a child
        # QLineEdit shown over whichever was clicked.
        self._cb_min_hit: Optional[QRectF] = None
        self._cb_max_hit: Optional[QRectF] = None
        self._inline_editor: Optional[QLineEdit] = None
        self._inline_which: Optional[str] = None
        self._inline_cancelled = False
        self._inline_min_editable = True
        self._inline_max_editable = True

        # zoom / pan state (only used when enable_zoom_pan is True)
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._drag_last: Optional[QPointF] = None
        self._drag_origin: Optional[QPointF] = None
        self._dragging = False

        if enable_zoom_pan:
            self._reset_btn = QPushButton("Reset", self)
            self._reset_btn.setFixedSize(52, 24)
            f = QFont("Consolas", 9)
            f.setBold(True)
            self._reset_btn.setFont(f)
            self._reset_btn.setStyleSheet(
                f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT_DIM};"
                f"border:1px solid {THEME.BORDER};border-radius:8px;}}"
                f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};color:{THEME.TEXT};}}")
            self._reset_btn.clicked.connect(self.reset_view)
        else:
            self._reset_btn = None

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def set_modules(self, modules: List[Module]):
        if self._include_lms:
            self._modules = list(modules)
        else:
            self._modules = [m for m in modules if m.mod_type != "LMS"]
        if self._modules:
            self._geo_bounds = (
                min(m.x - m.sx / 2 for m in self._modules),
                max(m.x + m.sx / 2 for m in self._modules),
                min(m.y - m.sy / 2 for m in self._modules),
                max(m.y + m.sy / 2 for m in self._modules),
            )
        self._layout_dirty = True
        self.update()

    def set_values(self, values: Dict[str, float]):
        self._values = values
        self.update()

    def set_range(self, vmin: float, vmax: float):
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def set_palette(self, idx_or_name):
        """Set palette by index or name."""
        if isinstance(idx_or_name, str):
            idx = PALETTE_NAMES.index(idx_or_name)
        else:
            idx = int(idx_or_name)
        self._palette_idx = idx % len(PALETTES)
        self.update()

    def cycle_palette(self):
        self._palette_idx = (self._palette_idx + 1) % len(PALETTES)
        self.update()

    def set_log_scale(self, on: bool):
        self._log_scale = on
        self.update()

    def is_log_scale(self) -> bool:
        return self._log_scale

    def apply_view(self, view: MapViewState) -> None:
        """Push a per-tab view state into the widget and repaint.

        ``view.values`` is skipped when ``None`` so several tabs can share
        the underlying data dict and differ only in display settings.
        """
        if view.values is not None:
            self._values = view.values
        self._vmin = view.vmin
        self._vmax = view.vmax
        if PALETTES:
            self._palette_idx = view.palette_idx % len(PALETTES)
        self._log_scale = view.log_scale
        self.update()

    def capture_view(self, include_values: bool = False) -> MapViewState:
        """Snapshot the widget's current view settings.

        ``include_values=True`` deep-copies the data dict; default skips
        it (callers usually keep the data dict separate per tab).
        """
        return MapViewState(
            values=dict(self._values) if include_values else None,
            vmin=self._vmin,
            vmax=self._vmax,
            palette_idx=self._palette_idx,
            log_scale=self._log_scale,
        )

    def auto_range(self) -> Tuple[float, float]:
        """Set vmin/vmax from current values (min..max, or min..min+1 if flat)."""
        vals = list(self._values.values())
        if vals:
            self._vmin = min(vals)
            self._vmax = max(vals)
            if self._vmin == self._vmax:
                self._vmax = self._vmin + 1.0
            self.update()
        return self._vmin, self._vmax

    def palette_idx(self) -> int:
        return self._palette_idx

    def palette_stops(self):
        return list(PALETTES.values())[self._palette_idx]

    def reset_view(self):
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._layout_dirty = True
        self.update()

    def value_to_t(self, v: float) -> float:
        """Map a raw value to [0, 1] using current scale (linear or log)."""
        vmin, vmax = self._vmin, self._vmax
        if self._log_scale:
            floor = max(vmin, 1e-9)
            ceil = max(vmax, floor * 10)
            v = max(v, floor)
            return (math.log10(v) - math.log10(floor)) / \
                   (math.log10(ceil) - math.log10(floor))
        return (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5

    # ------------------------------------------------------------------
    #  Layout
    # ------------------------------------------------------------------

    def _recompute_layout(self):
        self._rects.clear()
        if not self._modules:
            self._rect_names_rev = []
            self._layout_dirty = False
            return

        w, h = self.width(), self.height()
        margin, top, bot = self._margin, self._margin_top, self._margin_bottom
        pw, ph = w - 2 * margin, h - top - bot
        x0, x1, y0, y1 = self._geo_bounds
        base_scale = min(pw / max(x1 - x0, 1e-9), ph / max(y1 - y0, 1e-9))
        sc = base_scale * self._zoom
        dw, dh = (x1 - x0) * sc, (y1 - y0) * sc
        ox = margin + (pw - dw) / 2 + self._pan_x
        oy = top + (ph - dh) / 2 + self._pan_y

        # Record layout geometry (useful for subclass overlays)
        self._geo_x0 = x0
        self._geo_y1 = y1
        self._geo_scale = sc
        self._geo_ox = ox
        self._geo_oy = oy

        shrink = self._shrink
        for m in self._modules:
            mw, mh = m.sx * sc * shrink, m.sy * sc * shrink
            cx = ox + (m.x - x0) * sc
            cy = oy + (y1 - m.y) * sc
            self._rects[m.name] = QRectF(cx - mw / 2, cy - mh / 2, mw, mh)
        self._rect_names_rev = list(self._rects)[::-1]
        self._layout_dirty = False

    def geo_to_canvas(self, gx: float, gy: float) -> QPointF:
        """Convert geometry-space coords to widget canvas coords."""
        return QPointF(self._geo_ox + (gx - self._geo_x0) * self._geo_scale,
                       self._geo_oy + (self._geo_y1 - gy) * self._geo_scale)

    def resizeEvent(self, event):
        self._layout_dirty = True
        if self._reset_btn is not None:
            self._reset_btn.move(self.width() - self._reset_btn.width() - 6, 6)
        super().resizeEvent(event)

    # ------------------------------------------------------------------
    #  Painting
    # ------------------------------------------------------------------

    def paintEvent(self, event):
        if self._layout_dirty:
            self._recompute_layout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, self.BG_COLOR)

        if not self._rects:
            self._paint_empty(p, w, h)
            p.end()
            return

        self._paint_before_modules(p, w, h)
        self._paint_modules(p)
        self._paint_overlays(p, w, h)
        if self._show_colorbar:
            self._paint_colorbar(p, w, h)
        self._paint_after_colorbar(p, w, h)
        p.end()

    # -- hook: empty state (no modules loaded) --
    def _paint_empty(self, p: QPainter, w: int, h: int):
        pass

    # -- hook: before modules (title etc.) --
    def _paint_before_modules(self, p: QPainter, w: int, h: int):
        pass

    # -- hook: per-module fill (default: colormap by value) --
    def _paint_modules(self, p: QPainter):
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        vmin, vmax = self._vmin, self._vmax
        log_scale = self._log_scale
        if log_scale:
            log_lo = math.log10(max(vmin, 1e-9))
            log_hi = math.log10(max(vmax, vmin * 10, 1e-8))
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                p.fillRect(rect, no_data)
            else:
                if log_scale:
                    lv = math.log10(max(v, 1e-9))
                    t = (lv - log_lo) / (log_hi - log_lo) if log_hi > log_lo else 0.5
                else:
                    t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                p.fillRect(rect, cmap_qcolor(t, stops))

    # -- hook: after modules, before colorbar --
    def _paint_overlays(self, p: QPainter, w: int, h: int):
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(self.HOVER_COLOR, 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

    # -- hook: after colorbar (legend, extra labels) --
    def _paint_after_colorbar(self, p: QPainter, w: int, h: int):
        pass

    # -- hook: value format in colorbar min/max labels --
    def _fmt_value(self, v: float) -> str:
        if v == 0:
            return "0"
        return f"{v:.6g}"

    # -- hook: colorbar center text --
    def _colorbar_center_text(self) -> str:
        name = PALETTE_NAMES[self._palette_idx]
        if self._log_scale:
            name += "  [log]"
        return name

    # -- hook: maximum colour bar width (default 400) --
    CB_MAX_WIDTH = 400

    def _paint_colorbar(self, p: QPainter, w: int, h: int):
        stops = self.palette_stops()
        cb_w = min(self.CB_MAX_WIDTH, w - 80)
        cb_h = 20
        cb_x = (w - cb_w) / 2
        cb_y = h - 50
        self._cb_rect = QRectF(cb_x, cb_y, cb_w, cb_h)

        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in stops:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(self._cb_rect, QBrush(grad))
        p.setPen(QPen(self.CB_BORDER, 1.0))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(self._cb_rect)

        p.setPen(self.CB_TEXT)
        font = QFont("Consolas", 11)
        p.setFont(font)
        fm = QFontMetricsF(font)
        tick_h   = 5          # tick mark length (px)
        label_y  = cb_y + cb_h + tick_h + 1
        label_h  = fm.height() + 2
        min_str  = self._fmt_value(self._vmin)
        max_str  = self._fmt_value(self._vmax)

        # ── intermediate ticks ────────────────────────────────────────
        # Candidate positions evenly spaced in linear/log, rounded to int.
        # Duplicates, values equal to endpoints, and ticks whose labels
        # would overlap with a neighbour or with the endpoint labels are
        # all discarded.
        import math as _math
        N_CANDS = 9     # generate more candidates than we may show
        vmin, vmax = self._vmin, self._vmax
        if self._log_scale and vmin > 0 and vmax > vmin:
            raw_vals = [10 ** (
                _math.log10(vmin) + (i + 1) / (N_CANDS + 1) *
                (_math.log10(vmax) - _math.log10(vmin))
            ) for i in range(N_CANDS)]
        else:
            span = vmax - vmin if vmax > vmin else 1.0
            raw_vals = [vmin + (i + 1) / (N_CANDS + 1) * span
                        for i in range(N_CANDS)]

        # Round to integer and deduplicate while preserving order
        seen = set()
        int_vals = []
        for rv in raw_vals:
            iv = int(round(rv))
            if iv in seen or iv <= int(round(vmin)) or iv >= int(round(vmax)):
                continue
            seen.add(iv)
            int_vals.append(iv)

        # Convert values → pixel positions, then drop ticks whose label
        # would overlap a neighbour (min pixel gap = label_min_gap).
        def val_to_x(v):
            if self._log_scale and vmin > 0 and vmax > vmin:
                t = (_math.log10(max(v, 1e-9)) - _math.log10(max(vmin, 1e-9))) / \
                    (_math.log10(max(vmax, 1e-9)) - _math.log10(max(vmin, 1e-9)))
            else:
                t = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            return cb_x + t * cb_w

        label_min_gap = 8   # minimum pixel gap between label edges
        accepted = []   # list of (iv, tx, lbl_w)
        for iv in int_vals:
            lbl  = str(iv)
            lw   = fm.horizontalAdvance(lbl) + 4
            tx   = val_to_x(iv)
            # check against already-accepted ticks and endpoint labels
            ok = True
            # endpoint label widths (approx)
            end_w = 130
            if tx - lw / 2 < cb_x + end_w + label_min_gap:
                ok = False
            elif tx + lw / 2 > cb_x + cb_w - end_w - label_min_gap:
                ok = False
            else:
                for _, ptx, plw in accepted:
                    if abs(tx - ptx) < (lw + plw) / 2 + label_min_gap:
                        ok = False
                        break
            if ok:
                accepted.append((iv, tx, lw))

        p.setPen(QPen(self.CB_TEXT, 1.0))
        for iv, tx, lbl_w in accepted:
            p.drawLine(QPointF(tx, cb_y + cb_h),
                       QPointF(tx, cb_y + cb_h + tick_h))
            lbl = str(iv)
            p.drawText(QRectF(tx - lbl_w / 2, label_y, lbl_w, label_h),
                       Qt.AlignmentFlag.AlignCenter, lbl)

        # ── min / max / center text ───────────────────────────────────
        p.setPen(self.CB_TEXT)
        p.drawText(QRectF(cb_x, label_y, 130, label_h),
                   Qt.AlignmentFlag.AlignLeft, min_str)
        p.drawText(QRectF(cb_x + cb_w - 130, label_y, 130, label_h),
                   Qt.AlignmentFlag.AlignRight, max_str)
        # center text (palette name) sits above the bar
        center_font = QFont("Consolas", 10)
        p.setFont(center_font)
        p.drawText(QRectF(cb_x, cb_y - fm.height() - 2, cb_w, fm.height() + 2),
                   Qt.AlignmentFlag.AlignCenter, self._colorbar_center_text())
        p.setFont(font)  # restore

        if self._enable_inline_range_edit:
            # Faint pencil glyph next to each editable label hints
            # "click to edit".  Hit rects cover both the value text and
            # the pencil so the user can click anywhere on either.
            pencil = "✎"   # LOWER RIGHT PENCIL
            pencil_w = fm.horizontalAdvance(pencil) + 4
            min_w = fm.horizontalAdvance(min_str) + 4
            max_w = fm.horizontalAdvance(max_str) + 4
            p.setPen(QColor(THEME.TEXT_DIM))
            if self._inline_min_editable:
                self._cb_min_hit = QRectF(
                    cb_x - 2, label_y - 1, min_w + pencil_w + 4, label_h + 2)
                p.drawText(QRectF(cb_x + min_w, label_y, pencil_w, label_h),
                           Qt.AlignmentFlag.AlignLeft, pencil)
            else:
                self._cb_min_hit = None
            if self._inline_max_editable:
                self._cb_max_hit = QRectF(
                    cb_x + cb_w - max_w - pencil_w - 2, label_y - 1,
                    max_w + pencil_w + 4, label_h + 2)
                p.drawText(QRectF(cb_x + cb_w - max_w - pencil_w, label_y,
                                  pencil_w, label_h),
                           Qt.AlignmentFlag.AlignLeft, pencil)
            else:
                self._cb_max_hit = None
        else:
            self._cb_min_hit = None
            self._cb_max_hit = None

    # ------------------------------------------------------------------
    #  Mouse / hit-test
    # ------------------------------------------------------------------

    def _hit(self, pos) -> Optional[str]:
        for name in self._rect_names_rev:
            if self._rects[name].contains(pos):
                return name
        return None

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        if v is None:
            return name
        return f"{name}: {self._fmt_value(v)}"

    def mousePressEvent(self, e):
        if self._enable_zoom_pan and e.button() == Qt.MouseButton.MiddleButton:
            self.reset_view()
            return
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            self._drag_last = e.position()
            self._drag_origin = e.position()
            self._dragging = False

    def mouseReleaseEvent(self, e):
        if e.button() not in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            return
        if self._dragging:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        elif e.button() == Qt.MouseButton.LeftButton:
            self._handle_click(e.position())
        self._drag_last = None
        self._drag_origin = None
        self._dragging = False

    def _handle_click(self, pos):
        """Default: vmin/vmax label hit → inline edit; colour-bar body →
        paletteClicked; else → moduleClicked."""
        if self._check_inline_range_edit_click(pos):
            return
        if self._cb_rect and self._cb_rect.contains(pos):
            self.paletteClicked.emit()
            return
        name = self._hit(pos)
        self.moduleClicked.emit(name or "")

    def _check_inline_range_edit_click(self, pos) -> bool:
        """Hit-test ``pos`` against the colorbar vmin/vmax labels and open
        the inline editor if it lands on one.  Returns True if handled.

        Subclasses that override ``mousePressEvent`` (e.g. paint editors
        that intercept clicks on press rather than release) should call
        this at the top of their override to keep the inline-edit feature
        working — otherwise their override swallows clicks on the labels.
        """
        if not self._enable_inline_range_edit:
            return False
        if self._cb_min_hit is not None and self._cb_min_hit.contains(pos):
            self._show_inline_editor("min", self._cb_min_hit)
            return True
        if self._cb_max_hit is not None and self._cb_max_hit.contains(pos):
            self._show_inline_editor("max", self._cb_max_hit)
            return True
        return False

    def mouseMoveEvent(self, e):
        # zoom/pan drag
        if self._enable_zoom_pan and self._drag_last is not None:
            pos = e.position()
            if not self._dragging:
                dx = pos.x() - self._drag_origin.x()
                dy = pos.y() - self._drag_origin.y()
                if dx * dx + dy * dy > self._CLICK_THRESHOLD ** 2:
                    self._dragging = True
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
            if self._dragging:
                self._pan_x += pos.x() - self._drag_last.x()
                self._pan_y += pos.y() - self._drag_last.y()
                self._drag_last = pos
                self._layout_dirty = True
                self.update()
            return

        # hover
        pos = e.position()
        if self._enable_inline_range_edit and (
                (self._cb_min_hit is not None and self._cb_min_hit.contains(pos))
                or (self._cb_max_hit is not None and self._cb_max_hit.contains(pos))):
            self.setCursor(Qt.CursorShape.IBeamCursor)
        elif self._cb_rect and self._cb_rect.contains(pos):
            self.setCursor(Qt.CursorShape.PointingHandCursor)
        else:
            self.setCursor(Qt.CursorShape.ArrowCursor)

        found = self._hit(pos)
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                QToolTip.showText(e.globalPosition().toPoint(),
                                  self._tooltip_text(found), self)
                self.moduleHovered.emit(found)
            else:
                QToolTip.hideText()

    def wheelEvent(self, e):
        if not self._enable_zoom_pan:
            return
        factor = 1.15 if e.angleDelta().y() > 0 else 1.0 / 1.15
        new_zoom = max(0.5, min(self._zoom * factor, 20.0))
        if new_zoom == self._zoom:
            return
        pos = e.position()
        ratio = new_zoom / self._zoom
        self._pan_x = pos.x() + (self._pan_x - pos.x()) * ratio
        self._pan_y = pos.y() + (self._pan_y - pos.y()) * ratio
        self._zoom = new_zoom
        self._layout_dirty = True
        self.update()

    # ------------------------------------------------------------------
    #  Inline range edit (colorbar vmin / vmax labels)
    # ------------------------------------------------------------------

    def set_inline_edit_targets(self, min_editable: bool = True,
                                max_editable: bool = True):
        """Lock vmin and/or vmax against inline editing.

        When a side is locked, its pencil glyph is hidden and clicks on
        the label fall through to the colour-bar (palette cycle) handler.
        Used by ``ColorRangeController`` when ``min_fixed`` is set so the
        user can't bypass the lock by clicking the label.
        """
        self._inline_min_editable = bool(min_editable)
        self._inline_max_editable = bool(max_editable)
        self.update()

    def _show_inline_editor(self, which: str, rect: QRectF):
        """Open a child QLineEdit over the clicked vmin/vmax label."""
        if self._inline_editor is None:
            self._inline_editor = QLineEdit(self)
            self._inline_editor.setValidator(
                QDoubleValidator(-1e12, 1e12, 6, self._inline_editor))
            self._inline_editor.setStyleSheet(themed(
                f"QLineEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
                f"border:1px solid {THEME.ACCENT_STRONG};border-radius:3px;"
                f"padding:1px 4px;font:9pt Consolas;}}"))
            self._inline_editor.installEventFilter(self)
            self._inline_editor.editingFinished.connect(self._commit_inline_edit)
        self._inline_which = which
        self._inline_cancelled = False
        cur = self._vmin if which == "min" else self._vmax
        self._inline_editor.setText(self._fmt_value(cur))
        ed_w = max(80, int(rect.width()))
        ed_h = max(20, int(rect.height()) + 4)
        if which == "max":
            ed_x = int(rect.right()) - ed_w
        else:
            ed_x = int(rect.left())
        ed_y = int(rect.top()) - 2
        self._inline_editor.setGeometry(ed_x, ed_y, ed_w, ed_h)
        self._inline_editor.show()
        self._inline_editor.raise_()
        self._inline_editor.selectAll()
        self._inline_editor.setFocus()

    def _commit_inline_edit(self):
        # editingFinished fires on Enter and on focus-loss.  Esc handler
        # sets _inline_cancelled and hides the editor; bail in that case.
        if self._inline_editor is None or self._inline_which is None:
            return
        if self._inline_cancelled:
            self._inline_which = None
            return
        try:
            new_val = float(self._inline_editor.text())
        except ValueError:
            self._inline_editor.hide()
            self._inline_which = None
            return
        if self._inline_which == "min":
            new_min, new_max = new_val, self._vmax
        else:
            new_min, new_max = self._vmin, new_val
        if (math.isfinite(new_min) and math.isfinite(new_max)
                and new_max > new_min):
            self.set_range(new_min, new_max)
            self.rangeEdited.emit(new_min, new_max)
        self._inline_editor.hide()
        self._inline_which = None

    def eventFilter(self, obj, event):
        if obj is self._inline_editor and event.type() == QEvent.Type.KeyPress:
            if event.key() == Qt.Key.Key_Escape:
                self._inline_cancelled = True
                self._inline_editor.hide()
                self._inline_which = None
                return True
        return super().eventFilter(obj, event)

    def sizeHint(self):
        return QSize(680, 680)


# ===========================================================================
#  Reusable colormap-range control
# ===========================================================================
#
# ColorRangeControl wraps the "two min/max edit boxes + Auto button +
# optional Log toggle" pattern that several scripts in this directory
# duplicate.  Bind it to a HyCalMapWidget directly (simple case), to a
# MapViewState (per-tab profile), or to ``(state, widget)`` to drive
# both at once.
#
# Auto-button gestures
# --------------------
#   * Single click               → one-shot fit; pin state unchanged.
#   * Double click               → one-shot fit + enter persistent mode
#                                  (button highlights with ACCENT_STRONG).
#   * Click while pinned         → exit persistent mode.
#   * Double-click while pinned  → exit persistent mode (self-cancelling).
#   * Editing a range field      → exits persistent mode automatically.
#
# Auto-fit strategies
# -------------------
# The ``auto_fit`` argument selects how the Auto button computes a range
# from the current data dict.  Built-in presets:
#   * ``"minmax"``           — full ``min..max``.
#   * ``"minmax_nonzero"``   — ``min..max`` of values != 0 (use when
#                              zero is a sentinel — e.g. masked channels).
#   * ``"percentile"``       — ``np.percentile(values, lo, hi)`` with
#                              ``auto_fit_percentile=(lo, hi)`` (default
#                              ``(2, 98)``).
# Or pass a callable ``f(values: Dict[str, float]) -> (vmin, vmax)`` for
# anything custom.
#
# Migration cookbook
# ------------------
#   ctrl = ColorRangeControl(map_widget,
#                            auto_fit="minmax_nonzero",
#                            include_log=True)
#   layout.addWidget(ctrl)
#   # whenever the data on the map changes:
#   map_widget.set_values(new_values)
#   ctrl.notify_values_changed(new_values)   # re-fits if pinned


class _AutoButton(QPushButton):
    """QPushButton that distinguishes single-click from double-click.

    Qt fires ``clicked`` for both presses of a double-click.  We use a
    counter + ``QApplication.doubleClickInterval()`` timer to disambiguate:
    the first ``clicked`` starts the timer; if a second ``clicked`` arrives
    before the timer expires, it's a double-click; otherwise single.
    """

    oneshotRequested  = pyqtSignal()
    pinToggleRequested = pyqtSignal()

    def __init__(self, text: str = "Auto", parent: Optional[QWidget] = None):
        super().__init__(text, parent)
        self._click_count = 0
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(QApplication.doubleClickInterval())
        self._timer.timeout.connect(self._fire_pending)
        self.clicked.connect(self._on_clicked)

    def _on_clicked(self):
        self._click_count += 1
        if self._click_count == 1:
            self._timer.start()
        else:
            self._timer.stop()
            self._click_count = 0
            self.pinToggleRequested.emit()

    def _fire_pending(self):
        if self._click_count == 1:
            self._click_count = 0
            self.oneshotRequested.emit()
        else:
            self._click_count = 0


class ColorRangeController(QObject):
    """Headless controller for a HyCalMapWidget colormap range.

    Holds auto-fit / pin / log / min-fixed behaviour and exposes signals
    and methods so scripts can wire whatever GUI they want — a single
    "Auto" toolbar button, a menu action, the inline edits on the
    colorbar, or no GUI at all.

    Parameters
    ----------
    target
        ``HyCalMapWidget``                  — push range/log into the widget.
        ``MapViewState``                    — push into a per-tab profile.
        ``(MapViewState, HyCalMapWidget)``  — drive both.
    auto_fit
        ``"minmax"`` / ``"minmax_nonzero"`` / ``"percentile"`` or a
        callable ``f(values) -> (vmin, vmax)``.
    auto_fit_percentile
        ``(lo, hi)`` for the ``"percentile"`` preset.  Default ``(2, 98)``.
    min_fixed
        Lock the lower bound at this value (auto-fit and ``set_range``
        clamp vmin).  ``None`` (default) means both bounds are free.

    Signals
    -------
    rangeChanged(vmin, vmax) — any range change (auto-fit, set_range, or
                               relayed widget inline edit).
    autoPinned(on)           — persistent auto mode flipped.
    logToggled(on)           — log/linear flipped.

    Public API
    ----------
    notify_values_changed(values)  — call when the data dict changes
                                     (re-fits if pinned).
    auto_fit(values=None)          — one-shot programmatic fit.
    set_range(vmin, vmax)          — push a range; turns off pin.
    set_pinned(on) / toggle_pinned() / is_pinned()
    set_log(on)
    set_state(view)                — rebind to a different MapViewState.
    Properties: vmin, vmax, pinned, min_fixed, values.
    """

    rangeChanged = pyqtSignal(float, float)
    autoPinned   = pyqtSignal(bool)
    logToggled   = pyqtSignal(bool)

    _AUTO_FIT_PRESETS = ("minmax", "minmax_nonzero", "percentile")
    AutoFit = Union[str, Callable[[Dict[str, float]], Tuple[float, float]]]

    def __init__(self,
                 target,
                 *,
                 auto_fit: AutoFit = "minmax",
                 auto_fit_percentile: Tuple[float, float] = (2.0, 98.0),
                 min_fixed: Optional[float] = None,
                 parent: Optional[QObject] = None):
        super().__init__(parent)

        self._map: Optional[HyCalMapWidget] = None
        self._state: Optional[MapViewState] = None
        if isinstance(target, HyCalMapWidget):
            self._map = target
        elif isinstance(target, MapViewState):
            self._state = target
        elif (isinstance(target, tuple) and len(target) == 2 and
              isinstance(target[0], MapViewState) and
              isinstance(target[1], HyCalMapWidget)):
            self._state, self._map = target
        else:
            raise TypeError(
                "ColorRangeController target must be HyCalMapWidget, "
                "MapViewState, or (MapViewState, HyCalMapWidget)")

        if not (callable(auto_fit) or auto_fit in self._AUTO_FIT_PRESETS):
            raise ValueError(
                f"auto_fit must be callable or one of {self._AUTO_FIT_PRESETS}; "
                f"got {auto_fit!r}")
        self._auto_fit = auto_fit
        self._auto_pct = auto_fit_percentile
        self._min_fixed: Optional[float] = (
            float(min_fixed) if min_fixed is not None else None)

        self._pinned = False
        self._values: Dict[str, float] = {}

        # Relay widget inline edits — they always override pin (manual edit).
        if self._map is not None:
            self._map.rangeEdited.connect(self._on_widget_range_edited)

        # Honour min_fixed up front so any pre-existing target range obeys it.
        if self._min_fixed is not None:
            _, vmax = self._read_target_range()
            if vmax <= self._min_fixed:
                vmax = self._min_fixed + 1.0
            self._push_range(self._min_fixed, vmax)
            # Lock the widget's inline-edit min so it can't be overridden.
            if self._map is not None:
                self._map.set_inline_edit_targets(
                    min_editable=False, max_editable=True)

    # ---- target plumbing -------------------------------------------------

    def _read_target_range(self) -> Tuple[float, float]:
        if self._state is not None:
            return self._state.vmin, self._state.vmax
        if self._map is not None:
            return self._map._vmin, self._map._vmax
        return 0.0, 1.0

    def _push_range(self, vmin: float, vmax: float):
        if self._min_fixed is not None:
            vmin = self._min_fixed
        if self._state is not None:
            self._state.vmin = vmin
            self._state.vmax = vmax
        if self._map is not None:
            self._map.set_range(vmin, vmax)

    # ---- public API ------------------------------------------------------

    def notify_values_changed(self, values: Dict[str, float]):
        """Call when the data dict changes; re-fits if pinned."""
        self._values = values or {}
        if self._pinned:
            self._do_auto_fit_and_apply()

    def auto_fit(self, values: Optional[Dict[str, float]] = None):
        """One-shot fit.  Pin state unchanged.  ``values`` overrides cache."""
        if values is not None:
            self._values = values
        self._do_auto_fit_and_apply()

    def set_range(self, vmin: float, vmax: float):
        """Push a range; clamps vmin if min_fixed is set; turns off pin."""
        if self._min_fixed is not None:
            vmin = self._min_fixed
        if not (math.isfinite(vmin) and math.isfinite(vmax)) or vmax <= vmin:
            return
        if self._pinned:
            self._set_pinned(False)
        self._push_range(vmin, vmax)
        self.rangeChanged.emit(vmin, vmax)

    def set_pinned(self, on: bool):
        self._set_pinned(bool(on))

    def toggle_pinned(self):
        self._set_pinned(not self._pinned)

    def is_pinned(self) -> bool:
        return self._pinned

    def set_log(self, on: bool):
        on = bool(on)
        if self._state is not None:
            self._state.log_scale = on
        if self._map is not None:
            self._map.set_log_scale(on)
        self.logToggled.emit(on)

    def set_state(self, view: MapViewState):
        """Rebind to a different MapViewState (per-tab profile switch)."""
        self._state = view
        self.rangeChanged.emit(*self._read_target_range())

    @property
    def pinned(self) -> bool:
        return self._pinned

    @property
    def min_fixed(self) -> Optional[float]:
        return self._min_fixed

    @property
    def vmin(self) -> float:
        return self._read_target_range()[0]

    @property
    def vmax(self) -> float:
        return self._read_target_range()[1]

    @property
    def values(self) -> Dict[str, float]:
        return self._values

    # ---- internal --------------------------------------------------------

    def _on_widget_range_edited(self, vmin: float, vmax: float):
        """Inline edit on the widget's colorbar — always exits pin."""
        if self._pinned:
            self._set_pinned(False)
        if self._state is not None:
            self._state.vmin = vmin
            self._state.vmax = vmax
        self.rangeChanged.emit(vmin, vmax)

    def _set_pinned(self, on: bool):
        if self._pinned == on:
            return
        self._pinned = on
        self.autoPinned.emit(on)

    def _do_auto_fit_and_apply(self):
        vmin, vmax = self._compute_auto_fit()
        if self._min_fixed is not None:
            vmin = self._min_fixed
        if not (math.isfinite(vmin) and math.isfinite(vmax)) or vmax <= vmin:
            pad = max(abs(vmin) * 0.05, 1e-6)
            vmax = vmin + pad
        self._push_range(vmin, vmax)
        self.rangeChanged.emit(vmin, vmax)

    def _compute_auto_fit(self) -> Tuple[float, float]:
        values = self._values
        if callable(self._auto_fit):
            return tuple(self._auto_fit(values))
        if not values:
            return self._read_target_range()
        if self._auto_fit == "minmax":
            vals = [v for v in values.values()
                    if v is not None and not (isinstance(v, float)
                                              and math.isnan(v))]
            if not vals:
                return 0.0, 1.0
            return float(min(vals)), float(max(vals))
        if self._auto_fit == "minmax_nonzero":
            vals = [v for v in values.values()
                    if v is not None and v != 0.0 and
                    not (isinstance(v, float) and math.isnan(v))]
            if not vals:
                return 0.0, 1.0
            return float(min(vals)), float(max(vals))
        if self._auto_fit == "percentile":
            try:
                import numpy as np
            except ImportError:
                vals = list(values.values())
                return float(min(vals)), float(max(vals))
            arr = np.asarray(list(values.values()), dtype=float)
            arr = arr[np.isfinite(arr)]
            if arr.size == 0:
                return 0.0, 1.0
            lo, hi = self._auto_pct
            return (float(np.percentile(arr, lo)),
                    float(np.percentile(arr, hi)))
        return self._read_target_range()


class ColorRangeControl(QWidget):
    """Default min/max + Auto button + Log toggle widget.

    Thin shim over :class:`ColorRangeController` that provides the
    classic two-edit row UI.  Construct it for the simple case; for
    custom UIs (single button, inline-only, …) build a
    ``ColorRangeController`` directly and wire your own widgets.

    See :class:`ColorRangeController` for the auto-fit / pin /
    min_fixed / target-binding semantics; this widget passes those
    parameters straight through.

    Extra parameters
    ----------------
    include_log    — add an inline Log toggle button.
    orientation    — ``"horizontal"`` or ``"vertical"``.
    start_pinned   — start in persistent auto-fit mode.

    The underlying controller is exposed as ``self.controller`` for
    advanced wiring.
    """

    rangeChanged = pyqtSignal(float, float)
    autoPinned   = pyqtSignal(bool)
    logToggled   = pyqtSignal(bool)

    AutoFit = ColorRangeController.AutoFit

    def __init__(self,
                 target,
                 *,
                 auto_fit: AutoFit = "minmax",
                 auto_fit_percentile: Tuple[float, float] = (2.0, 98.0),
                 include_log: bool = False,
                 orientation: str = "horizontal",
                 start_pinned: bool = False,
                 min_fixed: Optional[float] = None,
                 parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._ctrl = ColorRangeController(
            target,
            auto_fit=auto_fit,
            auto_fit_percentile=auto_fit_percentile,
            min_fixed=min_fixed,
            parent=self,
        )
        self._build_ui(orientation, include_log)
        # Re-emit + sync on controller events.
        self._ctrl.rangeChanged.connect(self._on_ctrl_range_changed)
        self._ctrl.autoPinned.connect(self._on_ctrl_pin_changed)
        self._ctrl.logToggled.connect(self._on_ctrl_log_changed)
        # Initial sync.
        self._set_edits(self._ctrl.vmin, self._ctrl.vmax)
        if start_pinned:
            self._ctrl.set_pinned(True)

    # ---- accessors -------------------------------------------------------

    @property
    def controller(self) -> ColorRangeController:
        return self._ctrl

    # Pass-through public API for back-compat with the old widget.
    def notify_values_changed(self, values: Dict[str, float]):
        self._ctrl.notify_values_changed(values)

    def auto_fit(self, values: Optional[Dict[str, float]] = None):
        self._ctrl.auto_fit(values)

    def set_range(self, vmin: float, vmax: float):
        self._ctrl.set_range(vmin, vmax)

    def set_state(self, view: MapViewState):
        self._ctrl.set_state(view)
        # Pull updated state into the widgets.
        self._set_edits(self._ctrl.vmin, self._ctrl.vmax)
        if self._log_btn is not None:
            self._log_btn.blockSignals(True)
            self._log_btn.setChecked(view.log_scale)
            self._log_btn.blockSignals(False)
            self._update_log_btn_style()

    def is_pinned(self) -> bool:
        return self._ctrl.is_pinned()

    # ---- UI construction --------------------------------------------------

    def _build_ui(self, orientation: str, include_log: bool):
        edit_css = themed(
            f"QLineEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};border-radius:4px;"
            f"padding:2px 6px;}}")
        self._min_edit = QLineEdit()
        self._min_edit.setMaximumWidth(90)
        self._min_edit.setValidator(
            QDoubleValidator(-1e12, 1e12, 6, self._min_edit))
        self._min_edit.editingFinished.connect(self._on_edit)
        self._min_edit.setStyleSheet(edit_css)

        self._max_edit = QLineEdit()
        self._max_edit.setMaximumWidth(90)
        self._max_edit.setValidator(
            QDoubleValidator(-1e12, 1e12, 6, self._max_edit))
        self._max_edit.editingFinished.connect(self._on_edit)
        self._max_edit.setStyleSheet(edit_css)

        self._auto_btn = _AutoButton("Auto", self)
        self._auto_btn.setToolTip(
            "Click: auto-fit once   ·   Double-click: keep auto-fitting")
        self._auto_btn.oneshotRequested.connect(self._on_auto_oneshot)
        self._auto_btn.pinToggleRequested.connect(self._on_auto_double)
        self._update_auto_btn_style()

        self._log_btn: Optional[QPushButton] = None
        if include_log:
            self._log_btn = QPushButton("Log")
            self._log_btn.setCheckable(True)
            self._log_btn.toggled.connect(self._on_log_clicked)
            self._update_log_btn_style()

        self.setStyleSheet(themed(
            f"QLabel{{color:{THEME.TEXT};background:transparent;}}"))

        single_value = self._ctrl.min_fixed is not None
        if orientation == "vertical":
            outer = QVBoxLayout(self)
            outer.setContentsMargins(0, 0, 0, 0)
            outer.setSpacing(4)
            if not single_value:
                r1 = QHBoxLayout(); r1.setSpacing(4)
                r1.addWidget(QLabel("min:"))
                r1.addWidget(self._min_edit); r1.addStretch()
                outer.addLayout(r1)
            r2 = QHBoxLayout(); r2.setSpacing(4)
            r2.addWidget(QLabel("max:"))
            r2.addWidget(self._max_edit); r2.addStretch()
            r3 = QHBoxLayout(); r3.setSpacing(6)
            r3.addWidget(self._auto_btn)
            if self._log_btn is not None:
                r3.addWidget(self._log_btn)
            r3.addStretch()
            outer.addLayout(r2); outer.addLayout(r3)
        else:
            row = QHBoxLayout(self)
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(6)
            if single_value:
                row.addWidget(QLabel("Max:"))
                row.addWidget(self._max_edit)
            else:
                row.addWidget(QLabel("Range:"))
                row.addWidget(self._min_edit)
                row.addWidget(QLabel("–"))
                row.addWidget(self._max_edit)
            row.addWidget(self._auto_btn)
            if self._log_btn is not None:
                row.addWidget(self._log_btn)
            row.addStretch()

    def _update_auto_btn_style(self):
        if self._ctrl.is_pinned():
            self._auto_btn.setStyleSheet(themed(
                f"QPushButton{{background:{THEME.ACCENT_STRONG};"
                f"color:{THEME.TEXT};border:1px solid {THEME.ACCENT_STRONG};"
                f"padding:5px 14px;font:10pt;border-radius:6px;}}"))
        else:
            self._auto_btn.setStyleSheet(themed(
                f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT};"
                f"border:1px solid {THEME.BORDER};padding:5px 14px;"
                f"font:10pt;border-radius:6px;}}"
                f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}"))

    def _update_log_btn_style(self):
        if self._log_btn is None:
            return
        if self._log_btn.isChecked():
            self._log_btn.setStyleSheet(themed(
                f"QPushButton{{background:{THEME.ACCENT};color:{THEME.TEXT};"
                f"border:1px solid {THEME.ACCENT};padding:5px 14px;"
                f"font:10pt;border-radius:6px;}}"))
        else:
            self._log_btn.setStyleSheet(themed(
                f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT_DIM};"
                f"border:1px solid {THEME.BORDER};padding:5px 14px;"
                f"font:10pt;border-radius:6px;}}"
                f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};"
                f"color:{THEME.TEXT};}}"))

    def _set_edits(self, vmin: float, vmax: float):
        self._min_edit.blockSignals(True)
        self._max_edit.blockSignals(True)
        self._min_edit.setText(f"{vmin:.6g}")
        self._max_edit.setText(f"{vmax:.6g}")
        self._min_edit.blockSignals(False)
        self._max_edit.blockSignals(False)

    # ---- handlers --------------------------------------------------------

    def _on_edit(self):
        try:
            if self._ctrl.min_fixed is not None:
                vmin = self._ctrl.min_fixed
            else:
                vmin = float(self._min_edit.text())
            vmax = float(self._max_edit.text())
        except ValueError:
            return
        if vmax <= vmin:
            return
        self._ctrl.set_range(vmin, vmax)

    def _on_auto_oneshot(self):
        if self._ctrl.is_pinned():
            self._ctrl.set_pinned(False)
        else:
            self._ctrl.auto_fit()

    def _on_auto_double(self):
        if self._ctrl.is_pinned():
            self._ctrl.set_pinned(False)
        else:
            self._ctrl.auto_fit()
            self._ctrl.set_pinned(True)

    def _on_log_clicked(self, on: bool):
        self._update_log_btn_style()
        self._ctrl.set_log(on)

    # ---- controller relays -----------------------------------------------

    def _on_ctrl_range_changed(self, vmin: float, vmax: float):
        self._set_edits(vmin, vmax)
        self.rangeChanged.emit(vmin, vmax)

    def _on_ctrl_pin_changed(self, on: bool):
        self._update_auto_btn_style()
        self.autoPinned.emit(on)

    def _on_ctrl_log_changed(self, on: bool):
        if self._log_btn is not None and self._log_btn.isChecked() != on:
            self._log_btn.blockSignals(True)
            self._log_btn.setChecked(on)
            self._log_btn.blockSignals(False)
            self._update_log_btn_style()
        self.logToggled.emit(on)
