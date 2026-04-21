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
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import QWidget, QPushButton, QSizePolicy, QToolTip
from PyQt6.QtCore import Qt, QRectF, QPointF, QSize, pyqtSignal
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QLinearGradient, QPalette,
)


# ===========================================================================
#  Shared dark theme
# ===========================================================================

class THEME:
    """Central dark-theme colour palette shared by every ``scripts/`` GUI."""
    BG       = "#0d1117"   # window background
    PANEL    = "#161b22"   # secondary surfaces (text edits, tables)
    BUTTON   = "#21262d"   # raised controls
    BORDER   = "#30363d"
    TEXT     = "#c9d1d9"
    TEXT_DIM = "#8b949e"
    ACCENT   = "#58a6ff"


def apply_dark_palette(widget) -> None:
    """Install the shared dark colour palette on ``widget``.

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

PALETTES: Dict[str, List[Tuple[float, Tuple[int, int, int]]]] = {
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
    "rainbow": [
        (0.00, (30,   58,  95)), (0.25, (59,  130, 246)),
        (0.50, (45,  212, 160)), (0.75, (234, 179,   8)),
        (1.00, (245, 101, 101)),
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

    _CLICK_THRESHOLD = 4

    # Default (dark) theme.  set_light_theme(True) switches the instance to
    # the light palette; subclasses can override any of these before/after
    # super().__init__() to customise a single entry.
    BG_COLOR       = QColor("#0a0e14")
    NO_DATA_COLOR  = QColor("#1a1a2e")
    HOVER_COLOR    = QColor("#58a6ff")
    CB_BORDER      = QColor("#58a6ff")
    CB_TEXT        = QColor("#8b949e")
    EMPTY_TEXT     = QColor("#555555")

    _THEME_DARK = {
        "BG_COLOR":      QColor("#0a0e14"),
        "NO_DATA_COLOR": QColor("#1a1a2e"),
        "HOVER_COLOR":   QColor("#58a6ff"),
        "CB_BORDER":     QColor("#58a6ff"),
        "CB_TEXT":       QColor("#8b949e"),
        "EMPTY_TEXT":    QColor("#555555"),
    }
    _THEME_LIGHT = {
        "BG_COLOR":      QColor("#f6f8fa"),
        "NO_DATA_COLOR": QColor("#d0d7de"),
        "HOVER_COLOR":   QColor("#0969da"),
        "CB_BORDER":     QColor("#0969da"),
        "CB_TEXT":       QColor("#57606a"),
        "EMPTY_TEXT":    QColor("#8b949e"),
    }

    def __init__(self, parent=None, *,
                 shrink: float = 0.92,
                 margin: int = 12,
                 margin_top: int = 10,
                 margin_bottom: int = 50,
                 include_lms: bool = False,
                 show_colorbar: bool = True,
                 enable_zoom_pan: bool = False,
                 min_size: Tuple[int, int] = (400, 400)):
        super().__init__(parent)
        self._shrink = shrink
        self._margin = margin
        self._margin_top = margin_top
        self._margin_bottom = margin_bottom
        self._include_lms = include_lms
        self._show_colorbar = show_colorbar
        self._enable_zoom_pan = enable_zoom_pan

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
        self._light_theme = False

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
                "QPushButton{background:rgba(22,27,34,200);color:#8b949e;"
                "border:1px solid #30363d;border-radius:3px;}"
                "QPushButton:hover{background:rgba(33,38,45,220);color:#c9d1d9;}")
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

    def set_light_theme(self, on: bool):
        """Switch between dark (default) and light background schemes."""
        self._light_theme = bool(on)
        theme = self._THEME_LIGHT if on else self._THEME_DARK
        for k, v in theme.items():
            setattr(self, k, v)
        self._on_theme_changed()
        self.update()

    def is_light_theme(self) -> bool:
        return self._light_theme

    def _on_theme_changed(self):
        """Subclass hook — called after theme colors are swapped."""
        pass

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
        cb_h = 14
        cb_x = (w - cb_w) / 2
        cb_y = h - 40
        self._cb_rect = QRectF(cb_x, cb_y, cb_w, cb_h)

        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in stops:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(self._cb_rect, QBrush(grad))
        p.setPen(QPen(self.CB_BORDER, 1.0))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(self._cb_rect)

        p.setPen(self.CB_TEXT)
        p.setFont(QFont("Consolas", 9))
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, 120, 14),
                   Qt.AlignmentFlag.AlignLeft, self._fmt_value(self._vmin))
        p.drawText(QRectF(cb_x + cb_w - 120, cb_y + cb_h + 2, 120, 14),
                   Qt.AlignmentFlag.AlignRight, self._fmt_value(self._vmax))
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, cb_w, 14),
                   Qt.AlignmentFlag.AlignCenter, self._colorbar_center_text())

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
        """Default: colour-bar hit → paletteClicked, else → moduleClicked."""
        if self._cb_rect and self._cb_rect.contains(pos):
            self.paletteClicked.emit()
            return
        name = self._hit(pos)
        self.moduleClicked.emit(name or "")

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
        if self._cb_rect and self._cb_rect.contains(pos):
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

    def sizeHint(self):
        return QSize(680, 680)
