"""
HyCal geo-view widget with integrated scaler overlay.

Combines the snake-scan ``HyCalScanMapWidget`` (module colours, path preview,
limit box, motor crosshair) with the scaler-map colour-palette/mapping code
so that live scaler rates can be rendered underneath the scan state layer.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

from PyQt6.QtWidgets import QWidget, QSizePolicy
from PyQt6.QtCore import Qt, QRectF, QPointF, pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QPen, QFont

from scan_utils import (
    C, Module, module_to_ptrans,
    BEAM_CENTER_X, BEAM_CENTER_Y,
    PTRANS_X_MIN, PTRANS_X_MAX, PTRANS_Y_MIN, PTRANS_Y_MAX,
)

# ===========================================================================
#  Colour palettes
# ===========================================================================

PALETTES = {
    "viridis": [
        (0.00, (68, 1, 84)), (0.25, (59, 82, 139)),
        (0.50, (33, 145, 140)), (0.75, (94, 201, 98)),
        (1.00, (253, 231, 37)),
    ],
    "inferno": [
        (0.00, (0, 0, 4)), (0.25, (120, 28, 109)),
        (0.50, (229, 89, 52)), (0.75, (253, 198, 39)),
        (1.00, (252, 255, 164)),
    ],
    "coolwarm": [
        (0.00, (59, 76, 192)), (0.25, (141, 176, 254)),
        (0.50, (221, 221, 221)), (0.75, (245, 148, 114)),
        (1.00, (180, 4, 38)),
    ],
    "hot": [
        (0.00, (11, 0, 0)), (0.33, (230, 0, 0)),
        (0.66, (255, 210, 0)), (1.00, (255, 255, 255)),
    ],
    "rainbow": [
        (0.00, (30, 58, 95)), (0.25, (59, 130, 246)),
        (0.50, (45, 212, 160)), (0.75, (234, 179, 8)),
        (1.00, (245, 101, 101)),
    ],
}
PALETTE_NAMES = list(PALETTES.keys())


def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)


def _cmap_qcolor(t: float, stops) -> QColor:
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
#  HyCal geo-view widget
# ===========================================================================

class HyCalScanMapWidget(QWidget):
    moduleClicked = pyqtSignal(str)
    PAD = 8
    SHRINK = 0.90

    def __init__(self, all_modules, parent=None):
        super().__init__(parent)
        self._drawn = [m for m in all_modules if m.mod_type != "LMS"]
        self._mod_by_name = {m.name: m for m in all_modules}
        self._colors = {}           # scan state colours
        self._path_line = []
        self._dash_line = []
        self._marker_hx = self._marker_hy = None
        self._highlight = None
        self._hover_name = None
        self._lim_hx_min = PTRANS_X_MIN - BEAM_CENTER_X
        self._lim_hx_max = PTRANS_X_MAX - BEAM_CENTER_X
        self._lim_hy_min = BEAM_CENTER_Y - PTRANS_Y_MAX
        self._lim_hy_max = BEAM_CENTER_Y - PTRANS_Y_MIN
        self._scale = 1.0
        self._ox = self._oy = 0.0
        self._x_min = self._y_max = 0.0
        self._rects: Dict[str, QRectF] = {}
        # --- scaler overlay state ---
        self._scaler_values: Dict[str, float] = {}
        self._scaler_enabled = True
        self._scaler_vmin = 0.0
        self._scaler_vmax = 1000.0
        self._scaler_auto = True
        self._scaler_log = False
        self._scaler_auto_exclude = 5  # exclude top N hot channels from auto range
        self._palette_idx = 0
        # --- zoom / pan state ---
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._drag_last = None   # QPointF during button drag
        self._drag_origin = None # press position (for click vs drag)
        self._dragging = False
        self._layout_dirty = True  # recompute rects on next paint

        self.setMouseTracking(True)
        self.setMinimumSize(400, 400)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    # -- scan state helpers (unchanged) ----------------------------------

    def setModuleColors(self, c):  self._colors = c; self.update()
    def setPathPreview(self, p):   self._path_line = p; self.update()
    def setDashPreview(self, p):   self._dash_line = p; self.update()
    def setHighlight(self, n):     self._highlight = n; self.update()

    def setMarkerPosition(self, hx, hy):
        self._marker_hx = hx; self._marker_hy = hy; self.update()

    def modCenter(self, m):
        cx = self._ox + (m.x - self._x_min) * self._scale
        cy = self._oy + (self._y_max - m.y) * self._scale
        return QPointF(cx, cy)

    # -- scaler overlay methods ------------------------------------------

    def setScalerValues(self, vals: Dict[str, float]):
        self._scaler_values = vals
        if self._scaler_auto and vals:
            self._computeAutoRange()
        self.update()

    def setScalerEnabled(self, on: bool):
        self._scaler_enabled = on; self.update()

    def setScalerRange(self, vmin: float, vmax: float):
        self._scaler_vmin = vmin; self._scaler_vmax = vmax; self.update()

    def setScalerAutoRange(self, on: bool):
        self._scaler_auto = on
        if on and self._scaler_values:
            self._computeAutoRange()
        self.update()

    def setScalerAutoExclude(self, n: int):
        self._scaler_auto_exclude = max(0, n)
        if self._scaler_auto and self._scaler_values:
            self._computeAutoRange()
            self.update()

    def _computeAutoRange(self):
        """Auto-range excluding the top N hottest channels."""
        v = sorted(self._scaler_values.values())
        if not v: return
        n = self._scaler_auto_exclude
        trimmed = v[:max(len(v) - n, 1)]
        self._scaler_vmin = trimmed[0]
        self._scaler_vmax = trimmed[-1]
        if self._scaler_vmin == self._scaler_vmax:
            self._scaler_vmax = self._scaler_vmin + 1.0

    def setScalerLogScale(self, on: bool):
        self._scaler_log = on; self.update()

    def cyclePalette(self):
        self._palette_idx = (self._palette_idx + 1) % len(PALETTES)
        self.update()

    def scalerRange(self):
        return self._scaler_vmin, self._scaler_vmax

    def _value_to_t(self, v: float) -> float:
        vmin, vmax = self._scaler_vmin, self._scaler_vmax
        if self._scaler_log:
            floor = max(vmin, 1e-6)
            ceil = max(vmax, floor * 10)
            v = max(v, floor)
            return (math.log10(v) - math.log10(floor)) / \
                   (math.log10(ceil) - math.log10(floor))
        else:
            return (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5

    # -- layout ----------------------------------------------------------

    def _updateLayout(self):
        if not self._layout_dirty or not self._drawn: return
        self._layout_dirty = False
        x_min = min(m.x - m.sx / 2 for m in self._drawn)
        x_max = max(m.x + m.sx / 2 for m in self._drawn)
        y_min = min(m.y - m.sy / 2 for m in self._drawn)
        y_max = max(m.y + m.sy / 2 for m in self._drawn)
        w, h = self.width(), self.height()
        uw, uh = w - 2 * self.PAD, h - 2 * self.PAD
        base_scale = min(uw / max(x_max - x_min, 1), uh / max(y_max - y_min, 1))
        self._scale = base_scale * self._zoom
        dw = (x_max - x_min) * self._scale
        dh = (y_max - y_min) * self._scale
        self._ox = self.PAD + (uw - dw) / 2 + self._pan_x
        self._oy = self.PAD + (uh - dh) / 2 + self._pan_y
        self._x_min = x_min; self._y_max = y_max
        self._rects.clear()
        for m in self._drawn:
            cx = self._ox + (m.x - self._x_min) * self._scale
            cy = self._oy + (self._y_max - m.y) * self._scale
            hw = m.sx * self._scale * self.SHRINK / 2
            hh = m.sy * self._scale * self.SHRINK / 2
            self._rects[m.name] = QRectF(cx - hw, cy - hh, 2 * hw, 2 * hh)

    # -- painting --------------------------------------------------------

    # pre-built QColors for static colours used every frame
    _BG_COLOR = QColor("#0a0e14")
    _NO_DATA_COLOR = QColor("#15181d")
    _EXCLUDED_COLOR = QColor(C.MOD_EXCLUDED)

    def paintEvent(self, event):
        self._updateLayout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        p.fillRect(self.rect(), self._BG_COLOR)

        stops = list(PALETTES.values())[self._palette_idx]
        show_scaler = self._scaler_enabled and bool(self._scaler_values)
        colors = self._colors
        scaler_vals = self._scaler_values
        value_to_t = self._value_to_t
        no_data = self._NO_DATA_COLOR
        excluded = self._EXCLUDED_COLOR

        BORDER_W = 2.0
        # cache QColor objects for state colours to avoid re-creating per module
        _qcolor_cache: Dict[str, QColor] = {}

        for m in self._drawn:
            r = self._rects.get(m.name)
            if not r:
                continue
            sc_hex = colors.get(m.name)

            if show_scaler:
                sv = scaler_vals.get(m.name)
                if sv is not None:
                    p.fillRect(r, _cmap_qcolor(value_to_t(sv), stops))
                else:
                    p.fillRect(r, no_data)
                if sc_hex:
                    qc = _qcolor_cache.get(sc_hex)
                    if qc is None:
                        qc = QColor(sc_hex)
                        _qcolor_cache[sc_hex] = qc
                    p.setPen(QPen(qc, BORDER_W))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawRect(r)
            else:
                if sc_hex:
                    qc = _qcolor_cache.get(sc_hex)
                    if qc is None:
                        qc = QColor(sc_hex)
                        _qcolor_cache[sc_hex] = qc
                    p.fillRect(r, qc)
                else:
                    p.fillRect(r, excluded)

        # limit box
        bx0 = self._ox + (self._lim_hx_min - self._x_min) * self._scale
        by0 = self._oy + (self._y_max - self._lim_hy_max) * self._scale
        bx1 = self._ox + (self._lim_hx_max - self._x_min) * self._scale
        by1 = self._oy + (self._y_max - self._lim_hy_min) * self._scale
        p.setPen(QPen(QColor(C.RED), 1, Qt.PenStyle.DashLine))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(QRectF(bx0, by0, bx1 - bx0, by1 - by0))

        # path lines (white, visible over heat map)
        for pts, style, lw in [(self._path_line, Qt.PenStyle.SolidLine, 1.8),
                                (self._dash_line, Qt.PenStyle.DashLine, 1.5)]:
            if len(pts) >= 2:
                p.setPen(QPen(QColor(C.PATH_LINE), lw, style))
                for i in range(len(pts) - 1):
                    p.drawLine(pts[i], pts[i + 1])

        # highlight
        if self._highlight and self._highlight in self._rects:
            p.setPen(QPen(QColor(C.ACCENT), 2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._highlight])

        # beam position marker (red dot with white outline + crosshair)
        if self._marker_hx is not None:
            cx = self._ox + (self._marker_hx - self._x_min) * self._scale
            cy = self._oy + (self._y_max - self._marker_hy) * self._scale
            R, ARM = 3, 8
            p.setPen(QPen(QColor("#ffffff"), 2.0))
            p.setBrush(QColor(C.RED))
            p.drawEllipse(QPointF(cx, cy), R, R)
            p.setPen(QPen(QColor(C.RED), 1.5))
            p.drawLine(QPointF(cx - ARM, cy), QPointF(cx + ARM, cy))
            p.drawLine(QPointF(cx, cy - ARM), QPointF(cx, cy + ARM))

        p.end()

    # -- mouse / wheel events --------------------------------------------

    _CLICK_THRESHOLD = 4  # pixels — below this is a click, above is a drag

    def mousePressEvent(self, e):
        if e.button() == Qt.MouseButton.MiddleButton:
            self.resetView(); return
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            self._drag_last = e.position()
            self._drag_origin = e.position()
            self._dragging = False

    def mouseReleaseEvent(self, e):
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            if self._dragging:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            elif e.button() == Qt.MouseButton.LeftButton:
                # short click — select module
                pos = e.position()
                for name in reversed(list(self._rects)):
                    if self._rects[name].contains(pos):
                        self.moduleClicked.emit(name); break
            self._drag_last = None
            self._drag_origin = None
            self._dragging = False

    def mouseMoveEvent(self, e):
        # drag (left or right button held)
        if self._drag_last is not None:
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
        # hover tooltip
        pos = e.position()
        found = None
        for name in reversed(list(self._rects)):
            if self._rects[name].contains(pos): found = name; break
        if found != self._hover_name:
            self._hover_name = found
            if found:
                m = self._mod_by_name.get(found)
                if m:
                    px, py = module_to_ptrans(m.x, m.y)
                    tip = f"{m.name} ({m.mod_type})\nHyCal({m.x:.1f}, {m.y:.1f})\nptrans({px:.1f}, {py:.1f})"
                    sv = self._scaler_values.get(m.name)
                    if sv is not None:
                        tip += f"\nScaler: {sv:.1f}"
                    self.setToolTip(tip)
            else:
                self.setToolTip("")

    def wheelEvent(self, e):
        factor = 1.15 if e.angleDelta().y() > 0 else 1.0 / 1.15
        new_zoom = max(0.5, min(self._zoom * factor, 20.0))
        if new_zoom == self._zoom: return
        # zoom centred on cursor: adjust pan so the point under cursor stays put
        pos = e.position()
        ratio = new_zoom / self._zoom
        self._pan_x = pos.x() + (self._pan_x - pos.x()) * ratio
        self._pan_y = pos.y() + (self._pan_y - pos.y()) * ratio
        self._zoom = new_zoom
        self._layout_dirty = True
        self.update()

    def resetView(self):
        self._zoom = 1.0; self._pan_x = 0.0; self._pan_y = 0.0
        self._layout_dirty = True; self.update()

    def resizeEvent(self, e):
        self._layout_dirty = True; super().resizeEvent(e)

    def sizeHint(self):
        from PyQt6.QtCore import QSize
        return QSize(680, 680)
