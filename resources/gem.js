// gem.js — GEM detector visualization tab
//
// Left:  per-detector cluster occupancy heatmaps (2×2 grid)
// Right: tracking-efficiency cards + last-good-event ZX/ZY display
//        (HyCal-anchored 4-point line fits, see runGemEfficiency in
//         app_state.cpp).  No per-event refresh on the right panel —
//         the snapshot is server-side and only changes when a new event
//         passes the χ² + acceptance gates.

'use strict';

const GEM_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'];

let gemEffData = null;        // last /api/gem/efficiency response

// Theme-aware layout factories (read from the active THEME at call time).
function PL_GEM_OCC() {
    return {
        ...plotlyLayout(),
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor:  THEME.canvas,
        font: { color: THEME.text, size: 10 },
        margin: { l: 45, r: 10, t: 28, b: 32 },
        hovermode: 'closest',
        showlegend: false,
    };
}

function PL_GEM_EFF() {
    return {
        ...plotlyLayout(),
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor:  THEME.canvas,
        font: { color: THEME.text, size: 10 },
        margin: { l: 50, r: 12, t: 24, b: 36 },
        hovermode: 'closest',
        showlegend: false,
    };
}

// --- fetch + render ---------------------------------------------------------

function fetchGemAccum() {
    fetch('/api/gem/occupancy').then(r => r.json()).then(plotGemOccupancy).catch(() => {});
    fetch('/api/gem/efficiency').then(r => r.json()).then(updateGemEfficiency).catch(() => {});
}

// --- occupancy heatmap (left, 2x2 per-detector) ----------------------------

const GEM_OCC_IDS = ['gem-occ-0', 'gem-occ-1', 'gem-occ-2', 'gem-occ-3'];

function plotGemOccupancy(data) {
    if (!data || !data.enabled) {
        GEM_OCC_IDS.forEach(id => {
            const div = document.getElementById(id);
            if (div) div.innerHTML = '<div style="color:var(--dim);padding:20px;text-align:center">GEM not enabled</div>';
        });
        return;
    }

    const detectors = data.detectors || [];
    const total = data.total || 0;
    const scale = total > 0 ? 1.0 / total : 0;

    // Pre-compute per-detector z matrices and the global zmax so all four
    // heatmaps share one colour axis.  Sharing the scale lets the eye
    // compare absolute occupancy across GEMs, not just shape per GEM.
    const dets = GEM_OCC_IDS.map((_, detId) => detectors.find(d => d.id === detId));
    const grids = dets.map(det => {
        if (!det) return null;
        const nx = det.nx, ny = det.ny;
        const z = [];
        let local_max = 0;
        for (let iy = 0; iy < ny; iy++) {
            const row = [];
            for (let ix = 0; ix < nx; ix++) {
                const v = (det.bins[iy * nx + ix] || 0) * scale;
                row.push(v);
                if (v > local_max) local_max = v;
            }
            z.push(row);
        }
        return { det, z, local_max };
    });
    let zmax = 0;
    for (const g of grids) if (g && g.local_max > zmax) zmax = g.local_max;
    if (zmax <= 0) zmax = 1e-6;   // avoid Plotly auto-scaling to a flat plot

    // Compact per-heatmap layout: thin colourbar only on the right column
    // (cells 1 and 3), no axis titles, small title font.
    const compactMargin  = { l: 28, r: 8,  t: 18, b: 20 };
    const compactMarginR = { l: 28, r: 42, t: 18, b: 20 };

    GEM_OCC_IDS.forEach((divId, idx) => {
        const g = grids[idx];
        const onRightCol = (idx % 2) === 1;
        const showBar = onRightCol;
        const det = g && g.det;
        const titleText = det
            ? det.name + (total > 0 ? ` (${total})` : '')
            : 'GEM' + idx;
        const frameColor = GEM_COLORS[idx] || THEME.text;

        // Dashed detector frame outline — visible even when no events have
        // accumulated yet (heatmap is uniformly zero), so the active area is
        // always shown.  Plus an explicit axis range pinned to the detector
        // size so the frame doesn't shift around when events first arrive.
        const shapes = [];
        let xRange = null, yRange = null;
        if (det && det.x_size && det.y_size) {
            shapes.push({
                type: 'rect', xref: 'x', yref: 'y',
                x0: -det.x_size / 2, x1: det.x_size / 2,
                y0: -det.y_size / 2, y1: det.y_size / 2,
                line: { color: frameColor, width: 1.2, dash: 'dash' },
                fillcolor: 'rgba(0,0,0,0)',
            });
            const padX = det.x_size * 0.04, padY = det.y_size * 0.04;
            xRange = [-det.x_size / 2 - padX, det.x_size / 2 + padX];
            yRange = [-det.y_size / 2 - padY, det.y_size / 2 + padY];
        }

        const layout = Object.assign({}, PL_GEM_OCC(), {
            title: { text: titleText, font: { size: 11, color: THEME.text } },
            xaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                     ticks: 'outside', ticklen: 3,
                     range: xRange, autorange: xRange ? false : true },
            yaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                     ticks: 'outside', ticklen: 3,
                     range: yRange, autorange: yRange ? false : true },
            margin: showBar ? compactMarginR : compactMargin,
            shapes: shapes,
        });

        if (!g) {
            Plotly.react(divId,
                [{ x: [], y: [], z: [[]], type: 'heatmap' }],
                layout, { responsive: true, displayModeBar: false });
            return;
        }

        const nx = det.nx, ny = det.ny;
        const xStep = det.x_size / nx, yStep = det.y_size / ny;
        const x0 = -det.x_size / 2 + xStep / 2;
        const y0 = -det.y_size / 2 + yStep / 2;
        const xArr = Array.from({length: nx}, (_, i) => x0 + i * xStep);
        const yArr = Array.from({length: ny}, (_, i) => y0 + i * yStep);

        const trace = {
            x: xArr, y: yArr, z: g.z,
            type: 'heatmap',
            colorscale: 'Hot',
            zmin: 0, zmax: zmax,
            zauto: false,
            hovertemplate: det.name + '<br>x=%{x:.0f}<br>y=%{y:.0f}<br>rate=%{z:.4f}<extra></extra>',
            showscale: showBar,
        };
        if (showBar) {
            trace.colorbar = { thickness: 6, tickfont: { size: 8 }, tickformat: '.2f', len: 0.92 };
        }

        Plotly.react(divId, [trace], layout, { responsive: true, displayModeBar: false });
    });
}

// --- efficiency cards + snapshot view (right) ------------------------------

function updateGemEfficiency(data) {
    if (!data || !data.enabled) {
        const c = document.getElementById('gem-eff-cards');
        if (c) c.innerHTML = '<span style="color:var(--dim);grid-column:1/-1;align-self:center;text-align:center">GEM not enabled</span>';
        const info = document.getElementById('gem-eff-info');
        if (info) info.textContent = '';
        plotGemEffEmpty();
        return;
    }
    gemEffData = data;
    renderGemEffCards();
    renderGemEffSnapshot();
}

function renderGemEffCards() {
    if (!gemEffData) return;
    const root = document.getElementById('gem-eff-cards');
    if (!root) return;
    const counters = gemEffData.counters || [];
    const cfg = gemEffData.config || {};
    const minDen  = cfg.min_denom_for_eff || 0;
    const healthy = cfg.healthy || 90;
    const warning = cfg.warning || 70;
    root.innerHTML = '';
    counters.forEach(c => {
        let cls = 'gray', txt = '—', fillPct = 0;
        if (c.den >= minDen) {
            txt = c.eff_pct.toFixed(1) + '%';
            fillPct = Math.max(0, Math.min(100, c.eff_pct));
            cls = c.eff_pct >= healthy ? 'green'
                : c.eff_pct >= warning ? 'amber' : 'red';
        }
        const el = document.createElement('div');
        el.className = 'gem-eff-card ' + cls;
        // Translucent left-to-right fill behind the text — width tracks the
        // efficiency ratio so the box itself "shows" the value at a glance.
        el.style.setProperty('--fill-pct', fillPct + '%');
        const color = GEM_COLORS[c.id] || THEME.text;
        el.innerHTML =
            `<div class="name" style="color:${color}">${c.name || ('GEM' + c.id)}</div>` +
            `<div class="pct">${txt}</div>` +
            `<div class="cnt">${c.num} / ${c.den}</div>`;
        root.appendChild(el);
    });
}

function renderGemEffSnapshot() {
    const info = document.getElementById('gem-eff-info');
    if (!gemEffData) { plotGemEffEmpty(); return; }
    const snap = gemEffData.snapshot;
    if (!snap) {
        if (info) info.innerHTML = 'Waiting for matched event…';
        plotGemEffView(null);
        return;
    }
    const tests = snap.tests || [];
    if (info) {
        // χ²/dof range across all tested detectors (single value if only one).
        const chi2s = tests.filter(t => t && t.tested).map(t => +t.chi2_per_dof);
        let chi2Str = '—';
        if (chi2s.length === 1) chi2Str = chi2s[0].toFixed(2);
        else if (chi2s.length > 1) {
            const lo = Math.min(...chi2s), hi = Math.max(...chi2s);
            chi2Str = lo.toFixed(2) === hi.toFixed(2)
                ? lo.toFixed(2)
                : `${lo.toFixed(2)}–${hi.toFixed(2)}`;
        }
        // Per-detector flags, colored to match the line/star in the plot.
        const flags = tests.map((tt, i) => {
            if (!tt || !tt.tested) return '';
            const c = GEM_COLORS[i] || THEME.text;
            return `<span style="color:${c}">GEM${i}${tt.found ? '✓' : '✗'}</span>`;
        }).filter(Boolean).join(' ');
        info.innerHTML = `Event #${snap.event_id} &nbsp; χ²/dof=${chi2Str} &nbsp; ${flags}`;
    }
    plotGemEffView(snap);
}

// Empty wrapper used when /api/gem/efficiency hasn't been fetched yet.
function plotGemEffEmpty() {
    plotGemEffView(null);
}

// Compute lab-frame axis ranges from the detector geometry alone, so the GEM
// frames + HyCal z are always in view — even before any event arrives.
function gemEffViewRanges() {
    const dets = (gemEffData && gemEffData.detectors) || [];
    const hycalZ = (gemEffData && gemEffData.hycal_z) || 0;
    let xMin = +Infinity, xMax = -Infinity;
    let yMin = +Infinity, yMax = -Infinity;
    dets.forEach(d => {
        const pos = d.position || [0, 0, 0];
        if (d.x_size) {
            xMin = Math.min(xMin, pos[0] - d.x_size / 2);
            xMax = Math.max(xMax, pos[0] + d.x_size / 2);
        }
        if (d.y_size) {
            yMin = Math.min(yMin, pos[1] - d.y_size / 2);
            yMax = Math.max(yMax, pos[1] + d.y_size / 2);
        }
    });
    if (!isFinite(xMin)) { xMin = -300; xMax = 300; }
    if (!isFinite(yMin)) { yMin = -300; yMax = 300; }
    const xPad = (xMax - xMin) * 0.06;
    const yPad = (yMax - yMin) * 0.06;
    const zMax = (hycalZ > 0 ? hycalZ : 5800) * 1.05;
    return {
        xy: { x: [xMin - xPad, xMax + xPad], y: [yMin - yPad, yMax + yPad] },
        zy: { z: [-100, zMax],               y: [yMin - yPad, yMax + yPad] },
    };
}

// Always-on reference shapes: dashed GEM frame rectangles in XY, dashed
// vertical lines at each detector's z (and HyCal z) in ZY.
function gemEffViewShapes() {
    const dets = (gemEffData && gemEffData.detectors) || [];
    const hycalZ = (gemEffData && gemEffData.hycal_z) || 0;
    const shapesXY = [], shapesZY = [];
    dets.forEach(d => {
        const pos = d.position || [0, 0, 0];
        const c = GEM_COLORS[d.id] || THEME.text;
        if (d.x_size && d.y_size) {
            shapesXY.push({
                type: 'rect', xref: 'x', yref: 'y',
                x0: pos[0] - d.x_size / 2, x1: pos[0] + d.x_size / 2,
                y0: pos[1] - d.y_size / 2, y1: pos[1] + d.y_size / 2,
                line: { color: c, width: 1, dash: 'dash' },
                fillcolor: 'rgba(0,0,0,0)',
            });
        }
        if (pos[2]) {
            shapesZY.push({
                type: 'line', xref: 'x', yref: 'paper',
                x0: pos[2], x1: pos[2], y0: 0, y1: 1,
                line: { color: c, width: 1, dash: 'dash' },
            });
        }
    });
    if (hycalZ) {
        shapesZY.push({
            type: 'line', xref: 'x', yref: 'paper',
            x0: hycalZ, x1: hycalZ, y0: 0, y1: 1,
            line: { color: THEME.text, width: 1, dash: 'dash' },
        });
    }
    return { shapesXY, shapesZY };
}

// Render both XY and ZY views.  `snap` may be null — in that case the view
// only shows the GEM frame outlines (XY) and detector / HyCal z markers (ZY).
function plotGemEffView(snap) {
    const tracesXY = [], tracesZY = [];
    const hycalZ = (gemEffData && gemEffData.hycal_z) || 5800;

    if (snap) {
        // HyCal anchor — square marker
        tracesXY.push({
            x: [snap.hycal_lab[0]], y: [snap.hycal_lab[1]],
            mode: 'markers', type: 'scatter', name: 'HyCal',
            marker: { symbol: 'square', color: THEME.text, size: 11,
                      line: { color: THEME.selectBorder, width: 1 } },
            hovertemplate: 'HyCal<br>x=%{x:.1f}<br>y=%{y:.1f}<extra></extra>',
        });
        tracesZY.push({
            x: [snap.hycal_lab[2]], y: [snap.hycal_lab[1]],
            mode: 'markers', type: 'scatter', name: 'HyCal',
            marker: { symbol: 'square', color: THEME.text, size: 11,
                      line: { color: THEME.selectBorder, width: 1 } },
            hovertemplate: 'HyCal<br>z=%{x:.0f}<br>y=%{y:.1f}<extra></extra>',
        });

        const tests = snap.tests || [];
        // Dedupe GEM hits across tests by detector (seed scan is deterministic
        // so the candidate on a given detector is the same across tests).
        const hitByDet = {};
        tests.forEach(t => {
            if (!t || !t.tested) return;
            (t.fit_hits || []).forEach(h => {
                if (!(h.det in hitByDet)) hitByDet[h.det] = h;
            });
        });
        Object.values(hitByDet).forEach(h => {
            const c = GEM_COLORS[h.det] || THEME.text;
            tracesXY.push({
                x: [h.lx], y: [h.ly], mode: 'markers', type: 'scatter',
                name: 'GEM' + h.det,
                marker: { color: c, size: 8, line: { color: THEME.selectBorder, width: 1 } },
                hovertemplate: 'GEM' + h.det + '<br>x=%{x:.2f}<br>y=%{y:.2f}<extra></extra>',
            });
            tracesZY.push({
                x: [h.lz], y: [h.ly], mode: 'markers', type: 'scatter',
                name: 'GEM' + h.det,
                marker: { color: c, size: 8, line: { color: THEME.selectBorder, width: 1 } },
                hovertemplate: 'GEM' + h.det + '<br>z=%{x:.0f}<br>y=%{y:.2f}<extra></extra>',
            });
        });

        // Per-test overlays — fit line + prediction (+ found marker) per
        // tested detector R, colored to match GEM-R.
        const z0 = 0, z1 = hycalZ;
        tests.forEach((t, R) => {
            if (!t || !t.tested) return;
            const c = GEM_COLORS[R] || THEME.text;
            tracesXY.push({
                x: [t.fit.ax + t.fit.bx * z0, t.fit.ax + t.fit.bx * z1],
                y: [t.fit.ay + t.fit.by * z0, t.fit.ay + t.fit.by * z1],
                mode: 'lines', type: 'scatter', name: `Fit (test G${R})`,
                line: { color: c, width: 1, dash: 'dot' },
                opacity: 0.7, hoverinfo: 'skip',
            });
            tracesZY.push({
                x: [z0, z1],
                y: [t.fit.ay + t.fit.by * z0, t.fit.ay + t.fit.by * z1],
                mode: 'lines', type: 'scatter', name: `Fit (test G${R})`,
                line: { color: c, width: 1, dash: 'dot' },
                opacity: 0.7, hoverinfo: 'skip',
            });
            tracesXY.push({
                x: [t.predicted_lab[0]], y: [t.predicted_lab[1]],
                mode: 'markers', type: 'scatter', name: `Pred G${R}`,
                marker: { symbol: 'star', color: c, size: 12,
                          line: { color: THEME.selectBorder, width: 1 } },
                hovertemplate: `Pred GEM${R}<br>x=%{x:.2f}<br>y=%{y:.2f}<extra></extra>`,
            });
            tracesZY.push({
                x: [t.predicted_lab[2]], y: [t.predicted_lab[1]],
                mode: 'markers', type: 'scatter', name: `Pred G${R}`,
                marker: { symbol: 'star', color: c, size: 12,
                          line: { color: THEME.selectBorder, width: 1 } },
                hovertemplate: `Pred GEM${R}<br>z=%{x:.0f}<br>y=%{y:.2f}<extra></extra>`,
            });
            if (t.found) {
                tracesXY.push({
                    x: [t.found_lab[0]], y: [t.found_lab[1]],
                    mode: 'markers', type: 'scatter', name: `Found G${R}`,
                    marker: { symbol: 'circle-open', color: c, size: 14, line: { width: 2 } },
                    hovertemplate: `Found GEM${R}<br>x=%{x:.2f}<br>y=%{y:.2f}<extra></extra>`,
                });
                tracesZY.push({
                    x: [t.found_lab[2]], y: [t.found_lab[1]],
                    mode: 'markers', type: 'scatter', name: `Found G${R}`,
                    marker: { symbol: 'circle-open', color: c, size: 14, line: { width: 2 } },
                    hovertemplate: `Found GEM${R}<br>z=%{x:.0f}<br>y=%{y:.2f}<extra></extra>`,
                });
            }
        });
    }

    const ranges = gemEffViewRanges();
    const { shapesXY, shapesZY } = gemEffViewShapes();

    Plotly.react('gem-eff-xy', tracesXY, Object.assign({}, PL_GEM_EFF(), {
        title: { text: 'Front view (X–Y)', font: { size: 10, color: THEME.text } },
        xaxis: { title: 'x (mm)', gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 range: ranges.xy.x, scaleanchor: 'y', scaleratio: 1 },
        yaxis: { title: 'y (mm)', gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 range: ranges.xy.y },
        shapes: shapesXY,
    }), { responsive: true, displayModeBar: false });
    Plotly.react('gem-eff-zy', tracesZY, Object.assign({}, PL_GEM_EFF(), {
        title: { text: 'Side view (Z–Y)', font: { size: 10, color: THEME.text } },
        xaxis: { title: 'z (mm)', gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 range: ranges.zy.z },
        yaxis: { title: 'y (mm)', gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 range: ranges.zy.y },
        shapes: shapesZY,
    }), { responsive: true, displayModeBar: false });
}

// --- resize -----------------------------------------------------------------

function resizeGem() {
    GEM_OCC_IDS.forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
    ['gem-eff-xy', 'gem-eff-zy'].forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
}
