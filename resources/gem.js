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
let gemOccupancyData = null;  // last /api/gem/occupancy response — cached for theme flips

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
    return Promise.all([
        fetch('/api/gem/occupancy').then(r => r.json()).then(d => {
            gemOccupancyData = d;
            plotGemOccupancy(d);
        }).catch(() => {}),
        fetch('/api/gem/efficiency').then(r => r.json()).then(updateGemEfficiency).catch(() => {}),
    ]);
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

        // scaleanchor keeps mm in x and y at the same screen scale, so the
        // detector frame is drawn at its true geometric ratio (matches the
        // efficiency-grid panel on the right).
        const layout = Object.assign({}, PL_GEM_OCC(), {
            title: { text: titleText, font: { size: 11, color: THEME.text } },
            xaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                     ticks: 'outside', ticklen: 3,
                     range: xRange, autorange: xRange ? false : true,
                     constrain: 'domain' },
            yaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                     ticks: 'outside', ticklen: 3,
                     range: yRange, autorange: yRange ? false : true,
                     scaleanchor: 'x', scaleratio: 1, constrain: 'domain' },
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
        plotGemEffGrid(null);
        plotGemZTargetHist(null);
        return;
    }
    gemEffData = data;
    renderGemEffCards();
    renderGemEffSnapshot();
    plotGemEffGrid(data);
    plotGemZTargetHist(data.z_target_hist);
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
    if (info) {
        const chi2 = (typeof snap.chi2_per_dof === 'number')
            ? snap.chi2_per_dof.toFixed(2) : '—';
        // Per-detector ✓/✗ flag — ✓ = detector contributed to the fit (its hit
        // is within match_window of the seed line), ✗ = no in-window hit.
        // Missing detectors are dimmed so the present ones read at a glance.
        const flags = (snap.dets || []).map((d, i) => {
            const c = GEM_COLORS[i] || THEME.text;
            const ok = d && d.used_in_fit;
            const style = ok
                ? `color:${c}`
                : `color:${c};opacity:0.35`;
            return `<span style="${style}">GEM${i}${ok ? '✓' : '✗'}</span>`;
        }).join(' ');
        // Projected target z = closest approach of the fit line to the lab
        // z-axis (server-computed, only present when (bx²+by²) > 0).
        let zt = '';
        if (typeof snap.z_target_offset === 'number') {
            const sign = snap.z_target_offset >= 0 ? '+' : '−';
            zt = ` &nbsp; z<sub>t</sub>=${sign}${Math.abs(snap.z_target_offset).toFixed(1)} mm`;
        }
        info.innerHTML = `Event #${snap.event_id} &nbsp; χ²/dof=${chi2}${zt} &nbsp; ${flags}`;
    }
    plotGemEffView(snap);
    plotGemZTargetHist(gemEffData.z_target_hist);
}

// Empty wrapper used when /api/gem/efficiency hasn't been fetched yet.
function plotGemEffEmpty() {
    plotGemEffView(null);
}

// Compute lab-frame Z-Y axis ranges from the detector geometry alone, so
// the side view always shows every GEM plane + HyCal z, even before any
// event arrives.
function gemEffViewRanges() {
    const dets = (gemEffData && gemEffData.detectors) || [];
    const hycalZ = (gemEffData && gemEffData.hycal_z) || 0;
    let yMin = +Infinity, yMax = -Infinity;
    dets.forEach(d => {
        const pos = d.position || [0, 0, 0];
        if (d.y_size) {
            yMin = Math.min(yMin, pos[1] - d.y_size / 2);
            yMax = Math.max(yMax, pos[1] + d.y_size / 2);
        }
    });
    if (!isFinite(yMin)) { yMin = -300; yMax = 300; }
    const yPad = (yMax - yMin) * 0.06;
    const zMax = (hycalZ > 0 ? hycalZ : 5800) * 1.05;
    return { zy: { z: [-100, zMax], y: [yMin - yPad, yMax + yPad] } };
}

// Always-on reference shapes for the side view: dashed vertical lines at
// each detector's z and at HyCal z.
function gemEffViewShapes() {
    const dets = (gemEffData && gemEffData.detectors) || [];
    const hycalZ = (gemEffData && gemEffData.hycal_z) || 0;
    const shapesZY = [];
    dets.forEach(d => {
        const pos = d.position || [0, 0, 0];
        const c = GEM_COLORS[d.id] || THEME.text;
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
    return { shapesZY };
}

// Render the Z-Y side view of the latest matched event.  `snap` may be
// null — in that case the panel only shows the detector / HyCal z guides.
function plotGemEffView(snap) {
    const tracesZY = [];
    const hycalZ = (gemEffData && gemEffData.hycal_z) || 5800;

    if (snap) {
        // HyCal anchor — square marker
        tracesZY.push({
            x: [snap.hycal_lab[2]], y: [snap.hycal_lab[1]],
            mode: 'markers', type: 'scatter', name: 'HyCal',
            marker: { symbol: 'square', color: THEME.text, size: 11,
                      line: { color: THEME.selectBorder, width: 1 } },
            hovertemplate: 'HyCal<br>z=%{x:.0f}<br>y=%{y:.1f}<extra></extra>',
        });

        // Single fit line through the good track (HyCal + matched GEMs).
        // Dotted line from z=0 to HyCal z, drawn in theme text color.
        const fit = snap.fit || {};
        const z0 = 0, z1 = hycalZ;
        tracesZY.push({
            x: [z0, z1],
            y: [fit.ay + fit.by * z0, fit.ay + fit.by * z1],
            mode: 'lines', type: 'scatter', name: 'Fit',
            line: { color: THEME.text, width: 1.2, dash: 'dot' },
            opacity: 0.8, hoverinfo: 'skip',
        });

        // Per-detector overlays: filled circle at the hit, star at the
        // prediction.  Both are gated by hit_present so a no-hit (dimmed)
        // detector leaves the Z–Y plane blank, matching the dim/active
        // state shown by the GEM✓/✗ flags above the plot.  When a hit is
        // present we draw both markers so the residual is visible.
        (snap.dets || []).forEach(d => {
            if (!d.hit_present) return;
            const R = d.id;
            const c = GEM_COLORS[R] || THEME.text;
            if (d.hit_lab) {
                tracesZY.push({
                    x: [d.hit_lab[2]], y: [d.hit_lab[1]],
                    mode: 'markers', type: 'scatter', name: 'GEM' + R,
                    marker: { color: c, size: 8, line: { color: THEME.selectBorder, width: 1 } },
                    hovertemplate: 'GEM' + R + ' hit<br>z=%{x:.0f}<br>y=%{y:.2f}<extra></extra>',
                });
            }
            if (d.predicted_lab) {
                tracesZY.push({
                    x: [d.predicted_lab[2]], y: [d.predicted_lab[1]],
                    mode: 'markers', type: 'scatter', name: 'Pred G' + R,
                    marker: { symbol: 'star', color: c, size: 12,
                              line: { color: THEME.selectBorder, width: 1 } },
                    hovertemplate: `Pred GEM${R}<br>z=%{x:.0f}<br>y=%{y:.2f}<extra></extra>`,
                });
            }
        });
    }

    const ranges = gemEffViewRanges();
    const { shapesZY } = gemEffViewShapes();
    // Vertical dashed guide at the inferred vertex z.
    if (snap && typeof snap.z_target_lab === 'number') {
        shapesZY.push({
            type: 'line', xref: 'x', yref: 'paper',
            x0: snap.z_target_lab, x1: snap.z_target_lab, y0: 0, y1: 1,
            line: { color: THEME.text, width: 1.2, dash: 'dot' },
        });
    }

    Plotly.react('gem-eff-zy', tracesZY, Object.assign({}, PL_GEM_EFF(), {
        title: { text: 'Side view (Z–Y)', font: { size: 10, color: THEME.text } },
        xaxis: { title: 'z (mm)', gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 range: ranges.zy.z },
        yaxis: { title: 'y (mm)', gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 range: ranges.zy.y },
        shapes: shapesZY,
    }), { responsive: true, displayModeBar: false });
}

// --- per-detector efficiency-vs-position grid (left of the side view) -------
// Four heatmaps in a 2x2 layout, one per GEM, showing num/den efficiency over
// detector-local (x, y).  Bins with zero denominator are masked (rendered as
// the canvas color) so empty cells don't bias the eye.  Color scale is fixed
// to [0, 1] so cards comparing tiers stay consistent across detectors.
const GEM_EFF_GRID_IDS = ['gem-eff-grid-0', 'gem-eff-grid-1',
                          'gem-eff-grid-2', 'gem-eff-grid-3'];

function plotGemEffGrid(data) {
    const compactMargin  = { l: 28, r: 8,  t: 18, b: 20 };
    const compactMarginR = { l: 28, r: 42, t: 18, b: 20 };

    if (!data || !data.enabled) {
        GEM_EFF_GRID_IDS.forEach(id => {
            const div = document.getElementById(id);
            if (div) div.innerHTML = '<div style="color:var(--dim);padding:20px;text-align:center">GEM not enabled</div>';
        });
        return;
    }

    const detectors = data.detectors || [];
    const dets = GEM_EFF_GRID_IDS.map((_, detId) =>
        detectors.find(d => d.id === detId));

    GEM_EFF_GRID_IDS.forEach((divId, idx) => {
        const det = dets[idx];
        const onRightCol = (idx % 2) === 1;
        const showBar = onRightCol;
        const frameColor = GEM_COLORS[idx] || THEME.text;
        const detName = det && det.name ? det.name : ('GEM' + idx);

        // Always-on dashed detector frame so the active area is visible
        // even before any event arrives.  Range is pinned to the frame.
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

        // scaleanchor keeps mm in x and y at the same screen scale so the
        // detector frame is drawn at its true geometric ratio.  GEM frames
        // are ~1:2 portrait, so a wider cell will leave horizontal margin
        // around the data area — that's the correct trade-off for accuracy.
        const layout = Object.assign({}, PL_GEM_EFF(), {
            title: { text: detName, font: { size: 11, color: frameColor } },
            xaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                     ticks: 'outside', ticklen: 3,
                     range: xRange, autorange: xRange ? false : true,
                     constrain: 'domain' },
            yaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                     ticks: 'outside', ticklen: 3,
                     range: yRange, autorange: yRange ? false : true,
                     scaleanchor: 'x', scaleratio: 1, constrain: 'domain' },
            margin: showBar ? compactMarginR : compactMargin,
            shapes: shapes,
        });

        const grid = det && det.eff_grid;
        if (!grid || !grid.nx || !grid.ny || !grid.den || !grid.num) {
            Plotly.react(divId,
                [{ x: [], y: [], z: [[]], type: 'heatmap' }],
                layout, { responsive: true, displayModeBar: false });
            return;
        }

        const nx = grid.nx, ny = grid.ny;
        const xSize = grid.x_size || (det && det.x_size) || 0;
        const ySize = grid.y_size || (det && det.y_size) || 0;
        const xStep = xSize / nx, yStep = ySize / ny;
        const x0 = -xSize / 2 + xStep / 2;
        const y0 = -ySize / 2 + yStep / 2;
        const xArr = Array.from({length: nx}, (_, i) => x0 + i * xStep);
        const yArr = Array.from({length: ny}, (_, i) => y0 + i * yStep);

        // Per-bin eff = num/den.  null when den==0 so Plotly renders that
        // cell as transparent (instead of the lowest cmap color), letting
        // the canvas + dashed frame outline show through.
        const z = [];
        let totalDen = 0, totalNum = 0;
        for (let iy = 0; iy < ny; iy++) {
            const row = [];
            for (let ix = 0; ix < nx; ix++) {
                const k = iy * nx + ix;
                const den = grid.den[k] || 0;
                const num = grid.num[k] || 0;
                totalDen += den; totalNum += num;
                row.push(den > 0 ? (num / den) : null);
            }
            z.push(row);
        }
        const titleText = totalDen > 0
            ? `${detName}  (eff=${(100 * totalNum / totalDen).toFixed(1)}%, n=${totalDen})`
            : detName;
        layout.title.text = titleText;

        // Red → yellow → green gradient, fixed at [0, 1] so colors mean
        // the same thing across detectors and over time.  Tracks the
        // green/amber/red tiers used by the .gem-eff-card pills above.
        const trace = {
            x: xArr, y: yArr, z: z,
            type: 'heatmap',
            colorscale: [
                [0.0, '#d62728'],
                [0.7, '#ffbb33'],
                [0.9, '#2ca02c'],
                [1.0, '#2ca02c'],
            ],
            zmin: 0, zmax: 1,
            zauto: false,
            hoverongaps: false,
            hovertemplate: detName + '<br>x=%{x:.0f}<br>y=%{y:.0f}<br>eff=%{z:.2f}<extra></extra>',
            showscale: showBar,
        };
        if (showBar) {
            trace.colorbar = { thickness: 6, tickfont: { size: 8 },
                               tickformat: '.1f', len: 0.92 };
        }

        Plotly.react(divId, [trace], layout,
                     { responsive: true, displayModeBar: false });
    });
}

// --- projected-target-z histogram (right of the 2×2 efficiency cards) ------
// Shows DOCA-to-z-axis minus target_z, accumulated server-side; the current
// snapshot's value is drawn as a dotted vertical guide line.
function plotGemZTargetHist(hist) {
    const div = document.getElementById('gem-eff-zhist');
    if (!div) return;
    const layout = Object.assign({}, PL_GEM_EFF(), {
        title: { text: 'Projected vertex z − target z (mm)',
                 font: { size: 10, color: THEME.text } },
        margin: { l: 42, r: 8, t: 22, b: 32 },
        xaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 ticks: 'outside', ticklen: 3 },
        // Counts are non-negative — anchor the y-axis at 0 so the empty
        // histogram still shows a [0, …] range instead of dipping below.
        yaxis: { gridcolor: THEME.grid, zerolinecolor: THEME.border,
                 ticks: 'outside', ticklen: 3, rangemode: 'nonnegative' },
        bargap: 0.05,
    });
    if (!hist || !hist.bins || !hist.bins.length) {
        Plotly.react(div, [], layout, { responsive: true, displayModeBar: false });
        return;
    }
    const min = hist.min, step = hist.step;
    const x = hist.bins.map((_, i) => min + (i + 0.5) * step);
    const trace = {
        x: x, y: hist.bins, type: 'bar',
        marker: { color: THEME.accent || '#3aa0ff', line: { width: 0 } },
        hovertemplate: 'z=%{x:.1f} mm<br>count=%{y}<extra></extra>',
    };
    layout.xaxis.range = [hist.min, hist.max];
    const snap = gemEffData && gemEffData.snapshot;
    if (snap && typeof snap.z_target_offset === 'number') {
        layout.shapes = [{
            type: 'line', xref: 'x', yref: 'paper',
            x0: snap.z_target_offset, x1: snap.z_target_offset,
            y0: 0, y1: 1,
            line: { color: THEME.text, width: 1.2, dash: 'dot' },
        }];
    }
    Plotly.react(div, [trace], layout, { responsive: true, displayModeBar: false });
}

// --- resize -----------------------------------------------------------------

function resizeGem() {
    GEM_OCC_IDS.forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
    GEM_EFF_GRID_IDS.forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
    ['gem-eff-zy', 'gem-eff-zhist'].forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
}

// Theme flip — every GEM plot embeds THEME values in titles, frame outlines,
// fit lines, and marker/edge colors at draw time.  Replay both occupancy
// (from cached /api/gem/occupancy) and efficiency (from gemEffData) so the
// new theme reaches every text/marker, not just the chrome.
if (typeof onThemeChange === 'function') {
    onThemeChange(() => {
        if (gemOccupancyData) plotGemOccupancy(gemOccupancyData);
        if (gemEffData) {
            renderGemEffCards();
            renderGemEffSnapshot();
            plotGemEffGrid(gemEffData);
        } else {
            plotGemEffEmpty();
            plotGemEffGrid(null);
            plotGemZTargetHist(null);
        }
    });
}
