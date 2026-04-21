// gem.js — GEM detector visualization tab
//
// Left column:  per-event 2D cluster scatter (two planes stacked)
// Right column: accumulated cluster occupancy heatmaps
//
// Hit coordinates from the backend are centered: (0,0) = beam center.

'use strict';

// --- configuration ----------------------------------------------------------
const GEM_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'];
const GEM_PLANES = [
    { name: 'Plane 1 (upstream)',   dets: [0, 1], hitId: 'gem-plane-0' },
    { name: 'Plane 2 (downstream)', dets: [2, 3], hitId: 'gem-plane-1' },
];

const PL_GEM = {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: '#1a1a2e',
    font: { color: '#e0e0e0', size: 11 },
    margin: { l: 50, r: 20, t: 30, b: 40 },
    hovermode: 'closest',
    xaxis: { gridcolor: '#333', zerolinecolor: '#555' },
    yaxis: { gridcolor: '#333', zerolinecolor: '#555' },
};

const PL_GEM_OCC = {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: '#1a1a2e',
    font: { color: '#e0e0e0', size: 10 },
    margin: { l: 45, r: 10, t: 28, b: 32 },
    hovermode: 'closest',
    showlegend: false,
    xaxis: { gridcolor: '#333', zerolinecolor: '#555' },
    yaxis: { gridcolor: '#333', zerolinecolor: '#555' },
};

let gemConfig = null;

// --- helpers ----------------------------------------------------------------

function gemDetInfo(detId) {
    const def = { xSize: 614.4, ySize: 512.0, xOff: 0, yOff: 0 };
    if (!gemConfig || !gemConfig.layers) return def;
    const layer = gemConfig.layers.find(l => l.id === detId);
    if (!layer) return def;
    const pos = layer.position || [0, 0, 0];
    return {
        xSize: layer.x_size || layer.x_apvs * 128 * layer.x_pitch,
        ySize: layer.y_size || layer.y_apvs * 128 * layer.y_pitch,
        xOff:  pos[0] || 0,
        yOff:  pos[1] || 0,
    };
}

// --- fetch + render ---------------------------------------------------------

function fetchGemData() {
    const configReady = gemConfig
        ? Promise.resolve(gemConfig)
        : fetch('/api/gem/config').then(r => r.json()).then(cfg => { gemConfig = cfg; return cfg; });

    configReady.then(() => {
        fetch('/api/gem/hits').then(r => r.json()).then(plotGemHits).catch(() => {});
    });
}

function fetchGemAccum() {
    fetch('/api/gem/occupancy').then(r => r.json()).then(plotGemOccupancy).catch(() => {});
    fetch('/api/gem/hist').then(r => r.json()).then(plotGemHist).catch(() => {});
}

// --- event cluster scatter (left) -------------------------------------------

function plotGemHits(data) {
    if (!data || !data.enabled) {
        GEM_PLANES.forEach(plane => {
            const div = document.getElementById(plane.hitId);
            if (div) div.innerHTML = '<div style="color:var(--dim);padding:40px;text-align:center">GEM not enabled</div>';
        });
        return;
    }

    const detectors = data.detectors || [];

    GEM_PLANES.forEach((plane) => {
        const traces = [];
        const shapes = [];

        plane.dets.forEach((detId) => {
            const det = detectors.find(d => d.id === detId);
            if (!det) return;

            const hits = det.hits_2d || [];
            const color = GEM_COLORS[detId] || '#888';
            const detName = det.name || ('GEM' + detId);

            traces.push({
                x: hits.map(h => h.x),
                y: hits.map(h => h.y),
                mode: 'markers',
                type: 'scatter',
                name: detName,
                marker: {
                    color: color, size: 6, opacity: 0.8,
                    line: { width: 0.5, color: '#fff' },
                },
                hovertemplate: detName + '<br>x=%{x:.1f} mm<br>y=%{y:.1f} mm<extra></extra>',
            });

            // detector outline — offset to lab frame position
            const info = gemDetInfo(detId);
            shapes.push({
                type: 'rect',
                x0: info.xOff - info.xSize / 2, y0: info.yOff - info.ySize / 2,
                x1: info.xOff + info.xSize / 2, y1: info.yOff + info.ySize / 2,
                line: { color: color, width: 1.5, dash: 'dot' },
                fillcolor: 'rgba(0,0,0,0)',
            });
        });

        if (traces.length === 0) {
            traces.push({ x: [], y: [], mode: 'markers', type: 'scatter',
                          name: 'No data', marker: { size: 0 } });
        }

        const layout = Object.assign({}, PL_GEM, {
            title: { text: plane.name, font: { size: 13, color: '#e0e0e0' } },
            xaxis: {
                title: 'X (mm)', gridcolor: '#333', zerolinecolor: '#555',
                scaleanchor: 'y', scaleratio: 1,
            },
            yaxis: { title: 'Y (mm)', gridcolor: '#333', zerolinecolor: '#555' },
            shapes: shapes,
            showlegend: true,
            legend: { x: 0.01, y: 0.99, bgcolor: 'rgba(0,0,0,0.3)', font: { size: 10 } },
        });

        Plotly.react(plane.hitId, traces, layout, { responsive: true, displayModeBar: false });
    });
}

// --- occupancy heatmap (right, 2x2 per-detector) ---------------------------

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

    GEM_OCC_IDS.forEach((divId, detId) => {
        const det = detectors.find(d => d.id === detId);
        if (!det) {
            Plotly.react(divId,
                [{ x: [], y: [], z: [[]], type: 'heatmap' }],
                Object.assign({}, PL_GEM_OCC, { title: { text: 'GEM' + detId, font: { size: 12, color: '#e0e0e0' } } }),
                { responsive: true, displayModeBar: false });
            return;
        }

        const nx = det.nx, ny = det.ny;
        const xSize = det.x_size, ySize = det.y_size;
        const xStep = xSize / nx, yStep = ySize / ny;

        // build z matrix as rate (count / total_events)
        const z = [];
        const scale = total > 0 ? 1.0 / total : 0;
        for (let iy = 0; iy < ny; iy++) {
            const row = [];
            for (let ix = 0; ix < nx; ix++)
                row.push((det.bins[iy * nx + ix] || 0) * scale);
            z.push(row);
        }

        const x0 = -xSize / 2 + xStep / 2;
        const y0 = -ySize / 2 + yStep / 2;
        const xArr = Array.from({length: nx}, (_, i) => x0 + i * xStep);
        const yArr = Array.from({length: ny}, (_, i) => y0 + i * yStep);

        const traces = [{
            x: xArr, y: yArr, z: z,
            type: 'heatmap',
            colorscale: 'Hot', reversescale: true,
            hovertemplate: det.name + '<br>x=%{x:.0f} mm<br>y=%{y:.0f} mm<br>rate=%{z:.4f}<extra></extra>',
            colorbar: { thickness: 10, tickfont: { size: 9 }, tickformat: '.3f' },
        }];

        const layout = Object.assign({}, PL_GEM_OCC, {
            title: { text: det.name + (total > 0 ? ' (' + total + ' evts)' : ''),
                     font: { size: 12, color: '#e0e0e0' } },
            xaxis: { title: 'X (mm)', gridcolor: '#333', zerolinecolor: '#555' },
            yaxis: { title: 'Y (mm)', gridcolor: '#333', zerolinecolor: '#555' },
        });

        Plotly.react(divId, traces, layout, { responsive: true, displayModeBar: false });
    });
}

// --- GEM histograms (bottom right) ------------------------------------------

const GEM_HIST_IDS = ['gem-ncl-hist', 'gem-theta-hist'];
let currentGemNclHist = null, currentGemThetaHist = null;

function plotGemHist(data) {
    if (!data) { currentGemNclHist = currentGemThetaHist = null; return; }

    function plotOne(divId, hdata, title, xlabel, color, refKey) {
        if (!hdata || !hdata.bins || hdata.bins.length === 0) {
            Plotly.react(divId, [], Object.assign({}, PL_GEM_OCC, {
                title: { text: title, font: { size: 12, color: '#e0e0e0' } },
            }), { responsive: true, displayModeBar: false });
            return null;
        }
        const n = hdata.bins.length;
        const x = Array.from({length: n}, (_, i) => hdata.min + (i + 0.5) * hdata.step);
        const entries = hdata.bins.reduce((a, b) => a + b, 0);
        // store non-zero bins for copy
        const cx = [], cy = [];
        for (let i = 0; i < n; i++) { if (hdata.bins[i] > 0) { cx.push(x[i]); cy.push(hdata.bins[i]); } }
        Plotly.react(divId, [{
            x: x, y: hdata.bins, type: 'bar',
            marker: { color: color, line: { width: 0 } },
            hovertemplate: xlabel + '=%{x:.1f}<br>count=%{y}<extra></extra>',
        }], Object.assign({}, PL_GEM_OCC, {
            title: { text: title + '<br><span style="font-size:9px;color:#888">' + entries + ' entries</span>',
                     font: { size: 12, color: '#e0e0e0' } },
            xaxis: { title: xlabel, gridcolor: '#333', zerolinecolor: '#555' },
            yaxis: { title: 'Counts', gridcolor: '#333', zerolinecolor: '#555' },
            bargap: 0.05,
            shapes: refKey ? refShapes(refKey) : [],
        }), { responsive: true, displayModeBar: false });
        return { x: cx, y: cy };
    }

    currentGemNclHist = plotOne('gem-ncl-hist', data.nclusters, 'GEM Clusters / Event', 'N clusters', '#51cf66', 'gem_clusters');
    currentGemThetaHist = plotOne('gem-theta-hist', data.theta, 'GEM Hit Angle', 'θ (deg)', '#00b4d8', 'gem_theta');
}

// --- resize -----------------------------------------------------------------

function resizeGem() {
    GEM_PLANES.forEach(plane => {
        try { Plotly.Plots.resize(plane.hitId); } catch (e) {}
    });
    GEM_OCC_IDS.forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
    GEM_HIST_IDS.forEach(id => {
        try { Plotly.Plots.resize(id); } catch (e) {}
    });
}
