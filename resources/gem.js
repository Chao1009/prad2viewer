// gem.js — GEM detector visualization tab
//
// Two-plane display (2x1 grid):
//   Top:    Plane 1 — GEM0 (left, blue) + GEM1 (right, orange)
//   Bottom: Plane 2 — GEM2 (left, green) + GEM3 (right, red)
//
// Hits from different detectors are color-coded.
// Detector outlines shown as rectangles.

'use strict';

// --- configuration ----------------------------------------------------------
const GEM_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'];  // blue, orange, green, red
const GEM_PLANES = [
    { name: 'Plane 1 (upstream)', dets: [0, 1], plotId: 'gem-plane-0' },
    { name: 'Plane 2 (downstream)', dets: [2, 3], plotId: 'gem-plane-1' },
];

const PL_GEM = {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: '#1a1a2e',
    font: { color: '#e0e0e0', size: 11 },
    margin: { l: 50, r: 20, t: 30, b: 40 },
    hovermode: 'closest',
};

let gemConfig = null;  // cached from /api/gem/config

// --- fetch + render ---------------------------------------------------------

function fetchGemData() {
    // fetch config once, then hits
    const configReady = gemConfig
        ? Promise.resolve(gemConfig)
        : fetch('/api/gem/config').then(r => r.json()).then(cfg => { gemConfig = cfg; return cfg; });

    configReady.then(() => {
        fetch('/api/gem/hits').then(r => r.json()).then(data => {
            plotGemPlanes(data);
        }).catch(() => {});
    });
}

function gemDetSize(detId) {
    // compute detector size in mm from config
    if (!gemConfig || !gemConfig.layers) return { xSize: 614.4, ySize: 1228.8 };
    const layer = gemConfig.layers.find(l => l.id === detId);
    if (!layer) return { xSize: 614.4, ySize: 1228.8 };
    return {
        xSize: layer.x_apvs * 128 * layer.x_pitch,
        ySize: layer.y_apvs * 128 * layer.y_pitch,
    };
}

function plotGemPlanes(data) {
    if (!data || !data.enabled) {
        GEM_PLANES.forEach(plane => {
            const div = document.getElementById(plane.plotId);
            if (div) div.innerHTML = '<div style="color:var(--dim);padding:40px;text-align:center">GEM not enabled</div>';
        });
        return;
    }

    const detectors = data.detectors || [];

    GEM_PLANES.forEach((plane, pi) => {
        const traces = [];
        const shapes = [];

        plane.dets.forEach((detId, di) => {
            const det = detectors.find(d => d.id === detId);
            if (!det) return;

            const hits = det.hits_2d || [];
            const color = GEM_COLORS[detId] || '#888';
            const detName = det.name || ('GEM' + detId);

            // scatter trace for 2D hits
            traces.push({
                x: hits.map(h => h.x),
                y: hits.map(h => h.y),
                mode: 'markers',
                type: 'scatter',
                name: detName,
                marker: {
                    color: color,
                    size: 6,
                    opacity: 0.8,
                    line: { width: 0.5, color: '#fff' },
                },
                hovertemplate: detName + '<br>x=%{x:.1f} mm<br>y=%{y:.1f} mm<extra></extra>',
            });

            // detector outline as a shape
            const sz = gemDetSize(detId);
            shapes.push({
                type: 'rect',
                x0: 0, y0: -sz.ySize / 2,
                x1: sz.xSize, y1: sz.ySize / 2,
                line: { color: color, width: 1.5, dash: 'dot' },
                fillcolor: 'rgba(0,0,0,0)',
            });
        });

        // if no hits, add an empty trace so the plot still renders
        if (traces.length === 0) {
            traces.push({
                x: [], y: [], mode: 'markers', type: 'scatter',
                name: 'No data', marker: { size: 0 },
            });
        }

        const layout = Object.assign({}, PL_GEM, {
            title: { text: plane.name, font: { size: 13, color: '#e0e0e0' } },
            xaxis: {
                title: 'X (mm)',
                gridcolor: '#333',
                zerolinecolor: '#555',
                scaleanchor: 'y',
                scaleratio: 1,
            },
            yaxis: {
                title: 'Y (mm)',
                gridcolor: '#333',
                zerolinecolor: '#555',
            },
            shapes: shapes,
            showlegend: true,
            legend: { x: 0.01, y: 0.99, bgcolor: 'rgba(0,0,0,0.3)', font: { size: 10 } },
        });

        Plotly.react(plane.plotId, traces, layout, { responsive: true, displayModeBar: false });
    });
}

function resizeGem() {
    GEM_PLANES.forEach(plane => {
        try { Plotly.Plots.resize(plane.plotId); } catch (e) {}
    });
}
