// cut_dialog.js — Waveform-Tab peak filter dialog ("Cut Settings…")
// Called from init() in viewer.js.

function initCutDialog(){
    const cutBackdrop = document.getElementById('cut-backdrop');
    const cutDialog   = document.getElementById('cut-dialog');
    if (!cutBackdrop || !cutDialog) return;

    const $ = id => document.getElementById(id);

    // Build quality-bit checkbox lists from histConfig.quality_bits.
    // Always rebuilds — cheap, and avoids subtle bugs when the bit palette
    // changes (e.g. server config reloaded between opens).
    function buildBitList(containerId, set){
        const c = $(containerId);
        if (!c) return;
        const bits = (window.histConfig && histConfig.quality_bits) || [];
        c.innerHTML = '';
        if (!bits.length) {
            const empty = document.createElement('div');
            empty.style.cssText = 'color:var(--dim);font-size:11px;font-style:italic';
            empty.textContent = '(no bits exposed by server)';
            c.appendChild(empty);
            return;
        }
        bits.forEach(d => {
            const lbl = document.createElement('label');
            const cb  = document.createElement('input');
            cb.type = 'checkbox';
            cb.dataset.bit = d.name;
            cb.checked = !!(set && set.has(d.name));
            lbl.appendChild(cb);
            lbl.appendChild(document.createTextNode(d.label || d.name));
            c.appendChild(lbl);
        });
    }

    function openCutDialog(){
        cutBackdrop.classList.add('open');
        cutDialog.classList.add('open');
        $('cut-status-msg').textContent = '';

        const f = (window.histConfig && histConfig.filter) || {};
        const fillAxis = (axis, idMin, idMax) => {
            const r = f[axis] || {};
            $(idMin).value = r.min != null ? r.min : '';
            $(idMax).value = r.max != null ? r.max : '';
        };
        fillAxis('time',     'cut-time-min',     'cut-time-max');
        fillAxis('integral', 'cut-integral-min', 'cut-integral-max');
        fillAxis('height',   'cut-height-min',   'cut-height-max');

        const qb     = f.quality_bits || {};
        buildBitList('cut-accept-list', new Set(qb.accept || []));
        buildBitList('cut-reject-list', new Set(qb.reject || []));
    }

    function closeCutDialog(){
        cutBackdrop.classList.remove('open');
        cutDialog.classList.remove('open');
    }

    function readAxis(idMin, idMax){
        const a = $(idMin).value, b = $(idMax).value;
        const out = {};
        if (a !== '') out.min = parseFloat(a);
        if (b !== '') out.max = parseFloat(b);
        return Object.keys(out).length ? out : null;
    }

    function readBitNames(containerId){
        return Array.from($(containerId).querySelectorAll('input[type="checkbox"]'))
            .filter(cb => cb.checked).map(cb => cb.dataset.bit);
    }

    function buildFilter(){
        const f = {};
        const t = readAxis('cut-time-min',     'cut-time-max');     if (t) f.time     = t;
        const i = readAxis('cut-integral-min', 'cut-integral-max'); if (i) f.integral = i;
        const h = readAxis('cut-height-min',   'cut-height-max');   if (h) f.height   = h;
        const acc = readBitNames('cut-accept-list');
        const rej = readBitNames('cut-reject-list');
        if (acc.length || rej.length) f.quality_bits = {accept: acc, reject: rej};
        return f;
    }

    // Local redraw when the *show* toggle flips: pulls overlays from
    // histConfig.filter and the cut-show state.  Server isn't touched.
    // Bypasses the histogram refresh throttle (online mode) so the toggle
    // feels responsive — without this, showHistograms() may early-return
    // for ~1s and the overlays don't update.
    function redrawAll(){
        if (typeof lastHistModule !== 'undefined') lastHistModule = '';
        if (typeof selectedModule !== 'undefined' && selectedModule
            && typeof showWaveform === 'function') {
            showWaveform(selectedModule);
        } else if (typeof redrawGeo === 'function') {
            redrawGeo();
        }
    }

    // --- wiring -----------------------------------------------------------
    $('cut-settings-btn').onclick = openCutDialog;
    $('cut-dialog-close').onclick = closeCutDialog;
    $('cut-cancel').onclick       = closeCutDialog;
    cutBackdrop.onclick           = closeCutDialog;

    document.addEventListener('keydown', e => {
        if (e.key === 'Escape' && cutDialog.classList.contains('open'))
            closeCutDialog();
    });

    $('cut-reset').onclick = () => {
        ['cut-time-min','cut-time-max','cut-integral-min','cut-integral-max',
         'cut-height-min','cut-height-max'].forEach(id => $(id).value = '');
        ['cut-accept-list','cut-reject-list'].forEach(cid => {
            $(cid).querySelectorAll('input[type="checkbox"]').forEach(cb => cb.checked = false);
        });
    };

    $('cut-apply-btn').onclick = () => {
        const body = {
            filter:        buildFilter(),
            filter_enable: $('cut-apply').checked
        };
        $('cut-status-msg').textContent = 'Saving…';
        fetch('/api/hist_config', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(body)
        }).then(r => r.json()).then(d => {
            if (d.error) {
                $('cut-status-msg').textContent = 'Error: ' + d.error;
                return;
            }
            closeCutDialog();
            // hist_config_updated WS broadcast triggers config refresh on all clients;
            // call here too in case WS is down or this client is the one editing.
            if (typeof fetchConfigAndApply === 'function') fetchConfigAndApply();
        }).catch(() => {
            $('cut-status-msg').textContent = 'Request failed';
        });
    };

    // "apply" toggle: immediate server POST.  Just flips peak_filter.enable.
    $('cut-apply').onchange = function(){
        fetch('/api/hist_config', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({filter_enable: this.checked})
        }).then(() => {
            if (typeof fetchConfigAndApply === 'function') fetchConfigAndApply();
        }).catch(() => {});
    };

    // "show" toggle: client-side overlay only — no server roundtrip.
    $('cut-show').onchange = redrawAll;
}
