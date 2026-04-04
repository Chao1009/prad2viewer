// viewer.js — Orchestrator: init, tabs, event navigation, WebSocket, mode switching
// =========================================================================
// State
// =========================================================================
let modules=[], totalEvents=0, currentEvent=1;
let currentEventNumber=0, currentTriggerBits=0;  // DAQ event number + trigger from last loaded event
let triggerBitsDef=[];  // [{bit, mask, name, label}, ...]
let eventChannels={};
let selectedModule=null, hoveredModule=null;
const PC=['#00b4d8','#ff6b6b','#51cf66','#ffd43b','#cc5de8','#ff922b','#20c997','#f06595'];
const CRATE_NAME={0x80:'adchycal1',0x82:'adchycal2',0x84:'adchycal3',0x86:'adchycal4',0x88:'adchycal5',0x8a:'adchycal6',0x8c:'adchycal7',
    0x01:'PRadTS',0x04:'PRadROC_1',0x05:'PRadROC_2',0x06:'PRadROC_3',0x07:'PRadSRS_1',0x08:'PRadSRS_2'};
function crateName(r){return CRATE_NAME[r]||`ROC 0x${r.toString(16)}`;}
let histEnabled=false, histConfig={};
let mode='idle';    // 'idle', 'file', or 'online'
let etAvailable=false, fileAvailable=false;
let sourceCaps={type:'',has_waveforms:true,has_peaks:true,has_clusters:true,has_epics:true};
let ws=null;        // WebSocket connection (always connected)
let autoFollow=true; // auto-load latest event
let lastEventFetch=0, lastHistFetch=0, lastRingFetch=0, lastOccFetch=0, lastLmsFetch=0;
let refreshEventMs=200, refreshRingMs=500, refreshHistMs=2000, refreshLmsMs=2000;

// occupancy data (fetched once per file load when histograms enabled)
let occData={}, occTcutData={}, occTotal=0;

let activeTab='dq';  // 'dq' or 'cluster'

// =========================================================================
// Plotly shared config
// =========================================================================
const PL={paper_bgcolor:'#1a1a2e',plot_bgcolor:'#11112a',font:{family:'Consolas,monospace',size:10,color:'#aaa'},
    margin:{l:45,r:10,t:24,b:32},xaxis:{gridcolor:'#1a1a3a',zerolinecolor:'#222'},
    yaxis:{gridcolor:'#1a1a3a',zerolinecolor:'#222'}};
const PC2={responsive:true,displayModeBar:false};
const PC_EPICS={responsive:true,displayModeBar:true,
    modeBarButtonsToRemove:['sendDataToCloud','lasso2d','select2d'],
    displaylogo:false};

// ── Plotly plot registry ──────────────────────────────────────────────
// All Plotly divs register here with their tab and default layout.
// Provides unified init, resize-by-tab, and resize-all.
const plotRegistry=[];  // [{id, tab, layout, config}]

function registerPlot(id, tab, title, config){
    const layout={...PL};
    if(title) layout.title={text:title, font:{size:10,color:'#555'}};
    plotRegistry.push({id, tab, layout, config: config||PC2});
}

function initRegisteredPlots(){
    for(const p of plotRegistry)
        Plotly.newPlot(p.id, [], p.layout, p.config);
}

function resizePlotsForTab(tab){
    for(const p of plotRegistry)
        if(p.tab===tab) try{Plotly.Plots.resize(p.id);}catch(e){}
}

function resizeAllPlots(){
    for(const p of plotRegistry)
        try{Plotly.Plots.resize(p.id);}catch(e){}
}

function redrawGeo(){
    if(activeTab==='cluster') geoCluster();
    else if(activeTab==='lms') geoLms();
    else geoDq();
}

function geoHandleClick(cx,cy){
    const m=hitTest(cx,cy);
    if(!m){
        // click on empty canvas — deselect
        if(activeTab==='cluster'){
            selectedCluster=-1;
            document.getElementById('cl-select').value='all';
            geoCluster(); updateClusterTable(); showClusterDetail();
        } else if(activeTab==='lms'){
            lmsSelectedModule=-1;
            currentLmsData=null;
            Plotly.react('lms-plot',[],{...PL,title:{text:'LMS History',font:{size:10,color:'#555'}}},PC2);
            document.getElementById('lms-detail-header').innerHTML=
                '<span class="cl-info-text">Click a module to view LMS history</span>';
            updateLmsTable(); geoLms();
        } else {
            selectedModule=null;
            currentWaveform=null;
            currentHist={};
            document.getElementById('detail-header').innerHTML=
                '<div class="empty-msg">Click a module to view details</div>';
            Plotly.react('waveform-div',[],{...PL,xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC'}},PC2);
            Plotly.react('heighthist-div',[],{...PL,title:{text:'Height Histogram',font:{size:10,color:'#555'}}},PC2);
            Plotly.react('inthist-div',[],{...PL,title:{text:'Integral Histogram',font:{size:10,color:'#555'}}},PC2);
            Plotly.react('poshist-div',[],{...PL,title:{text:'Position Histogram',font:{size:10,color:'#555'}}},PC2);
            document.getElementById('peaks-tbody').innerHTML='';
            geoDq();
        }
        return;
    }
    if(activeTab==='cluster'){
        selectedModule=null;
        const idx=modules.indexOf(m);
        if(clusterData && clusterData.clusters && clusterData.clusters.length){
            const clusters=clusterData.clusters;
            let found=-1;
            for(let ci=0;ci<clusters.length;ci++){
                if(clusters[ci].modules&&clusters[ci].modules.includes(idx)){ found=ci; break; }
            }
            if(found<0) selectedCluster=-1;
            else selectedCluster=(selectedCluster===found)?-1:found;
        } else {
            selectedCluster=-1;
        }
        document.getElementById('cl-select').value=selectedCluster>=0?selectedCluster:'all';
        geoCluster(); updateClusterTable(); showClusterDetail();
    } else if(activeTab==='lms'){
        const idx=modules.indexOf(m);
        lmsSelectedModule=idx;
        fetchLmsHistory(idx, m.n);
        updateLmsTable();
        geoLms();
    } else {
        showWaveform(m);
    }
}

// =========================================================================
// Event loading (works for both file and online mode)
// =========================================================================
let eventRequestId = 0;  // increments on each fetch, stale responses ignored

// Build sample label: "Sample 100 (Evt. 99)"
function sampleLabel(){
    const evn=currentEventNumber?` (Evt. ${currentEventNumber})`:'';
    return `Sample ${currentEvent}${evn}`;
}

// Update status bar based on active tab
function decodeTriggerBits(bits){
    if(!bits) return '';
    const names=[];
    for(const d of triggerBitsDef){
        if(bits & (1<<d.bit)) names.push(d.name);
    }
    // include hex and decoded names
    let s=` trig=0x${bits.toString(16)}`;
    if(names.length) s+=` [${names.join('+')}]`;
    return s;
}

function updateStatusBar(){
    const modeTag = mode === 'online' ? ' [LIVE]' : '';
    const trig = decodeTriggerBits(currentTriggerBits);
    const label = sampleLabel();

    if(activeTab==='cluster'){
        const nc=clusterData?clusterData.clusters?clusterData.clusters.length:0:0;
        const nh=clusterData?clusterData.hits?Object.keys(clusterData.hits).length:0:0;
        document.getElementById('status-bar').textContent=`${label}: ${nc} clusters, ${nh} hit modules${trig}${modeTag}`;
    } else if(activeTab==='lms'){
        const lmsN=lmsSummaryData?lmsSummaryData.events||0:0;
        document.getElementById('status-bar').textContent=`${label} | LMS: ${lmsN} events${modeTag}`;
    } else {
        const nch = Object.keys(eventChannels).length;
        const npk = Object.values(eventChannels).reduce((s,c) => s + (c.pk||[]).length, 0);
        document.getElementById('status-bar').textContent=`${label}: ${nch} channels, ${npk} peaks${trig}${modeTag}`;
    }
}

function loadEventData(reqId, data) {
    if (reqId !== eventRequestId) return;  // stale response, discard
    if (data.error) {
        document.getElementById('status-bar').textContent = data.error;
        // refresh ring selector to remove stale entries
        if (mode === 'online') updateRingSelector();
        return;
    }
    currentEvent = data.event;
    currentEventNumber = data.event_number || 0;
    currentTriggerBits = data.trigger_bits || 0;
    eventChannels = data.channels || {};
    if(mode==='online') sampleCount++;
    updateStatusBar();
    updateHeaderStats();
    if(activeTab==='cluster'){
        clusterEvent=-1; // invalidate cache
        loadClusterData(currentEvent);
    } else if(activeTab==='lms'){
        // LMS geo doesn't change per event — no redraw needed
    } else if(activeTab==='gem'){
        fetchGemData();
    } else {
        geoDq();
    }
    if (activeTab==='dq' && selectedModule) showWaveform(selectedModule);
    updateGeoTooltip();

    // refresh histograms from server (on-demand accumulation in non-preprocessed mode)
    if(mode==='file' && !histEnabled){
        fetchOccupancy(); fetchClHist();
    }
}

function loadEvent(evnum) {
    currentEvent = evnum;
    const reqId = ++eventRequestId;
    if (mode === 'file') document.getElementById('ev-input').value = evnum;
    document.getElementById('status-bar').textContent = `Loading sample ${evnum}...`;
    fetch(`/api/event/${evnum}`).then(r => r.json()).then(d => loadEventData(reqId, d))
        .catch(err => { document.getElementById('status-bar').textContent = `Error: ${err}`; });
}

function loadLatestEvent() {
    const reqId = ++eventRequestId;
    fetch('/api/event/latest').then(r => r.json()).then(d => loadEventData(reqId, d))
        .catch(err => { document.getElementById('status-bar').textContent = `Error: ${err}`; });
}

// =========================================================================
// Online mode: WebSocket + ring buffer
// =========================================================================
function updateRingSelector() {
    fetch('/api/ring').then(r => r.json()).then(data => {
        const sel = document.getElementById('ring-select');
        const prev = sel.value;
        sel.innerHTML = '';
        const ring = data.ring || [];
        for (let i = ring.length - 1; i >= 0; i--) {
            const o = document.createElement('option');
            o.value = ring[i];
            o.textContent = `Sample ${ring[i]}` + (i === ring.length - 1 ? ' (latest)' : '');
            sel.appendChild(o);
        }
        // keep selection if not auto-following
        if (!autoFollow && prev && ring.includes(parseInt(prev))) sel.value = prev;
        else if (ring.length) sel.value = ring[ring.length - 1];
    });
}

function setEtStatus(connected, waiting, retries) {
    const el = document.getElementById('et-status');
    if (connected) {
        el.textContent = '● Connected';
        el.style.color = '#51cf66';
    } else if (waiting) {
        el.textContent = `● Waiting for ET (${retries||'...'})`;
        el.style.color = '#ffd43b';
    } else {
        el.textContent = '● Disconnected';
        el.style.color = '#f66';
    }
}

function updateFollowStatus() {
    const el = document.getElementById('follow-status');
    if (!el) return;
    if (autoFollow) {
        el.style.display = 'none';
    } else {
        el.textContent = `⏸ Paused at ${sampleLabel()} — click or press F to resume`;
        el.style.display = '';
    }
}

function connectWebSocket() {
    const proto = location.protocol === 'https:' ? 'wss:' : 'ws:';
    ws = new WebSocket(`${proto}//${location.host}`);

    ws.onopen = () => {};
    ws.onclose = () => {
        setTimeout(connectWebSocket, 2000);
    };
    ws.onmessage = (evt) => {
        try {
            const msg = JSON.parse(evt.data);
            if (msg.type === 'new_event') {
                setEtStatus(true);  // receiving events means ET is connected
                const now = Date.now();
                // throttle event display to ~5 Hz
                if (autoFollow && now - lastEventFetch > refreshEventMs) {
                    lastEventFetch = now;
                    loadLatestEvent();
                }
                // throttle ring selector update to ~2 Hz
                if (now - lastRingFetch > refreshRingMs) {
                    lastRingFetch = now;
                    updateRingSelector();
                }
                // throttle occupancy + cluster hist refresh to ~0.5 Hz
                if (now - lastOccFetch > refreshHistMs) {
                    lastOccFetch = now;
                    if(histEnabled) { fetchOccupancy(); fetchClHist(); }
                    if(activeTab==='physics') fetchPhysics();
                    if(activeTab==='gem') fetchGemAccum();
                }
            } else if (msg.type === 'status') {
                setEtStatus(msg.connected, msg.waiting, msg.retries);
            } else if (msg.type === 'hist_cleared') {
                occData={}; occTcutData={}; occTotal=0;
                initClHist(); plotClHist(); plotClStatHists();
                if (selectedModule) showHistograms(selectedModule);
                clearPhysicsFrontend();
                redrawGeo();
                if(activeTab==='gem') fetchGemAccum();
            } else if (msg.type === 'lms_event') {
                // throttle LMS refresh to ~0.5 Hz
                const now2 = Date.now();
                if (!lastLmsFetch) lastLmsFetch = 0;
                if (now2 - lastLmsFetch > refreshLmsMs) {
                    lastLmsFetch = now2;
                    if(activeTab==='lms') fetchLmsSummary();
                    // also refresh selected module's history
                    if(activeTab==='lms' && lmsSelectedModule>=0){
                        const name=lmsSummaryData&&lmsSummaryData.modules&&lmsSummaryData.modules[String(lmsSelectedModule)]
                            ?lmsSummaryData.modules[String(lmsSelectedModule)].name:'';
                        fetchLmsHistory(lmsSelectedModule, name);
                    }
                }
            } else if (msg.type === 'lms_cleared') {
                lmsSummaryData=null; lmsSelectedModule=-1; currentLmsData=null;
                if(activeTab==='lms'){ geoLms(); updateLmsTable(); }
            } else if (msg.type === 'epics_event') {
                const now3 = Date.now();
                if (now3 - lastEpicsFetch > refreshEpicsMs) {
                    lastEpicsFetch = now3;
                    if(activeTab==='epics'){
                        fetchEpicsChannels();
                        fetchEpicsLatest();
                        fetchAllEpicsSlots();
                    }
                }
            } else if (msg.type === 'epics_cleared') {
                clearEpicsFrontend();
            } else if (msg.type === 'mode_changed') {
                if (msg.mode && msg.mode !== mode) {
                    clearFrontend();
                    fetchConfigAndApply();
                }
            }
        } catch (e) {}
    };
}


// =========================================================================
// File browser (lazy-loading)
// =========================================================================

function openFileDialog() {
    const hdr = document.querySelector('.file-dialog-header span');
    const list = document.getElementById('file-list');
    const filter = document.getElementById('file-filter');
    const opts = document.querySelector('.file-dialog-options');

    document.getElementById('file-dialog').classList.add('open');
    document.getElementById('file-backdrop').classList.add('open');

    if (!g_dataDirEnabled) {
        hdr.textContent = 'Open EVIO File';
        filter.style.display = 'none';
        if (opts) opts.style.display = 'none';
        list.innerHTML = '<div style="padding:20px;color:var(--dim);text-align:center">'
            + 'No data folder configured.<br>Start with <code>--data-dir /path</code></div>';
        return;
    }

    hdr.textContent = `Open EVIO File — ${g_dataDir}`;
    filter.style.display = '';
    if (opts) opts.style.display = '';
    filter.value = '';
    list.innerHTML = '';

    fetchDirEntries('', list);
    filter.focus();
}

function closeFileDialog() {
    document.getElementById('file-dialog').classList.remove('open');
    document.getElementById('file-backdrop').classList.remove('open');
}

function fetchDirEntries(dir, container) {
    const url = dir ? `/api/files?dir=${encodeURIComponent(dir)}` : '/api/files';
    fetch(url).then(r => r.json()).then(data => {
        const entries = data.entries || [];
        if (!entries.length) {
            container.innerHTML = '<div style="padding:12px;color:var(--dim);text-align:center">Empty</div>';
            return;
        }
        renderEntries(entries, container);
    });
}

function renderEntries(entries, container) {
    const currentFile = g_currentFile || '';
    for (const e of entries) {
        if (e.type === 'dir') {
            const dirName = e.name.split('/').pop() || e.name;
            const row = document.createElement('div');
            row.className = 'file-folder';
            row.innerHTML = `<span class="folder-arrow">▸</span><span>${dirName}/</span>`
                + `<span class="folder-count">(${e.count})</span>`;
            const contents = document.createElement('div');
            contents.className = 'file-folder-contents';
            let loaded = false;
            row.onclick = () => {
                const open = contents.classList.toggle('open');
                row.querySelector('.folder-arrow').textContent = open ? '▾' : '▸';
                if (open && !loaded) {
                    loaded = true;
                    fetchDirEntries(e.name, contents);
                }
            };
            container.appendChild(row);
            container.appendChild(contents);
        } else {
            const fname = e.name.split('/').pop() || e.name;
            const isCurrent = e.name === currentFile;
            const row = document.createElement('div');
            row.className = 'file-item' + (isCurrent ? ' current' : '');
            row.innerHTML = `<span>${fname}</span><span class="fsize">${e.size_mb} MB</span>`;
            row.onclick = () => { closeFileDialog(); loadNewFile(e.name); };
            container.appendChild(row);
        }
    }
}

function filterFileList(text) {
    const filt = text.toLowerCase();
    const list = document.getElementById('file-list');
    // Filter visible folder/file rows by name
    for (const el of list.querySelectorAll('.file-folder')) {
        const name = el.textContent.toLowerCase();
        const contents = el.nextElementSibling;
        const match = !filt || name.includes(filt);
        el.style.display = match ? '' : 'none';
        if (contents && contents.classList.contains('file-folder-contents'))
            contents.style.display = match && contents.classList.contains('open') ? '' : 'none';
    }
    for (const el of list.querySelectorAll('.file-item')) {
        // only filter root-level items (not inside folders)
        if (el.parentElement === list) {
            const name = el.textContent.toLowerCase();
            el.style.display = (!filt || name.includes(filt)) ? '' : 'none';
        }
    }
}

let g_currentFile = '';
let g_histCheckbox = false;
let g_dataDirEnabled = false;
let g_dataDir = '';

function loadNewFile(relpath) {
    g_histCheckbox = document.getElementById('hist-checkbox').checked;
    document.getElementById('status-bar').textContent = `Loading ${relpath}...`;
    const histParam = g_histCheckbox ? '1' : '0';
    fetch(`/api/load?file=${encodeURIComponent(relpath)}&hist=${histParam}`).then(r => r.json()).then(data => {
        if (data.error) {
            document.getElementById('status-bar').textContent = data.error;
            return;
        }
        // show progress overlay and start polling
        showProgress(relpath);
    });
}

function showProgress(filename) {
    document.getElementById('progress-overlay').classList.add('active');
    document.getElementById('progress-title').textContent = `Loading ${filename.replace(/.*\//, '')}`;
    document.getElementById('progress-bar').style.width = '0%';
    document.getElementById('progress-text').textContent = 'Starting...';
    pollProgress();
}

function pollProgress() {
    fetch('/api/progress').then(r => r.json()).then(data => {
        if (!data.loading) {
            // done — hide overlay, reload config + first event
            document.getElementById('progress-overlay').classList.remove('active');
            fetchConfigAndApply();
            return;
        }
        const pct = data.total > 0 ? Math.round(100 * data.current / data.total) : 0;
        const phaseText = data.phase === 'indexing' ? 'Indexing events' : 'Building histograms';
        document.getElementById('progress-bar').style.width = data.total > 0 ? `${Math.min(pct, 100)}%` : '100%';
        document.getElementById('progress-text').textContent = data.total > 0
            ? `${phaseText}... ${data.current} / ${data.total}`
            : `${phaseText}...`;
        setTimeout(pollProgress, 300);
    }).catch(() => setTimeout(pollProgress, 1000));
}

function fetchOccupancy() {
    fetch('/api/occupancy').then(r => r.json()).then(data => {
        occData = data.occ || {};
        occTcutData = data.occ_tcut || {};
        occTotal = data.total || 0;
        // redraw if currently showing occupancy on DQ tab
        if (activeTab === 'dq' && document.getElementById('color-metric').value === 'occupancy') {
            syncDqRange();
            geoDq();
        }
        updateGeoTooltip();
    }).catch(() => {});
}


let sampleCount=0;
function updateHeaderStats(){
    const el=document.getElementById('header-stats');
    if(mode==='online'){
        el.textContent=`${sampleCount} samples`;
    } else {
        el.textContent=`${totalEvents} events`;
    }
}

// =========================================================================
// Draggable dividers
// =========================================================================
function setupDivider(divId, axis, getTarget, getContainer, getOffset, minA, minB, onResize){
    const div=document.getElementById(divId);
    let active=false;
    div.addEventListener('mousedown',e=>{
        active=true; div.classList.add('active');
        document.body.style.cursor=axis==='x'?'col-resize':'row-resize';
        document.body.style.userSelect='none'; e.preventDefault();
    });
    document.addEventListener('mousemove',e=>{
        if(!active)return;
        const container=getContainer(), rect=container.getBoundingClientRect();
        const pos=axis==='x'?e.clientX-rect.left-getOffset():e.clientY-rect.top-getOffset();
        const max_=(axis==='x'?rect.width:rect.height)-getOffset()-minB;
        const val=Math.max(minA,Math.min(max_,pos));
        const target=getTarget();
        target.style.flex='none';  // override flex:1 so explicit size takes effect
        if(axis==='x') target.style.width=val+'px'; else target.style.height=val+'px';
        onResize();
    });
    document.addEventListener('mouseup',()=>{
        if(!active)return; active=false; div.classList.remove('active');
        document.body.style.cursor=''; document.body.style.userSelect='';
    });
}

// =========================================================================
// Tab switching
// =========================================================================
function switchTab(tab){
    if(tab===activeTab) return;
    activeTab=tab;
    selectedModule=null;
    // clear notification dot on the tab being opened
    if(tab==='lms') document.getElementById('lms-dot').className='tab-dot';
    if(tab==='epics') document.getElementById('epics-dot').className='tab-dot';
    document.querySelectorAll('.tab').forEach(t=>{
        t.classList.toggle('active', t.dataset.tab===tab);
    });
    const fullTab=tab==='epics'||tab==='physics'||tab==='gem';
    document.getElementById('geo-panel').style.display        = fullTab ? 'none' : '';
    document.getElementById('div-main').style.display         = fullTab ? 'none' : '';
    document.getElementById('geo-toolbar-dq').style.display   = tab==='dq' ? 'flex' : 'none';
    document.getElementById('geo-toolbar-cl').style.display   = tab==='cluster' ? 'flex' : 'none';
    document.getElementById('geo-toolbar-lms').style.display  = tab==='lms' ? 'flex' : 'none';
    document.getElementById('detail-panel').style.display     = tab==='dq' ? 'flex' : 'none';
    document.getElementById('cluster-panel').style.display    = tab==='cluster' ? 'flex' : 'none';
    document.getElementById('lms-panel').style.display        = tab==='lms' ? 'flex' : 'none';
    document.getElementById('epics-outer').style.display      = tab==='epics' ? 'flex' : 'none';
    document.getElementById('physics-outer').style.display    = tab==='physics' ? 'flex' : 'none';
    document.getElementById('gem-outer').style.display        = tab==='gem' ? 'flex' : 'none';

    // --- per-tab actions: fetch data + resize after layout settles ---
    const tabActions = {
        dq:      { hasGeo: true },
        cluster: { hasGeo: true,
                   fetch(){ loadClusterData(currentEvent); },
                   after(){ plotClHist(); plotClStatHists(); } },
        lms:     { hasGeo: true,
                   fetch(){ fetchLmsSummary(); } },
        epics:   { fetch(){ fetchEpicsChannels(); fetchEpicsLatest(); fetchAllEpicsSlots(); } },
        physics: { fetch(){ fetchPhysics(); } },
        gem:     { fetch(){ fetchGemData(); fetchGemAccum(); },
                   after(){ resizeGem(); } },
    };
    const action = tabActions[tab] || tabActions.dq;

    if (action.fetch) action.fetch();

    // after layout settles: resize geo (for geo tabs) + registered Plotly plots + custom after()
    setTimeout(()=>{
        if (action.hasGeo) resizeGeo();
        resizePlotsForTab(tab);
        if (action.after) action.after();
    }, 50);
    updateStatusBar();
}

// Init
// =========================================================================
function init(){
    drawColorBar(); initGeo();
    document.getElementById('colorbar-canvas').onclick=()=>{
        paletteIdx=(paletteIdx+1)%PALETTE_NAMES.length;
        drawColorBar(); redrawGeo();
    };
    registerPlot('waveform-div', 'dq', null);
    registerPlot('heighthist-div','dq', 'Height Histogram');
    registerPlot('inthist-div',  'dq', 'Integral Histogram');
    registerPlot('poshist-div',  'dq', 'Position Histogram');

    // --- copy data buttons ---
    function setupCopyBtn(btnId, getData) {
        document.getElementById(btnId).onclick=()=>{
            const d=getData();
            if(!d) return;
            const text=`x: [${d.x.join(', ')}]\ny: [${d.y.join(', ')}]`;
            navigator.clipboard.writeText(text).then(()=>{
                const btn=document.getElementById(btnId);
                btn.textContent='✓'; setTimeout(()=>{btn.textContent='copy';},1000);
            });
        };
    }
    setupCopyBtn('btn-copy-wf', ()=>currentWaveform);
    setupCopyBtn('btn-copy-heighthist', ()=>currentHist['heighthist-div']);
    setupCopyBtn('btn-copy-inthist', ()=>currentHist['inthist-div']);
    setupCopyBtn('btn-copy-poshist', ()=>currentHist['poshist-div']);

    // --- dividers ---
    // 1. main vertical: geo ↔ detail
    setupDivider('div-main','x',
        ()=>document.getElementById('geo-panel'),
        ()=>document.querySelector('.main'),
        ()=>0, 300, 350, ()=>{resizeGeo();resizeAllPlots();});
    // 2-4. horizontal dividers between the 4 DQ cells
    // Each divider resizes the cell above it; offset = top of that cell relative to container
    function dqDivider(divId, cellId, prevCellIds){
        setupDivider(divId, 'y',
            ()=>document.getElementById(cellId),
            ()=>document.getElementById('detail-panel'),
            ()=>{
                // offset = header + all cells/dividers above this cell
                let off=document.getElementById('detail-header').offsetHeight;
                for(const id of prevCellIds){
                    const el=document.getElementById(id);
                    if(el) off+=el.offsetHeight;
                }
                return off;
            },
            60, 60, resizeAllPlots);
    }
    dqDivider('div-dq-1','dq-cell-heighthist',[]);
    dqDivider('div-dq-2','dq-cell-histrow',['dq-cell-heighthist','div-dq-1']);
    dqDivider('div-dq-3','dq-cell-waveform',['dq-cell-heighthist','div-dq-1','dq-cell-histrow','div-dq-2']);

    // --- tab switching ---
    document.querySelectorAll('.tab').forEach(t=>{
        t.onclick=()=>switchTab(t.dataset.tab);
    });

    // --- cluster controls ---
    document.getElementById('cl-select').onchange=e=>{
        selectedCluster=e.target.value==='all'?-1:parseInt(e.target.value);
        geoCluster(); updateClusterTable(); showClusterDetail();
    };
    document.getElementById('cl-log-scale').onchange=()=>{ if(activeTab==='cluster') geoCluster(); };
    document.getElementById('cl-colorbar-canvas').onclick=()=>{
        paletteIdx=(paletteIdx+1)%PALETTE_NAMES.length;
        drawColorBar(); redrawGeo();
    };
    document.getElementById('lms-colorbar-canvas').onclick=()=>{
        paletteIdx=(paletteIdx+1)%PALETTE_NAMES.length;
        drawColorBar(); redrawGeo();
    };

    registerPlot('cl-energy-hist',  'cluster', 'Cluster Energy');
    registerPlot('cl-nclust-hist', 'cluster', 'Clusters per Event');
    registerPlot('cl-nblocks-hist','cluster', 'Blocks per Cluster');
    registerPlot('gem-ncl-hist',   'gem',     'GEM Clusters / Event');
    registerPlot('gem-theta-hist', 'gem',     'GEM Hit Angle');
    setupCopyBtn('btn-copy-cl-hist', ()=>currentClHist);
    setupCopyBtn('btn-copy-nclust', ()=>currentNclustHist);
    setupCopyBtn('btn-copy-nblocks', ()=>currentNblocksHist);
    setupCopyBtn('btn-copy-gem-ncl', ()=>currentGemNclHist);
    setupCopyBtn('btn-copy-gem-theta', ()=>currentGemThetaHist);

    // cluster stat row column divider
    setupDivider('div-cl-stat','x',
        ()=>document.querySelector('.cl-stat-cell'),
        ()=>document.querySelector('.cl-stat-row'),
        ()=>0, 80, 80, ()=>{
            try{Plotly.Plots.resize('cl-nclust-hist');}catch(e){}
            try{Plotly.Plots.resize('cl-nblocks-hist');}catch(e){}
        });

    // waveform stacking controls — reset checkbox to match JS state (browser may restore old form state)
    document.getElementById('wf-stack').checked=false;
    document.getElementById('wf-stack').onchange=e=>{
        wfStackEnabled=e.target.checked;
        document.getElementById('wf-stack-count').style.display=wfStackEnabled?'':'none';
        document.getElementById('btn-wf-stack-reset').style.display=wfStackEnabled?'':'none';
        if(!wfStackEnabled){ wfStackTraces=[]; wfStackModKey=''; }
        if(selectedModule) showWaveform(selectedModule);
    };
    document.getElementById('btn-wf-stack-reset').onclick=()=>{
        wfStackTraces=[]; wfStackModKey='';
        if(selectedModule) showWaveform(selectedModule);
    };

    // histogram log-scale toggles
    document.getElementById('heighthist-logy').onchange=()=>{ if(selectedModule) showHistograms(selectedModule); };
    document.getElementById('inthist-logy').onchange=()=>{ if(selectedModule) showHistograms(selectedModule); };
    document.getElementById('clhist-logy').onchange=plotClHist;
    setupCopyBtn('btn-copy-lms', ()=>currentLmsData);

    // cluster panel divider: histogram ↔ table
    setupDivider('div-cl-ht','y',
        ()=>document.getElementById('cl-hist-panel'),
        ()=>document.getElementById('cluster-panel'),
        ()=>0,
        80, 80, ()=>{try{Plotly.Plots.resize('cl-energy-hist');}catch(e){}});

    registerPlot('lms-plot', 'lms', 'LMS History');
    registerPlot('physics-plot',       'physics', null, PC_EPICS);
    registerPlot('moller-xy-plot',     'physics', null, PC_EPICS);
    registerPlot('moller-energy-plot', 'physics', null, PC_EPICS);
    for(let i=0;i<EPICS_NUM_SLOTS;i++)
        registerPlot('epics-plot-'+i, 'epics', null, PC_EPICS);

    initRegisteredPlots();

    setupDivider('div-lms-ht','y',
        ()=>document.getElementById('lms-plot-panel'),
        ()=>document.getElementById('lms-panel'),
        ()=>0,
        80, 80, ()=>{try{Plotly.Plots.resize('lms-plot');}catch(e){}});
    document.getElementById('lms-color-metric').onchange=geoLms;
    document.getElementById('lms-log-scale').onchange=geoLms;

    // LMS range editors
    function lmsRangeGet(isMax){
        const mt=document.getElementById('lms-color-metric').value;
        return getGeoRange('lms',mt)[isMax?1:0];
    }
    function lmsRangeSet(isMax, v){
        const mt=document.getElementById('lms-color-metric').value;
        const r=getGeoRange('lms',mt);
        if(isMax) setGeoRange('lms',mt,r[0],v);
        else setGeoRange('lms',mt,v,r[1]);
    }
    setupRangeEdit('lms-range-min-btn','lms-range-min-edit','lms-range-min-show',
        ()=>lmsRangeGet(false), v=>lmsRangeSet(false,v), geoLms);
    setupRangeEdit('lms-range-max-btn','lms-range-max-edit','lms-range-max-show',
        ()=>lmsRangeGet(true), v=>lmsRangeSet(true,v), geoLms);
    document.getElementById('lms-ref-select').onchange=e=>{
        g_lmsRefIndex=parseInt(e.target.value);
        fetchLmsSummary();
        if(lmsSelectedModule>=0){
            const name=lmsSummaryData&&lmsSummaryData.modules&&lmsSummaryData.modules[String(lmsSelectedModule)]
                ?lmsSummaryData.modules[String(lmsSelectedModule)].name:'';
            fetchLmsHistory(lmsSelectedModule, name);
        }
    };

    // --- file mode nav ---
    document.getElementById('btn-prev').onclick=()=>{if(currentEvent>1)loadEvent(currentEvent-1);};
    document.getElementById('btn-next').onclick=()=>{if(currentEvent<totalEvents)loadEvent(currentEvent+1);};
    document.getElementById('ev-input').onchange=e=>{const v=parseInt(e.target.value);if(v>=1&&v<=totalEvents)loadEvent(v);};
    document.getElementById('color-metric').onchange=()=>{syncDqRange();geoDq();};
    document.getElementById('log-scale').onchange=geoDq;
    document.getElementById('time-cut').onchange=geoDq;

    // --- file browser ---
    document.getElementById('btn-open').onclick = openFileDialog;
    document.getElementById('file-dialog-close').onclick = closeFileDialog;
    document.getElementById('file-backdrop').onclick = closeFileDialog;
    document.getElementById('file-filter').oninput = e => filterFileList(e.target.value);
    document.addEventListener('keydown', e => {
        if (e.key === 'Escape') {
            if (document.getElementById('file-dialog').classList.contains('open'))
                closeFileDialog();
            if (document.getElementById('et-dialog').classList.contains('open'))
                closeEtDialog();
        }
    });

    // --- range editing ---
    function setupRangeEdit(btnId, editId, showId, getVal, setVal, onApply) {
        const btn=document.getElementById(btnId);
        const edit=document.getElementById(editId);
        const show=document.getElementById(showId);
        let editing=false;

        function startEdit() {
            editing=true; btn.classList.add('editing'); btn.textContent='✓';
            edit.classList.add('active'); show.style.display='none';
            edit.value=getVal()||'';
            edit.focus(); edit.select();
        }
        function applyEdit() {
            if(!editing) return;
            editing=false; btn.classList.remove('editing'); btn.textContent='✎';
            edit.classList.remove('active'); show.style.display='';
            const v=parseFloat(edit.value);
            setVal(isNaN(v)?null:v);
            onApply();
        }

        btn.addEventListener('mousedown', e => e.preventDefault());
        btn.onclick=()=>{ if(editing) applyEdit(); else startEdit(); };
        edit.addEventListener('keydown',e=>{
            if(e.key==='Enter') applyEdit();
            if(e.key==='Escape'){editing=false;btn.classList.remove('editing');btn.textContent='✎';edit.classList.remove('active');show.style.display='';}
        });
        edit.addEventListener('blur',()=>{ applyEdit(); });
    }
    // DQ range editors
    function dqRangeApply(){
        const mt=document.getElementById('color-metric').value;
        setGeoRange('dq', mt, rangeMin, rangeMax);
        updateRangeDisplay(); geoDq();
    }
    setupRangeEdit('range-min-btn','range-min-edit','range-min-show',
        ()=>rangeMin, v=>{rangeMin=v;}, dqRangeApply);
    setupRangeEdit('range-max-btn','range-max-edit','range-max-show',
        ()=>rangeMax, v=>{rangeMax=v;}, dqRangeApply);
    // Cluster range editors
    function clRangeApply(){ geoCluster(); }
    setupRangeEdit('cl-range-min-btn','cl-range-min-edit','cl-range-min-show',
        ()=>getGeoRange('cluster','energy')[0],
        v=>{const r=getGeoRange('cluster','energy');setGeoRange('cluster','energy',v,r[1]);},
        clRangeApply);
    setupRangeEdit('cl-range-max-btn','cl-range-max-edit','cl-range-max-show',
        ()=>getGeoRange('cluster','energy')[1],
        v=>{const r=getGeoRange('cluster','energy');setGeoRange('cluster','energy',r[0],v);},
        clRangeApply);

    // --- online mode nav ---
    document.getElementById('ring-select').onchange=e=>{
        autoFollow=false; updateFollowStatus();
        loadEvent(parseInt(e.target.value));
    };
    document.getElementById('ring-select').onfocus=()=>{ updateRingSelector(); };
    document.getElementById('follow-status').onclick=()=>{ autoFollow=true; updateFollowStatus(); loadLatestEvent(); };
    // per-tab clear buttons

    // Clear All — resets all tabs' data for new run
    document.getElementById('btn-clear-all').onclick=()=>{
        // always clear server-side, then frontend
        Promise.all([
            fetch('/api/hist/clear').then(r=>r.json()),
            fetch('/api/lms/clear').then(r=>r.json()),
            fetch('/api/epics/clear').then(r=>r.json()),
        ]).then(clearFrontend).catch(()=>{
            document.getElementById('status-bar').textContent='Error clearing data';
        });
    };

    // mode toggle button — opens ET dialog when going online
    document.getElementById('btn-mode-toggle').onclick=()=>{
        if(mode==='online'){
            fetch('/api/mode/file',{method:'POST'}).then(()=>fetchConfigAndApply());
        } else {
            openEtDialog();
        }
    };

    // ET connect dialog
    const etBackdrop=document.getElementById('et-backdrop');
    const etDialog=document.getElementById('et-dialog');
    document.getElementById('et-dialog-close').onclick=()=>closeEtDialog();
    document.getElementById('et-cancel').onclick=()=>closeEtDialog();
    etBackdrop.onclick=()=>closeEtDialog();
    document.getElementById('et-connect').onclick=()=>{
        const cfg={
            host:    document.getElementById('et-input-host').value,
            port:    parseInt(document.getElementById('et-input-port').value)||11111,
            et_file: document.getElementById('et-input-file').value,
            station: document.getElementById('et-input-station').value,
        };
        document.getElementById('et-status-msg').textContent='Connecting...';
        fetch('/api/mode/online',{
            method:'POST',
            headers:{'Content-Type':'application/json'},
            body:JSON.stringify(cfg),
        }).then(r=>r.json()).then(d=>{
            if(d.error){
                document.getElementById('et-status-msg').textContent='Error: '+d.error;
            } else {
                closeEtDialog();
                fetchConfigAndApply();
            }
        }).catch(()=>{
            document.getElementById('et-status-msg').textContent='Connection failed';
        });
    };

    // geo mouse
    const tip=document.getElementById('geo-tooltip');
    // build tooltip text for a module
    function tooltipText(m){
        let t=`${m.n}  (${m.t==='G'?'PbGlass':'PbWO₄'})\n${crateName(m.roc)}  slot ${m.sl}  ch ${m.ch}`;
        if(activeTab==='lms' && lmsSummaryData){
            const idx=modules.indexOf(m);
            const md=lmsSummaryData.modules?lmsSummaryData.modules[String(idx)]:null;
            if(md){
                const rmsPct=md.mean>0?(md.rms/md.mean*100).toFixed(1):'--';
                t+=`\nLMS Mean: ${md.mean.toFixed(1)}  RMS: ${md.rms.toFixed(2)}  (${rmsPct}%)`;
                t+=`\n${md.count} pts  ${md.warn?'⚠ WARNING':'OK'}`;
            } else { t+='\nNo LMS data'; }
        } else if(activeTab==='cluster' && clusterData){
            const idx=modules.indexOf(m);
            const energy=clusterData.hits?clusterData.hits[String(idx)]:0;
            if(energy) t+=`\nEnergy: ${energy.toFixed(1)} MeV`;
            const clusters=clusterData.clusters||[];
            for(let ci=0;ci<clusters.length;ci++){
                if(clusters[ci].modules&&clusters[ci].modules.includes(idx)){
                    t+=`\nCluster #${ci} (${clusters[ci].center}, ${clusters[ci].energy.toFixed(0)} MeV)`;
                    break;
                }
            }
        } else {
            const d=eventChannels[`${m.roc}_${m.sl}_${m.ch}`];
            const key=`${m.roc}_${m.sl}_${m.ch}`;
            if(d&&d.pk&&d.pk.length){
                const pks=peaksInCut(d.pk);
                const bp=tallest(pks);
                const tc=isTimeCut();
                if(bp) t+=`\nPed ${d.pm.toFixed(1)}  H ${bp.h.toFixed(0)}  Int ${bp.i.toFixed(0)}  T ${bp.t.toFixed(0)}ns  Pk ${pks.length}${tc?' (tcut)':''}`;
                else t+=`\nPed ${d.pm.toFixed(1)}  (no peaks${tc?' in time cut':''})`;
            }
            else if(d)t+=`\nPed ${d.pm.toFixed(1)}  (no peaks)`;
            if(occTotal>0){
                const tc=isTimeCut();
                const occ=tc?occTcutData:occData;
                const pct=100.0*(occ[key]||0)/occTotal;
                t+=`\nOcc ${pct.toFixed(1)}%  (${occTotal} evts${tc?' tcut':''})`;
            } else if(histEnabled===false){
                t+=`\nOcc: not computed (enable histograms)`;
            }
        }
        return t;
    }
    // update tooltip content for currently hovered module (called on data refresh)
    updateGeoTooltip=()=>{
        if(hoveredModule) tip.textContent=tooltipText(hoveredModule);
    };
    geoCanvas.addEventListener('mousemove',e=>{
        const r=geoCanvas.getBoundingClientRect(),m=hitTest(e.clientX-r.left,e.clientY-r.top);
        if(m!==hoveredModule){
            hoveredModule=m;
            renderGeoOutlines(_geoOutlineFn, _geoDecorateFn);  // outlines only — fills unchanged
        }
        if(m){
            tip.textContent=tooltipText(m);tip.style.display='block';
            tip.style.left=(e.clientX-r.left+14)+'px';tip.style.top=(e.clientY-r.top-8)+'px';
        }else tip.style.display='none';
    });
    // click is now handled via geoHandleClick (called from mouseup when drag threshold not exceeded)
    geoCanvas.addEventListener('mouseleave',()=>{
        hoveredModule=null;tip.style.display='none';
        renderGeoOutlines(_geoOutlineFn, _geoDecorateFn);
    });
    document.addEventListener('keydown',e=>{
        if(e.target.tagName==='INPUT'||e.target.tagName==='SELECT')return;
        if(mode==='file'){
            if(e.key==='ArrowLeft'&&currentEvent>1)loadEvent(currentEvent-1);
            if(e.key==='ArrowRight'&&currentEvent<totalEvents)loadEvent(currentEvent+1);
        } else {
            if(e.key==='ArrowLeft'||e.key==='ArrowRight'){
                // navigate ring buffer
                const sel=document.getElementById('ring-select');
                const opts=[...sel.options].map(o=>parseInt(o.value));
                const cur=parseInt(sel.value);
                const idx=opts.indexOf(cur);
                if(e.key==='ArrowRight'&&idx>0){sel.value=opts[idx-1];autoFollow=false;updateFollowStatus();loadEvent(opts[idx-1]);}
                if(e.key==='ArrowLeft'&&idx<opts.length-1){sel.value=opts[idx+1];autoFollow=false;updateFollowStatus();loadEvent(opts[idx+1]);}
            }
            if(e.key==='f'||e.key==='F'){autoFollow=true;updateFollowStatus();loadLatestEvent();}
        }
    });

    // always connect WebSocket (for mode_changed and clear notifications)
    connectWebSocket();

    // load config and init mode
    fetchConfigAndApply();
}
window.addEventListener('DOMContentLoaded',init);

// Reset all frontend state (used by Clear All and mode switching)
function clearFrontend(){
    occData={}; occTcutData={}; occTotal=0;
    eventChannels={}; currentWaveform=null; currentHist={};

    // reset waveform stacking state
    wfStackTraces=[]; wfStackModKey=''; wfStackEnabled=false;
    wfRequestId++;  // invalidate any in-flight waveform fetches
    lastHistModule='';
    document.getElementById('wf-stack').checked=false;
    document.getElementById('wf-stack-count').style.display='none';
    document.getElementById('btn-wf-stack-reset').style.display='none';

    // deselect module and blank all DQ plots
    selectedModule=null;
    document.getElementById('detail-header').innerHTML=
        '<div class="empty-msg">Click a module to view details</div>';
    Plotly.react('waveform-div',[],{...PL,xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC'}},PC2);
    Plotly.react('heighthist-div',[],{...PL,title:{text:'Height Histogram',font:{size:10,color:'#555'}}},PC2);
    Plotly.react('inthist-div',[],{...PL,title:{text:'Integral Histogram',font:{size:10,color:'#555'}}},PC2);
    Plotly.react('poshist-div',[],{...PL,title:{text:'Peak Position',font:{size:10,color:'#555'}}},PC2);
    document.getElementById('peaks-tbody').innerHTML='';

    // cluster tab
    initClHist(); plotClHist(); plotClStatHists();
    clusterData=null; clusterEvent=-1; selectedCluster=-1;
    currentNclustHist=null; currentNblocksHist=null;
    document.getElementById('cl-select').innerHTML='<option value="all">All</option>';
    document.getElementById('cl-detail-header').innerHTML=
        '<span class="cl-info-text">Click a module or select a cluster</span>';
    document.getElementById('cl-tbody').innerHTML='';

    // LMS tab
    lmsSummaryData=null; lmsSelectedModule=-1; currentLmsData=null;
    Plotly.react('lms-plot',[],{...PL,title:{text:'LMS History',font:{size:10,color:'#555'}}},PC2);
    document.getElementById('lms-detail-header').innerHTML=
        '<span class="cl-info-text">Click a module to view LMS history</span>';
    document.getElementById('lms-tbody').innerHTML='';
    document.getElementById('lms-ref-select').innerHTML='<option value="-1">None</option>';

    document.getElementById('ring-select').innerHTML='';

    // GEM, EPICS, Physics tabs
    currentGemNclHist=null; currentGemThetaHist=null;
    clearEpicsFrontend();
    clearPhysicsFrontend();

    sampleCount=0;
    updateHeaderStats();
    redrawGeo();
    document.getElementById('status-bar').textContent='All data cleared';
}

// fetch /api/config and reconfigure the UI
function fetchConfigAndApply(){
    fetch('/api/config').then(r=>r.json()).then(applyConfig);
}

// ET connection dialog
function openEtDialog(){
    // populate fields with current ET config from last /api/config
    const etc=window._etConfig||{};
    document.getElementById('et-input-host').value=etc.host||'localhost';
    document.getElementById('et-input-port').value=etc.port||11111;
    document.getElementById('et-input-file').value=etc.et_file||'/tmp/et_sys_prad2';
    document.getElementById('et-input-station').value=etc.station||'prad2_monitor';
    document.getElementById('et-status-msg').textContent='';
    document.getElementById('et-backdrop').classList.add('open');
    document.getElementById('et-dialog').classList.add('open');
}
function closeEtDialog(){
    document.getElementById('et-backdrop').classList.remove('open');
    document.getElementById('et-dialog').classList.remove('open');
}

function applyConfig(data){
    const crateRoc=data.crate_roc||{};
    const rawMods=data.modules||[],rawDaq=data.daq||[];
    if(rawMods.length&&rawMods[0].roc!==undefined){
        modules=rawMods.map(m=>({...m,t:m.t==='PbGlass'?'G':(m.t==='G'?'G':'W')}));
    }else{
        const dm={};for(const d of rawDaq)dm[d.name]=d;modules=[];
        for(const m of rawMods){const d=dm[m.n];if(!d)continue;
            modules.push({n:m.n,t:m.t==='PbGlass'?'G':'W',x:m.x,y:m.y,sx:m.sx,sy:m.sy,
                roc:crateRoc[String(d.crate)]||0,sl:d.slot,ch:d.channel});}
    }
    totalEvents=data.total_events||0;
    histEnabled=data.hist_enabled||false;
    histConfig=data.hist||{};
    refLines=data.ref_lines||{};
    triggerBitsDef=data.trigger_bits||[];
    // cluster histogram configs
    if(data.cluster_hist){
        clHistMin=data.cluster_hist.min||0;
        clHistMax=data.cluster_hist.max||3000;
        clHistStep=data.cluster_hist.step||10;
    }
    if(data.nclusters_hist){
        nclustMin=data.nclusters_hist.min||0;
        nclustMax=data.nclusters_hist.max||20;
        nclustStep=data.nclusters_hist.step||1;
    }
    if(data.nblocks_hist){
        nblocksMin=data.nblocks_hist.min||0;
        nblocksMax=data.nblocks_hist.max||40;
        nblocksStep=data.nblocks_hist.step||1;
    }
    initClHist();
    if(data.lms){
        g_lmsWarnThresh=data.lms.warn_threshold||0.1;
        const sel=document.getElementById('lms-ref-select');
        sel.innerHTML='<option value="-1">None</option>';
        if(data.lms.ref_channels){
            for(const rc of data.lms.ref_channels){
                const o=document.createElement('option');
                o.value=rc.index;
                o.textContent=rc.name;
                sel.appendChild(o);
            }
        }
    }
    if(data.color_ranges){
        for(const [k,v] of Object.entries(data.color_ranges)){
            if(Array.isArray(v) && v.length===2 && !geoRangeOverrides[k])
                geoRangeOverrides[k]=v;
        }
    }
    if(data.refresh_ms){
        refreshEventMs=data.refresh_ms.event||200;
        refreshRingMs=data.refresh_ms.ring||500;
        refreshHistMs=data.refresh_ms.histogram||2000;
        refreshLmsMs=data.refresh_ms.lms||2000;
    }
    initReport(data);
    initEpics(data);
    initPhysics(data);
    updateTimeCutLabel();
    mode=data.mode||'file';
    etAvailable=data.et_available||false;
    if(data.et_config) window._etConfig=data.et_config;
    fileAvailable=data.file_available||false;
    if(data.source) sourceCaps=data.source;
    const appTitle=mode==='online'?'PRad-II HyCal Monitor':'PRad-II HyCal Event Viewer';
    document.title=appTitle;
    document.getElementById('app-title').textContent=appTitle;
    g_currentFile=data.current_file||'';
    g_dataDirEnabled=data.data_dir_enabled||false;
    g_dataDir=data.data_dir||'';
    g_histCheckbox=histEnabled;

    const hcb=document.getElementById('hist-checkbox');
    if(hcb) hcb.checked=histEnabled;

    // show/hide mode-specific UI
    document.getElementById('nav-file').style.display   = mode!=='online'?'flex':'none';
    document.getElementById('nav-online').style.display = mode==='online'?'flex':'none';
    document.getElementById('btn-open').style.display='';

    // mode toggle button — visible whenever ET is available
    const toggleBtn=document.getElementById('btn-mode-toggle');
    if(etAvailable){
        toggleBtn.style.display='';
        toggleBtn.textContent=mode==='online'?'Go Offline':'Go Online';
    } else {
        toggleBtn.style.display='none';
    }

    // hide tabs based on data source capabilities
    document.querySelectorAll('.tab').forEach(t=>{
        const tab=t.dataset.tab;
        if(tab==='lms')    t.style.display=sourceCaps.has_waveforms?'':'none';
        if(tab==='epics')  t.style.display=sourceCaps.has_epics?'':'none';
    });

    // auto-switch to cluster tab if source has no waveform data
    if(!sourceCaps.has_waveforms && activeTab==='dq'){
        switchTab('cluster');
    }

    if(mode==='file'){
        document.getElementById('ev-total').textContent=`/ ${totalEvents}`;

        updateHeaderStats();
        if(histEnabled) { fetchOccupancy(); fetchClHist(); }
        fetchEpicsChannels(); fetchEpicsLatest();
        if(activeTab==='epics') fetchAllEpicsSlots();
        if(activeTab==='physics') fetchPhysics();
        syncDqRange();
        geoViewInit=false; resizeGeo();
        if(totalEvents>0)loadEvent(1);
    } else if(mode==='online'){
        setEtStatus(data.et_connected||false, !data.et_connected);
        syncDqRange();
        fetchOccupancy();
        fetchEpicsChannels(); fetchEpicsLatest();
        if(activeTab==='physics') fetchPhysics();
        resizeGeo();
        updateRingSelector();
        loadLatestEvent();
    } else {
        // idle mode
        syncDqRange();
        resizeGeo();
        updateHeaderStats();
    }
}
