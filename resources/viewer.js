// =========================================================================
// State
// =========================================================================
let modules=[], totalEvents=0, currentEvent=1;
let eventChannels={};
let selectedModule=null, hoveredModule=null;
let geoCanvas, geoCtx, geoWrap, scale=1, offsetX=0, offsetY=0, canvasW, canvasH;
const PC=['#00b4d8','#ff6b6b','#51cf66','#ffd43b','#cc5de8','#ff922b','#20c997','#f06595'];
const CRATE_NAME={0x80:'adchycal1',0x82:'adchycal2',0x84:'adchycal3',0x86:'adchycal4',0x88:'adchycal5',0x8a:'adchycal6',0x8c:'adchycal7'};
function crateName(r){return CRATE_NAME[r]||`ROC 0x${r.toString(16)}`;}
let histEnabled=false, histConfig={};
let mode='file';    // 'file' or 'online'
let ws=null;        // WebSocket connection (online mode)
let autoFollow=true; // auto-load latest event
let lastEventFetch=0, lastHistFetch=0;  // throttle timestamps

// occupancy data (fetched once per file load when histograms enabled)
let occData={}, occTcutData={}, occTotal=0;

// color range: null = auto from data, number = user-set or hist-synced
let rangeMin=null, rangeMax=null;

// =========================================================================
// Color scale
// =========================================================================
function viridis(t){
    t=Math.max(0,Math.min(1,t));
    return `rgb(${Math.round(255*Math.min(1,Math.max(0,-0.87+4.26*t-4.85*t*t+2.5*t*t*t)))},${Math.round(255*Math.min(1,Math.max(0,-0.03+0.77*t+1.32*t*t-1.87*t*t*t)))},${Math.round(255*Math.min(1,Math.max(0,0.33+1.74*t-4.26*t*t+3.17*t*t*t)))})`;
}
function drawColorBar(){const c=document.getElementById('colorbar-canvas'),x=c.getContext('2d');for(let i=0;i<c.width;i++){x.fillStyle=viridis(i/c.width);x.fillRect(i,0,1,c.height);}}

// sync color range from histogram config when a matching metric is selected
function syncRangeFromHist(){
    const mt=document.getElementById('color-metric').value;
    const h=histConfig;
    if(mt==='integral' && h.bin_min!==undefined){
        rangeMin=h.bin_min; rangeMax=h.bin_max;
    } else if(mt==='time' && h.pos_min!==undefined){
        rangeMin=h.pos_min; rangeMax=h.pos_max;
    } else if((mt==='occupancy' || mt==='occupancy_tcut') && occTotal>0){
        rangeMin=0; rangeMax=occTotal;
    } else {
        rangeMin=null; rangeMax=null;
    }
    updateRangeDisplay();
}

function updateRangeDisplay(){
    const minEl=document.getElementById('range-min-show');
    const maxEl=document.getElementById('range-max-show');
    minEl.textContent=(rangeMin!==null)?rangeMin.toFixed(0):'auto';
    maxEl.textContent=(rangeMax!==null)?rangeMax.toFixed(0):'auto';
}

// =========================================================================
// Geo canvas
// =========================================================================
function initGeo(){
    geoWrap=document.getElementById('geo-wrap');
    geoCanvas=document.getElementById('geo-canvas');
    geoCtx=geoCanvas.getContext('2d');
    resizeGeo();
    new ResizeObserver(resizeGeo).observe(geoWrap);
}
function resizeGeo(){
    canvasW=geoWrap.clientWidth; canvasH=geoWrap.clientHeight;
    if(canvasW<10||canvasH<10)return;
    geoCanvas.width=canvasW; geoCanvas.height=canvasH;
    if(modules.length)fitView(); drawGeo();
}
function fitView(){
    const m=15;let x0=1e9,x1=-1e9,y0=1e9,y1=-1e9;
    for(const d of modules){x0=Math.min(x0,d.x-d.sx/2);x1=Math.max(x1,d.x+d.sx/2);y0=Math.min(y0,d.y-d.sy/2);y1=Math.max(y1,d.y+d.sy/2);}
    scale=Math.min((canvasW-2*m)/(x1-x0),(canvasH-2*m)/(y1-y0));
    offsetX=canvasW/2-(x0+x1)/2*scale; offsetY=canvasH/2+(y0+y1)/2*scale;
}
function d2c(x,y){return[x*scale+offsetX,-y*scale+offsetY];}
function c2d(cx,cy){return[(cx-offsetX)/scale,-(cy-offsetY)/scale];}
function modVal(m){
    const key=`${m.roc}_${m.sl}_${m.ch}`;
    const mt=document.getElementById('color-metric').value;
    // occupancy metrics use pre-computed counts, not per-event data
    if(mt==='occupancy') return occTotal>0 ? (occData[key]||0) : null;
    if(mt==='occupancy_tcut') return occTotal>0 ? (occTcutData[key]||0) : null;
    const d=eventChannels[key];
    if(!d)return null;
    if(mt==='pedestal')return d.pm||0;
    if(!d.pk||!d.pk.length)return null;
    if(mt==='height')return Math.max(...d.pk.map(p=>p.h));
    if(mt==='time'){
        let best=d.pk[0];
        for(let i=1;i<d.pk.length;i++) if(d.pk[i].h>best.h) best=d.pk[i];
        return best.t;
    }
    return d.pk.reduce((s,p)=>s+p.i,0); // integral
}
function drawGeo(){
    if(!geoCtx)return;const ctx=geoCtx;ctx.clearRect(0,0,canvasW,canvasH);
    const useLog=document.getElementById('log-scale').checked;
    const vals=modules.map(modVal);
    const numVals=vals.filter(v=>v!==null);
    const autoMin=numVals.length?Math.min(...numVals):0;
    const autoMax=numVals.length?Math.max(...numVals):1;
    const vmin=rangeMin!==null?rangeMin:0;
    const vmax=rangeMax!==null?rangeMax:Math.max(autoMax,1);
    const span=vmax-vmin||1;

    // update display if auto
    if(rangeMin===null) document.getElementById('range-min-show').textContent=vmin.toFixed(0);
    if(rangeMax===null) document.getElementById('range-max-show').textContent=vmax.toFixed(0);

    for(let i=0;i<modules.length;i++){
        const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale,v=vals[i];
        let t=0;
        if(v!==null){
            const clamped=Math.max(vmin,Math.min(vmax,v));
            if(useLog) t=Math.log1p(clamped-vmin)/Math.log1p(span);
            else t=(clamped-vmin)/span;
        }
        ctx.fillStyle=(v!==null)?viridis(t):(m.t==='G'?'#1a1a2e':'#12122a');
        ctx.fillRect(cx-w/2,cy-h/2,w,h);
        const sel=selectedModule&&selectedModule.n===m.n,hov=hoveredModule&&hoveredModule.n===m.n;
        ctx.strokeStyle=sel?'#fff':hov?'#00b4d8':'#333';ctx.lineWidth=sel?2.5:hov?1.5:0.5;
        ctx.strokeRect(cx-w/2,cy-h/2,w,h);
    }
}
function hitTest(cx,cy){
    const[dx,dy]=c2d(cx,cy);
    for(let i=modules.length-1;i>=0;i--){const m=modules[i];if(Math.abs(dx-m.x)<=m.sx/2&&Math.abs(dy-m.y)<=m.sy/2)return m;}
    return null;
}

// =========================================================================
// Plotly shared config
// =========================================================================
const PL={paper_bgcolor:'#1a1a2e',plot_bgcolor:'#11112a',font:{family:'Consolas,monospace',size:10,color:'#aaa'},
    margin:{l:45,r:10,t:24,b:32},xaxis:{gridcolor:'#1a1a3a',zerolinecolor:'#222'},
    yaxis:{gridcolor:'#1a1a3a',zerolinecolor:'#222'}};
const PC2={responsive:true,displayModeBar:false};

function resizeAllPlots(){
    for(const id of['waveform-div','inthist-div','poshist-div'])
        try{Plotly.Plots.resize(id);}catch(e){}
}

// =========================================================================
// Waveform
// =========================================================================
function showWaveform(mod){
    selectedModule=mod;
    const key=`${mod.roc}_${mod.sl}_${mod.ch}`;
    const d=eventChannels[key];
    const pedInfo=d?` &nbsp; Ped: ${d.pm.toFixed(1)} ± ${d.pr.toFixed(1)}`:'';
    document.getElementById('detail-header').innerHTML=
        `<span class="mod-name">${mod.n}</span> <span class="mod-daq">${crateName(mod.roc)} &middot; slot ${mod.sl} &middot; ch ${mod.ch}${pedInfo}</span>`;

    if(!d||!d.s){
        Plotly.react('waveform-div',[],{...PL,title:{text:`${mod.n} — No data`,font:{size:11,color:'#555'}}},PC2);
        document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No data</td></tr>';
        showHistograms(mod); drawGeo(); return;
    }

    const samples=d.s, peaks=d.pk||[], x=samples.map((_,i)=>i);
    const traces=[
        {x,y:samples,type:'scatter',mode:'lines',name:'Waveform',line:{color:'#7777aa',width:1}},
        {x:[0,samples.length-1],y:[d.pm,d.pm],type:'scatter',mode:'lines',name:'Pedestal',line:{color:'#555',width:1,dash:'dash'}},
    ];
    const thr=d.pm+Math.max(5*d.pr,3);
    traces.push({x:[0,samples.length-1],y:[thr,thr],type:'scatter',mode:'lines',line:{color:'#333',width:1,dash:'dot'},showlegend:false});
    peaks.forEach((p,i)=>{
        const col=PC[i%PC.length],px=[],py=[];
        for(let j=p.l;j<=p.r;j++){px.push(j);py.push(samples[j]);}
        // hex to rgba for fill
        const r=parseInt(col.slice(1,3),16),g=parseInt(col.slice(3,5),16),b=parseInt(col.slice(5,7),16);
        const fill=`rgba(${r},${g},${b},0.18)`;
        // shaded integral region: pedestal baseline then waveform with fill between
        traces.push({x:px,y:px.map(()=>d.pm),type:'scatter',mode:'lines',
            line:{width:0},showlegend:false,hoverinfo:'skip'});
        traces.push({x:px,y:py,type:'scatter',mode:'lines',name:`Peak ${i}`,
            line:{color:col,width:2},fill:'tonexty',fillcolor:fill});
        traces.push({x:[p.p],y:[samples[p.p]],type:'scatter',mode:'markers',
            marker:{color:col,size:7,symbol:'diamond'},showlegend:false});
    });

    Plotly.react('waveform-div',traces,{...PL,
        title:{text:`${mod.n} — Event ${currentEvent}`,font:{size:11,color:'#ccc'}},
        xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC'},
        legend:{x:1,y:1,xanchor:'right',bgcolor:'rgba(0,0,0,0.6)',font:{size:9}},
    },PC2);

    // peaks table
    let rows='';
    peaks.forEach((p,i)=>{
        const col=PC[i%PC.length];
        rows+=`<tr style="border-left:3px solid ${col}"><td>${i}</td><td>${p.p}</td><td>${p.t.toFixed(0)}</td><td>${p.h.toFixed(1)}</td><td>${p.i.toFixed(0)}</td><td>${p.l}</td><td>${p.r}</td><td style="text-align:center">${p.o?'⚠':''}</td></tr>`;
    });
    if(!peaks.length) rows='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No peaks</td></tr>';
    document.getElementById('peaks-tbody').innerHTML=rows;

    showHistograms(mod);
    drawGeo();
}

// =========================================================================
// Histograms
// =========================================================================
function fetchAndPlotHist(divId, url, title, xTitle, binMin, binStep, barColor){
    if(!histEnabled){
        Plotly.react(divId,[],{...PL,title:{text:'--hist not enabled',font:{size:10,color:'#444'}}},PC2);
        return;
    }
    fetch(url).then(r=>r.json()).then(data=>{
        if(data.error||!data.bins||!data.bins.length){
            Plotly.react(divId,[],{...PL,title:{text:`${title} — No data`,font:{size:10,color:'#555'}}},PC2);
            return;
        }
        const x=data.bins.map((_,i)=>binMin+(i+0.5)*binStep);
        const entries=data.bins.reduce((a,b)=>a+b,0)+data.underflow+data.overflow;
        const stats=`${data.events} evts | Entries: ${entries}  Under: ${data.underflow}  Over: ${data.overflow}`;
        Plotly.react(divId,[{
            x,y:data.bins,type:'bar',marker:{color:barColor,line:{width:0}},
            hovertemplate:'%{x:.0f}: %{y}<extra></extra>',
        }],{...PL,
            title:{text:`${title}<br><span style="font-size:9px;color:#888">${stats}</span>`,font:{size:10,color:'#ccc'}},
            xaxis:{...PL.xaxis,title:xTitle},yaxis:{...PL.yaxis,title:'Counts'},bargap:0.05,
        },PC2);
    }).catch(()=>{
        Plotly.react(divId,[],{...PL,title:{text:'Fetch error',font:{size:10,color:'#f66'}}},PC2);
    });
}

function showHistograms(mod){
    // in online mode, throttle histogram fetches to ~1 Hz
    if (mode === 'online') {
        const now = Date.now();
        if (now - lastHistFetch < 1000) return;
        lastHistFetch = now;
    }
    const key=`${mod.roc}_${mod.sl}_${mod.ch}`;
    const h=histConfig;
    fetchAndPlotHist('inthist-div',`/api/hist/${key}`,
        `${mod.n} Integral [${h.time_min||170}-${h.time_max||190} ns]`,
        'Peak Integral', h.bin_min||0, h.bin_step||100, '#00b4d8');
    fetchAndPlotHist('poshist-div',`/api/poshist/${key}`,
        `${mod.n} Peak Position`,
        'Time (ns)', h.pos_min||0, h.pos_step||4, '#51cf66');
}

// =========================================================================
// Event loading (works for both file and online mode)
// =========================================================================
function loadEventData(data) {
    if (data.error) { document.getElementById('status-bar').textContent = data.error; return; }
    currentEvent = data.event;
    eventChannels = data.channels || {};
    const nch = Object.keys(eventChannels).length;
    const npk = Object.values(eventChannels).reduce((s,c) => s + (c.pk||[]).length, 0);
    const modeTag = mode === 'online' ? ' [LIVE]' : '';
    document.getElementById('status-bar').textContent = `Event ${currentEvent}: ${nch} channels, ${npk} peaks${modeTag}`;
    drawGeo();
    if (selectedModule) showWaveform(selectedModule);
}

function loadEvent(evnum) {
    currentEvent = evnum;
    if (mode === 'file') document.getElementById('ev-input').value = evnum;
    document.getElementById('status-bar').textContent = `Loading event ${evnum}...`;
    fetch(`/api/event/${evnum}`).then(r => r.json()).then(loadEventData)
        .catch(err => { document.getElementById('status-bar').textContent = `Error: ${err}`; });
}

function loadLatestEvent() {
    fetch('/api/event/latest').then(r => r.json()).then(loadEventData)
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
            o.textContent = `Event ${ring[i]}` + (i === ring.length - 1 ? ' (latest)' : '');
            sel.appendChild(o);
        }
        // keep selection if not auto-following
        if (!autoFollow && prev && ring.includes(parseInt(prev))) sel.value = prev;
        else if (ring.length) sel.value = ring[ring.length - 1];
    });
}

function setEtStatus(connected) {
    const el = document.getElementById('et-status');
    el.textContent = connected ? '● Connected' : '● Disconnected';
    el.style.color = connected ? '#51cf66' : '#f66';
}

function connectWebSocket() {
    const proto = location.protocol === 'https:' ? 'wss:' : 'ws:';
    ws = new WebSocket(`${proto}//${location.host}`);

    ws.onopen = () => { console.log('WS connected'); };
    ws.onclose = () => {
        console.log('WS disconnected, reconnecting in 2s');
        setTimeout(connectWebSocket, 2000);
    };
    ws.onmessage = (evt) => {
        try {
            const msg = JSON.parse(evt.data);
            if (msg.type === 'new_event') {
                const now = Date.now();
                // throttle event display to ~5 Hz
                if (autoFollow && now - lastEventFetch > 200) {
                    lastEventFetch = now;
                    loadLatestEvent();
                }
                // throttle ring selector update to ~2 Hz
                if (now - lastHistFetch > 500) {
                    lastHistFetch = now;
                    updateRingSelector();
                }
            } else if (msg.type === 'status') {
                setEtStatus(msg.connected);
            } else if (msg.type === 'hist_cleared') {
                if (selectedModule) showHistograms(selectedModule);
            }
        } catch (e) {}
    };
}

function clearHistograms() {
    fetch('/api/hist/clear').then(r => r.json()).then(data => {
        document.getElementById('status-bar').textContent = 'Histograms cleared';
        if (selectedModule) showHistograms(selectedModule);
    });
}

// =========================================================================
// File browser
// =========================================================================
let allFiles = [];

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

    fetch('/api/files').then(r => r.json()).then(data => {
        allFiles = data.files || [];
        renderFileList('');
        filter.focus();
    });
}

function closeFileDialog() {
    document.getElementById('file-dialog').classList.remove('open');
    document.getElementById('file-backdrop').classList.remove('open');
}

function renderFileList(filter) {
    const list = document.getElementById('file-list');
    const currentFile = (g_currentFile || '').replace(/.*\//, ''); // basename
    const filt = filter.toLowerCase();
    let html = '';
    for (const f of allFiles) {
        if (filt && !f.path.toLowerCase().includes(filt)) continue;
        const isCurrent = f.path === g_currentFile || f.path.endsWith('/' + currentFile);
        html += `<div class="file-item${isCurrent ? ' current' : ''}" data-path="${f.path}">
            <span>${f.path}</span><span class="fsize">${f.size_mb} MB</span></div>`;
    }
    if (!html) html = '<div style="padding:12px;color:var(--dim);text-align:center">No matching files</div>';
    list.innerHTML = html;

    // click handlers
    list.querySelectorAll('.file-item[data-path]').forEach(el => {
        el.onclick = () => {
            closeFileDialog();
            loadNewFile(el.dataset.path);
        };
    });
}

let g_currentFile = '';
let g_histCheckbox = false;
let g_dataDirEnabled = false;
let g_dataDir = '';  // tracks the "Process histograms" checkbox

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
            fetch('/api/config').then(r => r.json()).then(cfg => {
                totalEvents = cfg.total_events || 0;
                g_currentFile = cfg.current_file || '';
                histEnabled = cfg.hist_enabled || false;
                histConfig = cfg.hist || {};
                g_histCheckbox = histEnabled;
                const hcb = document.getElementById('hist-checkbox');
                if (hcb) hcb.checked = histEnabled;
                document.getElementById('ev-total').textContent = `/ ${totalEvents}`;
                updateHeaderInfo(cfg);
                if (histEnabled) fetchOccupancy();
                syncRangeFromHist();
                drawGeo();
                if (totalEvents > 0) loadEvent(1);
            });
            return;
        }
        const pct = data.total > 0 ? Math.round(100 * data.current / data.total) : 0;
        const phaseText = data.phase === 'indexing' ? 'Indexing events' : 'Building histograms';
        document.getElementById('progress-bar').style.width = `${Math.min(pct, 100)}%`;
        document.getElementById('progress-text').textContent =
            `${phaseText}... ${data.current}` + (data.total > 0 ? ` / ${data.total}` : '');
        setTimeout(pollProgress, 300);
    }).catch(() => setTimeout(pollProgress, 1000));
}

function fetchOccupancy() {
    fetch('/api/occupancy').then(r => r.json()).then(data => {
        occData = data.occ || {};
        occTcutData = data.occ_tcut || {};
        occTotal = data.total || 0;
        // redraw if currently showing occupancy
        const mt = document.getElementById('color-metric').value;
        if (mt === 'occupancy' || mt === 'occupancy_tcut') {
            syncRangeFromHist();
            drawGeo();
        }
    }).catch(() => {});
}

function updateHeaderInfo(cfg) {
    const hi = cfg.hist_enabled ? ` · hist [${(cfg.hist||{}).time_min||170}-${(cfg.hist||{}).time_max||190} ns]` : '';
    const fname = g_currentFile.replace(/.*\//, '') || '';
    document.getElementById('header-info').textContent =
        `${modules.length} mod · ${totalEvents} evts${hi}` + (fname ? ` · ${fname}` : '');
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
        if(axis==='x') target.style.width=val+'px'; else target.style.height=val+'px';
        onResize();
    });
    document.addEventListener('mouseup',()=>{
        if(!active)return; active=false; div.classList.remove('active');
        document.body.style.cursor=''; document.body.style.userSelect='';
    });
}

// =========================================================================
// Init
// =========================================================================
function init(){
    drawColorBar(); initGeo();
    Plotly.newPlot('waveform-div',[],{...PL,xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC'}},PC2);
    Plotly.newPlot('inthist-div',[],{...PL,title:{text:'Integral Histogram',font:{size:10,color:'#555'}}},PC2);
    Plotly.newPlot('poshist-div',[],{...PL,title:{text:'Position Histogram',font:{size:10,color:'#555'}}},PC2);

    // --- dividers ---
    // 1. main vertical: geo ↔ detail
    setupDivider('div-main','x',
        ()=>document.getElementById('geo-panel'),
        ()=>document.querySelector('.main'),
        ()=>0, 300, 350, ()=>{resizeGeo();resizeAllPlots();});
    // 2. horizontal: top-panel (histograms) ↔ bottom-panel
    setupDivider('div-tb','y',
        ()=>document.getElementById('top-panel'),
        ()=>document.getElementById('detail-panel'),
        ()=>document.getElementById('detail-header').offsetHeight,
        120, 120, resizeAllPlots);
    // 3. horizontal inside top-panel: inthist ↔ poshist
    setupDivider('div-hist','y',
        ()=>document.getElementById('inthist-div'),
        ()=>document.getElementById('top-panel'),
        ()=>0, 60, 60, resizeAllPlots);
    // 4. vertical inside bottom-panel: table ↔ waveform
    setupDivider('div-tw','x',
        ()=>document.getElementById('table-wrap'),
        ()=>document.getElementById('bottom-panel'),
        ()=>0, 150, 200, resizeAllPlots);

    // --- file mode nav ---
    document.getElementById('btn-prev').onclick=()=>{if(currentEvent>1)loadEvent(currentEvent-1);};
    document.getElementById('btn-next').onclick=()=>{if(currentEvent<totalEvents)loadEvent(currentEvent+1);};
    document.getElementById('ev-input').onchange=e=>{const v=parseInt(e.target.value);if(v>=1&&v<=totalEvents)loadEvent(v);};
    document.getElementById('color-metric').onchange=()=>{syncRangeFromHist();drawGeo();};
    document.getElementById('log-scale').onchange=drawGeo;

    // --- file browser ---
    document.getElementById('btn-open').onclick = openFileDialog;
    document.getElementById('file-dialog-close').onclick = closeFileDialog;
    document.getElementById('file-backdrop').onclick = closeFileDialog;
    document.getElementById('file-filter').oninput = e => renderFileList(e.target.value);
    document.addEventListener('keydown', e => {
        if (e.key === 'Escape' && document.getElementById('file-dialog').classList.contains('open'))
            closeFileDialog();
    });

    // --- range editing ---
    function setupRangeEdit(btnId, editId, showId, isMax) {
        const btn=document.getElementById(btnId);
        const edit=document.getElementById(editId);
        const show=document.getElementById(showId);
        let editing=false;

        function startEdit() {
            editing=true; btn.classList.add('editing'); btn.textContent='✓';
            edit.classList.add('active'); show.style.display='none';
            edit.value=(isMax?rangeMax:rangeMin)||'';
            edit.focus(); edit.select();
        }
        function applyEdit() {
            editing=false; btn.classList.remove('editing'); btn.textContent='✎';
            edit.classList.remove('active'); show.style.display='';
            const v=parseFloat(edit.value);
            if(isMax) rangeMax=isNaN(v)?null:v; else rangeMin=isNaN(v)?null:v;
            updateRangeDisplay(); drawGeo();
        }

        btn.onclick=()=>{ if(editing) applyEdit(); else startEdit(); };
        edit.addEventListener('keydown',e=>{ if(e.key==='Enter') applyEdit(); if(e.key==='Escape'){editing=false;btn.classList.remove('editing');btn.textContent='✎';edit.classList.remove('active');show.style.display='';} });
        edit.addEventListener('blur',()=>{ if(editing) applyEdit(); });
    }
    setupRangeEdit('range-min-btn','range-min-edit','range-min-show',false);
    setupRangeEdit('range-max-btn','range-max-edit','range-max-show',true);

    // --- online mode nav ---
    document.getElementById('ring-select').onchange=e=>{
        autoFollow=false;
        loadEvent(parseInt(e.target.value));
    };
    document.getElementById('btn-clear-hist').onclick=clearHistograms;

    // geo mouse
    const tip=document.getElementById('geo-tooltip');
    geoCanvas.addEventListener('mousemove',e=>{
        const r=geoCanvas.getBoundingClientRect(),m=hitTest(e.clientX-r.left,e.clientY-r.top);
        if(m!==hoveredModule){hoveredModule=m;drawGeo();}
        if(m){
            const d=eventChannels[`${m.roc}_${m.sl}_${m.ch}`];
            const key=`${m.roc}_${m.sl}_${m.ch}`;
            let t=`${m.n}  (${m.t==='G'?'PbGlass':'PbWO₄'})\n${crateName(m.roc)}  slot ${m.sl}  ch ${m.ch}`;
            if(d&&d.pk&&d.pk.length){
                let best=d.pk[0]; for(let j=1;j<d.pk.length;j++) if(d.pk[j].h>best.h) best=d.pk[j];
                t+=`\nPed ${d.pm.toFixed(1)}  H ${best.h.toFixed(0)}  Int ${best.i.toFixed(0)}  T ${best.t.toFixed(0)}ns`;
            }
            else if(d)t+=`\nPed ${d.pm.toFixed(1)}  (no peaks)`;
            if(occTotal>0){
                const o=occData[key]||0, ot=occTcutData[key]||0;
                t+=`\nOcc ${o}/${occTotal}  Occ(t) ${ot}/${occTotal}`;
            } else if(histEnabled===false){
                t+=`\nOcc: not computed (enable histograms)`;
            }
            tip.textContent=t;tip.style.display='block';
            tip.style.left=(e.clientX-r.left+14)+'px';tip.style.top=(e.clientY-r.top-8)+'px';
        }else tip.style.display='none';
    });
    geoCanvas.addEventListener('click',e=>{const r=geoCanvas.getBoundingClientRect(),m=hitTest(e.clientX-r.left,e.clientY-r.top);if(m)showWaveform(m);});
    geoCanvas.addEventListener('mouseleave',()=>{hoveredModule=null;tip.style.display='none';drawGeo();});
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
                if(e.key==='ArrowRight'&&idx>0){sel.value=opts[idx-1];autoFollow=false;loadEvent(opts[idx-1]);}
                if(e.key==='ArrowLeft'&&idx<opts.length-1){sel.value=opts[idx+1];autoFollow=false;loadEvent(opts[idx+1]);}
            }
            if(e.key==='f'||e.key==='F'){autoFollow=true;loadLatestEvent();}  // F = follow latest
        }
    });

    // load config and init mode
    fetch('/api/config').then(r=>r.json()).then(data=>{
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
        mode=data.mode||'file';
        g_currentFile=data.current_file||'';
        g_dataDirEnabled=data.data_dir_enabled||false;
        g_dataDir=data.data_dir||'';
        g_histCheckbox=histEnabled;

        // init histogram checkbox
        const hcb=document.getElementById('hist-checkbox');
        if(hcb) hcb.checked=histEnabled;

        // show/hide mode-specific UI
        document.getElementById('nav-file').style.display   = mode==='file'?'flex':'none';
        document.getElementById('nav-online').style.display = mode==='online'?'flex':'none';

        // always show file browser button
        document.getElementById('btn-open').style.display='';

        if(mode==='file'){
            document.getElementById('ev-total').textContent=`/ ${totalEvents}`;
            updateHeaderInfo(data);
            if(histEnabled) fetchOccupancy();
            syncRangeFromHist();
            resizeGeo();
            if(totalEvents>0)loadEvent(1);
        } else {
            setEtStatus(data.et_connected||false);
            document.getElementById('header-info').textContent=
                `${modules.length} modules · ONLINE · ring ${data.ring_buffer_size||20}`;
            resizeGeo();
            connectWebSocket();
            updateRingSelector();
            loadLatestEvent();
        }
    });
}
window.addEventListener('DOMContentLoaded',init);
