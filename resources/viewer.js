// =========================================================================
// State
// =========================================================================
let modules=[], totalEvents=0, currentEvent=1;
let eventChannels={};
let selectedModule=null, hoveredModule=null;
let geoCanvas, geoCtx, geoWrap, scale=1, offsetX=0, offsetY=0, canvasW, canvasH;
const PC=['#00b4d8','#ff6b6b','#51cf66','#ffd43b','#cc5de8','#ff922b','#20c997','#f06595'];
const CRATE_NAME={0x80:'adchycal1',0x82:'adchycal2',0x84:'adchycal3',0x86:'adchycal4',0x88:'adchycal5',0x8a:'adchycal6',0x8c:'adchycal7',
    0x01:'PRadTS',0x04:'PRadROC_1',0x05:'PRadROC_2',0x06:'PRadROC_3',0x07:'PRadSRS_1',0x08:'PRadSRS_2'};
function crateName(r){return CRATE_NAME[r]||`ROC 0x${r.toString(16)}`;}
let histEnabled=false, histConfig={};
let mode='file';    // 'file' or 'online'
let ws=null;        // WebSocket connection (online mode)
let autoFollow=true; // auto-load latest event
let lastEventFetch=0, lastHistFetch=0, lastRingFetch=0, lastOccFetch=0;  // throttle timestamps

// occupancy data (fetched once per file load when histograms enabled)
let occData={}, occTcutData={}, occTotal=0;
let currentWaveform=null;  // {x:[], y:[]} for copy button
let currentHist={};  // {divId: {x:[], y:[]}} for histogram copy

// color range: null = auto from data, number = user-set or hist-synced
let rangeMin=null, rangeMax=null;
let rangeUserEdited=false;
const rangeUserOverrides={};  // {metric: [min, max]} — persists user edits per metric

// --- clustering tab state ---
let activeTab='dq';  // 'dq' or 'cluster'
let clusterData=null;  // {hits:{}, clusters:[]}
let selectedCluster=-1;  // -1 = all
let clusterEvent=-1;  // event number for cached cluster data
let clRangeMin=null, clRangeMax=null;  // null = auto from data

// cluster energy histogram (accumulated on frontend)
let clHistBins=null, clHistEvents=0;
let clHistMin=0, clHistMax=3000, clHistStep=10;
let currentClHist=null;  // {x:[], y:[]} for copy button

// default ranges per metric (overridden by hist config)
const RANGE_DEFAULTS={
    integral:  [0, 10000],
    height:    [0, 1000],
    count:     [0, 10],
    time:      [0, 400],
    pedestal:  [0, 500],
    occupancy: [0, 100],
};

// =========================================================================
// Color scale — click the colorbar to cycle palettes
// =========================================================================
const PALETTES={
    viridis(t){
        return `rgb(${Math.round(255*Math.min(1,Math.max(0,-0.87+4.26*t-4.85*t*t+2.5*t*t*t)))},${Math.round(255*Math.min(1,Math.max(0,-0.03+0.77*t+1.32*t*t-1.87*t*t*t)))},${Math.round(255*Math.min(1,Math.max(0,0.33+1.74*t-4.26*t*t+3.17*t*t*t)))})`;
    },
    inferno(t){
        return `rgb(${Math.round(255*Math.min(1,Math.max(0,-0.02+2.16*t+4.79*t*t-8.13*t*t*t+2.17*t*t*t*t)))},${Math.round(255*Math.min(1,Math.max(0,-0.02-0.35*t+5.87*t*t-8.29*t*t*t+3.7*t*t*t*t)))},${Math.round(255*Math.min(1,Math.max(0,0.01+3.1*t-9.34*t*t+12.45*t*t*t-5.24*t*t*t*t)))})`;
    },
    coolwarm(t){
        const r=Math.round(255*Math.min(1,Math.max(0, 0.23+2.22*t-1.83*t*t)));
        const g=Math.round(255*Math.min(1,Math.max(0, 0.30+1.58*t-2.36*t*t+0.56*t*t*t)));
        const b=Math.round(255*Math.min(1,Math.max(0, 0.75-0.44*t-0.81*t*t+0.53*t*t*t)));
        return `rgb(${r},${g},${b})`;
    },
    hot(t){
        const r=Math.round(255*Math.min(1, t*2.8));
        const g=Math.round(255*Math.max(0, Math.min(1, (t-0.35)*2.8)));
        const b=Math.round(255*Math.max(0, Math.min(1, (t-0.7)*3.3)));
        return `rgb(${r},${g},${b})`;
    },
    jet(t){
        t=0.125+t*0.75;  // remap to skip dark blue/red ends
        const r=Math.round(255*Math.min(1, Math.max(0, 1.5-Math.abs(t-0.75)*4)));
        const g=Math.round(255*Math.min(1, Math.max(0, 1.5-Math.abs(t-0.5)*4)));
        const b=Math.round(255*Math.min(1, Math.max(0, 1.5-Math.abs(t-0.25)*4)));
        return `rgb(${r},${g},${b})`;
    },
    greyscale(t){
        const v=Math.round(255*t);
        return `rgb(${v},${v},${v})`;
    },
};
const PALETTE_NAMES=Object.keys(PALETTES);
let paletteIdx=0;
function colorScale(t){ t=Math.max(0,Math.min(1,t)); return PALETTES[PALETTE_NAMES[paletteIdx]](t); }
function drawColorBar(){
    const c=document.getElementById('colorbar-canvas'),x=c.getContext('2d');
    for(let i=0;i<c.width;i++){x.fillStyle=colorScale(i/c.width);x.fillRect(i,0,1,c.height);}
    c.title=PALETTE_NAMES[paletteIdx]+' (click to change)';
    // also draw cluster colorbar
    const c2=document.getElementById('cl-colorbar-canvas');
    if(c2){const x2=c2.getContext('2d');
        for(let i=0;i<c2.width;i++){x2.fillStyle=colorScale(i/c2.width);x2.fillRect(i,0,1,c2.height);}
        c2.title=PALETTE_NAMES[paletteIdx]+' (click to change)';
    }
}

// sync color range: hist config > defaults. Only called on metric change (not per-event).
function syncRangeFromHist(){
    const mt=document.getElementById('color-metric').value;
    // restore user edits for this metric if they exist
    if(rangeUserOverrides[mt]){
        rangeMin=rangeUserOverrides[mt][0];
        rangeMax=rangeUserOverrides[mt][1];
        rangeUserEdited=true;
        updateRangeDisplay();
        return;
    }
    rangeUserEdited=false;
    const h=histConfig;
    if(mt==='integral' && h.bin_min!==undefined){
        rangeMin=h.bin_min; rangeMax=h.bin_max;
    } else if(mt==='time' && h.pos_min!==undefined){
        rangeMin=h.pos_min; rangeMax=h.pos_max;
    } else {
        const def=RANGE_DEFAULTS[mt]||[null,null];
        rangeMin=def[0]; rangeMax=def[1];
    }
    updateRangeDisplay();
}

function updateTimeCutLabel(){
    const el=document.getElementById('tcut-label');
    if(!el) return;
    const h=histConfig;
    if(h.time_min!==undefined && h.time_max!==undefined)
        el.textContent=`Time Cut [${h.time_min}, ${h.time_max}] ns`;
    else
        el.textContent='Time Cut';
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
let geoViewInit=false;  // has fitView been called at least once?
let geoDragging=false, geoDragStartX=0, geoDragStartY=0, geoDragX=0, geoDragY=0;
let geoDragMoved=false;  // true if drag exceeded click threshold
const GEO_DRAG_THRESHOLD=4;  // pixels — movement below this is treated as click, not drag

function initGeo(){
    geoWrap=document.getElementById('geo-wrap');
    geoCanvas=document.getElementById('geo-canvas');
    geoCtx=geoCanvas.getContext('2d');
    resizeGeo();
    new ResizeObserver(resizeGeo).observe(geoWrap);

    // zoom with mouse wheel
    geoCanvas.addEventListener('wheel',e=>{
        e.preventDefault();
        const r=geoCanvas.getBoundingClientRect();
        const cx=e.clientX-r.left, cy=e.clientY-r.top;
        const factor=e.deltaY<0?1.15:1/1.15;
        offsetX=cx-(cx-offsetX)*factor;
        offsetY=cy-(cy-offsetY)*factor;
        scale*=factor;
        redrawGeo();
    },{passive:false});

    // pan with left-click drag (with threshold to distinguish from click)
    geoCanvas.addEventListener('mousedown',e=>{
        if(e.button===0){
            geoDragging=true; geoDragMoved=false;
            geoDragStartX=e.clientX; geoDragStartY=e.clientY;
            geoDragX=e.clientX; geoDragY=e.clientY;
            e.preventDefault();
        }
    });
    window.addEventListener('mousemove',e=>{
        if(!geoDragging) return;
        const dx=e.clientX-geoDragStartX, dy=e.clientY-geoDragStartY;
        if(!geoDragMoved && Math.abs(dx)<GEO_DRAG_THRESHOLD && Math.abs(dy)<GEO_DRAG_THRESHOLD) return;
        if(!geoDragMoved){ geoDragMoved=true; geoCanvas.style.cursor='grabbing'; }
        offsetX+=e.clientX-geoDragX;
        offsetY+=e.clientY-geoDragY;
        geoDragX=e.clientX; geoDragY=e.clientY;
        redrawGeo();
    });
    window.addEventListener('mouseup',e=>{
        if(!geoDragging) return;
        geoDragging=false;
        geoCanvas.style.cursor='';
        // if not dragged beyond threshold, treat as a click
        if(!geoDragMoved){
            const r=geoCanvas.getBoundingClientRect();
            geoHandleClick(e.clientX-r.left, e.clientY-r.top);
        }
    });

    // double-click to reset view
    geoCanvas.addEventListener('dblclick',e=>{
        e.preventDefault();
        if(modules.length) fitView();
        redrawGeo();
    });

    // reset view button
    document.getElementById('btn-reset-view').onclick=()=>{
        if(modules.length) fitView();
        redrawGeo();
    };
}

function resetGeoView(){
    if(modules.length) fitView();
    redrawGeo();
}
function resizeGeo(){
    canvasW=geoWrap.clientWidth; canvasH=geoWrap.clientHeight;
    if(canvasW<10||canvasH<10)return;
    geoCanvas.width=canvasW; geoCanvas.height=canvasH;
    if(modules.length && !geoViewInit){ fitView(); geoViewInit=true; }
    redrawGeo();
}
function fitView(){
    const m=15;let x0=1e9,x1=-1e9,y0=1e9,y1=-1e9;
    for(const d of modules){x0=Math.min(x0,d.x-d.sx/2);x1=Math.max(x1,d.x+d.sx/2);y0=Math.min(y0,d.y-d.sy/2);y1=Math.max(y1,d.y+d.sy/2);}
    scale=Math.min((canvasW-2*m)/(x1-x0),(canvasH-2*m)/(y1-y0));
    offsetX=canvasW/2-(x0+x1)/2*scale; offsetY=canvasH/2+(y0+y1)/2*scale;
}
function redrawGeo(){
    if(activeTab==='cluster') drawClusterGeo(); else drawGeo();
}
function geoHandleClick(cx,cy){
    const m=hitTest(cx,cy);
    if(!m) return;
    if(activeTab==='cluster'){
        selectedModule=null;
        if(!clusterData || !clusterData.clusters || !clusterData.clusters.length){
            selectedCluster=-1;
            drawClusterGeo(); updateClusterTable(); showClusterDetail();
            return;
        }
        const idx=modules.indexOf(m);
        const clusters=clusterData.clusters;
        let found=-1;
        for(let ci=0;ci<clusters.length;ci++){
            if(clusters[ci].modules&&clusters[ci].modules.includes(idx)){ found=ci; break; }
        }
        if(found<0) selectedCluster=-1;
        else selectedCluster=(selectedCluster===found)?-1:found;
        document.getElementById('cl-select').value=selectedCluster>=0?selectedCluster:'all';
        drawClusterGeo(); updateClusterTable(); showClusterDetail();
    } else {
        showWaveform(m);
    }
}
function d2c(x,y){return[x*scale+offsetX,-y*scale+offsetY];}
function c2d(cx,cy){return[(cx-offsetX)/scale,-(cy-offsetY)/scale];}
// check if time cut checkbox is active and config exists
function isTimeCut(){
    return document.getElementById('time-cut').checked
        && histConfig.time_min!==undefined && histConfig.time_max!==undefined;
}

// filter peaks by time cut if active
function peaksInCut(peaks){
    if(!peaks||!peaks.length) return [];
    if(!isTimeCut()) return peaks;
    const tmin=histConfig.time_min, tmax=histConfig.time_max;
    return peaks.filter(p=>p.t>=tmin && p.t<=tmax);
}

// tallest peak from a list
function tallest(peaks){
    if(!peaks||!peaks.length) return null;
    let best=peaks[0];
    for(let i=1;i<peaks.length;i++) if(peaks[i].h>best.h) best=peaks[i];
    return best;
}

function modVal(m){
    const key=`${m.roc}_${m.sl}_${m.ch}`;
    const mt=document.getElementById('color-metric').value;
    if(mt==='occupancy'){
        if(occTotal<=0) return null;
        const src=isTimeCut()?occTcutData:occData;
        return 100.0*(src[key]||0)/occTotal;
    }
    const d=eventChannels[key];
    if(!d)return null;
    if(mt==='pedestal')return d.pm||0;
    const pks=peaksInCut(d.pk);
    if(mt==='count') return pks.length;
    const bp=tallest(pks);
    if(!bp)return null;
    if(mt==='height')return bp.h;
    if(mt==='time')return bp.t;
    return bp.i;
}
function drawGeo(){
    if(!geoCtx)return;const ctx=geoCtx;ctx.clearRect(0,0,canvasW,canvasH);
    const useLog=document.getElementById('log-scale').checked;
    const vals=modules.map(modVal);
    const vmin=rangeMin!==null?rangeMin:0;
    const vmax=rangeMax!==null?rangeMax:100;
    const span=vmax-vmin||1;

    for(let i=0;i<modules.length;i++){
        const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale,v=vals[i];
        let t=0;
        if(v!==null){
            const clamped=Math.max(vmin,Math.min(vmax,v));
            if(useLog) t=Math.log1p(clamped-vmin)/Math.log1p(span);
            else t=(clamped-vmin)/span;
        }
        ctx.fillStyle=(v!==null)?colorScale(t):(m.t==='G'?'#1a1a2e':'#12122a');
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
    for(const id of['waveform-div','inthist-div','poshist-div','cl-energy-hist'])
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
        currentWaveform=null;
        Plotly.react('waveform-div',[],{...PL,title:{text:`${mod.n} — No data`,font:{size:11,color:'#555'}}},PC2);
        document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No data</td></tr>';
        showHistograms(mod); drawGeo(); return;
    }

    const samples=d.s, peaks=d.pk||[], x=samples.map((_,i)=>i);
    currentWaveform={x, y:Array.from(samples)};
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
        currentHist[divId]=null;
        Plotly.react(divId,[],{...PL,title:{text:'--hist not enabled',font:{size:10,color:'#444'}}},PC2);
        return;
    }
    fetch(url).then(r=>r.json()).then(data=>{
        if(data.error||!data.bins||!data.bins.length){
            currentHist[divId]=null;
            Plotly.react(divId,[],{...PL,title:{text:`${title} — No data`,font:{size:10,color:'#555'}}},PC2);
            return;
        }
        const x=data.bins.map((_,i)=>binMin+(i+0.5)*binStep);
        // store non-zero bins for copy
        const cx=[], cy=[];
        for(let i=0;i<data.bins.length;i++){if(data.bins[i]>0){cx.push(x[i]);cy.push(data.bins[i]);}}
        currentHist[divId]={x:cx,y:cy};

        const entries=data.bins.reduce((a,b)=>a+b,0)+data.underflow+data.overflow;
        const stats=`${data.events} evts | Entries: ${entries}  Under: ${data.underflow}  Over: ${data.overflow}`;
        const xMin=binMin, xMax=binMin+data.bins.length*binStep;
        Plotly.react(divId,[{
            x,y:data.bins,type:'bar',marker:{color:barColor,line:{width:0}},
            hovertemplate:'%{x:.0f}: %{y}<extra></extra>',
        }],{...PL,
            title:{text:`${title}<br><span style="font-size:9px;color:#888">${stats}</span>`,font:{size:10,color:'#ccc'}},
            xaxis:{...PL.xaxis,title:xTitle,range:[xMin,xMax]},yaxis:{...PL.yaxis,title:'Counts'},bargap:0.05,
        },PC2);
    }).catch(()=>{
        currentHist[divId]=null;
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
let eventRequestId = 0;  // increments on each fetch, stale responses ignored

function loadEventData(reqId, data) {
    if (reqId !== eventRequestId) return;  // stale response, discard
    if (data.error) {
        document.getElementById('status-bar').textContent = data.error;
        // refresh ring selector to remove stale entries
        if (mode === 'online') updateRingSelector();
        return;
    }
    currentEvent = data.event;
    eventChannels = data.channels || {};
    const nch = Object.keys(eventChannels).length;
    const npk = Object.values(eventChannels).reduce((s,c) => s + (c.pk||[]).length, 0);
    const modeTag = mode === 'online' ? ' [LIVE]' : '';
    document.getElementById('status-bar').textContent = `Event ${currentEvent}: ${nch} channels, ${npk} peaks${modeTag}`;
    if(activeTab==='cluster'){
        clusterEvent=-1; // invalidate cache
        loadClusterData(currentEvent);
    } else {
        drawGeo();
    }
    if (selectedModule) showWaveform(selectedModule);
}

function loadEvent(evnum) {
    currentEvent = evnum;
    const reqId = ++eventRequestId;
    if (mode === 'file') document.getElementById('ev-input').value = evnum;
    document.getElementById('status-bar').textContent = `Loading event ${evnum}...`;
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
            o.textContent = `Event ${ring[i]}` + (i === ring.length - 1 ? ' (latest)' : '');
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
        el.textContent = `⏸ Paused at event ${currentEvent} — click or press F to resume`;
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
                const now = Date.now();
                // throttle event display to ~5 Hz
                if (autoFollow && now - lastEventFetch > 200) {
                    lastEventFetch = now;
                    loadLatestEvent();
                }
                // throttle ring selector update to ~2 Hz
                if (now - lastRingFetch > 500) {
                    lastRingFetch = now;
                    updateRingSelector();
                }
                // throttle occupancy refresh to ~0.5 Hz
                if (now - lastOccFetch > 2000) {
                    lastOccFetch = now;
                    fetchOccupancy();
                }
            } else if (msg.type === 'status') {
                setEtStatus(msg.connected, msg.waiting, msg.retries);
            } else if (msg.type === 'hist_cleared') {
                occData={}; occTcutData={}; occTotal=0;
                if (selectedModule) showHistograms(selectedModule);
                drawGeo();
            }
        } catch (e) {}
    };
}

function clearHistograms() {
    fetch('/api/hist/clear').then(r => r.json()).then(data => {
        occData={}; occTcutData={}; occTotal=0;
        document.getElementById('status-bar').textContent = 'Entries cleared';
        if (selectedModule) showHistograms(selectedModule);
        drawGeo();
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
                if(cfg.cluster_hist){
                    clHistMin=cfg.cluster_hist.min||0;
                    clHistMax=cfg.cluster_hist.max||3000;
                    clHistStep=cfg.cluster_hist.step||10;
                }
                initClHist(); plotClHist();
                updateTimeCutLabel();
                g_histCheckbox = histEnabled;
                const hcb = document.getElementById('hist-checkbox');
                if (hcb) hcb.checked = histEnabled;
                document.getElementById('ev-total').textContent = `/ ${totalEvents}`;
                updateHeaderInfo(cfg);
                if (histEnabled) { fetchOccupancy(); fetchClHist(); }
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
        if (mt === 'occupancy') {
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
// Tab switching
// =========================================================================
function switchTab(tab){
    if(tab===activeTab) return;
    activeTab=tab;
    document.querySelectorAll('.tab').forEach(t=>{
        t.classList.toggle('active', t.dataset.tab===tab);
    });
    const isDQ = tab==='dq';
    document.getElementById('geo-toolbar-dq').style.display = isDQ ? 'flex' : 'none';
    document.getElementById('geo-toolbar-cl').style.display = isDQ ? 'none' : 'flex';
    document.getElementById('detail-panel').style.display = isDQ ? 'flex' : 'none';
    document.getElementById('cluster-panel').style.display = isDQ ? 'none' : 'flex';
    if(tab==='cluster') {
        selectedModule=null;  // clear DQ selection highlight
        loadClusterData(currentEvent);
        // Plotly can't size in hidden divs — resize now that panel is visible
        setTimeout(()=>{
            try{Plotly.Plots.resize('cl-energy-hist');}catch(e){}
            plotClHist();
        }, 50);
    } else {
        drawGeo();
    }
}

// =========================================================================
// Clustering
// =========================================================================
function loadClusterData(evnum){
    if(clusterEvent===evnum && clusterData) { drawClusterGeo(); updateClusterTable(); return; }
    document.getElementById('status-bar').textContent=`Loading clusters for event ${evnum}...`;
    fetch(`/api/clusters/${evnum}`).then(r=>r.json()).then(data=>{
        if(data.error){ document.getElementById('status-bar').textContent=data.error; return; }
        clusterData=data;
        clusterEvent=evnum;
        selectedCluster=-1;
        drawClusterGeo();
        updateClusterUI();
        // accumulate energy histogram (only if no prebuilt data)
        if(!histEnabled && data.clusters && data.clusters.length>0){
            fillClHist(data.clusters);
            plotClHist();
        }
        const nc=data.clusters?data.clusters.length:0;
        const nh=data.hits?Object.keys(data.hits).length:0;
        document.getElementById('status-bar').textContent=`Event ${evnum}: ${nc} clusters, ${nh} hit modules`;
    }).catch(err=>{ document.getElementById('status-bar').textContent=`Error: ${err}`; });
}

function updateClusterUI(){
    if(!clusterData) return;
    const clusters=clusterData.clusters||[];
    // populate cluster selector
    const sel=document.getElementById('cl-select');
    sel.innerHTML='<option value="all">All ('+clusters.length+')</option>';
    clusters.forEach((cl,i)=>{
        const o=document.createElement('option');
        o.value=i;
        o.textContent=`#${i} ${cl.center} — ${cl.energy.toFixed(0)} MeV`;
        sel.appendChild(o);
    });
    // summary
    const totalE=clusters.reduce((s,c)=>s+c.energy,0);
    document.getElementById('cl-summary').textContent=
        `${clusters.length} clusters, ${totalE.toFixed(0)} MeV total`;
    updateClusterTable();
}

function updateClusterTable(){
    if(!clusterData) return;
    const clusters=clusterData.clusters||[];
    const tbody=document.getElementById('cl-tbody');
    let rows='';
    clusters.forEach((cl,i)=>{
        const sel=selectedCluster===i;
        const col=PC[i%PC.length];
        rows+=`<tr class="cl-table-row${sel?' selected':''}" data-idx="${i}" style="border-left:3px solid ${col}">
            <td style="text-align:center">${i}</td>
            <td>${cl.center}</td>
            <td>${cl.energy.toFixed(1)}</td>
            <td>${cl.x.toFixed(1)}</td>
            <td>${cl.y.toFixed(1)}</td>
            <td style="text-align:center">${cl.nblocks}</td>
        </tr>`;
    });
    if(!clusters.length) rows='<tr><td colspan="6" style="text-align:center;color:var(--dim);padding:8px">No clusters</td></tr>';
    tbody.innerHTML=rows;
    // click handlers
    tbody.querySelectorAll('.cl-table-row').forEach(tr=>{
        tr.onclick=()=>{
            const idx=parseInt(tr.dataset.idx);
            selectedCluster=(selectedCluster===idx)?-1:idx;
            document.getElementById('cl-select').value=selectedCluster>=0?selectedCluster:'all';
            drawClusterGeo();
            updateClusterTable();
            showClusterDetail();
        };
    });
}

function showClusterDetail(){
    const hdr=document.getElementById('cl-detail-header');
    if(selectedCluster<0 || !clusterData || !clusterData.clusters[selectedCluster]){
        hdr.innerHTML='<span class="cl-info-text">Click a module or select a cluster</span>';
        return;
    }
    const cl=clusterData.clusters[selectedCluster];
    const col=PC[selectedCluster%PC.length];
    hdr.innerHTML=`<span class="mod-name" style="color:${col}">Cluster #${selectedCluster}</span>
        <span class="mod-daq">Center: ${cl.center} (ID ${cl.center_id}) &middot;
        ${cl.energy.toFixed(1)} MeV &middot; (${cl.x.toFixed(1)}, ${cl.y.toFixed(1)}) &middot;
        ${cl.nblocks} blocks, ${cl.npos} pos</span>`;
}

// build a set of module indices belonging to a cluster
function clusterModuleSet(clIdx){
    if(!clusterData||clIdx<0) return null;
    const cl=clusterData.clusters[clIdx];
    if(!cl) return null;
    return new Set(cl.modules);
}

function drawClusterGeo(){
    if(!geoCtx) return;
    const ctx=geoCtx;
    ctx.clearRect(0,0,canvasW,canvasH);
    if(!clusterData){ drawGeo(); return; }

    const hits=clusterData.hits||{};
    const clusters=clusterData.clusters||[];
    const useLog=document.getElementById('cl-log-scale').checked;

    // find energy range from hits (auto range)
    let autoMax=0;
    for(const k in hits) if(hits[k]>autoMax) autoMax=hits[k];
    if(autoMax<=0) autoMax=1;
    const emin=clRangeMin!==null?clRangeMin:0;
    const emax=clRangeMax!==null?clRangeMax:autoMax;
    document.getElementById('cl-range-min-show').textContent=emin.toFixed(0);
    document.getElementById('cl-range-max-show').textContent=emax.toFixed(0);

    // build module-index → cluster-index map
    const modCluster={};
    clusters.forEach((cl,ci)=>{
        (cl.modules||[]).forEach(mi=>{ modCluster[mi]=ci; });
    });

    // selected cluster module set
    const selSet=selectedCluster>=0?clusterModuleSet(selectedCluster):null;

    for(let i=0;i<modules.length;i++){
        const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
        const energy=hits[String(i)]||0;
        const ci=modCluster[i];

        // dim modules not in selected cluster
        let dimmed=false;
        if(selSet && !selSet.has(i)) dimmed=true;

        let fillColor;
        const span=emax-emin||1;
        if(energy>0 && !dimmed){
            const clamped=Math.max(emin,Math.min(emax,energy));
            let t=useLog?Math.log1p(clamped-emin)/Math.log1p(span):(clamped-emin)/span;
            t=Math.max(0,Math.min(1,t));
            fillColor=colorScale(t);
        } else {
            fillColor=dimmed?'#0a0a18':(m.t==='G'?'#1a1a2e':'#12122a');
        }
        ctx.fillStyle=fillColor;
        ctx.fillRect(cx-w/2,cy-h/2,w,h);

        // border: highlight cluster members
        let strokeColor='#333', lw=0.5;
        if(ci!==undefined && !dimmed){
            strokeColor=PC[ci%PC.length];
            lw=1.5;
        }
        if(selectedCluster>=0 && ci===selectedCluster){
            strokeColor=PC[ci%PC.length];
            lw=2.5;
        }
        if(hoveredModule&&hoveredModule.n===m.n){ strokeColor='#00b4d8'; lw=1.5; }
        if(selectedModule&&selectedModule.n===m.n){ strokeColor='#fff'; lw=2.5; }
        ctx.strokeStyle=strokeColor; ctx.lineWidth=lw;
        ctx.strokeRect(cx-w/2,cy-h/2,w,h);
    }

    // draw cluster center markers (crosshairs)
    clusters.forEach((cl,ci)=>{
        if(selSet && ci!==selectedCluster) return;
        const [cx,cy]=d2c(cl.x,cl.y);
        const col=PC[ci%PC.length];
        ctx.strokeStyle=col; ctx.lineWidth=2;
        const sz=6;
        ctx.beginPath(); ctx.moveTo(cx-sz,cy); ctx.lineTo(cx+sz,cy); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(cx,cy-sz); ctx.lineTo(cx,cy+sz); ctx.stroke();
    });
}

// =========================================================================
// =========================================================================
// Cluster energy histogram (accumulated)
// =========================================================================
function initClHist(){
    const nbins=Math.max(1,Math.ceil((clHistMax-clHistMin)/clHistStep));
    clHistBins=new Array(nbins).fill(0);
    clHistEvents=0;
    currentClHist=null;
}
function clearClHist(){ initClHist(); plotClHist(); }

function fetchClHist(){
    fetch('/api/cluster_hist').then(r=>r.json()).then(data=>{
        if(!data.bins||!data.bins.length) return;
        // use server config if available
        if(data.min!==undefined) clHistMin=data.min;
        if(data.max!==undefined) clHistMax=data.max;
        if(data.step!==undefined) clHistStep=data.step;
        clHistBins=data.bins;
        clHistEvents=data.events||0;
        plotClHist();
    }).catch(()=>{});
}

function fillClHist(clusters){
    if(!clHistBins) initClHist();
    if(!clusters||!clusters.length) return;
    for(const cl of clusters){
        const b=Math.floor((cl.energy-clHistMin)/clHistStep);
        if(b>=0 && b<clHistBins.length) clHistBins[b]++;
    }
    clHistEvents++;
}

function plotClHist(){
    const div='cl-energy-hist';
    if(!clHistBins||!clHistBins.length){
        currentClHist=null;
        Plotly.react(div,[],{...PL,title:{text:'Cluster Energy — No data',font:{size:10,color:'#555'}}},PC2);
        return;
    }
    const x=clHistBins.map((_,i)=>clHistMin+(i+0.5)*clHistStep);
    const entries=clHistBins.reduce((a,b)=>a+b,0);
    // store non-zero for copy
    const cx=[],cy=[];
    for(let i=0;i<clHistBins.length;i++){if(clHistBins[i]>0){cx.push(x[i]);cy.push(clHistBins[i]);}}
    currentClHist={x:cx,y:cy};

    Plotly.react(div,[{
        x,y:clHistBins,type:'bar',marker:{color:'#ff922b',line:{width:0}},
        hovertemplate:'%{x:.0f} MeV: %{y}<extra></extra>',
    }],{...PL,
        title:{text:`Cluster Energy<br><span style="font-size:9px;color:#888">${clHistEvents} evts | ${entries} clusters</span>`,font:{size:10,color:'#ccc'}},
        xaxis:{...PL.xaxis,title:'Energy (MeV)',range:[clHistMin,clHistMax]},
        yaxis:{...PL.yaxis,title:'Counts'},bargap:0.05,
    },PC2);
}

// Init
// =========================================================================
function init(){
    drawColorBar(); initGeo();
    document.getElementById('colorbar-canvas').onclick=()=>{
        paletteIdx=(paletteIdx+1)%PALETTE_NAMES.length;
        drawColorBar(); drawGeo();
    };
    Plotly.newPlot('waveform-div',[],{...PL,xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC'}},PC2);
    Plotly.newPlot('inthist-div',[],{...PL,title:{text:'Integral Histogram',font:{size:10,color:'#555'}}},PC2);
    Plotly.newPlot('poshist-div',[],{...PL,title:{text:'Position Histogram',font:{size:10,color:'#555'}}},PC2);

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
    setupCopyBtn('btn-copy-inthist', ()=>currentHist['inthist-div']);
    setupCopyBtn('btn-copy-poshist', ()=>currentHist['poshist-div']);

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
        ()=>document.getElementById('inthist-div').parentElement,
        ()=>document.getElementById('top-panel'),
        ()=>0, 60, 60, resizeAllPlots);
    // 4. vertical inside bottom-panel: table ↔ waveform
    setupDivider('div-tw','x',
        ()=>document.getElementById('table-wrap'),
        ()=>document.getElementById('bottom-panel'),
        ()=>0, 150, 200, resizeAllPlots);

    // --- tab switching ---
    document.querySelectorAll('.tab').forEach(t=>{
        t.onclick=()=>switchTab(t.dataset.tab);
    });

    // --- cluster controls ---
    document.getElementById('cl-select').onchange=e=>{
        selectedCluster=e.target.value==='all'?-1:parseInt(e.target.value);
        drawClusterGeo(); updateClusterTable(); showClusterDetail();
    };
    document.getElementById('cl-log-scale').onchange=()=>{ if(activeTab==='cluster') drawClusterGeo(); };
    document.getElementById('cl-colorbar-canvas').onclick=()=>{
        paletteIdx=(paletteIdx+1)%PALETTE_NAMES.length;
        drawColorBar(); if(activeTab==='cluster') drawClusterGeo(); else drawGeo();
    };

    // cluster energy histogram
    Plotly.newPlot('cl-energy-hist',[],{...PL,title:{text:'Cluster Energy',font:{size:10,color:'#555'}}},PC2);
    setupCopyBtn('btn-copy-cl-hist', ()=>currentClHist);

    // cluster panel divider: histogram ↔ table
    setupDivider('div-cl-ht','y',
        ()=>document.getElementById('cl-hist-panel'),
        ()=>document.getElementById('cluster-panel'),
        ()=>0,
        80, 80, ()=>{try{Plotly.Plots.resize('cl-energy-hist');}catch(e){}});

    // --- file mode nav ---
    document.getElementById('btn-prev').onclick=()=>{if(currentEvent>1)loadEvent(currentEvent-1);};
    document.getElementById('btn-next').onclick=()=>{if(currentEvent<totalEvents)loadEvent(currentEvent+1);};
    document.getElementById('ev-input').onchange=e=>{const v=parseInt(e.target.value);if(v>=1&&v<=totalEvents)loadEvent(v);};
    document.getElementById('color-metric').onchange=()=>{syncRangeFromHist();drawGeo();};
    document.getElementById('log-scale').onchange=drawGeo;
    document.getElementById('time-cut').onchange=drawGeo;

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
        rangeUserEdited=true;
        const mt=document.getElementById('color-metric').value;
        rangeUserOverrides[mt]=[rangeMin,rangeMax];
        updateRangeDisplay(); drawGeo();
    }
    setupRangeEdit('range-min-btn','range-min-edit','range-min-show',
        ()=>rangeMin, v=>{rangeMin=v;}, dqRangeApply);
    setupRangeEdit('range-max-btn','range-max-edit','range-max-show',
        ()=>rangeMax, v=>{rangeMax=v;}, dqRangeApply);
    // Cluster range editors
    function clRangeApply(){ drawClusterGeo(); }
    setupRangeEdit('cl-range-min-btn','cl-range-min-edit','cl-range-min-show',
        ()=>clRangeMin, v=>{clRangeMin=v;}, clRangeApply);
    setupRangeEdit('cl-range-max-btn','cl-range-max-edit','cl-range-max-show',
        ()=>clRangeMax, v=>{clRangeMax=v;}, clRangeApply);

    // --- online mode nav ---
    document.getElementById('ring-select').onchange=e=>{
        autoFollow=false; updateFollowStatus();
        loadEvent(parseInt(e.target.value));
    };
    document.getElementById('ring-select').onfocus=()=>{ updateRingSelector(); };
    document.getElementById('follow-status').onclick=()=>{ autoFollow=true; updateFollowStatus(); loadLatestEvent(); };
    document.getElementById('btn-clear-hist').onclick=clearHistograms;

    // geo mouse
    const tip=document.getElementById('geo-tooltip');
    geoCanvas.addEventListener('mousemove',e=>{
        const r=geoCanvas.getBoundingClientRect(),m=hitTest(e.clientX-r.left,e.clientY-r.top);
        if(m!==hoveredModule){
            hoveredModule=m;
            if(activeTab==='cluster') drawClusterGeo(); else drawGeo();
        }
        if(m){
            let t=`${m.n}  (${m.t==='G'?'PbGlass':'PbWO₄'})\n${crateName(m.roc)}  slot ${m.sl}  ch ${m.ch}`;
            if(activeTab==='cluster' && clusterData){
                const idx=modules.indexOf(m);
                const energy=clusterData.hits?clusterData.hits[String(idx)]:0;
                if(energy) t+=`\nEnergy: ${energy.toFixed(1)} MeV`;
                // find which cluster this module belongs to
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
            tip.textContent=t;tip.style.display='block';
            tip.style.left=(e.clientX-r.left+14)+'px';tip.style.top=(e.clientY-r.top-8)+'px';
        }else tip.style.display='none';
    });
    // click is now handled via geoHandleClick (called from mouseup when drag threshold not exceeded)
    geoCanvas.addEventListener('mouseleave',()=>{
        hoveredModule=null;tip.style.display='none';
        if(activeTab==='cluster') drawClusterGeo(); else drawGeo();
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
        // cluster energy histogram config
        if(data.cluster_hist){
            clHistMin=data.cluster_hist.min||0;
            clHistMax=data.cluster_hist.max||3000;
            clHistStep=data.cluster_hist.step||10;
        }
        initClHist();
        updateTimeCutLabel();
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
            if(histEnabled) { fetchOccupancy(); fetchClHist(); }
            syncRangeFromHist();
            geoViewInit=false; resizeGeo();
            if(totalEvents>0)loadEvent(1);
        } else {
            setEtStatus(data.et_connected||false);
            document.getElementById('header-info').textContent=
                `${modules.length} modules · ONLINE · ring ${data.ring_buffer_size||20}`;
            syncRangeFromHist();
            fetchOccupancy();
            resizeGeo();
            connectWebSocket();
            updateRingSelector();
            loadLatestEvent();
        }
    });
}
window.addEventListener('DOMContentLoaded',init);
