// =========================================================================
// State
// =========================================================================
let modules=[], totalEvents=0, currentEvent=1;
let currentEventNumber=0, currentTriggerBits=0;  // DAQ event number + trigger from last loaded event
let eventChannels={};
let selectedModule=null, hoveredModule=null;
let geoCanvas, geoCtx, geoWrap, scale=1, offsetX=0, offsetY=0, canvasW, canvasH;
let geoOutlineCanvas, geoOutlineCtx;  // static outline layer
const PC=['#00b4d8','#ff6b6b','#51cf66','#ffd43b','#cc5de8','#ff922b','#20c997','#f06595'];
const CRATE_NAME={0x80:'adchycal1',0x82:'adchycal2',0x84:'adchycal3',0x86:'adchycal4',0x88:'adchycal5',0x8a:'adchycal6',0x8c:'adchycal7',
    0x01:'PRadTS',0x04:'PRadROC_1',0x05:'PRadROC_2',0x06:'PRadROC_3',0x07:'PRadSRS_1',0x08:'PRadSRS_2'};
function crateName(r){return CRATE_NAME[r]||`ROC 0x${r.toString(16)}`;}
let histEnabled=false, histConfig={};
let mode='idle';    // 'idle', 'file', or 'online'
let etAvailable=false, fileAvailable=false;
let ws=null;        // WebSocket connection (always connected)
let autoFollow=true; // auto-load latest event
let lastEventFetch=0, lastHistFetch=0, lastRingFetch=0, lastOccFetch=0, lastLmsFetch=0;
let refreshEventMs=200, refreshRingMs=500, refreshHistMs=2000, refreshLmsMs=2000;

// occupancy data (fetched once per file load when histograms enabled)
let occData={}, occTcutData={}, occTotal=0;
let currentWaveform=null;  // {x:[], y:[]} for copy button
let wfStackEnabled=false;
let wfStackTraces=[];      // [{x,y},...] accumulated waveforms
let wfStackModKey='';      // module key for current stack (clear on module change)
let currentHist={};  // {divId: {x:[], y:[]}} for histogram copy

// color range: per-tab user overrides, keyed by "tab:metric"
// Each entry is [min, max] where null = auto
const geoRangeOverrides={};
function geoRangeKey(tab, metric){ return tab+':'+metric; }
function getGeoRange(tab, metric){
    const k=geoRangeKey(tab, metric);
    return geoRangeOverrides[k] || [null, null];
}
function setGeoRange(tab, metric, min, max){
    geoRangeOverrides[geoRangeKey(tab, metric)] = [min, max];
}

// convenience: current tab's range
function curRange(){
    const tab=activeTab;
    let metric='default';
    if(tab==='dq') metric=document.getElementById('color-metric').value;
    else if(tab==='cluster') metric='energy';
    else if(tab==='lms') metric=document.getElementById('lms-color-metric').value;
    return getGeoRange(tab, metric);
}

// --- clustering tab state ---
let activeTab='dq';  // 'dq' or 'cluster'
let clusterData=null;  // {hits:{}, clusters:[]}
let selectedCluster=-1;  // -1 = all
let clusterEvent=-1;  // event number for cached cluster data

// cluster energy histogram (accumulated on frontend)
let clHistBins=null, clHistEvents=0;
let clHistMin=0, clHistMax=3000, clHistStep=10;
let currentClHist=null;  // {x:[], y:[]} for copy button
let currentNclustHist=null, currentNblocksHist=null;

// cluster count histograms (configurable via config.json clustering section)
let nclustBins=null, nblocksBins=null;
let nclustMin=0, nclustMax=20, nclustStep=1;
let nblocksMin=0, nblocksMax=40, nblocksStep=1;

// DQ tab working range (set by syncDqRange, used by geoDq)
let rangeMin=null, rangeMax=null;
let updateGeoTooltip=()=>{};  // set in init(), called on data refresh

// Light theme flag — set by report.js captureGeoForTab for print-friendly rendering
let geoLightTheme=false;
function geoEmptyColor(type){ return geoLightTheme?(type==='G'?'#e8e8f0':'#dde'):(type==='G'?'#1a1a2e':'#12122a'); }
function geoDimColor(){ return geoLightTheme?'#d0d0dd':'#0a0a18'; }
function geoStrokeColor(){ return geoLightTheme?'#aaa':'#333'; }

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
    for(const cid of ['cl-colorbar-canvas','lms-colorbar-canvas']){
        const c2=document.getElementById(cid);
        if(c2){const x2=c2.getContext('2d');
            for(let i=0;i<c2.width;i++){x2.fillStyle=colorScale(i/c2.width);x2.fillRect(i,0,1,c2.height);}
            c2.title=PALETTE_NAMES[paletteIdx]+' (click to change)';
        }
    }
}

// load DQ color range for current metric from overrides.
function syncDqRange(){
    const mt=document.getElementById('color-metric').value;
    const r=getGeoRange('dq', mt);
    rangeMin=r[0]; rangeMax=r[1];
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
    geoOutlineCanvas=document.getElementById('geo-outline-canvas');
    geoOutlineCtx=geoOutlineCanvas.getContext('2d');
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
    geoOutlineCanvas.width=canvasW; geoOutlineCanvas.height=canvasH;
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
    if(activeTab==='cluster') geoCluster();
    else if(activeTab==='lms') geoLms();
    else geoDq();
}

// ── Tab-specific geo providers ───────────────────────────────────────
// Helper: map a value into a color scale with optional log scaling
function geoValueColor(val, vmin, vmax, useLog, emptyColor){
    if(val===null||val===undefined) return emptyColor;
    const span=vmax-vmin||1;
    const clamped=Math.max(vmin,Math.min(vmax,val));
    const t=useLog?Math.log1p(clamped-vmin)/Math.log1p(span):(clamped-vmin)/span;
    return colorScale(Math.max(0,Math.min(1,t)));
}

function geoDq(){
    const useLog=document.getElementById('log-scale').checked;
    const vals=modules.map(modVal);
    const vmin=rangeMin!==null?rangeMin:0;
    const vmax=rangeMax!==null?rangeMax:100;

    renderGeo(
        i => vals[i]!==null ? geoValueColor(vals[i],vmin,vmax,useLog) : geoEmptyColor(modules[i].t),
        i => (selectedModule&&selectedModule.n===modules[i].n) ? {color:'#fff',width:2.5} : null,
        null
    );
}

function geoCluster(){
    if(!clusterData){ renderGeo(i=>geoEmptyColor(modules[i].t),null,null); return; }
    const hits=clusterData.hits||{};
    const clusters=clusterData.clusters||[];
    const useLog=document.getElementById('cl-log-scale').checked;

    let autoMax=0;
    for(const k in hits) if(hits[k]>autoMax) autoMax=hits[k];
    if(autoMax<=0) autoMax=1;
    const clr=getGeoRange('cluster','energy');
    const emin=clr[0]!==null?clr[0]:0;
    const emax=clr[1]!==null?clr[1]:autoMax;
    document.getElementById('cl-range-min-show').textContent=emin.toFixed(0);
    document.getElementById('cl-range-max-show').textContent=emax.toFixed(0);

    const modCluster={};
    clusters.forEach((cl,ci)=>{ (cl.modules||[]).forEach(mi=>{ modCluster[mi]=ci; }); });
    const selSet=selectedCluster>=0?clusterModuleSet(selectedCluster):null;

    renderGeo(
        i => {
            const energy=hits[String(i)]||0;
            const dimmed=selSet&&!selSet.has(i);
            if(energy>0&&!dimmed) return geoValueColor(energy,emin,emax,useLog);
            return dimmed?geoDimColor():geoEmptyColor(modules[i].t);
        },
        i => {
            if(selectedModule&&selectedModule.n===modules[i].n) return {color:'#fff',width:2.5};
            const ci=modCluster[i];
            const dimmed=selSet&&!selSet.has(i);
            if(ci!==undefined&&!dimmed){
                const w=(selectedCluster>=0&&ci===selectedCluster)?2.5:1.5;
                return {color:PC[ci%PC.length],width:w};
            }
            return null;
        },
        ctx => {
            clusters.forEach((cl,ci)=>{
                if(selSet&&ci!==selectedCluster) return;
                const [cx,cy]=d2c(cl.x,cl.y);
                ctx.strokeStyle=PC[ci%PC.length]; ctx.lineWidth=2;
                const sz=6;
                ctx.beginPath();ctx.moveTo(cx-sz,cy);ctx.lineTo(cx+sz,cy);ctx.stroke();
                ctx.beginPath();ctx.moveTo(cx,cy-sz);ctx.lineTo(cx,cy+sz);ctx.stroke();
            });
        }
    );
}

function geoLms(){
    const metric=document.getElementById('lms-color-metric').value;
    const useLog=document.getElementById('lms-log-scale').checked;
    const mods=lmsSummaryData?lmsSummaryData.modules:{};

    const lmsVal=md=>{
        if(!md) return null;
        if(metric==='mean') return md.mean;
        if(metric==='rms_frac') return md.mean>0?md.rms/md.mean:0;
        return md.warn?1:0;
    };

    let autoMax=0;
    for(const k in mods){ const v=lmsVal(mods[k]); if(v>autoMax) autoMax=v; }
    if(autoMax<=0) autoMax=1;
    const lmsr=getGeoRange('lms',metric);
    const vmin=lmsr[0]!==null?lmsr[0]:0;
    const vmax=lmsr[1]!==null?lmsr[1]:autoMax;
    document.getElementById('lms-range-min-show').textContent=vmin.toFixed(metric==='rms_frac'?3:0);
    document.getElementById('lms-range-max-show').textContent=vmax.toFixed(metric==='rms_frac'?3:0);

    renderGeo(
        i => {
            const md=mods[String(i)];
            const val=lmsVal(md);
            if(val!==null&&val>0){
                if(metric==='warn') return md.warn?'#f66':'#51cf66';
                return geoValueColor(val,vmin,vmax,useLog);
            }
            return geoEmptyColor(modules[i].t);
        },
        i => {
            if(lmsSelectedModule===i) return {color:'#fff',width:2.5};
            const md=mods[String(i)];
            if(md&&md.warn) return {color:'#f66',width:1.5};
            return null;
        },
        null
    );
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
// ── Unified geo renderer ─────────────────────────────────────────────
// colorFn(i)    → fill color string for module i
// outlineFn(i)  → null (default stroke) or {color, width} for special highlight
// decorateFn(ctx) → draw tab-specific extras (e.g. cluster crosshairs) on overlay
//
// Fills go on geoCtx, outlines go on geoOutlineCtx.
// Hover highlight is handled automatically for all tabs.

let _geoOutlineFn=null, _geoDecorateFn=null;  // stored for hover-only redraws

function renderGeoFills(colorFn){
    if(!geoCtx)return;
    const ctx=geoCtx;
    if(geoLightTheme){ctx.fillStyle='#fff';ctx.fillRect(0,0,canvasW,canvasH);}
    else ctx.clearRect(0,0,canvasW,canvasH);
    for(let i=0;i<modules.length;i++){
        const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
        ctx.fillStyle=colorFn(i);
        ctx.fillRect(cx-w/2,cy-h/2,w,h);
    }
}

function renderGeoOutlines(outlineFn, decorateFn){
    if(!geoOutlineCtx)return;
    const ctx=geoOutlineCtx;
    ctx.clearRect(0,0,canvasW,canvasH);

    // 1. batch default outlines in one path
    ctx.beginPath(); ctx.strokeStyle=geoStrokeColor(); ctx.lineWidth=0.5;
    for(let i=0;i<modules.length;i++){
        const style=outlineFn?outlineFn(i):null;
        const hov=hoveredModule&&hoveredModule.n===modules[i].n;
        if(!style&&!hov){
            const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
            ctx.rect(cx-w/2,cy-h/2,w,h);
        }
    }
    ctx.stroke();

    // 2. special outlines (cluster members, warn, selection)
    for(let i=0;i<modules.length;i++){
        const style=outlineFn?outlineFn(i):null;
        if(style){
            const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
            ctx.strokeStyle=style.color; ctx.lineWidth=style.width;
            ctx.strokeRect(cx-w/2,cy-h/2,w,h);
        }
    }

    // 3. hover highlight (same for all tabs)
    if(hoveredModule){
        const m=hoveredModule,[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
        ctx.strokeStyle='#00b4d8'; ctx.lineWidth=1.5;
        ctx.strokeRect(cx-w/2,cy-h/2,w,h);
    }

    // 4. tab-specific decorations
    if(decorateFn) decorateFn(ctx);
}

function renderGeo(colorFn, outlineFn, decorateFn){
    _geoOutlineFn=outlineFn; _geoDecorateFn=decorateFn;
    renderGeoFills(colorFn);
    renderGeoOutlines(outlineFn, decorateFn);
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

    if(!d){
        currentWaveform=null;
        Plotly.react('waveform-div',[],{...PL,title:{text:`${mod.n} — No data`,font:{size:11,color:'#555'}}},PC2);
        document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No data</td></tr>';
        showHistograms(mod); geoDq(); return;
    }

    // if samples are already present (e.g. from ring buffer), use them directly;
    // otherwise fetch on demand from /api/waveform/<event>/<key>
    if(d.s){
        renderWaveform(mod, key, d, d.s);
    } else {
        fetch(`/api/waveform/${currentEvent}/${key}`).then(r=>r.json()).then(wf=>{
            if(wf.error){ renderWaveform(mod, key, d, null); return; }
            // cache samples in eventChannels so re-clicking doesn't re-fetch
            d.s=wf.s;
            if(wf.pk) d.pk=wf.pk;
            if(wf.pm!==undefined) d.pm=wf.pm;
            if(wf.pr!==undefined) d.pr=wf.pr;
            renderWaveform(mod, key, d, d.s);
        }).catch(()=>renderWaveform(mod, key, d, null));
    }

    showHistograms(mod); geoDq();
}

function renderWaveform(mod, key, d, samples){
    if(!samples){
        currentWaveform=null;
        Plotly.react('waveform-div',[],{...PL,title:{text:`${mod.n} — No samples`,font:{size:11,color:'#555'}}},PC2);
        document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No waveform data</td></tr>';
        return;
    }

    const peaks=d.pk||[], x=samples.map((_,i)=>i);
    currentWaveform={x, y:Array.from(samples)};

    // --- stacking mode ---
    if(wfStackEnabled){
        // reset stack if module changed
        if(key!==wfStackModKey){ wfStackTraces=[]; wfStackModKey=key; }

        // add current waveform to stack
        wfStackTraces.push({x:Array.from(x), y:Array.from(samples)});

        // enforce max stack size
        const maxStack=200;
        while(wfStackTraces.length>maxStack) wfStackTraces.shift();

        // build traces: all stacked waveforms in dim color, no peaks
        const traces=wfStackTraces.map(w=>({
            x:w.x, y:w.y, type:'scatter', mode:'lines',
            line:{color:'rgba(119,119,170,0.25)', width:1},
            showlegend:false, hoverinfo:'skip',
        }));
        // highlight the latest waveform
        if(wfStackTraces.length>0){
            const last=wfStackTraces[wfStackTraces.length-1];
            traces.push({x:last.x, y:last.y, type:'scatter', mode:'lines',
                name:'Latest', line:{color:'#7777aa', width:1.5}, showlegend:false});
        }

        document.getElementById('wf-stack-count').textContent=`${wfStackTraces.length}/${maxStack}`;
        Plotly.react('waveform-div',traces,{...PL,
            title:{text:`${mod.n} — Stacked (${wfStackTraces.length})`,font:{size:11,color:'#ccc'}},
            xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC'},
        },PC2);

        // skip peaks table in stack mode
        document.getElementById('peaks-tbody').innerHTML=
            '<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">Stack mode — peaks hidden</td></tr>';
        return;
    }

    // --- normal (single event) mode ---
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
}

// =========================================================================
// Histograms
// =========================================================================
function fetchAndPlotHist(divId, url, title, xTitle, binMin, binStep, barColor, logXId, logYId){
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
            xaxis:{...PL.xaxis,title:xTitle,range:[xMin,xMax],
                type:logXId&&document.getElementById(logXId).checked?'log':'linear'},
            yaxis:{...PL.yaxis,title:'Counts',
                type:logYId&&document.getElementById(logYId).checked?'log':'linear'},
            bargap:0.05,
        },PC2);
    }).catch(()=>{
        currentHist[divId]=null;
        Plotly.react(divId,[],{...PL,title:{text:'Fetch error',font:{size:10,color:'#f66'}}},PC2);
    });
}

let lastHistModule = '';
function showHistograms(mod){
    const key=`${mod.roc}_${mod.sl}_${mod.ch}`;
    // in online mode, throttle auto-refreshes of the same module to ~1 Hz
    if (mode === 'online' && key === lastHistModule) {
        const now = Date.now();
        if (now - lastHistFetch < refreshHistMs) return;
        lastHistFetch = now;
    }
    lastHistFetch = Date.now();
    lastHistModule = key;
    const h=histConfig;
    fetchAndPlotHist('inthist-div',`/api/hist/${key}`,
        `${mod.n} Integral [${h.time_min||170}-${h.time_max||190} ns]`,
        'Peak Integral', h.bin_min||0, h.bin_step||100, '#00b4d8', 'inthist-logx', 'inthist-logy');
    fetchAndPlotHist('poshist-div',`/api/poshist/${key}`,
        `${mod.n} Peak Position`,
        'Time (ns)', h.pos_min||0, h.pos_step||4, '#51cf66');
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
function updateStatusBar(){
    const modeTag = mode === 'online' ? ' [LIVE]' : '';
    const trig = currentTriggerBits ? ` trig=0x${currentTriggerBits.toString(16)}` : '';
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
                if(cfg.nclusters_hist){
                    nclustMin=cfg.nclusters_hist.min||0;
                    nclustMax=cfg.nclusters_hist.max||20;
                    nclustStep=cfg.nclusters_hist.step||1;
                }
                if(cfg.nblocks_hist){
                    nblocksMin=cfg.nblocks_hist.min||0;
                    nblocksMax=cfg.nblocks_hist.max||40;
                    nblocksStep=cfg.nblocks_hist.step||1;
                }
                initClHist(); plotClHist(); plotClStatHists();
                updateTimeCutLabel();
                g_histCheckbox = histEnabled;
                const hcb = document.getElementById('hist-checkbox');
                if (hcb) hcb.checked = histEnabled;
                document.getElementById('ev-total').textContent = `/ ${totalEvents}`;
    
                updateHeaderStats();
                if (histEnabled) { fetchOccupancy(); fetchClHist(); }
                fetchEpicsChannels(); fetchEpicsLatest();
                if(activeTab==='epics') fetchAllEpicsSlots();
                if(activeTab==='physics') fetchPhysics();
                syncDqRange();
                redrawGeo();
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

// =========================================================================
// Clustering
// =========================================================================
function loadClusterData(evnum){
    if(clusterEvent===evnum && clusterData) { geoCluster(); updateClusterTable(); return; }
    document.getElementById('status-bar').textContent=`Loading clusters for sample ${evnum}...`;
    fetch(`/api/clusters/${evnum}`).then(r=>{
        if(!r.ok) throw new Error('not available');
        return r.json();
    }).then(data=>{
        if(data.error){ document.getElementById('status-bar').textContent=data.error; return; }
        clusterData=data;
        clusterEvent=evnum;
        selectedCluster=-1;
        geoCluster();
        updateClusterUI();
        updateGeoTooltip();
        // accumulate energy histogram from per-event cluster data
        if(data.clusters && data.clusters.length>0){
            fillClHist(data.clusters);
            plotClHist(); plotClStatHists();
        }
        updateStatusBar();
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
            geoCluster();
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


// =========================================================================
// =========================================================================
// Cluster energy histogram (accumulated)
// =========================================================================
function initClHist(){
    const nbins=Math.max(1,Math.ceil((clHistMax-clHistMin)/clHistStep));
    clHistBins=new Array(nbins).fill(0);
    clHistEvents=0;
    currentClHist=null;
    nclustBins=new Array(Math.ceil((nclustMax-nclustMin)/nclustStep)).fill(0);
    nblocksBins=new Array(Math.ceil((nblocksMax-nblocksMin)/nblocksStep)).fill(0);
    currentNclustHist=null;
    currentNblocksHist=null;
}
function clearClHist(){ initClHist(); plotClHist(); plotClStatHists(); }

function fetchClHist(){
    fetch('/api/cluster_hist').then(r=>r.json()).then(data=>{
        if(!data.bins||!data.bins.length) return;
        if(data.min!==undefined) clHistMin=data.min;
        if(data.max!==undefined) clHistMax=data.max;
        if(data.step!==undefined) clHistStep=data.step;
        clHistBins=data.bins;
        clHistEvents=data.events||0;
        // nclusters/nblocks from server
        if(data.nclusters&&data.nclusters.bins&&data.nclusters.bins.length){
            nclustMin=data.nclusters.min||0; nclustMax=data.nclusters.max||20;
            nclustStep=data.nclusters.step||1; nclustBins=data.nclusters.bins;
        }
        if(data.nblocks&&data.nblocks.bins&&data.nblocks.bins.length){
            nblocksMin=data.nblocks.min||0; nblocksMax=data.nblocks.max||40;
            nblocksStep=data.nblocks.step||1; nblocksBins=data.nblocks.bins;
        }
        plotClHist(); plotClStatHists();
    }).catch(()=>{});
}

function fillClHist(clusters){
    if(!clHistBins) initClHist();
    if(!clusters||!clusters.length) return;
    // energy histogram
    for(const cl of clusters){
        const b=Math.floor((cl.energy-clHistMin)/clHistStep);
        if(b>=0 && b<clHistBins.length) clHistBins[b]++;
    }
    // number of clusters per event
    const nc=clusters.length;
    const nb1=Math.floor((nc-nclustMin)/nclustStep);
    if(nclustBins && nb1>=0 && nb1<nclustBins.length) nclustBins[nb1]++;
    // number of blocks per cluster
    for(const cl of clusters){
        const nbl=cl.nblocks||0;
        const nb2=Math.floor((nbl-nblocksMin)/nblocksStep);
        if(nblocksBins && nb2>=0 && nb2<nblocksBins.length) nblocksBins[nb2]++;
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
        xaxis:{...PL.xaxis,title:'Energy (MeV)',range:[clHistMin,clHistMax],
            type:document.getElementById('clhist-logx').checked?'log':'linear'},
        yaxis:{...PL.yaxis,title:'Counts',
            type:document.getElementById('clhist-logy').checked?'log':'linear'},
        bargap:0.05,
    },PC2);
}

function plotClStatHists(){
    // number of clusters histogram
    function plotStat(divId, bins, bmin, bstep, title, xTitle, color, copyVar){
        if(!bins||!bins.length){
            return null;
        }
        const x=bins.map((_,i)=>bmin+i*bstep);  // left edge = actual value for integer bins
        const entries=bins.reduce((a,b)=>a+b,0);
        const cx=[],cy=[];
        for(let i=0;i<bins.length;i++){if(bins[i]>0){cx.push(x[i]);cy.push(bins[i]);}}
        Plotly.react(divId,[{
            x,y:bins,type:'bar',marker:{color,line:{width:0}},
            hovertemplate:'%{x}: %{y}<extra></extra>',
        }],{...PL,
            title:{text:`${title}<br><span style="font-size:9px;color:#888">${entries} entries</span>`,font:{size:10,color:'#ccc'}},
            xaxis:{...PL.xaxis,title:xTitle,range:[bmin-0.5,bmin+bins.length*bstep-0.5]},
            yaxis:{...PL.yaxis,title:'Counts'},bargap:0.05,
        },PC2);
        return {x:cx,y:cy};
    }
    currentNclustHist=plotStat('cl-nclust-hist',nclustBins,nclustMin,nclustStep,
        'Clusters per Event','# Clusters','#00b4d8');
    currentNblocksHist=plotStat('cl-nblocks-hist',nblocksBins,nblocksMin,nblocksStep,
        'Blocks per Cluster','# Blocks','#51cf66');
}

// =========================================================================
// LMS monitoring
// =========================================================================
let g_lmsWarnThresh=0.1;
let g_lmsRefIndex=-1;  // -1 = None (no normalization)
let currentLmsData=null;  // {x:[], y:[]} for copy button
let lmsSummaryData=null;  // {modules:{idx:{name,mean,rms,count,warn}}, events}
let lmsSelectedModule=-1;

function fetchLmsSummary(){
    const refQ=g_lmsRefIndex>=0?`?ref=${g_lmsRefIndex}`:'';
    fetch(`/api/lms/summary${refQ}`).then(r=>r.json()).then(data=>{
        lmsSummaryData=data;
        geoLms();
        updateLmsTable();
        if(hoveredModule) updateGeoTooltip();
        // update tab dot if not currently on LMS tab
        if(activeTab!=='lms') updateLmsDot();
    }).catch(()=>{});
}

function updateLmsDot(){
    const dot=document.getElementById('lms-dot');
    if(!lmsSummaryData||!lmsSummaryData.modules){dot.className='tab-dot';return;}
    const hasWarn=Object.values(lmsSummaryData.modules).some(m=>m.warn);
    dot.className='tab-dot'+(hasWarn?' alert':'');
}

function fetchLmsHistory(modIdx, modName){
    const refQ=g_lmsRefIndex>=0?`?ref=${g_lmsRefIndex}`:'';
    fetch(`/api/lms/${modIdx}${refQ}`).then(r=>r.json()).then(data=>{
        if(!data.time||!data.time.length){
            currentLmsData=null;
            Plotly.react('lms-plot',[],{...PL,
                title:{text:`${modName} — No LMS data`,font:{size:10,color:'#555'}}},PC2);
            return;
        }
        currentLmsData={x:Array.from(data.time), y:Array.from(data.integral)};
        const vals=data.integral;
        const mean=vals.reduce((a,b)=>a+b,0)/vals.length;
        const warnHi=mean*(1+g_lmsWarnThresh);
        const warnLo=mean*(1-g_lmsWarnThresh);
        const tRange=[data.time[0],data.time[data.time.length-1]];

        Plotly.react('lms-plot',[
            {x:data.time, y:data.integral, type:'scatter', mode:'markers',
             marker:{color:'#ff922b',size:3}, name:'LMS integral'},
            {x:tRange, y:[mean,mean],
             type:'scatter', mode:'lines', line:{color:'#51cf66',width:1,dash:'dash'}, name:`Mean ${mean.toFixed(0)}`},
            {x:tRange, y:[warnHi,warnHi],
             type:'scatter', mode:'lines', line:{color:'#f66',width:1,dash:'dot'}, showlegend:false},
            {x:tRange, y:[warnLo,warnLo],
             type:'scatter', mode:'lines', line:{color:'#f66',width:1,dash:'dot'}, showlegend:false},
        ],{...PL,
            title:{text:`LMS — ${modName} (${data.events} pts)${g_lmsRefIndex>=0?' [ref corrected]':''}`,
                font:{size:10,color:'#ccc'}},
            xaxis:{...PL.xaxis,
                title:data.sync_unix
                    ?`Time (s) after ${new Date((data.sync_unix - data.sync_rel_sec)*1000).toISOString().replace('T',' ').slice(0,19)} UTC`
                    :'Time (s)'},
            yaxis:{...PL.yaxis,title:g_lmsRefIndex>=0?'Corrected Integral':'Integral'},
            legend:{x:1,y:1,xanchor:'right',bgcolor:'rgba(0,0,0,0.6)',font:{size:9}},
            margin:{...PL.margin,t:28,b:36},
        },PC2);

        document.getElementById('lms-info-text').innerHTML=
            `<span class="mod-name">${modName}</span> <span class="mod-daq">Mean: ${mean.toFixed(1)} | RMS: ${(Math.sqrt(vals.reduce((s,v)=>s+(v-mean)**2,0)/vals.length)).toFixed(1)} | ${data.events} pts</span>`;
    }).catch(()=>{});
}

function updateLmsTable(){
    const tbody=document.getElementById('lms-tbody');
    if(!lmsSummaryData||!lmsSummaryData.modules){
        tbody.innerHTML='<tr><td colspan="6" style="text-align:center;color:var(--dim);padding:8px">No LMS data</td></tr>';
        return;
    }
    // sort: warnings first, then by rms/mean descending
    const entries=Object.entries(lmsSummaryData.modules).map(([idx,m])=>({idx:parseInt(idx),...m}));
    entries.sort((a,b)=>{
        if(a.warn!==b.warn) return a.warn?-1:1;
        const ra=a.mean>0?a.rms/a.mean:0, rb=b.mean>0?b.rms/b.mean:0;
        return rb-ra;
    });
    let rows='';
    for(const e of entries){
        const rmsFrac=e.mean>0?(e.rms/e.mean*100).toFixed(1):'--';
        const sel=lmsSelectedModule===e.idx;
        rows+=`<tr class="cl-table-row${sel?' selected':''}" data-idx="${e.idx}">
            <td>${e.name}</td>
            <td>${e.mean.toFixed(1)}</td>
            <td>${e.rms.toFixed(2)}</td>
            <td>${rmsFrac}%</td>
            <td style="text-align:center">${e.count}</td>
            <td style="text-align:center">${e.warn?'<span class="lms-warn">WARN</span>':'<span class="lms-ok">OK</span>'}</td>
        </tr>`;
    }
    tbody.innerHTML=rows;
    tbody.querySelectorAll('.cl-table-row').forEach(tr=>{
        tr.onclick=()=>{
            const idx=parseInt(tr.dataset.idx);
            lmsSelectedModule=idx;
            const mod=modules.find(m=>modules.indexOf(m)===idx);
            const name=lmsSummaryData.modules[idx]?lmsSummaryData.modules[idx].name:'';
            fetchLmsHistory(idx, name);
            updateLmsTable();
            geoLms();
        };
    });
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

    // waveform stacking controls
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
    document.getElementById('inthist-logx').onchange=()=>{ if(selectedModule) showHistograms(selectedModule); };
    document.getElementById('inthist-logy').onchange=()=>{ if(selectedModule) showHistograms(selectedModule); };
    document.getElementById('clhist-logx').onchange=plotClHist;
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
    document.getElementById('file-filter').oninput = e => renderFileList(e.target.value);
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
            fetch('/api/mode/file',{method:'POST'});
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
    if(selectedModule) showHistograms(selectedModule);

    initClHist(); plotClHist(); plotClStatHists();
    clusterData=null; clusterEvent=-1; selectedCluster=-1;
    currentNclustHist=null; currentNblocksHist=null;
    document.getElementById('cl-select').innerHTML='<option value="all">All</option>';
    document.getElementById('cl-detail-header').innerHTML=
        '<span class="cl-info-text">Click a module or select a cluster</span>';
    document.getElementById('cl-tbody').innerHTML='';

    lmsSummaryData=null; lmsSelectedModule=-1; currentLmsData=null;
    Plotly.react('lms-plot',[],{...PL,title:{text:'LMS History',font:{size:10,color:'#555'}}},PC2);
    document.getElementById('lms-detail-header').innerHTML=
        '<span class="cl-info-text">Click a module to view LMS history</span>';
    document.getElementById('lms-tbody').innerHTML='';

    currentGemNclHist=null; currentGemThetaHist=null;
    clearEpicsFrontend();
    clearPhysicsFrontend();

    sampleCount=0;
    updateHeaderStats();
    redrawGeo();
    document.getElementById('status-bar').textContent='';
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
        toggleBtn.textContent=mode==='online'?'View Files':'Go Online';
    } else {
        toggleBtn.style.display='none';
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
        setEtStatus(data.et_connected||false);
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
