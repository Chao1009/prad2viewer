// geo.js — HyCal geometry canvas rendering (two-layer: fills + outlines)

let geoCanvas, geoCtx, geoWrap, scale=1, offsetX=0, offsetY=0, canvasW, canvasH;
let geoOutlineCanvas, geoOutlineCtx;  // static outline layer
let geoFocusPbWO4=false;  // dim LG modules to highlight PbWO4

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


// DQ tab working range (set by syncDqRange, used by geoDq)
let rangeMin=null, rangeMax=null;
let updateGeoTooltip=()=>{};  // set in init(), called on data refresh

// Light theme flag — set by report.js captureGeoForTab for print-friendly rendering
let geoLightTheme=false;
function geoEmptyColor(type){ return THEME.noData; }
function geoDimColor(){ return THEME.altBase; }
function geoStrokeColor(){ return THEME.border; }

// =========================================================================
// Color scale — click the colorbar to cycle palettes
// First two ("rainbow", "blue-yellow") match prad2hvmon's web monitor
// palettes (resources/monitor_geo_view.js) so the two viewers share a
// color language.  Cycling starts from "rainbow".
// =========================================================================
function _lerpHex(a,b,t){
    const ar=parseInt(a.slice(1,3),16),ag=parseInt(a.slice(3,5),16),ab=parseInt(a.slice(5,7),16);
    const br=parseInt(b.slice(1,3),16),bg=parseInt(b.slice(3,5),16),bb=parseInt(b.slice(5,7),16);
    return `rgb(${Math.round(ar+(br-ar)*t)},${Math.round(ag+(bg-ag)*t)},${Math.round(ab+(bb-ab)*t)})`;
}
const PALETTES={
    rainbow(t){
        if(t<0.25) return _lerpHex('#1e3a5f','#3b82f6',t*4);
        if(t<0.50) return _lerpHex('#3b82f6','#2dd4a0',(t-0.25)*4);
        if(t<0.75) return _lerpHex('#2dd4a0','#eab308',(t-0.50)*4);
        return _lerpHex('#eab308','#f56565',(t-0.75)*4);
    },
    'blue-yellow'(t){
        if(t<0.5) return _lerpHex('#0b1628','#3b9eff',t*2);
        return _lerpHex('#3b9eff','#eab308',(t-0.5)*2);
    },
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
    // Values are shown in the Cut-Settings dialog (cut_dialog.js).
    // Kept as a no-op so call sites in config.js / mode-switch don't break.
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

    // Focus PbWO4 checkbox
    document.getElementById('focus-pbwo4').onchange=function(){
        geoFocusPbWO4=this.checked; redrawGeo();
    };
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
        i => (selectedModule&&selectedModule.n===modules[i].n) ? {color:THEME.selectBorder,width:2.5} : null,
        null
    );
}

function d2c(x,y){return[x*scale+offsetX,-y*scale+offsetY];}
function c2d(cx,cy){return[(cx-offsetX)/scale,-(cy-offsetY)/scale];}
// Time cut on the geo color map follows the *show* checkbox (visual decision)
// — separate from the server-side "apply" which gates histogram fill.
function isTimeCut(){
    const cb = document.getElementById('cut-show');
    if (!cb || !cb.checked) return false;
    const f = histConfig && histConfig.waveform_filter
        && histConfig.waveform_filter.time;
    return !!(f && (f.min != null || f.max != null));
}

// filter peaks by time cut if active
function peaksInCut(peaks){
    if(!peaks||!peaks.length) return [];
    if(!isTimeCut()) return peaks;
    const f = histConfig.waveform_filter.time;
    const tmin = f.min != null ? f.min : -Infinity;
    const tmax = f.max != null ? f.max :  Infinity;
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
    // 'integral' picks the peak with the largest integral.  Time window
    // is honoured only when the show-cut checkbox is on AND the filter
    // has a time range, so the color map matches what the user explicitly
    // asked for (overlay off ⇒ show the full event).  The peaks coming
    // from the analyzer are already gated by min_peak_height, so no
    // separate threshold check is needed here.
    if(mt==='integral'){
        if(!d.pk||!d.pk.length) return null;
        const useTcut = isTimeCut();
        const f = useTcut ? histConfig.waveform_filter.time : null;
        const tmin = (f && f.min != null) ? f.min : undefined;
        const tmax = (f && f.max != null) ? f.max : undefined;
        let best=-1;
        for(const p of d.pk){
            if(tmin!==undefined && p.t<tmin) continue;
            if(tmax!==undefined && p.t>tmax) continue;
            if(p.i>best) best=p.i;
        }
        return best>=0?best:null;
    }
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
    if(geoLightTheme){ctx.fillStyle=THEME.bg;ctx.fillRect(0,0,canvasW,canvasH);}
    else ctx.clearRect(0,0,canvasW,canvasH);
    const focus=geoFocusPbWO4;
    for(let i=0;i<modules.length;i++){
        const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
        if(focus) ctx.globalAlpha=(m.t==='G')?0.15:1;
        ctx.fillStyle=colorFn(i);
        ctx.fillRect(cx-w/2,cy-h/2,w,h);
    }
    if(focus) ctx.globalAlpha=1;
}

function renderGeoOutlines(outlineFn, decorateFn){
    if(!geoOutlineCtx)return;
    const ctx=geoOutlineCtx;
    ctx.clearRect(0,0,canvasW,canvasH);
    const focus=geoFocusPbWO4;

    // 1. batch default outlines (split into two passes when focus mode dims LG)
    ctx.strokeStyle=geoStrokeColor(); ctx.lineWidth=0.5;
    for(const pass of focus?[['W',1],['G',0.15]]:[['*',1]]){
        ctx.globalAlpha=pass[1];
        ctx.beginPath();
        for(let i=0;i<modules.length;i++){
            if(pass[0]!=='*'&&modules[i].t!==pass[0]) continue;
            const style=outlineFn?outlineFn(i):null;
            const hov=hoveredModule&&hoveredModule.n===modules[i].n;
            if(!style&&!hov){
                const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
                ctx.rect(cx-w/2,cy-h/2,w,h);
            }
        }
        ctx.stroke();
    }
    ctx.globalAlpha=1;

    // 2. special outlines (cluster members, warn, selection)
    for(let i=0;i<modules.length;i++){
        const style=outlineFn?outlineFn(i):null;
        if(style){
            if(focus) ctx.globalAlpha=(modules[i].t==='G')?0.15:1;
            const m=modules[i],[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
            ctx.strokeStyle=style.color; ctx.lineWidth=style.width;
            ctx.strokeRect(cx-w/2,cy-h/2,w,h);
        }
    }
    if(focus) ctx.globalAlpha=1;

    // 3. hover highlight (same for all tabs)
    if(hoveredModule){
        const m=hoveredModule,[cx,cy]=d2c(m.x,m.y),w=m.sx*scale,h=m.sy*scale;
        ctx.strokeStyle=THEME.accent; ctx.lineWidth=1.5;
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
