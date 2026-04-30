// waveform.js — Waveform display, peak visualization, histogram fetching

let currentWaveform=null;  // {x:[], y:[]} for copy button
let wfStackEnabled=false;
let wfStackTraces=[];      // [{x,y},...] accumulated waveforms
let wfStackModKey='';      // module key for current stack (clear on module change)
let wfDaqEnabled=false;    // DAQ-mode toggle: render firmware Mode 1/2/3 annotations
let currentHist={};  // {divId: {x:[], y:[]}} for histogram copy
let lastHistModule = '';
let wfRequestId = 0;  // sequence guard for async waveform fetches

const NS_PER_SAMPLE = 4;  // FADC250: 250 MHz → 4 ns/sample

// Firmware quality bitmask (must match prad2dec/include/Fadc250Data.h).
const Q_GOOD             = 0;
const Q_PEAK_AT_BOUNDARY = 1 << 0;
const Q_NSB_TRUNCATED    = 1 << 1;
const Q_NSA_TRUNCATED    = 1 << 2;
const Q_VA_OUT_OF_RANGE  = 1 << 3;

function qualityLabel(q){
    if(!q) return 'OK';
    const flags=[];
    if(q & Q_PEAK_AT_BOUNDARY) flags.push('peakBnd');
    if(q & Q_NSB_TRUNCATED)    flags.push('NSBtrunc');
    if(q & Q_NSA_TRUNCATED)    flags.push('NSAtrunc');
    if(q & Q_VA_OUT_OF_RANGE)  flags.push('VaOOR');
    return flags.join('+');
}

// Default waveform window from config (in ns). Used for x-range on empty plots.
function wfWindowNs(){
    // ptw (programmable time window) in samples, default 100 samples = 400 ns
    const ptw = histConfig.ptw || 100;
    return ptw * NS_PER_SAMPLE;
}

// Build time-cut shapes for waveform plot (dimmed regions + dashed lines)
function timeCutShapes(xMax){
    const shapes = refShapes('waveform') || [];
    if (isTimeCut()) {
        const dim = {type:'rect', yref:'paper', y0:0, y1:1,
            fillcolor:THEME.overlayLight, line:{width:0}, layer:'above'};
        shapes.push({...dim, xref:'x', x0:0, x1:histConfig.time_min});
        shapes.push({...dim, xref:'x', x0:histConfig.time_max, x1:xMax});
        const edge = {type:'line', yref:'paper', y0:0, y1:1,
            line:{color:THEME.highlight, width:1, dash:'dash'}};
        shapes.push({...edge, x0:histConfig.time_min, x1:histConfig.time_min});
        shapes.push({...edge, x0:histConfig.time_max, x1:histConfig.time_max});
    }
    return shapes;
}

// Waveform plot layout with proper x-range and time cut
function wfLayout(title, xMax){
    return {...PL,
        title:{text:title, font:{size:11,color:THEME.textDim}},
        xaxis:{...PL.xaxis, title:'Time (ns)', range:[0, xMax], autorange:false},
        yaxis:{...PL.yaxis, title:'ADC'},
        shapes:timeCutShapes(xMax),
    };
}

// =========================================================================
// Waveform
// =========================================================================
function showWaveform(mod){
    selectedModule=mod;

    // no waveform data available for this source type
    if(!sourceCaps.has_waveforms){
        document.getElementById('detail-header').innerHTML=
            `<span class="mod-name">${mod.n}</span> <span class="mod-daq">No waveform data for this file type</span>`;
        showHistograms(mod); redrawGeo(); return;
    }

    const key=`${mod.roc}_${mod.sl}_${mod.ch}`;
    const d=eventChannels[key];
    const pedInfo=d?` &nbsp; Ped: ${d.pm.toFixed(1)} ± ${d.pr.toFixed(1)}`:'';
    document.getElementById('detail-header').innerHTML=
        `<span class="mod-name">${mod.n}</span> <span class="mod-daq">${crateName(mod.roc)} &middot; slot ${mod.sl} &middot; ch ${mod.ch}${pedInfo}</span>`;

    // reset stack when switching to a different module
    if(wfStackEnabled && key!==wfStackModKey){
        wfStackTraces=[]; wfStackModKey=key;
        document.getElementById('wf-stack-count').textContent='0/200';
    }

    if(!d){
        if(!wfStackEnabled){
            currentWaveform=null;
            Plotly.react('waveform-div',[], wfLayout(`${mod.n} — No data`, wfWindowNs()), PC2);
            document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No data</td></tr>';
        } else if(wfStackTraces.length===0){
            Plotly.react('waveform-div',[], wfLayout(`${mod.n} — Stacked (0)`, wfWindowNs()), PC2);
        }
        showHistograms(mod); redrawGeo(); return;
    }

    // if samples are already present (e.g. from ring buffer), use them directly;
    // otherwise fetch on demand from /api/waveform/<event>/<key>
    if(d.s){
        renderWaveform(mod, key, d, d.s);
    } else {
        // don't clear the plot while fetching — avoids flash in stacking mode
        const reqId = ++wfRequestId;
        fetch(`/api/waveform/${currentEvent}/${key}`).then(r=>r.json()).then(wf=>{
            if(reqId !== wfRequestId) return;  // stale response, discard
            if(wf.error){ if(!wfStackEnabled) renderWaveform(mod, key, d, null); return; }
            d.s=wf.s;
            if(wf.pk) d.pk=wf.pk;
            if(wf.pm!==undefined) d.pm=wf.pm;
            if(wf.pr!==undefined) d.pr=wf.pr;
            if(wf.daq) d.daq=wf.daq;
            renderWaveform(mod, key, d, d.s);
        }).catch(()=>{ if(reqId === wfRequestId && !wfStackEnabled) renderWaveform(mod, key, d, null); });
    }

    showHistograms(mod); redrawGeo();
}

function renderWaveform(mod, key, d, samples){
    if(!samples){
        if(wfStackEnabled) return;  // skip empty events, keep existing stack
        currentWaveform=null;
        Plotly.react('waveform-div',[], wfLayout(`${mod.n} — No samples`, wfWindowNs()), PC2);
        document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No waveform data</td></tr>';
        document.getElementById('peaks-tbody-daq').innerHTML='';
        return;
    }

    const peaks=d.pk||[], x=samples.map((_,i)=>i*NS_PER_SAMPLE);
    const tMax = (samples.length-1)*NS_PER_SAMPLE;
    currentWaveform={x, y:Array.from(samples)};

    // --- DAQ (firmware-mode) annotations ---
    // Mutually exclusive with stack mode — DAQ annotations don't make sense
    // when overlaying many events.  Stack toggle disables DAQ; DAQ disables stack.
    if(wfDaqEnabled && d.daq){
        renderWaveformDaq(mod, d, samples, x, tMax);
        return;
    }

    // --- stacking mode ---
    if(wfStackEnabled){
        wfStackTraces.push({x:Array.from(x), y:Array.from(samples)});

        const maxStack=200;
        while(wfStackTraces.length>maxStack) wfStackTraces.shift();

        const traces=wfStackTraces.map(w=>({
            x:w.x, y:w.y, type:'scatter', mode:'lines',
            line:{color:THEME.accent, width:1, opacity:0.25},
            showlegend:false, hoverinfo:'skip',
        }));
        if(wfStackTraces.length>0){
            const last=wfStackTraces[wfStackTraces.length-1];
            traces.push({x:last.x, y:last.y, type:'scatter', mode:'lines',
                name:'Latest', line:{color:THEME.accent, width:1.5}, showlegend:false});
        }

        let ylo=Infinity, yhi=-Infinity;
        for(const w of wfStackTraces) for(const v of w.y){ if(v<ylo) ylo=v; if(v>yhi) yhi=v; }
        const pad=(yhi-ylo)*0.05||5;

        document.getElementById('wf-stack-count').textContent=`${wfStackTraces.length}/${maxStack}`;
        const stackLayout = wfLayout(`${mod.n} — Stacked (${wfStackTraces.length})`, tMax);
        stackLayout.yaxis = {...stackLayout.yaxis, range:[ylo-pad,yhi+pad], autorange:false};
        Plotly.react('waveform-div', traces, stackLayout, PC2);

        document.getElementById('peaks-tbody').innerHTML=
            '<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">Stack mode — peaks hidden</td></tr>';
        return;
    }

    // --- normal (single event) mode ---
    const traces=[
        {x,y:samples,type:'scatter',mode:'lines',name:'Waveform',line:{color:THEME.accent,width:1}},
        {x:[0,tMax],y:[d.pm,d.pm],type:'scatter',mode:'lines',name:'Pedestal',line:{color:THEME.textMuted,width:1,dash:'dash'}},
    ];
    const thr=d.pm+Math.max(5*d.pr,3);
    traces.push({x:[0,tMax],y:[thr,thr],type:'scatter',mode:'lines',line:{color:THEME.grid,width:1,dash:'dot'},showlegend:false});
    peaks.forEach((p,i)=>{
        const col=PC[i%PC.length],px=[],py=[];
        for(let j=p.l;j<=p.r;j++){px.push(j*NS_PER_SAMPLE);py.push(samples[j]);}
        const r=parseInt(col.slice(1,3),16),g=parseInt(col.slice(3,5),16),b=parseInt(col.slice(5,7),16);
        const fill=`rgba(${r},${g},${b},0.18)`;
        traces.push({x:px,y:px.map(()=>d.pm),type:'scatter',mode:'lines',
            line:{width:0},showlegend:false,hoverinfo:'skip'});
        traces.push({x:px,y:py,type:'scatter',mode:'lines',name:`Peak ${i}`,
            line:{color:col,width:2},fill:'tonexty',fillcolor:fill});
        traces.push({x:[p.p*NS_PER_SAMPLE],y:[samples[p.p]],type:'scatter',mode:'markers',
            marker:{color:col,size:7,symbol:'diamond'},showlegend:false});
    });

    const layout = wfLayout(`${mod.n} — Event ${currentEvent}`, tMax);
    layout.legend = {x:1,y:1,xanchor:'right',bgcolor:THEME.overlay,font:{size:9}};
    Plotly.react('waveform-div', traces, layout, PC2);

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
// DAQ-mode (firmware Mode 1/2/3) renderer
//
// Annotates the waveform with the firmware-emulated TDC + windowing per the
// FADC250 User's Manual: TET line, Vnoise baseline, NSB/NSA brackets around
// each Tcross, Vp marker, vertical T line at the interpolated mid-amplitude
// time, and the Mode 2 integration polygon (Σ).
// =========================================================================
function renderWaveformDaq(mod, d, samples, x, tMax){
    const daq = d.daq || {};
    const pulses = daq.pk || [];
    const pedUsed = (daq.ped_used !== undefined) ? daq.ped_used : (d.pm || 0);
    const tet     = daq.tet || 0;
    const nsb     = daq.nsb || 0;
    const nsa     = daq.nsa || 0;

    // Update the inline read-out next to the toggle.
    const info = document.getElementById('wf-daq-info');
    if(info){
        info.textContent = `TET=${tet} · NSB=${nsb} · NSA=${nsa} · PED=${pedUsed.toFixed(1)}`;
    }
    // Pedestal-subtracted samples — what the firmware sees and what all
    // annotations are positioned against.
    const ys = samples.map(v => Math.max(0, v - pedUsed));

    const traces = [
        {x, y:ys, type:'scatter', mode:'lines', name:'Waveform (ped-sub)',
         line:{color:THEME.accent, width:1}},
        {x:[0,tMax], y:[0,0], type:'scatter', mode:'lines', name:'Vnoise',
         line:{color:THEME.textMuted, width:1, dash:'dash'}},
        {x:[0,tMax], y:[tet,tet], type:'scatter', mode:'lines', name:`TET = ${tet}`,
         line:{color:'#ff6b6b', width:1.4, dash:'dash'}},
    ];

    const shapes = timeCutShapes(tMax);
    const annotations = [];

    pulses.forEach((p, idx) => {
        const col = PC[idx % PC.length];
        const r = parseInt(col.slice(1,3),16),
              g = parseInt(col.slice(3,5),16),
              b = parseInt(col.slice(5,7),16);
        const fill = `rgba(${r},${g},${b},0.20)`;
        const tCrossNs = p.cross * NS_PER_SAMPLE;

        // Mode-2 integration polygon — fill under the curve over [wlo, whi].
        const wx=[], wy=[];
        for(let j=p.wlo; j<=p.whi; ++j){
            wx.push(j*NS_PER_SAMPLE);
            wy.push(ys[j]);
        }
        if(wx.length){
            traces.push({x:wx, y:wx.map(()=>0), type:'scatter', mode:'lines',
                line:{width:0}, showlegend:false, hoverinfo:'skip'});
            traces.push({x:wx, y:wy, type:'scatter', mode:'lines',
                name:`Σ${idx} = ${p.i.toFixed(0)}`,
                line:{color:col, width:2}, fill:'tonexty', fillcolor:fill});
        }

        // Vp marker (open circle at the peak sample, not Tcross).
        const peakSample = (p.vp_pos !== undefined) ? p.vp_pos : p.cross;
        const tPeakNs = peakSample * NS_PER_SAMPLE;
        traces.push({x:[tPeakNs], y:[p.vp], type:'scatter',
            mode:'markers', marker:{color:col, size:9, symbol:'circle-open',
            line:{color:col, width:2}}, showlegend:false, hoverinfo:'skip'});
        annotations.push({
            x:tPeakNs, y:p.vp, ax:30, ay:-20,
            text:`Vp${idx}=${p.vp.toFixed(0)}`,
            font:{color:col, size:10}, arrowcolor:col,
            showarrow:true, arrowhead:2,
        });

        // Vertical T line at the interpolated mid-amplitude time.
        shapes.push({type:'line', x0:p.t, x1:p.t, yref:'paper', y0:0, y1:1,
            line:{color:col, width:1.6}});
        annotations.push({
            x:p.t, y:1, yref:'paper', yanchor:'bottom',
            text:`T${idx}=${p.t.toFixed(2)} ns`,
            font:{color:col, size:10}, showarrow:false,
        });

        // NSB/NSA brackets (just below the peak, light shaded blocks).
        const yBracket = Math.max(0.05*p.vp, 5);
        const tNsbLo = (p.cross - nsb) * NS_PER_SAMPLE;
        const tNsaHi = (p.cross + nsa) * NS_PER_SAMPLE;
        // NSB span (orange)
        shapes.push({type:'line', x0:tNsbLo, x1:tCrossNs,
            y0:yBracket, y1:yBracket,
            line:{color:'#ffa94d', width:1.5}});
        // NSA span (cyan)
        shapes.push({type:'line', x0:tCrossNs, x1:tNsaHi,
            y0:yBracket, y1:yBracket,
            line:{color:'#22b8cf', width:1.5}});
        // Labels in ns, floored to the nearest sample (NSB and NSA are
        // integer-sample firmware register values; one sample = 4 ns at
        // 250 MHz).
        const nsbNs = Math.floor(nsb) * NS_PER_SAMPLE;
        const nsaNs = Math.floor(nsa) * NS_PER_SAMPLE;
        annotations.push({x:(tNsbLo+tCrossNs)/2, y:yBracket,
            yanchor:'bottom', yshift:2,
            text:`NSB=${nsbNs} ns`, font:{color:'#ffa94d', size:9},
            showarrow:false});
        annotations.push({x:(tCrossNs+tNsaHi)/2, y:yBracket,
            yanchor:'bottom', yshift:2,
            text:`NSA=${nsaNs} ns`, font:{color:'#22b8cf', size:9},
            showarrow:false});

        // Tcross dotted vertical (subtle).
        shapes.push({type:'line', x0:tCrossNs, x1:tCrossNs, yref:'paper',
            y0:0, y1:1, line:{color:'#ff6b6b', width:0.8, dash:'dot'}});
    });

    const title = pulses.length
        ? `${mod.n} — DAQ mode · ${pulses.length} pulse${pulses.length>1?'s':''}`
        : `${mod.n} — DAQ mode · 0 pulses (Vp ≤ TET)`;
    const layout = wfLayout(title, tMax);
    layout.shapes = shapes;
    layout.annotations = annotations;
    layout.yaxis = {...layout.yaxis, title:'ADC (ped-subtracted)'};
    layout.legend = {x:1, y:1, xanchor:'right', bgcolor:THEME.overlay,
                     font:{size:9}};
    Plotly.react('waveform-div', traces, layout, PC2);

    // DAQ peaks table.
    let rows='';
    pulses.forEach((p, i) => {
        const col = PC[i % PC.length];
        rows += `<tr style="border-left:3px solid ${col}">`
              + `<td>${p.n}</td><td>${p.cross}</td><td>${p.t.toFixed(2)}</td>`
              + `<td>${p.vp.toFixed(0)}</td><td>${p.i.toFixed(0)}</td>`
              + `<td>${p.wlo}..${p.cross}</td><td>${p.cross}..${p.whi}</td>`
              + `<td style="text-align:center" title="${qualityLabel(p.q)}">${p.q?qualityLabel(p.q):'OK'}</td>`
              + `</tr>`;
    });
    if(!pulses.length) rows = '<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No firmware pulses (Vp ≤ TET)</td></tr>';
    document.getElementById('peaks-tbody-daq').innerHTML = rows;
}

// =========================================================================
// Histograms
// =========================================================================
function fetchAndPlotHist(divId, url, title, xTitle, binMin, binStep, barColor, logYId, refKey, timeCut){
    fetch(url).then(r=>r.json()).then(data=>{
        if(data.error||!data.bins||!data.bins.length){
            currentHist[divId]=null;
            Plotly.react(divId,[],{...PL,title:{text:`${title} — No data`,font:{size:10,color:THEME.textMuted}}},PC2);
            return;
        }
        const x=data.bins.map((_,i)=>binMin+(i+0.5)*binStep);
        const cx=[], cy=[];
        for(let i=0;i<data.bins.length;i++){if(data.bins[i]>0){cx.push(x[i]);cy.push(data.bins[i]);}}
        currentHist[divId]={x:cx,y:cy};

        const entries=data.bins.reduce((a,b)=>a+b,0)+data.underflow+data.overflow;
        const stats=`${data.events} evts | Entries: ${entries}  Under: ${data.underflow}  Over: ${data.overflow}`;
        const xMin=binMin, xMax=binMin+data.bins.length*binStep;
        const shapes = refKey ? (refShapes(refKey)||[]) : [];
        if (timeCut && isTimeCut()) {
            const dim={type:'rect',yref:'paper',y0:0,y1:1,
                fillcolor:THEME.overlayLight,line:{width:0},layer:'below'};
            shapes.push({...dim, x0:xMin, x1:histConfig.time_min});
            shapes.push({...dim, x0:histConfig.time_max, x1:xMax});
            const edge={type:'line',yref:'paper',y0:0,y1:1,
                line:{color:THEME.highlight, width:1, dash:'dash'}};
            shapes.push({...edge, x0:histConfig.time_min, x1:histConfig.time_min});
            shapes.push({...edge, x0:histConfig.time_max, x1:histConfig.time_max});
        }
        Plotly.react(divId,[{
            x,y:data.bins,type:'bar',marker:{color:barColor,line:{width:0}},
            hovertemplate:'%{x:.0f}: %{y}<extra></extra>',
        }],{...PL,
            title:{text:`${title}<br><span style="font-size:9px;color:var(--theme-text-dim)">${stats}</span>`,font:{size:10,color:THEME.textDim}},
            xaxis:{...PL.xaxis,title:xTitle,range:[xMin,xMax]},
            yaxis:{...PL.yaxis,title:'Counts',
                type:logYId&&document.getElementById(logYId).checked?'log':'linear'},
            bargap:0.05,
            shapes,
        },PC2);
    }).catch(()=>{
        currentHist[divId]=null;
        Plotly.react(divId,[],{...PL,title:{text:'Fetch error',font:{size:10,color:THEME.danger}}},PC2);
    });
}

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
    fetchAndPlotHist('heighthist-div',`/api/heighthist/${key}`,
        `${mod.n} Peak Height [${h.time_min||170}-${h.time_max||190} ns]`,
        'Peak Height', h.height_min||0, h.height_step||10, '#e599f7', 'heighthist-logy', 'height_hist');
    fetchAndPlotHist('inthist-div',`/api/hist/${key}`,
        `${mod.n} Integral [${h.time_min||170}-${h.time_max||190} ns]`,
        'Peak Integral', h.bin_min||0, h.bin_step||100, '#00b4d8', 'inthist-logy', 'integral_hist');
    fetchAndPlotHist('poshist-div',`/api/poshist/${key}`,
        `${mod.n} Peak Position`,
        'Time (ns)', h.pos_min||0, h.pos_step||4, '#51cf66', null, 'time_hist', true);
}
