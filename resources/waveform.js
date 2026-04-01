// waveform.js — Waveform display, peak visualization, histogram fetching

let currentWaveform=null;  // {x:[], y:[]} for copy button
let wfStackEnabled=false;
let wfStackTraces=[];      // [{x,y},...] accumulated waveforms
let wfStackModKey='';      // module key for current stack (clear on module change)
const wfYminSteps=[0,150], wfYmaxSteps=[500,1000,2000,4000];
let wfYminIdx=1, wfYmaxIdx=0;
let wfYmin=wfYminSteps[1], wfYmax=wfYmaxSteps[0];
let currentHist={};  // {divId: {x:[], y:[]}} for histogram copy
let lastHistModule = '';
let wfRequestId = 0;  // sequence guard for async waveform fetches

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

    if(!d){
        if(!wfStackEnabled){
            currentWaveform=null;
            Plotly.react('waveform-div',[],{...PL,title:{text:`${mod.n} — No data`,font:{size:11,color:'#555'}}},PC2);
            document.getElementById('peaks-tbody').innerHTML='<tr><td colspan="8" style="text-align:center;color:var(--dim);padding:8px">No data</td></tr>';
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
            renderWaveform(mod, key, d, d.s);
        }).catch(()=>{ if(reqId === wfRequestId && !wfStackEnabled) renderWaveform(mod, key, d, null); });
    }

    showHistograms(mod); redrawGeo();
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
            xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC',range:[wfYmin,wfYmax]},
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
        xaxis:{...PL.xaxis,title:'Sample'},yaxis:{...PL.yaxis,title:'ADC',range:[wfYmin,wfYmax]},
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
