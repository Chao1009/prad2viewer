// lms.js — LMS gain monitoring tab: geo provider, history, summary table

let g_lmsWarnThresh=0.1;
let g_lmsRefIndex=-1;  // -1 = None (no normalization)
let currentLmsData=null;  // {x:[], y:[]} for copy button
let lmsSummaryData=null;  // {modules:{idx:{name,mean,rms,count,warn}}, events}
let lmsSelectedModule=-1;

function geoLms(){
    const metric=document.getElementById('lms-color-metric').value;
    const useLog=document.getElementById('lms-log-scale').checked;
    const mods=lmsSummaryData?lmsSummaryData.modules:{};

    const lmsVal=md=>{
        if(!md) return null;
        if(metric==='mean')     return md.mean;
        if(metric==='rms')      return md.rms;
        if(metric==='rms_frac') return md.mean>0?md.rms/md.mean:0;
        return md.warn?1:0;   // 'warn'
    };

    let autoMax=0;
    for(const k in mods){ const v=lmsVal(mods[k]); if(v>autoMax) autoMax=v; }
    if(autoMax<=0) autoMax=1;
    const lmsr=getGeoRange('lms',metric);
    const vmin=lmsr[0]!==null?lmsr[0]:0;
    const vmax=lmsr[1]!==null?lmsr[1]:autoMax;
    const dp=metric==='rms_frac'?3:metric==='rms'?2:0;
    document.getElementById('lms-range-min-show').textContent=vmin.toFixed(dp);
    document.getElementById('lms-range-max-show').textContent=vmax.toFixed(dp);

    renderGeo(
        i => {
            const md=mods[String(i)];
            const val=lmsVal(md);
            if(val!==null&&val>0){
                if(metric==='warn') return md.warn?THEME.danger:THEME.success;
                return geoValueColor(val,vmin,vmax,useLog);
            }
            return geoEmptyColor(modules[i].t);
        },
        i => {
            if(lmsSelectedModule===i) return {color:THEME.selectBorder,width:2.5};
            const md=mods[String(i)];
            if(md&&md.warn) return {color:THEME.danger,width:1.5};
            return null;
        },
        null
    );
}

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
                title:{text:`${modName} — No LMS data`,font:{size:10,color:THEME.textMuted}}},PC2);
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
             type:'scatter', mode:'lines', line:{color:THEME.success,width:1,dash:'dash'}, name:`Mean ${mean.toFixed(0)}`},
            {x:tRange, y:[warnHi,warnHi],
             type:'scatter', mode:'lines', line:{color:THEME.danger,width:1,dash:'dot'}, showlegend:false},
            {x:tRange, y:[warnLo,warnLo],
             type:'scatter', mode:'lines', line:{color:THEME.danger,width:1,dash:'dot'}, showlegend:false},
        ],{...PL,
            title:{text:`LMS — ${modName} (${data.events} pts)${g_lmsRefIndex>=0?' [ref corrected]':''}`,
                font:{size:10,color:THEME.text}},
            xaxis:{...PL.xaxis,
                title:data.sync_unix
                    ?`Time (s) after ${new Date((data.sync_unix - data.sync_rel_sec)*1000).toISOString().replace('T',' ').slice(0,19)} UTC`
                    :'Time (s)'},
            yaxis:{...PL.yaxis,title:g_lmsRefIndex>=0?'Corrected Integral':'Integral'},
            legend:{x:1,y:1,xanchor:'right',bgcolor:THEME.overlay,font:{size:9}},
            margin:{...PL.margin,t:28,b:36},
            shapes:refShapes('lms'),
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
