// epics.js — EPICS slow control monitoring tab
//
// Depends on globals from viewer.js: PL, PC_EPICS, activeTab

const EPICS_COLORS=['#00b4d8','#ff6b6b','#51cf66','#ffd43b','#cc5de8','#ff922b'];
const EPICS_MAX_PER_SLOT=6;
const EPICS_NUM_SLOTS=6;
let epicsChannels=[];
let epicsSlots=[[],[],[],[],[],[]];
let epicsWarnThresh=0.1, epicsAlertThresh=0.2, epicsMinAvgPts=10;
let epicsLatestData=null;
let epicsSlotData=new Array(EPICS_NUM_SLOTS).fill(null); // cached fetch results per slot
let lastEpicsFetch=0, refreshEpicsMs=2000;

// =========================================================================
// Data fetching
// =========================================================================

function fetchEpicsChannels(){
    fetch('/api/epics/channels').then(r=>r.json()).then(data=>{
        epicsChannels=data.channels||[];
    }).catch(()=>{});
}

function fetchAndPlotEpicsSlot(slot){
    const names=epicsSlots[slot];
    if(!names.length){
        epicsSlotData[slot]=null;
        Plotly.react('epics-plot-'+slot,[],{...PL},PC_EPICS);
        return;
    }
    // batch fetch: single request for all channels in this slot
    const query=names.map(n=>'ch='+encodeURIComponent(n)).join('&');
    fetch(`/api/epics/batch?${query}`).then(r=>r.json()).then(batch=>{
        // reshape batch response to match the per-channel format
        const results=(batch.channels||[]).map(ch=>({
            name:ch.name, time:batch.time||[], value:ch.value||[], count:ch.count||0
        }));
        epicsSlotData[slot]=results;
        const traces=[];
        results.forEach((data,i)=>{
            if(!data||!data.time||!data.time.length) return;
            traces.push({
                x:data.time,y:data.value,type:'scatter',mode:'lines+markers',
                name:data.name,
                line:{color:EPICS_COLORS[i%EPICS_COLORS.length],width:1.5},
                marker:{size:3,color:EPICS_COLORS[i%EPICS_COLORS.length]},
                hovertemplate:`${data.name}: %{y:.3f} (%{x:.1f}s)<extra></extra>`,
            });
        });
        Plotly.react('epics-plot-'+slot,traces,{...PL,
            xaxis:{...PL.xaxis,title:'Time (s)'},
            yaxis:{...PL.yaxis,autorange:true,automargin:true,
                rangemode:'normal',constraintoward:'center',
                range: (()=>{
                    // add 10% padding to y range
                    let ymin=Infinity,ymax=-Infinity;
                    for(const t of traces) for(const v of t.y){if(v<ymin)ymin=v;if(v>ymax)ymax=v;}
                    if(!isFinite(ymin)) return undefined;
                    const pad=Math.max(Math.abs(ymax-ymin)*0.1,Math.abs(ymax)*0.01,1e-6);
                    return [ymin-pad,ymax+pad];
                })()},
            showlegend:traces.length>1,
            legend:{font:{size:9,color:'#aaa'},bgcolor:'rgba(0,0,0,0)',x:0,y:1},
        },PC_EPICS);
    });
}

function fetchAllEpicsSlots(){
    for(let i=0;i<EPICS_NUM_SLOTS;i++) fetchAndPlotEpicsSlot(i);
}

function fetchEpicsLatest(){
    fetch('/api/epics/latest').then(r=>r.json()).then(data=>{
        epicsLatestData=data;
        updateEpicsTable();
        if(activeTab!=='epics') updateEpicsDot();
    }).catch(()=>{});
}

function updateEpicsDot(){
    const dot=document.getElementById('epics-dot');
    if(!epicsLatestData||!epicsLatestData.channels){dot.className='tab-dot';return;}
    let worst=0; // 0=ok, 1=warn, 2=alert
    for(const ch of epicsLatestData.channels){
        if(ch.count<epicsMinAvgPts||ch.mean===0) continue;
        const dev=Math.abs(ch.value-ch.mean)/Math.abs(ch.mean);
        if(dev>=epicsAlertThresh) worst=2;
        else if(dev>=epicsWarnThresh && worst<1) worst=1;
        if(worst===2) break;
    }
    dot.className='tab-dot'+(worst===2?' alert':worst===1?' warn':'');
}

// =========================================================================
// Slot management
// =========================================================================

function addEpicsChannel(slot,name){
    if(epicsSlots[slot].includes(name)) return;
    if(epicsSlots[slot].length>=EPICS_MAX_PER_SLOT) return;
    epicsSlots[slot].push(name);
    renderEpicsChips(slot);
    fetchAndPlotEpicsSlot(slot);
}

function removeEpicsChannel(slot,name){
    epicsSlots[slot]=epicsSlots[slot].filter(n=>n!==name);
    renderEpicsChips(slot);
    fetchAndPlotEpicsSlot(slot);
}

function renderEpicsChips(slot){
    const container=document.getElementById('epics-chips-'+slot);
    container.innerHTML=epicsSlots[slot].map((name,i)=>
        `<span class="epics-chip" style="background:${EPICS_COLORS[i%EPICS_COLORS.length]}33;color:${EPICS_COLORS[i%EPICS_COLORS.length]}">`+
        `${name}<span class="chip-x" data-slot="${slot}" data-name="${name}">&times;</span></span>`
    ).join('');
    container.querySelectorAll('.chip-x').forEach(x=>{
        x.onclick=()=>removeEpicsChannel(parseInt(x.dataset.slot),x.dataset.name);
    });
}

// =========================================================================
// Summary table
// =========================================================================

function updateEpicsTable(){
    const tbody=document.getElementById('epics-tbody');
    if(!epicsLatestData||!epicsLatestData.channels||!epicsLatestData.channels.length){
        tbody.innerHTML='<tr><td colspan="5" style="text-align:center;color:var(--dim)">Waiting for EPICS data...</td></tr>';
        return;
    }
    let html='';
    for(const ch of epicsLatestData.channels){
        let cls='', statusText='OK';
        if(ch.count>=epicsMinAvgPts && ch.mean!==0){
            const dev=Math.abs(ch.value-ch.mean)/Math.abs(ch.mean);
            if(dev>=epicsAlertThresh){ cls='epics-alert'; statusText='JUMPING'; }
            else if(dev>=epicsWarnThresh){ cls='epics-warn'; statusText='CHANGING'; }
        }else if(ch.count<epicsMinAvgPts){
            statusText='--';
        }
        const fmtVal=typeof ch.value==='number'?ch.value.toFixed(3):ch.value;
        const fmtMean=typeof ch.mean==='number'?ch.mean.toFixed(3):ch.mean;
        html+=`<tr class="epics-table-row" data-channel="${ch.name}" draggable="true">`;
        html+=`<td style="text-align:left">${ch.name}</td>`;
        html+=`<td class="${cls}">${fmtVal}</td>`;
        html+=`<td>${fmtMean}</td>`;
        html+=`<td>${ch.count}</td>`;
        html+=`<td class="${cls}">${statusText}</td></tr>`;
    }
    tbody.innerHTML=html;
    tbody.querySelectorAll('.epics-table-row').forEach(row=>{
        row.ondragstart=(e)=>{
            e.dataTransfer.setData('text/plain',row.dataset.channel);
            e.dataTransfer.effectAllowed='copy';
        };
        row.onclick=()=>{
            const name=row.dataset.channel;
            for(let s=0;s<EPICS_NUM_SLOTS;s++){
                if(epicsSlots[s].length<EPICS_MAX_PER_SLOT && !epicsSlots[s].includes(name)){
                    addEpicsChannel(s,name); return;
                }
            }
        };
    });
}

// =========================================================================
// Clear
// =========================================================================

function clearEpicsFrontend(){
    epicsLatestData=null;
    epicsSlots=[[],[],[],[],[],[]];
    for(let s=0;s<EPICS_NUM_SLOTS;s++){
        renderEpicsChips(s);
        try{Plotly.react('epics-plot-'+s,[],{...PL},PC_EPICS);}catch(e){}
    }
    document.getElementById('epics-tbody').innerHTML='';
}

// =========================================================================
// Init — called from viewer.js init() with config data
// =========================================================================

function initEpics(data){
    // config
    if(data&&data.epics){
        if(data.epics.warn_threshold!==undefined) epicsWarnThresh=data.epics.warn_threshold;
        if(data.epics.alert_threshold!==undefined) epicsAlertThresh=data.epics.alert_threshold;
        if(data.epics.min_avg_points!==undefined) epicsMinAvgPts=data.epics.min_avg_points;
        const cfgSlots=data.epics.slots||[];
        for(let s=0;s<Math.min(EPICS_NUM_SLOTS,cfgSlots.length);s++){
            const names=cfgSlots[s]||[];
            for(const n of names){
                if(epicsSlots[s].length<EPICS_MAX_PER_SLOT)
                    epicsSlots[s].push(n);
            }
        }
        for(let s=0;s<EPICS_NUM_SLOTS;s++) renderEpicsChips(s);
    }

    // copy buttons
    document.querySelectorAll('.epics-copy').forEach(btn=>{
        const slot=parseInt(btn.closest('.epics-slot').dataset.slot);
        btn.onclick=()=>{
            const results=epicsSlotData[slot];
            if(!results) return;
            let text='';
            for(const data of results){
                if(!data||!data.time||!data.time.length) continue;
                text+=`# ${data.name}\n`;
                text+=`time: [${data.time.join(', ')}]\n`;
                text+=`value: [${data.value.join(', ')}]\n\n`;
            }
            if(!text) return;
            navigator.clipboard.writeText(text).then(()=>{
                btn.textContent='\u2713'; setTimeout(()=>{btn.textContent='copy';},1000);
            });
        };
    });

    // search-as-you-type
    document.querySelectorAll('.epics-search').forEach(input=>{
        const slot=parseInt(input.dataset.slot);
        const dropdown=input.parentElement.querySelector('.epics-dropdown');
        input.oninput=()=>{
            const q=input.value.toLowerCase();
            if(!q){dropdown.classList.remove('open');return;}
            const matches=epicsChannels.filter(n=>n.toLowerCase().includes(q)).slice(0,20);
            if(!matches.length){dropdown.classList.remove('open');return;}
            dropdown.innerHTML=matches.map(n=>
                `<div class="epics-dropdown-item" data-name="${n}">${n}</div>`
            ).join('');
            dropdown.classList.add('open');
            dropdown.querySelectorAll('.epics-dropdown-item').forEach(item=>{
                item.onclick=()=>{
                    addEpicsChannel(slot,item.dataset.name);
                    input.value='';
                    dropdown.classList.remove('open');
                };
            });
        };
        input.onblur=()=>setTimeout(()=>dropdown.classList.remove('open'),200);
    });

    // drag-and-drop
    document.querySelectorAll('.epics-slot').forEach(slotEl=>{
        const slot=parseInt(slotEl.dataset.slot);
        slotEl.ondragover=(e)=>{e.preventDefault();e.dataTransfer.dropEffect='copy';slotEl.classList.add('drag-over');};
        slotEl.ondragleave=()=>slotEl.classList.remove('drag-over');
        slotEl.ondrop=(e)=>{
            e.preventDefault();
            slotEl.classList.remove('drag-over');
            const name=e.dataTransfer.getData('text/plain');
            if(name) addEpicsChannel(slot,name);
        };
    });
}
