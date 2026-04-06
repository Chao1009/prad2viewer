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

    // blank DQ plots but keep selected module
    Plotly.react('waveform-div',[], wfLayout(selectedModule?selectedModule.n:'', wfWindowNs()), PC2);
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
    filterActive=data.filter_active||false;
    filteredCount=data.filtered_count||totalEvents;
    histEnabled=data.hist_enabled||false;
    histConfig=data.hist||{};
    // populate time cut display
    if(document.getElementById('tcut-min-show'))
        document.getElementById('tcut-min-show').textContent=
            histConfig.time_min!==undefined ? histConfig.time_min : '?';
    if(document.getElementById('tcut-max-show'))
        document.getElementById('tcut-max-show').textContent=
            histConfig.time_max!==undefined ? histConfig.time_max : '?';
    refLines=data.ref_lines||{};
    triggerBitsDef=data.trigger_bits||[];
    triggerTypeDef=data.trigger_type||[];
    // load per-tab trigger filters from server config
    const tf=data.trigger_filter||{};
    for(const [tab, filt] of Object.entries(tf)){
        tabTrigFilter[tab]={accept:filt.trigger_accept||0, reject:filt.trigger_reject||0};
    }
    buildTriggerFilterUI();
    restoreTrigFilterFromTab();
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
    // LIVETIME — temporary
    if(mode==='online') startLivetimePolling(); else stopLivetimePolling();
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
        if(filterActive){
            document.getElementById('ev-total').textContent=`/ ${filteredCount} (${totalEvents} total)`;
            fetch('/api/filter/indices').then(r=>r.json()).then(idx=>{
                filteredIndices=idx;
                if(idx.length>0) loadEvent(idx[0]);
            }).catch(()=>{filteredIndices=null;});
        } else {
            document.getElementById('ev-total').textContent=`/ ${totalEvents}`;
            filteredIndices=null;
        }

        updateHeaderStats();
        if(histEnabled) { fetchOccupancy(); fetchClHist(); }
        fetchEpicsChannels(); fetchEpicsLatest();
        if(activeTab==='epics') fetchAllEpicsSlots();
        if(activeTab==='physics') fetchPhysics();
        syncDqRange();
        geoViewInit=false; resizeGeo();
        if(totalEvents>0 && !filterActive) loadEvent(1);
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
