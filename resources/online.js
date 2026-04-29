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
        el.style.color = THEME.success;
    } else if (waiting) {
        el.textContent = `● Waiting for ET (${retries||'...'})`;
        el.style.color = THEME.warn;
    } else {
        el.textContent = '● Disconnected';
        el.style.color = THEME.danger;
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

// LIVETIME — poll server for DAQ livetime (server shells out to caget).
// Thresholds + poll interval come from the server config in applyConfig().
// livetimeEnabled is set from cfg.livetime.enabled (server-side command !=  "").
let livetimePollMs=5000;
let livetimeHealthy=90, livetimeWarning=80;
let livetimeEnabled=false;
let livetimeTimer=null;
function pollLivetime(){
    fetch('/api/livetime').then(r=>r.json()).then(d=>{
        const el=document.getElementById('livetime-display');
        if(!el) return;
        el.style.display='';
        const ts=(d.livetime>=0)?d.livetime:null;
        const meas=(d.measured>=0)?d.measured:null;
        if(ts==null && meas==null){
            el.textContent='DAQ Livetime: N/A';
            el.style.color=THEME.textDim;
            return;
        }
        const parts=[];
        if(ts!=null)   parts.push(ts.toFixed(1)+'% (TS)');
        if(meas!=null) parts.push(meas.toFixed(1)+'% (Meas.)');
        el.textContent='DAQ Livetime: '+parts.join(' / ');
        // Color by the worse of the two so a sick channel still flags red.
        const worst=Math.min(ts??meas, meas??ts);
        el.style.color=worst>=livetimeHealthy?THEME.success
                      :worst>=livetimeWarning?THEME.warn:THEME.danger;
    }).catch(()=>{});
}
function startLivetimePolling(){
    if(livetimeTimer) return;
    if(!livetimeEnabled) return;          // server has no command configured
    pollLivetime();
    livetimeTimer=setInterval(pollLivetime,livetimePollMs);
}
function stopLivetimePolling(){
    if(livetimeTimer){clearInterval(livetimeTimer);livetimeTimer=null;}
    const el=document.getElementById('livetime-display');
    if(el) el.style.display='none';
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
                    if(activeTab==='cluster') fetchGemResiduals();
                }
                // gem_apv tab is per-event; loadEventData (called by
                // loadLatestEvent above) hooks into activeTab and refetches
                // once currentEvent has been updated.
            } else if (msg.type === 'status') {
                setEtStatus(msg.connected, msg.waiting, msg.retries);
            } else if (msg.type === 'hist_cleared') {
                occData={}; occTcutData={}; occTotal=0;
                initClHist(); plotClHist(); plotClStatHists();
                gemResidData=null; plotGemResiduals();
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
            } else if (msg.type === 'gem_threshold_updated') {
                // Another viewer (or this one) changed the GEM σ — keep
                // every open tab's input in sync immediately so users
                // see the same value across windows.  hits[] update
                // arrives naturally with the next event.  Guarded so
                // online.js loading before gem_apv.js can't TDZ-throw.
                if (typeof syncGemApvZsSigmaInput === 'function')
                    syncGemApvZsSigmaInput(msg.zs_sigma);
                if (typeof gemApvCalib !== 'undefined' && gemApvCalib)
                    gemApvCalib.zs_sigma = msg.zs_sigma;
            }
        } catch (e) {}
    };
}

