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

// LIVETIME — temporary: poll server for DAQ livetime
const LIVETIME_POLL_MS=5000;
let livetimeTimer=null;
function pollLivetime(){
    fetch('/api/livetime').then(r=>r.json()).then(d=>{
        const el=document.getElementById('livetime-display');
        if(!el) return;
        el.style.display='';
        if(d.livetime>=0){
            el.textContent='Livetime: '+d.livetime.toFixed(1)+'%';
            el.style.color='#51cf66';
        } else {
            el.textContent='Livetime: N/A';
            el.style.color='#888';
        }
    }).catch(()=>{});
}
function startLivetimePolling(){
    if(livetimeTimer) return;
    pollLivetime();
    livetimeTimer=setInterval(pollLivetime,LIVETIME_POLL_MS);
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

