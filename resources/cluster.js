// cluster.js — Clustering tab: geo provider, cluster data, energy histograms

// --- clustering tab state ---
let clusterData=null;  // {hits:{}, clusters:[]}
let selectedCluster=-1;  // -1 = all
let clusterEvent=-1;  // event number for cached cluster data

// per-event GEM hits (lab xy + projected to HyCal-local) for the geo overlay
let gemHits=null;     // {detectors:[{id,name,hits_2d:[{x,y,proj_x,proj_y,...}]}]}
let gemHitsEvent=-1;
// accumulated GEM↔HyCal residuals for the 4 small panels
let gemResidData=null;
// per-detector palette for the geo overlay
const GEM_DOT_COLORS=['#ff6b6b','#51cf66','#00b4d8','#ffa500'];

// cluster energy histogram (accumulated on frontend)
let clHistBins=null, clHistEvents=0;
let clHistMin=0, clHistMax=3000, clHistStep=10;
let currentClHist=null;  // {x:[], y:[]} for copy button
let currentNclustHist=null, currentNblocksHist=null;
// per-event sum of all module energies — sibling of the cluster energy hist.
let rawEnergyBins=null;
let rawEnergyMin=0, rawEnergyMax=6000, rawEnergyStep=20;
let currentRawEnergyHist=null;

// cluster count histograms (configurable via monitor_config.json hycal_hist section)
let nclustBins=null, nblocksBins=null;
let nclustMin=0.5, nclustMax=10.5, nclustStep=1;
let nblocksMin=0, nblocksMax=40, nblocksStep=1;
// Per-Ncl bucket arrays (parallel to nclustBins): when the user clicks
// a bar in cl-nclust-hist, we redraw cl-energy-hist and cl-nblocks-hist
// from the corresponding bucket instead of the unfiltered histograms.
// `selectedNcl` is the bucket index, or -1 for "show unfiltered".
let clEnergyBinsByNcl=null, nblocksBinsByNcl=null;
let selectedNcl=-1;

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
            if(selectedModule&&selectedModule.n===modules[i].n) return {color:THEME.selectBorder,width:2.5};
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
            // GEM 2D hits projected onto HyCal local plane.  Draw all hits;
            // any that fall outside the canvas just clip naturally.
            if(gemHits && gemHits.detectors){
                const ringColor=THEME.text;
                gemHits.detectors.forEach(det=>{
                    const fill=GEM_DOT_COLORS[det.id % GEM_DOT_COLORS.length];
                    (det.hits_2d||[]).forEach(h=>{
                        if(h.proj_x==null||h.proj_y==null) return;
                        const [px,py]=d2c(h.proj_x,h.proj_y);
                        ctx.fillStyle=fill;
                        ctx.beginPath();ctx.arc(px,py,4,0,2*Math.PI);ctx.fill();
                        ctx.strokeStyle=ringColor; ctx.lineWidth=1;
                        ctx.stroke();
                    });
                });
            }
        }
    );
}

function loadGemHits(evnum){
    if(gemHitsEvent===evnum && gemHits) { geoCluster(); return; }
    fetch('/api/gem/hits').then(r=>r.json()).then(data=>{
        gemHits=data;
        gemHitsEvent=evnum;
        geoCluster();
    }).catch(()=>{ gemHits=null; });
}

// build a set of module indices belonging to a cluster
function clusterModuleSet(clIdx){
    if(!clusterData||clIdx<0) return null;
    const cl=clusterData.clusters[clIdx];
    if(!cl) return null;
    return new Set(cl.modules);
}

function loadClusterData(evnum){
    if(clusterEvent===evnum && clusterData) {
        geoCluster(); updateClusterTable();
        loadGemHits(evnum);
        return;
    }
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
        // refresh cluster histograms from server (accumulated there)
        fetchClHist();
        fetchGemResiduals();
        loadGemHits(evnum);
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
    const hits=clusterData.hits||{};
    const moduleE=Object.values(hits).reduce((s,e)=>s+e,0);
    const clusterE=clusters.reduce((s,c)=>s+c.energy,0);
    const hi='color:var(--accent);font-weight:700';
    document.getElementById('cl-summary').innerHTML=
        `E Sum = <span style="${hi}">${moduleE.toFixed(0)} MeV</span>; `
        + `NCl = <span style="${hi}">${clusters.length}</span>, `
        + `ECl Sum = <span style="${hi}">${clusterE.toFixed(0)} MeV</span>`;
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
    const rebins=Math.max(1,Math.ceil((rawEnergyMax-rawEnergyMin)/rawEnergyStep));
    rawEnergyBins=new Array(rebins).fill(0);
    currentRawEnergyHist=null;
}

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
            // Use ?? so legitimate 0 / 0.5 don't get clobbered by ||.
            nclustMin=data.nclusters.min ?? 0.5;
            nclustMax=data.nclusters.max ?? 10.5;
            nclustStep=data.nclusters.step ?? 1;
            nclustBins=data.nclusters.bins;
        }
        if(data.nblocks&&data.nblocks.bins&&data.nblocks.bins.length){
            nblocksMin=data.nblocks.min ?? 0;
            nblocksMax=data.nblocks.max ?? 40;
            nblocksStep=data.nblocks.step ?? 1;
            nblocksBins=data.nblocks.bins;
        }
        if(data.raw_energy&&data.raw_energy.bins&&data.raw_energy.bins.length){
            rawEnergyMin=data.raw_energy.min ?? 0;
            rawEnergyMax=data.raw_energy.max ?? 6000;
            rawEnergyStep=data.raw_energy.step ?? 20;
            rawEnergyBins=data.raw_energy.bins;
        }
        // Per-Ncl bucket arrays (added 2026-04 — older servers omit them,
        // in which case selection just falls back to the unfiltered hist).
        clEnergyBinsByNcl = Array.isArray(data.bins_by_ncl)
            ? data.bins_by_ncl : null;
        nblocksBinsByNcl = (data.nblocks && Array.isArray(data.nblocks.bins_by_ncl))
            ? data.nblocks.bins_by_ncl : null;
        // If the selected bucket no longer exists (e.g. server reconfig),
        // fall back to unfiltered.
        if (selectedNcl >= 0 && (!clEnergyBinsByNcl
                || selectedNcl >= clEnergyBinsByNcl.length)) {
            selectedNcl = -1;
        }
        plotClHist(); plotClStatHists(); plotRawEnergyHist();
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

// Pick the energy hist to display: a per-Ncl bucket if one is selected
// and available, otherwise the unfiltered one.  Returns {bins, label}
// where `label` annotates the title (empty for unfiltered).
function selectedClEnergy(){
    if (selectedNcl >= 0 && clEnergyBinsByNcl
        && selectedNcl < clEnergyBinsByNcl.length) {
        const ncl = nclustValueAt(selectedNcl);
        return { bins: clEnergyBinsByNcl[selectedNcl],
                 label: ` · Ncl=${ncl}` };
    }
    return { bins: clHistBins, label: '' };
}

function selectedNblocks(){
    if (selectedNcl >= 0 && nblocksBinsByNcl
        && selectedNcl < nblocksBinsByNcl.length) {
        const ncl = nclustValueAt(selectedNcl);
        return { bins: nblocksBinsByNcl[selectedNcl],
                 label: ` · Ncl=${ncl}` };
    }
    return { bins: nblocksBins, label: '' };
}

// Convert a bucket index to its Ncl value (bin center).  Used in titles
// so the user sees "Ncl=2", not "bucket 1".
function nclustValueAt(bucketIdx){
    return Math.round(nclustMin + (bucketIdx + 0.5) * nclustStep);
}

function plotClHist(){
    const div='cl-energy-hist';
    const sel=selectedClEnergy();
    const bins=sel.bins;
    if(!bins||!bins.length){
        currentClHist=null;
        Plotly.react(div,[],{...PL,title:{text:'Cluster Energy — No data',font:{size:10,color:THEME.textMuted}}},PC2);
        return;
    }
    const x=bins.map((_,i)=>clHistMin+(i+0.5)*clHistStep);
    const entries=bins.reduce((a,b)=>a+b,0);
    // store non-zero for copy
    const cx=[],cy=[];
    for(let i=0;i<bins.length;i++){if(bins[i]>0){cx.push(x[i]);cy.push(bins[i]);}}
    currentClHist={x:cx,y:cy};

    Plotly.react(div,[{
        x,y:bins,type:'bar',marker:{color:'#ff922b',line:{width:0}},
        hovertemplate:'%{x:.0f} MeV: %{y}<extra></extra>',
    }],{...PL,
        title:{text:`Cluster Energy${sel.label}<br><span style="font-size:9px;color:var(--theme-text-dim)">${clHistEvents} evts | ${entries} clusters</span>`,font:{size:10,color:THEME.textDim}},
        xaxis:{...PL.xaxis,title:'Energy (MeV)',range:[clHistMin,clHistMax]},
        yaxis:{...PL.yaxis,title:'Counts',
            type:document.getElementById('clhist-logy').checked?'log':'linear'},
        bargap:0.05,
        shapes:refShapes('cluster_energy'),
    },PC2);
}

function plotRawEnergyHist(){
    const div='cl-rawe-hist';
    const bins=rawEnergyBins;
    if(!bins||!bins.length){
        currentRawEnergyHist=null;
        Plotly.react(div,[],{...PL,title:{text:'Raw Energy Sum — No data',font:{size:10,color:THEME.textMuted}}},PC2);
        return;
    }
    const x=bins.map((_,i)=>rawEnergyMin+(i+0.5)*rawEnergyStep);
    const entries=bins.reduce((a,b)=>a+b,0);
    const cx=[],cy=[];
    for(let i=0;i<bins.length;i++){if(bins[i]>0){cx.push(x[i]);cy.push(bins[i]);}}
    currentRawEnergyHist={x:cx,y:cy};

    Plotly.react(div,[{
        x,y:bins,type:'bar',marker:{color:'#ffa94d',line:{width:0}},
        hovertemplate:'%{x:.0f} MeV: %{y}<extra></extra>',
    }],{...PL,
        title:{text:`Raw Energy Sum<br><span style="font-size:9px;color:var(--theme-text-dim)">${entries} evts</span>`,font:{size:10,color:THEME.textDim}},
        xaxis:{...PL.xaxis,title:'ΣE (MeV)',range:[rawEnergyMin,rawEnergyMax]},
        yaxis:{...PL.yaxis,title:'Counts',
            type:document.getElementById('clrawe-logy').checked?'log':'linear'},
        bargap:0.05,
        shapes:refShapes('raw_energy'),
    },PC2);
}

function plotClStatHists(){
    // Generic bar histogram.  `selectedIdx` (optional) highlights the
    // chosen bar in HIGHLIGHT colour; titlePrefix/Suffix get prepended /
    // appended to the title.
    function plotStat(divId, bins, bmin, bstep, title, xTitle, baseColor,
                      refKey, opts={}){
        if(!bins||!bins.length) return null;
        // Bin centers — works whether the user picked an integer-aligned
        // range (e.g. 0.5..10.5/1 → 1, 2, …) or a normal one.
        const x=bins.map((_,i)=>bmin+(i+0.5)*bstep);
        const entries=bins.reduce((a,b)=>a+b,0);
        const cx=[],cy=[];
        for(let i=0;i<bins.length;i++){if(bins[i]>0){cx.push(x[i]);cy.push(bins[i]);}}
        // Per-bar colour: highlight the selected bin so the user can
        // see at a glance which slice the dependent hists are showing.
        const colors = (opts.selectedIdx>=0 && opts.selectedIdx<bins.length)
            ? bins.map((_,i)=>i===opts.selectedIdx?THEME.highlight:baseColor)
            : baseColor;
        const fullTitle = (opts.titleSuffix||'')
            ? `${title}${opts.titleSuffix}` : title;
        Plotly.react(divId,[{
            x,y:bins,type:'bar',marker:{color:colors,line:{width:0}},
            hovertemplate:(opts.hoverFmt||'%{x}: %{y}<extra></extra>'),
        }],{...PL,
            title:{text:`${fullTitle}<br><span style="font-size:9px;color:var(--theme-text-dim)">${entries} entries${opts.titleHint||''}</span>`,font:{size:10,color:THEME.textDim}},
            xaxis:{...PL.xaxis,title:xTitle,
                   range:[bmin, bmin+bins.length*bstep]},
            yaxis:{...PL.yaxis,title:'Counts'},bargap:0.05,
            shapes:refKey?refShapes(refKey):[],
        },PC2);
        return {x:cx,y:cy};
    }
    const nblocksSel = selectedNblocks();
    currentNclustHist=plotStat('cl-nclust-hist',nclustBins,nclustMin,nclustStep,
        'Clusters per Event','# Clusters','#00b4d8','cluster_number',
        { selectedIdx: selectedNcl,
          hoverFmt: '%{x:.0f} clusters: %{y}<extra></extra>',
          titleHint: selectedNcl>=0
            ? ` · click again or another bar to change` : ` · click a bar to filter` });
    currentNblocksHist=plotStat('cl-nblocks-hist',nblocksSel.bins,nblocksMin,nblocksStep,
        'Blocks per Cluster','# Blocks','#51cf66','cluster_size',
        { titleSuffix: nblocksSel.label });

    // Wire up the click handler on the Ncl histogram once.  Plotly's
    // graphDiv keeps event subscriptions across Plotly.react calls, so
    // we only need to bind on the first paint.  Cache the flag on the
    // div itself so repeat calls don't stack listeners.
    const nclDiv = document.getElementById('cl-nclust-hist');
    if (nclDiv && nclDiv.on && !nclDiv._nclickBound) {
        nclDiv.on('plotly_click', ev => {
            if (!ev || !ev.points || !ev.points.length) return;
            const idx = ev.points[0].pointIndex;
            // Toggle: clicking the already-selected bar deselects.
            selectedNcl = (selectedNcl === idx) ? -1 : idx;
            plotClHist();
            plotClStatHists();
        });
        nclDiv._nclickBound = true;
    }
}

// =========================================================================
// GEM↔HyCal residuals (4 small panels at the top of the cluster tab)
// =========================================================================

function fetchGemResiduals(){
    fetch('/api/gem/residuals').then(r=>r.json()).then(data=>{
        if(!data.enabled) return;
        gemResidData=data;
        plotGemResiduals();
    }).catch(()=>{});
}

// Mean/sigma from a binned 1D histogram using bin centers.
function _residStats(h){
    if(!h||!h.bins||!h.bins.length) return {n:0,mean:0,sigma:0};
    let n=0,sum=0,sumSq=0;
    for(let i=0;i<h.bins.length;i++){
        const x=h.min+(i+0.5)*h.step;
        const c=h.bins[i];
        n+=c; sum+=c*x; sumSq+=c*x*x;
    }
    if(n===0) return {n:0,mean:0,sigma:0};
    const mean=sum/n;
    const variance=sumSq/n-mean*mean;
    return {n,mean,sigma:Math.sqrt(Math.max(0,variance))};
}

function _residTrace(h, color, name){
    const x=[],y=[];
    for(let i=0;i<h.bins.length;i++){
        x.push(h.min+(i+0.5)*h.step);
        y.push(h.bins[i]);
    }
    return {x,y,type:'scatter',mode:'lines',line:{shape:'hvh',color,width:1.5},
        name,hovertemplate:`${name}=%{x:.2f} mm: %{y}<extra></extra>`};
}

function plotGemResiduals(){
    for(let d=0;d<4;d++){
        const div='gem-resid-'+d;
        // GEM_COLORS is defined in gem.js (loaded before cluster.js); use it
        // here so each title's detector tag matches the GEM tab's color code.
        const gemColor=(typeof GEM_COLORS!=='undefined' && GEM_COLORS[d]) || THEME.text;
        const det=gemResidData && gemResidData.detectors && gemResidData.detectors[d];
        if(!det || !det.dx_hist || !det.dx_hist.bins || !det.dx_hist.bins.length){
            Plotly.react(div,[],{...PL,
                title:{text:`<span style="color:${gemColor}">GEM${d}</span> — No data`,
                       font:{size:10,color:THEME.textDim},
                       x:0.5,xanchor:'center',xref:'paper'},
                margin:{l:35,r:8,t:24,b:24}},PC2);
            continue;
        }
        const events=gemResidData.events||0;
        const sx=_residStats(det.dx_hist), sy=_residStats(det.dy_hist);
        const meanN=events>0 ? (det.matched_hits||0)/events : 0;
        const fmt=v=>(v>=0?' ':'')+v.toFixed(2);
        const titleText=`<span style="color:${gemColor}">${det.name||'GEM'+d}</span>`;
        // Mean/σ panel anchored to the top-right of the plot area.  Larger
        // font than the title so the numbers are legible at a glance; one
        // colored line per axis (X = blue, Y = red — matches the trace).
        // Third row is ⟨N⟩ (avg matched hits per event for this detector).
        const statsAnnotation={
            xref:'paper',yref:'paper',
            x:0.98,y:0.97,xanchor:'right',yanchor:'top',
            align:'right',showarrow:false,
            text:
                `<span style="color:#4dabf7">μₓ=${fmt(sx.mean)} σₓ=${sx.sigma.toFixed(2)}</span><br>`
               +`<span style="color:#ff6b6b">μᵧ=${fmt(sy.mean)} σᵧ=${sy.sigma.toFixed(2)}</span><br>`
               +`⟨N⟩=${meanN.toFixed(2)}`,
            font:{size:12,color:THEME.text},
            bgcolor:'rgba(0,0,0,0)',
        };
        Plotly.react(div,[
            _residTrace(det.dx_hist,'#4dabf7','ΔX'),
            _residTrace(det.dy_hist,'#ff6b6b','ΔY'),
        ],{...PL,
            title:{text:titleText,font:{size:11,color:THEME.text},
                   x:0.5,xanchor:'center',xref:'paper'},
            xaxis:{...PL.xaxis,title:'Residual (mm)',
                range:[det.dx_hist.min, det.dx_hist.min+det.dx_hist.bins.length*det.dx_hist.step]},
            yaxis:{...PL.yaxis,title:'Counts'},
            margin:{l:38,r:8,t:24,b:30},
            showlegend:false,
            annotations:[statsAnnotation],
        },PC2);
    }
}

// Theme flip — plot titles/traces bake THEME at draw time (e.g. THEME.textDim
// in titles, THEME.highlight in selected bars, THEME.text in residual stats).
// Replay each plot function from cached *Bins/*Data so the new theme reaches
// every text/marker, not just the chrome.  geoCluster paints the SHARED geo
// canvas — viewer.js's master listener handles that via redrawGeo() so we
// don't clobber the active tab's geo here.
if (typeof onThemeChange === 'function') {
    onThemeChange(() => {
        plotClHist();
        plotRawEnergyHist();
        plotClStatHists();
        plotGemResiduals();
    });
}
