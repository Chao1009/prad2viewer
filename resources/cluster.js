// cluster.js — Clustering tab: geo provider, cluster data, energy histograms

// --- clustering tab state ---
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

// build a set of module indices belonging to a cluster
function clusterModuleSet(clIdx){
    if(!clusterData||clIdx<0) return null;
    const cl=clusterData.clusters[clIdx];
    if(!cl) return null;
    return new Set(cl.modules);
}

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
        // refresh cluster histograms from server (accumulated there)
        fetchClHist();
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
        xaxis:{...PL.xaxis,title:'Energy (MeV)',range:[clHistMin,clHistMax]},
        yaxis:{...PL.yaxis,title:'Counts',
            type:document.getElementById('clhist-logy').checked?'log':'linear'},
        bargap:0.05,
        shapes:refShapes('cluster_energy'),
    },PC2);
}

function plotClStatHists(){
    // number of clusters histogram
    function plotStat(divId, bins, bmin, bstep, title, xTitle, color, refKey){
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
            shapes:refKey?refShapes(refKey):[],
        },PC2);
        return {x:cx,y:cy};
    }
    currentNclustHist=plotStat('cl-nclust-hist',nclustBins,nclustMin,nclustStep,
        'Clusters per Event','# Clusters','#00b4d8','cluster_number');
    currentNblocksHist=plotStat('cl-nblocks-hist',nblocksBins,nblocksMin,nblocksStep,
        'Blocks per Cluster','# Blocks','#51cf66','cluster_size');
}
