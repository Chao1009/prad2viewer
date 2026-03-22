// report.js — Report generation: markdown + images
//
// Produces markdown text + PNG image attachments.
// Used by both "Download Report" (.md file) and "Post to Elog" (text + attachments).
//
// Depends on globals from viewer.js (accessed at runtime, not load time).

// =========================================================================
// Registry
// =========================================================================
const reportRegistry=[];
let elogConfig={url:'',logbook:'',author:'',tags:[]};
let reportAttachments=[];  // [{data (base64), filename, caption, type}]

function registerReportSection(section){
    reportRegistry.push(section);
    reportRegistry.sort((a,b)=>a.order-b.order);
}

// =========================================================================
// Capture helpers
// =========================================================================

// Capture the geo canvas for a given tab at fixed resolution with light theme.
async function captureGeoForTab(tab){
    const prev={tab:activeTab,w:geoCanvas.width,h:geoCanvas.height,
                s:scale,ox:offsetX,oy:offsetY};
    try{
        geoCanvas.width=1200; geoCanvas.height=900;
        canvasW=1200; canvasH=900;
        activeTab=tab;
        geoLightTheme=true;
        fitView(); redrawGeo();
        return geoCanvas.toDataURL('image/png');
    }finally{
        geoLightTheme=false;
        geoCanvas.width=prev.w; geoCanvas.height=prev.h;
        canvasW=prev.w; canvasH=prev.h;
        scale=prev.s; offsetX=prev.ox; offsetY=prev.oy;
        activeTab=prev.tab;
        redrawGeo();
    }
}

// Capture a geo view, register as attachment, return markdown image reference.
async function captureGeo(tab,caption,filename){
    const img=await captureGeoForTab(tab);
    addAttachment(img,filename,caption);
    return `![${caption}](${filename})\n\n`;
}

// Render a Plotly chart off-screen, capture as PNG data URL.
async function plotToImage(plotFn,w,h){
    const div=document.createElement('div');
    div.style.cssText='position:fixed;left:-9999px;width:'+w+'px;height:'+h+'px';
    document.body.appendChild(div);
    await plotFn(div);
    const img=await Plotly.toImage(div,{format:'png',width:w,height:h});
    Plotly.purge(div);
    div.remove();
    return img;
}

// Light-theme Plotly layout.
const RPL={paper_bgcolor:'#fff',plot_bgcolor:'#fff',
    font:{family:'Helvetica,Arial,sans-serif',size:11,color:'#222'},
    margin:{l:50,r:15,t:28,b:36},
    xaxis:{gridcolor:'#ddd',zerolinecolor:'#bbb',linecolor:'#999',mirror:true},
    yaxis:{gridcolor:'#ddd',zerolinecolor:'#bbb',linecolor:'#999',mirror:true}};

// Capture a bar histogram, register as attachment, return markdown reference.
// Returns empty string if no data.
async function captureHist(bins,binMin,binStep,title,xTitle,color,filename,w,h){
    if(!bins||!bins.some(b=>b>0)) return '';
    try{
        const img=await plotToImage(async div=>{
            const x=bins.map((_,i)=>binMin+(i+0.5)*binStep);
            const entries=bins.reduce((a,b)=>a+b,0);
            const titleText=entries?`${title} (${entries} entries)`:title;
            await Plotly.newPlot(div,[{x,y:bins,type:'bar',
                marker:{color,line:{width:0}}}],
                {...RPL,title:{text:titleText,font:{size:12,color:'#222'}},
                 xaxis:{...RPL.xaxis,title:xTitle},
                 yaxis:{...RPL.yaxis,title:'Counts'},bargap:0.05});
        },w,h);
        addAttachment(img,filename,title);
        return `![${title}](${filename})\n\n`;
    }catch(e){ return ''; }
}

function addAttachment(dataUrl,filename,caption){
    const b64=dataUrl.split(',')[1];
    if(b64) reportAttachments.push({data:b64,filename,caption,type:'image/png'});
}

// =========================================================================
// Markdown table helper
// =========================================================================
function mdTable(headers,rows,alignments){
    // headers: ['Col1','Col2',...]
    // alignments: optional array of 'l','r','c' per column (default 'l')
    const aligns=alignments||headers.map(()=>'l');
    const sepMap={l:':---',r:'---:',c:':---:'};
    let md='| '+headers.join(' | ')+' |\n';
    md+='| '+aligns.map(a=>sepMap[a]||':---').join(' | ')+' |\n';
    for(const row of rows)
        md+='| '+row.join(' | ')+' |\n';
    return md+'\n';
}

// =========================================================================
// Report sections
// =========================================================================

// Helper: temporarily set LMS metric dropdown, capture geo, return md with range info.
async function captureLmsGeo(metric,caption,filename){
    const sel=document.getElementById('lms-color-metric');
    const prev=sel.value;
    sel.value=metric;
    try{
        let md=await captureGeo('lms',caption,filename);
        const rMin=document.getElementById('lms-range-min-show').textContent;
        const rMax=document.getElementById('lms-range-max-show').textContent;
        const useLog=document.getElementById('lms-log-scale').checked;
        md+=`Color range: ${rMin} – ${rMax} | Scale: ${useLog?'log':'linear'}\n\n`;
        return md;
    }finally{ sel.value=prev; }
}

// Helper: find the ref index for a given ref channel name (e.g. 'LMS3').
function findRefIndex(name){
    const sel=document.getElementById('lms-ref-select');
    for(const o of sel.options)
        if(o.textContent===name) return parseInt(o.value);
    return -1;
}

// --- Occupancy ---
registerReportSection({id:'occupancy',title:'Occupancy',order:10,
    generate:async()=>{
        const prevMetric=document.getElementById('color-metric').value;
        document.getElementById('color-metric').value='occupancy';
        syncDqRange();
        let md='## Occupancy\n\n';
        if(occTotal>0) md+=`Total events: ${occTotal}\n\n`;
        // report color range and scale
        const useLog=document.getElementById('log-scale').checked;
        md+=`Color range: ${rangeMin??0} – ${rangeMax??'auto'} | Scale: ${useLog?'log':'linear'}\n\n`;
        md+=await captureGeo('dq','Occupancy','occupancy.png');
        document.getElementById('color-metric').value=prevMetric;
        syncDqRange();
        return md;
    }
});

// --- Clustering ---
registerReportSection({id:'cluster',title:'Clustering',order:20,
    generate:async()=>{
        if(!clHistBins||!clHistBins.some(b=>b>0)) return null;
        let md='## Clustering\n\n';
        md+=await captureHist(clHistBins,clHistMin,clHistStep,
            `Cluster Energy (${clHistEvents} evts)`,'Energy (MeV)','#ff922b',
            'cluster_energy.png',800,300);
        md+=await captureHist(nclustBins,nclustMin,nclustStep,
            'Clusters per Event','# Clusters','#00b4d8',
            'clusters_per_event.png',500,300);
        md+=await captureHist(nblocksBins,nblocksMin,nblocksStep,
            'Blocks per Cluster','# Blocks','#51cf66',
            'blocks_per_cluster.png',500,300);
        return md;
    }
});

// --- LMS Monitoring ---
registerReportSection({id:'lms',title:'LMS Monitoring',order:30,
    generate:async()=>{
        // fetch with LMS3 reference
        const lms3Ref=findRefIndex('LMS3');
        const refQ=lms3Ref>=0?`?ref=${lms3Ref}`:'';
        const d=await fetch(`/api/lms/summary${refQ}`).then(r=>r.json());
        lmsSummaryData=d;
        const lmsEvents=lmsSummaryData?lmsSummaryData.events||0:0;
        const trigBit=lmsSummaryData?lmsSummaryData.trigger_bit:'?';
        if(!lmsSummaryData||!lmsSummaryData.modules||
            !Object.keys(lmsSummaryData.modules).length)
            return `## LMS Monitoring\n\nLMS events received: ${lmsEvents} (trigger bit = ${trigBit})\n\n`;
        const allEntries=Object.entries(lmsSummaryData.modules)
            .map(([idx,m])=>({idx:parseInt(idx),...m}));

        // temporarily set ref so geo drawing uses the corrected data
        const prevRef=g_lmsRefIndex;
        g_lmsRefIndex=lms3Ref;

        let md='## LMS Monitoring\n\n';
        const refLabel=lms3Ref>=0?' (ref: LMS3)':'';
        md+=`LMS events: ${lmsSummaryData.events||0} | `;
        md+=`Modules: ${allEntries.length}${refLabel}\n\n`;

        // LMS Mean geo view — errors propagate (no silent catch)
        md+=await captureLmsGeo('mean','LMS Mean','lms_mean.png');
        // RMS/Mean geo view
        md+=await captureLmsGeo('rms_frac','RMS / Mean','lms_rms_frac.png');

        // restore ref
        g_lmsRefIndex=prevRef;

        // build table: all WARN + top N OK by rms/mean
        allEntries.sort((a,b)=>{
            if(a.warn!==b.warn) return a.warn?-1:1;
            const ra=a.mean>0?a.rms/a.mean:0, rb=b.mean>0?b.rms/b.mean:0;
            return rb-ra;
        });
        const warnEntries=allEntries.filter(e=>e.warn);
        const okEntries=allEntries.filter(e=>!e.warn);
        const topOk=okEntries.slice(0,5);
        const tableEntries=[...warnEntries,...topOk];

        const warnCount=warnEntries.length;
        md+=`### Status Summary\n\n`;
        md+=`Warnings: **${warnCount}** / ${allEntries.length} modules\n\n`;

        if(tableEntries.length){
            md+=mdTable(
                ['Module','Mean','RMS','RMS/Mean %','Count','Status'],
                tableEntries.map(e=>[
                    e.name, e.mean.toFixed(1), e.rms.toFixed(2),
                    (e.mean>0?(e.rms/e.mean*100).toFixed(1):'--')+'%',
                    e.count, e.warn?'**WARN**':'OK'
                ]),['l','r','r','r','r','l']
            );
            if(okEntries.length>5)
                md+=`*Showing ${warnCount} warnings + top 5 of ${okEntries.length} OK modules by RMS/Mean*\n\n`;
        }
        return md;
    }
});

// =========================================================================
// Report generation core
// =========================================================================

async function refreshDataForReport(){
    const fetches=[];
    fetches.push(fetch('/api/occupancy').then(r=>r.json()).then(d=>{
        occData=d.occ||{}; occTcutData=d.occ_tcut||{}; occTotal=d.total||0;
    }).catch(()=>{}));
    fetches.push(fetch('/api/cluster_hist').then(r=>r.json()).then(d=>{
        if(d.bins&&d.bins.length){
            if(d.min!==undefined) clHistMin=d.min;
            if(d.max!==undefined) clHistMax=d.max;
            if(d.step!==undefined) clHistStep=d.step;
            clHistBins=d.bins; clHistEvents=d.events||0;
        }
        if(d.nclusters&&d.nclusters.bins&&d.nclusters.bins.length){
            nclustMin=d.nclusters.min||0; nclustMax=d.nclusters.max||20;
            nclustStep=d.nclusters.step||1; nclustBins=d.nclusters.bins;
        }
        if(d.nblocks&&d.nblocks.bins&&d.nblocks.bins.length){
            nblocksMin=d.nblocks.min||0; nblocksMax=d.nblocks.max||40;
            nblocksStep=d.nblocks.step||1; nblocksBins=d.nblocks.bins;
        }
    }).catch(()=>{}));
    // LMS section fetches its own data with LMS3 ref
    await Promise.all(fetches);
}

// Generate the report. Returns {md, attachments} or null.
// reportBy: optional name string to include in header.
async function generateReport(reportBy){
    if(!modules.length){
        alert('No data loaded. Please load data before generating a report.');
        return null;
    }
    const statusBar=document.getElementById('status-bar');
    const prevStatus=statusBar.textContent;
    statusBar.textContent='Generating report...';
    try{
        await refreshDataForReport();
        reportAttachments=[];
        const ts=new Date().toLocaleString();
        const samples=mode==='online'?sampleCount:totalEvents;
        let header=`# PRad2 HyCal Monitor Report\n\n`;
        header+=`- **Generated:** ${ts}\n`;
        header+=`- **Samples:** ${samples}\n`;
        if(reportBy) header+=`- **Report by:** ${reportBy}\n`;
        let sectionsMd='';
        for(const entry of reportRegistry){
            try{
                const section=await entry.generate();
                if(section) sectionsMd+=section;
            }catch(err){
                sectionsMd+=`## ${entry.title}\n\n*Error: ${err.message}*\n\n`;
            }
        }
        // append LMS warn summary to header (data available after LMS section runs)
        if(lmsSummaryData&&lmsSummaryData.modules){
            const warns=Object.values(lmsSummaryData.modules)
                .filter(m=>m.warn).map(m=>m.name)
                .sort((a,b)=>{
                    // W modules first, then G
                    const ta=a.startsWith('W')?0:1, tb=b.startsWith('W')?0:1;
                    if(ta!==tb) return ta-tb;
                    return a.localeCompare(b,undefined,{numeric:true});
                });
            if(warns.length)
                header+=`- **Gain Monitoring Warnings (${warns.length}):** ${warns.join(', ')}\n`;
            else
                header+=`- **Gain Monitoring Warnings:** None\n`;
        }
        let md=header+`\n---\n\n`+sectionsMd;
        md+=`---\n*PRad2 HyCal Online Monitor — Report generated ${ts}*\n`;
        statusBar.textContent=prevStatus;
        return {md, attachments:reportAttachments};
    }catch(err){
        statusBar.textContent=`Report error: ${err.message}`;
        return null;
    }
}

// =========================================================================
// Download report (.md + images)
// =========================================================================

function downloadBlob(blob,filename){
    const a=document.createElement('a');
    a.href=URL.createObjectURL(blob);
    a.download=filename;
    a.style.display='none';
    document.body.appendChild(a);
    a.click();
    // delay removal and revoke to let browser finish saving
    setTimeout(()=>{a.remove();URL.revokeObjectURL(a.href);},30000);
}

function b64toBlob(b64,type){
    const bin=atob(b64);
    const arr=new Uint8Array(bin.length);
    for(let i=0;i<bin.length;i++) arr[i]=bin.charCodeAt(i);
    return new Blob([arr],{type});
}

async function downloadReport(){
    const reportBy=prompt('Report by (your name, optional):','');
    const report=await generateReport(reportBy||'');
    if(!report) return;
    const statusBar=document.getElementById('status-bar');
    const ts=new Date().toISOString().slice(0,19).replace(/[T:]/g,'-');

    // download .md file
    downloadBlob(new Blob([report.md],{type:'text/markdown'}),
        `prad2_report_${ts}.md`);

    // download each image with delay to avoid browser blocking
    for(let i=0;i<report.attachments.length;i++){
        const a=report.attachments[i];
        await new Promise(r=>setTimeout(r,500));
        downloadBlob(b64toBlob(a.data,a.type),a.filename);
    }
    statusBar.textContent=`Report saved (${1+report.attachments.length} files)`;
    setTimeout(()=>{statusBar.textContent='Ready';},3000);
}

// =========================================================================
// Elog posting
// =========================================================================

function escXml(s){
    return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;')
        .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

function buildElogXml(title,logbook,author,tags,body,attachments){
    const parts=['<?xml version="1.0" encoding="UTF-8"?>','<Logentry>',
        `  <created>${new Date().toISOString()}</created>`,
        `  <Author><username>${escXml(author)}</username></Author>`,
        `  <title>${escXml(title)}</title>`,
        `  <body type="text"><![CDATA[${body}]]></body>`,
        '  <Logbooks>'];
    for(const lb of logbook.split(','))
        parts.push(`    <logbook>${escXml(lb.trim())}</logbook>`);
    parts.push('  </Logbooks>');
    if(tags&&tags.length){
        parts.push('  <Tags>');
        for(const t of tags) parts.push(`    <tag>${escXml(t.trim())}</tag>`);
        parts.push('  </Tags>');
    }
    if(attachments&&attachments.length){
        parts.push('  <Attachments>');
        for(const a of attachments)
            parts.push('    <Attachment>',
                `      <caption>${escXml(a.caption)}</caption>`,
                `      <filename>${escXml(a.filename)}</filename>`,
                `      <type>${escXml(a.type)}</type>`,
                `      <data encoding="base64">${a.data}</data>`,
                '    </Attachment>');
        parts.push('  </Attachments>');
    }
    parts.push('</Logentry>');
    return parts.join('\n');
}

function showElogDialog(){
    document.getElementById('elog-backdrop').classList.add('open');
    document.getElementById('elog-dialog').classList.add('open');
    document.getElementById('elog-status').textContent='';
}
function hideElogDialog(){
    document.getElementById('elog-backdrop').classList.remove('open');
    document.getElementById('elog-dialog').classList.remove('open');
}

async function postToElog(){
    const reportBy=document.getElementById('elog-report-by').value.trim();
    const title=document.getElementById('elog-title').value.trim();
    const logbook=document.getElementById('elog-logbook').value.trim();
    const author=document.getElementById('elog-author').value.trim();
    const tagsStr=document.getElementById('elog-tags').value.trim();
    const tags=tagsStr?tagsStr.split(',').map(s=>s.trim()).filter(s=>s):[];
    const statusEl=document.getElementById('elog-status');

    if(!title||!logbook||!author){
        statusEl.textContent='Title, logbook, and posting account are required.';
        statusEl.style.color='#c00';
        return;
    }
    statusEl.textContent='Generating report...';
    statusEl.style.color='var(--dim)';

    const report=await generateReport(reportBy);
    if(!report){statusEl.textContent='Failed to generate report.';statusEl.style.color='#c00';return;}

    // strip image markdown references from body (images go as attachments)
    const body=report.md.replace(/!\[[^\]]*\]\([^)]+\)\n*/g,'');

    statusEl.textContent='Posting to elog...';
    const xml=buildElogXml(title,logbook,author,tags,body,report.attachments);

    try{
        const resp=await fetch('/api/elog/post',{
            method:'POST',headers:{'Content-Type':'application/xml'},body:xml});
        const result=await resp.json();
        if(result.ok){
            statusEl.textContent='Posted successfully! (HTTP '+result.status+')';
            statusEl.style.color='#080';
            setTimeout(hideElogDialog,2000);
        }else{
            statusEl.textContent='Post failed: HTTP '+result.status+(result.error?' - '+result.error:'');
            statusEl.style.color='#c00';
        }
    }catch(err){
        statusEl.textContent='Network error: '+err.message;
        statusEl.style.color='#c00';
    }
}

// =========================================================================
// Init — called from viewer.js init() with config data
// =========================================================================
function initReport(data){
    const reportBtn=document.getElementById('btn-report');
    const reportMenu=document.getElementById('report-menu');
    reportBtn.onclick=(e)=>{e.stopPropagation();reportMenu.classList.toggle('open');};
    document.addEventListener('click',()=>reportMenu.classList.remove('open'));
    reportMenu.onclick=(e)=>e.stopPropagation();
    document.getElementById('btn-report-pdf').onclick=()=>{
        reportMenu.classList.remove('open'); downloadReport();};
    document.getElementById('btn-report-elog').onclick=()=>{
        reportMenu.classList.remove('open'); showElogDialog();};

    document.getElementById('elog-dialog-close').onclick=hideElogDialog;
    document.getElementById('elog-backdrop').onclick=hideElogDialog;
    document.getElementById('elog-cancel').onclick=hideElogDialog;
    document.getElementById('elog-submit').onclick=postToElog;

    if(data&&data.elog&&data.elog.url){
        elogConfig=data.elog;
        document.getElementById('elog-logbook').value=data.elog.logbook||'';
        document.getElementById('elog-author').value=data.elog.author||'';
        document.getElementById('elog-tags').value=(data.elog.tags||[]).join(', ');
    } else {
        const eb=document.getElementById('btn-report-elog');
        eb.disabled=true;
        eb.title='Configure "elog" section in config.json to enable';
    }
}
