// physics.js — Physics tab: energy vs angle + Møller XY + Møller energy
//
// Depends on globals from viewer.js: PL, PC_EPICS, activeTab

let physicsData=null, mollerData=null;

function fetchEnergyAngle(){
    fetch('/api/physics/energy_angle').then(r=>r.json()).then(data=>{
        physicsData=data;
        plotEnergyAngle();
    }).catch(()=>{});
}

function fetchMoller(){
    fetch('/api/physics/moller').then(r=>r.json()).then(data=>{
        mollerData=data;
        plotMollerXY();
        plotMollerEnergy();
    }).catch(()=>{});
}

function fetchPhysics(){
    fetchEnergyAngle();
    fetchMoller();
}

// ep elastic scattering: E' = E / (1 + (E/Mp)*(1 - cos(theta)))
function elasticEp(beamE, thetaDeg){
    const Mp=938.272;
    const th=thetaDeg*Math.PI/180;
    return beamE/(1+(beamE/Mp)*(1-Math.cos(th)));
}

function plotEnergyAngle(){
    const div='physics-plot';
    if(!physicsData||!physicsData.bins||!physicsData.bins.length||!physicsData.nx){
        Plotly.react(div,[],{...PL,title:{text:'Energy vs Angle — No data',font:{size:12,color:'#888'}}},PC_EPICS);
        document.getElementById('physics-stats').textContent='';
        return;
    }
    const d=physicsData;
    const logZ=document.getElementById('physics-logz').checked;
    const showElastic=document.getElementById('physics-elastic').checked;

    const z=[];
    for(let iy=0;iy<d.ny;iy++){
        const row=d.bins.slice(iy*d.nx,(iy+1)*d.nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
    }
    const x=[];for(let i=0;i<d.nx;i++) x.push(d.angle_min+(i+0.5)*d.angle_step);
    const y=[];for(let i=0;i<d.ny;i++) y.push(d.energy_min+(i+0.5)*d.energy_step);

    const traces=[{
        z:z, x:x, y:y,
        type:'heatmap', colorscale:'Hot', reversescale:false,
        hovertemplate:'θ=%{x:.2f}° E=%{y:.0f} MeV: %{text}<extra></extra>',
        text:z.map((row,iy)=>row.map((v,ix)=>String(d.bins[iy*d.nx+ix]))),
        colorbar:{title:logZ?'log₁₀(counts)':'counts',titleside:'right',
            titlefont:{size:10,color:'#aaa'},tickfont:{size:9,color:'#aaa'}},
    }];

    if(showElastic && d.beam_energy>0){
        const ex=[],ey=[];
        for(let th=d.angle_min+0.1;th<=d.angle_max;th+=0.05){
            const e=elasticEp(d.beam_energy,th);
            if(e>=d.energy_min&&e<=d.energy_max){ex.push(th);ey.push(e);}
        }
        traces.push({x:ex,y:ey,mode:'lines',
            line:{color:'#00ff88',width:2,dash:'dot'},
            name:`ep elastic (${d.beam_energy} MeV)`,
            hovertemplate:'θ=%{x:.2f}° E=%{y:.0f} MeV<extra>ep elastic</extra>'});
    }

    Plotly.react(div,traces,{...PL,
        title:{text:`Energy vs Angle (${d.events} evts)`,font:{size:12,color:'#ccc'}},
        xaxis:{...PL.xaxis,title:'Scattering Angle (deg)'},
        yaxis:{...PL.yaxis,title:'Energy (MeV)'},
        margin:{l:55,r:80,t:30,b:40},
        showlegend:showElastic,
        legend:{x:0.7,y:0.95,font:{size:10,color:'#aaa'},bgcolor:'rgba(0,0,0,0)'},
        shapes:refShapes('energy_angle'),
    },PC_EPICS);

    // stats line
    const ml=mollerData;
    let stats=`${d.events} evts | beam: ${d.beam_energy||'?'} MeV`;
    if(ml) stats+=` | Møller: ${ml.moller_events}`;
    document.getElementById('physics-stats').textContent=stats;
}

function plotMollerXY(){
    const div='moller-xy-plot';
    const d=mollerData;
    if(!d||!d.xy_bins||!d.xy_bins.length||!d.xy_nx){
        Plotly.react(div,[],{...PL,title:{text:'Møller XY — No data',font:{size:12,color:'#888'}}},PC_EPICS);
        return;
    }
    const logZ=document.getElementById('physics-logz').checked;
    const z=[];
    for(let iy=0;iy<d.xy_ny;iy++){
        const row=d.xy_bins.slice(iy*d.xy_nx,(iy+1)*d.xy_nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
    }
    const x=[];for(let i=0;i<d.xy_nx;i++) x.push(d.xy_x_min+(i+0.5)*d.xy_x_step);
    const y=[];for(let i=0;i<d.xy_ny;i++) y.push(d.xy_y_min+(i+0.5)*d.xy_y_step);

    const cuts=d.cuts||{};
    const fmtA=v=>v!=null?v.toFixed(2):'?';
    const cutTxt=`θ∈[${fmtA(cuts.angle_min)},${fmtA(cuts.angle_max)}]° Esum±${((cuts.energy_tolerance||0.1)*100).toFixed(0)}%`;

    Plotly.react(div,[{
        z:z, x:x, y:y,
        type:'heatmap', colorscale:'Hot', reversescale:false,
        hovertemplate:'x=%{x:.1f} y=%{y:.1f} mm: %{text}<extra></extra>',
        text:z.map((row,iy)=>row.map((v,ix)=>String(d.xy_bins[iy*d.xy_nx+ix]))),
        colorbar:{title:logZ?'log₁₀':'counts',titleside:'right',
            titlefont:{size:10,color:'#aaa'},tickfont:{size:9,color:'#aaa'}},
    }],{...PL,
        title:{text:`Møller XY (${d.moller_events} evts) ${cutTxt}`,font:{size:11,color:'#ccc'}},
        xaxis:{...PL.xaxis,title:'X (mm)',scaleanchor:'y',scaleratio:1},
        yaxis:{...PL.yaxis,title:'Y (mm)'},
        margin:{l:50,r:70,t:30,b:35},
        shapes:refShapes('moller_xy'),
    },PC_EPICS);
}

function plotMollerEnergy(){
    const div='moller-energy-plot';
    const d=mollerData;
    if(!d||!d.energy_hist||!d.energy_hist.bins||!d.energy_hist.bins.length){
        Plotly.react(div,[],{...PL,title:{text:'Møller Energy — No data',font:{size:12,color:'#888'}}},PC_EPICS);
        return;
    }
    const h=d.energy_hist;
    const x=[];for(let i=0;i<h.bins.length;i++) x.push(h.min+(i+0.5)*h.step);

    Plotly.react(div,[{
        x:x, y:h.bins,
        type:'bar',
        marker:{color:'#00b4d8'},
        hovertemplate:'E=%{x:.0f} MeV: %{y}<extra></extra>',
    }],{...PL,
        title:{text:`Møller Cluster Energy (${d.moller_events} evts)`,font:{size:11,color:'#ccc'}},
        xaxis:{...PL.xaxis,title:'Energy (MeV)'},
        yaxis:{...PL.yaxis,title:'Counts'},
        margin:{l:50,r:20,t:28,b:35},
        bargap:0,
        shapes:refShapes('moller_energy'),
    },PC_EPICS);
}

function clearPhysicsFrontend(){
    physicsData=null; mollerData=null;
    Plotly.react('physics-plot',[],{...PL},PC_EPICS);
    Plotly.react('moller-xy-plot',[],{...PL},PC_EPICS);
    Plotly.react('moller-energy-plot',[],{...PL},PC_EPICS);
    document.getElementById('physics-stats').textContent='';
}

function resizePhysics(){
    try{Plotly.Plots.resize('physics-plot');}catch(e){}
    try{Plotly.Plots.resize('moller-xy-plot');}catch(e){}
    try{Plotly.Plots.resize('moller-energy-plot');}catch(e){}
}

function initPhysics(data){
    document.getElementById('physics-logz').onchange=()=>{plotEnergyAngle();plotMollerXY();plotMollerEnergy();};
    document.getElementById('physics-elastic').onchange=plotEnergyAngle;
}
