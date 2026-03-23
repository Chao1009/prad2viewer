// physics.js — Physics tab: energy vs angle 2D histogram
//
// Depends on globals from viewer.js: PL, PC_EPICS (reuse mode bar config), activeTab

let physicsData=null;

function fetchEnergyAngle(){
    fetch('/api/physics/energy_angle').then(r=>r.json()).then(data=>{
        physicsData=data;
        plotEnergyAngle();
    }).catch(()=>{});
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

    // reshape flat bins to 2D [ny][nx]
    const z=[];
    for(let iy=0;iy<d.ny;iy++){
        const row=d.bins.slice(iy*d.nx,(iy+1)*d.nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
    }

    // axis tick values at bin centers
    const x=[];for(let i=0;i<d.nx;i++) x.push(d.angle_min+(i+0.5)*d.angle_step);
    const y=[];for(let i=0;i<d.ny;i++) y.push(d.energy_min+(i+0.5)*d.energy_step);

    Plotly.react(div,[{
        z:z, x:x, y:y,
        type:'heatmap',
        colorscale:'Hot',
        reversescale:false,
        hovertemplate:'θ=%{x:.2f}° E=%{y:.0f} MeV: %{text}<extra></extra>',
        text:z.map((row,iy)=>row.map((v,ix)=>{
            const raw=d.bins[iy*d.nx+ix];
            return String(raw);
        })),
        colorbar:{title:logZ?'log₁₀(counts)':'counts',titleside:'right',
            titlefont:{size:10,color:'#aaa'},tickfont:{size:9,color:'#aaa'}},
    }],{...PL,
        title:{text:`Energy vs Angle (${d.events} evts)`,font:{size:12,color:'#ccc'}},
        xaxis:{...PL.xaxis,title:'Scattering Angle (deg)'},
        yaxis:{...PL.yaxis,title:'Energy (MeV)'},
        margin:{l:55,r:80,t:30,b:40},
    },PC_EPICS);

    document.getElementById('physics-stats').textContent=
        `${d.events} events | HyCal z: ${(d.hycal_z||5800)/1000}m`;
}

function clearPhysicsFrontend(){
    physicsData=null;
    Plotly.react('physics-plot',[],{...PL},PC_EPICS);
    document.getElementById('physics-stats').textContent='';
}

function initPhysics(data){
    document.getElementById('physics-logz').onchange=plotEnergyAngle;
}
