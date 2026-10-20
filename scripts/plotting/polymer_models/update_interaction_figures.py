#!/usr/bin/env python
"""Maintained electronic interaction plots and an offline rotatable 3D viewer."""
import hashlib
import html
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.data import covalent_radii
from plotly.offline import get_plotlyjs

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'results/polymer_models'


def main():
    subprocess.run([sys.executable,str(ROOT/'scripts/validation/review_polymer_interactions.py')],check=True)
    feed=json.loads((OUT/'interaction_checked.json').read_text())
    figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
    assets=OUT/'assets';assets.mkdir(exist_ok=True)
    js=assets/'plotly.min.js'
    if not js.exists():js.write_text(get_plotlyjs())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    for kind,field,title,xlabel in [
        ('folding','stacked_minus_open_deltaE_kJ_mol','Stacked versus open covalent conformers','Stacked − open electronic ΔE (kJ mol⁻¹)'),
        ('association','association_deltaE_kJ_mol','Separate-segment association','Electronic association ΔE (kJ mol⁻¹)')]:
        values=[r for r in feed[kind] if r[field] is not None]
        fig,ax=plt.subplots(figsize=(8.4,4.6),layout='constrained')
        if values:
            labels=[r['candidate_id'].replace('design_','').replace('exp_','')+' / '+
                    (f"n{r['repeats']}, stereo {r['stereoisomer']}" if kind=='folding' else r['orientation']) for r in values]
            ax.scatter([r[field] for r in values],range(len(values)),color='#0072B2',s=45)
            ax.set_yticks(range(len(values)),labels,fontsize=9);ax.invert_yaxis();ax.axvline(0,color='.6',lw=1)
        else:
            ax.text(.5,.5,'Awaiting matched, stationary DFT structures',ha='center',va='center',transform=ax.transAxes,color='.4')
            ax.set_xticks([]);ax.set_yticks([])
        ax.set_title(title);ax.set_xlabel(xlabel)
        fig.supxlabel('Neutral MeCN method controls. Electronic energies; entropy, Hessian minima and association BSSE remain unqualified.',fontsize=8)
        for suffix in ('pdf','svg','png'):fig.savefig(figdir/f'{kind}.{suffix}',dpi=350)
        plt.close(fig)
    models=[]
    for row in feed['models']:
        path=Path(row['final_xyz'] or row.get('input_xyz',''))
        if not path.is_file():
            ensemble=ROOT/feed['status']['campaign']/'seeds'/row['model_key']/'ensemble.json'
            if not ensemble.exists() or row['model_type']!='covalent_oligomer':continue
            sys.path.insert(0,str(ROOT/'scripts/pipeline'))
            from prepare_pi_stacking import choose
            try:preview,fallback=choose(json.loads(ensemble.read_text()),row['target'])
            except (ValueError,AssertionError):continue
            row=dict(row,pi_groups=preview['pi_groups'],seed_stage='biased starting geometry; DFT pending' if fallback else 'xTB conformer preview; DFT pending')
            path=Path(preview['xyz'])
        atoms=read(path);positions=atoms.positions;elements=atoms.get_chemical_symbols();bonds=[]
        for i in range(len(atoms)):
            for j in range(i):
                distance=np.linalg.norm(positions[i]-positions[j])
                if .3<distance<1.2*(covalent_radii[atoms.numbers[i]]+covalent_radii[atoms.numbers[j]]):bonds.append([i,j])
        if not row['final_xyz']:
            from sys import path as python_path
            python_path.insert(0,str(ROOT/'src'))
            from redox.build.polymer_contacts import contact_metrics
            metrics=contact_metrics(positions,row['pi_groups'])
        else:metrics=row['contact_metrics']
        stage=row['stage'] if row['final_xyz'] else (
            'DFT failed (preserved) · ' if row['stage']=='failed_preserved' else 'DFT pending · ')+row.get('seed_stage','prepared starting geometry')
        models.append(dict(key=row['key'],label=row['candidate_id']+' · '+str(row['repeats'])+' centers · '+str(row['stereoisomer'])+' · '+row['target'],
            positions=positions.tolist(),elements=elements,bonds=bonds,groups=row['pi_groups'],
            stage=stage,metrics=metrics))
    status=feed['status']
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Polymer interactions</title>
<script src="assets/plotly.min.js"></script><style>
body{font:16px/1.5 system-ui,sans-serif;max-width:1200px;margin:30px auto;padding:0 20px;color:#203040;background:#f8fafb}
a{color:#0072B2}select{font:inherit;max-width:100%;padding:10px}#viewer{height:600px;width:100%;background:white}#metrics{overflow-wrap:anywhere}.plots{display:grid;grid-template-columns:1fr 1fr;gap:20px}.plots img{width:100%}.note{padding:14px;background:#fff4df;border-left:4px solid #b87500} @media(max-width:700px){.plots{display:block}#viewer{height:450px}}
</style><h1>Polymer-site interactions</h1><p><a href="../README.md">All results</a> · <a href="README.md">Models and scope</a> · <a href="summary.html">Model construction and key results</a> · <a href="interaction_status.json">Current status</a></p>'''
    page+=f"<p><strong>{status['stationary_electronic_cases']}/{status['expected_DFT_cases']} stationary DFT electronic structures collected.</strong> Neutral NDI, PMDI and phthalimide controls; monomer, dimer, trimer and two-segment models.</p>"
    page+=f"<p>{status.get('prepared_DFT_inputs',0)}/27 DFT inputs ready. Large-model memory qualification: {html.escape(status['memory_pilot'])}.</p>"
    page+='<p class="note">These are method controls, not synthesis-shortlist admissions. Prepared and xTB views are previews. Electronic ΔE is distinct from interaction ΔG and does not establish an equilibrium stacking population or poor polymer performance.</p>'
    page+='<p class="note"><strong>Design objective: weak interaction between redox sites and polymer segments.</strong> Strong persistent pi stacking and undesired intersite electron transfer are risks. More negative association energies are concerns to investigate, not benefits to maximize. This neutral pilot has not computed electron-transfer coupling or rates. It covers monomers, dimers, trimers and separate-segment pairs; no tetramers.</p>'
    page+='<h2>Rotate the molecular structures</h2><select id="model" aria-label="Choose a molecular model"></select><p id="stage"></p><div id="viewer"></div><p id="metrics"></p>'
    page+='<p>Blue, orange and green rings identify the redox-core atoms. Atom colors follow C/N/O/H identity; connecting lines show inferred connectivity, not computed bond orders. Pose labels describe the starting arrangement.</p>'
    page+='<div class="plots"><a href="figures/folding.pdf"><img src="figures/folding.png" alt="Matched covalent stacked versus open electronic energies"></a><a href="figures/association.pdf"><img src="figures/association.png" alt="Electronic association of two intact segments"></a></div>'
    page+='<script>const models='+json.dumps(models).replace('</','<\/')+r''';
const select=document.getElementById('model'),colors=['#0072B2','#D55E00','#009E73'];
models.forEach((m,i)=>{const o=document.createElement('option');o.value=i;o.textContent=m.label;select.append(o)});
function render(){if(!models.length){document.getElementById('stage').textContent='Conformer preparation pending';return}const m=models[+select.value],traces=[];
const elementColor={C:'#444444',H:'#b8bfc5',N:'#3366cc',O:'#d84040',F:'#009E73',Cl:'#009E73'};
let xs=[],ys=[],zs=[];m.bonds.forEach(b=>{xs.push(m.positions[b[0]][0],m.positions[b[1]][0],null);ys.push(m.positions[b[0]][1],m.positions[b[1]][1],null);zs.push(m.positions[b[0]][2],m.positions[b[1]][2],null)});
traces.push({type:'scatter3d',mode:'lines',x:xs,y:ys,z:zs,line:{color:'#929aa2',width:4},hoverinfo:'skip'});
[...new Set(m.elements)].forEach(e=>{let indices=m.elements.map((x,i)=>x===e?i:-1).filter(i=>i>=0);traces.push({type:'scatter3d',mode:'markers',name:e,x:indices.map(i=>m.positions[i][0]),y:indices.map(i=>m.positions[i][1]),z:indices.map(i=>m.positions[i][2]),marker:{size:e==='H'?3:5,color:elementColor[e]||'#999'},text:indices.map(i=>e+' '+(i+1)),hoverinfo:'text'})});
m.groups.forEach((g,k)=>{traces.push({type:'scatter3d',mode:'markers',name:'Redox core '+(k+1),x:g.map(i=>m.positions[i][0]),y:g.map(i=>m.positions[i][1]),z:g.map(i=>m.positions[i][2]),marker:{color:colors[k%3],symbol:'circle-open',size:9,line:{width:2}},hoverinfo:'name'})});
document.getElementById('stage').textContent='Stage: '+m.stage;
document.getElementById('metrics').textContent=m.metrics.map(p=>'Sites '+p.site_i+'–'+p.site_j+': centroid '+p.centroid_distance_A.toFixed(2)+' Å; normal projection '+p.plane_separation_A.toFixed(2)+' Å; angle '+p.interplane_angle_deg.toFixed(1)+'°; slip '+p.lateral_offset_A.toFixed(2)+' Å').join(' | ');
Plotly.react('viewer',traces,{margin:{l:0,r:0,t:0,b:0},showlegend:false,scene:{aspectmode:'data',xaxis:{visible:false},yaxis:{visible:false},zaxis:{visible:false},camera:{eye:{x:1.5,y:1.5,z:1.1}}}},{responsive:true,displaylogo:false});}
select.addEventListener('change',render);select.value=String(Math.max(0,models.findIndex(m=>m.label.includes('2 centers')&&m.label.endsWith('stacked'))));render();</script></html>'''
    temp=OUT/'index.html.tmp';temp.write_text(page);temp.replace(OUT/'index.html')
    subprocess.run([sys.executable,str(ROOT/'scripts/plotting/polymer_models/oligomer_summary.py')],check=True)
    manifest=dict(source_sha256=hashlib.sha256((OUT/'interaction_checked.json').read_bytes()).hexdigest(),
                  producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),viewable_models=len(models),
                  browser_sha256=hashlib.sha256((OUT/'index.html').read_bytes()).hexdigest(),
                  plotly_asset_sha256=hashlib.sha256(js.read_bytes()).hexdigest(),
                  files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in figdir.iterdir() if p.suffix in ('.pdf','.png','.svg')})
    (figdir/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Polymer interaction plots and',len(models),'rotatable model views refreshed')


if __name__=='__main__':main()
