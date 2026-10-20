#!/usr/bin/env python
"""Source-bound model construction, checked geometries and electronic results.

CPU-only; uses exact molecular graphs and saved DFT coordinates. Never fabricate
an open/stacked comparison, free energy, population, coupling or error bar.
"""
import hashlib
import html
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign
from rdkit.Chem.Draw import rdMolDraw2D

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.build.geometry import _skeleton
from redox.build.polymer_contacts import pi_core_groups, plane
from redox.build.oligomer_topology import chain_topology
from redox.workflows.checkpoints import verify_state_checkpoint

OUT=ROOT/'results/polymer_models'
FIG=OUT/'figures'
COLORS=['#0072B2','#D55E00','#009E73']
NAMES={'exp_phthalimide_control':'Phthalimide','design_pmdi_ethyl':'PMDI–ethyl','design_ndi_ethyl':'NDI–ethyl'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mapped_molecule(row):
    """Declared full graph with actual XYZ coordinates, independent of atom order."""
    mol=Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    observed=Chem.MolFromXYZFile(row['final_xyz'])
    rdDetermineBonds.DetermineConnectivity(observed,charge=0)
    mapping=_skeleton(observed).GetSubstructMatch(_skeleton(mol))
    if len(mapping)!=mol.GetNumAtoms():raise ValueError('XYZ/graph mapping failed: '+row['key'])
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,j in enumerate(mapping):conf.SetAtomPosition(i,observed.GetConformer().GetAtomPosition(j))
    mol.AddConformer(conf)
    return mol,mapping


def review(feed):
    rows={r['key']:r for r in feed['models']};checked={};sources={}
    for row in rows.values():
        if row['E_smd_eV'] is None:continue
        sd=Path(row['final_xyz']).parent;raw=json.loads((sd/'result.json').read_text())
        record=load_record(sd,raw['charge'],raw['mult'])
        assert record and raw['solvent']=='acetonitrile' and raw['charge']==0 and raw['mult']==1
        assert abs(record['e_smd_eV']-row['E_smd_eV'])<1e-9
        campaign=ROOT/(row.get('recovery_campaign') or feed['status']['campaign'])
        receipt=campaign/'done'/f"{row['key']}.json"
        report=verify_state_checkpoint(receipt)
        assert abs(report['E_smd_eV']-row['E_smd_eV'])<1e-9
        checked[row['key']]=row
        for p in [Path(row['final_xyz']),sd/'result.json',receipt]:sources[str(p)]=sha(p)
    folding=[]
    for comparison in feed['folding']:
        states={r['target']:r for r in rows.values() if r['model_key']==comparison['model_key']}
        a,b=states['open'],states['stacked']
        item=dict(comparison)
        if a['key'] in checked and b['key'] in checked:
            value=(b['E_smd_eV']-a['E_smd_eV'])*EV_KJ
            distinct=a['stacked_pairs']==0 and b['stacked_pairs']>0
            if comparison['stacked_minus_open_deltaE_kJ_mol'] is not None:
                assert distinct and abs(value-comparison['stacked_minus_open_deltaE_kJ_mol'])<1e-7
            # Same contact class is not proof of an identical conformer.
            ma,_=mapped_molecule(a);mb,_=mapped_molecule(b)
            rms=float(rdMolAlign.GetBestRMS(Chem.RemoveHs(ma),Chem.RemoveHs(mb),maxMatches=10000))
            item.update(open_start_final_stacked_pairs=a['stacked_pairs'],stacked_start_final_stacked_pairs=b['stacked_pairs'],
                        heavy_atom_best_RMSD_A=rms,distinct_open_stack_endpoints=distinct,
                        sampled_conformer_deltaE_kJ_mol=value,
                        interpretation='Matched open/stacked electronic screen' if distinct else 'Both starts have stacked contacts; no distinct open/stacked comparison')
        folding.append(item)
    monomers={r['candidate_id']:r for r in checked.values() if r['repeats']==1}
    associations=[]
    for comparison in feed['association']:
        item=dict(comparison)
        row=next(r for r in rows.values() if r['model_key']==comparison['model_key'])
        if row['key'] in checked and row['candidate_id'] in monomers:
            mono=monomers[row['candidate_id']]
            mol=Chem.AddHs(Chem.MolFromSmiles(row['smiles']));one=Chem.AddHs(Chem.MolFromSmiles(mono['smiles']))
            assert sorted(a.GetAtomicNum() for a in mol.GetAtoms())==sorted([a.GetAtomicNum() for a in one.GetAtoms()]*2)
            assert Chem.GetFormalCharge(mol)==2*Chem.GetFormalCharge(one)
            value=(row['E_smd_eV']-2*mono['E_smd_eV'])*EV_KJ
            assert abs(value-comparison['association_deltaE_kJ_mol'])<1e-7
            item.update(final_stacked_pairs=row['stacked_pairs'],contact_metrics=row['contact_metrics'],
                        interpretation='Attractive electronic screen; BSSE, association entropy and coupling/rates unqualified')
        associations.append(item)
    return dict(status='checked_electronic_screen',completed_structures=len(checked),
                folding=folding,association=associations,source_sha256=sources,
                limits=['Neutral molecular controls; not whole beads','Electronic deltaE, not deltaG or a stacking population',
                        'Association includes cap/backbone contacts and monomer deformation; not a pure redox-core pair energy',
                        'No validated rejection cutoff or electron-transfer coupling/rate'])


def save(fig,name):
    for ext in ('png','pdf','svg'):
        fig.savefig(FIG/f'{name}.{ext}',dpi=300,facecolor='white',bbox_inches='tight')
    plt.close(fig)


def construction(feed):
    ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
    root=ET.Element('{'+ns+'}svg',width='1500',height='1040',viewBox='0 0 1500 1040')
    ET.SubElement(root,'{'+ns+'}rect',width='1500',height='1040',fill='white')
    def label(x,y,text,size=22,color='#15304A',bold=False,anchor='start'):
        node=ET.SubElement(root,'{'+ns+'}text',x=str(x),y=str(y),fill=color,
                           style=f"font:{'bold ' if bold else ''}{size}px DejaVu Sans,sans-serif",attrib={'text-anchor':anchor})
        node.text=text
    label(30,40,'How the oligomers are modeled',30,bold=True)
    label(30,75,'Exact H-capped polystyrene fragments with pendant NDI groups; all depicted states are neutral.',19)
    for n,x in [(1,20),(2,515),(3,1010)]:
        row=next(r for r in feed['models'] if r['candidate_id']=='design_ndi_ethyl' and r['model_type']=='covalent_oligomer' and r['repeats']==n)
        mol=Chem.MolFromSmiles(row['smiles']);groups=pi_core_groups(mol,row['fragment'])
        topology=chain_topology(row['smiles'],n)
        atom_colors={atom:matplotlib.colors.to_rgb(COLORS[i]) for i,group in enumerate(groups) for atom in group}
        atom_colors.update({atom:matplotlib.colors.to_rgb('#6A4C93') for atom in topology['backbone_atoms']})
        bond_colors={b.GetIdx():atom_colors[b.GetBeginAtomIdx()] for b in mol.GetBonds() if b.GetBeginAtomIdx() in atom_colors and b.GetEndAtomIdx() in atom_colors and atom_colors[b.GetBeginAtomIdx()]==atom_colors[b.GetEndAtomIdx()]}
        drawer=rdMolDraw2D.MolDraw2DSVG(480,365);drawer.drawOptions().bondLineWidth=2.
        drawer.DrawMolecule(mol,highlightAtoms=list(atom_colors),highlightAtomColors=atom_colors,
                            highlightBonds=list(bond_colors),highlightBondColors=bond_colors)
        drawer.FinishDrawing();nested=ET.fromstring(drawer.GetDrawingText());nested.set('x',str(x));nested.set('y','150');root.append(nested)
        label(x+15,125,('A  Monomer' if n==1 else 'B  Dimer' if n==2 else 'C  Trimer'),24,bold=True)
        label(x+15,560,f'{n} repeat unit'+('' if n==1 else 's')+f' · {n} redox site'+('' if n==1 else 's'),20)
        label(x+15,595,f"{row['atoms']} atoms including H; H omitted in drawing",17,color='#556575')
    label(35,650,'Trimer backbone detail: one open linear chain, not a triangle (schematic):',23,bold=True)
    xs=[110+225*i for i in range(6)]
    for i,x in enumerate(xs):
        if i<5:ET.SubElement(root,'{'+ns+'}line',x1=str(x+35),y1='755',x2=str(xs[i+1]-35),y2='755',stroke='#6A4C93',attrib={'stroke-width':'5'})
        label(x,764,['CH₃','CH','CH₂','CH','CH₂','CH₂'][i],25,color='#6A4C93',bold=True,anchor='middle')
        label(x,805,'C'+str(i+1),17,color='#556575',anchor='middle')
        if i%2==1:
            ET.SubElement(root,'{'+ns+'}line',x1=str(x),y1='730',x2=str(x),y2='698',stroke=COLORS[i//2],attrib={'stroke-width':'4'})
            label(x,682,'Ar'+str(i//2+1),22,color=COLORS[i//2],bold=True,anchor='middle')
    label(35,855,'Arᵢ = p-C₆H₄–CH₂–NDI. Three pendant groups; six backbone carbons, five C–C bonds, two ends.',20)
    label(35,900,'Separate-segment model: two intact monomer segments, in parallel and antiparallel starting poses.',20)
    label(35,945,'DFT optimization: r2SCAN-D4/SVP/SMD(MeCN) → energy: ωB97M-V/TZVPD/SMD(MeCN).',20)
    label(35,990,'Purple: open polymer backbone. H caps are implicit. Current pilot uses 1–3 repeats; no whole-bead models.',19)
    path=FIG/'model_construction.svg';ET.ElementTree(root).write(path,encoding='unicode')
    import cairosvg
    cairosvg.svg2pdf(url=str(path),write_to=str(FIG/'model_construction.pdf'))
    cairosvg.svg2png(url=str(path),write_to=str(FIG/'model_construction.png'),output_width=3000)


def geometry_panel(ax,row,title):
    atoms=read(row['final_xyz']);xyz=atoms.positions
    center,normal,axis,_=plane(xyz[row['pi_groups'][0]])
    basis=np.array([axis,np.cross(normal,axis),normal]);xyz=(xyz-center)@basis.T
    mol,mapping=mapped_molecule(row)
    group_color={i:COLORS[k] for k,group in enumerate(row['pi_groups']) for i in group}
    for bond in mol.GetBonds():
        a,b=mapping[bond.GetBeginAtomIdx()],mapping[bond.GetEndAtomIdx()]
        if atoms.numbers[a]==1 or atoms.numbers[b]==1:continue
        color=group_color.get(a,'#9099A2') if group_color.get(a)==group_color.get(b) else '#9099A2'
        ax.plot(*xyz[[a,b]].T,color=color,lw=2.0 if a in group_color and b in group_color else 1.2,alpha=.95)
    heavy=[i for i,n in enumerate(atoms.numbers) if n!=1]
    ax.scatter(*xyz[heavy].T,c=[group_color.get(i,'#88939F') for i in heavy],s=12,depthshade=False)
    span=max(np.ptp(xyz[heavy],axis=0));mid=(xyz[heavy].max(axis=0)+xyz[heavy].min(axis=0))/2
    for setter,value in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],mid):setter(value-span*.55,value+span*.55)
    ax.set_box_aspect((1,1,1));ax.view_init(elev=20,azim=-55);ax.set_axis_off();ax.set_title(title,loc='left',fontsize=12,fontweight='bold',pad=0)
    contacts=[m for m in row['contact_metrics'] if m['stack_geometry']]
    text='; '.join(f"{m['plane_separation_A']:.2f} Å / {m['interplane_angle_deg']:.1f}°" for m in contacts)
    ax.text2D(.04,.02,'Core normal projection / angle:\n'+text,transform=ax.transAxes,fontsize=10)


def geometry_examples(feed):
    wanted=[('design_ndi_ethyl__n2__stereo0__open','A  NDI dimer: open start → stacked'),
            ('design_ndi_ethyl__n3__stereo1__stacked','B  NDI trimer: two stacked contacts'),
            ('exp_phthalimide_control__segments__parallel__parallel','C  Two separate phthalimide segments')]
    fig=plt.figure(figsize=(14.5,5.0));rows={r['key']:r for r in feed['models']}
    for i,(key,title) in enumerate(wanted):
        row=rows[key];assert row['E_smd_eV'] is not None and row['structural_review']['issues']==[]
        geometry_panel(fig.add_subplot(1,3,i+1,projection='3d'),row,title)
    fig.suptitle('Representative checked oligomer geometries',fontsize=17,fontweight='bold',y=.99)
    fig.text(.02,.01,'Actual saved coordinates; H hidden. Colors distinguish redox cores; gray is backbone/caps. Neutral SMD(MeCN) stationary structures; Hessian minima untested.',fontsize=9)
    fig.subplots_adjust(left=.01,right=.99,top=.89,bottom=.12,wspace=.02)
    save(fig,'checked_geometries')


def results_figure(feed,reviewed):
    fig=plt.figure(figsize=(12.5,6.6));gs=fig.add_gridspec(1,2,width_ratios=[1.5,1.],left=.22,right=.95,top=.77,bottom=.30,wspace=.7)
    ax=fig.add_subplot(gs[0,0]);assoc=[r for r in reviewed['association'] if r['association_deltaE_kJ_mol'] is not None]
    labels=[NAMES[r['candidate_id']]+'\n'+r['orientation']+' start' for r in assoc]
    vals=[r['association_deltaE_kJ_mol'] for r in assoc]
    ax.scatter(vals,range(len(vals)),color='#D55E00',s=75,zorder=3)
    ax.set_yticks(range(len(vals)),labels);ax.invert_yaxis();ax.axvline(0,color='.6',lw=1)
    for i,v in enumerate(vals):ax.annotate(f'{v:.1f}',(v,i),xytext=(0,11),textcoords='offset points',ha='center',fontsize=11)
    ax.set_xlim(-65,5);ax.set_ylim(len(vals)-.5,-.7);ax.set_xlabel('Electronic association ΔE (kJ/mol)');ax.set_title('A  Separate segments',loc='left',fontweight='bold')
    ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    ax=fig.add_subplot(gs[0,1]);fold=[r for r in reviewed['folding'] if r['stacked_minus_open_deltaE_kJ_mol'] is not None]
    if fold:
        vals=[r['stacked_minus_open_deltaE_kJ_mol'] for r in fold]
        ax.scatter(vals,range(len(vals)),color='#0072B2',s=75)
        ax.set_yticks(range(len(vals)),[NAMES[r['candidate_id']]+'\ntrimer, stereo '+str(r['stereoisomer']) for r in fold])
        for i,v in enumerate(vals):ax.annotate(f'{v:.1f}',(v,i),xytext=(0,13),textcoords='offset points',ha='center',fontsize=11)
        ax.set_xlim(min(-12,min(vals)-3),max(2,max(vals)+3));ax.set_ylim(-.8,len(vals)-.2)
    ax.axvline(0,color='.6',lw=1);ax.set_xlabel('Stacked − open electronic ΔE (kJ/mol)');ax.set_title('B  Covalent folding',loc='left',fontweight='bold')
    fig.suptitle('Quantifying nearby-site interaction strength',fontsize=16,fontweight='bold',y=.98)
    status=feed['status'];fig.text(.04,.88,f"{status['stationary_electronic_cases']}/{status['expected_DFT_cases']} DFT cases checked · MeCN · no interaction ΔG or electron-transfer rates",fontsize=12)
    dimers=[r for r in reviewed['folding'] if r['repeats']==2]
    both=sum(r.get('open_start_final_stacked_pairs',0)>0 and r.get('stacked_start_final_stacked_pairs',0)>0 for r in dimers)
    fig.text(.04,.18,f'Stacked contacts occur for {both}/{len(dimers)} dimer controls. Geometry alone does not quantify interaction strength.',fontsize=11,fontweight='bold')
    fig.text(.04,.11,'Association includes backbone/cap contacts and deformation; covalent folding includes conformational effects. BSSE and entropy remain unqualified.',fontsize=10)
    fig.text(.04,.06,'These deterministic structure comparisons have no quantified error bars. Geometry and electronic ΔE do not establish equilibrium populations or transfer rates.',fontsize=9)
    save(fig,'key_results')


def decomposition_figure():
    path=OUT/'interaction_decomposition.json'
    if not path.exists():return None
    feed=json.loads(path.read_text());rows=[r for r in feed['pairs'] if r['status']=='checked_electronic_decomposition']
    if not rows:return None
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    labels=[NAMES[r['candidate_id']]+' / '+r['orientation']+' start' for r in rows]
    for ax,key,title,color in [(axes[0],'frozen_interaction_deltaE_kJ_mol','A  Frozen full-segment interaction','#0072B2'),
                               (axes[1],'segment_deformation_deltaE_kJ_mol','B  Segment deformation','#D55E00')]:
        values=[r[key] for r in rows];ax.scatter(values,range(len(rows)),s=65,color=color,zorder=3)
        ax.set_yticks(range(len(rows)),labels if ax is axes[0] else ['']*len(rows));ax.invert_yaxis();ax.axvline(0,color='.7',lw=1)
        for i,v in enumerate(values):ax.annotate(f'{v:.2f}',(v,i),xytext=(0,12),textcoords='offset points',ha='center',fontsize=11)
        ax.set_title(title,loc='left',fontweight='bold');ax.set_xlabel('Electronic ΔE (kJ/mol)');ax.grid(axis='x',alpha=.15);ax.set_ylim(len(rows)-.5,-.7)
    axes[0].set_xlim(-65,5);axes[1].set_xlim(-.05,max(.8,max(r['segment_deformation_deltaE_kJ_mol'] for r in rows)*1.4))
    fig.suptitle('Interaction and deformation are different contributions',fontsize=16,fontweight='bold')
    fig.supxlabel('Association ΔE = frozen interaction + deformation. Neutral full segments in SMD(MeCN); no BSSE correction, free energy or pure-core attribution.',fontsize=9)
    save(fig,'fragment_decomposition')
    return feed


def main():
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    inherited=sorted(os.sched_getaffinity(0))
    if len(inherited)<=2:
        # A collector child reuses its parent's reserved worker CPU allocation.
        os.sched_setaffinity(0,{inherited[0]})
    else:
        with allocation_guard():
            host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    FIG.mkdir(exist_ok=True);source=OUT/'interaction_checked.json';source_hash=sha(source)
    subprocess.run([sys.executable,str(ROOT/'scripts/validation/check_oligomer_topology.py'),'--out',str(OUT/'model_topology.json')],check=True)
    feed=json.loads(source.read_text());reviewed=review(feed)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    construction(feed);geometry_examples(feed);results_figure(feed,reviewed);decomposition=decomposition_figure()
    assert sha(source)==source_hash,'Input results changed while figures were built'
    assert all(sha(p)==h for p,h in reviewed['source_sha256'].items()),'Evidence changed during review'
    atomic_write_json(OUT/'screening_review.json',dict(reviewed,input_sha256=source_hash))
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Oligomer models and key results</title>
<style>body{font:17px/1.55 system-ui,sans-serif;max-width:1250px;margin:35px auto;padding:0 20px;color:#15304a}img{width:100%;height:auto}a{color:#0072b2}.note{padding:16px;background:#fff4df}</style>
<h1>Oligomer models and key results</h1><p><a href="index.html">Rotate all 3D models</a> · <a href="README.md">Pipeline and scope</a> · <a href="screening_review.json">Checked numerical comparisons</a></p>
<p class="note"><b>Design objective: weak intersite and segment interactions.</b> Strong persistent stacking and undesired electron transfer are concerns. Electronic attraction alone does not establish a stacking population or electron-transfer rate.</p>'''
    for name,title in [('model_construction','1. Exact model construction'),('checked_geometries','2. Checked DFT geometries'),('key_results','3. Completed electronic comparisons')]:
        page+=f'<h2>{html.escape(title)}</h2><a href="figures/{name}.pdf"><img src="figures/{name}.png" alt="{html.escape(title)}"></a><p><a href="figures/{name}.pdf">Vector PDF</a> · <a href="figures/{name}.svg">SVG</a></p>'
    if decomposition:
        page+='<h2>4. Interaction strength separated from deformation</h2><a href="figures/fragment_decomposition.pdf"><img src="figures/fragment_decomposition.png" alt="Frozen full-segment interaction and deformation energies"></a><p><a href="interaction_decomposition.json">Verified values and definitions</a> · <a href="figures/fragment_decomposition.pdf">Vector PDF</a></p>'
    page+='<p>Models are neutral H-capped polymer fragments. The separate-segment complexes contain two intact molecules; none depicts whole beads. No tetramers are included. Recovering calculations remain unqualified until the normal QC passes.</p></html>'
    (OUT/'summary.html').write_text(page)
    files=[FIG/(name+'.'+ext) for name in ['model_construction','checked_geometries','key_results'] for ext in ['pdf','png','svg']]+[OUT/'summary.html',OUT/'screening_review.json']
    if decomposition:files += [FIG/('fragment_decomposition.'+ext) for ext in ['pdf','png','svg']]
    atomic_write_json(FIG/'summary_manifest.json',dict(input_sha256=source_hash,producer_sha256=sha(Path(__file__)),
        topology_review_sha256=sha(OUT/'model_topology.json'),topology_checker_sha256=sha(ROOT/'src/redox/build/oligomer_topology.py'),
        topology_XYZ_checker_sha256=sha(ROOT/'scripts/validation/check_oligomer_topology.py'),constructor_sha256=sha(ROOT/'src/redox/build/oligomers.py'),
        decomposition_source_sha256=sha(OUT/'interaction_decomposition.json') if decomposition else None,
        source_sha256=reviewed['source_sha256'],outputs={str(p.relative_to(ROOT)):sha(p) for p in files}))
    print(json.dumps(dict(checked_structures=reviewed['completed_structures'],completed_associations=sum(r['association_deltaE_kJ_mol'] is not None for r in reviewed['association']),completed_folding=sum(r['stacked_minus_open_deltaE_kJ_mol'] is not None for r in reviewed['folding']),figures=3)))


if __name__=='__main__':main()
