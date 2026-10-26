#!/usr/bin/env python
"""Draw the checked polymer-to-local-model cut from actual fixed coordinates."""
import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as effects
from matplotlib.lines import Line2D
import numpy as np
from rdkit import Chem

from redox.build.polymer_contacts import pi_core_groups,plane
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT/'results/figures/polymer_models/local_motif_construction');a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    source=ROOT/'results/figures/polymer_models/spacing_conformers.json'
    prepared=ROOT/'results/polymer_models/campaigns/2026-10-26_local_motif_fragments/manifest.json'
    reviewed=ROOT/'results/reviews/2026-10-26/local_motif_construction/summary.json'
    review=json.loads(reviewed.read_text());api.verify_sources(review['source_sha256'])
    row=next(r for r in json.loads(prepared.read_text())['entries'] if r['key']=='occ_2_3__pose0__cap109')
    original=next(r for r in json.loads(source.read_text())['states'] if r['key']==row['parent_key'])
    mol=Chem.AddHs(Chem.MolFromSmiles(original['task']['smiles']));coords=np.asarray(original['coordinates_A'])
    groups=pi_core_groups(mol,original['task']['fragment']);centers=[];normals=[]
    for group in groups:
        center,normal,_,_=plane(coords[group]);centers.append(center);normals.append(normal)
    if normals[0]@normals[1]<0:normals[1]*=-1
    z=normals[0]+normals[1];z/=np.linalg.norm(z);delta=centers[1]-centers[0];x=delta-z*(delta@z)
    if np.linalg.norm(x)<1e-6:x=coords[groups[0][1]]-coords[groups[0][0]];x-=z*(x@z)
    x/=np.linalg.norm(x);y=np.cross(z,x);theta=np.pi/4
    rotation=np.column_stack([x,np.cos(theta)*y+np.sin(theta)*z,-np.sin(theta)*y+np.cos(theta)*z])
    if abs(np.linalg.det(rotation)-1)>1e-10:raise ValueError('The camera must be a proper rotation')
    center=np.mean(centers,axis=0);colors=['#0072B2','#CC79A7'];parent_colors={}
    for color,part in zip(colors,row['fragments']):
        for i in part['parent_indices']:parent_colors[i]=color
    pair_item=next(m for m in row['molecules'] if m['id']=='AB');pair=Chem.MolFromMolFile(pair_item['molecular_graph'],removeHs=False)
    data=[line.split() for line in Path(pair_item['geometry']).read_text().splitlines()[2:]]
    pair_coords=np.asarray([[float(v) for v in line[1:]] for line in data]);offset=len(row['fragments'][0]['parent_indices'])+1
    cap_indices=[row['fragments'][0]['cap_index'],row['fragments'][1]['cap_index']+offset]
    pair_colors={i:colors[0 if i<offset else 1] for i in range(pair.GetNumAtoms())}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    figure,axes=plt.subplots(1,2,figsize=(7.2,3.45),gridspec_kw={'width_ratios':[1.3,1]})
    for axis,graph,xyz,color_map,caps,title in zip(axes,[mol,pair],[coords,pair_coords],[parent_colors,pair_colors],
            [[],cap_indices],['a  Five-repeat phthalimide control','b  Two capped motifs (zoom)']):
        view=(xyz-center)@rotation;visible=[atom.GetIdx() for atom in graph.GetAtoms() if atom.GetAtomicNum()!=1 or atom.GetIdx() in caps]
        for bond in sorted(graph.GetBonds(),key=lambda b:np.mean(view[[b.GetBeginAtomIdx(),b.GetEndAtomIdx()],2])):
            i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            if i not in visible or j not in visible:continue
            color='#D55E00' if i in caps or j in caps else color_map.get(i,'#a0a0a0') if color_map.get(i)==color_map.get(j) else '#888888'
            axis.plot(view[[i,j],0],view[[i,j],1],color=color,lw=1.2,solid_capstyle='round',
                path_effects=[effects.Stroke(linewidth=3.2,foreground='white'),effects.Normal()])
        for atom in graph.GetAtoms():
            i=atom.GetIdx()
            if i not in visible:continue
            if atom.GetAtomicNum()==6:
                axis.scatter(*view[i,:2],s=7,color=color_map.get(i,'#a0a0a0'),linewidths=0,zorder=3)
            if atom.GetAtomicNum() not in (1,6) or i in caps:
                axis.text(*view[i,:2],atom.GetSymbol(),fontsize=7 if i not in caps else 9,ha='center',va='center',
                    color='#D55E00' if i in caps else color_map.get(i,'#555555'),
                    path_effects=[effects.withStroke(linewidth=2,foreground='white')])
        if graph is mol:
            for part in row['fragments']:
                i,j=part['cut_parent_bond'];axis.plot(view[[i,j],0],view[[i,j],1],color='#D55E00',lw=2.2,ls='--')
        extent=view[visible,:2];span=np.ptp(extent,axis=0);margin=max(span)*.12
        axis.set_xlim(extent[:,0].min()-margin,extent[:,0].max()+margin);axis.set_ylim(extent[:,1].min()-margin,extent[:,1].max()+margin)
        axis.set_aspect('equal');axis.axis('off');axis.set_title(title,loc='left',fontsize=8)
    figure.legend(handles=[Line2D([],[],color=colors[0],label='Motif A'),Line2D([],[],color=colors[1],label='Motif B'),
        Line2D([],[],color='#a0a0a0',label='Omitted polymer'),Line2D([],[],color='#D55E00',ls='--',label='Cut / H cap')],
        loc='lower center',ncol=4,frameon=False,fontsize=7)
    figure.subplots_adjust(left=.025,right=.99,top=.91,bottom=.16,wspace=.15)
    out=a.out.resolve();out.parent.mkdir(parents=True,exist_ok=True)
    for suffix in ('pdf','svg','png'):figure.savefig(out.with_suffix('.'+suffix),dpi=220,facecolor='white')
    plt.close(figure)
    record=dict(status='checked_frozen_local_construction_rendered',parent_key=row['parent_key'],cap_length_A=1.09,
        proper_rotation=rotation.tolist(),rotation_determinant=float(np.linalg.det(rotation)),parent_atoms=mol.GetNumAtoms(),local_pair_atoms=pair.GetNumAtoms(),
        core_centroid_distance_A=row['contact']['centroid_distance_A'],relative_positions_changed=False,
        caption='Coordinate wireframe of the complete phthalimide motifs and retained benzyl carbons. Replacing each polymer-aryl bond by a 1.09-Angstrom H cap gives two disconnected N-methylphthalimide fragments at unchanged parent positions; the right panel is enlarged. White depth outlines distinguish crossing bonds. Bond orders are retained in the source MOL files. This construction figure does not establish an interaction free energy or coupling.',
        source_sha256={str(p):api.digest(p) for p in [Path(__file__).resolve(),source,prepared,reviewed,Path(pair_item['geometry']),Path(pair_item['molecular_graph'])]},
        outputs_sha256={str(out.with_suffix('.'+s)):api.digest(out.with_suffix('.'+s)) for s in ('pdf','svg','png')})
    atomic_write_json(out.with_suffix('.json'),record)


if __name__=='__main__':main()
