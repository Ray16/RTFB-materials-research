#!/usr/bin/env python
"""Render the actual native PMDI pair with both core planes visible."""
import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from rdkit import Chem

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cached',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    if any(args.out.with_suffix('.'+e).exists() for e in ('pdf','png','json')):
        raise ValueError('Preserve existing figure outputs')
    helper_path=ROOT/'scripts/plotting/polymer_models/plot_neighbor_grafting.py'
    helper=api.load_module('neighbor_geometry_checked_coordinate_cache',helper_path)
    data=json.loads(args.cached.read_text());helper.validate_cached(data)
    geo=data['geometry'];mol=Chem.AddHs(Chem.MolFromSmiles(geo['smiles']))
    xyz=np.array(geo['coordinates_A']);heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()!=1]
    coords=(xyz-xyz[heavy].mean(axis=0))@np.array(geo['rotation'])
    core={i:color for g,color in zip(geo['groups'],('#D55E00','#CC79A7')) for i in g};backbone=set(geo['backbone'])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
    figure,axis=plt.subplots(figsize=(6.8,4.2))
    for bond in mol.GetBonds():
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if i not in heavy or j not in heavy:continue
        color,width=('#0072B2',2.2) if i in backbone and j in backbone else (core[i],1.8) if i in core and core.get(j)==core[i] else ('#A7A7A7',1.)
        axis.plot(*coords[[i,j],:2].T,color=color,lw=width,zorder=2)
    axis.scatter(*coords[heavy,:2].T,s=13,c=[core.get(i,'#0072B2' if i in backbone else '#009E73' if mol.GetAtomWithIdx(i).GetSymbol()=='Cl' else '#A7A7A7') for i in heavy],zorder=3)
    centers=np.array([coords[g].mean(axis=0) for g in geo['groups']])
    axis.plot(*centers[:,:2].T,'--',color='#333333',lw=1.,zorder=4)
    for site in geo['occupied_sites']:
        point=coords[geo['backbone'][2*site-1]]
        offset=(18,-20) if site==2 else (-22,-14)
        axis.annotate(str(site),xy=point[:2],xytext=offset,textcoords='offset points',color='#0072B2',
            fontsize=10,fontweight='bold',arrowprops=dict(arrowstyle='-',color='#0072B2',lw=.6),
            bbox=dict(facecolor='white',edgecolor='none',pad=.4))
    axis.set_aspect('equal');axis.axis('off');axis.margins(.12,.18)
    handles=[Line2D([0],[0],color='#0072B2',lw=2,label='Open backbone'),
        Line2D([0],[0],color='#D55E00',lw=2,label='PMDI core 1'),
        Line2D([0],[0],color='#CC79A7',lw=2,label='PMDI core 2')]
    figure.legend(handles=handles,loc='upper center',ncol=3,frameon=False,fontsize=9)
    figure.text(.5,.07,f"3D core-centroid distance: {geo['contact']['centroid_distance_A']:.2f} Å; plane angle: {geo['contact']['interplane_angle_deg']:.1f}°",ha='center',fontsize=9)
    figure.subplots_adjust(left=.04,right=.96,bottom=.15,top=.88)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','png'):figure.savefig(args.out.with_suffix('.'+ext),dpi=600,facecolor='white')
    plt.close(figure)
    atomic_write_json(args.out.with_suffix('.json'),dict(status='checked_PMDI_pair_coordinate_cache_rendered',
        source_sha256={str(p.resolve()):api.digest(p) for p in (args.cached,Path(__file__),helper_path)},
        geometry=geo,proper_rotation_only=True,hydrogens_hidden=True,
        isolated_motif_binding_or_coupling_inferred=False,conformer_convergence_established=False,
        output_sha256={str(args.out.with_suffix('.'+ext)):api.digest(args.out.with_suffix('.'+ext)) for ext in ('pdf','png')}))


if __name__=='__main__':main()
