#!/usr/bin/env python
"""Plot conditional grafting and cache its exact matched PMDI geometry."""
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

from redox.build.polymer_cases import _path
from redox.build.polymer_contacts import pi_core_groups, contact_metrics
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def geometry(report):
    task=report['states']['pair']['task'];xyz=Path(report['states']['pair']['state_dir'])/'opt.xyz'
    renderer_path=ROOT/'scripts/plotting/polymer_models/oligomer_summary.py'
    renderer=api.load_module('neighbor_pair_saved_molecular_mapping',renderer_path)
    mol,mapping=renderer.mapped_molecule(dict(task,final_xyz=str(xyz)))
    groups=pi_core_groups(mol,task['fragment'])
    if {frozenset(mapping[i] for i in g) for g in groups}!={frozenset(g) for g in task['pi_groups']}:
        raise ValueError('Conjugated cores do not map to the original XYZ atoms')
    mapped=Chem.MolFromSmiles(task['mapped_smiles']);canonical=Chem.MolFromSmiles(task['smiles'])
    match=canonical.GetSubstructMatch(mapped,useChirality=True)
    if len(match)!=mapped.GetNumAtoms():raise ValueError('Backbone stereochemical mapping failed')
    backbone=[match[i] for _,i in sorted(_path(mapped,task['repeat_units']).items())]
    coords=np.asarray(mol.GetConformer().GetPositions());heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()!=1]
    metric=contact_metrics(coords,groups)[0]
    original=report['paired_component_QC']['geometry_audit']['contact_metrics'][0]
    if any(abs(metric[k]-original[k])>1e-8 for k in ('centroid_distance_A','interplane_angle_deg','plane_separation_A')):
        raise ValueError('Drawing changed the independently checked pair geometry')
    normals=[]
    for group in groups:
        centered=coords[group]-coords[group].mean(axis=0)
        _,_,vectors=np.linalg.svd(centered,full_matrices=False);normals.append(vectors[-1])
    if np.dot(normals[0],normals[1])<0:normals[1]*=-1
    vertical=normals[0]+normals[1];vertical/=np.linalg.norm(vertical)
    horizontal=coords[groups[1]].mean(axis=0)-coords[groups[0]].mean(axis=0)
    horizontal-=vertical*np.dot(horizontal,vertical);horizontal/=np.linalg.norm(horizontal)
    upward=np.cross(vertical,horizontal)
    rotation=np.stack([horizontal,upward,vertical])
    return dict(smiles=task['smiles'],fragment=task['fragment'],occupied_sites=task['occupied_sites'],
        atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],coordinates_A=coords.tolist(),
        rotation=rotation.T.tolist(),groups=groups,backbone=backbone,contact=metric),{
            str(xyz):api.digest(xyz),str(renderer_path):api.digest(renderer_path)}


def validate_cached(data):
    if data['status']!='matched_neighbor_grafting_and_saved_geometry_rendered':
        raise ValueError('A checked finite-cycle figure cache is required')
    geo=data['geometry'];mol=Chem.AddHs(Chem.MolFromSmiles(geo['smiles']))
    if [a.GetAtomicNum() for a in mol.GetAtoms()]!=geo['atomic_numbers']:
        raise ValueError('Cached graph and coordinate atom order differ')
    coords=np.array(geo['coordinates_A']);rotation=np.array(geo['rotation'])
    if (coords.shape!=(mol.GetNumAtoms(),3) or not np.isfinite(coords).all()
            or not np.allclose(rotation@rotation.T,np.eye(3),rtol=0,atol=1e-10)
            or abs(np.linalg.det(rotation)-1)>1e-10):
        raise ValueError('Cached view must preserve finite coordinates, handedness and scale')
    groups=pi_core_groups(mol,geo['fragment']);metric=contact_metrics(coords,groups)[0]
    if ({frozenset(g) for g in groups}!={frozenset(g) for g in geo['groups']}
            or abs(metric['centroid_distance_A']-geo['contact']['centroid_distance_A'])>1e-8):
        raise ValueError('Cached conjugated cores or centroid distance differ')
    for cutoff,row in data['by_cutoff'].items():
        g=row['state_G_eV'];expected=(g['pair']+g['empty']-g['site2']-g['site3'])*data['conversion_kJ_mol_per_eV']
        if abs(expected-row['deltaDeltaG_kJ_mol'])>1e-7:raise ValueError('Cached neighboring-graft arithmetic differs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--review',type=Path);group.add_argument('--cached',type=Path)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    if any(args.out.with_suffix('.'+ext).exists() for ext in ('pdf','png','json')):
        raise ValueError('Preserve existing checked figures; choose a fresh output')
    if args.cached:
        data=json.loads(args.cached.read_text());validate_cached(data)
        sources={str(args.cached.resolve()):api.digest(args.cached)}
    else:
        report=json.loads(args.review.read_text());api.verify_sources(report['source_sha256'])
        if (report['status']!='exact_candidate_neighbor_grafting_cycle_independently_reconstructed'
                or report['conformer_convergence_established'] or report['synthesis_admission']):
            raise ValueError('A finite, independently reconstructed neighbor cycle is required')
        geo,sources=geometry(report);sources[str(args.review.resolve())]=api.digest(args.review)
        data=dict(status='matched_neighbor_grafting_and_saved_geometry_rendered',geometry=geo,
            by_cutoff=report['by_cutoff'],conversion_kJ_mol_per_eV=report['conversion_kJ_mol_per_eV'],
            primary_cutoff_cm=50,conformer_convergence_established=False,electronic_coupling_computed=False)
    sources[str(Path(__file__).resolve())]=api.digest(__file__)
    for name in ('polymer_contacts.py','polymer_cases.py'):
        p=ROOT/'src/redox/build'/name;sources[str(p)]=api.digest(p)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
    figure,(a,b)=plt.subplots(1,2,figsize=(7.2,2.7))
    primary=data['by_cutoff']['50']['models']['free_ions']
    for index,site in enumerate((2,3)):
        first=primary[f'graft_{site}_on_empty_kJ_mol'];conditional=primary[f'graft_{site}_given_{5-site}_kJ_mol']
        a.plot([index,index],[first,conditional],color='#666666',linewidth=1,zorder=1)
        a.scatter(index,first,c='#0072B2',marker='o',s=36,zorder=3)
        a.scatter(index,conditional,c='#D55E00',marker='s',s=36,zorder=3)
    a.set(xlim=(-.5,1.5),ylim=(-102,-94),xticks=[0,1],xticklabels=['Attach at site 2','Attach at site 3'],
        ylabel='Grafting ΔG° / kJ mol⁻¹')
    a.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#0072B2',markeredgecolor='#0072B2',label='Other site empty'),
        Line2D([0],[0],marker='s',color='none',markerfacecolor='#D55E00',markeredgecolor='#D55E00',label='Other site grafted')],
        loc='upper right',frameon=False,fontsize=7.5)
    a.set_title('a  Conditional grafting',loc='left',fontweight='bold',fontsize=9)
    cutoffs=[25,50,100];delta=[data['by_cutoff'][str(k)]['deltaDeltaG_kJ_mol'] for k in cutoffs]
    b.axhline(0,color='#777777',linewidth=.8,linestyle='--')
    b.plot(cutoffs,delta,'o-',color='#8064A2',linewidth=1.4,markersize=5)
    b.set(xticks=cutoffs,xlabel='qRRHO rotor cutoff / cm⁻¹',ylabel='ΔΔG / kJ mol⁻¹',ylim=(-2.5,.3))
    b.set_title('b  Rotor-cutoff sensitivity',loc='left',fontweight='bold',fontsize=9)
    for ax in (a,b):ax.spines[['top','right']].set_visible(False)
    figure.subplots_adjust(left=.1,right=.98,top=.88,bottom=.22,wspace=.45)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','png'):figure.savefig(args.out.with_suffix('.'+ext),dpi=600,facecolor='white')
    plt.close(figure);api.verify_sources(sources)
    data.update(source_sha256=sources,cached_replay=bool(args.cached),
        caption='Finite one-minimum PMDI conditional grafting in MeCN. Free-ion DBU reactions shown; common precursor/base/byproduct terms cancel from ΔΔG. The proper-rotation orthographic view uses the checked XYZ and covalent graph, with hydrogens omitted. Blue marks the open backbone, orange/magenta the diimide cores and green the remaining chloride. Numbers label occupied repeats. No isolated core binding, electronic coupling, ensemble convergence or experimental occupancy is inferred.',
        output_sha256={str(args.out.with_suffix('.'+ext)):api.digest(args.out.with_suffix('.'+ext)) for ext in ('pdf','png')})
    atomic_write_json(args.out.with_suffix('.json'),data)


if __name__=='__main__':main()
