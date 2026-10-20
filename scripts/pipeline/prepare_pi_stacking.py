#!/usr/bin/env python
"""Prepare finite, topology-checked MeCN xTB ensembles for polymer method controls."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
from rdkit import Chem

from redox.build.oligomers import decorate_sites
from redox.build.polymer_contacts import starting_ensemble, enantiomer_representatives, contact_metrics, plane
from redox.build.structures import to_xyz
from redox.core.common import ROOT, HARTREE_EV
from redox.core.io import atomic_write_json
from redox.validation.protonation import hydrogen_topology_status

BASE=ROOT/'results/polymer_models/campaigns/2026-10-19_pi_stacking'
XTB=Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb')
CONTROLS=('exp_phthalimide_control','design_pmdi_ethyl','design_ndi_ethyl')


def xtb(directory,smiles,groups):
    command=[str(XTB),'input.xyz','--gfn','2','--alpb','acetonitrile','--chrg','0','--uhf','0','--opt','tight','--iterations','1000']
    with (directory/'xtb.log').open('w') as log:
        proc=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=600)
    text=(directory/'xtb.log').read_text()
    values=re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh',text)
    output=directory/'xtbopt.xyz'
    if proc.returncode or not output.exists() or not values or 'GEOMETRY OPTIMIZATION CONVERGED' not in text:
        raise RuntimeError('xTB did not prove a converged geometry')
    topology=hydrogen_topology_status(smiles,output)
    if topology!='ok':raise RuntimeError('Declared graph changed: '+topology)
    from ase.io import read
    atoms=read(output)
    metrics=contact_metrics(atoms.positions,groups)
    return dict(method='GFN2-xTB/ALPB(MeCN); exploratory electronic energy',
                E_eV=float(values[-1])*HARTREE_EV,xyz=str(output),metrics=metrics,
                stacked_pairs=sum(m['stack_geometry'] for m in metrics),topology=topology,
                graph_and_optimization_checked=True,DFT_qualified=False)


def sample_case(case):
    allowed=sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0,{allowed[mp.current_process()._identity[0]%len(allowed)]})
    directory=BASE/'seeds'/case['key']
    if (directory/'ensemble.json').exists():
        existing=json.loads((directory/'ensemble.json').read_text())
        assert existing['case']==case,'Existing exact graph changed'
        assert all(Path(r['xyz']).exists() for r in existing['records'] if 'E_eV' in r)
        return existing
    directory.mkdir(parents=True,exist_ok=False)
    records=[];failures=[]
    modes=('unbiased',) if case['repeats']==1 else ('unbiased','parallel','slipped','antiparallel','open')
    for mode in modes:
        try:mol,groups,ranked=starting_ensemble(case['smiles'],case['fragment'],mode,count=6)
        except Exception as exc:
            failures.append(dict(mode=mode,stage='embedding',error=repr(exc)));continue
        assert len(groups)==case['repeats']
        for i,(_,cid) in enumerate(ranked[:2]):
            seed=directory/f'{mode}_{i}';seed.mkdir()
            one=Chem.Mol(mol);one.RemoveAllConformers();one.AddConformer(mol.GetConformer(cid),assignId=True)
            (seed/'input.xyz').write_text(to_xyz(one,'biased FF starting geometry; constraints released for xTB'))
            initial=contact_metrics(one.GetConformer().GetPositions(),groups)
            item=dict(mode=mode,seed=i,pi_groups=groups,initial_metrics=initial,
                      initial_xyz=str(seed/'input.xyz'),FF_energy_is_not_an_interaction_energy=True)
            try:item.update(xtb(seed,case['smiles'],groups))
            except Exception as exc:
                item.update(error=repr(exc),DFT_qualified=False)
                failures.append(dict(mode=mode,seed=i,stage='xTB',error=repr(exc)))
            atomic_write_json(seed/'record.json',item);records.append(item)
    if not records:
        atomic_write_json(directory/'preparation_failure.json',dict(case=case,failures=failures))
        raise RuntimeError('No starting geometries: '+case['key'])
    result=dict(case=case,records=records,failures=failures)
    atomic_write_json(directory/'ensemble.json',result)
    return result


def choose(ensemble,target):
    records=ensemble['records']
    successful=[r for r in records if 'E_eV' in r]
    if target=='monomer':
        assert successful,'No qualified xTB warm start for monomer'
        return min(successful,key=lambda r:r['E_eV']),False
    filtered=[r for r in successful if (r['stacked_pairs']>0 if target=='stacked' else r['stacked_pairs']==0)]
    if filtered:return min(filtered,key=lambda r:r['E_eV']),False
    # A geometrically distinct FF seed may still be tried by DFT. Its FF energy
    # is never presented as a relaxed minimum or an interaction result.
    candidates=[r for r in records if any(m['stack_geometry'] for m in r['initial_metrics'])== (target=='stacked')]
    if not candidates:raise ValueError('No distinct '+target+' starting geometry for '+ensemble['case']['key'])
    chosen=candidates[0]
    return dict(chosen,xyz=chosen['initial_xyz']),True


def add_dft_task(tasks,case,record,target,fallback):
    key=case['key']+'__'+target
    directory=BASE/'inputs'/key;directory.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(record['xyz'],directory/'input.xyz')
    count=int((directory/'input.xyz').read_text().splitlines()[0])
    tasks.append(dict(key=key,candidate_id=case['candidate_id'],model_key=case['key'],
        model_type=case['model_type'],repeats=case['repeats'],stereoisomer=case['stereoisomer'],
        target=target,smiles=case['smiles'],fragment=case['fragment'],pi_groups=record['pi_groups'],
        atoms=count,charge=0,mult=1,input_xyz=str(directory/'input.xyz'),
        seed_stage='biased_FF; DFT_relaxation_required' if fallback else 'converged_topology_checked_xTB_ALPB',
        role='interaction_method_control; not synthesis admission'))


def segment_pair(case,record,orientation):
    from ase.io import read,write
    atoms=read(record['xyz']);positions=atoms.positions
    group=record['pi_groups'][0];center,normal,axis,_=plane(positions[group])
    second_axis=np.cross(normal,axis);basis=np.array([axis,second_axis,normal])
    local=(positions-center)@basis.T
    other=local.copy()
    if orientation=='antiparallel':other[:,:2]*=-1
    other+=np.array([2.,0.,3.6])
    # Rotate both segments back into the same reference frame.
    first_xyz=local@basis+center;second_xyz=other@basis+center
    n=len(atoms);pair=atoms.copy();pair.positions=first_xyz
    other_atoms=atoms.copy();other_atoms.positions=second_xyz
    pair+=other_atoms
    distance=np.linalg.norm(first_xyz[:,None,:]-second_xyz[None,:,:],axis=2).min()
    if distance<.9:raise ValueError('Initial segment pair has an atomic overlap')
    key=case['candidate_id']+'__segments__'+orientation
    d=BASE/'seeds'/key;d.mkdir(parents=True,exist_ok=False)
    write(d/'input.xyz',pair,format='xyz')
    groups=[group,[i+n for i in group]]
    smiles=case['smiles']+'.'+case['smiles']
    entry=dict(mode=orientation,pi_groups=groups,initial_xyz=str(d/'input.xyz'),
               initial_metrics=contact_metrics(pair.positions,groups))
    try:entry.update(xtb(d,smiles,groups));fallback=False
    except Exception as exc:entry.update(xyz=str(d/'input.xyz'),error=repr(exc));fallback=True
    atomic_write_json(d/'record.json',entry)
    description=dict(key=key,candidate_id=case['candidate_id'],smiles=smiles,
                     fragment=case['fragment'],repeats=2,stereoisomer='two achiral segments',model_type='two_nonbonded_segments')
    return description,entry,fallback


def main():
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    assert not (BASE/'preparation.json').exists(),'Preserve the completed preparation'
    BASE.mkdir(parents=True,exist_ok=True)
    feed=json.loads((ROOT/'results/candidates/current.json').read_text())
    rows={r['id']:r for r in feed['original39']+feed['additional_candidates']}
    cases=[]
    for key in CONTROLS:
        row=rows[key]
        frag=Chem.MolFromSmiles(row['frag'])
        assert all(a.GetChiralTag()==Chem.ChiralType.CHI_UNSPECIFIED for a in frag.GetAtoms()),'Enantiomer reduction needs achiral pendants'
        for n in (1,2,3):
            variants=enantiomer_representatives(decorate_sites(row['frag'],n,range(1,n+1)))
            for i,smi in enumerate(variants):
                cases.append(dict(key=key+f'__n{n}__stereo{i}',candidate_id=key,smiles=smi,
                                  fragment=row['frag'],repeats=n,stereoisomer=i,model_type='covalent_oligomer'))
    atomic_write_json(BASE/'cases.json',dict(cases=cases,scope='Neutral method controls, not new synthesis nominations'))
    with allocation_guard():host_capacity(4,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(4))
    atomic_write_json(BASE/'preparer_owner.json',dict(pid=os.getpid(),host=os.uname().nodename,
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        cpus=sorted(os.sched_getaffinity(0)),started=time.time()))
    ensembles=[]
    with ProcessPoolExecutor(max_workers=4,mp_context=mp.get_context('spawn')) as pool:
        futures={pool.submit(sample_case,case):case for case in cases}
        for future in as_completed(futures):
            result=future.result();ensembles.append(result)
            atomic_write_json(BASE/'preparation_progress.json',dict(completed_cases=len(ensembles),expected_cases=len(cases),case_keys=[e['case']['key'] for e in ensembles]))
            print('Completed conformer case',result['case']['key'],flush=True)
    tasks=[]
    for ensemble in sorted(ensembles,key=lambda e:e['case']['key']):
        case=ensemble['case']
        for target in (('monomer',) if case['repeats']==1 else ('open','stacked')):
            record,fallback=choose(ensemble,target)
            add_dft_task(tasks,case,record,target,fallback)
        if case['repeats']==1:
            record,_=choose(ensemble,'monomer')
            for orientation in ('parallel','antiparallel'):
                description,pair,fallback=segment_pair(case,record,orientation)
                add_dft_task(tasks,description,pair,orientation,fallback)
    atomic_write_json(BASE/'preparation.json',dict(status='finite_DFT_inputs_prepared',created=time.time(),
        controls=list(CONTROLS),cases=len(cases),tasks=tasks,DFT_tasks=len(tasks),
        seed_records=sum(len(e['records']) for e in ensembles),
        conformer_population_qualified=False,all_energies_are_initial_electronic_screening=True,
        solvent='acetonitrile',redox_scope='Neutral resting-state control first; charged-state interactions remain pending',
        symmetry='Enantiomer representatives in achiral continuum; both trimer diastereomers retained',
        limitations=['No DFT Hessians/free energies yet','xTB energies are exploratory, not final stacking energetics',
                     'Covalent folding and separate-segment association are different comparisons',
                     'BSSE, statistical degeneracy and configurational entropy need later qualification']))
    print('Prepared',len(tasks),'finite DFT cases',flush=True)


if __name__=='__main__':main()
