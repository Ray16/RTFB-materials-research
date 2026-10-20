#!/usr/bin/env python
"""Resolve steric pair-pose failures without rerunning completed quantum states."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import time

import numpy as np
from ase.io import read,write
from ase.data import vdw_radii

from prepare_pi_stacking import BASE,choose,xtb
from redox.build.polymer_contacts import plane,contact_metrics
from redox.core.io import atomic_write_json
from redox.workflows.checkpoints import verify_state_checkpoint
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pair_task(task):
    allowed=sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0,{allowed[mp.current_process()._identity[0]%len(allowed)]})
    mono_key=task['candidate_id']+'__n1__stereo0__monomer'
    report=verify_state_checkpoint(BASE/'done'/f'{mono_key}.json')
    source=Path(report['state_dir'])/'opt.xyz';atoms=read(source)
    groups=report['task']['pi_groups'];group=groups[0]
    center,normal,axis,_=plane(atoms.positions[group]);basis=np.array([axis,np.cross(normal,axis),normal])
    first=(atoms.positions-center)@basis.T;other=first.copy()
    if task['target']=='antiparallel':other[:,:2]*=-1
    radii=vdw_radii[atoms.numbers]
    poses=[]
    for z in (3.6,3.8,4.,4.2,4.6):
        for slip in (1.5,2.,2.5,3.,3.5,4.,4.5,5.):
            for angle in (0.,90.,180.,270.):
                theta=np.deg2rad(angle);shift=np.array([slip*np.cos(theta),slip*np.sin(theta),z])
                moved=other+shift
                distances=np.linalg.norm(first[:,None,:]-moved[None,:,:],axis=2)
                ratio=float((distances/(radii[:,None]+radii[None,:])).min())
                if distances.min()<1.15 or ratio<.60:continue
                score=(z-3.6)**2+.15*(slip-2.)**2
                poses.append((score,z,slip,angle,ratio,moved))
    poses.sort(key=lambda x:x[0])
    if not poses:raise RuntimeError('No clash-free offset pair pose found')
    choices=[]
    for pose in poses:
        if choices and pose[3]==choices[0][3]:continue
        choices.append(pose)
        if len(choices)==2:break
    directory=BASE/'pair_pose_recovery'/task['key'];directory.mkdir(parents=True,exist_ok=False)
    n=len(atoms);pair_groups=[group,[i+n for i in group]]
    records=[]
    for i,(_,z,slip,angle,ratio,moved) in enumerate(choices):
        pair=atoms.copy();pair.positions=first@basis+center
        second=atoms.copy();second.positions=moved@basis+center;pair+=second
        d=directory/f'pose_{i}';d.mkdir();write(d/'input.xyz',pair,format='xyz')
        record=dict(starting_separation_A=z,starting_slip_A=slip,offset_direction_deg=angle,
                    minimum_vdw_distance_ratio=ratio,pi_groups=pair_groups,
                    monomer_DFT_xyz=str(source),monomer_DFT_xyz_sha256=sha(source),
                    initial_xyz=str(d/'input.xyz'),initial_metrics=contact_metrics(pair.positions,pair_groups))
        try:record.update(xtb(d,task['smiles'],pair_groups))
        except Exception as exc:record.update(error=repr(exc),xyz=str(d/'input.xyz'))
        atomic_write_json(d/'record.json',record);records.append(record)
    successful=[r for r in records if 'E_eV' in r]
    selected=min(successful,key=lambda r:r['E_eV']) if successful else records[0]
    atomic_write_json(directory/'selection.json',dict(selected=selected,candidates=records,
        method='Clash-filtered offset grid; two directions relaxed with unchanged GFN2-xTB/ALPB(MeCN)',
        actual_DFT_attempts=0))
    return task['key'],selected


def main():
    assert not (BASE/'preparation.json').exists(),'Preserve an existing complete preparation'
    owner=json.loads((BASE/'preparer_owner.json').read_text());proc=Path('/proc')/str(owner['pid'])
    assert not proc.exists() or proc.joinpath('stat').read_text().rsplit(')',1)[1].split()[19]!=owner['process_start'],'Original conformer producer still live'
    archive=BASE/'preserved_attempts/overlapping_segment_start';archive.mkdir(parents=True,exist_ok=True)
    shutil.copyfile('/tmp/redox_pi_stacking_preparation_v2.log',archive/'preparation.log')
    atomic_write_json(archive/'reason.json',dict(reason='Single rigid offset pose had intersegment overlap; preserve all 12 completed covalent ensembles and three completed DFT monomers. Use clash-filtered relative offsets for unattempted segment cases only.',original_owner=owner,original_owner_proven_terminal=True))
    plan=json.loads((BASE/'plan.json').read_text());pair_tasks=[t for t in plan['tasks'] if t['model_type']=='two_nonbonded_segments']
    with allocation_guard():host_capacity(4,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(4))
    atomic_write_json(BASE/'pair_preparer_owner.json',dict(pid=os.getpid(),host=os.uname().nodename,
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],cpus=sorted(os.sched_getaffinity(0)),started=time.time()))
    pairs={}
    with ProcessPoolExecutor(max_workers=4,mp_context=mp.get_context('spawn')) as pool:
        for future in as_completed([pool.submit(pair_task,task) for task in pair_tasks]):
            key,record=future.result();pairs[key]=record;print('Prepared clash-free segment pair',key,flush=True)
    tasks=[]
    for task in plan['tasks']:
        if task['model_type']=='two_nonbonded_segments':
            record=pairs[task['key']];fallback='E_eV' not in record
        else:
            ensemble=json.loads((BASE/'seeds'/task['model_key']/'ensemble.json').read_text())
            record,fallback=choose(ensemble,task['target'])
        directory=BASE/'resolved_inputs'/task['key'];directory.mkdir(parents=True,exist_ok=False)
        destination=directory/'input.xyz';shutil.copyfile(record['xyz'],destination)
        tasks.append(dict(task,input_xyz=str(destination),pi_groups=record['pi_groups'],
            seed_stage='prepared_clash_free_start; DFT_required' if fallback else 'converged_topology_checked_xTB_ALPB'))
    atomic_write_json(BASE/'preparation.json',dict(status='all_finite_DFT_inputs_prepared',created=time.time(),tasks=tasks,
        DFT_tasks=len(tasks),controls=sorted({t['candidate_id'] for t in tasks}),covalent_ensembles=12,
        source='Preserved covalent ensembles plus recovered clash-free segment pairs',
        source_sha256=sha(Path(__file__)),quantum_state_retries=0,
        solvent='acetonitrile',scope='Neutral electronic method controls; no interaction free energies'))
    print('All',len(tasks),'DFT input poses prepared without repeating completed quantum states',flush=True)


if __name__=='__main__':main()
