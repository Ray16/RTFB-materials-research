#!/usr/bin/env python
"""Freeze an exact-state panel and measure real UMA inference, including Nitro."""
import argparse
import gc
from importlib import metadata
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.qm.uma_nitro import state_key


def prepare(out):
    import numpy as np
    from ase.io import read, write
    from redox.qm import uma
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():
        host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    identity=uma.model_identity('uma-s-1p2p1')
    out.mkdir(parents=True,exist_ok=False)
    cases=[];witnesses={}
    selected=[(gid,st) for gid in ('design_pmdi_ethyl','design_ndi_ethyl') for st in ('neu','red1','red2')]
    selected += [('viologen','ox2'),('viologen','ox1')]
    for gid,st in selected:
        source=ROOT/'calcs/uma'/gid/st
        record=json.loads((source/'result.json').read_text())
        xyz=source/'relaxed.xyz'
        if not record['converged']:
            raise ValueError('A converged stored warm-start geometry is required')
        atoms=read(xyz);q=int(record['charge']);mult=int(record['chosen_mult'])
        state_key([int(z) for z in atoms.numbers],q,mult)
        fingerprinted=record.get('model_identity')==identity and record.get('geometry_sha256')==api.digest(xyz)
        cases.append(dict(id=gid+'__'+st,charge=q,mult=mult,source=str(xyz),
            stored_energy_eV=record['energy_eV'] if fingerprinted else None,
            stored_energy_scope='fingerprinted same-checkpoint reference' if fingerprinted else 'unfingerprinted historical value excluded from admission'))
        witnesses.update({str(p):api.digest(p) for p in (xyz,source/'result.json')})
    parent=ROOT/'results/polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume'
    checkpoint=api.load_module('Nitro_exact_polymer_checkpoint',ROOT/'scripts/validation/qualify_xc_checkpoint.py')
    frame,row,files=checkpoint.selected_checkpoint(parent)
    polymer=out/'polymer_checkpoint.xyz';polymer.write_text(frame);witnesses.update(files)
    from redox.validation.protonation import hydrogen_topology_status
    from redox.build.geometry import geometry_stereo_status
    task=api.read_plan(parent)['tasks'][0]
    if (hydrogen_topology_status(task['smiles'],polymer)!='ok'
            or geometry_stereo_status(task['smiles'],polymer) not in ('ok','not_applicable')):
        raise ValueError('Polymer checkpoint graph or stereochemistry differs')
    witnesses[str(polymer)]=api.digest(polymer)
    cases.append(dict(id='PMDI_five_repeat_two_grafts',charge=0,mult=1,source=str(polymer),
        stored_energy_eV=None,stored_energy_scope='Unqualified DFT warm checkpoint; UMA execution control only'))
    triplet=dict(cases[0],id='PMDI_fixed_geometry_triplet',mult=3,stored_energy_eV=None,
        stored_energy_scope='Spin-conditioning method control at the singlet geometry')
    cases.append(triplet)
    for case in cases:
        folder=out/'inputs'/case['id'];folder.mkdir(parents=True)
        atoms=read(case['source']);base=atoms.positions.copy();rng=np.random.default_rng(20261023)
        angle=.37;rotation=np.array([[math.cos(angle),-math.sin(angle),0],[math.sin(angle),math.cos(angle),0],[0,0,1]])
        frames=[base,base+rng.normal(0,.002,base.shape),base+np.array([.12,-.07,.03]),base@rotation.T]
        case['frames']=[];case['atoms']=len(atoms)
        for i,coordinates in enumerate(frames):
            atoms.positions=coordinates;atoms.info.update(charge=case['charge'],spin=case['mult'])
            target=folder/f'frame_{i}.xyz';write(target,atoms,format='extxyz')
            case['frames'].append(str(target));witnesses[str(target)]=api.digest(target)
    plan=dict(model='uma-s-1p2p1',task='omol',model_identity=identity,cases=cases,
        input_sha256=witnesses,random_seed=20261023,model_seed=41,warmup_calls=4,timed_calls=20,
        energy_tolerance_eV_per_atom=1e-4,minimum_energy_tolerance_eV=.002,
        force_RMSE_tolerance_eV_A=.002,force_max_tolerance_eV_A=.01,
        repeat_energy_tolerance_eV_per_atom=1e-5,repeat_force_tolerance_eV_A=1e-4,
        timing='CUDA synchronize around every actual energy/force call; ASE cache reset; startup measured separately',
        roles='Existing exact molecules and explicit polymer/spin method controls; no candidate or DFT thermochemistry admission',
        producer_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(),ROOT/'src/redox/qm/uma_nitro.py',Path(uma.__file__))})
    atomic_write_json(out/'plan.json',plan)
    print(json.dumps(dict(prepared_cases=len(cases),frames=4*len(cases),historical_energy_references=sum(c['stored_energy_eV'] is not None for c in cases))))


def run(campaign,mode,out,case_id=None):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):
        raise ValueError('Two disjoint CPUs and native caps=2 required')
    import numpy as np
    import torch
    from ase.io import read
    from redox.qm import uma
    from redox.qm.uma_nitro import make_calculator,check_actual_nitro
    torch.set_num_threads(2)
    plan=api.read_plan(campaign);api.verify_sources(plan['input_sha256']);api.verify_sources(plan['producer_sha256'])
    identity=uma.model_identity(plan['model'])
    if (identity['checkpoint_sha256']!=plan['model_identity']['checkpoint_sha256']
            or identity['reference_sha256']!=plan['model_identity']['reference_sha256']):
        raise ValueError('Weights or atomic references differ from the frozen panel')
    cases=[c for c in plan['cases'] if case_id is None or c['id']==case_id]
    if not cases:raise ValueError('Unknown frozen case')
    out.mkdir(parents=True,exist_ok=False)
    report=dict(status='UMA_point_evaluation_running',mode=mode,started=time.time(),
        host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),GPU=os.environ['CUDA_VISIBLE_DEVICES'],
        device_name=torch.cuda.get_device_name(),capability=list(torch.cuda.get_device_capability()),
        fairchem_version=metadata.version('fairchem-core'),torch_version=torch.__version__,
        model_identity=identity,plan_sha256=api.digest(campaign/'plan.json'),cases=[])
    atomic_write_json(out/'status.json',report)
    try:
        if mode=='nitro' and report['capability'][0] not in (8,9,10,11,12):
            raise RuntimeError('This GPU cannot run actual Nitro')
        for case in cases:
            begun=time.perf_counter()
            if mode=='legacy':
                if report['fairchem_version']!='2.21.0':raise ValueError('Legacy reference requires its original loader')
                calc=uma.make_calculator(plan['model'],'cuda')
            elif mode=='reference':
                from fairchem.core import FAIRChemCalculator,pretrained_mlip
                predictor=pretrained_mlip.get_predict_unit(plan['model'],device='cuda',inference_settings='batch',seed=plan['model_seed'])
                calc=FAIRChemCalculator(predictor,task_name=plan['task'])
            else:
                a=read(case['frames'][0])
                calc=make_calculator([int(z) for z in a.numbers],case['charge'],case['mult'],seed=plan['model_seed'])
            torch.cuda.synchronize();setup=time.perf_counter()-begun
            values=[]
            def evaluate(index):
                a=read(case['frames'][index]);a.info.update(charge=case['charge'],spin=case['mult']);a.calc=calc
                calc.reset();torch.cuda.synchronize();start=time.perf_counter()
                forces=np.asarray(a.get_forces(),dtype=float);energy=float(a.get_potential_energy())
                torch.cuda.synchronize();elapsed=time.perf_counter()-start
                if not np.isfinite(energy) or forces.shape!=(len(a),3) or not np.isfinite(forces).all():
                    raise ValueError('Invalid UMA energy/force output')
                if mode=='nitro':check_actual_nitro(calc.predictor)
                return dict(energy_eV=energy,forces_eV_A=forces.tolist(),seconds=elapsed)
            for index in range(4):values.append(evaluate(index))
            for i in range(plan['warmup_calls']):evaluate(i%4)
            seconds=[evaluate(i%4)['seconds'] for i in range(plan['timed_calls'])]
            repeat=evaluate(0)
            error_E=abs(repeat['energy_eV']-values[0]['energy_eV'])
            error_F=float(np.max(abs(np.array(repeat['forces_eV_A'])-values[0]['forces_eV_A'])))
            if (error_E>plan['repeat_energy_tolerance_eV_per_atom']*case['atoms'] or error_F>plan['repeat_force_tolerance_eV_A']):
                raise ValueError('Repeated fixed-geometry inference is not reproducible')
            row=dict(id=case['id'],charge=case['charge'],mult=case['mult'],atoms=case['atoms'],frames=values,
                setup_seconds=setup,median_step_seconds=float(np.median(seconds)),step_seconds=seconds,
                repeat_energy_error_eV=error_E,repeat_force_error_eV_A=error_F,
                actual_mode='nitro' if mode=='nitro' else mode)
            if mode=='nitro':row['execution']=calc.execution
            report['cases'].append(row);atomic_write_json(out/'status.json',report)
            del calc
            if mode=='reference':del predictor
            gc.collect();torch.cuda.empty_cache()
        api.verify_sources(plan['input_sha256']);api.verify_sources(plan['producer_sha256'])
        report['status']='UMA_point_evaluations_complete_cross_mode_review_pending'
    except BaseException as error:
        report.update(status='UMA_point_evaluation_failed_preserved',error=repr(error));raise
    finally:
        report['finished']=time.time();atomic_write_json(out/'status.json',report)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run'))
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path)
    p.add_argument('--mode',choices=('legacy','reference','nitro'));p.add_argument('--case')
    a=p.parse_args()
    if a.action=='prepare':prepare(a.campaign.resolve())
    else:
        if a.mode is None or a.out is None:p.error('run requires --mode and --out')
        run(a.campaign.resolve(),a.mode,a.out.resolve(),a.case)


if __name__=='__main__':main()
