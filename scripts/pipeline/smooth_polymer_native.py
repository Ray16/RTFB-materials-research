#!/usr/bin/env python
"""Model-addressed native thermochemistry for a checked smooth-surface state.

Uses the already controlled electronic-response operators. Records live outside
the primary native cache and cannot be admitted as legacy-model free energies.
Every displaced gradient and SCF is retained for independent reconstruction.
"""
import argparse
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

MODEL=ROOT/'scripts/pipeline/smooth_cds_control.py'
ADAPTER=ROOT/'scripts/pipeline/residual_response_solver.py'
NATIVE_SETTINGS=dict(xc='r2scan',basis='def2-svp',disp='d4',solvent='acetonitrile',
    SCF_conv_tol=1e-12,SCF_conv_tol_grad=1e-8,SCF_max_cycle=200,grid_level=3,
    gradient_grid_response=True,Hessian_grid_response=True,CPHF_grids_matched=True,
    gradient_limit_Ha_Bohr=.00045,Hessian_asymmetry_limit=1e-6,
    FD_steps_Bohr=[.0005,.001],FD_relative_limit=.01,minimum_frequency_limit_cm=-5)


def fd_columns(gradient,symbols):
    gradient=np.asarray(gradient,dtype=float)
    if gradient.shape!=(len(symbols),3) or not np.isfinite(gradient).all():
        raise ValueError('Finite atom-matched Cartesian gradient required')
    hetero=next((i for i,s in enumerate(symbols) if s=='O'),
        next((i for i,s in enumerate(symbols) if s=='Fe'),len(symbols)-1))
    required={(0,1),(hetero,2)}
    required.update(tuple(map(int,np.unravel_index(i,gradient.shape))) for i in np.argsort(abs(gradient).ravel())[-2:])
    return sorted(required)


def validate_partition(partition,energy,method_id):
    if (partition['experimental_method_id']!=method_id or partition['primary_candidate_admission']
            or not partition['derivative_checks_passed'] or partition['symmetry_source']!='xTB autodetection on exact geometry'):
        raise ValueError('Alternate native state has invalid method/scope/derivative/symmetry metadata')
    for cutoff in ('25','50','100'):
        nuclear=partition['nuclear']['nuclear_by_cutoff_eV'][cutoff]
        expected=energy+nuclear+partition['spin_partition_eV']
        if not math.isfinite(expected) or abs(expected-partition['G_by_cutoff_eV'][cutoff])>1e-8:
            raise ValueError('Alternate free energy does not close from its exact electronic/nuclear/spin terms')


def prepare(energy_case,review_path,out):
    original=api.read_plan(energy_case);review=json.loads(review_path.read_text())
    if (review['status']!='alternate_electronic_geometry_and_resource_checked_native_G_pending'
            or review['campaign']!=str(energy_case) or review['native_Hessian_or_G_qualified']):
        raise ValueError('Completed independent alternate electronic review required')
    api.verify_sources(original['producer_sha256']);api.verify_sources(original['input_sha256'])
    api.verify_sources(review['source_sha256'])
    model=api.load_module('native_smooth_declared_model',MODEL)
    controls=ROOT/'results/validation/smooth_model_qualification_2026-10-25'
    witness=api.load_module('native_smooth_model_controls',ROOT/'scripts/validation/review_smooth_model_panel.py').checked_controls(controls)
    identity=hashlib.sha256(json.dumps(model.MODEL,sort_keys=True).encode()).hexdigest()
    if (original['experimental_model']!=model.MODEL or original['experimental_method_id']!=review['portable_method_id']
            or original['experimental_method_id']!=identity):
        raise ValueError('Native and electronic solvent models differ')
    if len(original['tasks'])!=1 or (original['tasks'][0]['charge'],original['tasks'][0]['mult'])!=(0,1):
        raise ValueError('First large native pilot is one neutral singlet')
    out.mkdir(parents=True,exist_ok=False);xyz=out/'input.xyz';shutil.copyfile(energy_case/'input.xyz',xyz)
    inputs={**original['input_sha256'],**review['source_sha256'],**witness}
    for p in (energy_case/'plan.json',review_path,xyz):inputs[str(p)]=api.digest(p)
    producers=dict(original['producer_sha256'])
    for p in (Path(__file__).resolve(),MODEL,ADAPTER,ROOT/'src/redox/qm/hessian_memory.py',
              ROOT/'src/redox/qm/qrrho.py',ROOT/'scripts/validation/audit/native_partition.py'):
        producers[str(p)]=api.digest(p)
    plan=dict(original,input_sha256=inputs,producer_sha256=producers,
        execution_driver=str(Path(__file__).resolve()),raw_root=str(ROOT/'calcs/polymer_models'/out.parent.name/out.name),
        source_electronic_case=str(energy_case),electronic_review=str(review_path),
        maximum_host_RSS_GiB=80,maximum_peak_GPU_MiB=31744,geometry_only_execution=False,
        energy_only_execution=False,prospective_resource_pilot=True,native_pipeline_enabled=False,
        experimental_native_settings=NATIVE_SETTINGS,
        purpose='Separate smooth-model native resource/minimum/derivative/G pilot; primary model cache untouched',created=time.time())
    atomic_write_json(out/'plan.json',plan)
    print(json.dumps(dict(status='experimental_native_pilot_prepared',primary_candidate_admission=False)))


@contextmanager
def execution_budget():
    from gpu4pyscf.df.hessian import rhf,uhf
    from redox.qm.native_thermal import NATIVE_THERMAL
    from redox.qm.hessian_memory import execution_budget,stage_budget_bytes
    policy=execution_budget(NATIVE_THERMAL,24,False)
    originals=[]
    for module in (rhf,uhf):
        previous=module.get_avail_mem
        def bounded(*args,_old=previous,_name=module.__name__,**kwargs):
            return min(_old(*args,**kwargs),stage_budget_bytes(policy,_name,sys._getframe(1).f_code.co_name))
        originals.append((module,previous));module.get_avail_mem=bounded
    try:yield policy
    finally:
        for module,previous in originals:module.get_avail_mem=previous


def run(campaign):
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    verify_gate_owner();plan=api.read_plan(campaign)
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or set(os.environ.get(k) for k in caps)!={'2'}:
        raise ValueError('Two disjoint CPUs and native thread caps=2 required')
    host_capacity(2,minimum_ram_gib=112);gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    guard=json.loads((campaign/'host_guard/status.json').read_text())
    if guard['pid']!=os.getpid() or guard['maximum_RSS_GiB']!=80 or guard['minimum_free_GiB']!=32:
        raise ValueError('The exact same-PID native memory guard is required')
    api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if plan['experimental_native_settings']!=NATIVE_SETTINGS:raise ValueError('Native numerical settings changed')
    task=plan['tasks'][0];xyz=campaign/'input.xyz'
    source_review=json.loads(Path(plan['electronic_review']).read_text())
    if api.digest(xyz)!=plan['scored_geometry_sha256']:raise ValueError('Native geometry changed after energy review')
    raw=Path(plan['raw_root']);raw.mkdir(parents=True,exist_ok=False);shutil.copyfile(xyz,raw/'geometry.xyz')
    stage=campaign/'experimental_native';stage.mkdir(exist_ok=False)
    report=dict(status='experimental_native_resource_and_physical_pilot_running',started=time.time(),
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        raw_root=str(raw),experimental_method_id=plan['experimental_method_id'],primary_candidate_admission=False)
    path=stage/'status.json';atomic_write_json(path,report)
    def stop(signum,frame):raise KeyboardInterrupt('Owned alternate native stop; retain all components and SCF/FD checkpoints')
    previous=signal.signal(signal.SIGTERM,stop)
    try:
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft,qrrho
        from redox.qm.cds_finite import finer_cds_hessian
        from redox.qm.vibrations import project_rigid_modes,write_xtb_hessian
        atoms=read(xyz);symbols=atoms.get_chemical_symbols();coordinates=atoms.positions/.529177210903
        model=api.load_module('native_smooth_actual_model',MODEL)
        if plan['experimental_model']!=model.MODEL:raise ValueError('Declared smooth model changed')
        adapter=api.load_module('native_smooth_checked_response',ADAPTER)
        monitor=api.load_module('native_smooth_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        equations=api.load_module('native_smooth_independent_nuclear',ROOT/'scripts/validation/audit/native_partition.py')
        counter=0
        def solve(coords):
            nonlocal counter
            gc.collect();dft._free_gpu_pool();folder=raw/'SCF_evaluations'/f'{counter:04d}';counter+=1;folder.mkdir(parents=True)
            mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',basis='def2-svp',charge=0,spin=0,verbose=0,max_memory=6000)
            mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
            mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3;mf.chkfile=str(folder/'SCF.h5')
            mf.callback=lambda e:atomic_write_json(folder/'progress.json',dict(cycle=int(e['cycle'])+1,E_Ha=float(e['e_tot']),updated=time.time(),qualified=False))
            atomic_write_json(folder/'settings.json',dict(NATIVE_SETTINGS,charge=0,mult=1,
                coordinates_Bohr=coords.tolist(),experimental_method_id=plan['experimental_method_id']))
            mf.kernel()
            if not mf.converged:raise ValueError('Alternate native SCF failed')
            spin=dft._orbital_diag(mf,0)
            if abs(spin.get('spin_contam',math.inf))>.1:raise ValueError('Alternate native spin gate failed')
            obj=mf.nuc_grad_method();obj.grid_response=True;value=obj.kernel()
            gradient=np.asarray(value.get() if hasattr(value,'get') else value)
            if not np.isfinite(gradient).all():raise ValueError('Nonfinite native gradient')
            atomic_write_json(folder/'evaluation.json',dict(coordinates_Bohr=coords.tolist(),gradient_Ha_Bohr=gradient.tolist(),
                E_Ha=float(mf.e_tot),SCF_converged=True,spin=spin,experimental_method_id=plan['experimental_method_id'],qualified=False))
            return mf,gradient
        with monitor.MemoryMonitor() as memory,(stage/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
            with execution_budget() as policy,adapter.compact_smd_execution(raw/'execution',28),\
                    model.smooth_cds(campaign/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'),finer_cds_hessian(1e-5):
                mf,gradient=solve(coordinates);np.save(raw/'gradient.npy',gradient)
                if abs(gradient).max()>.00045:raise ValueError('Alternate native source is not stationary')
                columns=fd_columns(gradient,symbols)
                mf.cphf_grids=mf.grids.copy();hobj=mf.Hessian();hobj.grid_response=True;hobj.max_cycle=100
                value=hobj.kernel();hessian=np.asarray(value.get() if hasattr(value,'get') else value).transpose(0,2,1,3).reshape(3*len(atoms),3*len(atoms))
                np.save(raw/'hessian.npy',hessian)
                asym=float(abs(hessian-hessian.T).max())
                atomic_write_json(raw/'hessian_diagnostics.json',dict(maximum_asymmetry=asym,all_finite=bool(np.isfinite(hessian).all()),allocation_policy=policy,qualified=False))
                if not np.isfinite(hessian).all() or asym>1e-6:raise ValueError('Alternate native Hessian symmetry/finite gate failed')
                del mf,hobj,value;gc.collect();dft._free_gpu_pool();checks=[]
                for atom,axis in columns:
                    for step in (.0005,.001):
                        gradients=[]
                        for sign in (-1,1):
                            coords=coordinates.copy();coords[atom,axis]+=sign*step
                            mf,g=solve(coords);gradients.append(g);del mf;gc.collect();dft._free_gpu_pool()
                        finite=(gradients[1]-gradients[0]).ravel()/(2*step);analytic=hessian[:,3*atom+axis]
                        error=float(np.linalg.norm(finite-analytic)/max(np.linalg.norm(finite),1e-14))
                        checks.append(dict(atom=atom,axis=axis,step_Bohr=step,relative_error=error))
                        atomic_write_json(raw/'derivative_checks.json',dict(checks=checks,complete=False))
                        if not math.isfinite(error) or error>.01:raise ValueError('Alternate native sampled derivative gate failed')
                atomic_write_json(raw/'derivative_checks.json',dict(checks=checks,complete=True,columns=columns,steps_Bohr=[.0005,.001]))
        if not 0<memory.peak<=31744:raise ValueError('Alternate native resource ceiling failed')
        mass=qrrho.atomic_masses(atoms.numbers);projected,rank=project_rigid_modes(hessian,atoms.positions,mass)
        write_xtb_hessian(raw/'hessian',projected);xtb=qrrho.LIBRARY.parent.parent/'bin/xtb'
        with (raw/'symmetry_thermo.log').open('x') as stream:
            process=subprocess.run([str(xtb),'thermo','--sthr','50','--temp','298.15','geometry.xyz','hessian'],cwd=raw,stdout=stream,stderr=subprocess.STDOUT,timeout=120)
        if process.returncode:raise ValueError('Symmetry autodetection failed')
        symmetry=int(re.search(r'rotational number\s+(\d+)',(raw/'symmetry_thermo.log').read_text())[1])
        nuclear=qrrho.nuclear_thermals(hessian,atoms.positions,atoms.numbers,symmetry)
        if min(nuclear['frequencies_cm'])<-5:raise ValueError('Alternate native geometry is not a minimum')
        differences={str(k):abs(equations.nuclear_equations(np.array(nuclear['frequencies_cm']),atoms.positions,
            np.array(nuclear['atomic_masses_amu']),symmetry,k)-nuclear['nuclear_by_cutoff_eV'][str(k)]) for k in (25,50,100)}
        if max(differences.values())>1e-8:raise ValueError('Independent alternate nuclear equations failed')
        partition=dict(experimental_method_id=plan['experimental_method_id'],primary_candidate_admission=False,
            native_settings=NATIVE_SETTINGS,
            derivative_checks_passed=True,symmetry_source='xTB autodetection on exact geometry',symmetry_number=symmetry,
            nuclear=nuclear,spin_partition_eV=0.,G_by_cutoff_eV={k:source_review['E_smd_eV']+v for k,v in nuclear['nuclear_by_cutoff_eV'].items()},
            derivative_checks=checks,required_columns=columns,independent_nuclear_errors_eV=differences)
        validate_partition(partition,source_review['E_smd_eV'],plan['experimental_method_id'])
        atomic_write_json(raw/'experimental_partition.json',partition)
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        outputs={str(p):api.digest(p) for p in raw.rglob('*') if p.is_file() and p.suffix in ('.json','.xyz','.npy')}
        outputs[str(campaign/'model/model.json')]=api.digest(campaign/'model/model.json')
        report.update(status='experimental_native_finished_independent_review_pending',observed_peak_GPU_MiB=memory.peak,
            G_eV=partition['G_by_cutoff_eV']['50'],minimum_internal_frequency_cm=min(nuclear['frequencies_cm']),
            gradient_max_Ha_Bohr=float(abs(gradient).max()),maximum_sampled_column_error=max(c['relative_error'] for c in checks),
            output_sha256=outputs)
    except BaseException as error:
        report.update(status='experimental_native_failed_preserved',error=repr(error));raise
    finally:
        signal.signal(signal.SIGTERM,previous);report['finished']=time.time();atomic_write_json(path,report)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run'))
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--energy-case',type=Path);p.add_argument('--review',type=Path)
    a=p.parse_args()
    if a.action=='prepare':
        if not a.energy_case or not a.review:p.error('prepare requires energy-case and review')
        prepare(a.energy_case.resolve(),a.review.resolve(),a.campaign.resolve())
    else:run(a.campaign.resolve())


if __name__=='__main__':main()
