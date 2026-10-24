#!/usr/bin/env python
"""Qualified GPU geometry stage, independently checkpointed from active energies."""
import argparse
from contextlib import contextmanager,redirect_stdout,redirect_stderr
import json
import math
import os
from pathlib import Path
import shutil
import signal
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,verify_gate_owner,gpu_suitability


def prepare(source,pilot,out):
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    old=api.read_plan(source);point=api.read_plan(pilot)
    for p in (old,point):api.verify_sources(p['producer_sha256']);api.verify_sources(p['input_sha256'])
    status=pilot/'geometry/status.json';qualified=json.loads(status.read_text())
    if qualified['status']!='split_geometry_resource_and_reference_point_passed':
        raise ValueError('The exact geometry resource/reference pilot must pass')
    api.verify_sources(qualified['output_sha256'])
    if old['tasks']!=[point['task']] or old['quantum_settings']!=point['quantum_settings']:
        raise ValueError('Pilot task or geometry protocol differs')
    owner=api.load_module('geometry_stage_pilot_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if owner.owner_alive(qualified):raise RuntimeError('Resource pilot owner must finish before stage handoff')
    out.mkdir(parents=True,exist_ok=False);shutil.copytree(source/'controls',out/'controls')
    for stage in ('energy','thermal'):
        for name in ('tasks','done','claims','failures','workers','pilot','launches'):(out/stage/name).mkdir(parents=True)
    key=point['task']['key'];target=out/'live_inputs'/key;target.mkdir(parents=True)
    xyz=target/'input.xyz';xyz.write_bytes(Path(point['geometry_input']).read_bytes())
    ready=json.loads((source/'live_inputs'/key/'task.json').read_text());ready['input_xyz']=str(xyz)
    ready['input_sha256'][str(xyz)]=api.digest(xyz);ready['source_restart_checkpoint']=point['geometry_reference']
    atomic_write_json(target/'task.json',ready)
    inputs={**old['input_sha256'],**point['input_sha256'],str(status):api.digest(status),str(pilot/'plan.json'):api.digest(pilot/'plan.json'),str(target/'task.json'):api.digest(target/'task.json'),str(xyz):api.digest(xyz)}
    sources=dict(old['producer_sha256'],**point['producer_sha256'])
    for p in (Path(__file__).resolve(),ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py',ROOT/'scripts/pipeline/resume_polymer_forces_adaptive.py'):
        sources[str(p)]=api.digest(p)
    plan=dict(old,raw_root=str(ROOT/'calcs/polymer_models'/out.parent.name/out.name),
        producer_sha256=sources,input_sha256=inputs,geometry_resource_pilot=str(pilot),
        geometry_resource_policy=point['geometry'],energy_only_execution=True,native_pipeline_enabled=False,
        geometry_only_execution=True,externally_managed_stages=True,
        purpose='Geometry optimization and total-energy force validation; active electronic parent and native G are separate stages',prepared=time.time())
    atomic_write_json(out/'plan.json',plan)
    print(json.dumps(dict(status='separate_geometry_stage_prepared',checkpoint_evaluation=point['geometry_reference']['evaluation'])))


@contextmanager
def geometry_engine(directory):
    from redox.qm import dft
    from redox.qm.stable_geometry_scf import stabilize
    adapter=api.load_module('geometry_stage_memory_execution',ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py')
    forces=api.load_module('geometry_stage_conservative_forces',ROOT/'scripts/pipeline/resume_polymer_forces_adaptive.py')
    counter=0
    with adapter.compact_smd_execution(directory/'execution',budget_GiB=28),forces.adaptive_forces(directory/'CDS'):
        original=dft._build_smd_mf
        def build(*args,**kwargs):
            nonlocal counter
            mf=original(*args,**kwargs);folder=directory/'SCF_checkpoints'/f'{counter:04d}';counter+=1;folder.mkdir(parents=True)
            stabilize(mf,folder/'recovery.jsonl');mf.chkfile=str(folder/'SCF.h5')
            previous=getattr(mf,'callback',None)
            def callback(environment):
                if callable(previous):previous(environment)
                atomic_write_json(folder/'progress.json',dict(cycle=int(environment['cycle'])+1,
                    E_Ha=float(environment['e_tot']),updated=time.time(),qualified_result=False,
                    conv_tol=float(mf.conv_tol),xc=mf.xc,basis=repr(mf.mol.basis)))
            mf.callback=callback
            return mf
        dft._build_smd_mf=build
        try:yield forces
        finally:dft._build_smd_mf=original


def run(campaign):
    verify_gate_owner();caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):raise ValueError('Two allocated CPUs and matching native caps required')
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    policy=plan['geometry_resource_policy'];host_capacity(2,minimum_ram_gib=max(192,policy['minimum_host_RAM_GiB']))
    if gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))['memory_MiB']<policy['minimum_GPU_MiB']:raise ValueError('Wrong device memory class')
    task=plan['tasks'][0];key=task['key'];claim=campaign/'energy/claims'/key;claim.mkdir()
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    atomic_write_json(claim/'owner.json',owner)
    directory=campaign/'geometry';directory.mkdir(exist_ok=False)
    raw=Path(plan['raw_root']);sd=raw/'energy'/key/'rest';sd.mkdir(parents=True,exist_ok=False)
    report=dict(status='GPU_geometry_stage_running',task=task,owner=owner,state_dir=str(sd),started=time.time(),
        active_electronic_energy_computed=False,native_Hessian_qualified=False,G_computed=False)
    atomic_write_json(directory/'status.json',report)
    def stop(signum,frame):raise KeyboardInterrupt('Owned geometry stop; checkpoints preserved')
    before=signal.signal(signal.SIGTERM,stop)
    try:
        from redox.qm import dft
        from redox.workflows.geometry_checkpoints import checkpoint_callback
        monitor=api.load_module('geometry_stage_memory_monitor',ROOT/'scripts/pipeline/run_pi_stacking.py')
        audit=api.load_module('geometry_stage_structural_audit',ROOT/'scripts/validation/audit_polymer_geometries.py')
        xyz=campaign/'live_inputs'/key/'input.xyz'
        with monitor.MemoryMonitor() as memory,geometry_engine(raw), (directory/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
            with api.optimizer_settings(plan['quantum_settings'],sd/'geometry_checkpoints'):
                result=dft.dft_smd(xyz,task['charge'],task['mult'],backend='gpu',do_opt=True,
                    do_smd=False,do_gas=False,do_freq=False,opt_out=sd/'opt.xyz',
                    conv_params=dict(convergence_gmax=1e-4,convergence_grms=6e-5,convergence_energy=1e-7,
                        convergence_dmax=.002,convergence_drms=.001,callback=checkpoint_callback(sd/'geometry_checkpoints')))
            if not result['geometry_optimization_converged'] or not result['optimization_SCF_converged']:
                raise ValueError('Geometry optimization did not converge')
            result.update(id=key,state='rest',n_e=0,geometry_only_result=True)
            atomic_write_json(sd/'result.json',result)
            review=audit.audit(task,sd/'opt.xyz')
            if review['issues']:raise ValueError('Optimized graph/geometry audit failed: '+repr(review['issues']))
            atomic_write_json(directory/'geometry_review.json',review)
            forces=api.load_module('geometry_final_force_checks',ROOT/'scripts/pipeline/resume_polymer_forces.py')
            forcecheck=forces.final_energy_check(sd/'opt.xyz',sd/'total_energy_force_check')
            baseline=json.loads((sd/'total_energy_force_check/baseline.json').read_text())
            maximum=max(abs(x) for row in baseline['gradient_Ha_Bohr'] for x in row)
            if not math.isfinite(maximum) or maximum>.00045:raise ValueError('Independent final force is not stationary')
        if not 0<memory.peak<=policy['maximum_peak_GPU_MiB']:raise ValueError('Geometry device-memory bound exceeded')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        report.update(status='stationary_geometry_checked_active_parent_and_native_minimum_pending',
            observed_peak_GPU_MiB=memory.peak,gradient_max_Ha_Bohr=maximum,conservative_total_force_check=forcecheck,
            output_sha256={str(p):api.digest(p) for p in [sd/'opt.xyz',sd/'result.json',directory/'geometry_review.json',
                *list((sd/'total_energy_force_check').glob('*.json'))]})
    except BaseException as error:
        report.update(status='geometry_stage_failed_preserved',error=repr(error));raise
    finally:
        signal.signal(signal.SIGTERM,before);report['finished']=time.time();atomic_write_json(directory/'status.json',report)
        shutil.rmtree(claim)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run'))
    for n in ('source','pilot','campaign'):p.add_argument('--'+n,type=Path)
    a=p.parse_args()
    if a.campaign is None:p.error('--campaign required')
    if a.action=='prepare':
        if not a.source or not a.pilot:p.error('prepare requires source and pilot')
        prepare(a.source.resolve(),a.pilot.resolve(),a.campaign.resolve())
    else:run(a.campaign.resolve())


if __name__=='__main__':main()
