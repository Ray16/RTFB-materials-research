#!/usr/bin/env python
"""Resource/continuity pilots for GPU geometry and CPU active-energy stages.

Both use the existing quantum levels. No starting-point evaluation is admitted as
an optimized state, free energy or synthesis decision.
"""
import argparse
from contextlib import contextmanager,redirect_stdout,redirect_stderr
import json
import math
import os
from pathlib import Path
import shutil
import signal
import socket
import threading
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


CAPS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')


def prepare(source,reference,CPU_controls,out):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    parent=api.read_plan(source);old=api.read_plan(reference)
    for plan in (parent,old):
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if (parent['tasks']!=old['tasks'] or len(parent['tasks'])!=1
            or (parent['tasks'][0]['charge'],parent['tasks'][0]['mult'])!=(0,1)
            or any(parent[k]!=old[k] for k in ('quantum_settings','thermochemistry','solvent'))):
        raise ValueError('One exact neutral-singlet parent and matched reference required')
    owners=api.load_module('split_pilot_stopped_owners',ROOT/'scripts/pipeline/prepare_polymer_postpilot_continuation.py').verify_stopped(source)
    CPU=json.loads(CPU_controls.read_text())
    if CPU['status']!='two_spin_CPU_GPU_energy_and_gradient_equivalence_passed' or len(CPU['cases'])!=2:
        raise ValueError('Both CPU/GPU execution controls must pass')
    api.verify_sources(CPU['source_sha256'])
    controlroot=ROOT/'results/validation/pcm_switch_curvature_guarded_2026-10-23'
    checks=api.load_module('split_pilot_GPU_controls',ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py').checked_controls(controlroot,28)
    frame,row,witnesses=api.load_module('split_pilot_checkpoint',ROOT/'scripts/validation/qualify_xc_checkpoint.py').selected_checkpoint(source)
    refstatus=reference/'energy/pilot/status.json';ref=json.loads(refstatus.read_text())
    if ref['status']!='memory_pilot_passed':raise ValueError('Checked large-device reference required')
    api.verify_sources(ref['output_sha256'])
    out.mkdir(parents=True,exist_ok=False)
    geometry=out/'geometry_input.xyz';geometry.write_text(frame)
    key=parent['tasks'][0]['key'];energy=out/'electronic_input.xyz'
    original=reference/'live_inputs'/key/'input.xyz';shutil.copyfile(original,energy)
    inputs=dict(witnesses)
    for p in (source/'plan.json',reference/'plan.json',CPU_controls,refstatus,refstatus.parent/'energy_and_gradient.json',geometry,energy,original):inputs[str(p)]=api.digest(p)
    sources={**checks,**CPU['source_sha256']}
    for p in (Path(__file__).resolve(),ROOT/'scripts/pipeline/cpu_quantum_execution.py',ROOT/'scripts/pipeline/resume_polymer_forces_adaptive.py'):
        sources[str(p)]=api.digest(p)
    plan=dict(task=parent['tasks'][0],quantum_settings=parent['quantum_settings'],solvent=parent['solvent'],
        geometry_input=str(geometry),electronic_input=str(energy),geometry_reference=row,
        electronic_reference=str(refstatus.parent/'energy_and_gradient.json'),input_sha256=inputs,producer_sha256=sources,
        source_owners=owners,prepared=time.time(),
        geometry=dict(workspace_GiB=28,minimum_GPU_MiB=32768,maximum_peak_GPU_MiB=31744,
            minimum_host_RAM_GiB=144,energy_tolerance_Ha=1e-7,gradient_max_tolerance_Ha_Bohr=1e-6),
        electronic=dict(library_memory_MiB=192*1024,maximum_RSS_GiB=300,minimum_host_RAM_GiB=400,
            maximum_scratch_GiB=1024,minimum_free_scratch_GiB=64,energy_tolerance_eV=1e-5),
        optimized_geometry_qualified=False,G_qualified=False)
    atomic_write_json(out/'plan.json',plan)
    print(json.dumps(dict(status='split_resource_pilots_prepared',atoms=plan['task']['atoms'],checkpoint_evaluation=row['evaluation'])))


@contextmanager
def CPU_resources(scratch,settings,out):
    import psutil
    done=threading.Event();samples=[]
    def observe():
        while not done.is_set():
            sizes=0
            try:
                for p in scratch.rglob('*'):
                    try:
                        if p.is_file():sizes+=p.stat().st_size
                    except FileNotFoundError:pass  # A completed PySCF temporary can disappear.
            except Exception as error:
                atomic_write_json(out/'resource_monitor_error.json',dict(error=repr(error)))
                os.kill(os.getpid(),signal.SIGTERM);return
            row=dict(time=time.time(),RSS_bytes=psutil.Process().memory_info().rss,scratch_bytes=sizes,
                free_scratch_bytes=shutil.disk_usage(scratch).free,available_RAM_bytes=psutil.virtual_memory().available)
            samples.append(row);atomic_write_json(out/'resource_progress.json',row)
            if (row['RSS_bytes']>settings['maximum_RSS_GiB']*2**30
                    or sizes>settings['maximum_scratch_GiB']*2**30
                    or row['free_scratch_bytes']<settings['minimum_free_scratch_GiB']*2**30
                    or row['available_RAM_bytes']<32*2**30):
                atomic_write_json(out/'resource_stop.json',row);os.kill(os.getpid(),signal.SIGTERM);return
            done.wait(5)
    worker=threading.Thread(target=observe,daemon=True);worker.start()
    try:yield
    finally:
        done.set();worker.join(timeout=6)
        atomic_write_json(out/'resources.json',dict(samples=samples,checkpoints_are_unqualified=True))


def run(campaign,phase,threads,scratch=None):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,verify_gate_owner,gpu_suitability
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if any(os.environ.get(k)!=str(threads) for k in CAPS):raise ValueError('All native caps must match the allocated CPUs')
    if phase=='geometry':
        verify_gate_owner()
        if threads!=2 or len(os.sched_getaffinity(0))!=2:raise ValueError('GPU geometry pilot requires two disjoint CPUs')
        host_capacity(threads,minimum_ram_gib=plan[phase]['minimum_host_RAM_GiB'])
        if gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))['memory_MiB']<plan[phase]['minimum_GPU_MiB']:
            raise ValueError('Insufficient device memory class')
    else:
        with allocation_guard():
            capacity=host_capacity(threads,minimum_ram_gib=plan[phase]['minimum_host_RAM_GiB'])
            if plan[phase]['maximum_RSS_GiB']*2**30>.75*capacity['available_RAM_bytes']:
                raise RuntimeError('Available RAM cannot cover the declared CPU bound')
            os.sched_setaffinity(0,cpu_allocation(threads))
        if scratch is None or not scratch.is_dir() or os.environ.get('PYSCF_TMPDIR')!=str(scratch):
            raise ValueError('Use a fresh task-owned local scratch directory via PYSCF_TMPDIR')
        if shutil.disk_usage(scratch).free<(plan[phase]['maximum_scratch_GiB']+64)*2**30:
            raise RuntimeError('Insufficient local scratch capacity')
    out=campaign/phase;out.mkdir(exist_ok=False)
    report=dict(status='split_resource_point_running',phase=phase,started=time.time(),host=socket.gethostname(),
        pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),G_qualified=False,optimized_geometry_qualified=False)
    atomic_write_json(out/'status.json',report)
    def stop(signum,frame):raise KeyboardInterrupt('Owned resource stop; checkpoints preserved')
    before=signal.signal(signal.SIGTERM,stop)
    try:
        import numpy as np
        import pyscf.lib
        from redox.qm import dft
        pyscf.lib.num_threads(threads)
        if phase=='geometry':
            from ase.io import read
            from pyscf import gto
            engine=api.load_module('split_geometry_memory_observer',ROOT/'scripts/pipeline/run_pi_stacking.py')
            adapter=api.load_module('split_geometry_checked_adapter',ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py')
            forces=api.load_module('split_geometry_conservative_forces',ROOT/'scripts/pipeline/resume_polymer_forces_adaptive.py')
            atoms=read(plan['geometry_input'])
            with engine.MemoryMonitor() as memory, adapter.compact_smd_execution(out/'execution',budget_GiB=28), forces.adaptive_forces(out/'CDS'), (out/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
                mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),basis='def2-svp',charge=0,spin=0,verbose=0,max_memory=6000)
                mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3
                mf.conv_tol,mf.conv_tol_grad,mf.max_cycle=1e-12,1e-8,300;mf.chkfile=str(out/'SCF.h5')
                mf.kernel()
                if not mf.converged:raise ValueError('Geometry-level SCF did not converge')
                gradient=mf.nuc_grad_method();gradient.grid_response=True
                raw=gradient.kernel();values=np.asarray(raw.get() if hasattr(raw,'get') else raw)
                if not np.isfinite(values).all():raise ValueError('Nonfinite geometry gradient')
                result=dict(energy_Ha=float(mf.e_tot),gradient_max_Ha_Bohr=float(abs(values).max()))
                np.save(out/'gradient_Ha_Bohr.npy',values);atomic_write_json(out/'point.json',result)
            report['observed_peak_GPU_MiB']=memory.peak
            errors=dict(energy_Ha=abs(result['energy_Ha']-plan['geometry_reference']['energy_Ha']),
                gradient_max_Ha_Bohr=abs(result['gradient_max_Ha_Bohr']-plan['geometry_reference']['gradient_max_Ha_Bohr']))
            if (any(not math.isfinite(v) for v in errors.values())
                    or not 0<memory.peak<=plan[phase]['maximum_peak_GPU_MiB']
                    or errors['energy_Ha']>plan[phase]['energy_tolerance_Ha']
                    or errors['gradient_max_Ha_Bohr']>plan[phase]['gradient_max_tolerance_Ha_Bohr']):
                raise ValueError('Geometry resource or exact-checkpoint comparison failed: '+repr(errors))
        else:
            execution=api.load_module('split_CPU_energy_adapter',ROOT/'scripts/pipeline/cpu_quantum_execution.py')
            with CPU_resources(scratch,plan[phase],out),execution.cpu_execution(out/'execution',plan[phase]['library_memory_MiB']), (out/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
                result=dft.dft_smd(Path(plan['electronic_input']),0,1,backend='cpu',do_opt=False,do_freq=False)
                atomic_write_json(out/'point.json',result)
                if not result['converged_smd'] or not result['converged_gas']:raise ValueError('CPU active SCF did not converge')
            reference=json.loads(Path(plan['electronic_reference']).read_text())['result']
            errors={k:abs(result[k]-reference[k]) for k in ('e_smd_eV','e_gas_eV')}
            resources=json.loads((out/'resources.json').read_text())['samples']
            if not resources or (out/'resource_stop.json').exists() or (out/'resource_monitor_error.json').exists():
                raise ValueError('CPU resource observation was incomplete or interrupted')
            if any(not math.isfinite(v) or v>plan[phase]['energy_tolerance_eV'] for v in errors.values()):
                raise ValueError('CPU active energies differ from exact GPU reference: '+repr(errors))
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        report.update(status='split_'+phase+'_resource_and_reference_point_passed',differences=errors,
            output_sha256={str(p):api.digest(p) for p in out.rglob('*') if p.is_file() and p.name!='status.json' and p.suffix in ('.json','.npy')})
    except BaseException as error:
        report.update(status='split_resource_point_failed_preserved',error=repr(error));raise
    finally:
        signal.signal(signal.SIGTERM,before);report['finished']=time.time();atomic_write_json(out/'status.json',report)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','geometry','electronic'))
    for name in ('campaign','source','reference','CPU-controls','scratch'):p.add_argument('--'+name,type=Path)
    p.add_argument('--threads',type=int,default=2);a=p.parse_args()
    if a.campaign is None:p.error('--campaign required')
    if a.action=='prepare':
        if not all((a.source,a.reference,a.CPU_controls)):p.error('prepare requires source, reference and CPU controls')
        prepare(a.source.resolve(),a.reference.resolve(),a.CPU_controls.resolve(),a.campaign.resolve())
    else:run(a.campaign.resolve(),a.action,a.threads,a.scratch.resolve() if a.scratch else None)


if __name__=='__main__':main()
