#!/usr/bin/env python3
"""Complete one source-bound native state inside its private archive namespace."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.core.protocol import load_record,record_path
    from redox.qm.native_thermal import ensure_native_thermal,select_policy,load_native_thermal,free_energy_native,cache_path
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    verify_gate_owner();plan=json.loads(a.plan.read_text());case=a.plan.parent
    if any(sha(path)!=digest for path,digest in plan['source_sha256'].items()):raise ValueError('Prepared molecular source changed')
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or {os.environ.get(k) for k in caps}!={'2'}:
        raise ValueError('Two disjoint CPUs and six thread caps=2 required')
    host_capacity(2,minimum_ram_gib=plan['maximum_host_RSS_GiB']+32)
    gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    guard=json.loads((case/'host_guard/status.json').read_text())
    if guard['pid']!=os.getpid() or guard['maximum_RSS_GiB']!=plan['maximum_host_RSS_GiB']:
        raise ValueError('The prospective same-PID host guard is required')
    if not plan['route_gate']['proceed_to_full_redox'] or not plan['route_gate']['exact_route_binding_verified']:
        raise ValueError('Exact prospective route admission required')
    owner=dict(host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status=dict(status='molecular_native_running_unqualified',owner=owner,started=time.time(),
        plan_sha256=sha(a.plan),candidate_alias=plan['candidate_alias'],state=plan['state'],
        native_G_qualified=False,material_admission=False,prospective_resource_pilot=True)
    atomic_write_json(case/'status.json',status)
    def stop(signum,frame):raise KeyboardInterrupt('Owned molecular native stop; retain SCF/geometry checkpoints')
    signal.signal(signal.SIGTERM,stop)
    try:
        checkpoint=load(case.parent/'sources/checkpoints.py','molecular_native_checkpoint_adapter')
        from redox.qm import dft
        dft._build_smd_mf=checkpoint.recording_builder(dft._build_smd_mf,case/'SCF_evaluations',owner)
        from pyscf.geomopt import geometric_solver
        geometric_solver.kernel=checkpoint.recording_optimizer(geometric_solver.kernel,case/'optimizer_returns',owner)
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','molecular_native_memory_monitor')
        directory=Path(plan['target_directory']);gid,state=plan['source_id'],plan['state']
        select_policy(directory)
        with monitor.MemoryMonitor() as memory:
            try:kind,thermal=ensure_native_thermal(gid,state)
            except RuntimeError as error:
                if not plan['allow_bounded_native_minimum_repair'] or not any(s in str(error) for s in ('not stationary','not a minimum')):raise
                from redox.qm.native_thermal_repair import repair_minimum
                kind,thermal=repair_minimum(gid,state,str(error))
            energy=load_record(directory,plan['charge'],plan['mult'])
            thermal=load_native_thermal(directory,plan['charge'],plan['mult'])
            total,quality=free_energy_native(gid,state)
        if (energy is None or thermal is None or total is None or quality!='ok'
                or not 0<memory.peak<=plan['maximum_peak_GPU_MiB']):
            raise ValueError('Native thermochemistry/active-energy/resource admission failed')
        surface,_=cache_path(directory,plan['charge'],plan['mult'])
        by_cutoff={str(c):energy['e_smd_eV']+thermal['nuclear_by_cutoff_eV'][str(c)]
                   +thermal['g_spin_partition_eV'] for c in (25,50,100)}
        result=dict(status='molecular_native_state_finished_independent_review_pending',candidate_alias=plan['candidate_alias'],
            source_id=gid,state=state,charge=plan['charge'],mult=plan['mult'],product_inchikey=plan['product_inchikey'],
            G_eV=total,G_by_cutoff_eV=by_cutoff,kind=kind,active_E_smd_eV=energy['e_smd_eV'],
            native_record_path=str(surface),active_record_path=str(record_path(directory,plan['charge'],plan['mult'])),
            final_geometry_path=str(directory/'opt.xyz'),final_geometry_sha256=sha(directory/'opt.xyz'),
            native_G_qualified=False,observed_peak_GPU_MiB=memory.peak,
            global_conformer_convergence_proven=False,material_admission=False,
            source_sha256={str(p):sha(p) for p in (surface,record_path(directory,plan['charge'],plan['mult']),
                directory/'opt.xyz',directory/'result.json',directory/'thermal_policy.json')})
        atomic_write_json(case/'result.json',result)
        status.update(status=result['status'],result_sha256=sha(case/'result.json'),observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:
        status.update(status='molecular_native_failed_preserved',error=repr(error));raise
    finally:
        status['finished']=time.time();atomic_write_json(case/'status.json',status)


if __name__=='__main__':main()
