#!/usr/bin/env python3
"""Run the archived frozen-slow solution calculation on an admitted exact pair.

This worker executes inside the read-only original-path namespace. Scientific
equations and SCF controls come from the source-bound archived implementation.
Completed SCF/density chunks remain on shared storage if execution is stopped.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    verify_gate_owner()
    plan=json.loads(args.plan.read_text());case=args.plan.parent
    if any(sha(path)!=digest for path,digest in plan['source_sha256'].items()):
        raise ValueError('Exact solution-pair source or input changed')
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or {os.environ.get(k) for k in caps}!={'2'}:
        raise ValueError('Two disjoint CPUs and thread caps=2 required')
    host_capacity(2,minimum_ram_gib=plan['maximum_host_RSS_GiB']+32)
    gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    guard=json.loads((case/'host_guard/status.json').read_text())
    if guard['pid']!=os.getpid() or guard['maximum_RSS_GiB']!=plan['maximum_host_RSS_GiB']:
        raise ValueError('The declared same-PID resource guard must be active')
    if plan['route_gate']['route_status']!='proceed' or not plan['route_gate']['exact_route_binding_verified']:
        raise ValueError('An exact reviewed route is required before the property calculation')
    status=dict(status='solution_pair_running_unqualified',started=time.time(),
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        plan_sha256=sha(args.plan),candidate_id=plan['candidate_id'],prospective_resource_pilot=True,
        solution_descriptor_qualified=False,material_admission=False)
    atomic_write_json(case/'status.json',status)
    def stop(signum,frame):raise KeyboardInterrupt('Owned solution pair stopped; retain all SCF checkpoints')
    signal.signal(signal.SIGTERM,stop)
    try:
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','solution_pair_device_monitor')
        control=load(ROOT/'scripts/validation/solution_reorg_campaign.py','solution_pair_backend_control')
        control.OUT=case/'backend_control';control.OUT.mkdir(exist_ok=False)
        with monitor.MemoryMonitor() as memory:
            control.independent_reference('gpu')
            pair=load(ROOT/'scripts/validation/solution_reorg_pair.py','source_bound_solution_pair')
            pair.evaluate(plan['pair'],case/'pair_results','gpu',None,None)
        if not 0<memory.peak<=plan['maximum_peak_GPU_MiB']:
            raise ValueError('Solution-pair device-memory pilot failed')
        if any(sha(path)!=digest for path,digest in plan['source_sha256'].items()):
            raise ValueError('Solution-pair inputs changed during evaluation')
        status.update(status='solution_pair_finished_independent_review_pending',observed_peak_GPU_MiB=memory.peak,
            result_sha256=sha(case/'pair_results/result.json'),backend_control_sha256=sha(case/'backend_control/native_smd_gpu_reference.json'))
    except BaseException as error:
        status.update(status='solution_pair_failed_preserved',error=repr(error));raise
    finally:
        status['finished']=time.time();atomic_write_json(case/'status.json',status)


if __name__=='__main__':main()
