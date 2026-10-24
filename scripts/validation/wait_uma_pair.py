#!/usr/bin/env python
"""Launch one fixed-device FP32/Nitro timing/point comparison when an H100 is free.

Reservation contention consumes no validation attempt. A started scientific
failure is preserved for review. Unreachable owners are never presumed dead.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,free_slots


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--raw-out',type=Path,required=True);p.add_argument('--max-hours',type=float,default=48)
    p.add_argument('--quiet-seconds',type=float,default=30)
    a=p.parse_args()
    if not 0<a.max_hours<=48 or not 0<=a.quiet_seconds<=60:p.error('Use a finite window and 0–60 second sampled quiet placement')
    with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    a.out.mkdir(parents=True,exist_ok=False)
    policy_path=ROOT/'config/uma_execution.json';policy=json.loads(policy_path.read_text())
    script=ROOT/'scripts/validation/qualify_uma_current.py'
    watched={str(path):api.digest(path) for path in (Path(__file__).resolve(),script,ROOT/'src/redox/qm/uma_current.py',policy_path,a.campaign/'plan.json')}
    owner=dict(host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),deadline=time.time()+a.max_hours*3600,
        source_sha256=watched,scope='One serial fixed-device reference/Nitro validation baseline; not a throughput fleet')
    atomic_write_json(a.out/'owner.json',owner)
    status=dict(status='waiting_for_free_compatible_H100',actual_scientific_starts=0,gate_races=0,residency_interruptions=0)
    alive=api.load_module('UMA_pair_owner_reader',ROOT/'scripts/pipeline/run_polymer_spacing.py').owner_alive
    quiet_since={}
    while time.time()<owner['deadline']:
        api.verify_sources(watched)
        # Check imports before watching for a short free-GPU window. CPU/RAM
        # assignments are checked again atomically by the remote launcher.
        host='lambda13';cached=a.out/f'{host}_preflight.json'
        cache=json.loads(cached.read_text()) if cached.exists() else {}
        fresh=cached.exists() and time.time()-cached.stat().st_mtime<600
        if not (fresh and cache.get('ok') and cache.get('python')==policy['nitro_python']):
            helper=ROOT.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
            pre=subprocess.run([sys.executable,str(helper),host,'--python',policy['nitro_python'],
                '--module','torch','--module','ase','--module','rdkit','--module','psutil',
                '--library-path',policy['nitro_library_path'],'--path',str(a.campaign/'plan.json'),
                '--min-available-ram-mib','16384'],capture_output=True,text=True)
            cached.write_text(pre.stdout)
            if pre.returncode:
                status.update(status='H100_environment_preflight_failed_no_scientific_attempt',error=pre.stderr)
                break
        try:slots=free_slots(['lambda13'])
        except (OSError,subprocess.SubprocessError,RuntimeError) as error:
            status.update(status='availability_unknown_no_restart',error=repr(error));slots=[]
        now=time.time();quiet_since={slot:quiet_since.get(slot,now) for slot in slots}
        slots=[slot for slot in slots if now-quiet_since[slot]>=a.quiet_seconds]
        if not slots:
            status['sampled_quiet_seconds_required']=a.quiet_seconds
            status['checked']=time.time();atomic_write_json(a.out/'status.json',status);time.sleep(10);continue
        host,gpu=slots[-1]
        attempt=status['gate_races']+status['residency_interruptions']+1
        raw_target=a.raw_out/f'attempt_{attempt:03d}'
        log=a.out/f'worker_{attempt:03d}.log'
        argv=[str(ROOT.parent/'bin/gpu_reserve'),'run',str(gpu),'--label','redox-official-FP32-Nitro-pair','--','env',
            *[key+'=2' for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')],
            'PYTHONPATH='+str(ROOT/'src'),'LD_LIBRARY_PATH='+policy['nitro_library_path'],
            'FAIRCHEM_CACHE_DIR='+policy['fairchem_cache'],policy['nitro_python'],
            str(ROOT/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(a.out/f'guard_{attempt:03d}'),
            str(script),'pair','--campaign',str(a.campaign.resolve()),'--out',str(raw_target.resolve())]
        code='''import json,os,subprocess,time
from pathlib import Path
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():
 host_capacity(2,minimum_ram_gib=16);cpus=cpu_allocation(2)
 argv=['taskset','--cpu-list',','.join(map(str,cpus))]+%r
 with open(%r,'xb') as log:p=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,cwd=%r)
print(json.dumps(dict(host=os.uname().nodename,pid=p.pid,process_start=Path('/proc/%%s/stat'%%p.pid).read_text().rsplit(')',1)[1].split()[19],CPU_affinity=cpus,argv=argv,started=time.time())))
'''%(argv,str(log),str(ROOT))
        try:
            launch=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,
                shlex.join(['env','PYTHONPATH='+str(ROOT/'src'),'LD_LIBRARY_PATH='+policy['nitro_library_path'],policy['nitro_python'],'-c',code])],capture_output=True,text=True,timeout=35,check=True)
            handle=json.loads(launch.stdout)
            handle.update(preflight=str(cached),preflight_sha256=api.digest(cached))
        except (OSError,subprocess.SubprocessError,ValueError) as error:
            status.update(status='launch_observation_unresolved_preserve_do_not_retry',error=repr(error));break
        atomic_write_json(a.out/f'launch_{attempt:03d}.json',handle)
        status.update(status='gate_started_actual_prediction_pending',handle=handle,raw_pair_output=str(raw_target),checked=time.time())
        atomic_write_json(a.out/'status.json',status)
        running=True
        started_noted=False
        while time.time()<owner['deadline']:
            try:running=alive(handle)
            except (OSError,RuntimeError,subprocess.SubprocessError) as error:
                status.update(status='owned_gate_observation_unknown_no_restart',error=repr(error))
                atomic_write_json(a.out/'status.json',status);time.sleep(10);continue
            if not running:break
            if (raw_target/'reference/status.json').exists() and not started_noted:
                status.update(status='matched_pair_actually_started',actual_scientific_starts=status['actual_scientific_starts']+1)
                started_noted=True
            status['checked']=time.time();atomic_write_json(a.out/'status.json',status);time.sleep(10)
        if running:
            status['status']='finite_observer_ended_owned_validation_continues_no_restart';break
        interruption=a.out/f'guard_{attempt:03d}'/'interruption.json'
        guard_state=a.out/f'guard_{attempt:03d}'/'state.json'
        initial_conflict=json.loads(guard_state.read_text()).get('foreign_resident_PIDs',[]) if guard_state.exists() else []
        if interruption.exists() or initial_conflict:
            status.update(status='owned_residency_interruption_preserved_retry_resources',
                residency_interruptions=status['residency_interruptions']+1,
                last_interruption=json.loads(interruption.read_text()) if interruption.exists() else json.loads(guard_state.read_text()))
            quiet_since.pop((host,gpu),None)
            atomic_write_json(a.out/'status.json',status)
            if status['residency_interruptions']>=3:
                status['status']='residency_interruptions_require_placement_review';break
            time.sleep(10);continue
        # A remote owner's exit is authoritative, but NFS publication may lag it.
        for _ in range(5):
            if (raw_target/'status.json').exists():break
            time.sleep(2)
        if (raw_target/'status.json').exists():
            status.update(status='owned_pair_ended_read_retained_outcome',
                outcome=json.loads((raw_target/'status.json').read_text()))
            break
        if log.exists() and any(word in log.read_text() for word in ('BUSY','reserved')):
            status.update(status='gate_contention_no_scientific_attempt',gate_races=status['gate_races']+1);time.sleep(10);continue
        status['status']='owned_start_ended_without_pair_output_review_required';break
    else:status['status']='finite_resource_wait_ended_no_validation_start'
    status['finished']=time.time();atomic_write_json(a.out/'status.json',status)
    owner['finished']=time.time();atomic_write_json(a.out/'owner.json',owner)


if __name__=='__main__':main()
