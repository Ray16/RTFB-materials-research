#!/usr/bin/env python
"""Wait for one qualified GPU and start a reviewed gas counterpoise point.

The coordinator never retries an uncertain launch and never modifies a finished
scientific attempt. All numerical work stays in the source-bound point worker.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT=Path(__file__).resolve()
PROFILE=ROOT/'environments/cutensor22/runtime.json'
INVENTORY=ROOT/'scripts/pipeline/coordinate_native_inventory_resilient.py'
PREPARE=ROOT/'scripts/pipeline/prepare_local_counterpoise_point.py'
CPU_PREFLIGHT=ROOT/'scripts/fleet/preflight_cpu_dependencies.py'


def parse_memory(output):
    lines=output.strip().splitlines()
    if len(lines)!=1:raise ValueError('GPU memory response must contain one value')
    value=float(lines[0])
    if not 0<value<10**7:raise ValueError('Invalid GPU memory response')
    return value


def classify_launch_response(response):
    """Only an explicit non-start can be retried; unknown outcomes are held."""
    if not isinstance(response,dict):return 'preserve'
    status=response.get('status')
    if status=='gate_started_scientific_start_pending' and isinstance(response.get('owner'),dict):return 'started'
    if status in ('atomic_placement_unavailable','resource_capacity_unavailable'):return 'not_started'
    return 'preserve'


def device_memory(host,gpu):
    remote=['nvidia-smi',f'--id={gpu}','--query-gpu=memory.total','--format=csv,noheader,nounits']
    result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(remote)],
        capture_output=True,text=True,timeout=20,check=True)
    return parse_memory(result.stdout)


def preflight(case,host,profile):
    path=case/(host+'_preflight.json')
    if path.exists() and time.time()-path.stat().st_mtime<600:
        cached=json.loads(path.read_text())
        if cached.get('ok') and cached.get('python')==profile['python'] and cached.get('GPU_imports_blocked_in_CPU_preflight'):
            return path
    command=[sys.executable,str(CPU_PREFLIGHT),host,'--python',profile['python'],
        '--module','pyscf','--module','gpu4pyscf','--module','cupy','--module','ase','--module','scipy',
        '--path',str(case/'plan.json'),'--library-path',profile['library_path'],'--min-available-ram-mib','65536']
    result=subprocess.run(command,capture_output=True,text=True,timeout=60)
    data=json.loads(result.stdout);atomic_write_json(path,data)
    if (result.returncode or not data.get('ok') or data.get('python')!=profile['python']
            or not data.get('GPU_imports_blocked_in_CPU_preflight') or data.get('GPU_modules_loaded')):
        raise ValueError('CPU-safe node preflight failed')
    return path


def validate_geometry_review(geometry,case):
    if (geometry.get('status')!='prepared_motif_counterpoise_geometry_independently_checked_no_scientific_start'
            or geometry.get('campaign')!=str(case) or geometry.get('primary_candidate_admission')
            or geometry.get('solvent_counterpoise_computed')):
        raise ValueError('An independently checked exact prepared geometry is required')


def run(case,geometry_review,hosts,hours):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    profile=json.loads(PROFILE.read_text());plan=api.read_plan(case);campaign=json.loads((case.parent/'campaign.json').read_text())
    if Path(sys.executable)!=Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH'):
        raise ValueError('Use the declared interpreter and host loader environment')
    if not 0<hours<=48 or plan['minimum_GPU_MiB']!=32768 or plan['primary_candidate_admission'] or plan['solvent_counterpoise_computed']:
        raise ValueError('Invalid time, device class or scientific scope')
    if len(campaign['entries'])!=1 or campaign['entries'][0]['plan_sha256']!=api.digest(case/'plan.json'):
        raise ValueError('Prepared campaign does not match the exact point')
    api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    control=json.loads(Path(plan['control_review']).read_text());api.verify_sources(control['source_sha256'])
    if control['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':raise ValueError('Numerical control is not checked')
    geometry=json.loads(geometry_review.read_text());api.verify_sources(geometry['source_sha256'])
    validate_geometry_review(geometry,case)
    inventory=api.load_module('counterpoise_inventory_contract',INVENTORY)
    if len(hosts)!=len(set(hosts)) or not set(hosts)<=inventory.permitted_hosts(inventory.POLICY):
        raise ValueError('Hosts do not match the current GPU reservation policy')
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    with (case/'coordinator.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        status=case/'coordinator.json';intent=case/'dispatch_intent.json'
        if status.exists() or intent.exists() or (case/'launch.json').exists() or (case/'counterpoise').exists():
            raise ValueError('Existing ownership or attempt must be inspected before any restart')
        canonical_preflight=ROOT.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
        sources={str(p):api.digest(p) for p in [SCRIPT,INVENTORY,PREPARE,CPU_PREFLIGHT,canonical_preflight,
            inventory.POLICY,inventory.GATE,PROFILE,geometry_review,case/'plan.json',case.parent/'campaign.json']}
        record=dict(status='waiting_for_qualified_32GB_GPU',owner=owner,started=time.time(),deadline=time.time()+hours*3600,
            CPU_affinity=sorted(os.sched_getaffinity(0)),scientific_starts=0,source_sha256=sources,
            primary_candidate_admission=False,solvent_counterpoise_computed=False)
        failures={}
        try:
            while time.time()<record['deadline']:
                api.verify_sources(sources);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
                slots,checks=inventory.query_all(hosts)
                record.update(updated=time.time(),inventory=checks,available_slots=slots);atomic_write_json(status,record)
                for host,gpu in slots:
                    if time.time()<failures.get(host,0):continue
                    try:
                        memory=device_memory(host,gpu)
                        if memory<32768:continue
                        preflight_path=preflight(case,host,profile)
                    except (OSError,ValueError,subprocess.SubprocessError) as error:
                        failures[host]=time.time()+120
                        record.setdefault('preflight_failures',[]).append(dict(host=host,GPU=gpu,error=repr(error),checked=time.time()))
                        atomic_write_json(status,record);continue
                    atomic_write_json(intent,dict(host=host,GPU=gpu,owner=owner,checked=time.time(),outcome='requested'))
                    command=['env','OMP_NUM_THREADS=1','OMP_THREAD_LIMIT=1','MKL_NUM_THREADS=1',
                        'OPENBLAS_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','XTB_THREADS=1','PYTHONPATH='+str(ROOT/'src'),
                        profile['python'],str(PREPARE),'launch','--case',str(case),'--gpu',str(gpu),'--preflight',str(preflight_path)]
                    try:
                        response=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(command)],
                            capture_output=True,text=True,timeout=45,check=True)
                        launch=json.loads(response.stdout)
                    except (OSError,ValueError,subprocess.SubprocessError) as error:
                        record.update(status='submission_ownership_unverified_preserved',error=repr(error));return
                    action=classify_launch_response(launch)
                    if action=='started':
                        record.update(status='reserved_worker_started_scientific_confirmation_pending',launch=launch)
                        atomic_write_json(status,record)
                        for _ in range(12):
                            state_path=case/'counterpoise/status.json';guard_path=case/'residency_guard/state.json'
                            if state_path.exists() and guard_path.exists():
                                state=json.loads(state_path.read_text());guard=json.loads(guard_path.read_text())
                                if state.get('status')=='gas_counterpoise_control_running' and guard.get('own_residency_observed'):
                                    owners=api.load_module('counterpoise_point_owner_check',ROOT/'scripts/pipeline/run_polymer_spacing.py')
                                    if owners.owner_alive(state['owner']):
                                        record.update(status='scientific_counterpoise_started',scientific_starts=1,scientific_owner=state['owner']);return
                            time.sleep(5)
                        record['status']='launched_scientific_ownership_requires_review';return
                    if action=='preserve':
                        record.update(status='unexpected_launch_response_preserved',launch=launch);return
                    atomic_write_json(intent,dict(host=host,GPU=gpu,owner=owner,checked=time.time(),outcome='not_started'))
                time.sleep(15)
            record['status']='wait_window_ended_no_scientific_attempt'
        except BaseException as error:
            record.update(status='coordinator_failed_existing_attempt_preserved',error=repr(error));raise
        finally:
            record['finished']=time.time();atomic_write_json(status,record)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True);p.add_argument('--geometry-review',type=Path,required=True);p.add_argument('--hosts',required=True);p.add_argument('--max-hours',type=float,default=24);a=p.parse_args()
    run(a.case.resolve(),a.geometry_review.resolve(),a.hosts.split(','),a.max_hours)
