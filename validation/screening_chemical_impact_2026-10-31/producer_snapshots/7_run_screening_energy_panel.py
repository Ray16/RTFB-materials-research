#!/usr/bin/env python3
"""Dispatch the frozen benchmark controls across suitable policy-approved GPUs.

Two largest RKS/UKS inputs run first. Wider dispatch requires a separately
checked resource-pilot receipt. Scientific failures and uncertain starts remain
preserved; no claim is reclaimed by age, and no native G is inferred.
"""
import argparse
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def launch(row,host,gpu,interpreter,claim,try_number):
    namespace=Path(row['namespace']);root=namespace.parent
    output=root/f'placement_{try_number:03d}'
    # SSH gives the allocator a fresh node affinity, independently of this
    # coordinator's one-CPU binding. No GPU work runs in the SSH client.
    code='''import json,os,pathlib,subprocess,sys
from redox.workflows.screening_control import probe
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
spec=pathlib.Path(sys.argv[1]);output=pathlib.Path(sys.argv[2]);gpu=sys.argv[3]
allowed=os.sched_getaffinity(0)
with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
probe(spec)
os.sched_setaffinity(0,allowed)
argv=[sys.executable,"-u","-m","redox.workflows.screening_control","launch","--spec",str(spec),"--gpu",gpu,"--out",str(output)]
logpath=pathlib.Path(sys.argv[4]);logpath.parent.mkdir(parents=True,exist_ok=True)
with logpath.open("x") as log:child=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
owner=dict(host=os.uname().nodename.split(".")[0],pid=child.pid,process_start=pathlib.Path(f"/proc/{child.pid}/stat").read_text().rsplit(")",1)[1].split()[19])
receipt=dict(status="controller_spawned_scientific_start_unverified",owner=owner,placement=str(output),argv=argv)
atomic_write_json(pathlib.Path(sys.argv[5]),receipt);print(json.dumps(receipt))
'''
    log=root/'dispatch_logs'/f'try_{try_number:03d}.log'
    receipt=claim/f'try_{try_number:03d}.json'
    argv=['env','PYTHONPATH='+str(ROOT/'src'),'OMP_NUM_THREADS=1','OMP_THREAD_LIMIT=1',
        'OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','XTB_THREADS=1',
        'NUMPY_MADVISE_HUGEPAGE=0',interpreter,'-c',code,str(namespace),str(output),str(gpu),str(log),str(receipt)]
    result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(argv)],
                          capture_output=True,text=True,timeout=180,check=True)
    return json.loads(result.stdout)


def inspect(row,claim,alive):
    root=Path(row['namespace']).parent;case=root/'case';stage=case/'status.json'
    if stage.exists():
        data=json.loads(stage.read_text())
        if 'finished' in data:
            return ('completed_review_pending' if data['status']=='screening_energy_control_finished_review_pending'
                    else 'scientific_failure_preserved'),data
        live=alive(data['owner'])
        return ('running' if live is True else 'stopped_attempt_preserved' if live is False else 'ownership_unverified'),data
    guard=case/'residency_guard/owner.json'
    if guard.exists():
        data=json.loads(guard.read_text());live=alive(data)
        return ('guard_started' if live is True else 'guard_start_failed_preserved' if live is False else 'ownership_unverified'),data
    if not claim.exists():return 'pending',{}
    tries=sorted(claim.glob('try_*.json'))
    if not tries:return 'uncertain_submission_preserved',{}
    receipt=json.loads(tries[-1].read_text());live=alive(receipt['owner'])
    if live is not False:return ('submitted' if live else 'ownership_unverified'),receipt
    placement=Path(receipt['placement'])
    if (placement/'launch.json').exists():
        wrapper=json.loads((placement/'launch.json').read_text())
        if alive(dict(host=wrapper['host'],pid=wrapper['pid'],process_start=wrapper['process_start'])) is not False:
            return 'submitted',receipt
    # A source-frozen case without a scientific stage or guard can be retried
    # only after both exact controllers are gone. Preserve every placement.
    return 'prelaunch_refusal_no_scientific_attempt',receipt


def coordinate(manifest_path,output,interpreter,hours,resume):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    utilities=load(ROOT/'scripts/pipeline/run_namespace_energy_queue.py','screening_panel_capacity')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','screening_panel_owners')
    manifest=json.loads(manifest_path.read_text());rows=manifest['states']
    if len(rows)!=manifest['target_states'] or len({r['key'] for r in rows})!=len(rows):
        raise ValueError('Distinct complete prepared state roster required')
    for row in rows:
        if sha256_file(row['namespace'])!=row['namespace_sha256']:raise ValueError('Prepared namespace changed')
        verify(Path(row['namespace']))
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    output.mkdir(parents=True,exist_ok=True)
    pilot_keys={max([r for r in rows if (r['source_state']['mult']>1)==open_shell],
                    key=lambda r:(r['nao'],-abs(r['source_state']['charge'])))['key'] for open_shell in (False,True)}
    hosts_path=ROOT.parent/'.config/gpu_reserve/config.env'
    # The project utility reads the same authoritative gate policy.
    import re
    match=re.search(r'^ALLOWED_HOSTS="([^"]+)"',hosts_path.read_text(),re.M)
    if match is None:raise ValueError('Cannot read the system GPU policy')
    hosts=match.group(1).split()
    with (output/'controller.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state_path=output/'status.json'
        if state_path.exists():
            old=json.loads(state_path.read_text())
            if not resume or owners.alive(old['owner']) is not False:raise ValueError('Proven stopped owner and explicit resume required')
            history=output/'history';history.mkdir(exist_ok=True);atomic_write_json(history/(str(time.time_ns())+'.json'),old)
        source_paths=[Path(__file__),manifest_path,hosts_path,ROOT/'scripts/pipeline/run_namespace_energy_queue.py',
                      ROOT/'scripts/pipeline/observe_archived_native_namespace.py']
        hashes={str(p):sha256_file(p) for p in source_paths}
        owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
        status=dict(status='screening_energy_panel_running',owner=owner,started=time.time(),deadline=time.time()+hours*3600,
                    source_sha256=hashes,pilot_keys=sorted(pilot_keys),native_G_qualified=False)
        blocked_hosts=set();consecutive_failures={};accounted=set();cooldown={}
        try:
            while time.time()<status['deadline']:
                if any(sha256_file(p)!=h for p,h in hashes.items()):raise ValueError('Dispatcher or prepared manifest changed')
                qualification=manifest_path.parent/'resource_pilots.json'
                released=False
                if qualification.exists():
                    receipt=json.loads(qualification.read_text())
                    released=(receipt['status']=='largest_RKS_UKS_energy_control_resources_independently_checked'
                              and set(receipt['keys'])==pilot_keys and all(sha256_file(p)==h for p,h in receipt['source_sha256'].items()))
                    if not released:raise ValueError('Changed or unsupported pilot qualification receipt')
                entries=[];pending=[]
                for row in rows:
                    claim=Path(row['namespace']).parent/'dispatch_claim'
                    condition,data=inspect(row,claim,owners.alive)
                    entries.append(dict(key=row['key'],status=condition))
                    if condition in ('scientific_failure_preserved','stopped_attempt_preserved','guard_start_failed_preserved') and row['key'] not in accounted:
                        host=data.get('owner',data).get('host')
                        if host:
                            consecutive_failures[host]=consecutive_failures.get(host,0)+1
                            if consecutive_failures[host]>=3:blocked_hosts.add(host)
                        accounted.add(row['key'])
                    if condition=='completed_review_pending' and row['key'] not in accounted:
                        consecutive_failures[data['owner']['host']]=0;accounted.add(row['key'])
                    if condition in ('pending','prelaunch_refusal_no_scientific_attempt') and (released or row['key'] in pilot_keys) and time.time()>=cooldown.get(row['key'],0):pending.append(row)
                status.update(updated=time.time(),entries=entries,resource_pilots_released=released,blocked_hosts=sorted(blocked_hosts))
                atomic_write_json(state_path,status)
                if all(e['status'] in ('completed_review_pending','scientific_failure_preserved','stopped_attempt_preserved','guard_start_failed_preserved','uncertain_submission_preserved') for e in entries):
                    status['status']='screening_energy_panel_terminal_independent_reviews_required';break
                if pending:
                    slots,errors=utilities.discover_slots([h for h in hosts if h not in blocked_hosts]);status['availability_errors']=errors
                    for host,gpu in slots:
                        if not pending:break
                        row=pending[0]
                        try:capacity=utilities.capacity(host,gpu,interpreter)
                        except subprocess.SubprocessError:continue
                        if not utilities.fits(dict(maximum_host_RSS_GiB=40),capacity):continue
                        pending.pop(0);claim=Path(row['namespace']).parent/'dispatch_claim'
                        if not claim.exists():claim.mkdir();atomic_write_json(claim/'owner.json',dict(owner=owner,queue=str(output)))
                        else:
                            old=json.loads((claim/'owner.json').read_text())
                            if old['queue']!=str(output):raise ValueError('Foreign dispatch claim preserved')
                        number=len(list(claim.glob('try_*.json')))+1
                        try:launch(row,host,gpu,interpreter,claim,number)
                        except (OSError,subprocess.SubprocessError) as error:
                            atomic_write_json(claim/f'uncertain_{number:03d}.json',dict(host=host,GPU=gpu,error=repr(error),checked=time.time()))
                            blocked_hosts.add(host)
                        cooldown[row['key']]=time.time()+30
                time.sleep(5)
            else:status['status']='screening_energy_dispatch_window_ended_workers_preserved'
        except BaseException as error:
            status.update(status='screening_energy_dispatch_failed_workers_preserved',error=repr(error));raise
        finally:
            status['finished']=time.time();atomic_write_json(state_path,status)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--python',required=True);parser.add_argument('--hours',type=float,default=12);parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not 0<args.hours<=24:parser.error('Use a positive bounded 24-hour dispatch window')
    coordinate(args.manifest.resolve(),args.out.resolve(),args.python,args.hours,args.resume)
