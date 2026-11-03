#!/usr/bin/env python3
"""Drain the full corrected electronic roster with independent CPU collection.

Reuses the bounded shared dispatcher. Failed scientific jobs are preserved;
only proven closed prelaunch refusals can retry, with cooldown and a limit.
This controller does not qualify missing native/ensemble quantities.
"""
import argparse
from collections import Counter
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.sampling_energy import dispatch
from redox.workflows.prelaunch_accounting import prelaunch_accounting,launcher_versions

SOURCE_BYTES=Path(__file__).read_bytes()


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def ready(condition,tries,checked,cooldown,now):
    return condition=='pending' or (condition=='prelaunch_refusal_no_scientific_attempt' and tries<3 and now>=checked+cooldown)


def storage_peers(path,previous):
    """An append-only roster can add newly prepared concurrent calculations."""
    peers=json.loads(path.read_text());current={r['namespace']:r['namespace_sha256'] for r in peers}
    if len(current)!=len(peers) or any(current.get(p)!=h for p,h in previous.items()):
        raise ValueError('Storage peer history cannot be dropped or changed while a queue resumes')
    if any(sha256_file(p)!=h for p,h in current.items()):raise ValueError('Storage peer namespace changed')
    return peers,current


def run(manifest_path,reviews,interpreter,peers_path,output,hours,resume=False):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation,GATE
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    output.mkdir(parents=True,exist_ok=True);owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','sampling_queue_owners')
    observer=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','sampling_queue_status')
    collector=load(ROOT/'scripts/validation/collect_corrected_sampling_energy.py','sampling_queue_collection')
    manifest=json.loads(manifest_path.read_text());rows=manifest['states']
    if manifest['target_states']!=240 or len(rows)!=240 or len({r['key'] for r in rows})!=240:
        raise ValueError('The full 240-input electronic objective is required')
    policy=Path(os.environ.get('GPU_RESERVE_CONFIG',str(GATE.parent.parent/'.config/gpu_reserve/config.env')))
    paths=[manifest_path,*reviews,Path(__file__),policy,
        ROOT/'src/redox/workflows/sampling_energy.py',ROOT/'src/redox/workflows/charged_native_followup.py',
        ROOT/'src/redox/workflows/screening_control.py',ROOT/'scripts/fleet/run_screening_energy_panel.py',
        ROOT/'src/redox/workflows/prelaunch_accounting.py',ROOT/'src/redox/workflows/submission_refusal.py',
        ROOT/'scripts/validation/collect_corrected_sampling_energy.py',ROOT/'src/redox/validation/sampling_receipts.py',
        ROOT/'scripts/validation/review_corrected_sampling_energy.py',ROOT/'src/redox/validation/execution_review.py']
    hashes={str(p):sha256_file(p) for p in paths};snapshots={}
    for p in paths:
        content=SOURCE_BYTES if p==Path(__file__) else p.read_bytes();digest=hashlib.sha256(content).hexdigest()
        frozen=output/'sources'/digest/p.name
        if not frozen.exists():frozen.parent.mkdir(parents=True,exist_ok=True);frozen.write_bytes(content)
        if sha256_file(frozen)!=digest:raise ValueError('Queue source snapshot differs')
        snapshots[str(frozen)]=digest
    with (output/'controller.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);status_path=output/'status.json'
        if status_path.exists():
            previous=json.loads(status_path.read_text())
            if not resume or owners.alive(previous['owner']) is not False:
                raise ValueError('Explicit resume and a proven stopped exact queue owner are required')
            atomic_write_json(output/'history'/f'{time.time_ns()}.json',previous)
        else:previous={}
        owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
        state=dict(status='corrected_sampling_queue_active',owner=owner,CPU_affinity=sorted(os.sched_getaffinity(0)),
            started=time.time(),deadline=time.time()+hours*3600,source_sha256=snapshots,
            declared_inputs=240,native_G_qualified=False,all_scientific_failures_preserved=True,
            accounted_terminal=previous.get('accounted_terminal',[]),host_failures=previous.get('host_failures',{}),
            blocked_hosts=previous.get('blocked_hosts',[]),storage_peers=previous.get('storage_peers',{}))
        stop=threading.Event();old_signal=signal.signal(signal.SIGTERM,lambda *args:stop.set())
        try:
            while not stop.is_set() and time.time()<state['deadline']:
                if any(sha256_file(p)!=h for p,h in hashes.items()):
                    state['status']='queue_source_changed_workers_preserved';break
                peers,state['storage_peers']=storage_peers(peers_path,state['storage_peers'])
                peer_digest=sha256_file(peers_path);peer_snapshot=output/'peer_snapshots'/(peer_digest+'.json')
                if not peer_snapshot.exists():atomic_write_json(peer_snapshot,peers)
                state['storage_peer_roster_snapshot']=str(peer_snapshot)
                entries=[];eligible=set();new_terminal=[];last_retry=state.setdefault('prelaunch_retry_observed',{})
                for row in rows:
                    claim=Path(row['namespace']).parent/'dispatch_claim'
                    condition,data=observer.inspect(row,claim,owners.alive,
                        finished_status='corrected_sampling_energy_finished_independent_review_pending')
                    tries=len(list(claim.glob('try_*.json')));key=row['key']
                    failures=tries;accounting=None;accounting_error=None
                    if condition=='prelaunch_refusal_no_scientific_attempt':
                        try:
                            accounting=prelaunch_accounting(Path(row['namespace']),owners.alive,launcher_versions(output))
                            failures=accounting['prelaunch_failures']
                        except (OSError,ValueError,subprocess.SubprocessError) as error:
                            condition='ownership_unverified';accounting_error=repr(error);data=dict(error=accounting_error)
                    if condition=='prelaunch_refusal_no_scientific_attempt':last_retry.setdefault(key,time.time())
                    if ready(condition,failures,last_retry.get(key,0),60,time.time()):eligible.add(key)
                    if condition in ('completed_review_pending','scientific_failure_preserved','stopped_attempt_preserved','guard_start_failed_preserved') and key not in state['accounted_terminal']:
                        new_terminal.append((data.get('finished',time.time()),key,condition,data.get('owner',data).get('host')))
                    entries.append(dict(key=key,status=condition,submission_receipts=tries,prelaunch_failures=failures,
                                        resource_refusals=accounting['resource_refusals'] if accounting else None,
                                        prelaunch_accounting_error=accounting_error))
                for _,key,condition,host in sorted(new_terminal):
                    state['accounted_terminal'].append(key)
                    if host:
                        state['host_failures'][host]=0 if condition=='completed_review_pending' else state['host_failures'].get(host,0)+1
                        if state['host_failures'][host]>=3 and host not in state['blocked_hosts']:state['blocked_hosts'].append(host)
                state.update(entries=entries,counts=dict(Counter(r['status'] for r in entries)),checked=time.time())
                atomic_write_json(status_path,state)
                summary=collector.collect(manifest_path,manifest_path.parent/'checked_electronic',output/'collections'/f'{time.time_ns()}.json')
                atomic_write_json(output/'electronic_current.json',summary)
                state['checked_electronic_inputs']=summary['checked_electronic_inputs']
                if summary['independent_review_errors']:
                    state.update(status='independent_review_errors_workers_preserved',review_errors=summary['independent_review_errors']);break
                if summary['checked_electronic_inputs']==240:
                    state['status']='all_240_electronic_inputs_independently_checked_native_pending';break
                if stop.is_set():state['status']='queue_stop_requested_workers_preserved';break
                if any(sha256_file(p)!=h for p,h in hashes.items()):
                    state['status']='queue_source_changed_workers_preserved';break
                if eligible:
                    result=dispatch(manifest_path,interpreter,reviews,False,peers,state['blocked_hosts'],eligible)
                    for row in result['launches']:last_retry[row['key']]=time.time()
                    state['last_dispatch']=dict(launches=len(result['launches']),refusals=result['refusals'],checked=time.time())
                accepted={r['key'] for r in summary['states'] if r['status']=='electronic_independently_checked'}
                active=any(r['status'] in ('running','guard_started','submitted','ownership_unverified','uncertain_submission_preserved')
                    or (r['status']=='completed_review_pending' and r['key'] not in accepted) for r in entries)
                if not eligible and not active and not any(r['status']=='prelaunch_refusal_no_scientific_attempt' and r['prelaunch_failures']<3 for r in entries):
                    state['status']='electronic_queue_terminal_incomplete_reviews_or_recovery_required';break
                atomic_write_json(status_path,state);stop.wait(20)
            else:state['status']='queue_stop_requested_workers_preserved' if stop.is_set() else 'queue_window_finished_workers_preserved'
        except BaseException as error:
            state.update(status='queue_controller_failed_workers_preserved',error=repr(error));raise
        finally:
            signal.signal(signal.SIGTERM,old_signal);state['finished']=time.time();atomic_write_json(status_path,state)
    return state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('manifest','storage-peers','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--pilot-review',type=Path,action='append',required=True);p.add_argument('--python',required=True)
    p.add_argument('--hours',type=float,default=24);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if not 0<a.hours<=72:p.error('A positive bounded queue window up to 72 hours is required')
    result=run(a.manifest.resolve(),[p.resolve() for p in a.pilot_review],a.python,a.storage_peers.resolve(),a.out.resolve(),a.hours,a.resume)
    print(json.dumps(dict(status=result['status'],checked_electronic_inputs=result.get('checked_electronic_inputs',0))))
