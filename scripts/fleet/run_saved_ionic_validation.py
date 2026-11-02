#!/usr/bin/env python3
"""Wait for suitable capacity and dispatch one prepared saved-ionic-Hessian check.

Uses the maintained common probe, reservation, CPU/RAM and storage dispatcher.
An uncertain or scientific failure is terminal here; no calculation is retried.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify,load
from redox.workflows.charged_native_followup import dispatch_rows
from redox.workflows.sampling_native_queue import extend_peers


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--storage-peers',type=Path,required=True);parser.add_argument('--python',required=True)
    parser.add_argument('--hours',type=float,default=12)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not 0<args.hours<=72:raise ValueError('A positive bounded dispatch window of at most 72 hours is required')
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    namespace=args.spec.resolve();spec=verify(namespace);case=Path(spec['case'])
    plan=json.loads((case/'plan.json').read_text())
    if not plan.get('charged_saved_Hessian_validation') or not plan.get('saved_full_Hessian'):
        raise ValueError('Only a prepared source-bound ionic saved-matrix validation is admitted')
    output=args.out.resolve();output.mkdir(parents=True,exist_ok=args.resume)
    row=dict(key=plan['prepared_case']['key'],namespace=str(namespace),namespace_sha256=sha256_file(namespace))
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','saved_ionic_owners')
    remote=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','saved_ionic_common_inspection')
    sources=[Path(__file__).resolve(),namespace,ROOT/'src/redox/workflows/charged_native_followup.py',
        ROOT/'src/redox/workflows/screening_control.py',ROOT/'scripts/fleet/run_screening_energy_panel.py',
        ROOT/'scripts/pipeline/run_namespace_energy_queue.py',ROOT/'scripts/fleet/run_charged_polymer_geometry_panel.py']
    hashes={str(p):sha256_file(p) for p in sources}
    state=dict(status='saved_ionic_validation_controller_started',started=time.time(),deadline=time.time()+args.hours*3600,
        owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
                   process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        CPU_affinity=sorted(os.sched_getaffinity(0)),source_sha256=hashes,namespace=str(namespace),
        native_G_qualified=False,source_failure_preserved=True)
    with (output/'controller.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        status_path=output/'status.json'
        if status_path.exists():
            previous=json.loads(status_path.read_text())
            if not args.resume or owners.alive(previous['owner']) is not False:
                raise ValueError('Explicit resume requires the exact previous controller to be closed')
            history=output/'history';history.mkdir(exist_ok=True)
            atomic_write_json(history/(str(time.time_ns())+'.json'),previous)
        atomic_write_json(status_path,state)
        peers=extend_peers(args.storage_peers.resolve(),[row])
        while time.time()<state['deadline']:
            if any(sha256_file(p)!=h for p,h in hashes.items()):
                state['status']='source_changed_existing_attempt_preserved';break
            condition,detail=remote.inspect(row,namespace.parent/'dispatch_claim',owners.alive,
                finished_status='charged_polymer_native_component_finished_review_pending')
            state.update(checked=time.time(),condition=condition,attempt_observation=detail)
            if condition=='completed_review_pending':state['status']='saved_ionic_validation_finished_independent_review_pending';break
            if condition in ('scientific_failure_preserved','stopped_attempt_preserved','guard_start_failed_preserved','uncertain_submission_preserved'):
                state['status']='saved_ionic_attempt_preserved_for_explicit_review';break
            if condition in ('pending','prelaunch_refusal_no_scientific_attempt'):
                tries=list((namespace.parent/'dispatch_claim').glob('try_*.json'))
                if len(tries)>=3:state['status']='bounded_prelaunch_refusals_preserved';break
                peers=extend_peers(args.storage_peers.resolve(),[row])
                dispatch=dispatch_rows([row],output,args.python,sources,storage_peers=peers)
                state['last_dispatch']=str(output/'followup_dispatch_current.json')
                state['status']=('saved_ionic_validation_submitted_start_unverified' if dispatch['launches']
                                 else 'waiting_for_suitable_free_GPU_or_host_capacity')
            else:state['status']='saved_ionic_validation_attempt_observed'
            atomic_write_json(output/'status.json',state)
            time.sleep(30)
        else:state['status']='dispatch_window_expired_attempt_preserved'
        state['finished']=time.time();atomic_write_json(output/'status.json',state)


if __name__=='__main__':main()
