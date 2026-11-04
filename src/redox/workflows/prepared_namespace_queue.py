"""Bounded waiting for already admitted, immutable scientific namespaces.

Scientific admission stays with the caller. This shares process ownership,
resource-refusal accounting, source guards and the existing fleet dispatcher;
it never retries a started scientific calculation or certifies its outputs.
"""
import fcntl
import json
import os
from pathlib import Path
import signal
import threading
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify,load
from redox.workflows.charged_native_followup import dispatch_rows
from redox.workflows.sampling_native_queue import extend_peers
from redox.workflows.prelaunch_accounting import nonresource_prelaunch_failures,launcher_versions

TERMINAL_FAILURES={'scientific_failure_preserved','stopped_attempt_preserved',
    'guard_start_failed_preserved','uncertain_submission_preserved'}


def disposition(entries):
    conditions={row['condition'] for row in entries}
    if not entries:raise ValueError('An explicit nonempty scientific roster is required')
    if conditions=={'completed_review_pending'}:return 'finished_independent_review_pending'
    if conditions&TERMINAL_FAILURES and conditions<=TERMINAL_FAILURES|{'completed_review_pending'}:
        return 'attempt_preserved_for_explicit_review'
    return 'observe_or_dispatch'


def run(rows,output,peers_path,interpreter,hours,*,finished_status,validate,sources,resume=False,label='prepared_namespace'):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    from redox.validation.execution_review import freeze_source
    if not 0<hours<=72 or not finished_status:raise ValueError('Explicit completion and bounded positive dispatch window required')
    if len({r['key'] for r in rows})!=len(rows) or len({r['namespace'] for r in rows})!=len(rows):
        raise ValueError('Distinct prepared inputs required')
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    for row in rows:
        if sha256_file(row['namespace'])!=row['namespace_sha256']:raise ValueError('Prepared namespace changed')
        verify(Path(row['namespace']));validate(row)
    if not rows:raise ValueError('An explicit nonempty scientific roster is required')
    output=Path(output);output.mkdir(parents=True,exist_ok=resume)
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','prepared_queue_owners')
    remote=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','prepared_queue_inspection')
    paths=[*map(Path,sources),Path(__file__),*[Path(r['namespace']) for r in rows],
        ROOT/'src/redox/workflows/charged_native_followup.py',ROOT/'src/redox/workflows/screening_control.py',
        ROOT/'src/redox/workflows/molecular_native_namespace.py',
        ROOT/'scripts/fleet/run_screening_energy_panel.py',ROOT/'src/redox/workflows/prelaunch_accounting.py',
        ROOT/'src/redox/workflows/submission_refusal.py',ROOT/'scripts/pipeline/run_namespace_energy_queue.py',
        ROOT/'scripts/fleet/run_charged_polymer_geometry_panel.py']
    hashes={str(p):sha256_file(p) for p in paths}
    snapshots={str(freeze_source(output,Path(p))):h for p,h in hashes.items()}
    state=dict(status=label+'_controller_started',started=time.time(),deadline=time.time()+hours*3600,
        owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        CPU_affinity=sorted(os.sched_getaffinity(0)),source_sha256=snapshots,source_guard_sha256=hashes,
        namespace=rows[0]['namespace'],declared_namespaces=[r['namespace'] for r in rows],
        native_G_qualified=False,scientific_retries_performed=0,source_failure_preserved=True)
    stop=threading.Event()
    with (output/'controller.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        status_path=output/'status.json'
        if status_path.exists():
            old=json.loads(status_path.read_text())
            if not resume or owners.alive(old['owner']) is not False:
                raise ValueError('Explicit resume requires the exact previous controller to be closed')
            atomic_write_json(output/'history'/(str(time.time_ns())+'.json'),old)
        previous_handler=signal.signal(signal.SIGTERM,lambda *args:stop.set())
        try:
            atomic_write_json(status_path,state)
            while not stop.is_set() and time.time()<state['deadline']:
                if any(sha256_file(p)!=h for p,h in hashes.items()):
                    state['status']='source_changed_existing_attempt_preserved';break
                entries=[];eligible=[]
                for row in rows:
                    namespace=Path(row['namespace']);validate(row)
                    condition,detail=remote.inspect(row,namespace.parent/'dispatch_claim',owners.alive,finished_status=finished_status)
                    entry=dict(key=row['key'],namespace=row['namespace'],condition=condition,observation=detail)
                    if condition in ('pending','prelaunch_refusal_no_scientific_attempt'):
                        failures=nonresource_prelaunch_failures(namespace,owners.alive,launcher_versions(output))
                        entry['noncontention_prelaunch_failures']=failures
                        if failures<3:eligible.append(row)
                        else:entry['bounded_prelaunch_failure']=True
                    entries.append(entry)
                state.update(checked=time.time(),entries=entries,namespace=eligible[0]['namespace'] if eligible else None)
                outcome=disposition(entries)
                if outcome!='observe_or_dispatch':state['status']=label+'_'+outcome;break
                if any(r.get('bounded_prelaunch_failure') for r in entries):
                    state['status']='bounded_prelaunch_refusals_preserved';break
                atomic_write_json(status_path,state)
                if eligible and not stop.is_set():
                    peers=extend_peers(Path(peers_path),rows)
                    result=dispatch_rows(rows,output,interpreter,paths,storage_peers=peers,
                        eligible_keys={r['key'] for r in eligible})
                    state['last_dispatch']=str(output/'followup_dispatch_current.json')
                    state['status']=(label+'_submitted_start_unverified' if result['launches']
                        else 'waiting_for_suitable_free_GPU_or_host_capacity')
                else:state['status']=label+'_attempt_observed'
                atomic_write_json(status_path,state);stop.wait(30)
            else:
                state['status']='controller_stop_requested_workers_preserved' if stop.is_set() else 'dispatch_window_expired_attempt_preserved'
        except BaseException as error:
            state.update(status='prepared_queue_failed_workers_preserved',error=repr(error));raise
        finally:
            signal.signal(signal.SIGTERM,previous_handler)
            state['finished']=time.time();state['namespace']=None;atomic_write_json(status_path,state)
    return state
