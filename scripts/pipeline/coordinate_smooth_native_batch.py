#!/usr/bin/env python
"""Run a resumable, dependency-gated alternate-model native thermochemistry batch.

GPU/CPU placement reuses the validated reservation and memory-promise contract.
Each completed stage receives the separate alternate-model independent review;
no primary-energy normalization or native-G admission is performed.
"""
import argparse
import fcntl
import json
import inspect
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

SCRIPT = Path(__file__).resolve()
PROFILE = ROOT / 'environments/cutensor22/runtime.json'


def cpu_preflight(helper, batch, host, profile):
    body = inspect.getsource(helper.preflight)
    old = "helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'"
    if body.count(old) != 1: raise ValueError('CPU-only dependency template changed')
    body = body.replace(old, "helper = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'", 1)
    if body.count("'65536'") != 1:raise ValueError('Native RAM preflight template changed')
    body = body.replace("'65536'", "str(int((json.loads((batch/'batch.json').read_text())['maximum_RSS_GiB']+32)*1024))", 1)
    body = body.replace("if path.exists() and time.time() - path.stat().st_mtime < 600:",
        "if path.exists() and time.time() - path.stat().st_mtime < 600 and json.loads(path.read_text()).get('ok'):", 1)
    scope = {}; exec(compile(body, str(SCRIPT), 'exec'), dict(vars(helper)), scope)
    return scope['preflight'](batch, host, profile)


def native_placement(profile):
    base = api.load_module('smooth_native_reserved_placement', ROOT/'scripts/pipeline/coordinate_geometry_batch.py')
    body = inspect.getsource(base.launch_local)
    old = "manifest.get('stage') == 'electronic'"
    anchor = "        with (attempt / 'worker.log').open('x') as log:"
    if body.count(old) != 1 or body.count(anchor) != 1:
        raise ValueError('Deferred native placement template changed')
    body = body.replace(old, "manifest.get('stage') == 'native'", 1)
    capacity = "        capacity = host_capacity(2, minimum_ram_gib=bound + 32)"
    if body.count(capacity) != 1:raise ValueError('Native host-capacity template changed')
    body = body.replace(capacity, """        try:
            capacity = host_capacity(2, minimum_ram_gib=bound + 32)
        except RuntimeError as error:
            return dict(status='resource_capacity_unavailable', host=host, reason=str(error))""", 1)
    inserted = """        if plan.get('runtime_profile') != str(_PROFILE):
            raise ValueError('Unexpected native runtime profile')
        argv.insert(argv.index('env')+1, 'LD_LIBRARY_PATH='+_profile['library_path'])
        index=argv.index(manifest['worker'])
        argv[index:index]=[str(ROOT/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(_PROFILE),
            '--out', str(campaign/'tensor_runtime')]
"""
    scope = {}; namespace = dict(vars(base), _PROFILE=PROFILE, _profile=profile)
    exec(compile(body.replace(anchor, inserted+anchor, 1), str(SCRIPT), 'exec'), namespace, scope)
    base.launch_local = scope['launch_local']
    return base


def await_declaration(request_path, batch, hours):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard(): host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1))
    request = json.loads(request_path.read_text()); api.verify_sources(request['source_sha256'])
    if request['output_batch'] != str(batch) or request['primary_candidate_admission']:
        raise ValueError('Native declaration request has another destination or scope')
    review = Path(request['native_review']); reference = Path(request['native_reference'])
    record = dict(status='waiting_for_independently_checked_native_resource_pilot',
        owner=dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        started=time.time(), deadline=time.time()+hours*3600, request_sha256=api.digest(request_path),
        source_sha256=request['source_sha256'], primary_candidate_admission=False, scientific_starts=0)
    status = request_path.with_name(request_path.stem+'_status.json')
    try:
        while time.time() < record['deadline']:
            api.verify_sources(request['source_sha256'])
            if api.digest(request_path) != record['request_sha256']:raise ValueError('Declaration request changed')
            if review.exists():
                preparer = api.load_module('smooth_native_measured_declaration', ROOT/'scripts/pipeline/prepare_smooth_native_batch.py')
                manifest = preparer.declare([Path(p) for p in request['electronic_batches']], reference, review, batch)
                record.update(status='checked_pilot_declared_dependency_gated_native_cohort', declared=len(manifest['entries']))
                return True
            path = reference/'experimental_native/status.json'
            if path.exists():
                row = json.loads(path.read_text())
                if 'finished' in row and row['status'] != 'experimental_native_finished_independent_review_pending':
                    record.update(status='native_resource_pilot_failed_no_expansion', source_status=row['status']);return False
            record['updated']=time.time();atomic_write_json(status,record);time.sleep(15)
        record['status']='native_declaration_window_ended_no_expansion';return False
    except BaseException as error:
        record.update(status='native_declaration_failed_workers_untouched',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(status,record)


def coordinate(batch, hosts, hours, resume, placement, profile):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    manifest = json.loads((batch / 'batch.json').read_text())
    if (manifest['stage'] != 'native' or manifest['primary_candidate_admission']
            or manifest['worker'] != str(ROOT/'scripts/pipeline/smooth_polymer_native.py')):
        raise ValueError('Only the declared alternate-model native thermochemistry batch is supported')
    api.verify_sources(manifest['source_sha256'])
    preparer = api.load_module('smooth_native_dependency_preparer', Path(manifest['preparer']))
    reviewer = api.load_module('smooth_native_independent_reviewer', ROOT/'scripts/validation/review_smooth_polymer_native.py')
    sources = {**manifest['source_sha256'], str(SCRIPT):api.digest(SCRIPT),
               str(batch/'batch.json'):api.digest(batch/'batch.json')}
    with allocation_guard(): host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1))
    with (batch / 'dispatch.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = batch / 'dispatch.json'
        if path.exists():
            prior = json.loads(path.read_text())
            if not resume or placement.owner_live(prior['owner']):
                raise ValueError('Prior controller is live or explicit resume is missing')
            history = batch/'dispatch_history'; history.mkdir(exist_ok=True)
            (history/(str(time.time_ns())+'.json')).write_bytes(path.read_bytes())
        owner = dict(host=socket.gethostname().split('.')[0], pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
        record = dict(status='alternate_native_batch_in_progress', owner=owner,
            started=time.time(), deadline=time.time()+hours*3600, source_sha256=sources,
            CPU_affinity=sorted(os.sched_getaffinity(0)), native_G_in_primary_cache=False,
            primary_candidate_admission=False)
        failures, processed, observations = {}, set(), {}
        try:
            while time.time() < record['deadline']:
                api.verify_sources(sources)
                pending, rows = [], []
                for entry in manifest['entries']:
                    key, case = entry['key'], Path(entry['campaign'])
                    if entry['existing_checked_reference']:
                        resource = json.loads(Path(manifest['resource_review']).read_text())
                        api.verify_sources(resource['source_sha256'])
                        rows.append(dict(key=key, campaign=str(case), status='checked_identical_reference',
                            G_eV=resource['G_eV'], review=manifest['resource_review'],
                            native_G_in_primary_cache=False)); continue
                    try:
                        ready = preparer.prepare_ready(entry, manifest)
                    except (OSError, subprocess.SubprocessError) as error:
                        rows.append(dict(key=key, campaign=str(case),
                            status='upstream_ownership_or_IO_unverified_preserved', error=repr(error))); continue
                    if ready != 'prepared':
                        rows.append(dict(key=key, campaign=str(case), status=ready)); continue
                    state_path = case/'experimental_native/status.json'
                    state = json.loads(state_path.read_text()) if state_path.exists() else {}
                    launch_path = case/'launch.json'
                    handle = json.loads(launch_path.read_text()) if launch_path.exists() else observations.get(key)
                    item = dict(key=key, campaign=str(case))
                    if 'finished' in state:
                        passed = state['status']=='experimental_native_finished_independent_review_pending'
                        if passed:
                            guard = case/'residency_guard/owner.json'
                            try: live = placement.owner_live(state['owner'])
                            except (OSError, RuntimeError, subprocess.SubprocessError): live = None
                            if live is not False or not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):
                                item['status']='completed_waiting_for_closed_ownership'; rows.append(item); continue
                            report = case/'independent_native_review.json'
                            if not report.exists():
                                reviewer.review(case, report)
                            checked = json.loads(report.read_text()); api.verify_sources(checked['source_sha256']); preparer.validate_review(checked, manifest['experimental_method_id'])
                            if (checked['status']!='alternate_native_state_independently_reconstructed'
                                    or checked['experimental_method_id']!=manifest['experimental_method_id']
                                    or checked['primary_candidate_admission']):
                                raise ValueError('Completed native review has a different model or scope')
                            item.update(status='checked_native', G_eV=checked['G_eV'],
                                review=str(report), review_sha256=api.digest(report), native_G_in_primary_cache=False)
                        else:
                            item.update(status='failed_preserved', error=state.get('error'))
                        if key not in processed:
                            host=state['owner']['host']; failures[host]=0 if passed else failures.get(host,0)+1
                            processed.add(key)
                    elif handle and handle.get('owner'):
                        try: live=placement.owner_live(handle['owner']) if handle['owner']['pid']>0 else None
                        except (OSError,RuntimeError,subprocess.SubprocessError): live=None
                        if live is not False: item['status']='running' if live else 'submission_ownership_unverified_preserved'
                        else:
                            item['status']='stopped_attempt_preserved_review_required'
                            if key not in processed:
                                host=handle['owner']['host'];failures[host]=failures.get(host,0)+1;processed.add(key)
                    else:
                        intent=case/'dispatch_intent.json'
                        unresolved=intent.exists() and json.loads(intent.read_text()).get('outcome')!='not_started'
                        item['status']='submission_ownership_unverified_preserved' if unresolved else 'pending'
                        if not unresolved: pending.append(entry)
                    rows.append(item)
                record.update(updated=time.time(), entries=rows,
                    checked_native=sum(r['status'] in ('checked_native','checked_identical_reference') for r in rows),
                    running=sum(r['status']=='running' for r in rows),
                    preserved_failures=sum(r['status'] in ('failed_preserved','stopped_attempt_preserved_review_required') for r in rows),
                    consecutive_host_failures=failures)
                atomic_write_json(path,record)
                waiting={'running','upstream_pending','upstream_ownership_or_IO_unverified_preserved',
                         'submission_ownership_unverified_preserved','completed_waiting_for_closed_ownership'}
                if not pending and not any(r['status'] in waiting for r in rows):
                    record['status']='alternate_native_batch_terminal_primary_cache_separate';return
                slots=free_slots(hosts)
                for host in sorted({h for h,_ in slots}):
                    if not pending or failures.get(host,0)>=3: continue
                    helper=api.load_module('smooth_native_node_preflight',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
                    if not cpu_preflight(helper,batch,host,profile):continue
                    for _,gpu in [s for s in slots if s[0]==host]:
                        if not pending:break
                        entry=pending[0];case=Path(entry['campaign']);intent=case/'dispatch_intent.json'
                        atomic_write_json(intent,dict(host=host,GPU=gpu,owner=owner,checked=time.time()))
                        command=['env','OMP_NUM_THREADS=1','OMP_THREAD_LIMIT=1','OPENBLAS_NUM_THREADS=1',
                            'MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','PYTHONPATH='+str(ROOT/'src'),profile['python'],
                            str(SCRIPT),'launch-local','--batch',str(batch),'--campaign',str(case),'--gpu',str(gpu)]
                        try:
                            response=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(command)],
                                capture_output=True,text=True,timeout=45,check=True)
                            observation=json.loads(response.stdout)
                        except (OSError,ValueError,subprocess.SubprocessError) as error:
                            observations[entry['key']]=dict(owner=dict(host=host,pid=-1,process_start='unknown'),error=repr(error))
                            pending.pop(0);continue
                        with (batch/'dispatch_events.jsonl').open('a') as stream:
                            stream.write(json.dumps(dict(key=entry['key'],observation=observation,checked=time.time()))+'\n')
                        if observation['status']=='gate_started_scientific_start_pending':
                            observations[entry['key']]=observation;pending.pop(0)
                        else:
                            atomic_write_json(intent,dict(host=host,GPU=gpu,owner=owner,checked=time.time(),outcome='not_started'))
                            if observation['status']=='resource_capacity_unavailable':break
                time.sleep(15)
            record['status']='dispatch_window_ended_workers_and_outputs_preserved'
        except BaseException as error:
            record.update(status='alternate_native_controller_failed_workers_untouched',error=repr(error));raise
        finally:
            record['finished']=time.time();atomic_write_json(path,record)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('coordinate','launch-local'))
    parser.add_argument('--batch',type=Path,required=True)
    parser.add_argument('--request',type=Path)
    parser.add_argument('--campaign',type=Path);parser.add_argument('--gpu',type=int)
    parser.add_argument('--hosts',default='lambda1,lambda2,lambda11,lambda13')
    parser.add_argument('--max-hours',type=float,default=24);parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();profile=json.loads(PROFILE.read_text())
    if Path(sys.executable)!=Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH'):
        parser.error('Use the declared interpreter with the host loader environment')
    if not 0<args.max_hours<=48 or any(not h or h.startswith('-') for h in args.hosts.split(',')):
        parser.error('Use valid hosts and a finite window up to 48 hours')
    helper=api.load_module('smooth_native_measured_resource_placement',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    placement=native_placement(profile)
    if args.action=='coordinate' and not (args.batch/'batch.json').exists():
        if args.request is None:parser.error('A source-bound request is required before declaration')
        with args.request.with_suffix('.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            if not await_declaration(args.request.resolve(),args.batch.resolve(),args.max_hours):return
    if args.action=='launch-local':
        if args.campaign is None or args.gpu is None:parser.error('launch-local requires campaign and GPU')
        print(json.dumps(placement.launch_local(args.batch.resolve(),args.campaign.resolve(),args.gpu)))
    else:coordinate(args.batch.resolve(),args.hosts.split(','),args.max_hours,args.resume,placement,profile)


if __name__=='__main__':main()
