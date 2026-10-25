#!/usr/bin/env python
"""Run a resumable, dependency-gated alternate-model electronic batch.

GPU/CPU placement reuses the validated reservation and memory-promise contract.
Each completed stage receives the separate alternate-model independent review;
no primary-energy normalization or native-G admission is performed.
"""
import argparse
import fcntl
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

SCRIPT = Path(__file__).resolve()
PROFILE = ROOT / 'environments/cutensor22/runtime.json'


def coordinate(batch, hosts, hours, resume, placement, profile):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    manifest = json.loads((batch / 'batch.json').read_text())
    if (manifest['stage'] != 'electronic' or manifest['primary_candidate_admission']
            or manifest['native_G_computed'] or manifest['worker'] != str(ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py')):
        raise ValueError('Only the declared alternate-model electronic batch is supported')
    api.verify_sources(manifest['source_sha256'])
    preparer = api.load_module('smooth_electronic_dependency_preparer', Path(manifest['preparer']))
    reviewer = api.load_module('smooth_electronic_independent_reviewer', ROOT/'scripts/validation/review_smooth_polymer_energy.py')
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
        record = dict(status='alternate_electronic_batch_in_progress', owner=owner,
            started=time.time(), deadline=time.time()+hours*3600, source_sha256=sources,
            CPU_affinity=sorted(os.sched_getaffinity(0)), native_G_computed=False,
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
                            E_smd_eV=resource['E_smd_eV'], review=manifest['resource_review'],
                            native_G_computed=False)); continue
                    try:
                        ready = preparer.prepare_ready(entry, manifest)
                    except (OSError, subprocess.SubprocessError) as error:
                        rows.append(dict(key=key, campaign=str(case),
                            status='upstream_ownership_or_IO_unverified_preserved', error=repr(error))); continue
                    if ready != 'prepared':
                        rows.append(dict(key=key, campaign=str(case), status=ready)); continue
                    state_path = case/'experimental_energy/status.json'
                    state = json.loads(state_path.read_text()) if state_path.exists() else {}
                    launch_path = case/'launch.json'
                    handle = json.loads(launch_path.read_text()) if launch_path.exists() else observations.get(key)
                    item = dict(key=key, campaign=str(case))
                    if 'finished' in state:
                        passed = state['status']=='experimental_hybrid_energy_finished_independent_review_pending'
                        if passed:
                            guard = case/'residency_guard/owner.json'
                            try: live = placement.owner_live(state['owner'])
                            except (OSError, RuntimeError, subprocess.SubprocessError): live = None
                            if live is not False or not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):
                                item['status']='completed_waiting_for_closed_ownership'; rows.append(item); continue
                            report = case/'independent_energy_review.json'
                            if not report.exists():
                                reviewer.review(case, report)
                            checked = json.loads(report.read_text()); api.verify_sources(checked['source_sha256'])
                            if (checked['status']!='alternate_electronic_geometry_and_resource_checked_native_G_pending'
                                    or checked['portable_method_id']!=manifest['experimental_method_id']
                                    or checked['native_Hessian_or_G_qualified']):
                                raise ValueError('Completed energy review has a different model or scope')
                            item.update(status='checked_electronic', E_smd_eV=checked['E_smd_eV'],
                                review=str(report), review_sha256=api.digest(report), native_G_computed=False)
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
                    checked_electronic=sum(r['status'] in ('checked_electronic','checked_identical_reference') for r in rows),
                    running=sum(r['status']=='running' for r in rows),
                    preserved_failures=sum(r['status'] in ('failed_preserved','stopped_attempt_preserved_review_required') for r in rows),
                    consecutive_host_failures=failures)
                atomic_write_json(path,record)
                waiting={'running','upstream_pending','upstream_ownership_or_IO_unverified_preserved',
                         'submission_ownership_unverified_preserved','completed_waiting_for_closed_ownership'}
                if not pending and not any(r['status'] in waiting for r in rows):
                    record['status']='alternate_electronic_batch_terminal_native_G_separate';return
                slots=free_slots(hosts)
                for host in sorted({h for h,_ in slots}):
                    if not pending or failures.get(host,0)>=3: continue
                    helper=api.load_module('smooth_electronic_node_preflight',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
                    if not helper.preflight(batch,host,profile):continue
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
            record.update(status='alternate_electronic_controller_failed_workers_untouched',error=repr(error));raise
        finally:
            record['finished']=time.time();atomic_write_json(path,record)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('coordinate','launch-local'))
    parser.add_argument('--batch',type=Path,required=True)
    parser.add_argument('--campaign',type=Path);parser.add_argument('--gpu',type=int)
    parser.add_argument('--hosts',default='lambda1,lambda2,lambda11,lambda13')
    parser.add_argument('--max-hours',type=float,default=24);parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();profile=json.loads(PROFILE.read_text())
    if Path(sys.executable)!=Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH'):
        parser.error('Use the declared interpreter with the host loader environment')
    if not 0<args.max_hours<=48 or any(not h or h.startswith('-') for h in args.hosts.split(',')):
        parser.error('Use valid hosts and a finite window up to 48 hours')
    helper=api.load_module('smooth_electronic_measured_resource_placement',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    placement=helper.placement_module(profile)
    if args.action=='launch-local':
        if args.campaign is None or args.gpu is None:parser.error('launch-local requires campaign and GPU')
        print(json.dumps(placement.launch_local(args.batch.resolve(),args.campaign.resolve(),args.gpu)))
    else:coordinate(args.batch.resolve(),args.hosts.split(','),args.max_hours,args.resume,placement,profile)


if __name__=='__main__':main()
