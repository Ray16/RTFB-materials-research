#!/usr/bin/env python
"""Run all eligible local electronic inputs within measured resource classes."""
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


def coordinate(campaign,hours,hosts):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,free_slots
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    manifest=json.loads((campaign/'campaign.json').read_text());api.verify_sources(manifest['source_sha256'])
    profile=json.loads((ROOT/'environments/cutensor22/runtime.json').read_text())
    if str(Path(sys.executable))!=profile['python'] or os.environ.get('LD_LIBRARY_PATH'):
        raise ValueError('Use the profile interpreter with the host loader environment')
    sources={**manifest['source_sha256'],str(campaign/'campaign.json'):api.digest(campaign/'campaign.json'),
        str(Path(__file__).resolve()):api.digest(__file__),str(ROOT/'scripts/validation/review_local_motif_energy.py'):api.digest(ROOT/'scripts/validation/review_local_motif_energy.py')}
    reviewer=api.load_module('local_batch_independent_electronic_review',ROOT/'scripts/validation/review_local_motif_energy.py')
    owners=api.load_module('local_batch_exact_ownership',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    preflight=api.load_module('local_batch_CPU_only_preflight',ROOT/'scripts/pipeline/dispatch_smooth_energy_capacity.py')
    helper=api.load_module('local_batch_preflight_template',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    probe=preflight.capacity_preflight(helper)
    # Preflight reads the same shared resource fields as the larger queues.
    batch_path=campaign/'batch.json'
    if not batch_path.exists():atomic_write_json(batch_path,manifest)
    sentinel=manifest.get('far_pair_resource_sentinel')
    with (campaign/'dispatch.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        path=campaign/'dispatch.json'
        if path.exists():raise ValueError('Preserve and inspect the existing local batch owner')
        record=dict(status='local_frozen_energy_batch_in_progress',started=time.time(),deadline=time.time()+hours*3600,
            owner=dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
            source_sha256=sources,primary_candidate_admission=False,thermochemistry_computed=False,
            CPU_affinity=sorted(os.sched_getaffinity(0)))
        failures={};processed=set()
        try:
            while time.time()<record['deadline']:
                api.verify_sources(sources);rows=[];pending=[]
                for entry in manifest['entries']:
                    case=Path(entry['campaign']);key=entry['key']
                    if api.digest(case/'plan.json')!=entry['plan_sha256']:raise ValueError('Declared local plan changed')
                    row=dict(key=key,campaign=str(case));state_path=case/'electronic/status.json'
                    state=json.loads(state_path.read_text()) if state_path.exists() else {}
                    if 'finished' in state:
                        try:live=owners.owner_alive(state['owner'])
                        except (OSError,RuntimeError,subprocess.SubprocessError):row['status']='ownership_unverified_preserved';rows.append(row);continue
                        guard=case/'residency_guard/owner.json'
                        if live or not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):
                            row['status']='finished_waiting_for_closed_ownership';rows.append(row);continue
                        passed=state['status']=='frozen_local_energy_finished_independent_review_pending'
                        if passed:
                            report=case/'independent_energy_review.json'
                            if not report.exists():reviewer.review(case,report)
                            checked=json.loads(report.read_text());api.verify_sources(checked['source_sha256'])
                            row.update(status='checked_electronic',report=str(report),report_sha256=api.digest(report))
                        else:row.update(status='failed_preserved',error=state.get('error'))
                        if key not in processed:
                            host=state['owner']['host'];failures[host]=0 if passed else failures.get(host,0)+1;processed.add(key)
                    elif (case/'launch.json').exists():
                        launch=json.loads((case/'launch.json').read_text())
                        try:live=owners.owner_alive(launch['owner'])
                        except (OSError,RuntimeError,subprocess.SubprocessError):live=None
                        row['status']='running' if live else 'submission_ownership_unverified_preserved'
                    elif (case/'submission_intent.json').exists():row['status']='submission_ownership_unverified_preserved'
                    else:row['status']='pending';pending.append(entry)
                    rows.append(row)
                sentinel_row=next((r for r in rows if r['key']==sentinel),None)
                if sentinel and sentinel_row['status']=='failed_preserved':
                    record.update(status='far_pair_resource_sentinel_failed_other_pair_inputs_preserved',entries=rows);return
                if sentinel and sentinel_row['status']!='checked_electronic':pending=[e for e in pending if e['key']==sentinel]
                record.update(updated=time.time(),entries=rows,checked_electronic=sum(r['status']=='checked_electronic' for r in rows),
                    running=sum(r['status']=='running' for r in rows),preserved_failures=sum(r['status']=='failed_preserved' for r in rows),consecutive_host_failures=failures)
                atomic_write_json(path,record)
                if pending and not record['running'] and all(failures.get(h,0)>=3 for h in hosts):
                    record['status']='local_batch_stopped_after_three_actual_node_failures';return
                if all(r['status'] in ('checked_electronic','failed_preserved','submission_ownership_unverified_preserved') for r in rows):
                    record['status']='local_batch_terminal_with_source_failures_or_unresolved_submissions_preserved';return
                slots=free_slots(hosts)
                for host in sorted({h for h,_ in slots}):
                    if not pending or failures.get(host,0)>=3 or not probe(campaign,host,profile):continue
                    for _,gpu in (s for s in slots if s[0]==host):
                        if not pending:break
                        entry=pending[0];case=Path(entry['campaign']);intent=case/'submission_intent.json'
                        atomic_write_json(intent,dict(host=host,GPU=gpu,owner=record['owner'],checked=time.time()))
                        args=['env','-u','LD_LIBRARY_PATH','OMP_NUM_THREADS=1','OMP_THREAD_LIMIT=1','OPENBLAS_NUM_THREADS=1',
                            'MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','PYTHONPATH='+str(ROOT/'src'),profile['python'],manifest['dispatcher'],
                            'launch','--campaign',str(case),'--gpu',str(gpu),'--preflight',str(campaign/(host+'_preflight.json'))]
                        try:
                            result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(args)],capture_output=True,text=True,timeout=45,check=True)
                            launch=json.loads(result.stdout)
                        except (OSError,ValueError,subprocess.SubprocessError) as error:
                            atomic_write_json(intent,dict(host=host,GPU=gpu,owner=record['owner'],checked=time.time(),error=repr(error)));pending.pop(0);continue
                        if launch['status']=='gate_started_scientific_start_pending':pending.pop(0)
                        else:
                            # An explicit pre-spawn refusal leaves no scientific owner.
                            history=case/'placement_history';history.mkdir(exist_ok=True)
                            intent.rename(history/(str(time.time_ns())+'.json'))
                            if launch['status']=='resource_capacity_unavailable':break
                time.sleep(15)
            record['status']='local_dispatch_window_ended_existing_workers_preserved'
        except BaseException as error:record.update(status='local_batch_CPU_controller_failed_workers_untouched',error=repr(error));raise
        finally:record['finished']=time.time();atomic_write_json(path,record)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--max-hours',type=float,default=24);p.add_argument('--hosts',default='lambda5');a=p.parse_args()
    if not 0<a.max_hours<=48:p.error('Use a finite dispatch window up to 48 hours')
    coordinate(a.campaign.resolve(),a.max_hours,a.hosts.split(','))


if __name__=='__main__':main()
