#!/usr/bin/env python
"""CPU follow-up for explicitly declared local electronic resource pilots."""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--max-hours',type=float,default=24);a=p.parse_args()
    if not 0<a.max_hours<=48:p.error('Use a finite follow-up window up to 48 hours')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    campaign=a.campaign.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=campaign/'campaign.json';manifest=json.loads(source.read_text());reviewer_path=ROOT/'scripts/validation/review_local_motif_energy.py'
    reviewer=api.load_module('CPU_local_energy_independent_review',reviewer_path)
    owners=api.load_module('CPU_local_energy_stopped_ownership',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),source,reviewer_path)}
    row=dict(status='local_electronic_resource_pilots_CPU_followup_running',started=time.time(),deadline=time.time()+a.max_hours*3600,
        owner=dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        source_sha256=sources,CPU_affinity=sorted(os.sched_getaffinity(0)),quantum_jobs_submitted=0,checked={})
    try:
        while time.time()<row['deadline']:
            api.verify_sources(sources);pending=0
            for entry in manifest['entries']:
                if entry['key'] in row['checked']:continue
                case=Path(entry['campaign']);status=case/'electronic/status.json'
                if api.digest(case/'plan.json')!=entry['plan_sha256']:raise ValueError('Declared local pilot changed')
                if not status.exists():pending+=1;continue
                state=json.loads(status.read_text())
                if 'finished' not in state:pending+=1;continue
                try:live=owners.owner_alive(state['owner'])
                except (OSError,RuntimeError,subprocess.SubprocessError):pending+=1;continue
                if live:pending+=1;continue
                if state['status']!='frozen_local_energy_finished_independent_review_pending':
                    row['checked'][entry['key']]=dict(status='source_failure_preserved',source_status=state['status'],error=state.get('error'));continue
                guard=case/'residency_guard/owner.json'
                if not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):pending+=1;continue
                report=out/'checked'/(entry['key']+'.json')
                try:reviewer.review(case,report)
                except Exception as error:row['checked'][entry['key']]=dict(status='independent_review_failed_preserved',error=repr(error))
                else:row['checked'][entry['key']]=dict(status='independent_review_passed',report=str(report),report_sha256=api.digest(report))
            row.update(updated=time.time(),pending=pending,independent_reviews_passed=sum(r['status']=='independent_review_passed' for r in row['checked'].values()))
            atomic_write_json(out/'status.json',row)
            if not pending:row['status']='all_local_pilots_reviewed_or_failures_preserved';return
            time.sleep(15)
        row['status']='local_pilot_followup_window_ended_pending_preserved'
    except BaseException as error:row.update(status='CPU_local_pilot_observer_failed_workers_untouched',error=repr(error));raise
    finally:row['finished']=time.time();atomic_write_json(out/'status.json',row)


if __name__=='__main__':main()
