#!/usr/bin/env python
"""CPU review of a completed alternate-native pilot, without quantum launch."""
import argparse
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--max-hours',type=float,default=24);a=p.parse_args()
    if not 0<a.max_hours<=48:p.error('Use a bounded observation window')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    case,out=a.campaign.resolve(),a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    checker=ROOT/'scripts/validation/review_smooth_polymer_native.py'
    sources={str(q):api.digest(q) for q in (Path(__file__).resolve(),checker,case/'plan.json')}
    owner=api.load_module('alternate_native_observer_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    reviewer=api.load_module('alternate_native_observer_review',checker)
    record=dict(status='waiting_for_completed_alternate_native_pilot',started=time.time(),deadline=time.time()+a.max_hours*3600,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        CPU_affinity=sorted(os.sched_getaffinity(0)),source_sha256=sources,quantum_jobs_submitted=0,
        primary_candidate_admission=False)
    try:
        while time.time()<record['deadline']:
            api.verify_sources(sources);path=case/'experimental_native/status.json'
            state=json.loads(path.read_text()) if path.exists() else {}
            if 'finished' in state:
                if state['status']!='experimental_native_finished_independent_review_pending':
                    record.update(status='native_source_failure_preserved',source_status=state['status'],error=state.get('error'));return
                guard=case/'residency_guard/owner.json'
                if not owner.owner_alive(state['owner']) and guard.exists() and json.loads(guard.read_text()).get('guard_stopped'):
                    reviewer.review(case,out/'independent_review.json')
                    record.update(status='alternate_native_pilot_independently_reviewed',
                        review=str(out/'independent_review.json'),review_sha256=api.digest(out/'independent_review.json'));return
            record.update(updated=time.time(),observed_native_status=state.get('status','not_started'))
            atomic_write_json(out/'status.json',record);time.sleep(15)
        record['status']='observation_window_ended_outputs_preserved'
    except BaseException as error:
        record.update(status='CPU_observer_failed_scientific_work_untouched',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(out/'status.json',record)


if __name__=='__main__':main()
