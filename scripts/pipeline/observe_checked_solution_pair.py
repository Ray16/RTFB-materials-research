#!/usr/bin/env python3
"""Wait for closed pair ownership, then independently audit saved solution records."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.solution_pair_review import review,ReviewPending
from redox.workflows.solution_pair_namespace import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--spec',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--hours',type=float,default=24);args=parser.parse_args()
    if not 0<args.hours<=48:parser.error('Bounded observation window up to 48 hours required')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False);spec=verify(args.spec)
    definition=importlib.util.spec_from_file_location('solution_pair_owner_checks',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    sources={str(p):sha256_file(p) for p in [Path(__file__).resolve(),args.spec,
        ROOT/'src/redox/validation/solution_pair_review.py',ROOT/'src/redox/workflows/solution_pair_namespace.py',
        ROOT/'src/redox/validation/protonation.py',ROOT/'src/redox/build/geometry.py',
        ROOT/'src/redox/properties/solution_reorg.py']}
    record=dict(status='waiting_for_solution_pair_and_closed_ownership',source_sha256=sources,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        started=time.time(),deadline=time.time()+args.hours*3600,solution_relaxation_values_qualified=False)
    while time.time()<record['deadline']:
        if any(sha256_file(path)!=digest for path,digest in sources.items()):raise ValueError('Observer/review source changed')
        stage=Path(spec['case'])/'status.json'
        if stage.exists():
            raw=json.loads(stage.read_text());record['scientific_status']=raw['status']
            if 'finished' in raw:
                try:live=owners.alive(raw['owner'])
                except (OSError,ValueError,subprocess.SubprocessError):live=None
                guard=Path(spec['case'])/'residency_guard/owner.json'
                if live is False and guard.exists() and json.loads(guard.read_text()).get('guard_stopped'):
                    if raw['status']=='solution_pair_finished_independent_review_pending':
                        try:result=review(args.spec,args.out/'independent_review.json')
                        except (ReviewPending,OSError) as error:
                            record.update(status='solution_pair_review_observation_unverified',error=repr(error),updated=time.time())
                            atomic_write_json(args.out/'status.json',record);time.sleep(20);continue
                        except (ValueError,subprocess.SubprocessError) as error:
                            record.update(status='solution_pair_independent_review_failed_preserved',error=repr(error))
                        else:
                            record.update(status='solution_pair_independently_checked',review=str(args.out/'independent_review.json'),
                                review_sha256=sha256_file(args.out/'independent_review.json'),
                                solution_relaxation_values_qualified=result['summary']['endpoint_model_verified'])
                    else:record.update(status='solution_pair_failed_preserved',error=raw.get('error'))
                    record['finished']=time.time();atomic_write_json(args.out/'status.json',record);return
        record['updated']=time.time();atomic_write_json(args.out/'status.json',record);time.sleep(20)
    record.update(status='solution_pair_observation_window_ended_worker_untouched',finished=time.time())
    atomic_write_json(args.out/'status.json',record)


if __name__=='__main__':main()
