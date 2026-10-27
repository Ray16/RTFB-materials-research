#!/usr/bin/env python
"""Independently review the new AB checkpoint and repeat orbital stationarity."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

REVIEW=ROOT/'scripts/validation/review_counterpoise_cpu_remaining.py'
CHECK=ROOT/'scripts/validation/check_counterpoise_remaining_stationarity.py'


def alive(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Observer requires the producer host')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return False
    return fields[19]==owner['process_start'] and fields[0] not in ('Z','X')


def watch(raw,out,hours,threads):
    if out.exists():raise FileExistsError('Preserve prior reviews')
    launcher=json.loads((raw/'launcher.json').read_text())
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),REVIEW,CHECK,raw/'launch.json',raw/'launcher.json')}
    with allocation_guard():
        host_capacity(1);permitted=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,cpu_allocation(1))
    out.mkdir(parents=True);record=dict(status='waiting_for_refined_AB_checkpoint',started=time.time(),deadline=time.time()+hours*3600,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        producer_owner=launcher['owner'],source_sha256=sources,full_counterpoise_qualified=False)
    child=None
    try:
        while time.time()<record['deadline']:
            api.verify_sources(sources);state=json.loads((raw/'status.json').read_text());live=alive(launcher['owner'])
            if child is None:
                if live:
                    record['status']='waiting_for_refined_AB_checkpoint'
                elif state.get('status')!='CPU_remaining_term_converged_independent_review_pending':
                    record.update(status='refinement_ended_without_reviewable_result',producer_status=state);return
                else:
                    checkpoint_review=out/'review.json'
                    if not checkpoint_review.exists():
                        api.load_module('independent_AB_refinement_checkpoint',REVIEW).review(raw,checkpoint_review)
                    try:
                        with allocation_guard():
                            previous=sorted(os.sched_getaffinity(0))
                            try:
                                os.sched_setaffinity(0,permitted);host_capacity(threads,minimum_ram_gib=112);cpus=cpu_allocation(threads)
                            finally:os.sched_setaffinity(0,previous)
                            scratch=out/'scratch';scratch.mkdir(exist_ok=True)
                            argv=[sys.executable,str(CHECK),'--raw',str(raw),'--review',str(checkpoint_review),'--out',str(out/'stationarity'),'--threads',str(threads)]
                            env=dict(os.environ,PYSCF_TMPDIR=str(scratch))
                            for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]=str(threads)
                            with (out/'stationarity.log').open('xb') as log:
                                # The checker allocates its disjoint mask under the same guard after startup.
                                child=subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,env=env,preexec_fn=lambda:os.sched_setaffinity(0,set(cpus)))
                            handle=dict(host=os.uname().nodename,pid=child.pid,process_start=Path('/proc/'+str(child.pid)+'/stat').read_text().rsplit(')',1)[1].split()[19])
                            record.update(status='refined_AB_independent_stationarity_running',stationarity_owner=handle,stationarity_CPU_affinity=cpus)
                            atomic_write_json(out/'stationarity_launch.json',dict(owner=handle,argv=argv,CPU_affinity=cpus))
                    except RuntimeError as error:
                        if child is not None or (out/'stationarity.log').exists():raise
                        record.update(status='waiting_for_independent_stationarity_CPU_capacity',capacity_error=str(error))
            elif child.poll() is not None:
                path=out/'stationarity/summary.json';result=json.loads(path.read_text()) if path.exists() else {}
                if child.returncode==0 and result.get('passed') is True:
                    record.update(status='refined_AB_independent_stationarity_passed',result=str(path),result_sha256=api.digest(path));return
                record.update(status='refined_AB_independent_stationarity_failed_preserved',returncode=child.returncode,result=result);return
            record['checked']=time.time();atomic_write_json(out/'status.json',record);time.sleep(15)
        record['status']='finite_observer_ended_outputs_and_workers_preserved'
    except BaseException as error:
        record.update(status='refinement_observer_failed_preserved',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(out/'status.json',record)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=24);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if not 0<a.max_hours<=48 or not 1<=a.threads<=80:raise ValueError('Finite runtime and valid CPU count required')
    watch(a.raw.resolve(),a.out.resolve(),a.max_hours,a.threads)
