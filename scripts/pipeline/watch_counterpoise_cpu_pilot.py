#!/usr/bin/env python3
"""Watch one CPU pilot and review it only after its owner has finished."""
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

REVIEWER=ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py'


def owner_state(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify an owner on another host')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return 'stopped'
    if fields[19]!=owner['process_start'] or fields[0] in ('Z','X'):return 'stopped'
    return 'live'


def classify(status,owner):
    state=status.get('status')
    if state=='CPU_pilot_converged_independent_review_pending':
        return 'review_ready' if owner=='stopped' else 'waiting_for_owner_exit'
    if state=='CPU_pilot_interrupted_or_failed_preserved':return 'pilot_failed_preserved'
    if state in ('CPU_pilot_initializing_unqualified','CPU_pilot_SCF_running_unqualified'):
        return 'waiting_for_pilot' if owner=='live' else 'pilot_owner_stopped_without_terminal_status'
    raise ValueError('Unknown pilot state; do not promote')


def watch(raw,out,max_hours):
    launch=json.loads((out/'launch.json').read_text());sources=launch['source_sha256'];deadline=time.time()+max_hours*3600
    watch_owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status_path=out/'status.json'
    while time.time()<deadline:
        api.verify_sources(sources)
        status=json.loads((raw/'status.json').read_text());owner=status['owner'];state=classify(status,owner_state(owner))
        atomic_write_json(status_path,dict(status=state,checked=time.time(),owner=watch_owner,pilot_owner=owner,source_sha256=sources,
            review_complete=False,full_counterpoise_computed=False))
        if state=='review_ready':
            module=api.load_module('cpu_counterpoise_followup_review',REVIEWER)
            try:
                result=module.review(raw,out/'pilot_review.json')
            except BaseException as error:
                atomic_write_json(status_path,dict(status='review_failed_preserved',error=repr(error),checked=time.time(),
                    owner=watch_owner,pilot_owner=owner,source_sha256=sources,review_complete=False,full_counterpoise_computed=False));raise
            atomic_write_json(status_path,dict(status='pilot_checkpoint_review_complete_gradient_pending',checked=time.time(),
                owner=watch_owner,pilot_owner=owner,source_sha256=sources,review_complete=True,review_status=result['status'],
                full_counterpoise_computed=False));return
        if state in ('pilot_failed_preserved','pilot_owner_stopped_without_terminal_status'):return
        time.sleep(30)
    atomic_write_json(status_path,dict(status='watch_timeout_pilot_not_promoted',checked=time.time(),owner=watch_owner,source_sha256=sources,
        review_complete=False,full_counterpoise_computed=False))


def launch(raw,out,max_hours):
    if out.exists():raise FileExistsError('Preserve previous watcher output')
    if not raw.is_relative_to((ROOT/'calcs').resolve()) or not out.is_relative_to((ROOT/'calcs').resolve()):
        raise ValueError('Use durable project storage')
    status=json.loads((raw/'status.json').read_text())
    if classify(status,owner_state(status['owner'])) not in ('waiting_for_pilot','waiting_for_owner_exit','review_ready'):
        raise ValueError('Pilot is not eligible for monitoring')
    sources={str(p):api.digest(p) for p in (raw/'launch.json',REVIEWER,Path(__file__).resolve())}
    with allocation_guard():
        capacity=host_capacity(1);affinity=cpu_allocation(1);out.mkdir(parents=True)
        command=[sys.executable,str(Path(__file__).resolve()),'watch','--raw',str(raw),'--out',str(out),'--max-hours',str(max_hours)]
        atomic_write_json(out/'launch.json',dict(status='CPU_pilot_watcher_prepared',raw=str(raw),command=command,
            source_sha256=sources,cpu_affinity=affinity,capacity=capacity,max_hours=max_hours))
        env=dict(os.environ)
        for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='1'
        with (out/'watcher.log').open('x') as log:
            process=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(out/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','watch'))
    p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--max-hours',type=float,default=24);a=p.parse_args()
    if not 0<a.max_hours<=48:raise ValueError('Watcher duration must be between zero and 48 hours')
    raw=a.raw.resolve();out=a.out.resolve()
    if a.action=='launch':print(json.dumps(dict(pid=launch(raw,out,a.max_hours),out=str(out))))
    else:watch(raw,out,a.max_hours)
