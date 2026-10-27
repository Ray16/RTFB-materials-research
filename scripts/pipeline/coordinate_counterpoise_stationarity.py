#!/usr/bin/env python3
"""Run the CPU stationarity check after a source-bound pilot review completes."""
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

STATIONARITY=ROOT/'scripts/validation/check_cpu_counterpoise_stationarity.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/stationarity_water_check.json'


def classify(status, owner_live):
    state=status.get('status')
    if state=='pilot_checkpoint_review_complete_gradient_pending' and status.get('review_complete') is True:
        return 'ready'
    if state in ('pilot_failed_preserved','pilot_owner_stopped_without_terminal_status','review_failed_preserved','watch_timeout_pilot_not_promoted'):
        return 'prior_stage_failed_preserved'
    if state in ('waiting_for_pilot','waiting_for_owner_exit','review_ready'):
        return 'waiting' if owner_live else 'review_watcher_stopped_without_result'
    raise ValueError('Unknown review-watcher state')


def live_owner(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify remote watcher ownership')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return False
    return fields[19]==owner['process_start'] and fields[0] not in ('Z','X')


def water_control():
    data=json.loads(CONTROL.read_text());api.verify_sources(data['source_sha256'])
    if (data['status']!='water_counterpoise_stationarity_recomputed' or not data['passed']
            or data['orbital_gradient_norm_au']>data['orbital_gradient_threshold_au']
            or data['reference_energy_error_Ha']>1e-8):
        raise ValueError('Matching water stationarity control required')
    return data


def watch(raw,prior,out,max_hours,threads):
    launch=json.loads((out/'launch.json').read_text());sources=launch['source_sha256'];deadline=time.time()+max_hours*3600
    watcher_affinity=sorted(os.sched_getaffinity(0));original_affinity=launch['permitted_cpu_affinity']
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status_path=out/'status.json'
    while time.time()<deadline:
        api.verify_sources(sources);prior_status=json.loads((prior/'status.json').read_text())
        state=classify(prior_status,live_owner(prior_status['owner']))
        atomic_write_json(status_path,dict(status=state,checked=time.time(),owner=owner,source_sha256=sources,
            stationarity_complete=False,full_counterpoise_computed=False))
        if state=='ready':
            review=prior/'pilot_review.json';data=json.loads(review.read_text())
            if data['status']!='CPU_pilot_checkpoint_and_sampled_resources_reviewed_gradient_pending':
                raise ValueError('Independent checkpoint review is missing or incompatible')
            api.verify_sources(data['source_sha256']);water_control()
            with allocation_guard():
                try:
                    try:
                        os.sched_setaffinity(0,set(original_affinity))
                        capacity=host_capacity(threads,minimum_ram_gib=80)
                        choices=cpu_allocation(threads+len(watcher_affinity))
                        affinity=[c for c in choices if c not in watcher_affinity][:threads]
                        if len(affinity)!=threads:raise RuntimeError('Insufficient disjoint CPU allocation')
                    finally:
                        os.sched_setaffinity(0,set(watcher_affinity))
                except RuntimeError as error:
                    reasons=('Insufficient shared-node CPU capacity','Insufficient available host RAM',
                             'Requested CPU allocation is unavailable','Insufficient disjoint CPU allocation')
                    if not any(str(error).startswith(reason) for reason in reasons):raise
                    atomic_write_json(status_path,dict(status='waiting_for_CPU_capacity',checked=time.time(),owner=owner,
                        reason=str(error),source_sha256=sources,stationarity_complete=False,full_counterpoise_computed=False))
                    affinity=None
                if affinity is not None:
                    scratch=out/'scratch';scratch.mkdir(exist_ok=False)
                    destination=out/'stationarity';env=dict(os.environ,PYSCF_TMPDIR=str(scratch))
                    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]=str(threads)
                    command=[sys.executable,str(STATIONARITY),'--pilot-raw',str(raw),'--pilot-review',str(review),
                             '--pilot-out',str(destination),'--threads',str(threads)]
                    with (out/'stationarity.log').open('x') as log:
                        process=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                            env=env,cwd=ROOT,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
                    atomic_write_json(out/'stationarity_launch.json',dict(pid=process.pid,host=os.uname().nodename,
                        launched=time.time(),cpu_affinity=affinity,capacity=capacity,command=command,source_sha256=sources))
            if affinity is None:
                time.sleep(30)
                continue
            atomic_write_json(status_path,dict(status='stationarity_process_started_unqualified',checked=time.time(),owner=owner,
                stationarity_pid=process.pid,cpu_affinity=affinity,source_sha256=sources,stationarity_complete=False,
                full_counterpoise_computed=False))
            code=process.wait();result_path=destination/'summary.json'
            if code!=0 or not result_path.is_file():
                atomic_write_json(status_path,dict(status='stationarity_failed_preserved',returncode=code,checked=time.time(),
                    owner=owner,source_sha256=sources,stationarity_complete=False,full_counterpoise_computed=False));return
            result=json.loads(result_path.read_text())
            if result['status']!='CPU_pilot_energy_and_orbital_stationarity_checked' or not result['passed']:
                raise ValueError('Stationarity output is not qualified')
            atomic_write_json(status_path,dict(status='pilot_stationarity_complete_no_full_interaction',checked=time.time(),owner=owner,
                stationarity_complete=True,result_sha256=api.digest(result_path),source_sha256=sources,
                full_counterpoise_computed=False));return
        if state!='waiting':return
        time.sleep(30)
    atomic_write_json(status_path,dict(status='continuation_timeout_not_promoted',checked=time.time(),owner=owner,
        source_sha256=sources,stationarity_complete=False,full_counterpoise_computed=False))


def launch(raw,prior,out,max_hours,threads):
    if out.exists():raise FileExistsError('Preserve previous continuation')
    for path in (raw,prior,out):
        if not path.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project storage')
    water_control();prior_launch=prior/'launch.json';raw_launch=raw/'launch.json'
    sources={str(p):api.digest(p) for p in (STATIONARITY,CONTROL,prior_launch,raw_launch,Path(__file__).resolve())}
    with allocation_guard():
        capacity=host_capacity(1);permitted=sorted(os.sched_getaffinity(0));affinity=cpu_allocation(1);out.mkdir(parents=True)
        cmd=[sys.executable,str(Path(__file__).resolve()),'watch','--raw',str(raw),'--prior',str(prior),'--out',str(out),
             '--max-hours',str(max_hours),'--threads',str(threads)]
        atomic_write_json(out/'launch.json',dict(status='stationarity_continuation_prepared',raw=str(raw),prior=str(prior),
            permitted_cpu_affinity=permitted,cpu_affinity=affinity,source_sha256=sources,capacity=capacity,command=cmd))
        env=dict(os.environ)
        for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
        with (out/'controller.log').open('x') as log:
            process=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(out/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','watch'))
    for name in ('raw','prior','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=24);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if not 0<a.max_hours<=48 or not 1<=a.threads<=80:raise ValueError('Invalid duration or thread count')
    raw=a.raw.resolve();prior=a.prior.resolve();out=a.out.resolve()
    if a.action=='launch':print(json.dumps(dict(pid=launch(raw,prior,out,a.max_hours,a.threads),out=str(out))))
    else:
        try:watch(raw,prior,out,a.max_hours,a.threads)
        except BaseException as error:
            atomic_write_json(out/'status.json',dict(status='stationarity_continuation_failed_preserved',
                checked=time.time(),error=repr(error),stationarity_complete=False,full_counterpoise_computed=False))
            raise
