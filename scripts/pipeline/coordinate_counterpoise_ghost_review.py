#!/usr/bin/env python3
"""Review and check stationarity after a gated CPU ghost pilot finishes."""
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

REVIEWER=ROOT/'scripts/validation/review_counterpoise_cpu_ghost.py'
STATIONARITY=ROOT/'scripts/validation/check_cpu_counterpoise_ghost_stationarity.py'
WATER=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/ghost_stationarity_water_check_v2.json'
WATER_PRIOR=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/ghost_water_check.json'


def owner_live(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify remote owner')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return False
    return fields[19]==owner['process_start'] and fields[0] not in ('Z','X')


def classify_upstream(status,live):
    state=status.get('status')
    if state=='ghost_pilot_started_unqualified' and status.get('ghost_started') is True:return 'ghost_started'
    if state in ('upstream_failed_preserved','upstream_owner_stopped_without_result',
                 'ghost_continuation_failed_preserved','ghost_continuation_timeout_not_promoted'):
        return 'upstream_failed_preserved'
    if state in ('waiting','ready','waiting_for_CPU_capacity'):
        return 'waiting' if live else 'upstream_owner_stopped_without_result'
    raise ValueError('Unknown ghost continuation state')


def classify_ghost(status,live):
    state=status.get('status')
    if state=='CPU_ghost_pilot_converged_independent_review_pending':return 'waiting_for_owner_exit' if live else 'review_ready'
    if state=='CPU_ghost_pilot_interrupted_or_failed_preserved':return 'ghost_failed_preserved'
    if state in ('CPU_ghost_pilot_initializing_unqualified','CPU_ghost_pilot_SCF_running_unqualified'):
        return 'waiting_for_ghost' if live else 'ghost_owner_stopped_without_terminal_status'
    raise ValueError('Unknown ghost pilot state')


def water_control():
    data=json.loads(WATER.read_text());prior=json.loads(WATER_PRIOR.read_text());api.verify_sources(prior['source_sha256'])
    if (data['status']!='exact_ghost_stationarity_evaluator_water_control_passed'
            or data['evaluator_sha256']!=api.digest(STATIONARITY)
            or data['prior_control_sha256']!=api.digest(WATER_PRIOR) or not prior['passed']
            or data['energy_error_Ha']>1e-8 or data['checkpoint_energy_error_Ha']>1e-8
            or data['orbital_gradient_norm_au']>data['orbital_gradient_threshold_au']):
        raise ValueError('Matching ghost-stationarity water control required')
    return data


def watch(prior,ghost_raw,out,max_hours,threads):
    launch=json.loads((out/'launch.json').read_text());sources=launch['source_sha256'];deadline=time.time()+max_hours*3600
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    own_affinity=sorted(os.sched_getaffinity(0));permitted=launch['permitted_cpu_affinity'];status_path=out/'status.json'
    while time.time()<deadline:
        api.verify_sources(sources);upstream=json.loads((prior/'status.json').read_text())
        upstream_live=owner_live(upstream['owner']) if isinstance(upstream.get('owner'),dict) else False
        stage=classify_upstream(upstream,upstream_live)
        if stage=='ghost_started':
            if Path(upstream['ghost_raw']).resolve()!=ghost_raw:raise ValueError('Upstream launched a different ghost target')
            ghost_status_path=ghost_raw/'status.json'
            if not ghost_status_path.exists():
                pid=upstream['ghost_pid'];cmd=Path('/proc')/str(pid)/'cmdline'
                try:contents=cmd.read_bytes()
                except FileNotFoundError:stage='ghost_owner_stopped_without_status'
                else:
                    if str(ghost_raw).encode() not in contents:raise ValueError('Ghost PID does not match the intended raw path')
                    stage='waiting_for_ghost_startup'
            else:
                status=json.loads(ghost_status_path.read_text())
                if status['owner']['pid']!=upstream['ghost_pid']:raise ValueError('Ghost owner differs from launch')
                stage=classify_ghost(status,owner_live(status['owner']))
        atomic_write_json(status_path,dict(status=stage,checked=time.time(),owner=owner,source_sha256=sources,
            ghost_stationarity_complete=False,full_counterpoise_computed=False))
        if stage=='review_ready':
            reviewer=api.load_module('ghost_pilot_checkpoint_review',REVIEWER);result=reviewer.review(ghost_raw,out/'ghost_review.json')
            if result['status']!='CPU_ghost_pilot_checkpoint_and_sampled_resources_reviewed_gradient_pending':
                raise ValueError('Ghost checkpoint review was not accepted')
            api.verify_sources(result['source_sha256']);water_control()
            while time.time()<deadline:
                with allocation_guard():
                    try:
                        try:
                            os.sched_setaffinity(0,set(permitted));capacity=host_capacity(threads,minimum_ram_gib=80)
                            choices=cpu_allocation(threads+len(own_affinity));affinity=[c for c in choices if c not in own_affinity][:threads]
                            if len(affinity)!=threads:raise RuntimeError('Insufficient disjoint CPU allocation')
                        finally:os.sched_setaffinity(0,set(own_affinity))
                    except RuntimeError as error:
                        reasons=('Insufficient shared-node CPU capacity','Insufficient available host RAM',
                                 'Requested CPU allocation is unavailable','Insufficient disjoint CPU allocation')
                        if not any(str(error).startswith(reason) for reason in reasons):raise
                        atomic_write_json(status_path,dict(status='waiting_for_CPU_capacity',checked=time.time(),owner=owner,
                            reason=str(error),source_sha256=sources,ghost_stationarity_complete=False,full_counterpoise_computed=False))
                        affinity=None
                    if affinity is not None:
                        scratch=out/'scratch';scratch.mkdir(exist_ok=False);destination=out/'stationarity';env=dict(os.environ,PYSCF_TMPDIR=str(scratch))
                        for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[k]=str(threads)
                        cmd=[sys.executable,str(STATIONARITY),'--raw',str(ghost_raw),'--review',str(out/'ghost_review.json'),
                             '--out',str(destination),'--threads',str(threads)]
                        with (out/'stationarity.log').open('x') as log:
                            process=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                                env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
                        atomic_write_json(out/'stationarity_launch.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),
                            cpu_affinity=affinity,capacity=capacity,command=cmd,source_sha256=sources))
                if affinity is None:time.sleep(30);continue
                atomic_write_json(status_path,dict(status='ghost_stationarity_process_started_unqualified',checked=time.time(),owner=owner,
                    stationarity_pid=process.pid,cpu_affinity=affinity,source_sha256=sources,ghost_stationarity_complete=False,full_counterpoise_computed=False))
                code=process.wait();result_path=destination/'summary.json'
                if code!=0 or not result_path.is_file():
                    atomic_write_json(status_path,dict(status='ghost_stationarity_failed_preserved',returncode=code,checked=time.time(),owner=owner,
                        source_sha256=sources,ghost_stationarity_complete=False,full_counterpoise_computed=False));return
                checked=json.loads(result_path.read_text())
                if checked['status']!='CPU_ghost_energy_and_orbital_stationarity_checked' or not checked['passed']:
                    raise ValueError('Ghost stationarity output not qualified')
                atomic_write_json(status_path,dict(status='ghost_stationarity_complete_no_full_interaction',checked=time.time(),owner=owner,
                    ghost_stationarity_complete=True,result_sha256=api.digest(result_path),source_sha256=sources,full_counterpoise_computed=False));return
            break
        if stage not in ('waiting','waiting_for_ghost_startup','waiting_for_ghost','waiting_for_owner_exit'):return
        time.sleep(30)
    atomic_write_json(status_path,dict(status='ghost_review_timeout_not_promoted',checked=time.time(),owner=owner,
        source_sha256=sources,ghost_stationarity_complete=False,full_counterpoise_computed=False))


def launch(prior,ghost_raw,out,max_hours,threads):
    if out.exists():raise FileExistsError('Preserve prior ghost-review continuation')
    for path in (prior,ghost_raw,out):
        if not path.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project storage')
    water_control();sources={str(p):api.digest(p) for p in (REVIEWER,STATIONARITY,WATER,WATER_PRIOR,prior/'launch.json',Path(__file__).resolve())}
    with allocation_guard():
        capacity=host_capacity(1);permitted=sorted(os.sched_getaffinity(0));affinity=cpu_allocation(1);out.mkdir(parents=True)
        cmd=[sys.executable,str(Path(__file__).resolve()),'watch','--prior',str(prior),'--ghost-raw',str(ghost_raw),
             '--out',str(out),'--max-hours',str(max_hours),'--threads',str(threads)]
        atomic_write_json(out/'launch.json',dict(status='ghost_review_continuation_prepared',command=cmd,source_sha256=sources,
            permitted_cpu_affinity=permitted,cpu_affinity=affinity,capacity=capacity))
        env=dict(os.environ)
        for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
        with (out/'controller.log').open('x') as log:
            process=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(out/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','watch'))
    for name in ('prior','ghost-raw','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=24);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if not 0<a.max_hours<=48 or not 1<=a.threads<=80:raise ValueError('Invalid duration or thread count')
    prior=a.prior.resolve();ghost_raw=a.ghost_raw.resolve();out=a.out.resolve()
    if a.action=='launch':print(json.dumps(dict(pid=launch(prior,ghost_raw,out,a.max_hours,a.threads),out=str(out))))
    else:
        try:watch(prior,ghost_raw,out,a.max_hours,a.threads)
        except BaseException as error:
            atomic_write_json(out/'status.json',dict(status='ghost_review_continuation_failed_preserved',checked=time.time(),
                error=repr(error),ghost_stationarity_complete=False,full_counterpoise_computed=False));raise
