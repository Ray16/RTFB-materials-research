#!/usr/bin/env python3
"""Review three completed CPU terms, check gradients and assemble the gas sum."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import psutil

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

REVIEWER=ROOT/'scripts/validation/review_counterpoise_cpu_remaining.py'
STATIONARITY=ROOT/'scripts/validation/check_counterpoise_remaining_stationarity.py'
COLLECTOR=ROOT/'scripts/validation/collect_counterpoise_cpu_terms.py'
EVALUATOR=ROOT/'scripts/validation/check_counterpoise_remaining_terms.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/remaining_water_check.json'
ROLES=('AB','B_own','B_ghost');GIB=1024**3


def live_owner(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify remote owner')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return False
    return fields[19]==owner['process_start'] and fields[0] not in ('Z','X')


def classify_upstream(status,live):
    state=status.get('status')
    if state=='remaining_terms_completed_checkpoint_review_pending':return 'waiting_for_owner_exit' if live else 'ready'
    if state in ('remaining_term_failure_or_ownership_held','remaining_continuation_failed_preserved',
                 'remaining_continuation_timeout_preserved','upstream_failed_preserved','upstream_owner_stopped_without_result'):
        return 'upstream_failed_preserved'
    if state in ('waiting','remaining_terms_in_progress_unqualified'):
        return 'waiting' if live else 'upstream_owner_stopped_without_result'
    raise ValueError('Unknown remaining-terms continuation state')


def reserved_future_gib(running):
    total=0
    for process in running.values():
        try:rss=psutil.Process(process.pid).memory_info().rss
        except psutil.NoSuchProcess:rss=0
        total+=max(0,32*GIB-rss)
    return math.ceil(total/GIB)


def water_control():
    data=json.loads(CONTROL.read_text());api.verify_sources(data['source_sha256'])
    if (data['status']!='remaining_water_counterpoise_terms_stationary' or not data['passed']
            or len(data['terms'])!=3 or {x['role'] for x in data['terms']}!=set(ROLES)
            or not all(x['passed'] for x in data['terms']) or data['source_sha256'].get(str(EVALUATOR))!=api.digest(EVALUATOR)):
        raise ValueError('All source-bound remaining-term water controls required')
    return data


def watch(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity,prior,remaining_root,out,max_hours,threads):
    launch=json.loads((out/'launch.json').read_text());sources=launch['source_sha256'];deadline=time.time()+max_hours*3600
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]);own_affinity=sorted(os.sched_getaffinity(0));permitted=launch['permitted_cpu_affinity'];status_path=out/'status.json'
    while time.time()<deadline:
        api.verify_sources(sources);upstream=json.loads((prior/'status.json').read_text());live=live_owner(upstream['owner']) if isinstance(upstream.get('owner'),dict) else False
        stage=classify_upstream(upstream,live)
        atomic_write_json(status_path,dict(status=stage,checked=time.time(),owner=owner,source_sha256=sources,
            stationarity_complete=False,full_counterpoise_computed=False))
        if stage=='ready':break
        if stage not in ('waiting','waiting_for_owner_exit'):return
        time.sleep(30)
    else:
        atomic_write_json(status_path,dict(status='remaining_review_timeout_preserved',checked=time.time(),owner=owner,source_sha256=sources,
            stationarity_complete=False,full_counterpoise_computed=False));return
    water_control();reviewer=api.load_module('remaining_term_independent_reviewer',REVIEWER)
    for role in ROLES:
        folder=out/role;folder.mkdir();result=reviewer.review(remaining_root/role,folder/'review.json')
        if result['status']!='CPU_remaining_term_checkpoint_and_sampled_resources_reviewed_gradient_pending':
            raise ValueError('Term checkpoint review did not pass: '+role)
        api.verify_sources(result['source_sha256'])
    running={};completed={};pending=list(ROLES)
    while time.time()<deadline:
        api.verify_sources(sources)
        for role,process in list(running.items()):
            code=process.poll()
            if code is None:continue
            result_path=out/role/'stationarity/summary.json'
            if code!=0 or not result_path.is_file():
                atomic_write_json(status_path,dict(status='remaining_stationarity_failed_preserved',checked=time.time(),owner=owner,
                    failed_role=role,returncode=code,running_pids={k:v.pid for k,v in running.items() if k!=role},
                    source_sha256=sources,stationarity_complete=False,full_counterpoise_computed=False));return
            result=json.loads(result_path.read_text())
            if result['status']!='CPU_remaining_term_energy_and_orbital_stationarity_checked' or not result['passed'] or result['role']!=role:
                raise ValueError('Term stationarity output is not qualified: '+role)
            completed[role]=api.digest(result_path);del running[role]
        if not pending and not running:
            collector=api.load_module('remaining_five_term_collector',COLLECTOR)
            result=collector.collect(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity,remaining_root,out,out/'counterpoise.json')
            atomic_write_json(status_path,dict(status='five_CPU_gas_counterpoise_terms_stationary',checked=time.time(),owner=owner,
                result_sha256=api.digest(out/'counterpoise.json'),source_sha256=sources,stationarity_complete=True,
                full_counterpoise_computed=True,solvent_counterpoise_computed=False));return
        while pending:
            role=pending[0];required=80+reserved_future_gib(running)+32;affinity=None
            with allocation_guard():
                try:
                    try:
                        os.sched_setaffinity(0,set(permitted));capacity=host_capacity(threads,minimum_ram_gib=required)
                        choices=cpu_allocation(threads+len(own_affinity));affinity=[c for c in choices if c not in own_affinity][:threads]
                        if len(affinity)!=threads:raise RuntimeError('Insufficient disjoint CPU allocation')
                    finally:os.sched_setaffinity(0,set(own_affinity))
                except RuntimeError as error:
                    reasons=('Insufficient shared-node CPU capacity','Insufficient available host RAM',
                             'Requested CPU allocation is unavailable','Insufficient disjoint CPU allocation')
                    if not any(str(error).startswith(reason) for reason in reasons):raise
                if affinity is not None:
                    folder=out/role;scratch=folder/'scratch';scratch.mkdir(exist_ok=False);destination=folder/'stationarity';env=dict(os.environ,PYSCF_TMPDIR=str(scratch))
                    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]=str(threads)
                    cmd=[sys.executable,str(STATIONARITY),'--raw',str(remaining_root/role),'--review',str(folder/'review.json'),
                         '--out',str(destination),'--threads',str(threads)]
                    with (folder/'stationarity.log').open('x') as log:
                        process=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                            env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
                    atomic_write_json(folder/'stationarity_launch.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),
                        cpu_affinity=affinity,capacity=capacity,command=cmd,source_sha256=sources))
                    running[role]=process;pending.pop(0)
            if affinity is None:break
        atomic_write_json(status_path,dict(status='remaining_stationarity_in_progress_unqualified',checked=time.time(),owner=owner,
            pending=pending,running_pids={k:v.pid for k,v in running.items()},completed=completed,source_sha256=sources,
            stationarity_complete=False,full_counterpoise_computed=False));time.sleep(30)
    atomic_write_json(status_path,dict(status='remaining_review_timeout_preserved',checked=time.time(),owner=owner,
        running_pids={k:v.pid for k,v in running.items()},source_sha256=sources,stationarity_complete=False,full_counterpoise_computed=False))


def launch(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity,prior,remaining_root,out,max_hours,threads):
    if out.exists():raise FileExistsError('Preserve previous review continuation')
    for p in (own,own_stationarity,ghost,ghost_stationarity,prior,remaining_root,out):
        if not p.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project storage')
    water_control();sources={str(p):api.digest(p) for p in (REVIEWER,STATIONARITY,COLLECTOR,EVALUATOR,CONTROL,prior/'launch.json',own/'launch.json',case/'plan.json',geometry_review,Path(__file__).resolve())}
    with allocation_guard():
        capacity=host_capacity(1);permitted=sorted(os.sched_getaffinity(0));affinity=cpu_allocation(1);out.mkdir(parents=True)
        cmd=[sys.executable,str(Path(__file__).resolve()),'watch','--case',str(case),'--geometry-review',str(geometry_review),'--own',str(own),
            '--own-stationarity',str(own_stationarity),'--ghost',str(ghost),'--ghost-stationarity',str(ghost_stationarity),
            '--prior',str(prior),'--remaining-root',str(remaining_root),'--out',str(out),'--max-hours',str(max_hours),'--threads',str(threads)]
        atomic_write_json(out/'launch.json',dict(status='remaining_CPU_review_continuation_prepared',command=cmd,source_sha256=sources,
            permitted_cpu_affinity=permitted,cpu_affinity=affinity,capacity=capacity,per_term_sampled_RSS_guard_GiB=32,headroom_GiB=80))
        env=dict(os.environ)
        for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
        with (out/'controller.log').open('x') as log:
            process=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(out/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','watch'))
    for name in ('case','geometry-review','own','own-stationarity','ghost','ghost-stationarity','prior','remaining-root','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=48);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if not 0<a.max_hours<=72 or not 1<=a.threads<=80:raise ValueError('Invalid duration or thread count')
    args=[getattr(a,k).resolve() for k in ('case','geometry_review','own','own_stationarity','ghost','ghost_stationarity','prior','remaining_root','out')]
    if a.action=='launch':print(json.dumps(dict(pid=launch(*args,a.max_hours,a.threads),out=str(args[-1]))))
    else:
        try:watch(*args,a.max_hours,a.threads)
        except BaseException as error:
            atomic_write_json(args[-1]/'status.json',dict(status='remaining_review_continuation_failed_preserved',checked=time.time(),
                error=repr(error),stationarity_complete=False,full_counterpoise_computed=False));raise
