#!/usr/bin/env python3
"""Gated, resource-bounded CPU dispatch of the three remaining gas terms."""
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

WORKER=ROOT/'scripts/pipeline/counterpoise_cpu_remaining.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/remaining_water_check.json'
ROLES=('AB','B_own','B_ghost')
GIB=1024**3


def live_owner(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify remote owner')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return False
    return fields[19]==owner['process_start'] and fields[0] not in ('Z','X')


def classify_upstream(status,live):
    state=status.get('status')
    if state=='ghost_stationarity_complete_no_full_interaction' and status.get('ghost_stationarity_complete') is True:return 'ready'
    if state in ('ghost_review_continuation_failed_preserved','ghost_stationarity_failed_preserved',
                 'ghost_failed_preserved','ghost_owner_stopped_without_terminal_status',
                 'ghost_owner_stopped_without_status','upstream_failed_preserved',
                 'upstream_owner_stopped_without_result','ghost_review_timeout_not_promoted'):
        return 'upstream_failed_preserved'
    if state in ('waiting','waiting_for_ghost_startup','waiting_for_ghost','waiting_for_owner_exit',
                 'review_ready','waiting_for_CPU_capacity','ghost_stationarity_process_started_unqualified'):
        return 'waiting' if live else 'upstream_owner_stopped_without_result'
    raise ValueError('Unknown ghost qualification state')


def classify_term(raw):
    if not raw.exists():return dict(status='pending')
    status_path=raw/'status.json';launcher_path=raw/'launcher.json'
    if not status_path.exists():
        if not launcher_path.exists():return dict(status='held_missing_launch_record')
        launcher=json.loads(launcher_path.read_text());pid=launcher['pid'];path=Path('/proc')/str(pid)/'cmdline'
        try:command=path.read_bytes()
        except FileNotFoundError:return dict(status='held_owner_stopped_without_status',pid=pid)
        if str(raw).encode() not in command:return dict(status='held_pid_identity_mismatch',pid=pid)
        return dict(status='starting',pid=pid)
    status=json.loads(status_path.read_text());state=status.get('status');owner=status.get('owner')
    if state=='CPU_remaining_term_interrupted_or_failed_preserved':return dict(status='failed_preserved',owner=owner)
    if not isinstance(owner,dict):return dict(status='held_missing_owner')
    live=live_owner(owner)
    if state=='CPU_remaining_term_converged_independent_review_pending':
        return dict(status='finishing' if live else 'completed_review_pending',owner=owner)
    if state in ('CPU_remaining_term_initializing_unqualified','CPU_remaining_term_SCF_running_unqualified'):
        return dict(status='running' if live else 'held_owner_stopped_without_terminal_status',owner=owner)
    return dict(status='held_unknown_state',owner=owner)


def future_reservation_gib(states):
    """Reserve each live worker's unconsumed sampled 32-GiB ceiling."""
    total=0
    for state in states.values():
        if state['status'] not in ('starting','running','finishing'):continue
        owner=state.get('owner');pid=owner['pid'] if owner else state['pid']
        try:rss=psutil.Process(pid).memory_info().rss
        except psutil.NoSuchProcess:rss=0
        total+=max(0,32*GIB-rss)
    return math.ceil(total/GIB)


def water_control():
    data=json.loads(CONTROL.read_text());api.verify_sources(data['source_sha256'])
    if (data['status']!='remaining_water_counterpoise_terms_stationary' or not data['passed']
            or len(data['terms'])!=3 or {x['role'] for x in data['terms']}!=set(ROLES)
            or not all(x['passed'] for x in data['terms'])):raise ValueError('All remaining water controls required')
    return data


def watch(case,review,own,own_stationarity,ghost,prior,raw_root,out,max_hours,threads):
    launch=json.loads((out/'launch.json').read_text());sources=launch['source_sha256'];deadline=time.time()+max_hours*3600
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]);own_affinity=sorted(os.sched_getaffinity(0));permitted=launch['permitted_cpu_affinity'];status_path=out/'status.json'
    while time.time()<deadline:
        api.verify_sources(sources);upstream=json.loads((prior/'status.json').read_text());prior_live=live_owner(upstream['owner']) if isinstance(upstream.get('owner'),dict) else False
        stage=classify_upstream(upstream,prior_live)
        if stage!='ready':
            atomic_write_json(status_path,dict(status=stage,checked=time.time(),owner=owner,source_sha256=sources,terms_started=0,full_counterpoise_computed=False))
            if stage!='waiting':return
            time.sleep(30);continue
        ghost_stationarity=prior/'stationarity/summary.json'
        if api.digest(ghost_stationarity)!=upstream['result_sha256']:raise ValueError('Ghost stationarity record changed')
        water_control();worker=api.load_module('remaining_cpu_term_launcher',WORKER)
        states={role:classify_term(raw_root/role) for role in ROLES}
        if any(x['status'].startswith(('failed','held')) for x in states.values()):
            atomic_write_json(status_path,dict(status='remaining_term_failure_or_ownership_held',checked=time.time(),owner=owner,
                source_sha256=sources,terms=states,full_counterpoise_computed=False));return
        if all(x['status']=='completed_review_pending' for x in states.values()):
            atomic_write_json(status_path,dict(status='remaining_terms_completed_checkpoint_review_pending',checked=time.time(),owner=owner,
                source_sha256=sources,terms=states,full_counterpoise_computed=False));return
        for role in ROLES:
            if states[role]['status']!='pending':continue
            required=80+future_reservation_gib(states)+32
            try:
                pid=worker.launch(role,case,review,own,own_stationarity,ghost,ghost_stationarity,raw_root/role,threads,
                    exclude_cpus=own_affinity,allowed_cpus=permitted,minimum_ram_gib=required)
            except RuntimeError as error:
                reasons=('Insufficient shared-node CPU capacity','Insufficient available host RAM',
                         'Requested CPU allocation is unavailable','Insufficient disjoint CPU allocation')
                if (raw_root/role).exists() or not any(str(error).startswith(reason) for reason in reasons):raise
                break
            states[role]=dict(status='starting',pid=pid)
        atomic_write_json(status_path,dict(status='remaining_terms_in_progress_unqualified',checked=time.time(),owner=owner,
            source_sha256=sources,terms=states,full_counterpoise_computed=False));time.sleep(30)
    atomic_write_json(status_path,dict(status='remaining_continuation_timeout_preserved',checked=time.time(),owner=owner,
        source_sha256=sources,full_counterpoise_computed=False))


def launch(case,review,own,own_stationarity,ghost,prior,raw_root,out,max_hours,threads):
    if out.exists() or raw_root.exists():raise FileExistsError('Preserve previous continuation and term outputs')
    for p in (own,own_stationarity,ghost,prior,raw_root,out):
        if not p.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project storage')
    water_control();sources={str(p):api.digest(p) for p in (WORKER,CONTROL,prior/'launch.json',own/'launch.json',case/'plan.json',review,Path(__file__).resolve())}
    with allocation_guard():
        capacity=host_capacity(1);permitted=sorted(os.sched_getaffinity(0));affinity=cpu_allocation(1);out.mkdir(parents=True);raw_root.mkdir(parents=True)
        cmd=[sys.executable,str(Path(__file__).resolve()),'watch','--case',str(case),'--review',str(review),'--own',str(own),
            '--own-stationarity',str(own_stationarity),'--ghost',str(ghost),'--prior',str(prior),'--raw-root',str(raw_root),
            '--out',str(out),'--max-hours',str(max_hours),'--threads',str(threads)]
        atomic_write_json(out/'launch.json',dict(status='remaining_CPU_terms_continuation_prepared',command=cmd,source_sha256=sources,
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
    for name in ('case','review','own','own-stationarity','ghost','prior','raw-root','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=48);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if not 0<a.max_hours<=72 or not 1<=a.threads<=80:raise ValueError('Invalid duration or thread count')
    args=[getattr(a,k).resolve() for k in ('case','review','own','own_stationarity','ghost','prior','raw_root','out')]
    if a.action=='launch':print(json.dumps(dict(pid=launch(*args,a.max_hours,a.threads),out=str(args[-1]))))
    else:
        try:watch(*args,a.max_hours,a.threads)
        except BaseException as error:
            atomic_write_json(args[-1]/'status.json',dict(status='remaining_continuation_failed_preserved',checked=time.time(),
                error=repr(error),full_counterpoise_computed=False));raise
