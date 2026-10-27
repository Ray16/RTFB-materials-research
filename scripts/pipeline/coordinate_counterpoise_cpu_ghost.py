#!/usr/bin/env python3
"""Start one CPU ghost-basis resource pilot only after A-own stationarity passes."""
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

GHOST=ROOT/'scripts/pipeline/counterpoise_cpu_ghost_pilot.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/ghost_water_check.json'


def classify(status,owner_live):
    state=status.get('status')
    if state=='pilot_stationarity_complete_no_full_interaction' and status.get('stationarity_complete') is True:
        return 'ready'
    if state in ('stationarity_failed_preserved','stationarity_continuation_failed_preserved',
                 'continuation_timeout_not_promoted','prior_stage_failed_preserved','review_watcher_stopped_without_result'):
        return 'upstream_failed_preserved'
    if state in ('waiting','ready','waiting_for_CPU_capacity','stationarity_process_started_unqualified'):
        return 'waiting' if owner_live else 'upstream_owner_stopped_without_result'
    raise ValueError('Unknown stationarity continuation state')


def live_owner(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify remote owner')
    path=Path('/proc')/str(owner['pid'])/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return False
    return fields[19]==owner['process_start'] and fields[0] not in ('Z','X')


def water_control():
    data=json.loads(CONTROL.read_text());api.verify_sources(data['source_sha256'])
    if (data['status']!='water_ghost_energy_and_orbital_stationarity_checked' or not data['passed']
            or data['role']!='A_ghost' or data['orbital_gradient_norm_au']>data['orbital_gradient_threshold_au']
            or data['energy_error_Ha']>1e-8):raise ValueError('Qualified water ghost control required')
    return data


def watch(case,review,own,prior,ghost_raw,out,max_hours,threads):
    launch=json.loads((out/'launch.json').read_text());sources=launch['source_sha256'];deadline=time.time()+max_hours*3600
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    own_affinity=sorted(os.sched_getaffinity(0));status_path=out/'status.json'
    while time.time()<deadline:
        api.verify_sources(sources);status=json.loads((prior/'status.json').read_text());state=classify(status,live_owner(status['owner']))
        atomic_write_json(status_path,dict(status=state,checked=time.time(),owner=owner,source_sha256=sources,
            ghost_started=False,full_counterpoise_computed=False))
        if state=='ready':
            stationarity=prior/'stationarity/summary.json'
            if api.digest(stationarity)!=status['result_sha256']:raise ValueError('Stationarity result changed')
            water_control();module=api.load_module('cpu_ghost_gated_launcher',GHOST)
            try:
                pid=module.launch(case,review,own,stationarity,ghost_raw,threads,
                    exclude_cpus=own_affinity,allowed_cpus=launch['permitted_cpu_affinity'])
            except RuntimeError as error:
                reasons=('Insufficient shared-node CPU capacity','Insufficient available host RAM',
                         'Requested CPU allocation is unavailable','Insufficient disjoint CPU allocation')
                if ghost_raw.exists() or not any(str(error).startswith(reason) for reason in reasons):raise
                atomic_write_json(status_path,dict(status='waiting_for_CPU_capacity',checked=time.time(),owner=owner,
                    reason=str(error),source_sha256=sources,ghost_started=False,full_counterpoise_computed=False))
                time.sleep(30);continue
            atomic_write_json(status_path,dict(status='ghost_pilot_started_unqualified',checked=time.time(),owner=owner,
                ghost_pid=pid,ghost_raw=str(ghost_raw),source_sha256=sources,ghost_started=True,
                full_counterpoise_computed=False));return
        if state!='waiting':return
        time.sleep(30)
    atomic_write_json(status_path,dict(status='ghost_continuation_timeout_not_promoted',checked=time.time(),owner=owner,
        source_sha256=sources,ghost_started=False,full_counterpoise_computed=False))


def launch(case,review,own,prior,ghost_raw,out,max_hours,threads):
    if out.exists() or ghost_raw.exists():raise FileExistsError('Preserve prior controller and ghost attempts')
    for path in (own,prior,ghost_raw,out):
        if not path.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable shared storage')
    water_control();sources={str(p):api.digest(p) for p in (GHOST,CONTROL,prior/'launch.json',own/'launch.json',
        case/'plan.json',review,Path(__file__).resolve())}
    with allocation_guard():
        capacity=host_capacity(1);permitted=sorted(os.sched_getaffinity(0));affinity=cpu_allocation(1);out.mkdir(parents=True)
        cmd=[sys.executable,str(Path(__file__).resolve()),'watch','--case',str(case),'--review',str(review),'--own',str(own),
             '--prior',str(prior),'--ghost-raw',str(ghost_raw),'--out',str(out),'--max-hours',str(max_hours),'--threads',str(threads)]
        atomic_write_json(out/'launch.json',dict(status='CPU_ghost_continuation_prepared',command=cmd,source_sha256=sources,
            permitted_cpu_affinity=permitted,cpu_affinity=affinity,capacity=capacity))
        env=dict(os.environ)
        for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[k]='1'
        with (out/'controller.log').open('x') as log:
            process=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,
                env=env,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(out/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','watch'))
    for name in ('case','review','own','prior','ghost-raw','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=24);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if not 0<a.max_hours<=48 or not 1<=a.threads<=80:raise ValueError('Invalid duration or thread count')
    case=a.case.resolve();review=a.review.resolve();own=a.own.resolve();prior=a.prior.resolve();ghost_raw=a.ghost_raw.resolve();out=a.out.resolve()
    if a.action=='launch':print(json.dumps(dict(pid=launch(case,review,own,prior,ghost_raw,out,a.max_hours,a.threads),out=str(out))))
    else:
        try:watch(case,review,own,prior,ghost_raw,out,a.max_hours,a.threads)
        except BaseException as error:
            atomic_write_json(out/'status.json',dict(status='ghost_continuation_failed_preserved',checked=time.time(),
                error=repr(error),ghost_started=False,full_counterpoise_computed=False));raise
