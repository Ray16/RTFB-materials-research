#!/usr/bin/env python3
"""Stop only exact owned workers when registered GPU exclusivity is violated.

This independent CPU watchdog does not reserve/release GPUs or restart work.
It supplements the historical in-process monitor, whose SIGTERM can be ignored
inside a scientific runtime. Escalation always rechecks UID, PID start and argv.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import signal
import time

import psutil

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def process_start(pid):
    return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]


def exact_owned(pid,start,plan):
    process=psutil.Process(pid)
    if process.uids().real!=os.getuid() or process_start(pid)!=start or str(plan) not in process.cmdline():
        raise ValueError('A reused PID or unrelated command must never be signaled')
    return process


def stop_owned(pid,start,plan,grace=10):
    process=exact_owned(pid,start,plan);process.send_signal(signal.SIGTERM)
    try:process.wait(timeout=grace);return dict(stopped=True,forced=False)
    except psutil.TimeoutExpired:pass
    try:
        process=exact_owned(pid,start,plan);process.kill();process.wait(timeout=10)
        return dict(stopped=True,forced=True)
    except (psutil.NoSuchProcess,FileNotFoundError):return dict(stopped=True,forced=False)
    except psutil.TimeoutExpired:return dict(stopped=False,forced=True,claims_preserved=True)


def run(roster,helper,output,hours):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    output.mkdir(parents=True,exist_ok=False)
    definition=importlib.util.spec_from_file_location('frozen_exclusivity_snapshot',helper)
    monitor=importlib.util.module_from_spec(definition);definition.loader.exec_module(monitor)
    host=os.uname().nodename.split('.')[0];deadline=time.time()+hours*3600;known={}
    owner=dict(host=host,pid=os.getpid(),process_start=process_start(os.getpid()),CPU_affinity=sorted(os.sched_getaffinity(0)),
        deadline=deadline,foreign_processes_signaled=False,source_sha256={str(p):sha256_file(p) for p in (Path(__file__),helper)})
    atomic_write_json(output/'owner.json',owner)
    try:
        while time.time()<deadline:
            checked=[];events=[];unresolved=[]
            try:
                rows=json.loads(roster.read_text());current={r['namespace']:r['namespace_sha256'] for r in rows}
                if len(rows)!=len(current) or any(current.get(k)!=v for k,v in known.items()):
                    raise ValueError('Exclusivity roster is append-only')
                devices,residents=monitor.snapshot()
                for name,digest in current.items():
                    namespace=Path(name)
                    if sha256_file(namespace)!=digest:raise ValueError('Registered namespace changed')
                    spec=json.loads(namespace.read_text());case=Path(spec['case']);path=case/'residency_guard/owner.json'
                    if not path.exists():continue
                    identity=json.loads(path.read_text())
                    if identity['host'].split('.')[0]!=host:continue
                    pid,start=identity['pid'],identity['process_start']
                    try:
                        exact_owned(pid,start,case/'plan.json')
                        other=monitor.conflicts(pid,identity['GPU'],devices,residents)
                        checked.append(dict(namespace=name,pid=pid,process_start=start,GPU=identity['GPU'],foreign_PIDs=other))
                        if not other:continue
                        # Snapshot queries can outlive a worker; revalidate its
                        # exact identity before writing a negative QC marker.
                        exact_owned(pid,start,case/'plan.json')
                        event=dict(checked=time.time(),host=host,namespace=name,namespace_sha256=digest,
                            pid=pid,process_start=start,GPU=identity['GPU'],foreign_resident_PIDs=other,
                            foreign_processes_signaled=False,watchdog=str(output/'owner.json'))
                        marker=case/'residency_guard/external_interruption.json'
                        if not marker.exists():atomic_write_json(marker,event)
                        event.update(stop_owned(pid,start,case/'plan.json'));events.append(event)
                    except (psutil.NoSuchProcess,FileNotFoundError):continue
                    except (OSError,ValueError,psutil.AccessDenied) as error:
                        unresolved.append(dict(namespace=name,error=repr(error)))
                known=current
            except Exception as error:unresolved.append(dict(error=repr(error),no_restart_inferred=True))
            state=dict(checked=time.time(),host=host,live_workers=checked,events=events,unresolved=unresolved,
                foreign_processes_signaled=False,claims_or_reservations_cleared=False)
            atomic_write_json(output/'state.json',state)
            if events:
                with (output/'events.jsonl').open('a') as stream:
                    for event in events:stream.write(json.dumps(event)+'\n')
            time.sleep(5)
    finally:
        owner.update(finished=time.time(),status='namespace_exclusivity_watchdog_stopped');atomic_write_json(output/'owner.json',owner)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('roster','helper','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--hours',type=float,default=24);a=p.parse_args()
    if not 0<a.hours<=72:p.error('Positive bounded watchdog window up to 72 hours required')
    run(a.roster.resolve(),a.helper.resolve(),a.out.resolve(),a.hours)
