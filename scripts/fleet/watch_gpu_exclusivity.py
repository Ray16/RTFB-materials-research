#!/usr/bin/env python
"""Audit this project's GPU workers and stop only owned conflicting processes.

The shared reservation gate remains authoritative for allocation. This separate
CPU observer detects delayed foreign CUDA residency after launch. It never
reserves/releases GPUs, clears claims, restarts work or signals foreign jobs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import time

import psutil

from redox.core.io import atomic_write_json


def process_start(pid):
    return Path(f"/proc/{pid}/stat").read_text().rsplit(")",1)[1].split()[19]


def conflicts(pid, gpu, devices, residents):
    if gpu not in devices:
        raise ValueError("Cannot identify the assigned physical GPU")
    return sorted({p for p, uuid in residents if uuid == devices[gpu] and p != pid})


def stop_owned(pid, start, expected_campaign):
    p = psutil.Process(pid)
    if (p.uids().real != os.getuid() or process_start(pid) != start
            or str(expected_campaign) not in p.cmdline()):
        raise ValueError("Do not signal a reused PID or an unrelated command")
    descendants = [(q.pid, process_start(q.pid)) for q in p.children(recursive=True)]
    p.send_signal(signal.SIGTERM)
    try:
        p.wait(timeout=10)
    except psutil.TimeoutExpired:
        return dict(status="owned_worker_termination_unconfirmed_preserve_claims", stopped=False)
    for child, child_start in descendants:
        try:
            if process_start(child) == child_start:
                psutil.Process(child).send_signal(signal.SIGTERM)
        except (FileNotFoundError, psutil.NoSuchProcess):
            pass
    return dict(status="owned_conflicting_worker_stopped_checkpoints_preserved", stopped=True,
                descendants_checked=len(descendants))


def snapshot():
    def query(arguments):
        return list(csv.reader(subprocess.check_output(["nvidia-smi", *arguments,
            "--format=csv,noheader,nounits"],text=True,timeout=15).splitlines(),skipinitialspace=True))
    devices = {int(i): uuid.strip() for i, uuid in query(["--query-gpu=index,uuid"])}
    residents = [(int(pid), uuid.strip()) for pid, uuid in query(["--query-compute-apps=pid,gpu_uuid"])]
    return devices, residents


def audit(root, out):
    host = socket.gethostname().split('.')[0]
    devices, residents = snapshot()
    checked, events, unresolved = [], [], []
    for path in sorted(root.glob("*/*/workers/*.json")):
        try:
            row = json.loads(path.read_text())
        except (OSError,ValueError) as error:
            unresolved.append(dict(manifest=str(path),error=repr(error)))
            continue
        if row.get("host", "").split('.')[0] != host:
            continue
        pid, start = row.get("pid"), row.get("process_start")
        if not isinstance(pid,int) or not start:
            continue
        try:
            p = psutil.Process(pid)
            if process_start(pid) != start or p.status() == psutil.STATUS_ZOMBIE:
                continue
            campaign = path.parents[2].resolve()
            if str(campaign) not in p.cmdline():
                raise ValueError("Worker PID does not run the recorded campaign")
            gpu = int(row["gpu"])
            other = conflicts(pid, gpu, devices, residents)
            checked.append(dict(manifest=str(path),pid=pid,process_start=start,gpu=gpu,
                GPU_residency_observed=any(r==pid for r,_ in residents),conflicting_PIDs=other))
            if other:
                event = dict(checked=time.time(),host=host,manifest=str(path),pid=pid,
                    process_start=start,gpu=gpu,foreign_resident_PIDs=other,
                    foreign_processes_signaled=False,claims_or_GPU_reservations_cleared=False)
                event.update(stop_owned(pid,start,campaign))
                events.append(event)
                with (out / "events.jsonl").open("a") as f:
                    f.write(json.dumps(event)+"\n")
        except (FileNotFoundError, psutil.NoSuchProcess):
            continue
        except (OSError, ValueError, psutil.AccessDenied) as error:
            unresolved.append(dict(manifest=str(path),error=repr(error)))
    return dict(checked=time.time(),host=host,live_workers=checked,events=events,
                unresolved=unresolved,scientific_outputs_qualified=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work-root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--deadline",type=float,required=True)
    p.add_argument("--poll-seconds",type=float,default=10)
    args = p.parse_args()
    if not 0 < args.deadline-time.time() <= 48*3600 or not 2 <= args.poll_seconds <= 60:
        p.error("Use a finite deadline within 48 hours and a 2–60 second poll")
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():
        host_capacity(1,minimum_ram_gib=1)
        os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False)
    owner=dict(host=socket.gethostname(),pid=os.getpid(),process_start=process_start(os.getpid()),
        CPU_affinity=sorted(os.sched_getaffinity(0)),deadline=args.deadline,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    atomic_write_json(args.out/'owner.json',owner)
    def stop(signum, frame):
        raise KeyboardInterrupt(f"observer signal {signum}; scientific workers unchanged")
    signal.signal(signal.SIGTERM,stop)
    try:
        while time.time() < args.deadline:
            try:
                state=audit(args.work_root.resolve(),args.out)
            except (OSError,subprocess.SubprocessError,ValueError) as error:
                state=dict(checked=time.time(),status="residency_observation_unresolved_no_restart",error=repr(error))
            atomic_write_json(args.out/'state.json',state)
            time.sleep(args.poll_seconds)
    finally:
        owner.update(finished=time.time(),status="finite_exclusivity_observer_stopped")
        atomic_write_json(args.out/'owner.json',owner)


if __name__=='__main__':
    main()
