#!/usr/bin/env python
"""Enforce a tighter host-memory bound on an already running owned worker.

No scientific inputs, solver settings, reservations or installed files change.
Every signal rechecks UID, PID start time, host and exact campaign/program argv.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import socket
import time

import psutil

from redox.core.io import atomic_write_json


def owned_process(owner, campaign, program):
    if owner['host'].split('.')[0] != socket.gethostname().split('.')[0]:
        raise ValueError('The worker belongs to another host')
    process = psutil.Process(owner['pid'])
    start = Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()[19]
    if (start != owner['process_start'] or process.uids().real != os.getuid()
            or str(campaign) not in process.cmdline() or str(program) not in process.cmdline()):
        raise ValueError('Refuse a reused PID or a different worker')
    return process


def enforce(owner, campaign, program, maximum_rss_bytes, minimum_free_bytes):
    process = owned_process(owner, campaign, program)
    rss = process.memory_info().rss
    available = psutil.virtual_memory().available
    result = dict(checked=time.time(), RSS_bytes=rss, free_RAM_bytes=available,
                  maximum_RSS_bytes=maximum_rss_bytes, minimum_free_RAM_bytes=minimum_free_bytes,
                  target_PID=process.pid, status='owned_worker_within_memory_bound')
    if rss > maximum_rss_bytes or available < minimum_free_bytes:
        # Recheck immediately before signaling; a stale manifest is insufficient.
        owned_process(owner, campaign, program).send_signal(signal.SIGTERM)
        result['status'] = 'owned_worker_memory_stop_requested'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('owner', 'campaign', 'program', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--maximum-rss-gib', type=float, required=True)
    parser.add_argument('--minimum-free-gib', type=float, default=32)
    parser.add_argument('--max-hours', type=float, default=48)
    args = parser.parse_args()
    if not (0 < args.max_hours <= 48 and args.maximum_rss_gib > 0 and args.minimum_free_gib > 0):
        parser.error('Use positive memory bounds and a finite window up to 48 hours')
    owner = json.loads(args.owner.read_text())
    campaign, program = args.campaign.resolve(), args.program.resolve()
    owned_process(owner, campaign, program)
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=1)
        os.sched_setaffinity(0, cpu_allocation(1))
    args.out.mkdir(parents=True, exist_ok=False)
    record = dict(host=socket.gethostname(), pid=os.getpid(),
                  process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                  CPU_affinity=sorted(os.sched_getaffinity(0)), target=owner,
                  campaign=str(campaign), program=str(program), started=time.time(),
                  deadline=time.time()+args.max_hours*3600,
                  maximum_RSS_bytes=int(args.maximum_rss_gib*2**30),
                  minimum_free_RAM_bytes=int(args.minimum_free_gib*2**30),
                  foreign_processes_signaled=False, scientific_settings_changed=False)
    atomic_write_json(args.out/'owner.json', record)
    try:
        while time.time() < record['deadline']:
            try:
                state = enforce(owner, campaign, program, record['maximum_RSS_bytes'],
                                record['minimum_free_RAM_bytes'])
            except (FileNotFoundError, psutil.NoSuchProcess):
                record['status'] = 'target_ended_memory_guard_finished'
                break
            atomic_write_json(args.out/'state.json', state)
            if state['status'] != 'owned_worker_within_memory_bound':
                record['status'] = state['status']
                break
            time.sleep(5)
        else:
            owned_process(owner, campaign, program).send_signal(signal.SIGTERM)
            record['status'] = 'finite_window_owned_stop_requested'
    finally:
        record['finished'] = time.time()
        atomic_write_json(args.out/'owner.json', record)


if __name__ == '__main__':
    main()
