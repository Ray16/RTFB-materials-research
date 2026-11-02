#!/usr/bin/env python3
"""Measure shared-fork versus actual private memory without quantum/GPU work.

A small initialized buffer is inherited by one owned child. The child first
shares it, then writes its own copy. RSS and proportional resident memory are
measured in each state; the child is joined before the control returns.
"""
import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import time

import psutil

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def totals(pid, child_pid):
    rows = []
    for identity in (pid, child_pid):
        p = psutil.Process(identity); memory = p.memory_full_info()
        rows.append(dict(pid=identity, process_start=p.create_time(), RSS_bytes=memory.rss,
                         PSS_bytes=memory.pss, USS_bytes=memory.uss))
    return dict(processes=rows, tree_RSS_bytes=sum(r['RSS_bytes'] for r in rows),
                tree_PSS_bytes=sum(r['PSS_bytes'] for r in rows),
                tree_USS_bytes=sum(r['USS_bytes'] for r in rows))


def run():
    size = 96 * 2**20; payload = bytearray(size)
    for i in range(0, size, 4096):
        payload[i] = 1
    parent, child = mp.get_context('fork').Pipe()
    def work():
        parent.close(); child.send('shared')
        if child.recv() != 'copy':
            raise ValueError('Owned control expected its explicit copy command')
        for i in range(0, size, 4096):
            payload[i] = 2
        child.send('private'); child.recv(); child.close()
    process = mp.get_context('fork').Process(target=work); process.start(); child.close()
    try:
        if not parent.poll(10) or parent.recv() != 'shared':
            raise RuntimeError('Shared child did not reach its observation barrier')
        shared = totals(os.getpid(), process.pid); parent.send('copy')
        if not parent.poll(10) or parent.recv() != 'private':
            raise RuntimeError('Private child did not reach its observation barrier')
        private = totals(os.getpid(), process.pid); parent.send('done'); process.join(10)
        if process.is_alive() or process.exitcode != 0:
            raise RuntimeError('Owned memory-control child did not finish normally')
    finally:
        if process.is_alive():
            process.terminate(); process.join(10)
        parent.close()
    delta = private['tree_PSS_bytes'] - shared['tree_PSS_bytes']
    if (not shared['tree_RSS_bytes'] > 1.6 * shared['tree_PSS_bytes']
            or not .85 * size <= delta <= 1.15 * size
            or abs(private['tree_RSS_bytes'] - shared['tree_RSS_bytes']) > .1 * size):
        raise ValueError('Controlled fork/private-memory contrast did not reproduce')
    return dict(status='shared_fork_RSS_double_count_and_private_allocation_independently_observed',
        payload_bytes=size, shared=shared, private=private, private_PSS_increase_bytes=delta,
        child_returncode=process.exitcode, GPU_work=False, failed_quantum_job_cause_proven=False,
        quantum_G_qualified=False, checked=time.time(), producer_sha256=sha256_file(__file__))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    if a.out.exists():
        raise FileExistsError('Preserve earlier memory-accounting observations')
    result = run(); atomic_write_json(a.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'payload_bytes', 'private_PSS_increase_bytes')}))
