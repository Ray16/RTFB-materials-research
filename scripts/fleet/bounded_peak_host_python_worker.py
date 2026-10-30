#!/usr/bin/env python3
"""Add an effective same-process RSS bound without changing a scientific worker.

The ordinary inner scientific guard remains active. This outer guard samples
the process tree more frequently and records Linux's self-RSS high-water mark.
Only this owned worker is signalled on a resource breach.
"""
import argparse
import os
from pathlib import Path
import runpy
import signal
import sys
import threading
import time

import psutil

from redox.core.io import atomic_write_json


def high_water_bytes(status):
    for line in status.splitlines():
        if line.startswith('VmHWM:'):
            fields = line.split()
            if len(fields) != 3 or fields[2] != 'kB':
                raise ValueError('Unexpected Linux RSS high-water units')
            return int(fields[1])*1024
    raise ValueError('Linux RSS high-water mark unavailable')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--maximum-rss-gib', type=float, required=True)
    p.add_argument('--minimum-free-gib', type=float, default=32)
    p.add_argument('--poll-seconds', type=float, default=.1)
    p.add_argument('--grace-seconds', type=float, default=5)
    p.add_argument('program', type=Path)
    p.add_argument('arguments', nargs=argparse.REMAINDER)
    a = p.parse_args()
    if (not 0 < a.maximum_rss_gib <= 80 or a.minimum_free_gib != 32
            or not .05 <= a.poll_seconds <= 1 or not 0 < a.grace_seconds <= 30):
        p.error('A positive effective bound up to 80 GiB, 32 GiB free floor and 0.05–1 s sampling are required')
    a.out.mkdir(parents=True, exist_ok=False)
    process = psutil.Process()
    record = dict(status='effective_memory_guard_running', pid=os.getpid(),
                  process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                  host=os.uname().nodename, maximum_RSS_GiB=a.maximum_rss_gib,
                  minimum_free_GiB=a.minimum_free_gib, CPU_affinity=process.cpu_affinity(),
                  peak_tree_RSS_bytes=0, self_VmHWM_bytes=0, started=time.time(),
                  polling_seconds=a.poll_seconds, signal_target='this owned process only')
    stopped, breached = threading.Event(), threading.Event()

    def sample():
        own = process.memory_info().rss
        children = 0
        for child in process.children(recursive=True):
            try:
                children += child.memory_info().rss
            except psutil.NoSuchProcess:
                continue
        hwm = high_water_bytes(Path('/proc/self/status').read_text())
        free = psutil.virtual_memory().available
        record.update(peak_tree_RSS_bytes=max(record['peak_tree_RSS_bytes'], own+children),
                      self_VmHWM_bytes=max(record['self_VmHWM_bytes'], hwm),
                      current_tree_RSS_bytes=own+children, free_RAM_bytes=free, checked=time.time())
        return max(own+children, hwm) > a.maximum_rss_gib*2**30 or free < a.minimum_free_gib*2**30

    def watch():
        requested = None
        while not stopped.is_set():
            try:
                if sample():
                    record['status'] = 'effective_memory_limit_breached'
                    breached.set()
                atomic_write_json(a.out/'status.json', record)
            except Exception as error:
                record.update(status='effective_memory_observation_failed', error=repr(error))
                breached.set()
                try:
                    atomic_write_json(a.out/'status.json', record)
                except OSError:
                    pass
            if breached.is_set():
                if requested is None:
                    requested = time.monotonic()
                    record['graceful_stop_requested'] = time.time()
                    os.kill(os.getpid(), signal.SIGTERM)
                elif time.monotonic()-requested >= a.grace_seconds:
                    record.update(status='effective_memory_guard_owned_forced_stop',
                                  checkpoint_grace_seconds=a.grace_seconds)
                    try:
                        atomic_write_json(a.out/'status.json', record)
                    finally:
                        os.kill(os.getpid(), signal.SIGKILL)
                    return
            stopped.wait(a.poll_seconds)

    def interrupt(signum, frame):
        raise KeyboardInterrupt('Effective host guard stopped its owned process; preserve checkpoints')

    previous_signal = signal.signal(signal.SIGTERM, interrupt)
    previous_argv = sys.argv
    thread = threading.Thread(target=watch, daemon=True)
    try:
        thread.start()
        sys.argv = [str(a.program), *a.arguments]
        runpy.run_path(str(a.program), run_name='__main__')
        if sample() or breached.is_set():
            raise RuntimeError('Effective host-memory profile failed its observed-peak check')
        record['status'] = 'effective_memory_guard_finished_profile_review_pending'
    except BaseException as error:
        record.update(status='effective_memory_guard_failed_preserved', error=repr(error))
        raise
    finally:
        stopped.set()
        thread.join(timeout=2)
        signal.signal(signal.SIGTERM, previous_signal)
        sys.argv = previous_argv
        record['finished'] = time.time()
        atomic_write_json(a.out/'status.json', record)


if __name__ == '__main__':
    main()
