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

PSS_PROFILE = 'corrected_hybrid_energy_pss'


def high_water_bytes(status):
    for line in status.splitlines():
        if line.startswith('VmHWM:'):
            fields = line.split()
            if len(fields) != 3 or fields[2] != 'kB':
                raise ValueError('Unexpected Linux RSS high-water units')
            return int(fields[1])*1024
    raise ValueError('Linux RSS high-water mark unavailable')


def validate_resources(bound, free_floor, poll, grace, profile='native'):
    allowed = (0 < bound <= 80) if profile == 'native' else (
        bound in (160, 240) if profile in ('corrected_hybrid_energy', PSS_PROFILE) else False)
    if not allowed or free_floor != 32 or not .05 <= poll <= 1 or not 0 < grace <= 30:
        raise ValueError('Explicit supported RSS profile, 32 GiB free floor and 0.05–1 s sampling required')


def tree_memory(process, proportional=False):
    """Retain raw RSS; optionally count shared descendant pages proportionally.

    Main-process RSS/high-water limits remain independent. PSS is queried only
    for the separately declared hybrid profile; native/RSS behavior is unchanged.
    A disappeared owned child has released its memory. Other observation errors
    propagate so the caller stops rather than dropping an unobserved allocation.
    """
    own = process.memory_info().rss
    children = []; readers = []; total = own
    for child in process.children(recursive=True):
        try:
            rss = child.memory_info().rss; total += rss
            row=dict(pid=child.pid,process_start=child.create_time(),RSS_bytes=rss)
            children.append(row);readers.append((child,row))
        except psutil.NoSuchProcess:
            continue
    value = None
    if proportional:
        # With no descendants, using full self RSS is conservative and avoids
        # scanning a large page table just to count the main process alone.
        value = own if not children else process.memory_full_info().pss
        for child,row in readers:
            try:
                if not child.is_running():raise psutil.NoSuchProcess(child.pid)
                info = child.memory_full_info()
                value += info.pss; row['PSS_bytes'] = info.pss
            except psutil.NoSuchProcess:
                row['exited_before_PSS_observation'] = True
        if value < 0:
            raise ValueError('Nonnegative proportional resident memory required')
    return dict(self_RSS_bytes=own, tree_RSS_bytes=total, tree_PSS_bytes=value, children=children)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--maximum-rss-gib', type=float, required=True)
    p.add_argument('--minimum-free-gib', type=float, default=32)
    p.add_argument('--poll-seconds', type=float, default=.1)
    p.add_argument('--grace-seconds', type=float, default=5)
    p.add_argument('--resource-profile', choices=('native', 'corrected_hybrid_energy', PSS_PROFILE), default='native')
    p.add_argument('program', type=Path)
    p.add_argument('arguments', nargs=argparse.REMAINDER)
    a = p.parse_args()
    try:
        validate_resources(a.maximum_rss_gib, a.minimum_free_gib, a.poll_seconds, a.grace_seconds, a.resource_profile)
    except ValueError as error:
        p.error(str(error))
    a.out.mkdir(parents=True, exist_ok=False)
    process = psutil.Process()
    record = dict(status='effective_memory_guard_running', pid=os.getpid(),
                  process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                  host=os.uname().nodename, maximum_RSS_GiB=a.maximum_rss_gib,
                  minimum_free_GiB=a.minimum_free_gib, CPU_affinity=process.cpu_affinity(),
                  resource_profile=a.resource_profile,
                  memory_accounting='max_self_HWM_tree_PSS' if a.resource_profile==PSS_PROFILE else 'max_self_HWM_tree_RSS',
                  peak_tree_PSS_bytes=0, peak_guard_usage_bytes=0,
                  peak_tree_RSS_bytes=0, self_VmHWM_bytes=0, started=time.time(),
                  polling_seconds=a.poll_seconds, signal_target='this owned process only')
    stopped, breached = threading.Event(), threading.Event()

    def sample():
        started=time.monotonic();observed=tree_memory(process,a.resource_profile==PSS_PROFILE)
        hwm = high_water_bytes(Path('/proc/self/status').read_text())
        free = psutil.virtual_memory().available
        usage=max(hwm,observed['tree_PSS_bytes'] if a.resource_profile==PSS_PROFILE else observed['tree_RSS_bytes'])
        if observed['tree_RSS_bytes']>record['peak_tree_RSS_bytes']:
            record['peak_tree_RSS_sample']=observed
        record.update(peak_tree_RSS_bytes=max(record['peak_tree_RSS_bytes'], observed['tree_RSS_bytes']),
                      self_VmHWM_bytes=max(record['self_VmHWM_bytes'], hwm),
                      peak_tree_PSS_bytes=max(record['peak_tree_PSS_bytes'],observed['tree_PSS_bytes'] or 0),
                      peak_guard_usage_bytes=max(record['peak_guard_usage_bytes'],usage),
                      current_tree_RSS_bytes=observed['tree_RSS_bytes'],
                      current_tree_PSS_bytes=observed['tree_PSS_bytes'],
                      last_children=observed['children'],free_RAM_bytes=free,checked=time.time(),
                      observation_seconds=time.monotonic()-started)
        return usage > a.maximum_rss_gib*2**30 or free < a.minimum_free_gib*2**30

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
