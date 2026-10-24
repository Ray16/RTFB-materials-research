#!/usr/bin/env python
"""Run a checkpoint-aware Python worker with an owned-process host RAM bound."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--maximum-rss-gib', type=float, required=True)
    parser.add_argument('--minimum-free-gib', type=float, default=32)
    parser.add_argument('program', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.maximum_rss_gib <= 0 or args.minimum_free_gib <= 0:
        parser.error('Positive resource bounds required')
    args.out.mkdir(parents=True, exist_ok=False)
    process = psutil.Process()
    report = dict(pid=os.getpid(), process_start=process.create_time(), host=os.uname().nodename,
        CPU_affinity=process.cpu_affinity(), started=time.time(), maximum_RSS_GiB=args.maximum_rss_gib,
        minimum_free_GiB=args.minimum_free_gib, peak_RSS_bytes=0, status='running',
        signal_target='this owned process only', program=str(args.program.resolve()))
    stop = threading.Event()

    def observe():
        while not stop.is_set():
            try:
                rss = process.memory_info().rss
                free = psutil.virtual_memory().available
                report['peak_RSS_bytes'] = max(report['peak_RSS_bytes'], rss)
                report.update(RSS_bytes=rss, free_RAM_bytes=free, checked=time.time())
                breached = rss > args.maximum_rss_gib*2**30 or free < args.minimum_free_gib*2**30
                if breached:
                    report['status'] = 'resource_bound_exceeded_own_stop_requested'
                atomic_write_json(args.out/'status.json', report)
                if breached:
                    os.kill(os.getpid(), signal.SIGTERM)
                    return
            except Exception as error:
                report.update(status='resource_observation_failed_own_stop_requested', error=repr(error))
                atomic_write_json(args.out/'status.json', report)
                os.kill(os.getpid(), signal.SIGTERM)
                return
            stop.wait(5)

    def interrupt(signum, frame):
        raise KeyboardInterrupt('Host RAM guard stopped its owned worker')

    previous_handler = signal.signal(signal.SIGTERM, interrupt)
    previous_argv = sys.argv
    thread = threading.Thread(target=observe, daemon=True)
    try:
        thread.start()
        sys.argv = [str(args.program), *args.arguments]
        runpy.run_path(str(args.program), run_name='__main__')
        report['status'] = 'program_returned_review_required'
    except BaseException as error:
        report.update(status='program_interrupted_or_failed_preserved', error=repr(error))
        raise
    finally:
        stop.set()
        thread.join(timeout=6)
        signal.signal(signal.SIGTERM, previous_handler)
        sys.argv = previous_argv
        report['finished'] = time.time()
        atomic_write_json(args.out/'status.json', report)


if __name__ == '__main__':
    main()
