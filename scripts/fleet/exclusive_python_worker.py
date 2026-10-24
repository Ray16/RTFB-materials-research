#!/usr/bin/env python
"""Run a Python GPU program in the reserved process with residency checks.

Using runpy keeps the gate parent, scientific PID, CPU affinity and thread caps
unchanged. A confirmed foreign resident signals only this process. The guard's
files must live outside any scientific output tree being hash-frozen.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import runpy
import signal
import socket
import sys
import threading
import time

from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--guard-out',type=Path,required=True)
    p.add_argument('program',type=Path)
    p.add_argument('arguments',nargs=argparse.REMAINDER)
    args=p.parse_args()
    verify_gate_owner()
    helper=Path(__file__).with_name('watch_gpu_exclusivity.py')
    spec=importlib.util.spec_from_file_location('reserved_process_residency',helper)
    monitor=importlib.util.module_from_spec(spec);spec.loader.exec_module(monitor)
    args.guard_out.mkdir(parents=True,exist_ok=False)
    pid=os.getpid();gpu=int(os.environ['CUDA_VISIBLE_DEVICES'])
    report=dict(host=socket.gethostname(),pid=pid,process_start=monitor.process_start(pid),
        GPU=gpu,CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),
        program=str(args.program.resolve()),foreign_processes_signaled=False)
    atomic_write_json(args.guard_out/'owner.json',report)

    def check():
        devices,residents=monitor.snapshot()
        other=monitor.conflicts(pid,gpu,devices,residents)
        state=dict(checked=time.time(),foreign_resident_PIDs=other,
            own_residency_observed=any(p==pid for p,_ in residents),GPU=gpu)
        atomic_write_json(args.guard_out/'state.json',state)
        return other

    if check():
        raise RuntimeError('Foreign GPU residency found before scientific program start')
    stop=threading.Event()
    def observe():
        while not stop.wait(5):
            try:
                other=check()
            except Exception as error:
                atomic_write_json(args.guard_out/'observation_error.json',dict(
                    checked=time.time(),error=repr(error),no_restart_inferred=True))
                continue
            if other:
                atomic_write_json(args.guard_out/'interruption.json',dict(
                    checked=time.time(),pid=pid,foreign_resident_PIDs=other,
                    signal_target='own reserved process only',foreign_processes_signaled=False))
                os.kill(pid,signal.SIGTERM)
                return
    thread=threading.Thread(target=observe,daemon=True);thread.start()
    prior=sys.argv
    try:
        sys.argv=[str(args.program),*args.arguments]
        runpy.run_path(str(args.program),run_name='__main__')
    finally:
        sys.argv=prior;stop.set();thread.join(timeout=2)
        report.update(finished=time.time(),guard_stopped=True)
        atomic_write_json(args.guard_out/'owner.json',report)


if __name__=='__main__':
    main()
