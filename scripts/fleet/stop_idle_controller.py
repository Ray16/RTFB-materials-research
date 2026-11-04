#!/usr/bin/env python3
"""Stop one exact owned CPU controller only at a child-free boundary."""
import argparse
import json
import os
from pathlib import Path
import signal
import time

import psutil

from redox.core.io import atomic_write_json


def stop(pid,start,program,output,seconds):
    if output.exists():raise FileExistsError('Preserve earlier controller ownership checks')
    target=psutil.Process(pid)
    def identity():
        actual=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]
        if actual!=start or str(program) not in target.cmdline() or target.uids().real!=os.getuid():
            raise ValueError('Controller PID/start, exact program or user differs')
    identity();record=dict(host=os.uname().nodename,pid=pid,process_start=start,argv=target.cmdline(),
        started=time.time(),only_controller_signalled=True,scientific_workers_signalled=False)
    stopped=False
    try:
        deadline=time.monotonic()+seconds
        while time.monotonic()<deadline:
            identity();os.kill(pid,signal.SIGSTOP);stopped=True
            children=target.children(recursive=True)
            if not children:
                # SIGSTOP closes the check/spawn race; an owned detached GPU
                # worker is not a child and retains its own gate/checkpoints.
                identity();record['children_at_stop']=[]
                os.kill(pid,signal.SIGTERM);os.kill(pid,signal.SIGCONT);stopped=False
                try:target.wait(timeout=10)
                except psutil.TimeoutExpired:raise RuntimeError('Controller did not stop; no source mutation authorized')
                record.update(status='exact_idle_controller_stopped',finished=time.time());atomic_write_json(output,record);return record
            record['last_children']=[]
            for child in children:
                try:record['last_children'].append(dict(pid=child.pid,argv=child.cmdline()))
                except (psutil.NoSuchProcess,psutil.ZombieProcess):continue
            os.kill(pid,signal.SIGCONT);stopped=False;time.sleep(.5)
        raise RuntimeError('No child-free boundary observed; controller continues unchanged')
    finally:
        if stopped:
            try:identity();os.kill(pid,signal.SIGCONT)
            except (ProcessLookupError,FileNotFoundError,psutil.NoSuchProcess):pass


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pid',type=int,required=True)
    p.add_argument('--start',required=True);p.add_argument('--program',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--seconds',type=float,default=60);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    print(json.dumps(stop(a.pid,a.start,a.program,a.out,a.seconds)))
