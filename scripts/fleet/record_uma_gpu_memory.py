#!/usr/bin/env python
"""Sample the owned UMA process GPU memory without importing CUDA libraries."""
import argparse
import csv
import math
import os
from pathlib import Path
import runpy
import subprocess
import sys
import threading
import time

from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner


def own_memory(output, uuid, pid):
    values = []
    for row in csv.reader(output.splitlines(), skipinitialspace=True):
        if len(row) != 3:
            raise ValueError('Malformed GPU process memory response')
        device, process, memory = [v.strip() for v in row]
        if device == uuid and int(process) == pid:
            value = float(memory)
            if not math.isfinite(value) or value < 0:
                raise ValueError('Invalid GPU process memory')
            values.append(value)
    return sum(values) if values else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('program', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    verify_gate_owner()
    gpu = os.environ['CUDA_VISIBLE_DEVICES']
    uuid = subprocess.check_output(['nvidia-smi','--id='+gpu,'--query-gpu=uuid',
                                    '--format=csv,noheader'],text=True,timeout=10).strip()
    if not uuid or '\n' in uuid:
        raise ValueError('Cannot identify the reserved GPU')
    args.out.mkdir(parents=True, exist_ok=False)
    pid = os.getpid()
    record = dict(status='GPU_memory_observation_running',pid=pid,GPU=gpu,uuid=uuid,
                  started=time.time(),samples=0,own_samples=0,peak_process_memory_MiB=None,
                  sampling_errors=[])
    stop = threading.Event()

    def observe():
        while not stop.is_set():
            try:
                raw = subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_gpu_memory',
                                               '--format=csv,noheader,nounits'],text=True,timeout=10)
                memory = own_memory(raw,uuid,pid)
                record['samples'] += 1
                if memory is not None:
                    record['own_samples'] += 1
                    record['peak_process_memory_MiB'] = max(record['peak_process_memory_MiB'] or 0,memory)
            except (OSError,ValueError,subprocess.SubprocessError) as error:
                record['sampling_errors'].append(repr(error))
            record['checked'] = time.time()
            atomic_write_json(args.out/'status.json',record)
            stop.wait(2)

    thread = threading.Thread(target=observe,daemon=True)
    previous = sys.argv
    try:
        thread.start()
        sys.argv = [str(args.program),*args.arguments]
        runpy.run_path(str(args.program),run_name='__main__')
        record['status'] = 'GPU_memory_observation_program_returned_review_required'
    except BaseException as error:
        record.update(status='GPU_memory_observation_program_failed_preserved',error=repr(error))
        raise
    finally:
        stop.set()
        thread.join(timeout=12)
        sys.argv = previous
        record['finished'] = time.time()
        atomic_write_json(args.out/'status.json',record)


if __name__ == '__main__':
    main()
