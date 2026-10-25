#!/usr/bin/env python
"""Launch exactly two small full-spin controls after local runtime preflight."""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import GATE, allocation_guard, cpu_allocation, free_slots, host_capacity
    host = socket.gethostname().split('.')[0]
    out = args.out.resolve()
    profile_path = ROOT / 'environments/cutensor22/runtime.json'
    profile = json.loads(profile_path.read_text())
    preflight = json.loads((out / (host + '_preflight.json')).read_text())
    if not preflight['ok'] or preflight['python'] != profile['python']:
        raise ValueError('Actual runtime interpreter must have a successful local preflight')
    path = out / 'launch.json'
    if path.exists():
        raise ValueError('Prior launch ownership must be reviewed; no automatic duplicate')
    qualifier = ROOT / 'scripts/validation/qualify_residual_response.py'
    sources = {str(p): api.digest(p) for p in (Path(__file__).resolve(), qualifier,
        ROOT / 'scripts/pipeline/residual_response_solver.py', profile_path)}
    with allocation_guard():
        capacity = host_capacity(4, minimum_ram_gib=96)
        helpers = api.load_module('response_control_existing_RAM_promises', ROOT / 'scripts/pipeline/coordinate_geometry_batch.py')
        if capacity['available_RAM_bytes'] - helpers.remaining_promises() < 64 * 2**30:
            raise ValueError('Existing RAM promises leave no room for both controls')
        cpus = cpu_allocation(4)
        slots = free_slots([host])
        if len(slots) < 2:
            raise ValueError('Two exclusive GPUs required for the two declared controls')
        selected = slots[:2]
        for _, gpu in selected:
            size = float(subprocess.check_output(['nvidia-smi', '--id='+str(gpu),
                '--query-gpu=memory.total', '--format=csv,noheader,nounits'], text=True).strip())
            if size < 15360:
                raise ValueError('Small control requires at least the qualified 16 GB class')
        record = dict(source_sha256=sources, host=host, started=time.time(), workers=[],
            scope='Exactly one full singlet and one full doublet replay, 14 GiB device ceiling each')
        atomic_write_json(path, record)
        for index, (_, gpu) in enumerate(selected):
            mult = index + 1
            caps = [f'{key}=2' for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS',
                'OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')]
            argv = ['taskset','-c',','.join(map(str,cpus[2*index:2*index+2])),str(GATE),'run',str(gpu),
                '--label','redox-residual-response-full-control','--','env',*caps,
                'PYTHONPATH='+str(ROOT/'src'),'LD_LIBRARY_PATH='+profile['library_path'],profile['python'],
                str(ROOT/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(out/'guards'/str(mult)),
                str(ROOT/'scripts/fleet/bounded_host_python_worker.py'),'--out',str(out/'host_guards'/str(mult)),
                '--maximum-rss-gib','16','--minimum-free-gib','32',
                str(ROOT/'scripts/fleet/tensor_runtime_worker.py'),'--profile',str(profile_path),
                '--out',str(out/'runtime'/str(mult)),str(qualifier),'--mult',str(mult),
                '--workspace-GiB','28','--out',str(out/'controls'/str(mult))]
            with (out/f'worker_{mult}.log').open('x') as log:
                child = subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,start_new_session=True)
            record['workers'].append(dict(mult=mult,GPU=gpu,CPU_affinity=cpus[2*index:2*index+2],
                gate_pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv))
            atomic_write_json(path,record)
            time.sleep(.2)
    print(json.dumps(dict(host=host,submitted_controls=len(record['workers']))))


if __name__ == '__main__':
    main()
