"""Isolated UMA state queue with proton topology and measured memory gates."""
import argparse
import csv
import gc
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time

from redox.core.io import atomic_write_json
from redox.qm import uma
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.resources import verify_gate_owner

base = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--tasks', required=True)
args = parser.parse_args()
plan = json.loads((base / 'uma_plan.json').read_text())
verify_gate_owner()
for path, digest in plan['source_sha256'].items():
    assert uma.file_sha256(path) == digest, path
import torch
assert torch.cuda.device_count() == 1
assert uma.model_identity(uma.DEFAULT_MODEL) == plan['model_identity']
torch.set_num_threads(1)
torch.cuda.set_per_process_memory_fraction(.75, 0)
device = torch.cuda.get_device_properties(0)
uma.LIBRARY, uma.OUT = base / 'isolated_library', base / 'uma'
with (uma.LIBRARY / 'manifest.csv').open() as stream:
    rows = list(csv.DictReader(stream))
uma.read_manifest = lambda: rows
owner = dict(host=socket.gethostname(), pid=os.getpid(),
             process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
             CPU_affinity=sorted(os.sched_getaffinity(0)),
             physical_gpu=os.environ['CUDA_VISIBLE_DEVICES'])
consecutive = 0


def interrupted(signum, frame):
    raise KeyboardInterrupt(f'signal {signum}')


signal.signal(signal.SIGTERM, interrupted)
for index in map(int, args.tasks.split(',')):
    task = plan['tasks'][index]
    row = task['row']
    out = base / 'tasks' / task['key']
    if (out / 'status.json').exists():
        raise RuntimeError('Existing attempt preserved; review it before resuming')
    for path, digest in {**plan['source_sha256'], **task['seed_sha256']}.items():
        assert uma.file_sha256(path) == digest, path
    out.mkdir(parents=True, exist_ok=False)
    claim = base / 'claims' / task['key']
    claim.mkdir()
    atomic_write_json(claim / 'owner.json', owner)
    report = dict(owner, status='computing', task=task, started=time.time(), attempts=1,
                  device=device.name, total_memory_bytes=device.total_memory,
                  allocation_fraction=.75, native_thread_caps={k:os.environ.get(k) for k in
                    ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
                  solution_properties_validated=False)
    atomic_write_json(out / 'status.json', report)
    stop, samples = threading.Event(), []

    def sample():
        while not stop.is_set():
            try:
                value = subprocess.check_output(['nvidia-smi', '--id='+owner['physical_gpu'],
                    '--query-gpu=memory.used', '--format=csv,noheader,nounits'], text=True, timeout=5)
                samples.append(dict(time=time.time(), used_MiB=float(value.strip())))
            except (OSError, ValueError, subprocess.SubprocessError):
                pass
            stop.wait(.5)

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    try:
        torch.cuda.reset_peak_memory_stats(0)
        sys.argv = ['redox.qm.uma', '--only', row['id']+':'+row['state'], '--device', 'cuda']
        uma.main()
        directory = uma.OUT / row['id'] / row['state']
        result = json.loads((directory / 'result.json').read_text())
        seeds = list(map(Path, task['seed_sha256']))
        assert not uma.cache_errors(result, row, directory/'relaxed.xyz', seeds, plan['model_identity'], .02)
        topology = hydrogen_topology_status(row['smiles'], directory/'relaxed.xyz')
        report['hydrogen_topology'] = topology
        assert topology == 'ok', topology
        for path, digest in task['seed_sha256'].items():
            assert uma.file_sha256(path) == digest
        torch.cuda.synchronize()
        report.update(peak_allocated_bytes=torch.cuda.max_memory_allocated(0),
                      peak_reserved_bytes=torch.cuda.max_memory_reserved(0),
                      peak_device_used_MiB=max((s['used_MiB'] for s in samples), default=None),
                      result=str(directory/'result.json'), result_sha256=uma.file_sha256(directory/'result.json'),
                      geometry_sha256=uma.file_sha256(directory/'relaxed.xyz'))
        assert samples and report['peak_device_used_MiB']*1024**2 < .8*device.total_memory
        assert report['peak_reserved_bytes'] < .8*device.total_memory
        report['status'] = 'gas_UMA_state_and_hydrogen_topology_verified'
        consecutive = 0
    except Exception as exc:
        report.update(status='failed_preserving_outputs', error=repr(exc))
        consecutive += 1
    except BaseException as exc:
        report.update(status='interrupted_preserving_outputs', error=repr(exc))
        raise
    finally:
        stop.set()
        sampler.join(timeout=6)
        report.update(memory_samples=samples, finished=time.time())
        atomic_write_json(out/'status.json', report)
        gc.collect()
        torch.cuda.empty_cache()
        shutil.rmtree(claim)
    if consecutive >= 3:
        raise RuntimeError('Three consecutive actual failures; stop worker for review')
