#!/usr/bin/env python3
"""Freeze live process/start, CPU/thread and GPU observations without restarting."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json

REMOTE = '''import json,os,pathlib,subprocess,sys
owners=json.loads(sys.argv[1]);result=[]
names=("OMP_NUM_THREADS","OMP_THREAD_LIMIT","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS","XTB_THREADS","CUDA_VISIBLE_DEVICES")
for owner in owners:
 p=pathlib.Path("/proc")/str(owner["pid"]);live=False
 try:live=p.joinpath("stat").read_text().rsplit(")",1)[1].split()[19]==owner["process_start"]
 except FileNotFoundError:pass
 row=dict(owner=owner,live=live)
 if live:
  values=dict(x.split("=",1) for x in p.joinpath("environ").read_text().split("\\0") if "=" in x)
  row.update(CPU_affinity=sorted(os.sched_getaffinity(owner["pid"])),environment={k:values.get(k) for k in names})
 result.append(row)
gpu=subprocess.run(["nvidia-smi","--query-compute-apps=pid,gpu_uuid,used_memory","--format=csv,noheader,nounits"],capture_output=True,text=True,timeout=20,check=True)
print(json.dumps(dict(processes=result,GPU_residency=gpu.stdout.splitlines())))
'''


def observe(manifest, output, extra_namespaces=()):
    if output.exists():
        raise FileExistsError('Preserve earlier live observations')
    data = json.loads(manifest.read_text()); records = []; owners = {}
    namespaces = [Path(r['namespace']) for r in data['states']] + list(extra_namespaces)
    for namespace in namespaces:
        spec = json.loads(namespace.read_text()); case = Path(spec['case']); path = case / 'status.json'
        if not path.exists():
            continue
        stage = json.loads(path.read_text()); owner = stage['owner']; host = owner['host'].split('.')[0]
        owners.setdefault(host, []).append(owner)
        progresses = []
        for p in (case / 'raw/SCF_checkpoints').glob('*/*/progress.json'):
            progresses.append(dict(path=str(p), snapshot=json.loads(p.read_text())))
        records.append(dict(namespace=str(namespace), namespace_sha256=sha256_file(namespace),
            stage_snapshot=stage, SCF_iterations=progresses,
            SCF_kernels_entered=len(list((case / 'raw/SCF_checkpoints').glob('*/*/settings.json'))),
            occupied_checkpoints_written=len(list((case / 'raw/SCF_checkpoints').glob('*/*/occupied.h5'))),
            native_G_qualified=False))
    def remote(item):
        host, entries = item
        try:
            run = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', host,
                shlex.join(['/usr/bin/python3', '-c', REMOTE, json.dumps(entries)])],
                capture_output=True, text=True, check=True, timeout=40)
            return host, json.loads(run.stdout)
        except (OSError, subprocess.SubprocessError) as error:
            return host, dict(status='observation_unverified_no_restart', error=repr(error))
    with ThreadPoolExecutor(max_workers=max(1, len(owners))) as pool:
        live = dict(pool.map(remote, owners.items()))
    result = dict(status='sampling_process_and_progress_observation_not_physical_admission',
        checked=time.time(), manifest=str(manifest), manifest_sha256=sha256_file(manifest),
        records=records, live_hosts=live, native_G_qualified=False, jobs_restarted=False,
        producer_sha256=sha256_file(__file__))
    atomic_write_json(output, result); return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--extra-namespace', type=Path, action='append', default=[]); a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    result = observe(a.manifest.resolve(), a.out.resolve(), [p.resolve() for p in a.extra_namespace])
    print(json.dumps(dict(status=result['status'], observed_scientific_records=len(result['records']),
        live_processes=sum(p['live'] for host in result['live_hosts'].values() for p in host.get('processes', [])),
        SCF_kernels_entered=sum(r['SCF_kernels_entered'] for r in result['records']))))
