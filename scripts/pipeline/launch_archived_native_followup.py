#!/usr/bin/env python3
"""Preflight and launch explicit archived-native placements, preserving receipts.

GPU leases and CPU/RAM allocation are acquired by the remote namespace launcher.
One host's dependency error stops further placements on that host. This script
does not infer a scientific start or completion from a launcher receipt.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


CAPS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')


def launch_host(host, entries, destination, interpreter):
    rows = []
    for entry in entries:
        case = destination/entry['key']
        case.mkdir()
        prefix = ['env', *(key+'=1' for key in CAPS), 'PYTHONPATH='+str(ROOT/'src')]
        if entry.get('namespace_bin_directory'):
            directory = Path(entry['namespace_bin_directory'])
            if not directory.is_absolute() or '..' in directory.parts:
                raise ValueError('Absolute declared namespace-binary directory required')
            prefix += ['PATH='+str(directory)+':/usr/local/bin:/usr/bin:/bin']
        prefix += [interpreter]
        probe = prefix+['-m', 'redox.workflows.archived_native_namespace', 'probe', '--spec', entry['spec']]
        response = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host,
                                   shlex.join(probe)], capture_output=True, text=True, timeout=240)
        (case/'preflight_stdout.txt').write_text(response.stdout)
        (case/'preflight_stderr.txt').write_text(response.stderr)
        if response.returncode:
            rows.append(dict(key=entry['key'], host=host, status='node_preflight_failed_no_launch',
                             returncode=response.returncode))
            break
        preflight = json.loads(response.stdout)
        if preflight['host'] != host or preflight['status'] != 'original_path_namespace_hashes_match':
            raise ValueError('Preflight node or source identity differs')
        argv = prefix+['-m', 'redox.workflows.archived_native_namespace', 'launch', '--spec', entry['spec'],
                       '--gpu', str(entry['gpu']), '--out', str(case/'placement')]
        code = ('import json,os,subprocess,pathlib; '
                'p=pathlib.Path('+repr(str(case/'controller.log'))+'); '
                'f=p.open("x"); '
                'child=subprocess.Popen('+repr(argv)+',stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); '
                'stat=pathlib.Path(f"/proc/{child.pid}/stat").read_text().rsplit(")",1)[1].split()[19]; '
                'print(json.dumps({"pid":child.pid,"process_start":stat,"host":os.uname().nodename}))')
        response = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host,
                                   shlex.join(['/usr/bin/python3', '-c', code])],
                                  capture_output=True, text=True, timeout=30, check=True)
        owner = json.loads(response.stdout)
        row = dict(key=entry['key'], host=host, GPU=entry['gpu'], owner=owner,
                   spec=entry['spec'], spec_sha256=sha256_file(entry['spec']),
                   argv=argv, started=time.time(), status='controller_spawned_scientific_start_unverified')
        atomic_write_json(case/'controller.json', row)
        rows.append(row)
        print(json.dumps(dict(key=entry['key'], host=host, GPU=entry['gpu'], controller=owner['pid'])), flush=True)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--placements', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--python', required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    rows = json.loads(a.placements.read_text())['entries']
    if len({r['key'] for r in rows}) != len(rows) or len({(r['host'], r['gpu']) for r in rows}) != len(rows):
        raise ValueError('Distinct exact cases and GPU placements required')
    policy = ROOT.parent/'.config/gpu_reserve/config.env'
    import re
    permitted = re.search(r'^ALLOWED_HOSTS="([^"]+)"', policy.read_text(), re.M)[1].split()
    if any(r['host'] not in permitted for r in rows):
        raise ValueError('Requested host is outside the current system policy')
    a.out.mkdir(parents=True, exist_ok=False)
    source = {str(path): sha256_file(path) for path in
              (Path(__file__).resolve(), a.placements.resolve(), policy,
               ROOT/'src/redox/workflows/archived_native_namespace.py')}
    hosts = sorted({r['host'] for r in rows})
    with ThreadPoolExecutor(max_workers=len(hosts)) as pool:
        futures = [pool.submit(launch_host, h, [r for r in rows if r['host'] == h],
                               a.out.resolve(), a.python) for h in hosts]
        result = [row for future in futures for row in future.result()]
    if any(sha256_file(path) != digest for path, digest in source.items()):
        raise ValueError('Placement inputs changed during submission; inspect retained receipts')
    atomic_write_json(a.out/'summary.json', dict(status='explicit_native_placements_submitted_review_pending',
                                               entries=result, source_sha256=source,
                                               actual_scientific_starts_verified=False,
                                               completed_native_free_energies_verified=False))


if __name__ == '__main__':
    main()
