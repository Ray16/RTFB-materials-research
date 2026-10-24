#!/usr/bin/env python
"""Launch the two prepared restoration controls on exclusive local GPUs.

Use the system reservation policy, shared CPU allocation and same-process
residency guard. Each prepared control has one actual attempt; interrupted work
is preserved and is not automatically retried.
"""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.polymer_campaign import digest, verify_sources
from redox.workflows.resources import (GATE, allocation_guard, cpu_allocation,
    free_slots, gpu_suitability, host_capacity)


CAPS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')


def launch(restoration, campaign, out, preflight):
    installation = json.loads((restoration / 'status.json').read_text())
    if installation['status'] != 'fresh_environment_installed_GPU_qualification_pending':
        raise ValueError('Dependency-checked installation required')
    check = json.loads(preflight.read_text())
    host = socket.gethostname().split('.')[0]
    if not check['ok'] or check['host'].split('.')[0] != host:
        raise ValueError('Successful preflight of the actual local node required')
    prefix, source = Path(installation['prefix']), Path(installation['source_checkout'])
    if Path(check['python']).resolve() != (prefix / 'bin/python').resolve():
        raise ValueError('Preflight used another interpreter')
    plan = json.loads((campaign / 'plan.json').read_text())
    verify_sources(plan['producer_sha256'])
    verify_sources(plan['input_sha256'])
    expected = {'formaldehyde', 'methyl_radical'}
    if {t['case'] for t in plan['tasks']} != expected or len(plan['tasks']) != 2:
        raise ValueError('Exactly the two predeclared baseline controls are required')
    if check.get('errors'):
        raise ValueError('Preflight contains unresolved errors')
    launches = out / 'launches'
    launches.mkdir(parents=True, exist_ok=True)
    pending = [t for t in plan['tasks'] if not (launches / (t['case'] + '.json')).exists()]
    records = []
    with allocation_guard():
        capacity = host_capacity(2 * len(pending), minimum_ram_gib=32)
        # The small controls reserve 16 GiB host RAM each. Existing large native
        # promises are included rather than treating free GPUs as free RAM.
        promises = 0
        workroot = ROOT / 'results/polymer_models/campaigns/2026-10-23_crest_refinement'
        for path in workroot.glob('*/thermal/workers/*.json'):
            row = json.loads(path.read_text())
            if row.get('host', '').split('.')[0] != host or row.get('status') != 'computing':
                continue
            proc = Path('/proc') / str(row['pid']) / 'stat'
            if not proc.exists() or proc.read_text().rsplit(')', 1)[1].split()[19] != row['process_start']:
                continue
            scientific_plan = json.loads((path.parents[2] / 'plan.json').read_text())
            promises += scientific_plan.get('host_required_bytes', 256 * 2**30)
        if promises + len(pending) * 16 * 2**30 > .75 * capacity['available_RAM_bytes']:
            raise RuntimeError('Shared host RAM cannot cover existing promises and controls')
        slots = []
        for name, gpu in free_slots([host]):
            try:
                device = gpu_suitability(gpu)
            except (RuntimeError, subprocess.SubprocessError):
                continue
            slots.append((gpu, device))
        if len(slots) < len(pending):
            raise RuntimeError('Insufficient exclusive suitable devices for the prepared controls')
        for task, (gpu, device) in zip(pending, slots):
            cpus = cpu_allocation(2)
            case = task['case']
            runner = ROOT / 'scripts/fleet/exclusive_python_worker.py'
            worker = ROOT / 'scripts/validation/qualify_restored_quantum.py'
            argv = ['taskset', '--cpu-list', ','.join(map(str, cpus)), str(GATE),
                'run', str(gpu), '--label', 'redox-fresh-env-' + case, '--', 'env',
                *[k + '=2' for k in CAPS], 'PYTHONNOUSERSITE=1',
                'PYTHONPATH=' + str(source / 'src'), 'LD_LIBRARY_PATH=' + str(prefix / 'lib'),
                str(prefix / 'bin/python'), '-s', str(runner), '--guard-out', str(out / 'guards' / case),
                str(worker), 'worker', '--restoration', str(restoration),
                '--campaign', str(campaign), '--case', case]
            with (out / (case + '.log')).open('x') as log:
                process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=log,
                    stderr=subprocess.STDOUT, start_new_session=True)
            report = dict(host=host, gate_pid=process.pid,
                gate_process_start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
                GPU=gpu, CPU_affinity=cpus, device=device, argv=argv, started=time.time(),
                host_capacity=capacity, existing_host_promises_bytes=promises,
                preflight_sha256=digest(preflight), producer_sha256={str(p):digest(p) for p in
                    (Path(__file__).resolve(), runner, runner.with_name('watch_gpu_exclusivity.py'), worker)})
            atomic_write_json(launches / (case + '.json'), report)
            records.append(report)
            # Allocation observation sees the newly bound process before the
            # next worker can be assigned an overlapping CPU set.
            time.sleep(.5)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('restoration', 'campaign', 'out', 'preflight'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    result = launch(args.restoration.resolve(), args.campaign.resolve(),
        args.out.resolve(), args.preflight.resolve())
    print(json.dumps(dict(launched=len(result), starts_require_verification=True)))


if __name__ == '__main__':
    main()
