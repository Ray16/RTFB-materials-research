#!/usr/bin/env python
"""Place one alternate-model electronic pilot with the existing fleet guards.

Reuse only resource placement, not primary energy admission or normalization.
An uncertain RPC or an existing scientific attempt is preserved for review.
"""
import argparse
import fcntl
import inspect
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
PROFILE = ROOT / 'environments/cutensor22/runtime.json'
WORKER = ROOT / 'scripts/pipeline/smooth_polymer_energy.py'


def placement_module(profile):
    source = ROOT / 'scripts/pipeline/coordinate_geometry_batch.py'
    module = api.load_module('experimental_energy_resource_contract', source)
    body = inspect.getsource(module.launch_local)
    anchor = "        with (attempt / 'worker.log').open('x') as log:"
    if body.count(anchor) != 1:
        raise ValueError('The validated resource-placement template changed')
    inserted = """        if plan.get('runtime_profile') != str(_PROFILE):
            raise ValueError('Unexpected experimental energy runtime')
        argv.insert(argv.index('env')+1, 'LD_LIBRARY_PATH='+_profile['library_path'])
        index=argv.index(manifest['worker'])
        argv[index:index]=[str(ROOT/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(_PROFILE),
            '--out', str(campaign/'tensor_runtime')]
"""
    scope = dict(vars(module), _PROFILE=PROFILE, _profile=profile); local = {}
    exec(compile(body.replace(anchor, inserted+anchor, 1), str(SCRIPT), 'exec'), scope, local)
    module.launch_local = local['launch_local']
    return module


def preflight(batch, host, profile):
    path = batch / (host + '_preflight.json')
    if path.exists() and time.time() - path.stat().st_mtime < 600:
        data = json.loads(path.read_text())
        return data['ok'] and data['python'] == profile['python']
    helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
    argv = [sys.executable, str(helper), host, '--python', profile['python'],
        '--module', 'pyscf', '--module', 'gpu4pyscf', '--module', 'cupy', '--module', 'ase',
        '--module', 'scipy', '--path', str(batch / 'batch.json'), '--library-path', profile['library_path'],
        '--min-available-ram-mib', '65536']
    result = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    data = json.loads(result.stdout); atomic_write_json(path, data)
    return result.returncode == 0 and data['ok']


def wait(batch, hosts, hours, placement, profile):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    manifest = json.loads((batch / 'batch.json').read_text())
    if (len(manifest['entries']) != 1 or manifest['worker'] != str(WORKER)
            or manifest['primary_candidate_admission'] or manifest['stage'] != 'electronic'):
        raise ValueError('Exactly one declared alternate-model electronic pilot is supported')
    api.verify_sources(manifest['source_sha256'])
    entry = manifest['entries'][0]; campaign = Path(entry['campaign'])
    if api.digest(campaign / 'plan.json') != entry['plan_sha256']:
        raise ValueError('Experimental electronic plan changed')
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    with (batch / 'dispatch.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = batch / 'dispatch.json'
        if status.exists():
            raise ValueError('Preserve the existing dispatcher and scientific ownership; inspect before resuming')
        owner = dict(host=socket.gethostname().split('.')[0], pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
        sources = {**manifest['source_sha256'], str(batch / 'batch.json'): api.digest(batch / 'batch.json'),
            str(SCRIPT): api.digest(SCRIPT)}
        record = dict(status='experimental_energy_pilot_waiting_for_free_capacity', owner=owner,
            started=time.time(), deadline=time.time()+hours*3600, scientific_starts=0,
            source_sha256=sources, CPU_affinity=sorted(os.sched_getaffinity(0)))
        try:
            while time.time() < record['deadline']:
                api.verify_sources(sources)
                path = campaign / 'experimental_energy/status.json'
                if path.exists() or (campaign / 'residency_guard').exists():
                    record['status'] = 'existing_scientific_attempt_preserved'; return
                for host, gpu in free_slots(hosts):
                    if not preflight(batch, host, profile): continue
                    intent = campaign / 'dispatch_intent.json'
                    atomic_write_json(intent, dict(host=host, GPU=gpu, owner=owner, checked=time.time()))
                    command = ['env', 'OMP_NUM_THREADS=1', 'OMP_THREAD_LIMIT=1', 'MKL_NUM_THREADS=1',
                        'OPENBLAS_NUM_THREADS=1', 'NUMEXPR_NUM_THREADS=1', 'PYTHONPATH='+str(ROOT/'src'),
                        profile['python'], str(SCRIPT), 'launch-local', '--batch', str(batch),
                        '--campaign', str(campaign), '--gpu', str(gpu)]
                    try:
                        result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8',
                            host, shlex.join(command)], capture_output=True, text=True, timeout=45, check=True)
                        observation = json.loads(result.stdout)
                    except (OSError, ValueError, subprocess.SubprocessError) as error:
                        record.update(status='submission_ownership_unverified_preserved', error=repr(error)); return
                    if observation['status'] == 'gate_started_scientific_start_pending':
                        record['launch'] = observation
                        for _ in range(12):
                            if path.exists():
                                record.update(status='experimental_electronic_scientific_stage_started', scientific_starts=1)
                                return
                            if not placement.owner_live(observation['owner']): break
                            time.sleep(5)
                        if (campaign/'residency_guard').exists() or placement.owner_live(observation['owner']):
                            record['status'] = 'scientific_start_needs_review_preserved'; return
                        record['status'] = 'gate_exited_before_scientific_stage_review_required'; return
                    else:
                        atomic_write_json(intent, dict(host=host, GPU=gpu, owner=owner,
                            checked=time.time(), outcome='not_started'))
                record['updated'] = time.time(); atomic_write_json(status, record)
                time.sleep(15)
            record['status'] = 'placement_window_ended_no_extra_scientific_attempt'
        except BaseException as error:
            record.update(status='dispatcher_failed_existing_outputs_preserved', error=repr(error)); raise
        finally:
            record['finished'] = time.time(); atomic_write_json(status, record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('wait', 'launch-local'))
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--campaign', type=Path)
    parser.add_argument('--gpu', type=int)
    parser.add_argument('--hosts', default='lambda1,lambda2,lambda11,lambda13')
    parser.add_argument('--max-hours', type=float, default=24)
    args = parser.parse_args(); profile = json.loads(PROFILE.read_text())
    if Path(sys.executable) != Path(profile['python']): parser.error('Use the declared profile interpreter')
    if os.environ.get('LD_LIBRARY_PATH'): parser.error('Keep the CUDA loader path confined to scientific workers')
    if not 0 < args.max_hours <= 48: parser.error('Use a finite window up to 48 hours')
    if any(not h or h.startswith('-') for h in args.hosts.split(',')): parser.error('Invalid hosts')
    placement = placement_module(profile)
    if args.action == 'launch-local':
        if args.campaign is None or args.gpu is None: parser.error('launch-local requires campaign and GPU')
        print(json.dumps(placement.launch_local(args.batch.resolve(), args.campaign.resolve(), args.gpu)))
    else: wait(args.batch.resolve(), args.hosts.split(','), args.max_hours, placement, profile)


if __name__ == '__main__': main()
