#!/usr/bin/env python
"""Finite multi-node dispatch of independent, checkpointed geometry campaigns.

Use every free suitable GPU within measured host-RAM and disjoint CPU budgets.
Actual scientific failures are preserved, never automatically retried. Atomic
gate contention is tracked separately. Existing running workers survive exit.
"""
import argparse
import fcntl
import json
import math
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
PREFLIGHT = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'


def remaining_promises():
    import psutil
    total = 0
    for process in psutil.process_iter():
        try:
            argv = process.cmdline()
            if not argv or 'python' not in Path(argv[0]).name or '--maximum-rss-gib' not in argv:
                continue
            bound = float(argv[argv.index('--maximum-rss-gib') + 1]) * 2**30
            if not math.isfinite(bound) or bound <= 0:
                raise RuntimeError('Invalid existing host resource promise')
            total += max(0, bound - process.memory_info().rss)
        except psutil.NoSuchProcess:
            continue
        except psutil.AccessDenied as error:
            raise RuntimeError('Cannot verify existing host resource promises') from error
    return total


def can_fit(available, promised, bound=16 * 2**30, floor=32 * 2**30):
    if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0
           for v in (available, promised, bound, floor)):
        raise ValueError('Finite nonnegative host resource quantities are required')
    return available - promised >= bound + floor


def launch_local(batch, campaign, gpu):
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation, free_slots, gpu_suitability, GATE
    manifest = json.loads((batch / 'batch.json').read_text())
    api.verify_sources(manifest['source_sha256'])
    entries = [e for e in manifest['entries'] if e['campaign'] == str(campaign)]
    if len(entries) != 1:
        raise ValueError('Campaign is not in the exact declared batch')
    if entries[0]['plan_sha256'] is None and manifest.get('stage') == 'electronic':
        preparer = api.load_module('deferred_electronics_launch_validation', Path(manifest['preparer']))
        preparer.validate_prepared(entries[0], manifest)
    elif api.digest(campaign / 'plan.json') != entries[0]['plan_sha256']:
        raise ValueError('Campaign is not in the exact declared batch')
    plan = api.read_plan(campaign)
    if (plan['solvent'] != 'acetonitrile' or len(plan['tasks']) != 1
            or (plan['tasks'][0]['charge'], plan['tasks'][0]['mult']) != (0, 1)):
        raise ValueError('This execution class requires one neutral singlet in acetonitrile')
    host = socket.gethostname().split('.')[0]
    preflight = json.loads((batch / (host + '_preflight.json')).read_text())
    if not preflight['ok'] or preflight['python'] != sys.executable:
        raise ValueError('Node/interpreter preflight must pass')
    if (campaign / 'geometry').exists() or (campaign / 'residency_guard').exists():
        raise ValueError('Preserve an existing scientific attempt')
    with allocation_guard():
        bound = manifest['maximum_RSS_GiB']
        capacity = host_capacity(2, minimum_ram_gib=bound + 32)
        promised = remaining_promises()
        if not can_fit(capacity['available_RAM_bytes'], promised, bound=bound * 2**30):
            return dict(status='resource_capacity_unavailable', host=host)
        try:
            cpus = cpu_allocation(2)
            gpu_suitability(gpu)
        except RuntimeError as error:
            return dict(status='resource_capacity_unavailable', host=host, reason=str(error))
        if (host, gpu) not in free_slots([host]):
            return dict(status='atomic_placement_unavailable', host=host)
        attempts = campaign / 'dispatch_attempts'; attempts.mkdir(exist_ok=True)
        attempt = attempts / f'{time.time_ns()}'; attempt.mkdir()
        caps = [f'{key}=2' for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
            'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')]
        argv = ['taskset', '--cpu-list', ','.join(map(str, cpus)), str(GATE), 'run', str(gpu),
            '--label', 'redox-conformer-geometry-batch', '--', 'env', *caps,
            'PYTHONPATH=' + str(ROOT / 'src'), sys.executable,
            str(ROOT / 'scripts/fleet/exclusive_python_worker.py'), '--guard-out', str(campaign / 'residency_guard'),
            str(ROOT / 'scripts/fleet/bounded_host_python_worker.py'), '--out', str(campaign / 'host_guard'),
            '--maximum-rss-gib', str(bound), '--minimum-free-gib', '32', manifest['worker'], 'run', '--campaign', str(campaign)]
        with (attempt / 'worker.log').open('x') as log:
            process = subprocess.Popen(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, start_new_session=True)
        time.sleep(.1)
        owner = dict(host=host, pid=process.pid,
            process_start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()[19])
        result = dict(status='gate_started_scientific_start_pending', owner=owner, GPU=gpu,
            CPU_affinity=cpus, argv=argv, started=time.time(), attempt=str(attempt),
            available_RAM_bytes=capacity['available_RAM_bytes'], remaining_other_promises_bytes=promised,
            plan_sha256=api.digest(campaign / 'plan.json'))
        atomic_write_json(attempt / 'launch.json', result)
        atomic_write_json(campaign / 'launch.json', result)
        return result


def owner_live(owner):
    module = api.load_module('geometry_dispatch_owner_check', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    return module.owner_alive(owner)


def preflight(batch, host):
    path = batch / (host + '_preflight.json')
    if path.exists() and time.time() - path.stat().st_mtime < 600:
        return json.loads(path.read_text())['ok']
    argv = [sys.executable, str(PREFLIGHT), host, '--python', sys.executable,
        '--module', 'pyscf', '--module', 'gpu4pyscf', '--module', 'cupy', '--module', 'rdkit',
        '--module', 'ase', '--module', 'scipy', '--path', str(batch / 'batch.json'),
        '--bin', str(ROOT.parent / 'bin/gpu_reserve'), '--min-available-ram-mib', '49152']
    result = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    row = json.loads(result.stdout)
    atomic_write_json(path, row)
    return result.returncode == 0 and row['ok']


def coordinate(batch, hosts, hours, resume, stage_root=None):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    manifest = json.loads((batch / 'batch.json').read_text())
    api.verify_sources(manifest['source_sha256'])
    sources = {str(p): api.digest(p) for p in (SCRIPT, batch / 'batch.json')}
    if stage_root:
        for p in (ROOT / 'scripts/pipeline/prepare_crest_refinement.py', stage_root / 'request.json'):
            sources[str(p)] = api.digest(p)
    stage = manifest.get('stage', 'geometry')
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    with (batch / 'dispatch.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status_path = batch / 'dispatch.json'
        if status_path.exists():
            prior = json.loads(status_path.read_text())
            if not resume or owner_live(prior['owner']):
                raise ValueError('Prior dispatcher is live or explicit resume is missing')
            history = batch / 'dispatch_history'; history.mkdir(exist_ok=True)
            (history / f'{time.time_ns()}.json').write_bytes(status_path.read_bytes())
        owner = dict(host=socket.gethostname().split('.')[0], pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
        record = dict(owner=owner, started=time.time(), deadline=time.time() + hours * 3600,
            hosts=hosts, source_sha256=sources, CPU_affinity=sorted(os.sched_getaffinity(0)))
        failures_by_host = {}; processed = set(); assignments = {}
        try:
            while time.time() < record['deadline']:
                api.verify_sources(sources); api.verify_sources(manifest['source_sha256'])
                if stage_root:
                    stager = api.load_module('continuing_source_conformer_staging', ROOT / 'scripts/pipeline/prepare_crest_refinement.py')
                    record['source_staging'] = stager.stage(stage_root)
                pending = []; rows = []
                for entry in manifest['entries']:
                    campaign = Path(entry['campaign']); key = entry['key']
                    if stage == 'electronic' and not (campaign / 'plan.json').exists():
                        preparer = api.load_module('dependency_checked_conformer_SP_preparation', Path(manifest['preparer']))
                        ready = preparer.prepare_ready(entry, manifest)
                        if ready != 'prepared':
                            rows.append(dict(key=key, campaign=str(campaign), status=ready)); continue
                    if entry['plan_sha256'] is not None and api.digest(campaign / 'plan.json') != entry['plan_sha256']:
                        raise ValueError('A declared campaign plan changed')
                    path = (campaign / 'geometry/status.json' if stage == 'geometry'
                            else campaign / 'energy/tasks' / key / 'status.json')
                    row = json.loads(path.read_text()) if path.exists() else {}
                    launch_path = campaign / 'launch.json'
                    handle = json.loads(launch_path.read_text()) if launch_path.exists() else assignments.get(key)
                    if 'finished' in row:
                        passed = row['status'] == ('stationary_geometry_checked_active_parent_and_native_minimum_pending'
                            if stage == 'geometry' else 'state_records_verified')
                        if passed:
                            api.verify_sources(row['output_sha256'])
                            if stage == 'electronic':
                                admitted, invalid = api.admitted_records(campaign, 'energy')
                                if invalid or len(admitted) != 1:
                                    raise ValueError('Canonical conformer energy receipt failed admission')
                                if not owner_live(row['owner']):
                                    normalizer = api.load_module('conformer_energy_normal_handoff', ROOT / 'scripts/pipeline/normalize_polymer_energy_handoff.py')
                                    normalizer.normalize(campaign)
                        state = ('checked_geometry' if stage == 'geometry' else 'checked_electronic') if passed else 'failed_preserved'
                        if key not in processed and handle and handle.get('owner'):
                            host = handle['owner']['host']
                            failures_by_host[host] = 0 if passed else failures_by_host.get(host, 0) + 1
                            processed.add(key)
                    elif handle and handle.get('owner'):
                        if handle['owner']['pid'] < 1:
                            live = None
                        else:
                            try:
                                live = owner_live(handle['owner'])
                            except (OSError, RuntimeError, subprocess.SubprocessError):
                                live = None
                        if live is not False:
                            state = 'running' if live else 'ownership_unverified_preserved'
                        elif (campaign / 'residency_guard').exists() or path.exists():
                            state = 'stopped_scientific_attempt_preserved'
                            if key not in processed:
                                host = handle['owner']['host']; failures_by_host[host] = failures_by_host.get(host, 0) + 1; processed.add(key)
                        elif len(list((campaign / 'dispatch_attempts').glob('*'))) >= 4:
                            state = 'placement_attempts_require_review'
                        else:
                            state = 'pending_gate_contention'; pending.append(entry)
                    else:
                        intent = campaign / 'dispatch_intent.json'
                        unresolved = intent.exists() and json.loads(intent.read_text()).get('outcome') != 'not_started'
                        if unresolved:
                            state = 'ownership_unverified_preserved'
                        else:
                            state = 'pending'; pending.append(entry)
                    rows.append(dict(key=key, campaign=str(campaign), status=state))
                record.update(status=stage + '_batch_in_progress', updated=time.time(), entries=rows,
                    checked_geometries=sum(r['status'] == 'checked_geometry' for r in rows),
                    checked_electronic=sum(r['status'] == 'checked_electronic' for r in rows),
                    running=sum(r['status'] == 'running' for r in rows),
                    preserved_failures=sum(r['status'] in ('failed_preserved', 'stopped_scientific_attempt_preserved') for r in rows),
                    consecutive_host_failures=failures_by_host, scientific_energies_or_G_supplied=False)
                atomic_write_json(status_path, record)
                if not pending and not any(r['status'] in ('running', 'ownership_unverified_preserved', 'upstream_pending') for r in rows):
                    record['status'] = stage + '_batch_terminal_physical_followup_separate'; return
                slots = free_slots(hosts)
                for host in sorted({h for h, _ in slots}):
                    if not pending or failures_by_host.get(host, 0) >= 3 or not preflight(batch, host):
                        continue
                    for _, gpu in [s for s in slots if s[0] == host]:
                        if not pending:
                            break
                        entry = pending[0]
                        command = ['env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1',
                            'PYTHONPATH=' + str(ROOT / 'src'), sys.executable, str(SCRIPT), 'launch-local',
                            '--batch', str(batch), '--campaign', entry['campaign'], '--gpu', str(gpu)]
                        # Save intent before RPC; an unknown response keeps this
                        # task preserved until its remote launch record resolves.
                        intent = Path(entry['campaign']) / 'dispatch_intent.json'
                        atomic_write_json(intent, dict(host=host, GPU=gpu, owner=owner, checked=time.time()))
                        try:
                            result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8',
                                host, shlex.join(command)], capture_output=True, text=True, timeout=45, check=True)
                            observation = json.loads(result.stdout)
                        except (OSError, ValueError, subprocess.SubprocessError) as error:
                            assignments[entry['key']] = dict(owner=dict(host=host, pid=-1, process_start='unknown'), error=repr(error))
                            pending.pop(0); continue
                        with (batch / 'dispatch_events.jsonl').open('a') as stream:
                            stream.write(json.dumps(dict(key=entry['key'], observation=observation, checked=time.time())) + '\n')
                        if observation['status'] == 'gate_started_scientific_start_pending':
                            assignments[entry['key']] = observation; pending.pop(0)
                        elif observation['status'] == 'resource_capacity_unavailable':
                            atomic_write_json(intent, dict(host=host, GPU=gpu, owner=owner,
                                checked=time.time(), outcome='not_started'))
                            break
                        else:
                            atomic_write_json(intent, dict(host=host, GPU=gpu, owner=owner,
                                checked=time.time(), outcome='not_started'))
                time.sleep(15)
            record['status'] = 'dispatch_window_ended_running_outputs_preserved'
        except BaseException as error:
            record.update(status='dispatch_stopped_running_workers_preserved', error=repr(error)); raise
        finally:
            record['finished'] = time.time()
            atomic_write_json(status_path, record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('coordinate', 'launch-local'))
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--campaign', type=Path)
    parser.add_argument('--gpu', type=int)
    parser.add_argument('--hosts', default='lambda1,lambda2,lambda11')
    parser.add_argument('--max-hours', type=float, default=48)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--stage-root', type=Path)
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48:
        parser.error('Use a finite window up to 48 hours')
    if any(not host or host.startswith('-') for host in args.hosts.split(',')):
        parser.error('Invalid host list')
    if args.action == 'launch-local':
        if args.campaign is None or args.gpu is None:
            parser.error('launch-local requires --campaign and --gpu')
        print(json.dumps(launch_local(args.batch.resolve(), args.campaign.resolve(), args.gpu)))
    else:
        coordinate(args.batch.resolve(), args.hosts.split(','), args.max_hours, args.resume,
            args.stage_root.resolve() if args.stage_root else None)


if __name__ == '__main__':
    main()
