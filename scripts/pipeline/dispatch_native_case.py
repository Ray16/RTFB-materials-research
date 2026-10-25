#!/usr/bin/env python
"""Wait for suitable free capacity and launch one exact native campaign once.

Only atomic-placement failures may be retried; existing scientific attempts and
unknown submission ownership are preserved. Bulk geometry workers continue.
"""
import argparse
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


def launch_local(campaign, gpu):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots, gpu_suitability, GATE
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    driver = Path(plan['execution_driver'])
    if plan['producer_sha256'].get(str(driver)) != api.digest(driver):
        raise ValueError('The exact execution driver is not source-bound')
    if (campaign / 'streamed_native').exists() or (campaign / 'residency_guard').exists():
        raise ValueError('Preserve the existing scientific attempt')
    host = socket.gethostname().split('.')[0]
    runtime_definition = json.loads(Path(plan['runtime_profile']).read_text()) if plan.get('runtime_profile') else {}
    runtime_python = runtime_definition.get('python', sys.executable)
    preflight = json.loads((campaign / 'native_dispatch' / (host + '_preflight.json')).read_text())
    if not preflight['ok'] or preflight['python'] != runtime_python:
        raise ValueError('The node/interpreter preflight must pass')
    with allocation_guard():
        capacity = host_capacity(2, minimum_ram_gib=128)
        helpers = api.load_module('native_dispatch_remaining_promises', ROOT / 'scripts/pipeline/coordinate_geometry_batch.py')
        promised = helpers.remaining_promises()
        if capacity['available_RAM_bytes'] - promised < 128 * 2**30:
            return dict(status='waiting_for_host_capacity')
        cpus = cpu_allocation(2)
        gpu_suitability(gpu)
        if (host, gpu) not in free_slots([host]):
            return dict(status='placement_race_no_scientific_start')
        attempt = campaign / 'native_dispatch' / f'attempt_{time.time_ns()}'
        attempt.mkdir()
        caps = [f'{k}=2' for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
            'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')]
        argv = ['taskset', '--cpu-list', ','.join(map(str, cpus)), str(GATE), 'run', str(gpu),
            '--label', 'redox-source-bound-native-case', '--', 'env', *caps, 'PYTHONPATH=' + str(ROOT / 'src'),
            runtime_python, str(ROOT / 'scripts/fleet/exclusive_python_worker.py'), '--guard-out', str(campaign / 'residency_guard'),
            str(ROOT / 'scripts/fleet/bounded_host_python_worker.py'), '--out', str(campaign / 'host_guard'),
            '--maximum-rss-gib', str(plan['maximum_host_RSS_GiB']), '--minimum-free-gib', '32',
            str(ROOT / 'scripts/fleet/checkpoint_scf_worker.py'), '--out', str(campaign / 'SCF_checkpoints'),
            str(driver), 'run', '--campaign', str(campaign)]
        if plan.get('runtime_profile'):
            profile_path = Path(plan['runtime_profile'])
            if plan['input_sha256'].get(str(profile_path)) != api.digest(profile_path):
                raise ValueError('The contraction runtime definition is not source-bound')
            profile = json.loads(profile_path.read_text())
            library_path = profile.get('library_path', profile['loader_directory'])
            if os.environ.get('LD_LIBRARY_PATH'):
                library_path += os.pathsep + os.environ['LD_LIBRARY_PATH']
            env_index = argv.index('env') + 1
            argv.insert(env_index, 'LD_LIBRARY_PATH=' + library_path)
            index = argv.index(str(driver))
            argv[index:index] = [str(ROOT / 'scripts/fleet/tensor_runtime_worker.py'),
                '--profile', str(profile_path), '--out', str(campaign / 'tensor_runtime')]
        with (attempt / 'worker.log').open('x') as log:
            process = subprocess.Popen(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, start_new_session=True)
        time.sleep(.1)
        result = dict(status='gate_started_scientific_start_pending', host=host, gate_pid=process.pid,
            process_start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
            GPU=gpu, CPU_affinity=cpus, argv=argv, attempt=str(attempt), started=time.time(),
            available_RAM_bytes=capacity['available_RAM_bytes'], remaining_other_promises_bytes=promised,
            plan_sha256=api.digest(campaign / 'plan.json'))
        atomic_write_json(attempt / 'launch.json', result)
        atomic_write_json(campaign / 'launch.json', result)
        return result


def wait(campaign, hosts, hours):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    out = campaign / 'native_dispatch'; out.mkdir(exist_ok=False)
    record = dict(status='waiting_for_suitable_free_GPU_and_host_capacity', started=time.time(),
        deadline=time.time() + hours * 3600, host=socket.gethostname(), pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), scientific_starts=0, placement_attempts=0,
        source_sha256={str(p): api.digest(p) for p in (SCRIPT, campaign / 'plan.json', ROOT / 'scripts/pipeline/coordinate_geometry_batch.py')})
    try:
        while time.time() < record['deadline']:
            api.verify_sources(record['source_sha256'])
            if (campaign / 'streamed_native').exists() or (campaign / 'residency_guard').exists():
                record['status'] = 'scientific_attempt_exists_preserve_no_duplicate'; return
            for host, gpu in free_slots(hosts):
                preflight_path = out / (host + '_preflight.json')
                if not preflight_path.exists() or time.time() - preflight_path.stat().st_mtime > 600:
                    helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
                    plan = api.read_plan(campaign)
                    definition = json.loads(Path(plan['runtime_profile']).read_text()) if plan.get('runtime_profile') else {}
                    command = [sys.executable, str(helper), host, '--python', definition.get('python', sys.executable),
                        '--module', 'pyscf', '--module', 'gpu4pyscf', '--module', 'cupy', '--module', 'ase',
                        '--module', 'scipy', '--path', str(campaign / 'plan.json'), '--min-available-ram-mib', '131072']
                    if definition.get('loader_directory'):
                        command += ['--library-path', definition.get('library_path', definition['loader_directory'])]
                    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
                    preflight = json.loads(result.stdout); atomic_write_json(preflight_path, preflight)
                if not json.loads(preflight_path.read_text())['ok']:
                    continue
                command = ['env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'PYTHONPATH=' + str(ROOT / 'src'),
                    sys.executable, str(SCRIPT), 'launch-local', '--campaign', str(campaign), '--gpu', str(gpu)]
                try:
                    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host, shlex.join(command)],
                        capture_output=True, text=True, timeout=45, check=True)
                    observation = json.loads(result.stdout)
                except (OSError, ValueError, subprocess.SubprocessError) as error:
                    record.update(status='submission_ownership_unverified_preserved', error=repr(error)); return
                if observation['status'] == 'gate_started_scientific_start_pending':
                    record['placement_attempts'] += 1
                    record['launch'] = observation
                    checker = api.load_module('native_dispatch_gate_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
                    gate_owner = dict(host=host, pid=observation['gate_pid'], process_start=observation['process_start'])
                    for _ in range(12):
                        if (campaign / 'streamed_native/status.json').exists():
                            record.update(status='native_scientific_stage_started', scientific_starts=1); return
                        if not checker.owner_alive(gate_owner):
                            break
                        time.sleep(5)
                    if (campaign / 'residency_guard').exists() or checker.owner_alive(gate_owner):
                        record['status'] = 'scientific_start_or_ownership_needs_review_preserved'; return
                    if record['placement_attempts'] >= 4:
                        record['status'] = 'bounded_placement_attempts_need_review'; return
            record['updated'] = time.time(); atomic_write_json(out / 'status.json', record)
            time.sleep(5)
        record['status'] = 'placement_window_ended_no_extra_scientific_attempt'
    except BaseException as error:
        record.update(status='native_dispatch_failed_existing_workers_preserved', error=repr(error)); raise
    finally:
        record['finished'] = time.time(); atomic_write_json(out / 'status.json', record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('wait', 'launch-local'))
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--gpu', type=int)
    parser.add_argument('--hosts', default='lambda1,lambda2,lambda11')
    parser.add_argument('--max-hours', type=float, default=48)
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48 or any(not h or h.startswith('-') for h in args.hosts.split(',')):
        parser.error('Use valid hosts and a finite window up to 48 hours')
    if args.action == 'launch-local':
        if args.gpu is None:
            parser.error('--gpu is required')
        print(json.dumps(launch_local(args.campaign.resolve(), args.gpu)))
    else:
        wait(args.campaign.resolve(), args.hosts.split(','), args.max_hours)


if __name__ == '__main__':
    main()
