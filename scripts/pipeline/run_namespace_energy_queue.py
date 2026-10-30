#!/usr/bin/env python3
"""Run prepared electronic namespace states as guarded GPU capacity becomes free.

The queue launches the maintained per-state adapter. Independent physical review
is a separate collector; a written energy never substitutes for native G.
Uncertain submissions and stopped scientific attempts remain preserved.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def module(path, name):
    specification = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(result)
    return result


def inspect_state(spec, alive):
    case = Path(spec['attempt'])
    path = case/'experimental_energy/status.json'
    if path.exists():
        state = json.loads(path.read_text())
        if 'finished' in state:
            return ('finished_review_pending' if state['status'] ==
                    'experimental_hybrid_energy_finished_independent_review_pending' else 'failed_preserved'), state
        try:
            live = alive(state['owner'])
        except (OSError, RuntimeError, ValueError, subprocess.SubprocessError):
            live = None
        return ('running' if live else 'stopped_attempt_preserved' if live is False
                else 'ownership_unverified_preserved'), state
    guard = case/'residency_guard/owner.json'
    if guard.exists():
        return 'guard_present_scientific_start_pending', json.loads(guard.read_text())
    claim = case/'dispatch_claim'
    if claim.exists():
        owner = claim/'owner.json'
        return 'submission_claim_preserved', json.loads(owner.read_text()) if owner.exists() else {}
    return 'pending', {}


def capacity(host, gpu, interpreter):
    code = ('import json,os,psutil,subprocess; '
            'from redox.workflows.archived_native_namespace import remaining_host_promises; '
            'm=psutil.virtual_memory(); '
            'g=subprocess.run(["nvidia-smi","--id",'+repr(str(gpu))+','
            '"--query-gpu=memory.total","--format=csv,noheader,nounits"],'
            'capture_output=True,text=True,timeout=20,check=True); '
            'print(json.dumps(dict(available=m.available,promised=remaining_host_promises(),'
            'GPU_MiB=float(g.stdout.strip()),load=os.getloadavg()[0],cpus=os.cpu_count())))')
    argv = ['env', 'PYTHONPATH='+str(ROOT/'src'), 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1',
            interpreter, '-c', code]
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host, shlex.join(argv)],
                            capture_output=True, text=True, timeout=35, check=True)
    return json.loads(result.stdout)


def fits(spec, resources):
    return (resources['GPU_MiB'] >= 31*1024 and resources['load']+2 <= resources['cpus']
            and resources['available']-resources['promised'] >=
            (spec['maximum_host_RSS_GiB']+32)*2**30)


def discover_slots(hosts):
    """A failed node observation contributes no free GPUs and preserves peers."""
    from redox.workflows.resources import GATE

    def observe(host):
        try:
            result = subprocess.run([str(GATE), 'list', '--hosts', host],
                                    capture_output=True, text=True, timeout=15, check=True)
            rows = []
            for line in result.stdout.splitlines():
                fields = line.split()
                if len(fields) != 2 or fields[0] != host or not fields[1].isdigit():
                    raise ValueError('Unexpected reservation-gate availability response')
                rows.append((host, int(fields[1])))
            return rows, None
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            return [], dict(host=host, error=repr(error), availability='unverified_no_submission')

    with ThreadPoolExecutor(max_workers=max(1, len(hosts))) as pool:
        observations = list(pool.map(observe, hosts))
    return [slot for rows, _ in observations for slot in rows], [error for _, error in observations if error]


def release_gate_contention(spec, output, alive):
    """Release only this queue's proven gate refusal, never an uncertain start."""
    from redox.validation.namespace_sources import NamespaceSources
    case = Path(spec['attempt'])
    claim = case/'dispatch_claim'
    if not (claim/'owner.json').exists() or not (claim/'launch.json').exists():
        return False
    owner = json.loads((claim/'owner.json').read_text())
    if owner.get('queue') != str(output):
        return False
    receipts = json.loads((claim/'launch.json').read_text())['receipts']
    if len(receipts) != 1 or receipts[0]['status'] != 'controller_spawned_scientific_start_unverified':
        return False
    receipt = receipts[0]
    placement = Path(receipt['argv'][receipt['argv'].index('--out')+1])
    if not (placement/'completion.json').exists():
        return False
    result = json.loads((placement/'completion.json').read_text())
    refusal = f'gpu_reserve: gpu{result["GPU"]} BUSY/reserved -> not running'
    if (result.get('returncode') in (None, 0) or refusal not in (placement/'command.log').read_text().splitlines()
            or (case/'residency_guard').exists() or (case/'experimental_energy').exists()
            or NamespaceSources(spec).path(spec['original_raw_root']).exists()):
        return False
    owners = [receipt['owner'], dict(host=result['host'], pid=result['pid'], process_start=result['process_start'])]
    try:
        if any(alive(value) for value in owners):
            return False
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError):
        return False
    history = output/'gate_contention_history'/spec['key']/str(time.time_ns())
    shutil.copytree(claim, history)
    atomic_write_json(history/'resolution.json', dict(status='gate_refused_no_scientific_attempt',
                      stopped_owners=owners, placement=str(placement), checked=time.time(),
                      no_guard_or_scientific_outputs=True, inference_attempt_consumed=False))
    for child in claim.iterdir():
        child.unlink()
    claim.rmdir()
    return True


def coordinate(manifest_path, output, interpreter, hours, resume):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    from redox.validation.namespace_sources import read_checked_review
    launcher = module(ROOT/'scripts/pipeline/launch_archived_native_followup.py', 'namespace_explicit_launcher')
    owner_checks = module(ROOT/'scripts/pipeline/observe_archived_native_namespace.py', 'namespace_owner_checks')
    config = json.loads(manifest_path.read_text())
    paths = [Path(p) for p in config['specs']]
    specs = [json.loads(p.read_text()) for p in paths]
    if (len(paths) != len(set(paths)) or len({s['key'] for s in specs}) != len(specs)
            or any(s['native_G_computed'] for s in specs)):
        raise ValueError('Distinct prepared electronic namespace states required')
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    output.mkdir(parents=True, exist_ok=True)
    with (output/'controller.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = output/'status.json'
        if status.exists():
            previous = json.loads(status.read_text())
            if not resume or owner_checks.alive(previous['owner']):
                raise ValueError('Existing controller must be proven stopped and explicitly resumed')
            history = output/'history'
            history.mkdir(exist_ok=True)
            (history/(str(time.time_ns())+'.json')).write_bytes(status.read_bytes())
        producers = [Path(__file__).resolve(), manifest_path, ROOT/'scripts/pipeline/launch_archived_native_followup.py',
                     ROOT/'src/redox/workflows/derived_energy_namespace.py',
                     ROOT/'scripts/pipeline/run_measured_energy_followup.py',
                     ROOT/'src/redox/workflows/archived_energy_namespace.py']
        hashes = {str(path): sha256_file(path) for path in [*producers, *paths]}
        owner = dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
                     process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
        record = dict(status='prepared_energy_queue_running', owner=owner, started=time.time(),
                      deadline=time.time()+hours*3600, source_sha256=hashes, native_G_computed=False,
                      declared_upstream_states=config.get('declared_upstream_states'), prepared_states=len(paths))
        checked, stopped_hosts, accounted = {}, set(), set()
        failures = {}
        placement_index = len(list((output/'placements').glob('*'))) if (output/'placements').exists() else 0
        try:
            while time.time() < record['deadline']:
                if any(sha256_file(path) != digest for path, digest in hashes.items()):
                    raise ValueError('Queue source, manifest or prepared namespace changed')
                pending, rows, unsettled_hosts, newly_terminal = [], [], set(), []
                for path, spec in zip(paths, specs):
                    key = spec['key']
                    if key in checked:
                        rows.append(checked[key])
                        continue
                    release_gate_contention(spec, output, owner_checks.alive)
                    state, raw = inspect_state(spec, owner_checks.alive)
                    row = dict(key=key, spec=str(path), status=state)
                    if state == 'finished_review_pending' and (Path(spec['attempt'])/'independent_energy_review.json').exists():
                        report, binding = read_checked_review(path, 'energy')
                        row.update(status='electronic_independently_checked', E_smd_eV=report['E_smd_eV'], binding=binding)
                        checked[key] = row
                        claim = Path(spec['attempt'])/'dispatch_claim'
                        if (claim/'owner.json').exists() and json.loads((claim/'owner.json').read_text())['queue'] == str(output):
                            saved = output/'closed_claims'/key
                            saved.mkdir(parents=True, exist_ok=True)
                            for child in claim.iterdir():
                                target = saved/child.name
                                if target.exists() and sha256_file(target) != sha256_file(child):
                                    raise ValueError('Preserved claim history differs')
                                if not target.exists():
                                    shutil.copyfile(child, target)
                            for child in claim.iterdir():
                                child.unlink()
                            claim.rmdir()
                    elif state == 'pending':
                        pending.append((path, spec))
                    elif state in ('guard_present_scientific_start_pending', 'submission_claim_preserved'):
                        if 'host' in raw and not (Path(spec['attempt'])/'host_guard/status.json').exists():
                            unsettled_hosts.add(raw['host'].split('.')[0])
                    if key not in accounted and state in ('finished_review_pending', 'failed_preserved'):
                        newly_terminal.append((raw['finished'], key, state, raw['owner']['host'].split('.')[0]))
                    rows.append(row)
                for _, key, state, host in sorted(newly_terminal):
                    failures[host] = 0 if state == 'finished_review_pending' else failures.get(host, 0)+1
                    if failures[host] >= 3:
                        stopped_hosts.add(host)
                    accounted.add(key)
                record.update(updated=time.time(), entries=rows, independently_checked=len(checked),
                              pending=len(pending), stopped_hosts=sorted(stopped_hosts),
                              consecutive_scientific_host_failures=failures)
                atomic_write_json(status, record)
                if not pending and all(r['status'] in ('electronic_independently_checked', 'failed_preserved',
                                                       'stopped_attempt_preserved') for r in rows):
                    record['status'] = 'energy_queue_terminal_native_G_separate'
                    return
                slots, discovery_errors = discover_slots(config['hosts'])
                record['availability_observation_errors'] = discovery_errors
                atomic_write_json(status, record)
                for host, gpu in slots:
                    if not pending:
                        break
                    if host in stopped_hosts or host in unsettled_hosts:
                        continue
                    try:
                        resources = capacity(host, gpu, interpreter)
                    except (OSError, ValueError, subprocess.SubprocessError):
                        stopped_hosts.add(host)
                        continue
                    selected = next(((path, spec) for path, spec in pending if fits(spec, resources)), None)
                    if selected is None:
                        continue
                    path, spec = selected
                    claim = Path(spec['attempt'])/'dispatch_claim'
                    try:
                        claim.mkdir()
                    except FileExistsError:
                        continue
                    atomic_write_json(claim/'owner.json', dict(owner=owner, queue=str(output),
                                      namespace=str(path), namespace_sha256=hashes[str(path)],
                                      host=host, GPU=gpu, resource_preflight=resources, claimed=time.time()))
                    placement_index += 1
                    destination = output/'placements'/f'{placement_index:04d}'
                    destination.mkdir(parents=True)
                    entry = dict(key=spec['key'], spec=str(path), host=host, gpu=gpu,
                                 namespace_bin_directory=str(ROOT.parent/'redox_polymer_continuations/runtime/bin'))
                    try:
                        receipts = launcher.launch_host(host, [entry], destination, interpreter, 'measured-energy')
                        atomic_write_json(claim/'launch.json', dict(receipts=receipts))
                    except (OSError, ValueError, subprocess.SubprocessError) as error:
                        atomic_write_json(claim/'uncertain_submission.json', dict(error=repr(error), preserve=True))
                        stopped_hosts.add(host)
                    pending.remove(selected)
                    # The next task on this host must see the bounded worker's
                    # actual reservation promise. Until then, reserve the host
                    # for this unresolved handoff rather than guessing from RSS.
                    unsettled_hosts.add(host)
                time.sleep(20)
            record['status'] = 'energy_queue_window_ended_existing_workers_untouched'
        except BaseException as error:
            record.update(status='energy_queue_stopped_existing_workers_preserved', error=repr(error))
            raise
        finally:
            record['finished'] = time.time()
            atomic_write_json(status, record)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--python', required=True)
    p.add_argument('--hours', type=float, default=24)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args()
    if not 0 < args.hours <= 48:
        p.error('Use a finite observation window up to 48 hours')
    coordinate(args.manifest.resolve(), args.out.resolve(), args.python, args.hours, args.resume)
