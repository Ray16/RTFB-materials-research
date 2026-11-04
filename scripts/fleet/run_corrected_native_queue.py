#!/usr/bin/env python3
"""Drain checked corrected sampling energies into native thermochemistry.

Wait for the largest independently reconstructed pilot before wider admission.
Preserve every failed or uncertain attempt; reuse the maintained shared GPU,
CPU, RAM and storage dispatcher rather than implementing another launcher.
"""
import argparse
from collections import Counter
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.sampling_native_receipts import ReceiptCache
from redox.workflows.molecular_native_namespace import load
from redox.workflows.sampling_native_queue import collect, prepare_ready, roster, extend_peers
from redox.workflows.prelaunch_accounting import prelaunch_accounting,launcher_versions


def terminal_incomplete(prepared_count, entries, eligible, checked_keys):
    if prepared_count != 240 or eligible:
        return False
    active = {'running', 'guard_started', 'submitted', 'ownership_unverified', 'uncertain_submission_preserved'}
    return not any(row['status'] in active
        or (row['status'] == 'completed_review_pending' and row['key'] not in checked_keys)
        or (row['status'] == 'prelaunch_refusal_no_scientific_attempt'
            and row.get('prelaunch_failures',row['submission_receipts']) < 3)
        for row in entries)


def run(manifest_path, energy_stream, pilot_manifest, resource_review, archive,
        peers_path, output, interpreter, hours, resume=False, seed_reviews=()):
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation, GATE
    from redox.workflows.sampling_native import dispatch
    from redox.validation.execution_review import freeze_source
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    roster(manifest_path)
    output.mkdir(parents=True, exist_ok=True)
    owners = load(ROOT / 'scripts/pipeline/observe_archived_native_namespace.py', 'native_queue_owners')
    observer = load(ROOT / 'scripts/fleet/run_screening_energy_panel.py', 'native_queue_observer')
    shared = load(ROOT / 'scripts/fleet/run_corrected_sampling_queue.py', 'native_queue_shared_policy')
    policy = Path(os.environ.get('GPU_RESERVE_CONFIG', str(GATE.parent.parent / '.config/gpu_reserve/config.env')))
    paths = [manifest_path, pilot_manifest, resource_review, policy, Path(__file__),
             ROOT / 'scripts/pipeline/observe_archived_native_namespace.py',
             ROOT / 'scripts/fleet/run_screening_energy_panel.py',
             ROOT / 'scripts/fleet/run_corrected_sampling_queue.py',
             ROOT / 'src/redox/workflows/charged_native_followup.py',
             ROOT / 'src/redox/workflows/screening_control.py',
             ROOT / 'src/redox/workflows/prelaunch_accounting.py',
             ROOT / 'src/redox/workflows/submission_refusal.py',
             ROOT / 'src/redox/workflows/sampling_native.py',
             ROOT / 'src/redox/workflows/sampling_native_queue.py',
             ROOT / 'src/redox/validation/sampling_native_receipts.py',
             ROOT / 'src/redox/validation/sampling_receipts.py',
             ROOT / 'src/redox/validation/execution_review.py',
             ROOT / 'scripts/validation/review_corrected_sampling_native.py']
    hashes = {str(p): sha256_file(p) for p in paths}
    snapshots = {str(freeze_source(output, p)): hashes[str(p)] for p in paths}
    with (output / 'controller.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status_path = output / 'status.json'
        previous = {}
        if status_path.exists():
            previous = json.loads(status_path.read_text())
            if not resume or owners.alive(previous['owner']) is not False:
                raise ValueError('Explicit resume and proven stopped controller required')
            atomic_write_json(output / 'history' / (str(time.time_ns()) + '.json'), previous)
        prepared_path = output / 'prepared_current.json'
        if prepared_path.exists():
            prepared = json.loads(prepared_path.read_text())['states']
        else:
            prepared = json.loads(pilot_manifest.read_text())['states']
            atomic_write_json(prepared_path, dict(declared_sampling_inputs=240, states=prepared))
        reviews = output / 'reviews'
        reviews.mkdir(exist_ok=True)
        by_key = {r['key']: r for r in prepared}
        cache = ReceiptCache()
        for path in seed_reviews:
            report = json.loads(path.read_text())
            row = by_key[report['key']]
            cache.read(path, Path(row['namespace']), row['key'], roster(manifest_path)['experimental_method_id'])
            target = reviews / (row['key'] + '.json')
            if target.exists() and target.read_bytes() != path.read_bytes():
                raise ValueError('Existing native review cannot be replaced')
            if not target.exists():
                target.write_bytes(path.read_bytes())
        owner = dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
                     process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
        state = dict(status='corrected_native_queue_active', owner=owner,
                     CPU_affinity=sorted(os.sched_getaffinity(0)), started=time.time(), deadline=time.time() + hours * 3600,
                     declared_inputs=240, source_sha256=snapshots, occupied_ensemble_cycles_qualified=False,
                     accounted_terminal=previous.get('accounted_terminal', []),
                     host_failures=previous.get('host_failures', {}), blocked_hosts=previous.get('blocked_hosts', []),
                     storage_peers=previous.get('storage_peers', {}),
                     prelaunch_retry_observed=previous.get('prelaunch_retry_observed', {}))
        stop = threading.Event()
        old_signal = signal.signal(signal.SIGTERM, lambda *args: stop.set())
        try:
            while not stop.is_set() and time.time() < state['deadline']:
                if any(sha256_file(p) != h for p, h in hashes.items()):
                    state['status'] = 'native_queue_source_changed_workers_preserved'
                    break
                peers = extend_peers(peers_path, prepared)
                peers, state['storage_peers'] = shared.storage_peers(peers_path, state['storage_peers'])
                result = collect(manifest_path, prepared, reviews, cache)
                atomic_write_json(output / 'native_current.json', result)
                state.update(checked=time.time(), checked_native_inputs=result['checked_native_inputs'])
                if result['independent_review_errors']:
                    state.update(status='native_review_errors_workers_preserved', review_errors=result['independent_review_errors'])
                    break
                checked = {r['key']: r for r in result['states'] if r['native_G_qualified']}
                if len(checked) == 240:
                    state['status'] = 'all_240_native_inputs_checked_ensemble_analysis_pending'
                    break
                largest = checked.get('all_same__occ_2_3__conf001')
                state['largest_native_pilot_checked'] = largest is not None
                entries, eligible = [], []
                terminal = []
                for row in prepared:
                    claim = Path(row['namespace']).parent / 'dispatch_claim'
                    condition, data = observer.inspect(row, claim, owners.alive,
                        finished_status='corrected_sampling_native_finished_independent_review_pending')
                    tries = len(list(claim.glob('try_*.json')))
                    key = row['key']
                    failures,accounting,accounting_error=tries,None,None
                    if condition == 'prelaunch_refusal_no_scientific_attempt':
                        try:
                            accounting=prelaunch_accounting(Path(row['namespace']),owners.alive,launcher_versions(output))
                            failures=accounting['prelaunch_failures']
                        except (OSError,ValueError,subprocess.SubprocessError) as error:
                            condition='ownership_unverified';accounting_error=repr(error);data=dict(error=accounting_error)
                    if condition == 'prelaunch_refusal_no_scientific_attempt':
                        state['prelaunch_retry_observed'].setdefault(key, time.time())
                    if shared.ready(condition, failures, state['prelaunch_retry_observed'].get(key, 0), 60, time.time()):
                        # Initial pilots were already explicitly authorized; all
                        # other rows can only exist after prepare's pilot gate.
                        eligible.append(row)
                    if condition in ('completed_review_pending', 'scientific_failure_preserved',
                                     'stopped_attempt_preserved', 'guard_start_failed_preserved') and key not in state['accounted_terminal']:
                        terminal.append((data.get('finished', time.time()), key, condition, data.get('owner', data).get('host')))
                    entries.append(dict(key=key,status=condition,submission_receipts=tries,prelaunch_failures=failures,
                                        resource_refusals=accounting['resource_refusals'] if accounting else None,
                                        prelaunch_accounting_error=accounting_error))
                for _, key, condition, host in sorted(terminal):
                    state['accounted_terminal'].append(key)
                    if host:
                        state['host_failures'][host] = 0 if condition == 'completed_review_pending' else state['host_failures'].get(host, 0) + 1
                        if state['host_failures'][host] >= 3 and host not in state['blocked_hosts']:
                            state['blocked_hosts'].append(host)
                state.update(prepared_inputs=len(prepared), counts=dict(Counter(r['status'] for r in entries)), entries=entries)
                state['phase'] = 'reviewing_and_dispatching_prepared_inputs'
                atomic_write_json(status_path, state)
                if eligible and not stop.is_set():
                    if any(sha256_file(p) != h for p, h in hashes.items()):
                        state['status'] = 'native_queue_source_changed_workers_preserved'
                        break
                    # Preserve all rows for storage accounting; only send the
                    # eligible subset to the existing scientifically gated API.
                    manifest = output / 'dispatch_manifests' / (str(time.time_ns()) + '.json')
                    atomic_write_json(manifest, dict(declared_sampling_inputs=240, states=eligible))
                    launched = dispatch(manifest, interpreter, peers, state['blocked_hosts'])
                    state['last_dispatch'] = dict(checked=time.time(), launches=len(launched['launches']), refusals=launched['refusals'])
                    for row in launched['launches']:
                        state['prelaunch_retry_observed'][row['key']] = time.time()
                largest_state = next((r['status'] for r in entries if r['key'] == 'all_same__occ_2_3__conf001'), None)
                if largest is None and largest_state in ('scientific_failure_preserved', 'stopped_attempt_preserved', 'guard_start_failed_preserved'):
                    state['status'] = 'largest_native_pilot_terminal_followups_held'
                    break
                if terminal_incomplete(len(prepared), entries, eligible, set(checked)):
                    state['status'] = 'native_queue_terminal_incomplete_recovery_required'
                    break
                # Start already prepared work before spending time validating
                # and copying the next growing group of checked energy inputs.
                if largest is not None and not stop.is_set():
                    state['phase'] = 'preparing_more_checked_inputs_after_dispatch'
                    atomic_write_json(status_path, state)
                    if any(sha256_file(p) != h for p, h in hashes.items()):
                        state['status'] = 'native_queue_source_changed_workers_preserved'
                        break
                    stream = json.loads(energy_stream.read_text())
                    prepared = prepare_ready(manifest_path, stream, prepared, Path(largest['review']),
                                             resource_review, archive, output)
                    peers = extend_peers(peers_path, prepared)
                atomic_write_json(status_path, state)
                stop.wait(30)
            else:
                state['status'] = 'native_queue_stop_requested_workers_preserved' if stop.is_set() else 'native_queue_window_finished_workers_preserved'
        except BaseException as error:
            state.update(status='native_queue_failed_workers_preserved', error=repr(error))
            raise
        finally:
            signal.signal(signal.SIGTERM, old_signal)
            state['finished'] = time.time()
            atomic_write_json(status_path, state)
    return state


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'energy-stream', 'pilot-manifest', 'resource-review', 'archive', 'storage-peers', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--python', required=True)
    p.add_argument('--seed-review', type=Path, action='append', default=[])
    p.add_argument('--hours', type=float, default=24)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    if not 0 < a.hours <= 72:
        p.error('A positive bounded window up to 72 hours is required')
    result = run(a.manifest.resolve(), a.energy_stream.resolve(), a.pilot_manifest.resolve(),
                 a.resource_review.resolve(), a.archive.resolve(), a.storage_peers.resolve(),
                 a.out.resolve(), a.python, a.hours, a.resume, [p.resolve() for p in a.seed_review])
    print(json.dumps(dict(status=result['status'], checked_native_inputs=result.get('checked_native_inputs', 0))))
