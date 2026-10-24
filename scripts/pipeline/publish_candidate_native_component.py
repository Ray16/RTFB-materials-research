#!/usr/bin/env python
"""Independently review a completed component before its matched-cycle handoff.

Optional bounded observation starts no quantum jobs and never retries a failed
state. Publication preserves the prior registry and requires ordinary energy,
native, exact-model and independent nuclear/geometry admission.
"""
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import socket
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

REVIEWER = ROOT / 'scripts/validation/review_polymer_native_component.py'
COORDINATOR = ROOT / 'scripts/pipeline/coordinate_candidate_baseline.py'


def component_binding(review, energy, thermal):
    if (review['status'] != 'available_exact_polymer_native_components_independently_reconstructed'
            or review['checked_components'] != 1 or len(review['states']) != 1
            or not review['normal_admission_rechecked'] or review['QM_recomputed']):
        raise ValueError('Exactly one independently reconstructed native component is required')
    row = review['states'][0]
    if (row['key'] != thermal['task']['key'] or energy['task'] != thermal['task']
            or row['state_dir'] != thermal['state_dir']
            or row['occupied_sites'] != thermal['task']['occupied_sites']
            or not all(math.isfinite(v) for v in (row['G_eV'], thermal['G_eV'], row['E_eV'], energy['E_smd_eV']))
            or abs(row['G_eV'] - thermal['G_eV']) > 1e-8
            or abs(row['E_eV'] - energy['E_smd_eV']) > 1e-8):
        raise ValueError('Independent review and exact electronic/native component differ')
    return thermal['task']['model_key']


def publish(base, campaign, review_path):
    review = json.loads(review_path.read_text())
    api.verify_sources(review['source_sha256'])
    energy, invalid_e = api.admitted_records(campaign, 'energy')
    thermal, invalid_g = api.admitted_records(campaign, 'thermal')
    if invalid_e or invalid_g or len(energy) != 1 or len(thermal) != 1:
        raise ValueError('Exactly one normally admitted electronic/native component is required')
    key = component_binding(review, next(iter(energy.values())), next(iter(thermal.values())))
    declaration = json.loads((base / 'baseline_plan.json').read_text())
    api.verify_sources(declaration['source_sha256'])
    cases = [case for case in declaration['cases'] if case['key'] == key]
    if len(cases) != 1:
        raise ValueError('The native component is not in the declared matched cycle')
    coordinator = api.load_module('native_publication_exact_model_check', COORDINATOR)
    plan = api.read_plan(campaign)
    coordinator.same_model(api.read_plan(Path(cases[0]['campaign'])), plan)
    api.verify_sources(plan['producer_sha256'])
    api.verify_sources(plan['input_sha256'])
    if not plan.get('recovery_scope'):
        raise ValueError('A recoverable model scope is required by the cycle collector')
    registry_path = base / 'baseline_recoveries.json'
    with registry_path.with_suffix('.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        registry = json.loads(registry_path.read_text())
        history = base / 'baseline_observer_history'
        history.mkdir(exist_ok=True)
        archived = history / f'registry_before_checked_native_publication_{time.time_ns()}.json'
        archived.write_bytes(registry_path.read_bytes())
        registry.setdefault('overrides', {})[key] = dict(campaign=str(campaign),
            independent_review=str(review_path), independent_review_sha256=api.digest(review_path),
            reason='Completed exact native component independently reconstructed; ordinary physical gates and exact electronic parent retained')
        registry.get('external_continuations', {}).pop(key, None)
        registry.update(previous_registry=str(archived), updated=time.time())
        atomic_write_json(registry_path, registry)
    return dict(status='checked_native_component_published_to_matched_cycle', key=key,
        campaign=str(campaign), independent_review=str(review_path),
        registry=str(registry_path), registry_sha256=api.digest(registry_path),
        prior_registry=str(archived), conformer_convergence_established=False,
        synthesis_admission=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('base', 'campaign', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--wait-hours', type=float, default=0)
    args = parser.parse_args()
    if not math.isfinite(args.wait_hours) or not 0 <= args.wait_hours <= 48:
        parser.error('Use a bounded observation window of zero to 48 hours')
    base, campaign, out = args.base.resolve(), args.campaign.resolve(), args.out.resolve()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    out.mkdir(parents=True, exist_ok=False)
    sources = {str(p): api.digest(p) for p in (Path(__file__).resolve(), REVIEWER, COORDINATOR, campaign / 'plan.json')}
    owner = dict(host=socket.gethostname(), pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), started=time.time(),
        deadline=time.time() + args.wait_hours * 3600, source_sha256=sources,
        quantum_jobs_submitted=0, status='waiting_for_completed_native_component')
    try:
        while True:
            api.verify_sources(sources)
            path = campaign / 'streamed_native/status.json'
            row = json.loads(path.read_text()) if path.exists() else {}
            status = row.get('status', '')
            if status == 'streamed_native_case_passed_independent_review_pending' and 'finished' in row:
                identity = row['owner']
                if identity['host'].split('.')[0] != socket.gethostname().split('.')[0]:
                    raise ValueError('Observe on the producer host so owner termination is directly verified')
                proc = Path('/proc') / str(identity['pid']) / 'stat'
                fields = proc.read_text().rsplit(')', 1)[1].split() if proc.exists() else []
                if not fields or fields[19] != identity['process_start'] or fields[0] == 'Z':
                    break
            if 'failed' in status or ('finished' in row and status != 'streamed_native_case_passed_independent_review_pending'):
                raise ValueError('Native calculation ended without qualification; outputs and registry preserved')
            owner.update(updated=time.time(), observed_native_status=status)
            atomic_write_json(out / 'observer.json', owner)
            if time.time() >= owner['deadline']:
                owner['status'] = 'observation_window_ended_outputs_preserved'
                return
            time.sleep(min(15, max(.1, owner['deadline'] - time.time())))
        reviewer = api.load_module('independent_native_component_publication_review', REVIEWER)
        reviewer.review(campaign, out / 'component')
        result = publish(base, campaign, out / 'component/summary.json')
        atomic_write_json(out / 'publication.json', result)
        owner['status'] = result['status']
    except BaseException as error:
        owner.update(status='publication_stopped_registry_preserved_unless_publication_record_exists', error=repr(error))
        raise
    finally:
        owner['finished'] = time.time()
        atomic_write_json(out / 'observer.json', owner)


if __name__ == '__main__':
    main()
