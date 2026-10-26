#!/usr/bin/env python
"""Review a finished fixed-device UMA pair without starting or retrying work."""
import argparse
import json
import os
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.validation.uma_parity import compare_points
from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def independently_review(plan_path, pair, guard):
    plan = json.loads(plan_path.read_text())
    ref_path, nitro_path = pair / 'reference/status.json', pair / 'nitro/status.json'
    reference, candidate = json.loads(ref_path.read_text()), json.loads(nitro_path.read_text())
    if any(record.get('plan_sha256') != api.digest(plan_path) for record in (reference, candidate)):
        raise ValueError('Pair records do not match the frozen input panel')
    if (reference['mode'] != 'reference' or candidate['mode'] != 'nitro'
            or reference['host'] != candidate['host'] or reference['GPU'] != candidate['GPU']
            or reference['capability'] != candidate['capability'] or candidate['capability'][0] not in (8,9,10,11,12)):
        raise ValueError('Pair does not use the same supported device and declared modes')
    frames = sum(len(case['frames']) for case in plan['cases'])
    if reference.get('explicit_reference_mode_checks', 0) < frames:
        raise ValueError('Reference-mode checks are incomplete')
    if any(row.get('actual_mode') != 'reference' for row in reference['cases']):
        raise ValueError('A reference case did not use the declared FP32 mode')
    if candidate.get('diagnostic_only', False):
        raise ValueError('Diagnostic behavior overrides cannot qualify the normal backend')
    for row in candidate['cases']:
        execution = row.get('execution', {})
        if (row.get('actual_mode') != 'nitro' or execution.get('actual_mode') != 'nitro'
                or 'uma_nitro' not in execution.get('executor', '')):
            raise ValueError('A candidate case did not use the actual Nitro executor')
    owner_path, state_path = guard / 'owner.json', guard / 'state.json'
    owner, state = json.loads(owner_path.read_text()), json.loads(state_path.read_text())
    if (not owner.get('guard_stopped') or state.get('foreign_resident_PIDs')
            or not state.get('own_residency_observed') or str(owner.get('GPU')) != str(candidate['GPU'])
            or owner.get('host', '').split('.')[0] != candidate['host'].split('.')[0]
            or (guard / 'interruption.json').exists() or (guard / 'observation_error.json').exists()):
        raise ValueError('Exclusive GPU residency or completed guard is not established')
    api.verify_sources(plan['input_sha256'])
    api.verify_sources(plan['producer_sha256'])
    result = compare_points(plan, reference, candidate)
    result.update(status='independent_fixed_device_pair_passed' if result['all_energy_and_force_checks_pass']
                  else 'independent_fixed_device_pair_failed',
                  actual_nitro_and_reference_modes_checked=True, same_device_checked=True,
                  exclusive_residency_checked=True, production_optimizer_qualified=False,
                  source_sha256={str(p): api.digest(p) for p in (plan_path, ref_path, nitro_path, owner_path,
                      state_path, Path(__file__).resolve(), ROOT/'src/redox/validation/uma_parity.py',
                      ROOT/'src/redox/qm/uma_nitro.py', ROOT/'scripts/validation/qualify_uma_current.py')},
                  scientific_scope='Fixed-geometry energy/force parity only; not optimizer, solution DFT or experimental accuracy')
    return result


def observe(waiter, plan, out, hours):
    if not 0 < hours <= 48 or out.exists():
        raise ValueError('Use a new observer directory and a finite window up to 48 hours')
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    out.mkdir(parents=True)
    owner = json.loads((waiter / 'owner.json').read_text())
    if owner['host'].split('.')[0] != os.uname().nodename.split('.')[0]:
        raise ValueError('Observer must verify the waiter on its own host')
    proof = {str(p): api.digest(p) for p in (Path(__file__).resolve(), plan,
             ROOT/'src/redox/validation/uma_parity.py', ROOT/'src/redox/qm/uma_nitro.py',
             ROOT/'scripts/validation/qualify_uma_current.py')}
    record = dict(status='waiting_for_fixed_device_pair', started=time.time(), deadline=time.time()+hours*3600,
                  owner=dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
                    process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
                  waiter=owner, source_sha256=proof, scientific_jobs_submitted=0, production_optimizer_qualified=False)
    status_path = out / 'status.json'
    atomic_write_json(status_path, record)
    try:
        while time.time() < record['deadline']:
            api.verify_sources(proof)
            status = json.loads((waiter / 'status.json').read_text())
            if status.get('status') == 'owned_pair_ended_read_retained_outcome':
                outcome = status.get('outcome', {})
                if outcome.get('status') != 'matched_device_UMA_pair_finished_independent_review_pending':
                    record.update(status='pair_finished_without_reviewable_result', outcome=outcome)
                    return
                pair = Path(outcome['comparison']).parent
                if pair != Path(status['raw_pair_output']).resolve():
                    raise ValueError('Reported pair output differs from the waiter record')
                attempt = pair.name.removeprefix('attempt_')
                if not attempt.isdigit():
                    raise ValueError('Unexpected pair attempt identifier')
                review = independently_review(plan, pair, waiter / ('guard_' + attempt))
                atomic_write_json(out / 'independent_review.json', review)
                record.update(status=review['status'], review=str(out/'independent_review.json'),
                              review_sha256=api.digest(out/'independent_review.json'))
                return
            proc = Path('/proc') / str(owner['pid']) / 'stat'
            live = proc.exists() and proc.read_text().rsplit(')',1)[1].split()[19] == owner['process_start']
            if not live:
                record.update(status='waiter_not_live_no_reviewable_result', waiter_status=status.get('status'))
                return
            record.update(checked=time.time(), waiter_status=status.get('status'),
                          actual_scientific_starts=status.get('actual_scientific_starts'))
            atomic_write_json(status_path, record)
            time.sleep(15)
        record['status'] = 'observer_window_ended_without_reviewable_result'
    except BaseException as error:
        record.update(status='independent_review_failed_preserved', error=repr(error))
        raise
    finally:
        record['finished'] = time.time()
        atomic_write_json(status_path, record)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--waiter', type=Path, required=True)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--max-hours', type=float, default=48)
    args = p.parse_args()
    observe(args.waiter.resolve(), args.plan.resolve(), args.out.resolve(), args.max_hours)


if __name__ == '__main__':
    main()
