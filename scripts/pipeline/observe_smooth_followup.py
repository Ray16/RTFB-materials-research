#!/usr/bin/env python
"""Review completed alternate-model stages without launching quantum work."""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def observe(batches, energies, out, hours):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    geometry_script = ROOT / 'scripts/validation/review_smooth_geometry.py'
    energy_script = ROOT / 'scripts/validation/review_smooth_polymer_energy.py'
    geometry = api.load_module('smooth_followup_geometry_review', geometry_script)
    electronic = api.load_module('smooth_followup_energy_review', energy_script)
    owners = api.load_module('smooth_followup_ownership', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    sources = {str(p):api.digest(p) for p in (Path(__file__).resolve(), geometry_script, energy_script)}
    tasks = []
    for batch in batches:
        path = batch / 'batch.json'; manifest = json.loads(path.read_text())
        if manifest['stage'] != 'geometry' or manifest['primary_candidate_admission']:
            raise ValueError('Only declared alternate-model geometry batches are supported')
        api.verify_sources(manifest['source_sha256']); sources[str(path)] = api.digest(path)
        for entry in manifest['entries']:
            case = Path(entry['campaign'])
            if api.digest(case / 'plan.json') != entry['plan_sha256']:
                raise ValueError('A declared geometry plan changed')
            tasks.append(dict(key=batch.name+'__'+entry['key'], campaign=case, kind='geometry'))
    for case in energies:
        plan = api.read_plan(case)
        if plan['primary_candidate_admission'] or plan['native_G_computed']:
            raise ValueError('Only separate electronic controls without native G are supported')
        sources[str(case / 'plan.json')] = api.digest(case / 'plan.json')
        tasks.append(dict(key=case.parent.name+'__'+case.name, campaign=case, kind='energy'))
    if len({t['key'] for t in tasks}) != len(tasks):
        raise ValueError('Each scientific stage must have one observer identity')
    out.mkdir(parents=True, exist_ok=False)
    record = dict(status='CPU_independent_stage_reviews_in_progress', started=time.time(),
        deadline=time.time()+hours*3600, owner=dict(host=socket.gethostname().split('.')[0],
            pid=os.getpid(), process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        CPU_affinity=sorted(os.sched_getaffinity(0)), source_sha256=sources,
        quantum_jobs_submitted=0, native_G_qualification_inferred=False, reviewed={})
    try:
        while time.time() < record['deadline']:
            api.verify_sources(sources)
            pending = 0
            for task in tasks:
                key, case, kind = task['key'], task['campaign'], task['kind']
                if key in record['reviewed']: continue
                path = case / ('geometry/status.json' if kind == 'geometry' else 'experimental_energy/status.json')
                if not path.exists(): pending += 1; continue
                row = json.loads(path.read_text())
                if 'finished' not in row: pending += 1; continue
                expected = ('experimental_smooth_geometry_checked_energy_and_native_QC_pending'
                    if kind == 'geometry' else 'experimental_hybrid_energy_finished_independent_review_pending')
                if row['status'] != expected:
                    record['reviewed'][key] = dict(status='source_failure_preserved', source_status=row['status'],
                        campaign=str(case), native_G_qualified=False)
                    continue
                try:
                    live = owners.owner_alive(row['owner'])
                except (OSError, RuntimeError, subprocess.SubprocessError):
                    pending += 1; continue
                guard = case / 'residency_guard/owner.json'
                if live or not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):
                    pending += 1; continue
                destination = out / 'checked' / (key+'.json')
                try:
                    if kind == 'geometry': geometry.review(case, destination)
                    else: electronic.review(case, destination)
                except Exception as error:
                    record['reviewed'][key] = dict(status='independent_review_failed_preserved',
                        error=repr(error), campaign=str(case), native_G_qualified=False)
                else:
                    record['reviewed'][key] = dict(status='independent_stage_review_passed',
                        campaign=str(case), report=str(destination), report_sha256=api.digest(destination),
                        native_G_qualified=False)
                atomic_write_json(out / 'status.json', record)
            record.update(updated=time.time(), declared_stages=len(tasks), pending_stages=pending,
                independent_reviews_passed=sum(r['status']=='independent_stage_review_passed' for r in record['reviewed'].values()))
            atomic_write_json(out / 'status.json', record)
            if not pending:
                record['status'] = 'all_declared_stages_reviewed_or_source_failures_preserved'; return
            time.sleep(15)
        record['status'] = 'observer_window_ended_pending_stages_preserved'
    except BaseException as error:
        record.update(status='CPU_observer_failed_scientific_workers_untouched', error=repr(error)); raise
    finally:
        record['finished'] = time.time(); atomic_write_json(out / 'status.json', record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry-batch', type=Path, action='append', default=[])
    parser.add_argument('--energy-case', type=Path, action='append', default=[])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-hours', type=float, default=24)
    args = parser.parse_args()
    if not args.geometry_batch and not args.energy_case: parser.error('Declare at least one exact stage')
    if not 0 < args.max_hours <= 48: parser.error('Use a finite observation window up to 48 hours')
    observe([p.resolve() for p in args.geometry_batch], [p.resolve() for p in args.energy_case],
            args.out.resolve(), args.max_hours)


if __name__ == '__main__': main()
