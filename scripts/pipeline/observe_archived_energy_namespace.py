#!/usr/bin/env python3
"""Independently review completed archived hybrid-energy namespace attempts."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.archived_energy_namespace import energy_mount
from redox.workflows.host_local_review import command_for_host
from redox.validation.namespace_sources import read_checked_review
from observe_archived_native_namespace import alive


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--specs', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--hours', type=float, default=24)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    paths = json.loads(a.specs.read_text())['specs']
    if len(paths) != len(set(paths)):
        raise ValueError('Distinct attempt identities required')
    a.out.mkdir(parents=True, exist_ok=False)
    proof = {path: sha256_file(path) for path in paths}
    for path in (Path(__file__).resolve(), ROOT/'src/redox/workflows/archived_energy_namespace.py',
                 ROOT/'src/redox/workflows/archived_native_namespace.py',
                 ROOT/'src/redox/workflows/host_local_review.py',
                 ROOT/'src/redox/validation/namespace_sources.py',
                 ROOT/'scripts/pipeline/observe_archived_native_namespace.py'):
        proof[str(path)] = sha256_file(path)
    record = dict(status='waiting_for_energy_and_closed_ownership', started=time.time(),
                  deadline=time.time()+a.hours*3600, source_sha256=proof,
                  owner=dict(host=os.uname().nodename, pid=os.getpid(),
                             process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
                  native_G_admitted=False)
    done = {}
    reviews_returned = set()
    while time.time() < record['deadline']:
        if any(sha256_file(path) != digest for path, digest in proof.items()):
            raise ValueError('Energy namespace changed')
        rows = []
        for path in paths:
            if path in done:
                rows.append(done[path])
                continue
            spec = json.loads(Path(path).read_text())
            case = Path(spec['attempt'])
            state_path = case/'experimental_energy/status.json'
            row = dict(spec=path, key=spec['key'], status='energy_stage_not_yet_written')
            if state_path.exists():
                state = json.loads(state_path.read_text())
                row['status'] = 'energy_stage_running_or_owner_unverified'
                if 'finished' in state:
                    if state['status'] != 'experimental_hybrid_energy_finished_independent_review_pending':
                        row.update(status='energy_failure_preserved', error=state.get('error'))
                        done[path] = row
                    else:
                        try:
                            live = alive(state['owner'])
                        except (OSError, ValueError, subprocess.SubprocessError):
                            live = None
                        guard = case/'residency_guard/owner.json'
                        if live is False and guard.exists() and json.loads(guard.read_text()).get('guard_stopped'):
                            prior = Path(spec['prior_root'])
                            profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
                            target = case/'independent_energy_review.json'
                            command = energy_mount(spec)+['env', 'PYTHONPATH='+str(prior/'src'), profile['python'],
                                str(prior/'scripts/validation/review_smooth_polymer_energy.py'),
                                '--campaign', spec['original_campaign'], '--out',
                                str(Path(spec['original_campaign'])/'independent_energy_review.json')]
                            command = command_for_host(command, state['owner']['host'])
                            if not target.exists() and path not in reviews_returned:
                                with (a.out/(spec['key']+'_review.log')).open('x') as log:
                                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
                                if result.returncode:
                                    row.update(status='energy_independent_review_failed_preserved', returncode=result.returncode)
                                    done[path] = row
                                    rows.append(row)
                                    continue
                                reviews_returned.add(path)
                            try:
                                report, binding = read_checked_review(path, 'energy')
                            except FileNotFoundError:
                                row['status'] = 'completed_energy_review_file_visibility_pending'
                            else:
                                row.update(status='energy_independently_reviewed', **binding,
                                           E_smd_eV=report['E_smd_eV'])
                                done[path] = row
            rows.append(row)
        record.update(updated=time.time(), entries=rows, terminal=len(done), expected=len(paths))
        atomic_write_json(a.out/'status.json', record)
        if len(done) == len(paths):
            record.update(status='energy_attempts_terminal_native_followup_separate', finished=time.time())
            atomic_write_json(a.out/'status.json', record)
            return
        time.sleep(20)
    record.update(status='energy_observation_window_ended_workers_untouched', finished=time.time())
    atomic_write_json(a.out/'status.json', record)


if __name__ == '__main__':
    main()
