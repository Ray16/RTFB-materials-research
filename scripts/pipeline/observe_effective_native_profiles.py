#!/usr/bin/env python3
"""Collect effective-budget native pilots with physical and resource gates."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import read_checked_review
from redox.workflows.effective_native_namespace import mount, review, verify
from redox.workflows.host_local_review import command_for_host


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profiles', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--hours', type=float, default=24)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    a.out.mkdir(parents=True, exist_ok=False)
    paths = [Path(path) for path in json.loads(a.profiles.read_text())['profiles']]
    if len(paths) != len(set(paths)):
        raise ValueError('Distinct execution profiles required')
    proof = {str(path): sha256_file(path) for path in [a.profiles.resolve(), *paths, Path(__file__).resolve(),
        ROOT/'src/redox/workflows/effective_native_namespace.py',
        ROOT/'scripts/validation/review_effective_native.py', ROOT/'src/redox/validation/namespace_sources.py']}
    copy_root = paths[0].parent/'qc_sources'/str(time.time_ns())
    copy_root.mkdir(parents=True)
    copied = copy_root/'review_effective_native.py'
    shutil.copyfile(ROOT/'scripts/validation/review_effective_native.py', copied)
    proof[str(copied)] = sha256_file(copied)
    specification = importlib.util.spec_from_file_location('effective_profile_owner_checks',
        ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(owners)
    record = dict(status='waiting_for_effective_native_physical_and_resource_checks', started=time.time(),
                  owner=dict(host=os.uname().nodename, pid=os.getpid(),
                    process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
                  deadline=time.time()+a.hours*3600, source_sha256=proof, profiles_qualified=False)
    record['review_source_copy'] = str(copied)
    done = {}
    while time.time() < record['deadline']:
        if any(sha256_file(path) != digest for path, digest in proof.items()):
            raise ValueError('Effective native observer/profile source changed')
        rows = []
        for path in paths:
            if str(path) in done:
                rows.append(done[str(path)])
                continue
            profile, spec = verify(path)
            case = Path(spec['attempt'])
            stage = case/'experimental_native/status.json'
            row = dict(profile=str(path), key=spec['key'], status='scientific_stage_not_yet_written')
            if stage.exists():
                state = json.loads(stage.read_text())
                row['status'] = 'native_running_or_ownership_unverified'
                if 'finished' in state:
                    if state['status'] != 'experimental_native_finished_independent_review_pending':
                        row.update(status='failed_native_attempt_preserved', error=state.get('error'))
                        done[str(path)] = row
                    else:
                        try:
                            live = owners.alive(state['owner'])
                        except (OSError, ValueError, subprocess.SubprocessError):
                            live = None
                        guard = case/'residency_guard/owner.json'
                        if live is False and guard.exists() and json.loads(guard.read_text()).get('guard_stopped'):
                            report_path = case/'independent_native_review.json'
                            if not report_path.exists():
                                runtime = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
                                command = mount(profile, spec)+['env', 'PYTHONPATH='+spec['prior_root']+'/src',
                                    runtime['python'], str(copied), '--campaign', spec['original_campaign'],
                                    '--profile', str(path), '--out', str(Path(spec['original_campaign'])/'independent_native_review.json')]
                                command = command_for_host(command, state['owner']['host'])
                                with (a.out/(spec['key']+'_physical_review.log')).open('x') as log:
                                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
                                if result.returncode:
                                    row.update(status='effective_native_independent_review_failed_preserved',
                                               returncode=result.returncode)
                                    done[str(path)] = row
                                    rows.append(row)
                                    continue
                            native, binding = read_checked_review(Path(profile['native_namespace']), 'native')
                            resource_path = case/'effective_resource_review.json'
                            if not resource_path.exists():
                                review(path, resource_path)
                            row.update(status='native_and_effective_resource_checks_passed', binding=binding,
                                       resource_review=str(resource_path), resource_review_sha256=sha256_file(resource_path),
                                       G_eV=native['G_by_cutoff_eV']['50'])
                            done[str(path)] = row
            rows.append(row)
        record.update(updated=time.time(), entries=rows, terminal=len(done), expected=len(paths))
        atomic_write_json(a.out/'status.json', record)
        if len(done) == len(paths):
            record.update(status='effective_native_pilots_terminal_outcomes_separate', finished=time.time())
            atomic_write_json(a.out/'status.json', record)
            return
        time.sleep(20)
    record.update(status='effective_native_observation_window_ended_workers_untouched', finished=time.time())
    atomic_write_json(a.out/'status.json', record)


if __name__ == '__main__':
    main()
