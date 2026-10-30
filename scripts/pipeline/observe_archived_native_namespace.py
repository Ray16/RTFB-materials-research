#!/usr/bin/env python3
"""Wait for retained native attempts, then run their unchanged independent review."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.archived_native_namespace import mount_command, verify_spec


def alive(owner):
    host = owner['host'].split('.')[0]
    path = f'/proc/{int(owner["pid"])}/stat'
    if host == os.uname().nodename.split('.')[0]:
        p = Path(path)
        return p.exists() and p.read_text().rsplit(')', 1)[1].split()[19] == owner['process_start']
    code = ('from pathlib import Path; p=Path('+repr(path)+'); '
            'print(p.exists() and p.read_text().rsplit(")",1)[1].split()[19]=='+repr(owner['process_start'])+')')
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host,
                             shlex.join(['/usr/bin/python3', '-c', code])],
                            capture_output=True, text=True, timeout=20, check=True)
    if result.stdout.strip() not in ('True', 'False'):
        raise ValueError('Unexpected owner observation')
    return result.stdout.strip() == 'True'


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
    a.out.mkdir(parents=True, exist_ok=False)
    specs = json.loads(a.specs.read_text())['specs']
    if len(set(specs)) != len(specs):
        raise ValueError('Distinct namespace attempts required')
    hashes = {path: sha256_file(path) for path in specs}
    frozen = {str(Path(__file__).resolve()): sha256_file(Path(__file__)),
              str(ROOT/'src/redox/workflows/archived_native_namespace.py'):
                  sha256_file(ROOT/'src/redox/workflows/archived_native_namespace.py')}
    result = dict(status='waiting_for_finished_native_states_and_closed_ownership', started=time.time(),
                  owner=dict(host=os.uname().nodename, pid=os.getpid(),
                             process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
                  deadline=time.time()+a.hours*3600, source_sha256=frozen, specs_sha256=hashes,
                  completed_native_G_inferred_from_launch=False)
    done = {}
    while time.time() < result['deadline']:
        if any(sha256_file(path) != digest for path, digest in {**frozen, **hashes}.items()):
            raise ValueError('Observer source or namespace changed')
        rows = []
        for path in specs:
            if path in done:
                rows.append(done[path])
                continue
            spec = verify_spec(Path(path))
            actual = Path(spec['attempt'])
            stage = actual/'experimental_native/status.json'
            row = dict(spec=path, key=spec['key'], status='scientific_stage_not_yet_written')
            if stage.exists():
                state = json.loads(stage.read_text())
                row.update(scientific_status=state['status'], status='native_stage_running_or_owner_unverified')
                if 'finished' in state:
                    if state['status'] != 'experimental_native_finished_independent_review_pending':
                        row.update(status='scientific_failure_preserved', error=state.get('error'))
                        done[path] = row
                    else:
                        try:
                            live = alive(state['owner'])
                        except (OSError, subprocess.SubprocessError, ValueError):
                            live = None
                        guard = actual/'residency_guard/owner.json'
                        if live is False and guard.exists() and json.loads(guard.read_text()).get('guard_stopped'):
                            prior = Path(spec['prior_root'])
                            profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
                            review = actual/'independent_native_review.json'
                            if review.exists():
                                raise FileExistsError('Preserve existing independent review; inspect before reuse')
                            command = mount_command(spec)+['env', 'PYTHONPATH='+str(prior/'src'),
                                profile['python'], str(prior/'scripts/validation/review_smooth_polymer_native.py'),
                                '--campaign', spec['original_campaign'], '--out',
                                str(Path(spec['original_campaign'])/'independent_native_review.json')]
                            with (a.out/(spec['key']+'_review.log')).open('x') as log:
                                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
                            if completed.returncode:
                                row.update(status='independent_review_failed_preserved', returncode=completed.returncode)
                            else:
                                report = json.loads(review.read_text())
                                if report['status'] != 'alternate_native_state_independently_reconstructed':
                                    raise ValueError('Unexpected independent review scope')
                                row.update(status='native_state_independently_checked_under_declared_namespace',
                                           review=str(review), review_sha256=sha256_file(review), G_eV=report['G_eV'])
                            done[path] = row
            rows.append(row)
        result.update(updated=time.time(), entries=rows, terminal=len(done), expected=len(specs))
        atomic_write_json(a.out/'status.json', result)
        if len(done) == len(specs):
            result.update(status='all_declared_attempts_terminal_physical_outcomes_separate', finished=time.time())
            atomic_write_json(a.out/'status.json', result)
            return
        time.sleep(20)
    result.update(status='observation_window_ended_workers_untouched', finished=time.time())
    atomic_write_json(a.out/'status.json', result)


if __name__ == '__main__':
    main()
