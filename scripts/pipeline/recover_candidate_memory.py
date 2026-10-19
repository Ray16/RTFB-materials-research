#!/usr/bin/env python
"""One bounded, claimed recovery of a documented memory failure on a larger GPU."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from redox.validation.spins import atomic_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--task', type=int, required=True)
    a = ap.parse_args()
    manifest = json.loads(a.manifest.read_text())
    task = manifest['tasks'][a.task]
    assert task['kind'] == 'memory_recovery'
    assert socket.gethostname().split('.')[0] == task['host']
    assert sorted(os.sched_getaffinity(0)) == [task['controller_cpu']]
    for path, digest in manifest['producer_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    out = a.manifest.parent / task['id']
    out.mkdir()
    fleet = Path(manifest['original_fleet'])
    claim = fleet / 'claims' / task['id']
    claim.mkdir()  # Never reclaim another worker's claim.
    child = None
    start = Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]
    report = dict(task=task, status='preparing', host=socket.gethostname(), pid=os.getpid(),
                  process_start=start, started=time.time(), attempts=1, retries=0)
    def save():
        atomic_json(out / 'status.json', report)
    def stop(sig, frame):
        raise KeyboardInterrupt(f'signal {sig}')
    signal.signal(signal.SIGTERM, stop)
    atomic_json(claim / 'owner.json', dict(host=socket.gethostname(), pid=os.getpid(),
                process_start=start, worker='H100_memory_recovery', manifest=str(a.manifest)))
    save()
    try:
        # Move only explicitly inspected incomplete thermal work, holding both writer locks.
        with (ROOT / 'logs/pipeline_locks' / (task['id'] + '.lock')).open('a+') as ml, \
             (ROOT / 'logs/pipeline_locks' / (task['id'] + '_native_thermal.lock')).open('a+') as tl:
            fcntl.flock(ml, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(tl, fcntl.LOCK_EX | fcntl.LOCK_NB)
            inventories = []
            for partial in task['partials']:
                source = Path(partial)
                assert source.is_dir() and not source.with_suffix('.json').exists()
                sums = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in source.rglob('*') if p.is_file()}
                assert set(sums) == {'producer_source.py'}
                target = out / 'archived_partial' / source.parent.parent.name / source.name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(source, target)
                assert sums == {str(p.relative_to(target)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in target.rglob('*') if p.is_file()}
                inventories.append(dict(source=str(source), archive=str(target), sha256=sums,
                                        reason='Inspected original memory failure; retry on H100'))
            atomic_json(out / 'partial_archive.json', inventories)
        cpus = ','.join(map(str, task['compute_cpus']))
        argv = ['taskset', '-c', cpus, sys.executable, '-u',
                str(ROOT / 'scripts/pipeline/run_motifs.py'), '--config',
                str(ROOT / 'config/custom_candidates.py'), '--only', task['id'],
                '--out', str(out / 'pipeline'), '--execute', '--gpu', str(task['gpu']),
                '--threads', '2', '--cpus', cpus, '--nconf', '3',
                '--thermal-definition', 'native', '--skip-review']
        env = dict(os.environ, PYTHONPATH=str(ROOT / 'src'))
        report.update(status='running', argv=argv)
        save()
        contention = 0
        deadline = json.loads((fleet / 'queue.json').read_text())['stop_new_work_at_unix']
        while True:
            with (out / f'pipeline_{contention:02d}.log').open('w') as output:
                child = subprocess.Popen(argv, cwd=ROOT, env=env, stdout=output,
                                         stderr=subprocess.STDOUT, start_new_session=True)
                report['pipeline_pid'] = child.pid
                save()
                rc = child.wait(timeout=12 * 3600)
                child = None
            if rc == 0:
                break
            pm = json.loads((out / 'pipeline/run_manifest.json').read_text())
            last = Path(pm['steps'][-1]['log'])
            if 'BUSY/reserved -> not running' not in last.read_text():
                break
            report.update(status='waiting_for_reservation', resource_contentions=contention+1)
            save()
            if time.time() >= deadline:
                raise RuntimeError('Reservation contention deadline; no scientific attempt consumed')
            if '--resume' not in argv:
                argv.append('--resume')
            contention += 1
            time.sleep(15)
        report['returncode'] = rc
        if rc:
            raise RuntimeError('Targeted recovery failed; preserve all outputs without another automatic attempt')
        from redox.core.common import read_manifest
        from redox.validation.readiness import inspect_states
        pm = json.loads((out / 'pipeline/run_manifest.json').read_text())
        qc = inspect_states([r for r in read_manifest() if r['id'] == task['id']])
        if pm.get('status') != 'selected_outputs_verified' or not qc or any(r['status'] != 'ready' for r in qc):
            raise RuntimeError('Required selected outputs did not pass verification')
        report.update(status='selected_records_verified', readiness=qc,
                      scientific_accuracy_validated=False, finished=time.time())
        atomic_json(fleet / 'done' / (task['id'] + '.json'), report)
        save()
    except BaseException as exc:
        report.update(status='failed_preserved', error=repr(exc), finished=time.time())
        save()
        raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        shutil.rmtree(claim)


if __name__ == '__main__':
    main()
