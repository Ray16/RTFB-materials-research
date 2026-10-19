"""One new exact graft/state per reservation; all failures remain reviewable."""
import argparse
from contextlib import ExitStack
import fcntl
import hashlib
import json
import os
from pathlib import Path
import socket
import signal
import psutil
import subprocess
import sys
import time
import traceback

from redox.core.common import ROOT, DFT, read_manifest
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--task', type=int, required=True)
    a = ap.parse_args()
    base = Path(__file__).resolve().parent
    plan = json.loads((base / 'plan.json').read_text())
    task = plan['tasks'][a.task]
    destination = base / 'tasks' / task['key']
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(status='preparing', task=task, host=socket.gethostname(), pid=os.getpid(),
                  process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                  cpu_affinity=sorted(os.sched_getaffinity(0)), started=time.time(),
                  experimental_accuracy_validated=False, primary_table_promotion=False)

    def save():
        atomic_write_json(destination / 'status.json', report)

    save()
    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    def stop_children():
        children=psutil.Process().children(recursive=True)
        for child in reversed(children):
            try:child.terminate()
            except psutil.NoSuchProcess:pass
        _,alive=psutil.wait_procs(children,timeout=10)
        for child in alive:
            try:child.kill()
            except psutil.NoSuchProcess:pass
        psutil.wait_procs(alive,timeout=10)
    try:
        verify_gate_owner()
        assert len(os.sched_getaffinity(0)) == 2
        for name in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
                     'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
            assert os.environ.get(name) == '2', name
        for path, expected in {**plan['producer_sha256'], **task['input_sha256']}.items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, path
        row = next(r for r in read_manifest() if r['id'] == task['id'] and r['state'] == task['state'])
        assert row == task['manifest_row'], 'Selected manifest input changed'
        state_dir = DFT / task['id'] / task['state']
        assert not state_dir.exists(), 'Existing state requires a separate resume/adoption review'
        locks = ROOT / 'logs/pipeline_locks'
        with ExitStack() as stack:
            # Shared molecule guards exclude the existing exclusive repair/motif paths.
            # Other tasks in this panel own different state directories.
            for suffix in ('', '_native_thermal'):
                lock = stack.enter_context((locks / f"{task['id']}{suffix}.lock").open('a+'))
                fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
            lock = stack.enter_context((locks / f"{task['key']}_cata_solution.lock").open('a+'))
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            stack.callback(stop_children)
            selection = task['id'] + ':' + task['state']
            for stage, argv in (
                ('dft', [sys.executable, '-u', '-m', 'redox.qm.dft', '--only', selection,
                         '--backend', 'gpu', '--nthreads', '2', '--nconf', '3', '--no-freq']),
                ('active_sp', [sys.executable, '-u', '-m', 'redox.qm.sp', '--only', selection,
                               '--backend', 'gpu']),
            ):
                report.update(status='computing', stage=stage, argv=argv)
                save()
                with (destination / f'{stage}.log').open('w') as log:
                    subprocess.run(argv, cwd=ROOT, check=True, stdout=log, stderr=subprocess.STDOUT)
            report.update(stage='native_thermal')
            save()
            from redox.qm.native_thermal import ensure_native_thermal, select_policy
            select_policy(state_dir)
            _, thermal = ensure_native_thermal(task['id'], task['state'])
            from redox.validation.readiness import inspect_states
            qc = inspect_states([row])
            assert len(qc) == 1 and qc[0]['status'] == 'ready', qc
            for path, expected in {**plan['producer_sha256'], **task['input_sha256']}.items():
                assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected,path
            report.update(status='state_records_verified', readiness=qc,
                          ground_spin_selection_complete=False,
                          thermal_eV=thermal['g_thermal_total_eV'],
                          output_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                         for p in state_dir.rglob('*') if p.is_file()
                                         and p.suffix in ('.json', '.xyz', '.npy')})
    except BaseException as exc:
        report.update(status='failed_preserved', error=repr(exc), traceback=traceback.format_exc())
        raise
    finally:
        report['finished'] = time.time()
        save()


if __name__ == '__main__':
    main()
