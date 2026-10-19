"""One private vinyl-monomer solution state; all failures remain reviewable."""
import argparse
from contextlib import ExitStack, redirect_stdout, redirect_stderr
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
        row = task['manifest_row']
        assert row == task['manifest_row'], 'Selected manifest input changed'
        private_root = base/'workspace'
        state_dir = private_root/'calcs/dft'/task['id']/task['state']
        assert not state_dir.exists(), 'Existing state requires a separate resume/adoption review'
        locks = ROOT / 'logs/pipeline_locks'
        with ExitStack() as stack:
            # Shared molecule guards exclude the existing exclusive repair/motif paths.
            # Other tasks in this panel own different state directories.
            for suffix in ('', '_native_thermal'):
                lock = stack.enter_context((locks / f"{task['id']}{suffix}.lock").open('a+'))
                fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
            lock = stack.enter_context((locks / f"{task['key']}_aryl_solution.lock").open('a+'))
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            stack.callback(stop_children)
            selection = task['id'] + ':' + task['state']
            import redox.qm.dft as dft
            import pyscf.lib
            pyscf.lib.num_threads(2)
            dft.ROOT = private_root  # Process-local paths; shared configuration stays canonical.
            original_smd = dft.dft_smd
            def tighter_smd(*args, **kwargs):
                if kwargs.get('do_opt',True):
                    kwargs['conv_params']=plan['optimizer_stopping_tolerances']
                return original_smd(*args,**kwargs)
            dft.dft_smd=tighter_smd
            stack.callback(setattr,dft,'dft_smd',original_smd)
            for stage in ('dft','active_sp'):
                report.update(status='computing',stage=stage)
                save()
                with (destination/f'{stage}.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
                    if stage=='dft':
                        assert dft.run_batch([row],False,do_opt=True,do_freq=False,backend='gpu',nconf=3)==0
                    else:
                        from redox.qm.sp import ensure
                        raw=json.loads((state_dir/'result.json').read_text())
                        ensure(task['id'],task['state'],state_dir,raw['charge'],raw['mult'],backend='gpu')
                if stage == 'dft':
                    from redox.validation.protonation import hydrogen_topology_status
                    topology = hydrogen_topology_status(row['smiles'], state_dir / 'opt.xyz')
                    report['hydrogen_topology'] = topology
                    save()
                    if topology != 'ok':
                        raise RuntimeError('Declared acid/phenol graph changed during solution relaxation: '+topology)
            report.update(stage='native_thermal')
            save()
            import redox.qm.native_thermal as native
            native.read_manifest=lambda:[row]
            native.select_policy(state_dir)
            with (destination/'native_thermal.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
                _,thermal=native.ensure_native_thermal(task['id'],task['state'],state_dir=state_dir,memory_cap_gib=24)
            from redox.validation.readiness import inspect_states
            qc = inspect_states([row],root=private_root/'calcs/dft')
            assert len(qc) == 1 and qc[0]['status'] == 'ready', qc
            for path, expected in {**plan['producer_sha256'], **task['input_sha256']}.items():
                assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected,path
            report.update(status='state_records_verified', readiness=qc,
                          ground_spin_selection_complete=False,
                          scope='private_vinyl_monomer_pilot',
                          proposed_spin_surface_only=True,
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
