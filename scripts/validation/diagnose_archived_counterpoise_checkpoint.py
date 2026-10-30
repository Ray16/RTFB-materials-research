#!/usr/bin/env python3
"""Recompute one saved CPU orbital gradient without admitting an SCF result."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.workflows.archived_native_namespace import remaining_host_promises
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('archive', 'prior-root', 'checkpoint-root', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--maximum-rss-gib', type=float, default=32)
    a = p.parse_args()
    if a.threads < 1 or a.maximum_rss_gib != 32:
        raise ValueError('Positive CPU count and the retained 32 GiB diagnostic bound required')
    if a.out.exists():
        raise FileExistsError('Preserve existing diagnostic evidence')
    reader = ArchivedSources(a.prior_root, a.archive, relative_to_prior=True)
    original_raw = a.checkpoint_root
    raw = reader.path(original_raw)
    launch = json.loads((raw/'launch.json').read_text())
    status = json.loads((raw/'status.json').read_text())
    if launch['role'] != 'AB' or status['role'] != 'AB':
        raise ValueError('This diagnostic is the retained AB checkpoint only')
    owner = status['owner']
    owner_path = Path('/proc')/str(owner['pid'])/'stat'
    if owner['host'] != os.uname().nodename:
        raise ValueError('Establish remote producer termination explicitly before using this entry point')
    if owner_path.exists() and owner_path.read_text().rsplit(')', 1)[1].split()[19] == owner['process_start']:
        raise ValueError('Saved producer is still live')
    reader.verify_all(launch['source_sha256'])
    plan_path = Path(launch['case'])/'plan.json'
    plan = json.loads(reader.path(plan_path).read_text())
    geometry = Path(plan['input_geometry'])
    evaluator = a.prior_root/'scripts/validation/check_counterpoise_remaining_terms.py'
    inputs = [original_raw/'launch.json', original_raw/'status.json', original_raw/'AB/SCF.h5',
              original_raw/'grid_coordinates_AU.npy', original_raw/'grid_weights.npy',
              original_raw/'AB/progress.json', plan_path, geometry, evaluator,
              a.prior_root/'scripts/validation/local_counterpoise.py',
              a.prior_root/'scripts/validation/run_cpu_counterpoise_reference.py']
    sources = {str(path): sha256_file(reader.path(path)) for path in inputs}
    reader.verify_all(sources)
    reference = json.loads(reader.path(original_raw/'AB/progress.json').read_text())['E_Ha']
    with allocation_guard():
        capacity = host_capacity(a.threads, minimum_ram_gib=80)
        promises = remaining_host_promises()
        if '--maximum-rss-gib' in sys.argv:
            promises -= max(0, int(32*2**30)-psutil.Process().memory_info().rss)
        if capacity['available_RAM_bytes']-promises < 64*2**30:
            raise RuntimeError('Insufficient RAM after accounting for live worker promises')
        cpus = cpu_allocation(a.threads)
        os.sched_setaffinity(0, cpus)
        a.out.mkdir(parents=True)
        state = dict(status='saved_AB_checkpoint_orbital_diagnostic_running', started=time.time(),
                     owner=dict(pid=os.getpid(), host=os.uname().nodename,
                                process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
                     source_sha256=sources, CPU_affinity=cpus, threads=a.threads,
                     producer_sha256=sha256_file(Path(__file__)), source_bindings=len(reader.receipts),
                     full_counterpoise_qualified=False, SCF_convergence_inferred=False)
        atomic_write_json(a.out/'status.json', state)
    program = ('import importlib.util,json,pathlib; from pyscf import lib; lib.num_threads('+str(a.threads)+'); '
               'p=pathlib.Path('+repr(str(evaluator))+'); '
               's=importlib.util.spec_from_file_location("saved_CP_orbital_diagnostic",p); '
               'm=importlib.util.module_from_spec(s); s.loader.exec_module(m); '
               'g=json.loads(pathlib.Path('+repr(str(geometry))+').read_text()); '
               'r=m.evaluate(g,"AB",pathlib.Path('+repr(str(original_raw/'AB/SCF.h5'))+'),'
               'pathlib.Path('+repr(str(original_raw/'grid_coordinates_AU.npy'))+'),'
               'pathlib.Path('+repr(str(original_raw/'grid_weights.npy'))+'),'+repr(reference)+'); '
               'r.update(full_counterpoise_qualified=False,SCF_convergence_inferred=False); '
               'pathlib.Path('+repr(str(a.out.resolve()/'diagnostic.json'))+').write_text(json.dumps(r,indent=2)+"\\n")')
    profile = json.loads((a.archive/'environments/cutensor22/runtime.json').read_text())
    cmd = ['bwrap', '--die-with-parent', '--ro-bind', '/', '/', '--ro-bind', str(a.archive), str(a.prior_root),
           '--bind', str(a.out.resolve()), str(a.out.resolve()), '--proc', '/proc', '--dev', '/dev',
           '--tmpfs', '/tmp', '--chdir', str(a.prior_root), '--', profile['python'], '-c', program]
    env = dict(os.environ, PYTHONPATH=str(a.prior_root/'src'), PYTHONDONTWRITEBYTECODE='1', PYSCF_TMPDIR='/tmp')
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        env[key] = str(a.threads)
    atomic_write_json(a.out/'command.json', dict(argv=cmd, CPU_affinity=cpus, source_sha256=sources))
    with (a.out/'worker.log').open('x') as log:
        child = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                 env=env, start_new_session=True)
        peak = 0
        while child.poll() is None:
            processes = psutil.Process(child.pid).children(recursive=True)+[psutil.Process(child.pid)]
            total = 0
            for process in processes:
                try:
                    total += process.memory_info().rss
                except psutil.NoSuchProcess:
                    pass
            peak = max(peak, total)
            if total > 32*2**30 or psutil.virtual_memory().available < 32*2**30:
                os.killpg(child.pid, signal.SIGTERM)
                state.update(status='owned_diagnostic_resource_stop_requested')
            time.sleep(2)
    reader.recheck()
    state.update(returncode=child.returncode, peak_RSS_bytes=peak, finished=time.time())
    if child.returncode == 0:
        result = json.loads((a.out/'diagnostic.json').read_text())
        state.update(status='saved_checkpoint_diagnostic_finished_not_SCF_admission',
                     orbital_gradient_norm_au=result['orbital_gradient_norm_au'],
                     orbital_gradient_threshold_au=result['orbital_gradient_threshold_au'])
    else:
        state['status'] = 'checkpoint_diagnostic_failed_preserved'
    atomic_write_json(a.out/'status.json', state)
    print(json.dumps(state))


if __name__ == '__main__':
    main()
