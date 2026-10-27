#!/usr/bin/env python
"""Resume source-bound smooth geometry batches with failure-isolated inventory.

Scientific workers and acceptance criteria remain in the original geometry queue.
Shared discovery is cached briefly; each actual placement rechecks the gate.
"""
import fcntl
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows import resources

PROFILE = ROOT / 'environments/cutensor22/runtime.json'
ORIGINAL = ROOT / 'scripts/pipeline/coordinate_geometry_batch.py'
INVENTORY = ROOT / 'scripts/pipeline/coordinate_native_inventory_resilient.py'
CPU_PREFLIGHT = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'
CACHE = ROOT / 'calcs/validation/smooth_geometry_inventory_cache'


def cache_valid(record, hosts, policy_hash, now):
    return (record.get('hosts') == sorted(hosts) and record.get('policy_sha256') == policy_hash
            and isinstance(record.get('checked'), (int, float))
            and 0 <= now - record['checked'] <= 30
            and isinstance(record.get('slots'), list) and isinstance(record.get('checks'), list))


def shared_inventory(hosts, inventory):
    CACHE.mkdir(parents=True, exist_ok=True)
    with (CACHE / 'lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = CACHE / 'inventory.json'
        try:
            row = json.loads(path.read_text())
        except (OSError, ValueError):
            row = {}
        policy_hash = api.digest(inventory.POLICY)
        if not cache_valid(row, hosts, policy_hash, time.time()):
            slots, checks = inventory.query_all(hosts)
            row = dict(hosts=sorted(hosts), policy_sha256=policy_hash, checked=time.time(),
                       slots=slots, checks=checks)
            atomic_write_json(path, row)
        return [tuple(s) for s in row['slots']], row['checks']


def make_preflight(profile):
    def preflight(batch, host):
        path = batch / (host + '_preflight.json')
        try:
            row = json.loads(path.read_text())
            if (time.time() - path.stat().st_mtime < 600 and row.get('ok')
                    and row.get('python') == profile['python']
                    and row.get('GPU_imports_blocked_in_CPU_preflight') is True
                    and row.get('GPU_modules_loaded') == []):
                return True
        except (OSError, ValueError):
            pass
        argv = [sys.executable, str(CPU_PREFLIGHT), host, '--python', profile['python'],
                '--module', 'pyscf', '--module', 'gpu4pyscf', '--module', 'cupy',
                '--module', 'rdkit', '--module', 'ase', '--module', 'scipy',
                '--path', str(batch / 'batch.json'), '--bin', str(ROOT.parent / 'bin/gpu_reserve'),
                '--library-path', profile['library_path'], '--min-available-ram-mib', '65536']
        try:
            result = subprocess.run(argv, capture_output=True, text=True, timeout=90)
            row = json.loads(result.stdout)
            valid = (result.returncode == 0 and row.get('ok') is True
                     and row.get('python') == profile['python']
                     and row.get('GPU_imports_blocked_in_CPU_preflight') is True
                     and row.get('GPU_modules_loaded') == [])
            if not valid:
                row['ok'] = False
            atomic_write_json(path, row)
            return valid
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            atomic_write_json(path, dict(ok=False, host=host, error=repr(error), checked=time.time()))
            return False
    return preflight


def adapt(module, profile):
    module.SCRIPT = Path(__file__).resolve()
    module.preflight = make_preflight(profile)
    body = inspect.getsource(module.launch_local)
    anchor = "        with (attempt / 'worker.log').open('x') as log:"
    if body.count(anchor) != 1:
        raise ValueError('Atomic geometry placement template changed')
    inserted = """        if plan.get('runtime_profile') != str(_PROFILE):
            raise ValueError('Task runtime differs from the model batch')
        argv.insert(argv.index('env')+1, 'LD_LIBRARY_PATH='+_profile['library_path'])
        index=argv.index(manifest['worker'])
        argv[index:index]=[str(ROOT/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(_PROFILE),
            '--out', str(campaign/'tensor_runtime')]
"""
    namespace = dict(vars(module), _profile=profile, _PROFILE=PROFILE)
    local = {}
    exec(compile(body.replace(anchor, inserted + anchor, 1), __file__, 'exec'), namespace, local)
    module.launch_local = local['launch_local']
    body = inspect.getsource(module.coordinate)
    old = 'stationary_geometry_checked_active_parent_and_native_minimum_pending'
    if body.count(old) != 1:
        raise ValueError('Geometry result contract changed')
    body = body.replace(old, 'experimental_smooth_geometry_checked_energy_and_native_QC_pending', 1)
    old = "                            assignments[entry['key']] = dict(owner=dict(host=host, pid=-1, process_start='unknown'), error=repr(error))\n"
    new = ("                            assignments[entry['key']] = dict(owner=dict(host=host, pid=-1, process_start='unknown'), error=repr(error))\n"
           "                            atomic_write_json(Path(entry['campaign']) / 'placement_error.json', dict(status='placement_response_unverified_preserved', owner=owner, host=host, GPU=gpu, error=repr(error), checked=time.time()))\n")
    if body.count(old) != 1:
        raise ValueError('Placement error contract changed')
    local = {}
    exec(compile(body.replace(old, new, 1), __file__, 'exec'), dict(vars(module)), local)
    module.coordinate = local['coordinate']


def main():
    os.chdir(ROOT)
    profile = json.loads(PROFILE.read_text())
    if Path(sys.executable) != Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH'):
        raise ValueError('Use the declared interpreter and host loader environment')
    batch = Path(sys.argv[sys.argv.index('--batch') + 1]).resolve()
    manifest = json.loads((batch / 'batch.json').read_text())
    if (manifest['stage'] != 'geometry' or manifest['worker'] != str(ROOT / 'scripts/pipeline/run_smooth_geometry.py')
            or manifest['primary_candidate_admission']):
        raise ValueError('Only declared alternate-model geometry batches are supported')
    api.verify_sources(manifest['source_sha256'])
    inventory = api.load_module('geometry_per_host_inventory', INVENTORY)
    hosts = (sys.argv[sys.argv.index('--hosts') + 1].split(',') if '--hosts' in sys.argv
             else [os.uname().nodename.split('.')[0]])
    if len(hosts) != len(set(hosts)) or not set(hosts) <= inventory.permitted_hosts(inventory.POLICY):
        raise ValueError('Use distinct policy-permitted hosts')
    proof = {str(p): api.digest(p) for p in (Path(__file__).resolve(), ORIGINAL, INVENTORY, CPU_PREFLIGHT,
                                             PROFILE, inventory.GATE, inventory.POLICY)}
    module = api.load_module('resilient_smooth_geometry_queue', ORIGINAL)
    adapt(module, profile)
    out = batch / 'resilient_geometry_inventory'
    out.mkdir(exist_ok=True)
    owner = dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    original = resources.free_slots
    def free_slots(requested=None):
        api.verify_sources(proof)
        chosen = hosts if requested is None else list(requested)
        if not set(chosen) <= inventory.permitted_hosts(inventory.POLICY):
            raise ValueError('Current policy disallows requested host')
        slots, checks = (inventory.query_all(chosen) if len(chosen) == 1 else shared_inventory(chosen, inventory))
        if len(chosen) > 1:
            atomic_write_json(out / 'status.json', dict(owner=owner, checked=time.time(), checks=checks,
                              available_slots=slots, source_sha256=proof))
        return slots
    resources.free_slots = free_slots
    try:
        module.main()
    finally:
        resources.free_slots = original


if __name__ == '__main__':
    main()
