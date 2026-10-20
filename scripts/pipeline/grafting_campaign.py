#!/usr/bin/env python
"""Finite, reserved GPU campaign for balanced Merrifield-site grafting in MeCN."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stderr, redirect_stdout
import csv
import gc
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time
import traceback

from redox.core.common import ROOT, HARTREE_EV, apply_thermal_policy, free_energy
from redox.core.io import atomic_write_json

BASE = ROOT / 'results/grafting/campaigns/2026-10-19_mecn_reference'
SCRIPT = Path(__file__).resolve()
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
           'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes(directory):
    return {str(p): digest(p) for p in sorted(Path(directory).rglob('*'))
            if p.is_file() and p.suffix in ('.json', '.xyz', '.npy')}


def verify(mapping):
    for name, expected in mapping.items():
        if digest(name) != expected:
            raise RuntimeError('Source or result changed: ' + name)


def build_input(spec):
    from rdkit import Chem
    from redox.build.structures import conformer_ensemble, to_xyz
    # Each preparer has one disjoint CPU and single-thread native libraries.
    allowed = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {allowed[mp.current_process()._identity[0] % len(allowed)]})
    mol = Chem.MolFromSmiles(spec['smiles'])
    if mol is None or Chem.GetFormalCharge(mol) != spec['charge']:
        raise ValueError('Invalid declared species: ' + spec['id'])
    molh = Chem.AddHs(mol)
    electrons = sum(a.GetAtomicNum() for a in molh.GetAtoms()) - spec['charge']
    assert (electrons - (spec['mult'] - 1)) % 2 == 0
    directory = BASE / 'workspace/library' / spec['id']
    directory.mkdir(parents=True, exist_ok=False)
    conformers = directory / 'conformers'
    conformers.mkdir()
    if spec['id'] == 'chloride':
        xyz = '1\nCl- isolated atom\nCl 0.0 0.0 0.0\n'
        (conformers/'conf_00.xyz').write_text(xyz)
        count = 1
    elif spec['id'] == 'hydrogen_chloride':
        # Non-axis-aligned diatomic supports the native finite-gradient checks.
        length = 1.275 / math.sqrt(3)
        xyz = f'2\nHCl neutral reference\nCl 0.0 0.0 0.0\nH {length} {length} {length}\n'
        (conformers/'conf_00.xyz').write_text(xyz)
        count = 1
    else:
        molh, cids, ff, energies = conformer_ensemble(mol, net_charge=spec['charge'], top_k=3)
        for i, (cid, energy) in enumerate(zip(cids, energies)):
            one = Chem.Mol(molh)
            one.RemoveAllConformers()
            one.AddConformer(molh.GetConformer(cid), assignId=True)
            (conformers/f'conf_{i:02d}.xyz').write_text(to_xyz(one, f'{ff} seed E={energy}'))
        count = len(cids)
    shutil.copyfile(conformers/'conf_00.xyz', directory/'rest.xyz')
    row = dict(id=spec['id'], name=spec['id'], family='grafting_component', state='rest',
               charge=spec['charge'], mult_hint=spec['mult'], n_e=0, n_conf=count,
               solv_preopt=int(spec['charge'] != 0), smiles=spec['smiles'])
    return dict(key=spec['id'], id=spec['id'], state='rest', role=spec['role'],
                atoms=molh.GetNumAtoms(), manifest_row=row, input_sha256=hashes(directory),
                atom=molh.GetNumAtoms() == 1,
                preopt='none' if molh.GetNumAtoms() <= 2 else 'alpb')


def prepare():
    from redox.core.common import read_result
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    assert not (BASE/'plan.json').exists(), 'Preserve the existing campaign and inspect it'
    BASE.mkdir(parents=True, exist_ok=True)
    for name in ('claims', 'done', 'failures', 'tasks', 'workers', 'workspace/library'):
        (BASE/name).mkdir(parents=True, exist_ok=True)
    actual = ROOT/'calcs/grafting'/BASE.name
    actual.mkdir(parents=True, exist_ok=False)
    (BASE/'workspace/calcs').mkdir()
    (BASE/'workspace/calcs/dft').symlink_to(actual, target_is_directory=True)
    inventory = json.loads((ROOT/'results/grafting/species.json').read_text())
    reactions = json.loads((ROOT/'results/grafting/reactions.json').read_text())
    specs = inventory['shared_references'] + inventory['precursors']
    specs = [dict(s, role=s.get('role', 'shared_reference')) for s in specs]
    cached = {}
    for r in reactions['reactions']:
        raw = read_result(r['product_id'], r['product_state'])
        G, qc = free_energy(raw)
        if G is not None and qc == 'ok' and raw.get('thermal_definition') == 'native_solution_dft_qrrho_spin':
            state_dir = ROOT/'calcs/dft'/r['product_id']/r['product_state']
            cached[r['id']] = dict(G_eV=G, state_dir=str(state_dir), source_sha256=hashes(state_dir))
        else:
            specs.append(dict(id='product_'+r['id'], smiles=r['product_smiles'],
                              charge=r['product_charge'], mult=1, role='missing_product'))
    atomic_write_json(BASE/'reused_products.json', cached)
    # Four CPU preparers; no GPU calculation occurs during graph embedding.
    with allocation_guard():
        host_capacity(4, minimum_ram_gib=8)
        cpus = cpu_allocation(4)
        os.sched_setaffinity(0, cpus)
    with ProcessPoolExecutor(max_workers=4, mp_context=mp.get_context('spawn')) as pool:
        tasks = list(pool.map(build_input, specs))
    # Shared references first, then the initial representative panel, then by size.
    reactions_by_precursor = {r['precursor_id']: r['id'] for r in reactions['reactions']}
    pilots = set(json.loads((ROOT/'results/polymer_models/status.json').read_text())['pilot_candidates'])
    tasks.sort(key=lambda t: (t['role'] != 'shared_reference',
                             reactions_by_precursor.get(t['id']) not in pilots, t['atoms'], t['key']))
    producer_files = list((ROOT/'src/redox').rglob('*.py')) + list((ROOT/'config').glob('*.json'))
    producer_files += list((ROOT/'config').glob('*.py')) + [SCRIPT, ROOT/'scripts/pipeline/collect_grafting_results.py']
    inputs = [ROOT/'results/grafting/reactions.json', ROOT/'results/grafting/species.json', BASE/'reused_products.json']
    plan = dict(schema=1, created=time.time(), solvent='acetonitrile', temperature_K=298.15,
                tasks=tasks, reused_products=len(cached), reaction_count=len(reactions['reactions']),
                raw_root=str(actual), interpreter=sys.executable, attempts_per_item=1,
                worker_consecutive_failure_limit=3, source_snapshot=reactions['source_snapshot'],
                producer_sha256={str(p): digest(p) for p in sorted(set(producer_files))},
                input_sha256={str(p): digest(p) for p in inputs},
                optimizer_stopping_tolerances=dict(convergence_gmax=5e-5, convergence_grms=3e-5,
                    convergence_energy=1e-7, convergence_dmax=.001, convergence_drms=.0005),
                spin_scope='Closed-shell precursor/reference and resting-product surfaces; not a general ground-spin proof',
                base='unspecified; NH/OH reference retains HCl explicitly')
    atomic_write_json(BASE/'plan.json', plan)
    print(json.dumps(dict(prepared_tasks=len(tasks), reused_products=len(cached), raw_root=str(actual))))


def monatomic(sd, row):
    """Atomic translation/spin uses the qualified kernel; no fictitious Hessian."""
    from redox.qm.sp import ensure
    from redox.core.protocol import load_record
    from redox.qm.qrrho import atomic_masses, kernel_thermal
    raw = dict(id=row['id'], state='rest', charge=-1, mult=1, optimized=False,
               geometry_definition='isolated atom; internal optimization inapplicable')
    atomic_write_json(sd/'result.json', raw)
    ensure(row['id'], 'rest', sd, -1, 1, backend='gpu')
    rec = load_record(sd, -1, 1)
    if rec is None or not math.isfinite(rec['e_smd_eV']):
        raise RuntimeError('Atomic reference SP failed protocol/convergence admission')
    if abs(rec.get('smd_spin_contam', math.inf)) > .1:
        raise RuntimeError('Atomic spin QC failed')
    mass = float(atomic_masses([17])[0])
    thermal = kernel_thermal([], [0., 0., 0.], 0., mass, 1, atom=True)
    T = 298.15
    entropy = 2.2868*(5*math.log10(T)+3*math.log10(mass))-2.3135
    analytic = (2.5*1.98726*T-T*entropy)/1000/627.50947428*HARTREE_EV
    g = thermal['G_Ha']*HARTREE_EV
    if not math.isclose(g, analytic, abs_tol=1e-10):
        raise RuntimeError('Independent monatomic partition check failed')
    atomic_write_json(sd/'monatomic_thermal.json', dict(
        definition='monatomic ideal translation and singlet spin; 1 atm convention',
        temperature_K=T, mass_amu=mass, internal_modes=0, G_thermal_eV=g,
        analytic_G_thermal_eV=analytic, kernel=thermal,
        kernel_source_sha256=digest(ROOT/'src/redox/qm/qrrho.py')))
    return rec['e_smd_eV']+g


def run_state(task, plan):
    from redox.workflows.resources import verify_gate_owner
    import redox.qm.dft as dft
    import pyscf.lib
    from redox.qm.sp import ensure
    import redox.qm.native_thermal as native
    from redox.validation.readiness import inspect_states
    verify_gate_owner()
    verify(task['input_sha256'])
    dest = BASE/'tasks'/task['key']
    dest.mkdir(exist_ok=False)
    sd = BASE/'workspace/calcs/dft'/task['id']/task['state']
    row = task['manifest_row']
    report = dict(status='computing', task=task, started=time.time(), host=socket.gethostname(),
                  pid=os.getpid(), stage='dft', primary_candidate_promotion=False)
    def save(): atomic_write_json(dest/'status.json', report)
    save()
    original = dft.dft_smd
    previous_root = dft.ROOT
    previous_manifest = native.read_manifest
    def tighter(*args, **kwargs):
        if kwargs.get('do_opt', True):
            kwargs['conv_params'] = plan['optimizer_stopping_tolerances']
        return original(*args, **kwargs)
    try:
        pyscf.lib.num_threads(2)
        dft.ROOT = BASE/'workspace'
        dft.dft_smd = tighter
        native.read_manifest = lambda: [row]
        with (dest/'calculation.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
            if task['atom']:
                sd.mkdir(parents=True, exist_ok=False)
                shutil.copy(BASE/'workspace/library'/task['id']/'rest.xyz', sd/'opt.xyz')
                G = monatomic(sd, row)
                readiness = 'qualified_monatomic_reference'
            else:
                failed = dft.run_batch([row], False, do_opt=True, do_freq=False, backend='gpu',
                                       nconf=3, preopt=task['preopt'])
                if failed: raise RuntimeError('All conformer optimizations failed')
                report['stage'] = 'active_sp'; save()
                raw = json.loads((sd/'result.json').read_text())
                ensure(task['id'], task['state'], sd, raw['charge'], raw['mult'], backend='gpu')
                from redox.validation.protonation import hydrogen_topology_status
                topology = hydrogen_topology_status(row['smiles'], sd/'opt.xyz')
                if topology != 'ok': raise RuntimeError('Declared hydrogen/connectivity changed: '+topology)
                report['stage'] = 'native_thermal'; save()
                native.select_policy(sd)
                native.ensure_native_thermal(task['id'], task['state'], state_dir=sd, memory_cap_gib=24)
                readiness = inspect_states([row], root=BASE/'workspace/calcs/dft')
                if readiness[0]['status'] != 'ready': raise RuntimeError(str(readiness))
                from redox.core.protocol import load_record
                rec = load_record(sd, raw['charge'], raw['mult'])
                G, qc = free_energy(apply_thermal_policy(dict(raw, e_smd_eV=rec['e_smd_eV']), sd))
                if G is None or qc != 'ok': raise RuntimeError('Component G not admitted: '+qc)
        verify(plan['producer_sha256']); verify(plan['input_sha256'])
        report.update(status='state_records_verified', G_eV=G, readiness=readiness,
                      state_dir=str(sd.resolve()), output_sha256=hashes(sd))
    except BaseException as exc:
        report.update(status='failed_preserved', error=repr(exc), traceback=traceback.format_exc())
        raise
    finally:
        dft.dft_smd = original
        dft.ROOT = previous_root
        native.read_manifest = previous_manifest
        report['finished'] = time.time(); save()


def worker():
    from redox.workflows.resources import verify_gate_owner, gpu_suitability
    from redox.workflows.checkpoints import verify_state_checkpoint
    from redox.workflows.processes import stop_owned_children
    verify_gate_owner()
    plan = json.loads((BASE/'plan.json').read_text())
    verify(plan['producer_sha256']); verify(plan['input_sha256'])
    from redox.qm.native_thermal import implementation_controls
    implementation_controls()
    assert len(os.sched_getaffinity(0)) == 2
    assert all(os.environ.get(k) == '2' for k in THREADS)
    device = gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    large = device['memory_MiB'] >= 80*1024
    owner = dict(host=socket.gethostname(), pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                 cpus=sorted(os.sched_getaffinity(0)), gpu=int(os.environ['CUDA_VISIBLE_DEVICES']))
    report = dict(owner, status='claiming', started=time.time(), items=[])
    path = BASE/'workers'/f"{owner['host']}_gpu{owner['gpu']}.json"
    def stop(signum, frame): raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, stop)
    consecutive = 0
    try:
        for task in plan['tasks']:
            key = task['key']
            # The large/unknown molecular memory profile stays on H100s.
            if not large and task['atoms'] > 48: continue
            done = BASE/'done'/f'{key}.json'
            if done.exists():
                verify_state_checkpoint(done); continue
            if (BASE/'failures'/f'{key}.json').exists(): continue
            claim = BASE/'claims'/key
            try: claim.mkdir()
            except FileExistsError: continue
            atomic_write_json(claim/'owner.json', owner)
            try:
                if done.exists(): verify_state_checkpoint(done); continue
                if (BASE/'tasks'/key).exists():
                    raise RuntimeError('Unreceipted prior attempt requires owner review: '+key)
                report.update(status='computing', active_key=key)
                atomic_write_json(path, report)
                print('Actual calculation start', key, flush=True)
                try:
                    run_state(task, plan)
                except Exception as exc:
                    if not (BASE/'tasks'/key/'status.json').exists():
                        raise RuntimeError('Pre-execution failure; stop before consuming more inputs') from exc
                    atomic_write_json(BASE/'failures'/f'{key}.json', dict(key=key, owner=owner,
                        attempts=1, error=repr(exc), state_status=str(BASE/'tasks'/key/'status.json'), finished=time.time()))
                    report['items'].append(dict(key=key, status='failed_preserved'))
                    consecutive += 1
                else:
                    status_path = BASE/'tasks'/key/'status.json'
                    atomic_write_json(done, dict(key=key, owner=owner, status_path=str(status_path),
                                                status_sha256=digest(status_path), finished=time.time()))
                    verify_state_checkpoint(done)
                    report['items'].append(dict(key=key, status='state_records_verified'))
                    consecutive = 0
            finally:
                stop_owned_children()
                gc.collect()
                from redox.qm.dft import _free_gpu_pool
                _free_gpu_pool()
                (claim/'owner.json').unlink(); claim.rmdir()
                atomic_write_json(path, report)
            # The collector is read-only with respect to scientific records.
            env = dict(os.environ, **{k: '1' for k in THREADS})
            result = subprocess.run([sys.executable, str(ROOT/'scripts/pipeline/collect_grafting_results.py')],
                                    env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            with (BASE/'collector.log').open('a') as log: log.write(result.stdout)
            if result.returncode: print('Collector needs review; records preserved', flush=True)
            if consecutive >= 3: raise RuntimeError('Three consecutive actual failures; worker stopped')
        report.update(status='no_unclaimed_unattempted_suitable_tasks', active_key=None)
    except BaseException as exc:
        report.update(status='stopped_preserving_outputs', error=repr(exc))
        raise
    finally:
        report['finished'] = time.time(); atomic_write_json(path, report)


def launch():
    from redox.workflows.resources import allocation_guard, cpu_allocation, free_slots, gpu_suitability, host_capacity
    plan = json.loads((BASE/'plan.json').read_text())
    verify(plan['producer_sha256']); verify(plan['input_sha256'])
    host = socket.gethostname().split('.')[0]
    preflight = json.loads((BASE/f'{host}_preflight.json').read_text())
    assert preflight['ok'] and preflight['python'] == sys.executable
    receipt = BASE/f'launch_{host}.json'
    assert not receipt.exists(), 'Inspect the existing launch before any resume'
    with allocation_guard():
        suitable = []
        for _, gpu in free_slots([host]):
            try: suitable.append((gpu, gpu_suitability(gpu)))
            except RuntimeError: continue
        if not suitable: raise RuntimeError('No free suitable GPUs')
        capacity = host_capacity(2*len(suitable), minimum_ram_gib=32*len(suitable))
        cpus = cpu_allocation(2*len(suitable))
        # Reserve through the shared gate directly. Each queue process is its
        # direct child, retaining ownership checks across sequential tasks.
        import site
        libraries = [str(p) for b in site.getsitepackages() for p in Path(b).glob('nvidia/*/lib')]
        env = dict(os.environ, PYTHONPATH=str(ROOT/'src'), REDOX_HESSIAN_MEMORY_GIB='24')
        env['PATH'] = str(Path(sys.executable).parent)+os.pathsep+env.get('PATH', '')
        env['LD_LIBRARY_PATH'] = ':'.join(libraries+[env.get('LD_LIBRARY_PATH', '')])
        env.update({k:'2' for k in THREADS})
        report = dict(host=host, status='submitting', created=time.time(), resources=capacity,
                      plan=str(BASE/'plan.json'), plan_sha256=digest(BASE/'plan.json'), workers=[])
        atomic_write_json(receipt, report)
        for i, (gpu, device) in enumerate(suitable):
            assigned = cpus[2*i:2*i+2]
            reservation = BASE/f'reserved_{host}_gpu{gpu}'
            reservation.mkdir()
            # The same gated dependency probe used by the existing workflow.
            argv = ['taskset', '--cpu-list', ','.join(map(str, assigned)),
                    '/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve', 'run', str(gpu),
                    '--label', 'redox-grafting', '--', sys.executable, '-m', 'redox.workflows.gpu',
                    '--out', str(reservation), '--probe', '--', sys.executable, '-u', str(SCRIPT), 'worker']
            log = BASE/f'controller_{host}_gpu{gpu}.log'
            with log.open('w') as stream:
                proc = subprocess.Popen(argv, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                    stdout=stream, stderr=subprocess.STDOUT, start_new_session=True,
                    preexec_fn=lambda assigned=assigned: os.sched_setaffinity(0, assigned))
            report['workers'].append(dict(gpu=gpu, device=device, cpus=assigned, gate_pid=proc.pid,
                process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
                argv=argv, log=str(log), submitted=time.time()))
            atomic_write_json(receipt, report)
            time.sleep(.5)
        report['status'] = 'submitted_actual_start_verification_pending'
        atomic_write_json(receipt, report)
        print(host, len(suitable), 'reserved queue workers submitted')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'worker', 'launch'])
    args = parser.parse_args()
    globals()[args.action]()
