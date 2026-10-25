#!/usr/bin/env python
"""Search away from legacy-CDS derivative-scale artifacts, retaining the energy.

Broader finite-difference CDS directions are used only inside optimization.
Final forces use the original adaptive derivative, original total-energy tests
and a stronger all-axis CDS scale-consistency check. No gate is relaxed.
"""
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

GEOMETRY = ROOT / 'scripts/pipeline/optimize_polymer_geometry.py'
GC = ROOT / 'scripts/pipeline/smd_gc_workspace.py'


def provisional(row):
    result = dict(row)
    if result.get('status') == 'stationary_geometry_checked_active_parent_and_native_minimum_pending':
        result['status'] = 'base_geometry_checks_passed_extended_CDS_checks_pending'
    return result


def prepare(source, out):
    old = api.read_plan(source)
    api.verify_sources(old['producer_sha256']); api.verify_sources(old['input_sha256'])
    status_path = source / 'geometry/status.json'; failed = json.loads(status_path.read_text())
    alive = api.load_module('force_repair_source_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py').owner_alive
    if (failed['status'] != 'geometry_stage_failed_preserved' or 'force fails energy differences' not in failed.get('error', '')
            or alive(failed['owner']) or len(old['tasks']) != 1):
        raise ValueError('A finished exact geometry with failed independent force validation is required')
    state = Path(failed['state_dir']); result_path = state / 'result.json'; result = json.loads(result_path.read_text())
    if not result.get('geometry_optimization_converged') or not result.get('optimization_SCF_converged'):
        raise ValueError('The source must contain a converged optimizer geometry')
    seed = state / 'geometry_checkpoints/optimizer_returned.xyz'
    if not seed.exists():
        seed = state / 'opt.xyz'
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(source / 'controls', out / 'controls')
    shutil.copytree(source / 'live_inputs', out / 'live_inputs')
    task = old['tasks'][0]; ready_path = out / 'live_inputs' / task['key'] / 'task.json'; xyz = ready_path.parent / 'input.xyz'
    shutil.copyfile(seed, xyz)
    audit = api.load_module('force_repair_input_graph_audit', ROOT / 'scripts/validation/audit_polymer_geometries.py').audit(task, xyz)
    if audit['issues']:
        raise ValueError('Saved optimizer geometry has a graph/stereo issue')
    atomic_write_json(ready_path.parent / 'recovery_start_review.json', audit)
    ready = json.loads(ready_path.read_text()); ready['input_xyz'] = str(xyz); ready['input_sha256'][str(xyz)] = api.digest(xyz)
    ready['source_restart_checkpoint'] = dict(source=str(seed), source_qualified=False)
    atomic_write_json(ready_path, ready)
    for stage in ('energy', 'thermal'):
        for name in ('tasks', 'done', 'claims', 'failures', 'workers', 'pilot', 'launches'):
            (out / stage / name).mkdir(parents=True)
    inputs = dict(old['input_sha256'])
    for p in (source / 'plan.json', status_path, result_path, seed, state / 'total_energy_force_check/summary.json', xyz, ready_path):
        inputs[str(p)] = api.digest(p)
    producers = dict(old['producer_sha256'])
    for p in (Path(__file__).resolve(), GEOMETRY, GC):
        producers[str(p)] = api.digest(p)
    plan = dict(old, raw_root=str(ROOT / 'calcs/polymer_models' / out.parent.name / out.name),
        input_sha256=inputs, producer_sha256=producers, recovery_from=str(source),
        optimizer_direction_settings=dict(CDS_step_Bohr=.0005, energy_function_changed=False,
            force_used_for_final_admission=False, convergence_gmax=1e-5, convergence_grms=6e-6,
            convergence_energy=1e-8, convergence_dmax=.0005, convergence_drms=.0003),
        extended_CDS_gate=dict(steps_Bohr=[1e-4,5e-4,1e-3], maximum_gradient_scale_difference_Ha_Bohr=1e-6),
        purpose='Same-energy geometry repair; independent adaptive forces and original total-energy gates plus full-axis CDS consistency')
    atomic_write_json(out / 'plan.json', plan)


@contextmanager
def search_directions(path):
    import numpy as np
    from gpu4pyscf.solvent.grad import smd as grad
    from gpu4pyscf.solvent import smd
    from redox.qm.cds_finite import energy_gradient
    previous = grad.get_cds; counter = 0
    def evaluate(solvent):
        nonlocal counter
        mol, copy = solvent.mol.copy(), solvent.copy()
        positions = mol.atom_coords(unit='Bohr')
        def energy(coords):
            mol.set_geom_(coords, unit='Bohr'); copy.reset(mol)
            return smd.get_cds_legacy(copy)[0]
        value = energy_gradient(positions, energy, .0005)
        if not np.isfinite(value).all():
            raise ValueError('Nonfinite CDS search direction')
        counter += 1
        atomic_write_json(path, dict(evaluations=counter, step_Bohr=.0005,
            energy_function_changed=False, used_for_physical_admission=False, updated=time.time()))
        return value
    grad.get_cds = evaluate
    try:
        yield
    finally:
        grad.get_cds = previous


def run(campaign):
    import numpy as np
    from ase.io import read
    from pyscf import gto
    from gpu4pyscf.solvent.smd import get_cds_legacy
    from types import SimpleNamespace
    from redox.qm import dft
    from redox.qm.cds_finite import energy_gradient
    from redox.workflows.resources import host_capacity
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    guard = json.loads((campaign / 'host_guard/status.json').read_text())
    if guard['pid'] != os.getpid() or guard['maximum_RSS_GiB'] != 16 or guard['minimum_free_GiB'] != 32:
        raise ValueError('The geometry-only same-PID host guard is required')
    module = api.load_module('energy_consistent_geometry_original_physical_gates', GEOMETRY)
    module.host_capacity = lambda threads, **kwargs: host_capacity(threads, minimum_ram_gib=48)
    original_settings = api.optimizer_settings
    @contextmanager
    def settings(*args, **kwargs):
        with original_settings(*args, **kwargs):
            previous = dft.dft_smd
            def tight(*a, **kw):
                values = dict(kw.get('conv_params') or {})
                values.update({k:v for k,v in plan['optimizer_direction_settings'].items() if k.startswith('convergence_')})
                kw['conv_params'] = values
                return previous(*a, **kw)
            dft.dft_smd = tight
            try:
                with search_directions(campaign / 'optimizer_direction_progress.json'):
                    yield
            finally:
                dft.dft_smd = previous
    original_writer = module.atomic_write_json
    def write(path, row):
        original_writer(path, provisional(row) if Path(path) == campaign / 'geometry/status.json' else row)
    module.atomic_write_json = write; api.optimizer_settings = settings
    gc = api.load_module('energy_consistent_geometry_GC', GC)
    metrics = dict(collections=0, unreachable_objects_collected=0)
    try:
        with gc.collect_before_build(metrics):
            module.run(campaign)
        status_path = campaign / 'geometry/status.json'; record = json.loads(status_path.read_text())
        state = Path(record['state_dir']); baseline = json.loads((state / 'total_energy_force_check/baseline.json').read_text())
        coords = np.array(baseline['coordinates_Bohr']); symbols = read(state / 'opt.xyz').get_chemical_symbols()
        mol = gto.M(atom=list(zip(symbols,coords)), unit='Bohr', basis='def2-svp', charge=0, spin=0, verbose=0)
        obj = SimpleNamespace(mol=mol, solvent_descriptors=None, solvent='acetonitrile')
        def energy(x):
            obj.mol.set_geom_(x, unit='Bohr'); return get_cds_legacy(obj)[0]
        values = {str(h): energy_gradient(coords,energy,h) for h in plan['extended_CDS_gate']['steps_Bohr']}
        maximum = max(float(np.max(abs(values['0.0001'] - value))) for value in values.values())
        proof = campaign / 'CDS_full_axis_consistency.json'
        atomic_write_json(proof, dict(maximum_scale_difference_Ha_Bohr=maximum,
            steps_Bohr=plan['extended_CDS_gate']['steps_Bohr'], all_Cartesian_components_checked=True,
            passed=maximum<=1e-6, source_geometry_sha256=api.digest(state/'opt.xyz')))
        if maximum > 1e-6:
            raise ValueError('Legacy CDS is not derivative-scale consistent at the returned geometry')
        record['status'] = 'stationary_geometry_checked_active_parent_and_native_minimum_pending'
        record['extended_CDS_consistency'] = json.loads(proof.read_text())
        record['output_sha256'][str(proof)] = api.digest(proof)
        record['finished'] = time.time()
        atomic_write_json(status_path, record)
    except BaseException as error:
        path = campaign / 'geometry/status.json'
        if path.exists():
            record = json.loads(path.read_text()); record.update(status='geometry_repair_failed_preserved', error=repr(error), finished=time.time()); atomic_write_json(path,record)
        raise
    finally:
        api.optimizer_settings = original_settings
        atomic_write_json(campaign / 'geometry_GC.json', metrics)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','run'))
    parser.add_argument('--source', type=Path)
    parser.add_argument('--campaign', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        if not args.source: parser.error('--source required')
        prepare(args.source.resolve(),args.campaign.resolve())
    else:
        run(args.campaign.resolve())


if __name__ == '__main__':
    main()
