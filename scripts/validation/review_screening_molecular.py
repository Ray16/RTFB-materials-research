#!/usr/bin/env python3
"""Recompute full finite-difference checks from closed molecular control records."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def review(namespace):
    spec = verify(namespace)
    case = Path(spec['case'])
    used = [namespace]
    def read(path):
        used.append(path)
        return json.loads(path.read_text())
    plan = read(case / 'plan.json')
    model = read(Path(spec['sources']) / 'molecule.json')
    state = model['states'][plan['molecular_control']]
    atoms = state.get('atoms', model['atoms'])
    coordinates = np.asarray([a[1] for a in atoms], dtype=np.float64)
    stage = read(case / 'status.json')
    launch = read(case.parent / 'placement/launch.json')
    completion = read(case.parent / 'placement/completion.json')
    host = read(case / 'host_guard/status.json')
    runtime = read(case / 'tensor_runtime/status.json')
    owner = read(case / 'residency_guard/owner.json')
    residency = read(case / 'residency_guard/state.json')
    definition = importlib.util.spec_from_file_location('molecular_control_owner', ROOT / 'scripts/pipeline/observe_archived_native_namespace.py')
    checker = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(checker)
    if (stage['status'] != 'screening_molecular_control_finished_review_pending'
            or stage['native_G_qualified'] or checker.alive(stage['owner']) is not False
            or completion['returncode'] != 0 or completion['budget_stop'] is not None
            or host['status'] != 'effective_memory_guard_finished_profile_review_pending'
            or runtime['status'] != 'task_returned_under_checked_cutensor_runtime'
            or not owner['guard_stopped'] or not residency['own_residency_observed']
            or residency['foreign_resident_PIDs'] or host['free_RAM_bytes'] < 32 * 2**30
            or host['maximum_RSS_GiB'] != plan['maximum_host_RSS_GiB']
            or max(host['peak_tree_RSS_bytes'], host['self_VmHWM_bytes']) > plan['maximum_host_RSS_GiB'] * 2**30
            or not stage['owner']['pid'] == host['pid'] == owner['pid'] == runtime['pid']
            or not stage['owner']['process_start'] == host['process_start'] == owner['process_start']
            or not host['CPU_affinity'] == owner['CPU_affinity'] == runtime['CPU_affinity'] == launch['CPU_affinity']
            or len(host['CPU_affinity']) != 2):
        raise ValueError('Closed molecular control and normal resource/ownership checks required')
    values = {}
    checks = []
    for variant in ('installed', 'corrected'):
        folder = case / variant
        points = [read(p) for p in sorted(folder.glob('point_*/evaluation.json'))]
        expected_points = 1 + 2 * coordinates.size * len(model['FD_steps_Bohr'])
        if len(points) != expected_points:
            raise ValueError('All Cartesian columns at both steps must be present')
        for point in points:
            if (not point['SCF_converged'] or point['variant'] != variant
                    or point['native_G_qualified'] or abs(point['spin']['spin_contam']) > .1
                    or not np.isfinite(point['energy_Ha'])
                    or np.asarray(point['gradient_Ha_Bohr']).shape != coordinates.shape
                    or not np.isfinite(point['gradient_Ha_Bohr']).all()):
                raise ValueError('Invalid independent SCF record')
        center = points[0]
        if not np.array_equal(center['coordinates_Bohr'], coordinates):
            raise ValueError('Center differs from fixed prospective control geometry')
        hessian_path = folder / 'Hessian.npy'
        used.append(hessian_path)
        matrix = np.load(hessian_path, allow_pickle=False)
        if matrix.shape != (coordinates.size, coordinates.size) or not np.isfinite(matrix).all():
            raise ValueError('Invalid analytical Hessian shape/values')
        asymmetry = float(abs(matrix - matrix.T).max())
        if asymmetry > model['Hessian_asymmetry_limit_Ha_Bohr2']:
            raise ValueError('Unchanged Hessian symmetry gate failed')
        report = read(folder / 'derivatives.json')
        if not report['complete'] or len(report['checks']) != coordinates.size * len(model['FD_steps_Bohr']):
            raise ValueError('Incomplete producer finite-difference panel')
        observed = set()
        for column in report['checks']:
            atom, axis, step = column['atom'], column['axis'], column['step_Bohr']
            key = (atom, axis, step)
            if key in observed or atom not in range(len(atoms)) or axis not in range(3) or step not in model['FD_steps_Bohr']:
                raise ValueError('Duplicate or undeclared finite-difference coordinate')
            observed.add(key)
            displaced = []
            for sign in (-1, 1):
                target = coordinates.copy()
                target[atom, axis] += sign * step
                matches = [p for p in points if np.array_equal(p['coordinates_Bohr'], target)]
                if len(matches) != 1:
                    raise ValueError('One independent exact displaced SCF required')
                displaced.append(matches[0])
            finite = (np.asarray(displaced[1]['gradient_Ha_Bohr']) - np.asarray(displaced[0]['gradient_Ha_Bohr'])).ravel() / (2 * step)
            error = float(np.linalg.norm(finite - matrix[:, 3 * atom + axis]) / max(np.linalg.norm(finite), 1e-14))
            energy_gradient = (displaced[1]['energy_Ha'] - displaced[0]['energy_Ha']) / (2 * step)
            gradient_error = energy_gradient - center['gradient_Ha_Bohr'][atom][axis]
            if (abs(error - column['relative_Hessian_column_error']) > 1e-14
                    or abs(gradient_error - column['energy_gradient_error_Ha_Bohr']) > 1e-14
                    or error > model['FD_relative_limit']
                    or abs(gradient_error) > model['energy_gradient_absolute_limit_Ha_Bohr']):
                raise ValueError('Reconstructed energy/gradient/Hessian derivative gate failed')
            checks.append(dict(variant=variant, atom=atom, axis=axis, step_Bohr=step,
                               relative_Hessian_error=error, energy_gradient_error_Ha_Bohr=gradient_error))
        values[variant] = dict(energy=center['energy_Ha'], gradient=np.asarray(center['gradient_Ha_Bohr']), Hessian=matrix)
    a, b = [values[v] for v in ('installed', 'corrected')]
    differences = dict(energy_Ha=abs(b['energy'] - a['energy']),
        gradient_Ha_Bohr=float(abs(b['gradient'] - a['gradient']).max()),
        Hessian_Ha_Bohr2=float(abs(b['Hessian'] - a['Hessian']).max()))
    if any(differences[k] > limit for k, limit in model['installed_corrected_limits'].items()):
        raise ValueError('Installed/corrected small-control agreement failed')
    verify(namespace)
    return dict(status='small_molecule_screening_energy_force_Hessian_control_independently_checked',
        namespace=str(namespace), control=plan['molecular_control'], charge=state['charge'], mult=state['mult'],
        device=launch['device'], atoms=len(atoms), derivative_checks=checks, differences=differences,
        maximum_relative_Hessian_error=max(c['relative_Hessian_error'] for c in checks),
        maximum_energy_gradient_error_Ha_Bohr=max(abs(c['energy_gradient_error_Ha_Bohr']) for c in checks),
        polymer_implementation_qualified=False, native_G_qualified=False, production_default_changed=False,
        source_sha256={str(p): sha256_file(p) for p in used}, producer_sha256=sha256_file(__file__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=2)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise FileExistsError('Preserve earlier molecular control reviews')
    result = review(args.spec.resolve())
    atomic_write_json(args.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'control', 'maximum_relative_Hessian_error', 'differences')}))
