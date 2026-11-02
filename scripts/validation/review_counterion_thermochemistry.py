#!/usr/bin/env python3
"""Reconstruct a complete chloride microstate component, keeping association pending."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
from ase.io import read as read_atoms

from redox.core.common import ROOT, HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify, load
from redox.validation.execution_review import review_closure, freeze_source, SOURCE_BYTES as CLOSURE_SOURCE
from redox.validation import primary_checkpoints
from redox.validation.scf_frame import freeze_reviewer, check_spin_contamination
from redox.workflows.ionic_stationary_cycles import screening_signature

SOURCE_BYTES = Path(__file__).read_bytes()


def review(namespace):
    spec = verify(namespace)
    case, sources = Path(spec['case']), Path(spec['sources'])
    used = [namespace]
    def read(path):
        used.append(Path(path))
        return json.loads(Path(path).read_text())
    plan = read(case/'plan.json')
    if (plan['role'] != 'primary_counterion_thermochemistry' or (plan['charge'], plan['mult']) != (1, 1)
            or plan['geometry_repair_permitted'] or plan['parent_formal_core_charge'] != 2):
        raise ValueError('Exact primary oxidized-core chloride microstate required')
    owners = load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py', 'chloride_thermo_owner')
    stage = review_closure(case, read, owners.alive,
        expected_status='counterion_thermochemistry_finished_independent_review_pending',
        host_GiB=40, GPU_MiB=31744, worker=sources/'worker.py')
    used.extend([freeze_reviewer(case), freeze_source(case, Path(__file__), SOURCE_BYTES),
        freeze_source(case, ROOT/'src/redox/validation/execution_review.py', CLOSURE_SOURCE),
        freeze_source(case, Path(primary_checkpoints.__file__), primary_checkpoints.SOURCE_BYTES)])
    geometry = case/'state/opt.xyz'
    if sha256_file(geometry) != plan['geometry_sha256'] or sha256_file(geometry) != sha256_file(sources/'reference_geometry.xyz'):
        raise ValueError('Native and checked geometry differ')
    used.extend([geometry, sources/'reference_geometry.xyz'])
    geometry_review = read(sources/'geometry_review.json')
    if any(sha256_file(p) != h for p, h in geometry_review['source_sha256'].items()):
        raise ValueError('Source geometry review backing changed')
    from redox.validation.protonation import hydrogen_topology_status
    if hydrogen_topology_status(plan['smiles'], geometry) != 'ok':
        raise ValueError('Exact disconnected chloride/covalent graft connectivity differs')
    energy = read(case/'electronic.json')
    protocol = plan['electronic_protocol']
    check_spin_contamination(energy.get('smd_spin_contam'))
    if (energy['energy_protocol'] != protocol or energy['optimized'] or not energy['converged_smd']
            or (energy['charge'], energy['mult']) != (1, 1)
            or abs(energy['e_smd_eV'] - HARTREE_EV*energy['e_smd_Ha']) > 1e-8):
        raise ValueError('Matched primary electronic energy context differs')
    electronic = primary_checkpoints.checkpoint_series(case/'electronic_SCF_evaluations', geometry,
        dict(charge=1, mult=1, xc=protocol['xc'], basis=protocol['basis'], disp=protocol['disp'],
             conv_tol=1e-9, max_cycle=200, solvent_class='SMD'), read, used,
        fixed_geometry=True, require_all_converged=False)
    if not electronic['converged_energies_Ha'] or abs(min(electronic['converged_energies_Ha']) - energy['e_smd_Ha']) > 1e-9:
        raise ValueError('Selected electronic energy is not the lowest converged checkpoint')
    result = read(case/'result.json')
    native_path = Path(result['native_record_path'])
    native = read(native_path)
    matrix_path = native_path.with_suffix('')/'hessian_cartesian.npy'
    matrix = np.load(matrix_path, allow_pickle=False)
    used.append(matrix_path)
    atoms = read_atoms(geometry)
    if (native['identity']['protocol'] != plan['native_protocol']
            or (native['identity']['charge'], native['identity']['mult']) != (1, 1)
            or native['status'] != 'derivative_checked_solution_thermal' or not native['stationary_minimum']
            or not native['converged'] or not native['grid_response'] or not 0 <= native['gradient_max'] <= .00045
            or matrix.shape != (3*len(atoms), 3*len(atoms)) or matrix.dtype != np.float64
            or not np.isfinite(matrix).all() or abs(matrix - matrix.T).max() > 1e-6):
        raise ValueError('Actual native state, stationarity or full Hessian gate failed')
    check_spin_contamination(native['smd_spin_diagnostics'].get('spin_contam'))
    series = primary_checkpoints.checkpoint_series(case/'SCF_evaluations', geometry,
        dict(charge=1, mult=1, basis='def2-svp', xc='r2scan', disp='d4', conv_tol=1e-12,
             conv_tol_grad=1e-8, max_cycle=200, grid_level=3, solvent_class='SMD'),
        read, used, expected_count=9, native_diagnostics=native['smd_spin_diagnostics'])
    fd = load(ROOT/'scripts/validation/review_molecular_native.py', 'chloride_native_FD')
    used.append(freeze_source(case, Path(fd.__file__)))
    checks, paths = fd.finite_difference_review(case, geometry, native, matrix)
    used.extend(paths)
    audit = load(ROOT/'scripts/validation/audit/native_partition.py', 'chloride_independent_partition')
    used.append(freeze_source(case, Path(audit.__file__)))
    partition = audit.read_partition(native_path, geometry)
    if partition is None or partition['protocol'] != plan['thermal_protocol']:
        raise ValueError('Independent modes, masses, symmetry, nuclear or spin partition failed')
    used.extend(Path(p) for p in partition['evidence_sha256'])
    projection = read(Path(partition['path']))
    values = {cutoff: energy['e_smd_eV'] + partition['nuclear_by_cutoff_eV'][cutoff]
              + projection['g_spin_partition_eV'] for cutoff in ('25','50','100')}
    if (any(not math.isfinite(v) or abs(v - result['G_by_cutoff_eV'][c]) > 1e-8 for c, v in values.items())
            or any(not math.isfinite(result['nuclear_by_cutoff_eV'][c])
                   or abs(result['nuclear_by_cutoff_eV'][c] - partition['nuclear_by_cutoff_eV'][c]) > 1e-8 for c in values)
            or result['electronic_eV'] != energy['e_smd_eV']):
        raise ValueError('Own electronic/nuclear/spin composite does not reconstruct')
    implementation = case/'screening_implementation/implementation.json'
    kernel = implementation.with_name('screening_kernel.cu')
    signature = screening_signature(read(implementation), sources, kernel)
    used.append(kernel)
    verify(namespace)
    host = read(case/'host_guard/status.json')
    return dict(status='primary_counterion_thermochemistry_independently_checked', namespace=str(namespace),
        key=plan['source_id'], smiles=plan['smiles'], charge=1, mult=1, formal_core_charge=2,
        dimensions=plan['dimensions'], geometry_sha256=plan['geometry_sha256'],
        electronic_protocol=protocol, native_protocol=plan['native_protocol'], thermal_protocol=plan['thermal_protocol'],
        temperature_K=298.15, standard_state=plan['thermal_protocol']['pressure'],
        reference='absolute composite microstate G within the primary molecular model; separated omitted iodide spectator',
        electronic_eV=energy['e_smd_eV'], nuclear_by_cutoff_eV=partition['nuclear_by_cutoff_eV'],
        spin_partition_eV=projection['g_spin_partition_eV'], G_by_cutoff_eV=values,
        minimum_frequency_cm=min(projection['frequencies_cm']), maximum_column_error=max(r['relative_error'] for r in checks),
        actual_electronic_SCF_frames_checked=electronic['actual_SCF_frames_checked'],
        actual_native_SCF_frames_checked=series['actual_SCF_frames_checked'], numerical_implementation=signature,
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'], peak_host_RSS_bytes=max(host['peak_tree_RSS_bytes'], host['self_VmHWM_bytes']),
        native_G_qualified=True, association_free_energy_computed=False,
        conformer_ensemble_qualified=False, experimental_populations_inferred=False, material_admission=False,
        source_sha256={str(p):sha256_file(p) for p in used}, producer_sha256=hashlib.sha256(SOURCE_BYTES).hexdigest())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spec', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    if a.out.exists():
        raise FileExistsError('Preserve earlier independent microstate reviews')
    result = review(a.spec.resolve()); atomic_write_json(a.out, result)
    print(json.dumps({k:result[k] for k in ('status','key','minimum_frequency_cm','maximum_column_error')}))
