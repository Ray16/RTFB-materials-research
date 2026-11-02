#!/usr/bin/env python3
"""Review the actual closed gas/SMD SCF determinants of one corrected input."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np
from ase.io import read as read_atoms
from pyscf.lib import param

from redox.core.common import ROOT, HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.workflows.molecular_native_namespace import verify

SOURCE_BYTES = Path(__file__).read_bytes()


def load(path, name):
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition); definition.loader.exec_module(module)
    return module


def review(namespace):
    from redox.validation.execution_review import review_closure, freeze_source, SOURCE_BYTES as closure_bytes
    from redox.validation.scf_frame import checkpoint_frame, freeze_reviewer
    from redox.workflows.sampling_energy import profile
    spec = verify(namespace); case = Path(spec['case']); sources = Path(spec['sources']); used = [namespace]
    def read(path):
        used.append(path); return json.loads(path.read_text())
    plan = read(case / 'plan.json'); task = plan['tasks'][0]; declared = profile(task)
    if (plan['role'] != 'corrected_sampling_energy' or plan['electronic_protocol'] != ACTIVE_SP
            or any(plan.get(k) != v for k, v in declared.items()) or plan['native_G_computed']
            or plan['primary_candidate_admission'] or plan['historical_energy_or_thermal_values_transferred']):
        raise ValueError('Exact corrected same-model electronic-only scope required')
    owners = load(ROOT / 'scripts/pipeline/observe_archived_native_namespace.py', 'sampling_energy_ownership')
    stage = review_closure(case, read, owners.alive,
        expected_status='corrected_sampling_energy_finished_independent_review_pending',
        host_GiB=plan['maximum_host_RSS_GiB'], GPU_MiB=31744, worker=sources / 'worker.py')
    used += [freeze_source(case, Path(__file__), SOURCE_BYTES), freeze_reviewer(case),
             freeze_source(case, ROOT / 'src/redox/validation/execution_review.py', closure_bytes)]
    original = read(case / 'experimental_energy/status.json')
    if (original['status'] != 'experimental_hybrid_energy_finished_independent_review_pending'
            or original['owner'] != stage['owner']
            or original['observed_peak_GPU_MiB'] != stage['observed_peak_GPU_MiB']):
        raise ValueError('Original and corrected execution ownership/resources differ')
    for path, digest in original['output_sha256'].items():
        if sha256_file(path) != digest:
            raise ValueError('Raw original electronic output changed')
        used.append(Path(path))
    raw = Path(plan['raw_root']); result = read(raw / 'experimental_energy.json')
    driver = load(Path(spec['archive']) / 'scripts/pipeline/smooth_polymer_energy.py', 'sampling_energy_numerical_gates')
    driver.check_energy(result, task, plan['experimental_method_id'])
    if result['energy_protocol'] != plan['electronic_protocol']:
        raise ValueError('Actual hybrid electronic protocol differs from its declaration')
    atoms = read_atoms(case / 'input.xyz'); coordinates = atoms.positions / param.BOHR
    saved = load(sources / 'occupied.py', 'sampling_energy_frozen_occupied_format')
    groups = {'gas': [], 'smd': []}; counts = {'gas': 0, 'smd': 0}; errors = []; orthogonal = []
    for checkpoint in sorted((raw / 'SCF_checkpoints').glob('*/final/occupied.h5')):
        folder = checkpoint.parent.parent
        attempts = sorted(p for p in folder.glob('*/settings.json') if p.parent.name != 'final')
        if not attempts:
            raise ValueError('Final occupied determinant lacks its actual solver settings')
        settings = read(attempts[-1]); solvent = settings['solvent']
        name = 'gas' if solvent is None else 'smd' if solvent == 'acetonitrile' else None
        if (name is None or settings['charge'] != 0 or settings['mult'] != 1
                or settings['xc'] != ACTIVE_SP['xc'] or settings['disp'] != ACTIVE_SP['disp']
                or settings['nlc'] != ACTIVE_SP['nlc'] or settings['conv_tol'] != 1e-9
                or settings['max_cycle'] != 200 or settings['grid_level'] != 3):
            raise ValueError('Actual solvent/SCF/functional settings differ')
        frame = checkpoint_frame(checkpoint, charge=0, mult=1, basis=ACTIVE_SP['basis'],
            coordinates_Bohr=coordinates, symbols=atoms.get_chemical_symbols())
        proof = saved.review(checkpoint); used.append(checkpoint)
        if proof['basis_functions'] != task['active_basis_functions']:
            raise ValueError('Actual checkpoint basis size differs from the declared graph')
        counts[name] += 1; errors.append(frame['maximum_coordinate_error_Bohr'])
        orthogonal.append(proof['maximum_orthonormality_error'])
        if proof['converged']:
            groups[name].append(proof['energy_Ha'])
    for name, energies in groups.items():
        if (not energies or abs(min(energies) - result['e_' + name + '_eV'] / HARTREE_EV) > 1e-9):
            raise ValueError('Lowest converged gas/SMD energy lacks actual occupied-checkpoint evidence')
    implementation = read(case / 'screening_implementation/implementation.json')
    if (implementation['status'] != 'run_local_screening_context_closed_qualification_separate'
            or implementation['installed_library_files_modified'] or implementation['physics_model_changed']):
        raise ValueError('Closed run-local corrected screening required')
    used += list((case / 'screening_implementation').glob('*.py'))
    used += list((case / 'screening_implementation').glob('*.cu'))
    verify(namespace)
    host = read(case / 'host_guard/status.json')
    return dict(status='corrected_sampling_energy_independently_checked', key=task['key'],
        namespace=str(namespace), task=task, portable_method_id=plan['experimental_method_id'],
        experimental_model=plan['experimental_model'], electronic_protocol=plan['electronic_protocol'],
        E_smd_eV=result['e_smd_eV'], E_gas_eV=result['e_gas_eV'],
        actual_checkpoint_counts=counts, actual_checkpoint_frames_verified=True,
        maximum_checkpoint_coordinate_error_Bohr=max(errors), occupied_wavefunctions_checked=True,
        maximum_occupied_orthonormality_error=max(orthogonal),
        virtual_coefficients_retained=False, response_restart_supported=False,
        resources_checked=True, peak_GPU_MiB=stage['observed_peak_GPU_MiB'],
        peak_host_RSS_bytes=max(host['self_VmHWM_bytes'], host['peak_tree_RSS_bytes']),
        native_G_qualified=False, primary_candidate_admission=False,
        source_sha256={str(p): sha256_file(p) for p in used},
        reviewer_sha256=hashlib.sha256(SOURCE_BYTES).hexdigest())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--namespace', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    if a.out.exists():
        raise FileExistsError('Preserve earlier independent energy reviews')
    result = review(a.namespace.resolve()); atomic_write_json(a.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'key', 'actual_checkpoint_counts', 'peak_GPU_MiB')}))
