#!/usr/bin/env python3
"""CPU reconstruction of sealed polymer native records after checkout relocation.

No SCF, optimizer or GPU computation runs. The original completed records and
source hashes remain immutable. This replays numerical/physical evidence; it
does not supply missing conformers or qualify the current production runtime.
"""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import re

import numpy as np
from ase.io import read

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources


def configuration(task):
    names = [task[k] for k in ('configuration_label', 'stereoisomer') if k in task]
    if not names or any(name != names[0] for name in names):
        raise ValueError('An explicit, consistent backbone configuration is required')
    return names[0]


def spin_valid(point):
    value = point['spin'].get('spin_contam', math.inf)
    return type(value) in (float, int) and math.isfinite(value) and abs(value) <= .1


def projected_modes(matrix, xyz, mass):
    """Independent rigid-space removal used in the original scientific audit."""
    centered = xyz - np.average(xyz, axis=0, weights=mass)
    vectors = [(np.sqrt(mass)[:, None]*axis).ravel() for axis in np.eye(3)]
    vectors += [(np.sqrt(mass)[:, None]*np.cross(centered, axis)).ravel() for axis in np.eye(3)]
    rigid = []
    for vector in vectors:
        residual = vector.copy()
        for _ in range(2):
            for unit in rigid:
                residual -= unit*np.dot(unit, residual)
        if np.linalg.norm(residual) > 1e-10*max(1., np.linalg.norm(vector)):
            rigid.append(residual/np.linalg.norm(residual))
    basis = np.array(rigid).T
    values, vectors = np.linalg.eigh(np.eye(matrix.shape[0])-basis@basis.T)
    internal = vectors[:, values > .5]
    roots = np.repeat(np.sqrt(mass*(1.660539040e-27/9.10938356e-31)), 3)
    weighted = (matrix+matrix.T)/(2*np.outer(roots, roots))
    values = np.linalg.eigvalsh(internal.T@weighted@internal)
    return np.sign(values)*np.sqrt(abs(values))*219474.63067, len(rigid)


def finite_difference_checks(raw, native, matrix, gradient, symbols, coords, method, read_json):
    hetero = next((i for i, symbol in enumerate(symbols) if symbol == 'O'),
                  next((i for i, symbol in enumerate(symbols) if symbol == 'Fe'), len(symbols)-1))
    required = {(0, 1), (hetero, 2)}
    required.update(tuple(map(int, np.unravel_index(i, gradient.shape)))
                    for i in np.argsort(abs(gradient).ravel())[-2:])
    if {tuple(v) for v in native['required_columns']} != required:
        raise ValueError('Hessian-column coverage differs')
    baseline = read_json(raw/'SCF_evaluations/0000/evaluation.json')
    if (not baseline['SCF_converged'] or not spin_valid(baseline)
            or baseline['experimental_method_id'] != method
            or not np.allclose(baseline['coordinates_Bohr'], coords, rtol=0, atol=1e-10)
            or not np.allclose(baseline['gradient_Ha_Bohr'], gradient, rtol=0, atol=1e-13)):
        raise ValueError('Baseline SCF/gradient/geometry binding differs')
    index, checks = 1, []
    for atom, axis in sorted(required):
        for step in (.0005, .001):
            points = [read_json(raw/f'SCF_evaluations/{i:04d}/evaluation.json') for i in (index, index+1)]
            index += 2
            vectors = []
            for sign, point in zip((-1, 1), points):
                expected = coords.copy()
                expected[atom, axis] += sign*step
                vector = np.asarray(point['gradient_Ha_Bohr'], dtype=float)
                if (vector.shape != gradient.shape or not np.isfinite(vector).all()
                        or not np.allclose(point['coordinates_Bohr'], expected, rtol=0, atol=1e-12)
                        or not point['SCF_converged'] or not spin_valid(point)
                        or point['experimental_method_id'] != method):
                    raise ValueError('Displaced SCF/gradient/model binding differs')
                vectors.append(vector)
            finite = (vectors[1]-vectors[0]).ravel()/(2*step)
            error = float(np.linalg.norm(finite-matrix[:, 3*atom+axis])/max(np.linalg.norm(finite), 1e-14))
            stored = next(v for v in native['derivative_checks']
                          if (v['atom'], v['axis'], v['step_Bohr']) == (atom, axis, step))
            if (not math.isfinite(error) or not math.isfinite(stored['relative_error'])
                    or error > .01 or abs(error-stored['relative_error']) > 1e-10):
                raise ValueError('Sampled Hessian derivative check differs or fails')
            checks.append(error)
    return checks, required, index


def replay(item, sources, nuclear_equations):
    review_path = sources.verify(item['review'], item['review_sha256'])
    report = json.loads(review_path.read_text())
    if (report['status'] != 'alternate_native_state_independently_reconstructed'
            or report['primary_candidate_admission'] is not False
            or report['task']['key'] != item['key']
            or (report['task']['charge'], report['task']['mult']) != (0, 1)):
        raise ValueError('Only sealed checked alternate-model native states may replay')
    bindings = report['source_sha256']
    sources.verify_all(bindings)

    def bound(path):
        path = str(path)
        if path not in bindings:
            raise ValueError('Unbound raw evidence: '+path)
        return sources.verify(path, bindings[path])

    def read_json(path):
        return json.loads(bound(path).read_text())

    campaign = Path(report['campaign'])
    plan = read_json(campaign/'plan.json')
    state = read_json(campaign/'experimental_native/status.json')
    if (state['status'] != 'experimental_native_finished_independent_review_pending'
            or plan['tasks'] != [report['task']]
            or plan['experimental_method_id'] != report['experimental_method_id']
            or plan['experimental_model'] != report['experimental_model']
            or plan['experimental_model']['solvent'] != 'acetonitrile'):
        raise ValueError('Completed state, exact graph or solvent method differs')
    sources.verify_all(plan['producer_sha256'])
    sources.verify_all(plan['input_sha256'])
    sources.verify_all(state['output_sha256'])
    raw = Path(state['raw_root'])
    geometry = bound(raw/'geometry.xyz')
    if sha256_file(geometry) != plan['scored_geometry_sha256']:
        raise ValueError('Native and scored electronic geometries differ')
    native = read_json(raw/'experimental_partition.json')
    parent = read_json(Path(plan['electronic_review']))
    sources.verify_all(parent['source_sha256'])
    if (not math.isfinite(parent['E_smd_eV']) or abs(native['spin_partition_eV']) > 1e-12
            or not math.isfinite(native['spin_partition_eV'])):
        raise ValueError('Finite neutral-singlet electronic and spin terms required')
    if (native['experimental_method_id'] != report['experimental_method_id']
            or native['primary_candidate_admission'] is not False
            or not native['derivative_checks_passed']
            or native['symmetry_source'] != 'xTB autodetection on exact geometry'
            or native['native_settings'] != plan['experimental_native_settings']):
        raise ValueError('Native scope, derivative, symmetry or settings differ')
    matrix = np.load(bound(raw/'hessian.npy'), allow_pickle=False)
    gradient = np.load(bound(raw/'gradient.npy'), allow_pickle=False)
    atoms = read(geometry)
    count = len(atoms)
    if (matrix.shape != (3*count, 3*count) or gradient.shape != (count, 3)
            or not np.isfinite(matrix).all() or not np.isfinite(gradient).all()
            or abs(matrix-matrix.T).max() > 1e-6 or abs(gradient).max() > .00045):
        raise ValueError('Hessian or geometry stationarity fails')
    checks, required, index = finite_difference_checks(
        raw, native, matrix, gradient, atoms.get_chemical_symbols(),
        atoms.positions/.529177210903, report['experimental_method_id'], read_json)
    if len(list(sources.path(raw/'SCF_evaluations').glob('*/evaluation.json'))) != index:
        raise ValueError('Displaced-gradient inventory differs')
    nuclear = native['nuclear']
    masses = np.asarray(nuclear['atomic_masses_amu'], dtype=float)
    if masses.shape != (count,) or not np.isfinite(masses).all() or min(masses) <= 0:
        raise ValueError('Finite positive atom-matched masses required')
    frequencies, rank = projected_modes(matrix, atoms.positions, masses)
    error = float(np.max(abs(frequencies-np.asarray(nuclear['frequencies_cm']))))
    if (rank != 6 or len(frequencies) != 3*count-6 or not math.isfinite(error)
            or not np.isfinite(frequencies).all() or error > 1e-5 or min(frequencies) < -5):
        raise ValueError('Internal frequencies or minimum gate fails')
    sigma = int(re.search(r'rotational number\s+(\d+)', bound(raw/'symmetry_thermo.log').read_text())[1])
    if sigma < 1 or sigma != native['symmetry_number']:
        raise ValueError('Symmetry identity differs')
    errors, free = {}, {}
    for rotor in ('25', '50', '100'):
        value = float(nuclear_equations(frequencies, atoms.positions, masses, sigma, int(rotor)))
        errors[rotor] = abs(value-nuclear['nuclear_by_cutoff_eV'][rotor])
        free[rotor] = float(parent['E_smd_eV']+value+native['spin_partition_eV'])
        if (not all(math.isfinite(v) for v in (value, errors[rotor], native['G_by_cutoff_eV'][rotor],
                                               report['G_by_cutoff_eV'][rotor])) or errors[rotor] > 1e-8
                or abs(free[rotor]-native['G_by_cutoff_eV'][rotor]) > 1e-8
                or abs(free[rotor]-report['G_by_cutoff_eV'][rotor]) > 1e-8):
            raise ValueError('Independent nuclear or free-energy sum differs')
    guard = campaign/'residency_guard'
    ownership = read_json(guard/'owner.json')
    residency = read_json(guard/'state.json')
    host = read_json(campaign/'host_guard/status.json')
    if (not ownership['guard_stopped'] or not residency['own_residency_observed']
            or residency['foreign_resident_PIDs'] or sources.path(guard/'interruption.json').exists()
            or host['status'] != 'program_returned_review_required'
            or not 0 < host['peak_RSS_bytes'] <= 80*2**30
            or not 0 < state['observed_peak_GPU_MiB'] <= 31744):
        raise ValueError('Sealed execution-resource or exclusive-residency checks fail')
    return dict(key=item['key'], configuration=configuration(report['task']),
                occupied_sites=report['task']['occupied_sites'], graph_sha256=report['task']['graph_sha256'],
                backbone_stereo_fingerprint=report['task']['backbone_stereo_fingerprint'],
                method_id=report['experimental_method_id'], G_by_cutoff_eV=free,
                E_eV=float(parent['E_smd_eV']), review=item['review'],
                review_sha256=item['review_sha256'], atoms=count,
                sampled_columns=len(required), maximum_column_error=max(checks),
                minimum_frequency_cm=float(min(frequencies)), frequency_error_cm=error,
                nuclear_errors_eV=errors, status='relocated_native_physical_replay_passed')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('archive', 'prior-root', 'manifest', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        raise FileExistsError('Preserve previous replay evidence')
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    manifest = json.loads(a.manifest.read_text())
    relative_base = manifest.get('relative_paths_base')
    if relative_base is not None and relative_base != str(a.prior_root):
        raise ValueError('Relative source paths require the exact declared prior checkout directory')
    reader = ArchivedSources(a.prior_root, a.archive, relative_to_prior=relative_base is not None)
    producer_paths = [Path(__file__), ROOT/'src/redox/validation/archived_sources.py',
                      ROOT/'src/redox/core/hashing.py', ROOT/'src/redox/core/io.py']
    producer_hashes = {str(path.resolve()): sha256_file(path) for path in producer_paths}
    manifest_hash = sha256_file(a.manifest)
    equations_path = ROOT/'scripts/validation/audit/native_partition.py'
    original_equations = str(a.prior_root/'scripts/validation/audit/native_partition.py')
    if sha256_file(equations_path) != sha256_file(reader.path(original_equations)):
        raise ValueError('Maintained independent nuclear equations differ from the archived source')
    spec = importlib.util.spec_from_file_location('archived_nuclear_equations', equations_path)
    equations = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(equations)
    keys = [item['key'] for item in manifest['reviews']]
    if not keys or len(keys) != len(set(keys)):
        raise ValueError('Nonempty distinct conformer keys required')
    rows = []
    for item in manifest['reviews']:
        rows.append(replay(item, reader, equations.nuclear_equations))
        atomic_write_json(a.out.with_name(a.out.stem+'_progress.json'),
                          dict(status='partial_replay_final_source_recheck_pending',
                               completed_states=len(rows), expected_states=len(keys), records=rows,
                               producer_sha256=producer_hashes, manifest_sha256=manifest_hash))
        print(json.dumps(dict(checked=len(rows), key=item['key'])), flush=True)
    reader.recheck()
    if (sha256_file(a.manifest) != manifest_hash
            or any(sha256_file(path) != value for path, value in producer_hashes.items())):
        raise ValueError('Replay implementation or declared inputs changed during execution')
    result = dict(status='archived_polymer_native_physical_records_replayed', states=len(rows), records=rows,
                  prior_root=str(a.prior_root), archive_root=str(a.archive.resolve()),
                  sources_verified=len(reader.receipts), relocated_sources=reader.receipts,
                  manifest_sha256=manifest_hash, producer_sha256=producer_hashes,
                  equations_sha256=sha256_file(equations_path), CPU_affinity=sorted(os.sched_getaffinity(0)),
                  quantum_calculations_run=False, current_production_runtime_qualified=False,
                  missing_conformers_supplied=False, global_ensemble_convergence_proven=False)
    atomic_write_json(a.out, result)


if __name__ == '__main__':
    main()
