#!/usr/bin/env python3
"""Replay cached Fc and small solvent-surface controls from packaged files.

Electronic/component arithmetic needs only the standard library. Full stored
force/Hessian/mode/nuclear replay additionally requires the pinned NumPy wheel.
No SCF, CUDA, xTB, production redox import or original calculation path is used.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local(root, name):
    relative = Path(name)
    if relative.is_absolute() or not relative.parts or '..' in relative.parts:
        raise ValueError('Package members must have contained relative paths')
    target = root / relative
    if any((root / Path(*relative.parts[:i])).is_symlink() for i in range(1, len(relative.parts)+1)):
        raise ValueError('Package members cannot be symlinks')
    target = target.resolve(strict=True)
    if not target.is_relative_to(root.resolve()) or not target.is_file():
        raise ValueError('Package member escapes its declared root')
    return target


def read_json(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON value: '+value)
    return json.loads(Path(path).read_text(), parse_constant=invalid)


def module(root, path, name):
    spec = importlib.util.spec_from_file_location(name, local(root, path))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def verify_inventory(root):
    manifest = read_json(local(root, 'data/manifest.json'))
    if (manifest['schema'] != 1 or manifest['primary_candidate_admission'] is not False
            or manifest['polymer_free_energy_qualification'] is not False):
        raise ValueError('Cached controls cannot supply primary candidate or polymer G admission')
    for path, expected in manifest['files_sha256'].items():
        if digest(local(root, path)) != expected:
            raise ValueError('Packaged file differs from its checksum: '+path)
    return manifest


def check_quantum(root, manifest):
    import numpy as np
    equations = module(root, 'nuclear_equations.py', 'cached_control_nuclear_equations')
    modes = module(root, 'projected_modes.py', 'cached_control_projected_modes')
    expected = {(name, arm) for name in ('formaldehyde', 'methyl', 'ammonia', 'chloromethyl')
                for arm in ('legacy', 'smooth')}
    rows = manifest['quantum_states']
    if len(rows) != 8 or {(r['species'], r['arm']) for r in rows} != expected:
        raise ValueError('The declared four-species/two-model raw control panel is incomplete')
    checked = []
    for row in rows:
        files = row['files']; get = lambda key: local(root, files[key])
        energy = read_json(get('energy')); nuclear = read_json(get('nuclear'))
        hessian = np.load(get('hessian'), allow_pickle=False)
        gradient = np.load(get('gradient'), allow_pickle=False)
        lines = get('geometry').read_text().splitlines()
        if int(lines[0]) != 4 or len(lines[2:]) != 4:
            raise ValueError('The cached control must have exactly four atoms')
        xyz = np.array([[float(v) for v in line.split()[1:]] for line in lines[2:]])
        mult = 1 if row['species'] in ('formaldehyde', 'ammonia') else 2
        if (energy['charge'] != 0 or energy['mult'] != mult or energy['solvent'] != 'acetonitrile'
                or energy['opt_xc'] != 'r2scan' or energy['opt_basis'] != 'def2-svp'
                or energy['sp_xc'] != 'wb97m-v' or energy['sp_basis'] != 'def2-tzvpd'
                or not energy['geometry_optimization_converged'] or not energy['optimization_SCF_converged']
                or not energy['converged_smd'] or not energy['converged_gas']
                or energy['primary_candidate_admission']):
            raise ValueError('Cached quantum state identity/protocol/convergence differs')
        profile = energy['CDS_profile']
        if profile != ('legacy_SMD_CDS' if row['arm'] == 'legacy' else manifest['experimental_model']):
            raise ValueError('Legacy and alternate surface models were mixed')
        if energy['energy_geom_sha1'] != hashlib.sha1(get('geometry').read_bytes()).hexdigest()[:16]:
            raise ValueError('Cached electronic energy belongs to another geometry')
        for medium in ('smd', 'gas'):
            if (not math.isfinite(energy[f'e_{medium}_eV'])
                    or abs(energy[f'e_{medium}_Ha']*27.211386245988-energy[f'e_{medium}_eV']) > 1e-8
                    or abs(energy[f'{medium}_spin_contam']) > .1):
                raise ValueError('Cached electronic units or spin gate failed')
        if (hessian.shape != (12, 12) or gradient.shape != (4, 3)
                or not np.isfinite(hessian).all() or not np.isfinite(gradient).all()
                or np.max(abs(hessian-hessian.T)) > 1e-6 or np.max(abs(gradient)) > .00045):
            raise ValueError('Cached Hessian/stationarity checks failed')
        evaluations = row['evaluations']
        if len(evaluations) != 49:
            raise ValueError('Complete baseline and force-displacement records are required')
        records = [read_json(local(root, p)) for p in evaluations]
        if any(not r['SCF_converged'] or abs(r['spin']['spin_contam']) > .1
               or r['arm'] != row['arm'] for r in records):
            raise ValueError('A stored displaced SCF/spin/model check failed')
        coordinates = np.array(records[0]['coordinates_Bohr'])
        if (not np.allclose(coordinates, xyz/.529177210903, rtol=0, atol=1e-10)
                or not np.allclose(records[0]['gradient_Ha_Bohr'], gradient, rtol=0, atol=1e-13)):
            raise ValueError('Baseline force/coordinates do not match their geometry')
        force_errors, hessian_errors = [], []
        index = 1
        for atom in range(4):
            for axis in range(3):
                for step in (.0005, .001):
                    a, b = records[index:index+2]; index += 2
                    delta = np.zeros_like(coordinates); delta[atom, axis] = step
                    if (not np.allclose(np.array(a['coordinates_Bohr'])-coordinates, delta, rtol=0, atol=1e-12)
                            or not np.allclose(np.array(b['coordinates_Bohr'])-coordinates, -delta, rtol=0, atol=1e-12)):
                        raise ValueError('A stored force displacement belongs to another coordinate or step')
                    force_errors.append(abs((a['E_Ha']-b['E_Ha'])/(2*step)-gradient[atom, axis]))
                    if (atom, axis) in {(0, 1), (1, 2)}:
                        finite = (np.array(a['gradient_Ha_Bohr'])-np.array(b['gradient_Ha_Bohr'])).ravel()/(2*step)
                        analytic = hessian[:, 3*atom+axis]
                        hessian_errors.append(float(np.linalg.norm(finite-analytic)/max(np.linalg.norm(analytic), 1e-14)))
        if max(force_errors) > 1e-6 or max(hessian_errors) > .01:
            raise ValueError('Stored full-energy force or sampled-Hessian arithmetic failed')
        mass = np.array(nuclear['atomic_masses_amu'])
        if mass.shape != (4,) or not np.isfinite(mass).all() or np.min(mass) <= 0:
            raise ValueError('Invalid cached isotope-average masses')
        frequencies, rank = modes.projected_modes(hessian, xyz, mass)
        frequency_error = float(np.max(abs(frequencies-np.array(nuclear['frequencies_cm']))))
        if rank != 6 or len(frequencies) != 6 or frequency_error > 1e-5 or min(frequencies) < -5:
            raise ValueError('Cached mode projection or internal-frequency gate failed')
        if nuclear['imposed_comparison_symmetry_number'] != 1 or nuclear['primary_native_thermal_admission']:
            raise ValueError('Comparison symmetry/scope changed')
        errors = {}
        for cutoff in (25, 50, 100):
            value = equations.nuclear_equations(frequencies, xyz, mass, 1, cutoff)
            errors[str(cutoff)] = abs(value-nuclear['nuclear_by_cutoff_eV'][str(cutoff)])
        if not all(math.isfinite(v) and v <= 1e-8 for v in errors.values()):
            raise ValueError('Cached nuclear partition does not reconstruct')
        reference = row['reference_metrics']
        if (abs(max(force_errors)-reference['maximum_force_error_Ha_Bohr']) > 1e-12
                or abs(max(hessian_errors)-reference['maximum_Hessian_column_relative_error']) > 1e-10):
            raise ValueError('Stored derivative maxima differ from the independent source review')
        checked.append(dict(species=row['species'], arm=row['arm'], frequency_max_error_cm=frequency_error,
            maximum_force_error_Ha_Bohr=float(max(force_errors)),
            maximum_Hessian_column_relative_error=max(hessian_errors), nuclear_errors_eV=errors))
    return checked


def replay(root, electronic_only=False):
    root = Path(root).resolve()
    manifest = verify_inventory(root)
    claims = module(root, 'claims.py', 'cached_control_claims')
    def load(name):
        return read_json(local(root, manifest['evidence'][name]))
    fc = claims.audit_fc_functional_controls(load)
    surface = claims.audit_surface_controls(load)
    if fc is None or surface is None:
        raise ValueError('Both frozen electronic and surface evidence sets are required')
    geometries = {}
    expected = {(method, state) for method in ('wb97mv', 'b3lyp', 'pbe0') for state in ('neu', 'ox')}
    if (len(manifest['electronic_states']) != 6
            or {(r['method'], r['state']) for r in manifest['electronic_states']} != expected):
        raise ValueError('Exactly six Fc electronic states are required')
    for row in manifest['electronic_states']:
        if row['energy_evidence'] != f"Fc_functional_{row['method']}_{row['state']}.json":
            raise ValueError('Fc method/state identity does not match its stored evidence')
        geometry = local(root, row['geometry']); energy = load(row['energy_evidence'])
        value = hashlib.sha1(geometry.read_bytes()).hexdigest()[:16]
        if value != energy['energy_geom_sha1']:
            raise ValueError('Stored Fc electronic state belongs to another geometry')
        state = row['state']
        if state in geometries and geometry.read_bytes() != geometries[state]:
            raise ValueError('The functional controls do not share the declared fixed geometry')
        geometries[state] = geometry.read_bytes()
    quantum = [] if electronic_only else check_quantum(root, manifest)
    return dict(status='packaged_electronic_and_component_arithmetic_reconstructed' if electronic_only
                else 'packaged_control_forces_Hessians_modes_nuclear_and_Fc_reconstructed',
        files_verified=len(manifest['files_sha256']), Fc_functional_controls=fc,
        surface_component_controls=surface, quantum_states_reconstructed=len(quantum), states=quantum,
        quantum_calculations_run=False, CUDA_or_xTB_loaded=False,
        original_physical_calculations_repeated=False, polymer_G_qualified=False,
        primary_candidate_admission=False, experimental_accuracy_inferred=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--electronic-only', action='store_true')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = replay(args.package, args.electronic_only)
    if args.out:
        with args.out.open('x') as stream:
            json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'states'}))


if __name__ == '__main__': main()
