"""Independent projected-partition reader; no production thermochemistry imports.

Reconstruct mode exclusion, frequencies, nuclear partitioning and spin terms from
the backing Hessian. The native-surface admission remains the caller's job.
"""
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[3]
LIBRARY = Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/lib/libxtb.so.6.7.1')
LIBRARY_SHA = 'c405797f34ed96470f23158406e6d8a8e102681d3c3cbbaea16692f514f22853'
PROTOCOL = dict(schema=1, model='projected internal-mode xTB6.7.1 qRRHO',
    projection='mass-weighted SVD internal subspace; symmetric full-precision Cartesian Hessian',
    masses='installed xTB NIST isotopic-average atomic masses',
    rigid_mode_exclusion='geometric subspace; no frequency threshold',
    nuclear_kernel='installed GNU Fortran thermodyn; explicit internal modes',
    kernel_isolation='dedicated CPU subprocess; no PySCF/D4 Python imports',
    library_sha256=LIBRARY_SHA, frequency_scale=1., temperature_K=298.15,
    rotor_cutoffs_cm=[25, 50, 100], primary_rotor_cutoff_cm=50,
    imaginary_mode_policy='absolute value only within the existing -5 cm^-1 minimum gate',
    zero_internal_frequency_policy='refuse; no artificial floor',
    symmetry='xTB symmetry detection for the exact input geometry',
    pressure='same translation convention as installed xTB thermodyn')
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def nuclear_equations(frequencies, coordinates, mass, symmetry, cutoff):
    """Independent Python equations from the pinned upstream thermodyn source."""
    temperature, gas_constant, boltzmann = 298.15, 1.98726, 3.166808578545117e-6
    beta = 1 / (boltzmann * temperature)
    omega = np.abs(frequencies) / 219474.63067
    centered = coordinates - np.average(coordinates, axis=0, weights=mass)
    inertia = sum(m * (np.eye(3) * np.dot(x, x) - np.outer(x, x)) for m, x in zip(mass, centered))
    inertia += np.array([[1., 2., 4.], [2., 3., 5.], [4., 5., 6.]]) * 1e-10
    moments = np.linalg.eigvalsh(inertia)
    rotational = np.array([16.8576522 / v if v >= 3e-4 else 0. for v in moments])[::-1]
    average_SI = 1e-47 * 1.66053 * np.mean(moments)
    average_au = average_SI / 9.10938356e-31 * (1 / .52917726)**2 * 1e20
    mu = .5 / (omega + 1e-14)
    reduced = mu * average_au / (mu + average_au)
    harmonic_weight = 1 / (1 + (cutoff / np.abs(frequencies))**4)
    x = omega * beta
    sv_ho = x / np.expm1(x) - np.log(-np.expm1(-x))
    sv_fr = .5 + .5 * np.log(2 * np.pi * reduced / beta)
    entropy_vib = gas_constant * np.sum(harmonic_weight * sv_ho + (1 - harmonic_weight) * sv_fr)
    enthalpy_vib = np.sum(omega / np.expm1(x)) * gas_constant / boltzmann
    if len(mass) == 1:
        entropy_rot, enthalpy_rot = 0., 0.
    elif rotational[-1] < 1e-10:
        entropy_rot = gas_constant * (math.log(1 / (beta * rotational[0] / 219474.63067 * symmetry)) + 1)
        enthalpy_rot = gas_constant * temperature
    else:
        entropy_rot = gas_constant * (.5 * (math.log(np.pi / np.prod(rotational / 219474.63067) / beta**3)
                                                    - 2 * math.log(symmetry)) + 1.5)
        enthalpy_rot = 1.5 * gas_constant * temperature
    entropy_trans = 2.2868 * (5 * math.log10(temperature) + 3 * math.log10(sum(mass))) - 2.3135
    enthalpy_trans = 2.5 * gas_constant * temperature
    return (.5 * sum(omega) + ((enthalpy_vib + enthalpy_rot + enthalpy_trans)
            - temperature * (entropy_vib + entropy_rot + entropy_trans)) / 1000 / 627.50947428) * 27.211386245988


def read_partition(native_path, geometry):
    """Return independently reconstructed values/evidence, or fail closed."""
    try:
        native_path, geometry = Path(native_path), Path(geometry)
        proof_path = ROOT/'config/projected_thermal_validation.json'
        proof = json.loads(proof_path.read_text())
        seal = proof.pop('certificate_content_sha256')
        if (hashlib.sha256(json.dumps(proof, sort_keys=True, separators=(',', ':')).encode()).hexdigest() != seal
                or proof['status'] != 'independently_qualified' or proof['protocol'] != PROTOCOL
                or proof['library_sha256'] != LIBRARY_SHA or sha(LIBRARY) != LIBRARY_SHA
                or proof['kernel_worker_source_sha256'] != sha(ROOT/'src/redox/qm/xtb_kernel.py')
                or proof['qrrho_source_sha256'] != sha(ROOT/'src/redox/qm/qrrho.py')):
            return None
        evidence = {str(proof_path): sha(proof_path)}
        for path, digest in proof['evidence_sha256'].items():
            if sha(ROOT/path) != digest:return None
            evidence[str(ROOT/path)] = digest
        native = json.loads(native_path.read_text())
        matrix_path = native_path.with_suffix('')/'hessian_cartesian.npy'
        identity = dict(geometry_sha256=sha(geometry), raw_native_sha256=sha(native_path),
                        Cartesian_Hessian_sha256=sha(matrix_path), charge=native['identity']['charge'],
                        mult=native['identity']['mult'], protocol=PROTOCOL)
        if identity['geometry_sha256'] != native['identity']['geometry_sha256']:return None
        address = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
        target = native_path.parent.parent/'projected_thermal'/(address+'.json')
        record = json.loads(target.read_text())
        if (record['status'] != 'projected_internal_qrrho_checked' or record['identity'] != identity
                or record['implementation_certificate_sha256'] != seal):return None
        execution,owner=record.get('kernel_execution') or {},record.get('owner') or {}
        if (execution.get('library_sha256')!=LIBRARY_SHA
                or execution.get('source_sha256')!=proof['kernel_worker_source_sha256']
                or execution.get('pid')==owner.get('pid') or execution.get('host')!=owner.get('host')
                or set(execution.get('thread_limits',{}).values())!={'1'}
                or not execution.get('cpu_affinity')
                or not set(execution['cpu_affinity'])<=set(owner.get('cpu_affinity',[]))):return None
        resolved={}
        for path, digest in record['artifacts_sha256'].items():
            local=target.with_suffix('')/Path(path).name
            candidate=local if local.is_file() else Path(path)
            if sha(candidate) != digest:return None
            evidence[str(candidate)] = digest
            resolved[str(candidate)] = digest
        evidence.update({str(p): sha(p) for p in (native_path, matrix_path, geometry, target)})
        files = {Path(p).name: Path(p) for p in resolved}
        required = {'projected_hessian.npy', 'hessian', 'mol.xyz', 'qrrho_producer.py', 'record_producer.py',
                    'thermo_CLI_25.log', 'thermo_CLI_50.log', 'thermo_CLI_100.log'}
        if not required <= set(files) or files['mol.xyz'].read_bytes() != geometry.read_bytes():return None
        mass_path = next(ROOT/p for p in proof['evidence_sha256'] if p.endswith('src_param_atomicmass.f90'))
        array = mass_path.read_text().split('= amutoau * [', 1)[1].split(']', 1)[0]
        table = np.array([float(v) for v in re.findall(r'([0-9.]+)_wp', array)])
        if len(table) != 118:return None
        atoms = read(geometry)
        mass = table[atoms.numbers - 1]
        if abs(mass - record['atomic_masses_amu']).max() > 1e-8:return None
        centered = atoms.positions - np.average(atoms.positions, axis=0, weights=mass)
        candidates = [(np.sqrt(mass)[:, None] * axis).ravel() for axis in np.eye(3)]
        candidates += [(np.sqrt(mass)[:, None] * np.cross(centered, axis)).ravel() for axis in np.eye(3)]
        rigid = []
        for vector in candidates:
            v = vector.copy()
            for _ in range(2):
                for u in rigid:v -= u * np.dot(u, v)
            if np.linalg.norm(v) > 1e-10 * max(1., np.linalg.norm(vector)):rigid.append(v / np.linalg.norm(v))
        basis = np.array(rigid).T
        eigen, vectors = np.linalg.eigh(np.eye(3 * len(atoms)) - basis @ basis.T)
        internal = vectors[:, eigen > .5]
        roots = np.repeat(np.sqrt(mass * (1.660539040e-27 / 9.10938356e-31)), 3)
        hessian = np.load(matrix_path)
        if not np.isfinite(hessian).all() or abs(hessian - hessian.T).max() > 1e-6:return None
        weighted = (hessian + hessian.T) / (2 * np.outer(roots, roots))
        force = np.linalg.eigvalsh(internal.T @ weighted @ internal)
        frequencies = np.sign(force) * np.sqrt(abs(force)) * 219474.63067
        if (len(frequencies) != record['internal_mode_count'] or len(rigid) != record['rigid_rank']
                or len(frequencies) != 3 * len(atoms) - len(rigid)
                or np.max(abs(frequencies - record['frequencies_cm']), initial=0) > 1e-5
                or (len(frequencies) and min(frequencies) < -5) or np.any(frequencies == 0)
                or np.where(abs(frequencies) <= 1)[0].tolist() != record['retained_soft_mode_indices']):return None
        thermal = {}
        for cutoff in (25, 50, 100):
            text = files[f'thermo_CLI_{cutoff}.log'].read_text()
            symmetry = int(re.search(r'rotational number\s+(\d+)', text)[1])
            if symmetry != record['symmetry_number']:return None
            value = nuclear_equations(frequencies, atoms.positions, mass, symmetry, cutoff)
            if not math.isfinite(value) or abs(value - record['nuclear_by_cutoff_eV'][str(cutoff)]) > 1e-8:return None
            control = next(c for c in record['CLI_comparisons'] if c['cutoff_cm'] == cutoff)
            cli = float(re.search(r'G\(RRHO\) contrib\.\s+([-+0-9.Ee]+)', text)[1]) * 27.211386245988
            if (int(re.search(r'# frequencies\s+(\d+)', text)[1]) != control['CLI_retained_modes']
                    or control['reference'] != 'exact serialized Cartesian Hessian and CLI mode filter'
                    or abs(cli - control['CLI_eV']) > 1e-8
                    or abs(abs(cli - control['explicit_filtered_kernel_eV']) - control['absolute_difference_eV']) > 1e-12
                    or not 0 <= control['absolute_difference_eV'] <= 1e-6):return None
            thermal[str(cutoff)] = value
        spin = -8.617333262145e-5 * 298.15 * math.log(identity['mult'])
        total = thermal['50'] + spin
        if (abs(record['g_nuclear_qrrho_eV'] - thermal['50']) > 1e-8
                or abs(record['g_spin_partition_eV'] - spin) > 1e-8
                or abs(record['g_thermal_total_eV'] - total) > 1e-8):return None
        if any(sha(p) != digest for p, digest in evidence.items()):return None
        return dict(g_thermal_total_eV=float(total), nuclear_by_cutoff_eV={k:float(v) for k,v in thermal.items()},
                    path=str(target), sha256=sha(target), evidence_sha256=evidence,
                    protocol=PROTOCOL, address=address)
    except (OSError, ValueError, TypeError, KeyError, AttributeError, StopIteration, np.linalg.LinAlgError):
        return None
