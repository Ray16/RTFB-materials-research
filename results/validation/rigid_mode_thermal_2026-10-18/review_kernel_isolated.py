"""Independent projection, analytic qRRHO equations and exact-input CLI review."""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
with allocation_guard():
    resources = host_capacity(1)
    cpus = cpu_allocation(1)
    os.sched_setaffinity(0, cpus)
import numpy as np
from ase.io import read
from redox.core.io import atomic_write_json

BASE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
HARTREE_EV = 27.211386245988


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
    if rotational[-1] < 1e-10:
        entropy_rot = gas_constant * (math.log(1 / (beta * rotational[0] / 219474.63067 * symmetry)) + 1)
        enthalpy_rot = gas_constant * temperature
    else:
        entropy_rot = gas_constant * (.5 * (math.log(np.pi / np.prod(rotational / 219474.63067) / beta**3)
                                                    - 2 * math.log(symmetry)) + 1.5)
        enthalpy_rot = 1.5 * gas_constant * temperature
    entropy_trans = 2.2868 * (5 * math.log10(temperature) + 3 * math.log10(sum(mass))) - 2.3135
    enthalpy_trans = 2.5 * gas_constant * temperature
    return (.5 * sum(omega) + ((enthalpy_vib + enthalpy_rot + enthalpy_trans)
            - temperature * (entropy_vib + entropy_rot + entropy_trans)) / 1000 / 627.50947428) * HARTREE_EV


def main():
    plan_path = BASE / 'kernel_qualification_plan_isolated.json'
    source = BASE / 'kernel_qualification_isolated/summary.json'
    plan, data = json.loads(plan_path.read_text()), json.loads(source.read_text())
    assert data['status'] == 'bounded_qualification_complete'
    assert len(data['states']) == len(plan['tasks'])
    checks, states, excluded = [], [], []
    evidence = {str(source): sha(source), str(plan_path): sha(plan_path)}
    mass_source = BASE / 'xtb_source/src_param_atomicmass.f90'
    array = mass_source.read_text().split('= amutoau * [', 1)[1].split(']', 1)[0]
    expected_masses = np.array([float(v) for v in re.findall(r'([0-9.]+)_wp', array)])
    assert len(expected_masses) == 118
    evidence[str(mass_source)] = sha(mass_source)

    def check(key, name, passed):checks.append(dict(key=key, check=name, passed=bool(passed)))

    for task, record in zip(plan['tasks'], data['states'], strict=True):
        key = task['id'] + ':' + task['state'] + ':' + Path(task['native_record']).stem
        check(key, 'source_identity', record['native_record'] == task['native_record'])
        if record['status'] != 'explicit_internal_kernel_qualified':
            excluded.append(dict(key=key, error=record.get('error')))
            continue
        start = len(checks)
        for path, digest in {**record['source_sha256'], **record['artifacts_sha256']}.items():
            check(key, 'hash:' + path, sha(path) == digest)
            evidence[path] = digest
        atoms = read(task['geometry'])
        mass = expected_masses[atoms.numbers - 1]
        partition = record['partition']
        check(key, 'upstream_atomic_masses', abs(mass - partition['atomic_masses_amu']).max() < 1e-8)
        xyz = atoms.positions - np.average(atoms.positions, axis=0, weights=mass)
        candidates = [(np.sqrt(mass)[:, None] * axis).ravel() for axis in np.eye(3)]
        candidates += [(np.sqrt(mass)[:, None] * np.cross(xyz, axis)).ravel() for axis in np.eye(3)]
        rigid = []
        for candidate in candidates:
            v = candidate.copy()
            for _ in range(2):
                for u in rigid:v -= np.dot(u, v) * u
            if np.linalg.norm(v) > 1e-10 * max(1., np.linalg.norm(candidate)):
                rigid.append(v / np.linalg.norm(v))
        rigid = np.array(rigid).T
        projector = np.eye(3 * len(atoms)) - rigid @ rigid.T
        eigen, vectors = np.linalg.eigh(projector)
        internal = vectors[:, eigen > .5]
        matrix = np.load(task['hessian'])
        repeated = np.repeat(np.sqrt(mass * (1.660539040e-27 / 9.10938356e-31)), 3)
        weighted = (matrix + matrix.T) / (2 * np.outer(repeated, repeated))
        force = np.linalg.eigvalsh(internal.T @ weighted @ internal)
        frequencies = np.sign(force) * np.sqrt(abs(force)) * 219474.63067
        check(key, 'independent_projected_spectrum', len(frequencies) == partition['internal_mode_count']
              and abs(frequencies - partition['frequencies_cm']).max() < 1e-5)
        check(key, 'all_internal_modes_included', len(frequencies) == 3 * len(atoms) - len(rigid.T)
              and len(rigid.T) == partition['rigid_rank'])
        check(key, 'minimum_and_nonzero', min(frequencies) >= -5 and np.all(frequencies != 0))
        check(key, 'soft_modes_retained', np.where(abs(frequencies) <= 1)[0].tolist()
              == partition['retained_soft_mode_indices'])
        deviations = []
        for cutoff in (25, 50, 100):
            value = nuclear_equations(frequencies, atoms.positions, mass, record['symmetry_number'], cutoff)
            error = abs(value - partition['nuclear_by_cutoff_eV'][str(cutoff)])
            check(key, 'independent_qRRHO:' + str(cutoff), error < 1e-8)
            deviations.append(error)
            control = next(c for c in record['CLI_comparisons'] if c['cutoff_cm'] == cutoff)
            log_path = next(Path(p) for p in record['artifacts_sha256'] if Path(p).name == f'thermo_{cutoff}.log')
            text = log_path.read_text()
            match = re.search(r'G\(RRHO\) contrib\.\s+([-+0-9.Ee]+)', text)
            count = re.search(r'# frequencies\s+(\d+)', text)
            symmetry = re.search(r'rotational number\s+(\d+)', text)
            cli = float(match[1]) * HARTREE_EV
            check(key, 'raw_CLI_control:' + str(cutoff), int(count[1]) == control['CLI_retained_modes']
                  and int(symmetry[1]) == record['symmetry_number']
                  and abs(cli - control['CLI_eV']) < 1e-8
                  and abs(abs(cli - control['explicit_filtered_kernel_eV']) - control['absolute_difference_eV']) < 1e-12
                  and control['absolute_difference_eV'] <= 1e-6)
        states.append(dict(key=key, passed=all(c['passed'] for c in checks[start:]),
                           maximum_independent_qRRHO_difference_eV=max(deviations),
                           soft_modes=len(partition['retained_soft_mode_indices'])))
    for path, digest in evidence.items():check('final_inventory', path, sha(path) == digest)
    out = BASE / 'kernel_isolated_independent_review'
    out.mkdir(exist_ok=False)
    summary = dict(status='explicit_internal_qRRHO_reviewed', states=len(states), excluded=excluded,
                   checks=len(checks), failed_checks=sum(not c['passed'] for c in checks),
                   all_states_passed=all(s['passed'] for s in states), protocol=data['protocol'],
                   maximum_qRRHO_equation_difference_eV=max(s['maximum_independent_qRRHO_difference_eV'] for s in states),
                   soft_mode_records=sum(s['soft_modes'] > 0 for s in states),
                   evidence_sha256=evidence, resources=resources, cpu_affinity=cpus,
                   review_source_sha256=sha(Path(__file__)), production_admission=False, finished=time.time())
    atomic_write_json(out / 'summary.json', summary)
    atomic_write_json(out / 'checks.json', checks)
    atomic_write_json(out / 'states.json', states)
    print({k: v for k, v in summary.items() if k not in ('evidence_sha256', 'resources', 'protocol')})
    assert summary['failed_checks'] == 0 and not excluded


if __name__ == '__main__':main()
