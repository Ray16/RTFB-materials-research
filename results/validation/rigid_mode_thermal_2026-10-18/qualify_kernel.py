"""Compare explicit-mode kernel with freshly projected CLI controls, on CPU."""
import hashlib
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time

for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
            'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
    os.environ[key] = '2'
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
with allocation_guard():
    resources = host_capacity(2)
    cpus = cpu_allocation(2)
    os.sched_setaffinity(0, cpus)

import numpy as np
from ase.io import read
from redox.core.io import atomic_write_json
from redox.qm import qrrho as Q
from redox.qm.vibrations import project_rigid_modes, write_xtb_hessian, parse_xtb_thermal_eV
from redox.qm.native_thermal import XTB
from pyscf.dispersion import dftd4 as caller_D4_interface

BASE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    plan_path = args.manifest
    plan = json.loads(plan_path.read_text())
    for path, digest in plan['producer_sha256'].items():
        assert sha(path) == digest, path
    out = args.out
    out.mkdir(exist_ok=False)
    result = dict(status='running', protocol=Q.PROTOCOL, resources=resources, cpu_affinity=cpus,
                  host=os.uname().nodename, pid=os.getpid(), started=time.time(), states=[],
                  process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                  plan_sha256=sha(plan_path), production_admission=False)
    atomic_write_json(out / 'summary.json', result)
    for task in plan['tasks']:
        work = out / task['id'] / task['state'] / Path(task['native_record']).stem
        work.mkdir(parents=True, exist_ok=False)
        item = dict(id=task['id'], state=task['state'], native_record=task['native_record'])
        result['states'].append(item)
        try:
            for path, digest in task['source_sha256'].items():
                assert sha(path) == digest, path
            atoms = read(task['geometry'])
            masses = Q.atomic_masses(atoms.numbers)
            matrix = np.load(task['hessian'])
            projected, rank = project_rigid_modes(matrix, atoms.positions, masses)
            np.save(work / 'projected_hessian.npy', projected)
            write_xtb_hessian(work / 'hessian', projected)
            (work / 'mol.xyz').write_bytes(Path(task['geometry']).read_bytes())
            frequencies, internal_rank = Q.internal_vibrations(matrix, atoms.positions, masses)
            assert rank == internal_rank
            comparisons, symmetry_number = [], None
            rot, average, linear = Q.inertia_properties(atoms.positions, masses)
            numbers = np.fromstring(' '.join(line for line in (work/'hessian').read_text().splitlines()
                                             if not line.startswith('$')), sep=' ')
            serialized = numbers.reshape(projected.shape).T
            roots = np.repeat(np.sqrt(masses * Q.AMUTOAU), 3)
            eigenvalues = np.linalg.eigvalsh(serialized / np.outer(roots, roots))
            serialized_frequencies = np.sign(eigenvalues) * np.sqrt(abs(eigenvalues)) * Q.AUTORCM
            for cutoff in (25, 50, 100):
                log = work / f'thermo_{cutoff}.log'
                with log.open('w') as stream:
                    process = subprocess.run([str(XTB), 'thermo', '--sthr', str(cutoff),
                                              '--temp', '298.15', 'mol.xyz', 'hessian'],
                                             cwd=work, stdout=stream, stderr=subprocess.STDOUT, timeout=120)
                assert process.returncode == 0
                text = log.read_text()
                number = int(re.search(r'rotational number\s+(\d+)', text)[1])
                if symmetry_number is None:symmetry_number = number
                assert number == symmetry_number
                cli_count = int(re.search(r'# frequencies\s+(\d+)', text)[1])
                retained = abs(serialized_frequencies[abs(serialized_frequencies) > 1.])
                assert cli_count == len(retained)
                partial = Q.kernel_thermal(retained, rot, average, sum(masses), number,
                                           cutoff_cm=cutoff, linear=linear)['G_Ha'] * Q.HARTREE_EV
                cli = parse_xtb_thermal_eV(text)
                difference = abs(partial - cli)
                comparisons.append(dict(cutoff_cm=cutoff, CLI_retained_modes=cli_count,
                                         reference='exact serialized Cartesian Hessian and CLI mode filter',
                                         explicit_filtered_kernel_eV=partial, CLI_eV=cli,
                                         absolute_difference_eV=difference))
                assert difference <= plan['CLI_reproduction_bound_eV'], (difference, log)
            partition = Q.nuclear_thermals(matrix, atoms.positions, atoms.numbers, symmetry_number)
            assert partition['internal_mode_count'] == 3 * len(atoms) - rank
            item.update(status='explicit_internal_kernel_qualified', partition=partition,
                        symmetry_number=symmetry_number, CLI_comparisons=comparisons,
                        source_sha256=task['source_sha256'],
                        artifacts_sha256={str(p): sha(p) for p in work.iterdir() if p.is_file()})
            for path, digest in task['source_sha256'].items():
                assert sha(path) == digest, path
        except Exception as exc:
            item.update(status='failed_preserved', error=repr(exc))
        atomic_write_json(work / 'result.json', item)
        atomic_write_json(out / 'summary.json', result)
        print(item['id'], item['state'], item['status'],
              item.get('partition', {}).get('retained_soft_mode_indices'), item.get('error', ''), flush=True)
    for path, digest in plan['producer_sha256'].items():
        assert sha(path) == digest, path
    result.update(status='bounded_qualification_complete', finished=time.time(),
                  thermodynamic_kernel_worker=Q.kernel_execution(),
                  caller_D4_library_preloaded=True,
                  qualified=sum(r['status'] == 'explicit_internal_kernel_qualified' for r in result['states']),
                  failed=sum(r['status'] == 'failed_preserved' for r in result['states']))
    atomic_write_json(out / 'summary.json', result)


if __name__ == '__main__':
    main()
