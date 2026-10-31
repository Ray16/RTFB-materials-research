#!/usr/bin/env python3
"""Check the isolated barrier correction on closed/open-shell MeCN derivatives.

Both derivatives and their displaced SCF calculations are retained. These fixed
geometries are method controls, never minimum/free-energy or material evidence.
"""
import argparse
from contextlib import nullcontext
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    case = args.plan.parent
    sources = case.parent / 'sources'
    plan = json.loads(args.plan.read_text())
    if plan['role'] != 'screening_molecular_control':
        raise ValueError('Declared molecular implementation control required')
    def check_sources():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest() != h
               for p, h in plan['source_sha256'].items()):
            raise ValueError('Frozen control source changed')
    check_sources()
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in
            ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
             'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')):
        raise ValueError('Two CPUs and matching native caps required')
    record = dict(status='screening_molecular_control_running', started=time.time(),
                  owner=dict(host=os.uname().nodename, pid=os.getpid(),
                             process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
                  native_G_qualified=False, production_default_changed=False)
    atomic_write_json(case / 'status.json', record)
    try:
        import cupy as cp
        import numpy as np
        from pyscf import gto
        from redox.qm import dft
        from redox.qm.cds_finite import finer_cds_hessian
        helper = load(sources / 'helper.py', 'frozen_molecular_screening_helper')
        model = json.loads((sources / 'molecule.json').read_text())
        state = model['states'][plan['molecular_control']]
        atoms = state.get('atoms', model['atoms'])
        symbols = [a[0] for a in atoms]
        coordinates = np.asarray([a[1] for a in atoms], dtype=np.float64)
        results = {}
        for variant in ('installed', 'corrected'):
            folder = case / variant
            folder.mkdir()
            counter = 0
            def solve(coords):
                nonlocal counter
                gc.collect()
                cp.get_default_memory_pool().free_all_blocks()
                point = folder / f'point_{counter:03d}'
                counter += 1
                point.mkdir()
                mol = gto.M(atom=list(zip(symbols, coords.tolist())), unit='Bohr',
                            basis=model['basis'], charge=state['charge'],
                            spin=state['mult'] - 1, verbose=0, max_memory=2000)
                mf = dft._build_smd_mf(mol, model['xc'], model['solvent'], 'gpu', disp=model['disp'])
                mf.conv_tol = model['SCF_conv_tol']
                mf.conv_tol_grad = model['SCF_conv_tol_grad']
                mf.max_cycle = model['SCF_max_cycle']
                mf.grids.level = model['grid_level']
                mf.chkfile = str(point / 'SCF.h5')
                mf.callback = lambda e: atomic_write_json(point / 'progress.json',
                    dict(cycle=int(e['cycle']) + 1, E_Ha=float(e['e_tot']), updated=time.time()))
                mf.kernel()
                if not mf.converged:
                    raise ValueError('Control SCF did not converge')
                grad = mf.nuc_grad_method()
                grad.grid_response = True
                value = grad.kernel()
                gradient = np.asarray(value.get() if hasattr(value, 'get') else value)
                if not np.isfinite(gradient).all() or not np.isfinite(mf.e_tot):
                    raise ValueError('Nonfinite converged control evaluation')
                spin = dft._orbital_diag(mf, state['mult'] - 1)
                if abs(spin['spin_contam']) > .1:
                    raise ValueError('Unchanged spin contamination gate failed')
                row = dict(coordinates_Bohr=coords.tolist(), energy_Ha=float(mf.e_tot),
                           gradient_Ha_Bohr=gradient.tolist(), SCF_converged=True,
                           spin=spin, variant=variant, native_G_qualified=False)
                atomic_write_json(point / 'evaluation.json', row)
                return mf, gradient, row
            context = nullcontext() if variant == 'installed' else helper.corrected_screening(
                sources / 'kernel.cu', folder / 'implementation')
            with context, finer_cds_hessian(1e-5):
                mf, gradient, center = solve(coordinates)
                mf.cphf_grids = mf.grids.copy()
                hobj = mf.Hessian()
                hobj.grid_response = True
                hobj.max_cycle = 100
                value = hobj.kernel()
                raw = np.asarray(value.get() if hasattr(value, 'get') else value)
                matrix = raw.transpose(0, 2, 1, 3).reshape(coordinates.size, coordinates.size)
                np.save(folder / 'Hessian.npy', matrix, allow_pickle=False)
                if not np.isfinite(matrix).all():
                    raise ValueError('Nonfinite analytical Hessian')
                asymmetry = float(abs(matrix - matrix.T).max())
                checks = []
                del hobj, mf
                # All Cartesian columns, not a sampled subset, at both steps.
                for atom in range(len(symbols)):
                    for axis in range(3):
                        for step in model['FD_steps_Bohr']:
                            displaced = []
                            for sign in (-1, 1):
                                coords = coordinates.copy()
                                coords[atom, axis] += sign * step
                                mf, _, row = solve(coords)
                                displaced.append(row)
                                del mf
                            finite = (np.asarray(displaced[1]['gradient_Ha_Bohr']) -
                                      np.asarray(displaced[0]['gradient_Ha_Bohr'])).ravel() / (2 * step)
                            analytic = matrix[:, 3 * atom + axis]
                            relative = float(np.linalg.norm(finite - analytic) / max(np.linalg.norm(finite), 1e-14))
                            energy_gradient = (displaced[1]['energy_Ha'] - displaced[0]['energy_Ha']) / (2 * step)
                            checks.append(dict(atom=atom, axis=axis, step_Bohr=step,
                                relative_Hessian_column_error=relative,
                                energy_gradient_error_Ha_Bohr=energy_gradient - gradient[atom, axis],
                                displaced=displaced))
                            atomic_write_json(folder / 'derivatives.json', dict(checks=checks, complete=False))
                results[variant] = dict(center=center, asymmetry_Ha_Bohr2=asymmetry, checks=checks)
                atomic_write_json(folder / 'derivatives.json', dict(results[variant], complete=True))
                atomic_write_json(case / 'comparison.json', dict(results=results, native_G_qualified=False))
        a, b = [results[v] for v in ('installed', 'corrected')]
        differences = dict(energy_Ha=abs(b['center']['energy_Ha'] - a['center']['energy_Ha']),
            gradient_Ha_Bohr=float(abs(np.asarray(b['center']['gradient_Ha_Bohr']) - np.asarray(a['center']['gradient_Ha_Bohr'])).max()),
            Hessian_Ha_Bohr2=float(abs(np.load(case / 'corrected/Hessian.npy') - np.load(case / 'installed/Hessian.npy')).max()))
        passed = (all(differences[k] <= v for k, v in model['installed_corrected_limits'].items())
                  and all(r['asymmetry_Ha_Bohr2'] <= model['Hessian_asymmetry_limit_Ha_Bohr2'] for r in results.values())
                  and all(c['relative_Hessian_column_error'] <= model['FD_relative_limit']
                          and abs(c['energy_gradient_error_Ha_Bohr']) <= model['energy_gradient_absolute_limit_Ha_Bohr']
                          for r in results.values() for c in r['checks']))
        atomic_write_json(case / 'comparison.json', dict(results=results, differences=differences,
                          declared_gates_passed=passed, native_G_qualified=False))
        check_sources()
        record.update(status='screening_molecular_control_finished_review_pending',
                      declared_gates_passed=passed, differences=differences)
    except BaseException as error:
        record.update(status='screening_molecular_control_failed_preserved', error=repr(error))
        raise
    finally:
        record['finished'] = time.time()
        atomic_write_json(case / 'status.json', record)


if __name__ == '__main__':
    main()
