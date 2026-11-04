#!/usr/bin/env python3
"""Qualify a corrected full polymer Hessian against retained derivative evidence.

No active electronic parent or total grafting/occupancy G is assembled here.
The corrected implementation and all response checkpoints have separate identity.
"""
import argparse
import gc
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import time


def load(path, name):
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    case = args.plan.parent
    sources = case.parent / 'sources'
    plan = json.loads(args.plan.read_text())
    if (plan['role'] != 'screening_polymer_Hessian_control'
            or plan['screening_variant'] != 'corrected' or plan['native_G_qualified']):
        raise ValueError('Separate corrected-Hessian implementation control required')
    def check_sources():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest() != h
               for p, h in plan['source_sha256'].items()):
            raise ValueError('Hessian control source changed')
    check_sources()
    if (len(os.sched_getaffinity(0)) != 2 or os.environ.get('NUMPY_MADVISE_HUGEPAGE') != '0'
            or any(os.environ.get(k) != '2' for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT',
                'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'))):
        raise ValueError('Two CPUs, six matching caps and import-time NumPy profile required')
    owner = dict(host=os.uname().nodename, pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    record = dict(status='screening_polymer_Hessian_control_running', owner=owner,
                  started=time.time(), native_G_qualified=False, production_default_changed=False)
    atomic_write_json(case / 'status.json', record)
    try:
        import cupy as cp
        import numpy as np
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft, qrrho
        from redox.qm.cds_finite import finer_cds_hessian
        helper = load(sources / 'helper.py', 'corrected_polymer_screening_helper')
        responses = load(sources / 'response_helper.py', 'corrected_polymer_response_checkpoints')
        checkpoints = load(sources / 'checkpoints.py', 'corrected_polymer_SCF_checkpoints')
        restart = load(sources / 'restart_helper.py', 'corrected_polymer_exact_restart')
        response_restart = (load(sources / 'response_resume_helper.py', 'corrected_polymer_checked_response_reuse')
                            if plan.get('resume_response') else None)
        native = load(ROOT / 'scripts/pipeline/smooth_polymer_native.py', 'unchanged_polymer_Hessian_settings')
        adapter = load(ROOT / 'scripts/pipeline/residual_response_solver.py', 'unchanged_polymer_response_operator')
        model = load(ROOT / 'scripts/pipeline/smooth_cds_control.py', 'unchanged_polymer_MeCN_model')
        monitor = load(ROOT / 'scripts/pipeline/run_pi_stacking.py', 'unchanged_polymer_memory_monitor')
        symbols = plan['symbols']
        center = json.loads(Path(plan['points'][0]['evaluation']).read_text())
        coordinates = np.asarray(center['coordinates_Bohr'], dtype=np.float64)
        if plan['numerical_settings'] != native.NATIVE_SETTINGS:
            raise ValueError('The original native numerical settings changed')
        builder = checkpoints.recording_builder(dft._build_smd_mf, case / 'SCF_evaluations', owner)
        counter = 0
        def solve(coords):
            nonlocal counter
            gc.collect()
            cp.get_default_memory_pool().free_all_blocks()
            folder = case / 'points' / f'{counter:04d}'
            counter += 1
            folder.mkdir(parents=True)
            mol = gto.M(atom=list(zip(symbols, coords.tolist())), unit='Bohr',
                        basis='def2-svp', charge=0, spin=0, verbose=0, max_memory=6000)
            mf = builder(mol, 'r2scan', 'acetonitrile', 'gpu', disp='d4')
            mf.conv_tol = 1e-12
            mf.conv_tol_grad = 1e-8
            mf.max_cycle = 200
            mf.grids.level = 3
            mf.callback = lambda e: atomic_write_json(folder / 'progress.json',
                dict(cycle=int(e['cycle']) + 1, E_Ha=float(e['e_tot']), updated=time.time()))
            mf.kernel()
            if not mf.converged:
                raise ValueError('Corrected control SCF did not converge')
            spin = dft._orbital_diag(mf, 0)
            if abs(spin.get('spin_contam', math.inf)) > .1:
                raise ValueError('Unchanged spin contamination gate failed')
            obj = mf.nuc_grad_method()
            obj.grid_response = True
            value = obj.kernel()
            gradient = np.asarray(value.get() if hasattr(value, 'get') else value)
            if not np.isfinite(gradient).all():
                raise ValueError('Nonfinite corrected polymer gradient')
            row = dict(coordinates_Bohr=coords.tolist(), gradient_Ha_Bohr=gradient.tolist(),
                       E_Ha=float(mf.e_tot), spin=spin, SCF_converged=True,
                       SCF_checkpoint=str(mf.chkfile),SCF_state_restored=False,native_G_qualified=False)
            atomic_write_json(folder / 'evaluation.json', row)
            return mf, gradient, row, folder
        with helper.corrected_screening(sources / 'kernel.cu', case / 'screening_implementation'), \
                monitor.MemoryMonitor() as memory, native.execution_budget() as policy, \
                adapter.compact_smd_execution(case / 'execution', 28), \
                model.smooth_cds(case / 'model', ROOT / 'results/validation/smooth_CDS_2026-10-25/pilot'), \
                finer_cds_hessian(1e-5):
            saved_components=None
            if plan.get('resume_components'):
                frame,saved_components=restart.read_frame(plan['resume_components']['source_case'])
                if frame!=plan['resume_components']['frame']:raise ValueError('Restart frame changed')
                mol=gto.M(atom=list(zip(symbols,coordinates.tolist())),unit='Bohr',basis='def2-svp',charge=0,spin=0,verbose=0,max_memory=6000)
                mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
                mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3
                restart.restore_mf(mf,frame,saved_components,case/'restart_energy.json')
                obj=mf.nuc_grad_method();obj.grid_response=True;value=obj.kernel()
                gradient=np.asarray(value.get() if hasattr(value,'get') else value)
                center_folder=case/'points/0000';center_folder.mkdir(parents=True);counter=1
                atomic_write_json(center_folder/'evaluation.json',dict(coordinates_Bohr=coordinates.tolist(),
                    gradient_Ha_Bohr=gradient.tolist(),E_Ha=float(mf.e_tot),SCF_converged=True,
                    SCF_state_restored=True,source_frame=frame,native_G_qualified=False))
            else:mf, gradient, _, center_folder = solve(coordinates)
            maximum_gradient = float(abs(gradient).max())
            if maximum_gradient > native.NATIVE_SETTINGS['gradient_limit_Ha_Bohr']:
                raise ValueError('Unchanged stationarity gate failed at the retained geometry')
            columns = native.fd_columns(gradient, symbols)
            required = (plan['failed_column']['atom'], plan['failed_column']['axis'])
            columns = sorted(set(columns) | {required})
            atomic_write_json(case / 'columns.json', dict(columns=columns, steps_Bohr=[.0005, .001]))
            mf.cphf_grids = mf.grids.copy()
            if saved_components is None:
                np.savez(case / 'grids.npz', coords=cp.asnumpy(mf.grids.coords), weights=cp.asnumpy(mf.grids.weights))
            else:shutil.copyfile(Path(plan['resume_components']['source_case'])/'grids.npz',case/'grids.npz')
            identity = dict(center_evaluation=str(center_folder / 'evaluation.json'),
                center_evaluation_sha256=hashlib.sha256((center_folder / 'evaluation.json').read_bytes()).hexdigest(),
                grids_sha256=hashlib.sha256((case / 'grids.npz').read_bytes()).hexdigest(),
                source_namespace=plan['original_namespace'], screening_variant='corrected')
            hobj = mf.Hessian()
            hobj.grid_response = True
            hobj.max_cycle = 100
            from gpu4pyscf.hessian import rhf as response_rhf
            from contextlib import nullcontext
            reuse=nullcontext() if saved_components is None else restart.reuse_components(hobj,saved_components,case/'component_reuse')
            response_reuse=(nullcontext() if response_restart is None else response_restart.checked_blocks(adapter,
                plan['resume_response']['blocks'],case/'response_reuse',expected_identity=plan['resume_response']['identity'],
                source_h1=plan['resume_response']['full_H1'],current_h1=case/'solution_H1/identity.json'))
            with responses.response_memory_cap(response_rhf, plan['response_available_memory_GiB'] * 2**30), reuse, \
                    restart.capture_solution_h1(hobj,case/'solution_H1',identity), \
                    response_reuse, \
                    responses.trace_solver(adapter, case / 'trace', checkpoint_root=case / 'response_blocks',
                                           checkpoint_identity=identity):
                value = hobj.kernel()
            tensor = np.asarray(value.get() if hasattr(value, 'get') else value)
            matrix = tensor.transpose(0, 2, 1, 3).reshape(coordinates.size, coordinates.size)
            np.save(case / 'Hessian.npy', matrix, allow_pickle=False)
            asymmetry = float(abs(matrix - matrix.T).max())
            atomic_write_json(case / 'Hessian_diagnostics.json', dict(maximum_asymmetry=asymmetry,
                              all_finite=bool(np.isfinite(matrix).all()), allocation_policy=policy))
            if not np.isfinite(matrix).all() or asymmetry > 1e-6:
                raise ValueError('Unchanged full-Hessian finite/symmetry gate failed')
            del mf, hobj, value, tensor
            gc.collect()
            cp.get_default_memory_pool().free_all_blocks()
            checks = []
            for atom, axis in columns:
                for step in (.0005, .001):
                    displaced = []
                    for sign in (-1, 1):
                        coords = coordinates.copy()
                        coords[atom, axis] += sign * step
                        mf, value, _, folder = solve(coords)
                        displaced.append((value, str(folder / 'evaluation.json')))
                        del mf
                    finite = (displaced[1][0] - displaced[0][0]).ravel() / (2 * step)
                    error = float(np.linalg.norm(finite - matrix[:, 3 * atom + axis]) / max(np.linalg.norm(finite), 1e-14))
                    checks.append(dict(atom=atom, axis=axis, step_Bohr=step, relative_error=error,
                                       minus=displaced[0][1], plus=displaced[1][1]))
                    atomic_write_json(case / 'derivatives.json', dict(checks=checks, complete=False))
                    if not math.isfinite(error) or error > .01:
                        raise ValueError('Unchanged sampled derivative gate failed')
            atomic_write_json(case / 'derivatives.json', dict(checks=checks, complete=True))
        if not 0 < memory.peak <= 31744:
            raise ValueError('Corrected full-Hessian GPU ceiling failed')
        reference = plan['checked_native_reference']
        partition = json.loads(Path(reference['partition']).read_text())
        atoms = read(reference['geometry'])
        nuclear = qrrho.nuclear_thermals(matrix, atoms.positions, atoms.numbers, partition['symmetry_number'])
        if min(nuclear['frequencies_cm']) < -5:
            raise ValueError('Unchanged internal-mode minimum gate failed')
        original = np.load(reference['Hessian'], allow_pickle=False)
        comparison = dict(maximum_Hessian_difference_Ha_Bohr2=float(abs(matrix - original).max()),
            relative_Hessian_difference=float(np.linalg.norm(matrix - original) / np.linalg.norm(original)),
            nuclear_difference_eV={k:nuclear['nuclear_by_cutoff_eV'][k] - partition['nuclear']['nuclear_by_cutoff_eV'][k]
                                   for k in ('25', '50', '100')},
            source_geometric_symmetry_number=partition['symmetry_number'], native_G_qualified=False)
        atomic_write_json(case / 'nuclear.json', nuclear)
        atomic_write_json(case / 'reference_comparison.json', comparison)
        check_sources()
        record.update(status='screening_polymer_Hessian_control_finished_review_pending',
            observed_peak_GPU_MiB=memory.peak, maximum_gradient_Ha_Bohr=maximum_gradient,
            maximum_sampled_column_error=max(c['relative_error'] for c in checks),
            minimum_frequency_cm=min(nuclear['frequencies_cm']), comparison=comparison)
    except BaseException as error:
        record.update(status='screening_polymer_Hessian_control_failed_preserved', error=repr(error))
        raise
    finally:
        record['finished'] = time.time()
        atomic_write_json(case / 'status.json', record)


if __name__ == '__main__':
    main()
