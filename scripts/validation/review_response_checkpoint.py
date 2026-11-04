#!/usr/bin/env python3
"""Review a completed response execution control without admitting a Hessian/G.

Checks source/owner/resource closure and saved FP64 arrays. Residual norms are
source-bound solver records, not a new evaluation of the GPU response equation.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def review(namespace, placement):
    spec = verify(namespace)
    case = Path(spec['case'])
    used = []
    def read(path):
        used.append(path)
        return json.loads(path.read_text())
    plan = read(case / 'plan.json')
    stage = read(case / 'status.json')
    completion = read(placement / 'completion.json')
    launch = read(placement / 'launch.json')
    owner = read(case / 'residency_guard/owner.json')
    residency = read(case / 'residency_guard/state.json')
    host = read(case / 'host_guard/status.json')
    runtime = read(case / 'tensor_runtime/status.json')
    energy = read(case / 'energy_reproduction.json')
    trace = read(case / 'trace/trace.json')
    solver = read(case / 'execution/response_solver.json')
    summary = read(case / 'response_summary.json')
    definition = importlib.util.spec_from_file_location('response_closed_owner', ROOT / 'scripts/pipeline/observe_archived_native_namespace.py')
    checker = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(checker)
    if (plan['role'] != 'orbital_response_diagnostic' or not plan['checkpoint_response']
            or stage['status'] != 'orbital_response_diagnostic_finished_no_native_admission'
            or checker.alive(stage['owner']) is not False
            or completion['returncode'] != 0 or completion['budget_stop'] is not None
            or completion['namespace_sha256'] != sha256_file(namespace)
            or not owner['guard_stopped'] or residency['foreign_resident_PIDs']
            or not residency['own_residency_observed']
            or host['status'] != 'effective_memory_guard_finished_profile_review_pending'
            or runtime['status'] != 'task_returned_under_checked_cutensor_runtime'
            or host['maximum_RSS_GiB'] != launch['effective_maximum_host_RSS_GiB']
            or max(host['peak_tree_RSS_bytes'], host['self_VmHWM_bytes']) > host['maximum_RSS_GiB'] * 2**30
            or host['free_RAM_bytes'] < 32 * 2**30
            or not 0 < stage['observed_peak_GPU_MiB'] <= plan['maximum_peak_GPU_MiB']
            or runtime['program'] != str(Path(spec['sources']) / 'worker.py')):
        raise ValueError('Completed exact response control and normal guards required')
    if (not stage['owner']['pid'] == owner['pid'] == host['pid'] == runtime['pid']
            or not stage['owner']['process_start'] == owner['process_start'] == host['process_start']
            or not owner['CPU_affinity'] == host['CPU_affinity'] == runtime['CPU_affinity'] == launch['CPU_affinity']
            or len(owner['CPU_affinity']) != 2):
        raise ValueError('Response ownership or allocation chain differs')
    natom = len(plan['checkpoint']['identity']['coordinates_Bohr'])
    norb, nmo = plan['checkpoint']['identity']['mo_coeff']['shape']
    nocc = summary['arrays']['mo1']['shape'][-1]
    calls = solver['calls']
    if (solver['status'] != 'response_equations_solved_full_physical_QC_separate'
            or solver['equations_changed'] or solver['level_shift_added']
            or trace['asynchronous_C_watchdog_used'] or trace['solver_replaced']
            or len(calls) != len(trace['calls']) or len(calls) != math.ceil(natom / 16)
            or sum(c['shape'][0] for c in calls) != 3 * natom):
        raise ValueError('Declared full response coverage and unchanged equation required')
    def array_record(record, expected):
        path = Path(record['path'])
        if not path.is_relative_to(case) or sha256_file(path) != record['file_sha256']:
            raise ValueError('Missing or changed response checkpoint')
        array = np.load(path, mmap_mode='r', allow_pickle=False)
        if (array.dtype != np.float64 or array.shape != tuple(expected)
                or list(array.shape) != record['shape'] or not np.isfinite(array).all()
                or hashlib.sha256(memoryview(array).cast('B')).hexdigest() != record['array_sha256']):
            raise ValueError('Response array shape, type, values or payload hash differs')
        used.append(path)
    rows = []
    for i, (call, timing) in enumerate(zip(calls, trace['calls'])):
        meta = read(case / f'response_blocks/block_{i:03d}.json')
        identity = meta['identity']
        if (timing['status'] != 'returned' or timing['index'] != i
                or timing['shape'] != call['shape'] or call['shape'][1] != nmo * nocc
                or call['shape'][0] != 3 * min(16, natom - 16 * i)
                or call['precision'] != 'float64' or call['solver_info'] != 0
                or call['level_shift_added'] or not math.isfinite(call['true_residual_norm'])
                or not 0 <= call['true_residual_norm'] <= call['solver_absolute_target']
                or call['solver_absolute_target'] != .01 * call['requested_absolute_tolerance']
                or identity['source_namespace'] != plan['source_namespace']
                or identity['checkpoint'] != plan['checkpoint']
                or sha256_file(Path(identity['recreated_grids']['path'])) != identity['recreated_grids']['sha256']
                or meta['kind'] != 'raw_Krylov_solution_before_occupied_block_reconstruction'
                or meta['native_G_qualified'] or len(meta['rhs_sha256']) != 64
                or meta['checkpoint'] != timing['checkpoint']):
            raise ValueError('Response block frame, shape or source-bound solver residual failed')
        source = Path(spec['archive']) / 'scripts/pipeline/residual_response_solver.py'
        if sha256_file(source) != meta['original_solver_source_sha256']:
            raise ValueError('Original solver implementation changed')
        array_record(meta['checkpoint'], call['shape'])
        rows.append(dict(index=i, shape=call['shape'], producer_residual_norm=call['true_residual_norm'],
                         producer_absolute_target=call['solver_absolute_target']))
    for name, shape in [('mo1', (natom, 3, norb, nocc)), ('e1', (natom, 3, nocc, nocc))]:
        array_record(summary['arrays'][name]['checkpoint'], shape)
    difference = energy['recreated_energy_Ha'] - energy['saved_energy_Ha']
    if (not math.isfinite(difference) or abs(difference) > 1e-8
            or difference != energy['difference_Ha'] or energy['original_grid_arrays_available']
            or summary['original_grid_arrays_available'] or summary['native_G_qualified']
            or stage['native_G_qualified']):
        raise ValueError('Fixed-density energy or explicitly limited restart scope differs')
    used.extend([namespace, case / 'recreated_grids.npz'])
    hashes = {str(p): sha256_file(p) for p in used}
    verify(namespace)
    return dict(status='checkpointed_response_execution_independently_checked_no_native_admission',
                namespace=str(namespace), completed_blocks=len(rows), Cartesian_RHS=3 * natom,
                blocks=rows, energy_reproduction_difference_Ha=difference,
                residuals_independently_reevaluated=False, residual_scope='source-bound producer original-equation checks',
                original_grid_arrays_available=False, exact_Hessian_restart_qualified=False,
                corrected_screening_implementation_used=False, native_G_qualified=False,
                material_admission=False, peak_GPU_MiB=stage['observed_peak_GPU_MiB'],
                peak_host_RSS_bytes=max(host['peak_tree_RSS_bytes'], host['self_VmHWM_bytes']),
                source_sha256=hashes, producer_sha256=sha256_file(__file__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('spec', 'placement', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise FileExistsError('Preserve earlier response reviews')
    result = review(args.spec.resolve(), args.placement.resolve())
    atomic_write_json(args.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'completed_blocks', 'Cartesian_RHS', 'peak_GPU_MiB', 'native_G_qualified')}))
