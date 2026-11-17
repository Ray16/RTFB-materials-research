"""Solve the original matrix-free response equation with a checked residual.

Only the linear solver changes. The supplied response operator, orbital-energy
denominators, electronic state and precision remain untouched. SciPy GMRES uses
host vectors; evaluations of a supplied CuPy operator remain on the GPU.
"""
from contextlib import contextmanager
import math
from pathlib import Path
import time

import numpy as np


def residual_krylov(aop, b, x0=None, tol=1e-10, max_cycle=30, lindep=1e-14,
                    dot=None, callback=None, verbose=None, *, record=None):
    from scipy.sparse.linalg import LinearOperator, gmres
    if (not math.isfinite(tol) or tol <= 0 or type(max_cycle) is not int
            or max_cycle < 1 or dot is not None or callback is not None):
        raise ValueError('Positive finite residual tolerance and supported solver arguments required')
    on_gpu = type(b).__module__.split('.')[0] == 'cupy'
    if on_gpu:
        import cupy as cp
        host = cp.asnumpy
        device = cp.asarray
    else:
        host = np.asarray
        device = np.asarray
    rhs = np.asarray(host(b))
    if rhs.ndim not in (1, 2) or rhs.dtype != np.float64 or not np.isfinite(rhs).all():
        raise ValueError('Finite FP64 one- or multiple-right-hand-side response required')
    shape = rhs.shape
    operator_shape = shape if rhs.ndim == 2 else (1, rhs.size)
    evaluations = 0
    limit = 2 * max_cycle + 10
    def apply(vector):
        nonlocal evaluations
        evaluations += 1
        if evaluations > limit:
            raise RuntimeError('Bounded response operator budget exhausted')
        value = np.asarray(host(aop(device(vector.reshape(operator_shape)))))
        if value.shape != operator_shape or not np.isfinite(value).all():
            raise ValueError('Response operator returned invalid dimensions or nonfinite values')
        return vector + value.ravel()
    initial = None if x0 is None else np.asarray(host(x0)).ravel().copy()
    if initial is not None and (initial.size != rhs.size or not np.isfinite(initial).all()):
        raise ValueError('Invalid response initial guess')
    # Stronger than the supplied absolute norm threshold. No relative tolerance
    # can silently weaken the criterion for a large perturbation block.
    target = tol * .01
    restart = min(30, max_cycle, rhs.size)
    history = []
    solution, info = gmres(LinearOperator((rhs.size, rhs.size), matvec=apply, dtype=np.float64),
        rhs.ravel(), x0=initial, rtol=0., atol=target, restart=restart,
        maxiter=math.ceil(max_cycle / restart), callback=lambda value: history.append(float(value)),
        callback_type='pr_norm')
    residual = apply(solution) - rhs.ravel()
    norm = float(np.linalg.norm(residual))
    row = dict(algorithm='SciPy restarted GMRES, unchanged I+A equation',
        requested_absolute_tolerance=tol, solver_absolute_target=target,
        true_residual_norm=norm, maximum_absolute_residual=float(np.max(abs(residual))),
        rhs_norm=float(np.linalg.norm(rhs)), solver_info=int(info), operator_evaluations=evaluations,
        iterations=len(history), residual_history=history, shape=list(shape),
        precision='float64', level_shift_added=False, checked=time.time())
    if record is not None:
        record(row)
    if info != 0 or not math.isfinite(norm) or norm > target:
        raise RuntimeError('Response failed independent original-equation residual check')
    return device(solution.reshape(shape))


@contextmanager
def checked_response_solver(directory):
    from gpu4pyscf.hessian import rhf, uhf
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    calls = []
    proof = dict(status='residual_response_solver_active_unqualified', started=time.time(),
        equations_changed=False, level_shift_added=False, calls=calls,
        source_sha256={str(p): api.digest(p) for p in
            (Path(__file__).resolve(), Path(rhf.__file__), Path(uhf.__file__))})
    def record(row):
        calls.append(row)
        atomic_write_json(directory / 'response_solver.json', proof)
    def solve(*args, **kwargs):
        return residual_krylov(*args, **kwargs, record=record)
    originals = [(module, module.krylov) for module in (rhf, uhf)]
    for module, _ in originals:
        module.krylov = solve
    try:
        yield proof
        if not calls:
            raise ValueError('No response solve exercised')
        proof['status'] = 'response_equations_solved_full_physical_QC_separate'
    except BaseException as error:
        proof.update(status='response_solver_failed_preserved', error=repr(error))
        raise
    finally:
        for module, previous in originals:
            module.krylov = previous
        proof['finished'] = time.time()
        atomic_write_json(directory / 'response_solver.json', proof)


@contextmanager
def compact_smd_execution(directory, budget_GiB, occupied_cap=32):
    from redox.core.common import ROOT
    from redox.workflows import polymer_campaign as api
    parent = api.load_module('residual_solver_unchanged_bounded_operators',
        ROOT / 'scripts/pipeline/df_occupied_workspace.py')
    with parent.compact_smd_execution(directory, budget_GiB, occupied_cap) as proof:
        with checked_response_solver(directory):
            yield proof
