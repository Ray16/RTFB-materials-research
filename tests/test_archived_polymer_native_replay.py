import copy
import importlib.util
from pathlib import Path

import numpy as np
import pytest

SOURCE = Path(__file__).resolve().parents[1]/'scripts/validation/replay_archived_polymer_native.py'
spec = importlib.util.spec_from_file_location('polymer_replay', SOURCE)
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def harmonic_points():
    coords = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    gradient = np.arange(1., 10.).reshape(3, 3)*1e-6
    symbols = ['C', 'O', 'N']
    required = {(0, 1), (1, 2), (2, 1), (2, 2)}
    points, checks, index = {}, [], 1

    def point(xyz, grad):
        return dict(coordinates_Bohr=xyz, gradient_Ha_Bohr=grad,
                    SCF_converged=True, spin=dict(spin_contam=0.), experimental_method_id='method')

    points[0] = point(coords.copy(), gradient.copy())
    for atom, axis in sorted(required):
        for step in (.0005, .001):
            for sign in (-1, 1):
                xyz, grad = coords.copy(), gradient.copy()
                xyz[atom, axis] += sign*step
                grad[atom, axis] += sign*step
                points[index] = point(xyz, grad)
                index += 1
            checks.append(dict(atom=atom, axis=axis, step_Bohr=step, relative_error=0.))
    native = dict(required_columns=sorted(required), derivative_checks=checks)
    reader = lambda path: points[int(path.parent.name)]
    return native, points, gradient, symbols, coords, reader


def test_raw_displaced_gradients_reconstruct_analytic_columns():
    native, _, gradient, symbols, coords, reader = harmonic_points()
    checks, columns, total = replay.finite_difference_checks(
        Path('/raw'), native, np.eye(9), gradient, symbols, coords, 'method', reader)
    assert max(checks) < 1e-10
    assert len(columns) == 4 and total == 17


@pytest.mark.parametrize('failure', ['geometry', 'gradient', 'hessian', 'method', 'missing_column', 'nan_spin'])
def test_rejects_physically_inconsistent_displacements(failure):
    native, points, gradient, symbols, coords, reader = harmonic_points()
    matrix = np.eye(9)
    if failure == 'geometry': points[1]['coordinates_Bohr'][0, 0] += .1
    if failure == 'gradient': points[1]['gradient_Ha_Bohr'][0, 0] = np.nan
    if failure == 'hessian': matrix *= 1.02
    if failure == 'method': points[1]['experimental_method_id'] = 'other'
    if failure == 'missing_column': native['required_columns'].pop()
    if failure == 'nan_spin': points[1]['spin']['spin_contam'] = float('nan')
    with pytest.raises(ValueError):
        replay.finite_difference_checks(Path('/raw'), native, matrix, gradient,
                                        symbols, coords, 'method', reader)


def test_rigid_motion_projection_is_translation_and_rotation_invariant():
    coords = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    masses = np.full(3, 12.)
    first, rank = replay.projected_modes(np.eye(9), coords, masses)
    rotation = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
    second, other_rank = replay.projected_modes(np.eye(9), coords@rotation+100., masses)
    assert rank == other_rank == 6
    assert len(first) == 3
    assert np.all(first > 0)
    assert np.allclose(first, second, rtol=0, atol=1e-8)


def test_explicit_configuration_aliases_cannot_contradict_each_other():
    assert replay.configuration(dict(stereoisomer='all_same')) == 'all_same'
    with pytest.raises(ValueError):
        replay.configuration(dict(configuration_label='alternating', stereoisomer='all_same'))
    with pytest.raises(ValueError):
        replay.configuration(dict(key='all_same__occ_2__conf001'))
