"""Physical invariants of the Cartesian thermochemistry boundary."""
import numpy as np
import pytest

from redox.qm.vibrations import project_rigid_modes


@pytest.mark.parametrize('xyz,rank', [
    ([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]], 6),
    ([[-1., 0., 0.], [0., 0., 0.], [2., 0., 0.]], 5),
    ([[0., 0., 0.]], 3),
])
def test_rigid_displacements_are_null_and_projection_is_idempotent(xyz, rank):
    xyz = np.asarray(xyz)
    mass = np.arange(len(xyz)) + 1.
    rng = np.random.default_rng(42)
    arbitrary = rng.normal(size=(3 * len(xyz),) * 2)
    projected, actual = project_rigid_modes(arbitrary, xyz, mass)
    assert actual == rank
    for axis in np.eye(3):
        translation = np.tile(axis, len(xyz))
        rotation = np.cross(axis, xyz - np.average(xyz, axis=0, weights=mass)).ravel()
        assert np.linalg.norm(projected @ translation) < 1e-12
        assert np.linalg.norm(projected @ rotation) < 1e-12
    again, _ = project_rigid_modes(projected, xyz, mass)
    assert np.allclose(projected, again, atol=1e-12)
    assert np.allclose(projected, projected.T, atol=1e-12)


def test_stretch_and_instability_survive_projection():
    xyz = np.array([[-1., 0., 0.], [1., 0., 0.]])
    stretching = np.array([-1., 0., 0., 1., 0., 0.])
    for sign in (-1., 1.):
        physical = sign * np.outer(stretching, stretching)
        projected, rank = project_rigid_modes(physical, xyz, [1., 1.])
        assert rank == 5
        assert np.allclose(projected, physical, atol=1e-12)


def test_translation_and_coordinate_units_do_not_change_result():
    xyz = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    h = np.arange(81).reshape(9, 9)
    projected, _ = project_rigid_modes(h, xyz, [12., 1., 16.])
    moved, _ = project_rigid_modes(h, xyz * 100 + [4., -9., 6.], [12., 1., 16.])
    assert np.allclose(projected, moved, atol=1e-10)


@pytest.mark.parametrize('mass', [[0., 1.], [-1., 1.], [np.nan, 1.]])
def test_invalid_masses_are_refused(mass):
    with pytest.raises(ValueError):
        project_rigid_modes(np.eye(6), [[0., 0., 0.], [1., 0., 0.]], mass)
