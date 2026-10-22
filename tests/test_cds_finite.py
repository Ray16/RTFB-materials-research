"""Conservative nuclear derivatives preserve signs, coordinates and ordering."""

import numpy as np
import pytest

from redox.qm.cds_finite import energy_gradient, gradient_hessian


def test_energy_gradient_and_hessian_reproduce_a_known_coupled_potential():
    rng = np.random.default_rng(33)
    a = rng.normal(size=(6, 6))
    a = a + a.T
    x = rng.normal(size=(2, 3))
    saved = x.copy()
    energy = lambda y: 0.5 * y.ravel() @ a @ y.ravel()
    gradient = lambda y: (a @ y.ravel()).reshape(2, 3)
    assert np.max(abs(energy_gradient(x, energy) - gradient(x))) < 1e-9
    tensor = gradient_hessian(x, gradient)
    assert np.max(abs(tensor.transpose(0, 2, 1, 3).reshape(6, 6) - a)) < 1e-9
    assert np.array_equal(x, saved)


@pytest.mark.parametrize("step", [0, -1, float("nan")])
def test_bad_energy_derivative_step_is_rejected(step):
    with pytest.raises(ValueError):
        energy_gradient(np.zeros((1, 3)), lambda x: 1.0, step)
