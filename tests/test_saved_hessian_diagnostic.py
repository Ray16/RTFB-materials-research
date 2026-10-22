"""Component extraction keeps Hessian coordinate ordering and does not mask skew."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def diagnostic():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts/validation/diagnose_saved_hessian.py"
    )
    spec = importlib.util.spec_from_file_location("saved_H_diagnostic", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_finite_gradient_component_recovers_full_known_tensor(diagnostic):
    rng = np.random.default_rng(67)
    a = rng.normal(size=(6, 6))
    hessian = a + a.T
    x = rng.normal(size=(2, 3))
    result = diagnostic.finite_gradient_hessian(
        x, lambda y: (hessian @ y.ravel()).reshape(2, 3), 1e-4
    )
    assert np.max(abs(result - hessian)) < 1e-10


def test_summary_exposes_asymmetry_instead_of_symmetrizing_it(diagnostic):
    matrix = np.eye(6)
    matrix[0, 3] = 0.01
    report = diagnostic.matrix_summary(matrix)
    assert report["maximum_asymmetry"] == 0.01 and report["worst_indices"] == [0, 3]
    assert matrix[3, 0] == 0
