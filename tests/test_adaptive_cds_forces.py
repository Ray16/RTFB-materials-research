import importlib.util
from pathlib import Path

import numpy as np
import pytest


def load():
    path = Path(__file__).resolve().parents[1] / "scripts/pipeline/resume_polymer_forces_adaptive.py"
    spec = importlib.util.spec_from_file_location("adaptive_cds_driver_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_adaptive_derivative_refines_and_keeps_tolerance():
    module = load()
    point = np.array([[.1, .2, .3]])

    def energy(coords):
        # Rapidly varying smooth term forces refinement of x only.
        return float(np.sum(coords ** 2) + 1e-6 * np.sin(10000 * (coords[0, 0] - .1)))

    result, evidence = module.adaptive_gradient(point, energy)
    np.testing.assert_allclose(result, [[.21, .4, .6]], atol=1e-6, rtol=0)
    assert len(evidence[0]["attempts"]) > 2
    assert all(c["final_step_difference"] <= 1e-6 for c in evidence)
    assert len(evidence[1]["attempts"]) == 2


def test_adaptive_derivative_rejects_nonfinite_and_persistent_jump():
    module = load()
    with pytest.raises(ValueError, match="Nonfinite"):
        module.adaptive_gradient(np.zeros((1, 3)), lambda coords: np.nan)
    with pytest.raises(ValueError, match="remains step inconsistent"):
        module.adaptive_gradient(np.zeros((1, 3)), lambda coords: float(coords[0, 0] > 0))
