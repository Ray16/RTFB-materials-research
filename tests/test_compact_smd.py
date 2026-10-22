"""The in-place solvent matrix must retain the original equation and inputs."""

import numpy as np
import pytest
from types import SimpleNamespace
from redox.qm.compact_smd import compact_build


def build(self):
    S, D, A, f_epsilon = self.S, self.D, self.A, self.factor
    DA = D * A
    DAS = np.dot(DA, S)
    K = S - f_epsilon / (2.0 * np.pi) * DAS
    self.K = K.copy()
    K = None
    return self


def test_reused_matrix_matches_reference_and_keeps_physical_matrices_unchanged():
    rng = np.random.default_rng(91)
    S, D = rng.normal(size=(2, 17, 17))
    A = rng.uniform(0.2, 1.0, 17)
    state = SimpleNamespace(S=S.copy(), D=D.copy(), A=A, factor=0.9)
    expected = build(state).K.copy()
    updated, source = compact_build(build)
    actual = updated(state).K
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)
    np.testing.assert_array_equal(state.S, S)
    np.testing.assert_array_equal(state.D, D)
    assert "DAS *= " in source


def unrelated(self):
    K = self.S + self.D
    return K


def test_unreviewed_library_expression_is_refused():
    with pytest.raises(ValueError, match="differs"):
        compact_build(unrelated)
