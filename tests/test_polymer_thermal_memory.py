"""Check strided solvent contractions and fail-closed full-matrix replay gates."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def adapter():
    path = Path(__file__).parents[1] / "scripts/pipeline/resume_polymer_thermal.py"
    spec = importlib.util.spec_from_file_location("thermal_memory_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("atoms,points", [(1, 1), (4, 37), (11, 129)])
def test_atom_blocks_reproduce_strided_analytic_area_contraction(adapter, atoms, points):
    rng = np.random.default_rng(196)
    tensor = rng.normal(size=(points, atoms, atoms, 3, 3)).transpose(1, 2, 3, 4, 0)
    left, right = rng.normal(size=(2, points))
    expected = tensor @ (left * right)
    actual = adapter.area_contraction(tensor, left, right)
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("problem", ["failed", "changed_evidence", "matrix_error"])
def test_replay_gate_rejects_failures_changed_outputs_and_bad_matrix(adapter, tmp_path, monkeypatch, problem):
    monkeypatch.setattr(adapter, "CAMPAIGN", tmp_path)
    evidence = tmp_path / "matrix.npy"
    evidence.write_bytes(b"immutable matrix")
    expected = adapter.api.digest(evidence)
    for mult in (1, 2):
        path = tmp_path / "controls" / str(mult) / "status.json"
        path.parent.mkdir(parents=True)
        row = dict(status="full_native_replay_passed", native_gates_passed=True,
                   evidence_sha256={str(evidence): expected},
                   errors={name: 1e-10 for name in ("relative_matrix_error", "maximum_matrix_difference", "maximum_thermal_difference_eV")})
        if mult == 2 and problem == "failed":
            row["status"] = "failed_preserved"
        if mult == 2 and problem == "matrix_error":
            row["errors"]["relative_matrix_error"] = 1e-3
        path.write_text(json.dumps(row))
    if problem == "changed_evidence":
        evidence.write_bytes(b"changed")
    with pytest.raises((RuntimeError, ValueError)):
        adapter.checked_replays()
