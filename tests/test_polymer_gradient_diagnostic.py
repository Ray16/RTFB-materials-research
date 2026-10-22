"""A derivative diagnostic can reuse only identical, converged evaluations."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def diagnostic():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts/validation/polymer_gradient_diagnostic.py"
    )
    spec = importlib.util.spec_from_file_location("gradient_replay_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_probe_coordinates_are_explicit_and_validated(diagnostic):
    assert diagnostic.parse_probe("7:1") == (7, 1)
    for value in ("-1:0", "0:3"):
        with pytest.raises(Exception):
            diagnostic.parse_probe(value)


@pytest.mark.parametrize("change", ["coordinates", "convergence", "energy"])
def test_changed_or_unconverged_records_are_refused(diagnostic, tmp_path, change):
    record = dict(E_Ha=-2.0, SCF_converged=True, coordinates_Bohr=[[0.0, 0.0, 0.0]])
    if change == "coordinates":
        record["coordinates_Bohr"][0][1] = 0.01
    if change == "convergence":
        record["SCF_converged"] = False
    if change == "energy":
        record["E_Ha"] = float("nan")
    (tmp_path / "point.json").write_text(json.dumps(record))
    with pytest.raises(ValueError):
        diagnostic.reusable_record(tmp_path, "point", np.zeros((1, 3)))


def test_valid_record_retains_its_original_bytes_and_hash(diagnostic, tmp_path):
    record = dict(E_Ha=-2.0, SCF_converged=True, coordinates_Bohr=[[0.0, 0.0, 0.0]])
    path = tmp_path / "point.json"
    path.write_text(json.dumps(record))
    original = path.read_bytes()
    replay = diagnostic.reusable_record(tmp_path, "point", np.zeros((1, 3)))
    assert replay["reused_sha256"] == diagnostic.file_hash(path)
    assert path.read_bytes() == original
