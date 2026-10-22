"""A cached derivative repair must preserve its exact component identity."""

import json

import numpy as np
import pytest

from redox.qm import hessian_resume as resume
from redox.workflows import polymer_campaign as api


@pytest.fixture
def spec(tmp_path):
    base = np.eye(6)
    coarse = np.zeros((6, 6))
    coarse[0, 1] = 0.001
    arrays = {
        "source": base + coarse,
        "coarse": coarse,
        "fine": np.zeros((6, 6)),
        "refined": base,
    }
    for name, array in arrays.items():
        np.save(tmp_path / (name + ".npy"), array)
    reference = tmp_path / "controls.json"
    reference.write_text(
        json.dumps(
            dict(
                status="two_full_recomputation_replays_passed",
                controls=[
                    dict(
                        mult=m,
                        input_sha256={
                            str(tmp_path / "source.npy"): api.digest(
                                tmp_path / "source.npy"
                            )
                        },
                        maximum_matrix_difference=0.0,
                        relative_matrix_error=0.0,
                        maximum_thermal_difference_eV=0.0,
                    )
                    for m in (1, 2)
                ],
            )
        )
    )
    result = dict(
        source_hessian=str(tmp_path / "source.npy"),
        coarse_CDS=str(tmp_path / "coarse.npy"),
        fine_CDS=str(tmp_path / "fine.npy"),
        refined_hessian=str(tmp_path / "refined.npy"),
        reconstruction_replay=str(reference),
    )
    result["input_sha256"] = {str(p): api.digest(p) for p in tmp_path.iterdir()}
    return result


def test_additive_reconstruction_removes_only_the_declared_component(spec):
    matrix, indices = resume.reconstructed_matrix(spec)
    assert np.array_equal(matrix, np.eye(6)) and indices == [0, 1]
    assert np.load(spec["source_hessian"])[0, 1] == 0.001


def test_changed_component_is_not_accepted_from_a_done_marker(spec):
    np.save(spec["fine_CDS"], np.eye(6))
    with pytest.raises(RuntimeError, match="changed"):
        resume.reconstructed_matrix(spec)


def test_oversized_replay_error_cannot_be_hidden_by_good_target_symmetry(spec):
    from pathlib import Path

    path = Path(spec["reconstruction_replay"])
    value = json.loads(path.read_text())
    value["controls"][1]["maximum_matrix_difference"] = 0.01
    path.write_text(json.dumps(value))
    spec["input_sha256"][str(path)] = api.digest(path)
    with pytest.raises(ValueError, match="fixed tolerance"):
        resume.reconstructed_matrix(spec)
