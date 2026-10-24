"""Independent replay must reject inconsistent or incomplete force evidence."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fine_force_review",
    ROOT / "scripts/validation/review_polymer_gradient_profile.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture
def harmonic():
    plan = dict(acceptance=dict(maximum_total_derivative_error_Ha_Bohr=1e-6,
        finest_pair_derivative_difference_Ha_Bohr=1e-6, finest_steps_Bohr=[5e-5, 1e-5]))
    records = dict(baseline=dict(coordinates_Bohr=[[0., 0., 0.]], gradient_Ha_Bohr=[[0.2, 0., 0.]]))
    differences = []
    for h in (1e-4, 5e-5, 1e-5):
        for sign in (-1, 1):
            x = sign*h
            records[f"atom0_axis0_h{h}_sign{sign}"] = dict(
                coordinates_Bohr=[[x, 0., 0.]], E_Ha=0.2*x+0.5*x*x, SCF_converged=True)
        differences.append(dict(atom=0, axis=0, step_Bohr=h, analytic=0.2, numerical=0.2))
    return plan, records, differences


def test_harmonic_energy_pairs_reconstruct_the_analytic_force(harmonic):
    result = m.compare(*harmonic)
    assert result["derivative_criteria_passed"]
    assert result["maximum_total_derivative_error_Ha_Bohr"] < 1e-12


def test_an_inconsistent_analytic_force_fails_the_same_tolerance(harmonic):
    plan, records, differences = harmonic
    records["baseline"]["gradient_Ha_Bohr"][0][0] += 1e-5
    for r in differences:
        r["analytic"] += 1e-5
    assert not m.compare(plan, records, differences)["derivative_criteria_passed"]


@pytest.mark.parametrize("change", ["unconverged", "other_coordinate", "duplicate"])
def test_invalid_pairs_cannot_be_admitted_by_a_good_reported_derivative(harmonic, change):
    plan, records, differences = harmonic
    if change == "unconverged":
        records["atom0_axis0_h0.0001_sign1"]["SCF_converged"] = False
    elif change == "other_coordinate":
        records["atom0_axis0_h0.0001_sign1"]["coordinates_Bohr"][0][2] = 1e-4
    else:
        differences.append(dict(differences[0]))
    with pytest.raises(ValueError):
        m.compare(plan, records, differences)


def test_a_live_diagnostic_is_not_reviewed_as_complete(tmp_path):
    (tmp_path / "plan.json").write_text("{}")
    (tmp_path / "status.json").write_text(json.dumps(dict(status="fine_step_diagnostic_running")))
    with pytest.raises(ValueError, match="not complete"):
        m.review(tmp_path)
