"""The scoped numerical refinement must not become an unexplained gate bypass."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("saved_polymer_energy_qualification",
    ROOT / "scripts/pipeline/qualify_saved_polymer_energy.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture
def evidence():
    coarse = dict(passed=False, acceptance_Ha_Bohr=1e-6, checks=[
        dict(atom=a, axis=1, step_Bohr=h, analytic=0., numerical=(1e-8 if h==.0005 else 2e-5),
             absolute_error=(1e-8 if h==.0005 else 2e-5)) for a in (29, 100) for h in (.0005, .001)])
    fine = dict(diagnostic_criteria_passed=True, derivative_criteria_passed=True,
        checks=[dict(atom=a, axis=1, step_Bohr=h, analytic=0., numeric=1e-7, error_Ha_Bohr=1e-7)
                for a in (29, 100) for h in (.0001, .00005, .00001)])
    surface = dict(maximum_remaining_difference_Ha_Bohr=1e-8, maximum_library_energy_difference_Ha=0.)
    return coarse, fine, surface


def test_four_step_refinement_keeps_the_original_failure_and_tolerance(evidence):
    result = m.qualification_checks(*evidence)
    assert result["passed"] and result["original_coarse_check_passed"] is False
    assert result["acceptance_Ha_Bohr"] == 1e-6
    assert len(result["checks"]) == 8
    assert evidence[0]["passed"] is False and len(evidence[0]["checks"]) == 4


@pytest.mark.parametrize("bad", ["unexplained", "nonfinite", "fine_failure", "missing_atom", "missing_step", "force_error"])
def test_incomplete_or_unexplained_refinement_is_refused(evidence, bad):
    coarse, fine, surface = evidence
    if bad == "unexplained":
        surface["maximum_remaining_difference_Ha_Bohr"] = 2e-5
    elif bad == "nonfinite":
        surface["maximum_remaining_difference_Ha_Bohr"] = float("nan")
    elif bad == "fine_failure":
        fine["diagnostic_criteria_passed"] = False
    elif bad == "missing_atom":
        fine["checks"] = [r for r in fine["checks"] if r["atom"] != 100]
    elif bad == "missing_step":
        fine["checks"].pop()
    else:
        fine["checks"][0]["error_Ha_Bohr"] = 1.1e-6
    with pytest.raises(ValueError):
        m.qualification_checks(coarse, fine, surface)
