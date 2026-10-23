"""Alternative-basis controls cannot be promoted into primary certificates."""
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("package_auxiliary", ROOT / "scripts/release/review_package_auxiliary.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def surface():
    return dict(primary_native_certificate_claimed=False, identity=dict(protocol=dict(basis="def2-svpd")),
        status="derivative_checked_solution_thermal", stationary_minimum=True, grid_response=True,
        gradient_max=1e-5, maximum_asymmetry=1e-10,
        derivative_checks=[dict(relative_error=1e-4) for _ in range(4)])


def test_alternative_control_retains_its_separate_scope():
    MODULE.validate_alternative_surface(surface())


@pytest.mark.parametrize("change", [
    lambda x: x.update(primary_native_certificate_claimed=True),
    lambda x: x["identity"]["protocol"].update(basis="def2-svp"),
    lambda x: x.update(gradient_max=.001),
    lambda x: x.update(derivative_checks=[dict(relative_error=.1)]*4),
])
def test_alternative_control_refuses_scope_or_gate_changes(change):
    value = surface()
    change(value)
    with pytest.raises(ValueError, match="Alternative basis"):
        MODULE.validate_alternative_surface(value)
