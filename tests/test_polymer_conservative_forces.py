import importlib.util
import json
from pathlib import Path

import pytest


def driver():
    path = Path(__file__).resolve().parents[1] / "scripts/pipeline/resume_polymer_forces.py"
    spec = importlib.util.spec_from_file_location("checked_force_continuation", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_reject_unstable_or_mismatched_seed_derivative(tmp_path):
    module = driver()
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\n\nH 0 0 0\n")
    report = dict(geometry_sha256=module.api.digest(xyz), steps_Bohr=[1e-4, 5e-5],
                  maximum_step_disagreement=1e-8, full_energy_difference_checks=[dict(absolute_error=2e-8)])
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    assert module.checked_diagnostic(path, xyz) == report
    for key, value in [("maximum_step_disagreement", 2e-6), ("geometry_sha256", "other"),
                       ("steps_Bohr", [1e-3, 1e-4]),
                       ("full_energy_difference_checks", [dict(absolute_error=2e-6)])]:
        path.write_text(json.dumps(dict(report, **{key: value})))
        with pytest.raises(ValueError):
            module.checked_diagnostic(path, xyz)
