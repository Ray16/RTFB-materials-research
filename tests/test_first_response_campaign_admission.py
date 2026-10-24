"""Native preparation requires the exact geometry's checked response evidence."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("first_response_campaign_admission",
    ROOT / "scripts/pipeline/first_response_polymer_campaign.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture
def proof(tmp_path):
    folder = tmp_path / "proof"
    folder.mkdir()
    state = tmp_path / "state"
    state.mkdir()
    (state / "opt.xyz").write_text("exact saved geometry")
    (folder / "input.xyz").write_bytes((state / "opt.xyz").read_bytes())
    plan = dict(input_xyz=str(folder / "input.xyz"), task=dict(key="exact_case"),
                producer_sha256={}, input_sha256={})
    status = dict(status="response_resource_and_sampled_derivatives_passed", output_sha256={})
    (folder / "plan.json").write_text(json.dumps(plan))
    (folder / "status.json").write_text(json.dumps(status))
    review = dict(status="fixed_density_response_matrices_independently_reconstructed",
        maximum_relative_error=1e-5, qualified_native_Hessian=False, qualified_G=False,
        observed_peak_GPU_MiB=31000,
        source_sha256={str(folder / n): hashlib.sha256((folder/n).read_bytes()).hexdigest()
                       for n in ("plan.json", "status.json")})
    (folder / "independent_review.json").write_text(json.dumps(review))
    report = dict(state_dir=str(state), task=dict(key="exact_case"))
    return folder, report


def test_response_qualification_is_bound_to_original_geometry_and_state(proof):
    folder, report = proof
    assert len(m.response_evidence(folder, report)) == 3


@pytest.mark.parametrize("change", ["geometry", "state", "promoted_scope", "large_error", "mutated_status"])
def test_other_or_unqualified_evidence_is_not_a_native_resource_pass(proof, change):
    folder, report = proof
    if change == "geometry":
        (folder / "input.xyz").write_text("another geometry")
    elif change == "state":
        report["task"]["key"] = "other_case"
    elif change == "mutated_status":
        (folder / "status.json").write_text('{}')
    else:
        p = folder / "independent_review.json"
        data = json.loads(p.read_text())
        data["qualified_G" if change == "promoted_scope" else "maximum_relative_error"] = True if change == "promoted_scope" else .02
        p.write_text(json.dumps(data))
    with pytest.raises((ValueError, KeyError)):
        m.response_evidence(folder, report)
