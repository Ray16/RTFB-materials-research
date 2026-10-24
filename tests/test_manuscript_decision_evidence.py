"""The paper's independent reconstruction must reject misleading promotions."""
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("manuscript_claims", ROOT / "manuscript/check_claims.py")
CLAIMS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLAIMS)
EPOCH = ROOT / "manuscript/evidence/2026-10-23_reference_origin_review"


def load(name):
    return json.loads((EPOCH / name).read_text())


def test_checked_followup_reconstructs():
    result = CLAIMS.audit_decision_followup(load)
    assert result["finite_native_G"] == 12
    assert result["capacity_compositions"] == 21
    assert result["all_axis_eligible"] == 0


def deployment_load(name):
    return json.loads((ROOT / 'manuscript/evidence/2026-10-23_deployment_revision' / name).read_text())


def test_fresh_deployment_claims_reconstruct():
    assert CLAIMS.audit_quantum_deployment(deployment_load)['controls'] == 2


@pytest.mark.parametrize('filename,mutate,error', [
    ('quantum_environment.json', lambda x: x.update(full_portability_claimed=True), 'qualification scope'),
    ('quantum_environment_states.json', lambda x: x['states'][0]['fresh'].update(G_eV=0), 'free-energy reconstruction'),
    ('quantum_environment_states.json', lambda x: x['states'][0]['fresh']['thermal']['frequencies_cm'].pop(), 'mode inventory'),
    ('quantum_environment_states.json', lambda x: x['states'].pop(), 'Both deployment'),
    ('quantum_environment_compatibility.json', lambda x: x['third_party_changes'][0].update(name='numpy'), 'qualification scope'),
])
def test_fresh_deployment_cannot_be_promoted_or_misbound(filename, mutate, error):
    def changed(name):
        value = deployment_load(name)
        if name == filename:
            mutate(value)
        return value
    with pytest.raises(ValueError, match=error):
        CLAIMS.audit_quantum_deployment(changed)


@pytest.mark.parametrize("filename,mutate,error", [
    ("pmdi_capacity.json", lambda x: x["exact_graphs_at_handle_density_0_6"][1].update(
        two_e_upper_bound_mAh_g=153.9), "Capacity upper bound"),
    ("all_axes_decisions.json", lambda x: x["candidates"][0].update(
        all_required_axes_pass=True, eligible_for_Pareto=True), "bypassed all-axis"),
    ("finite_ensemble.json", lambda x: x["cycles"][0]["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"].update(
        {"50": -10.98}), "ensemble occupancy cycle"),
    ("pmdi_precursor_review.json", lambda x: x.update(
        exact_polymer_grafting_demonstrated=True), "completed attachment route"),
])
def test_reconstruction_rejects_promoted_or_mixed_evidence(filename, mutate, error):
    def changed(name):
        value = load(name)
        if name == filename:
            mutate(value)
        return value
    with pytest.raises(ValueError, match=error):
        CLAIMS.audit_decision_followup(changed)


def manuscript_builder():
    spec = importlib.util.spec_from_file_location("paper_build", ROOT / "manuscript/build.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_accepts_line_wrapped_draft_labels():
    manuscript_builder().check_draft_labels("Working draft\nDBU\nremain sub-\nmission\nrequirements")


def test_scientific_conclusion_does_not_need_editorial_task_language():
    manuscript_builder().check_draft_labels("Working draft\nDBU\nFinite conformer sampling limits material inference.")


def test_build_still_requires_pending_labels():
    with pytest.raises(RuntimeError, match="Draft labels"):
        manuscript_builder().check_draft_labels("Working draft\nDBU\nready for publication")


def test_paper_reconstruction_rejects_computed_reference_promoted_to_experiment():
    def changed(name):
        value = load(name)
        if name == "basicity_reference_origin.json":
            value["independent_experimental_pKaH_MeCN"] = 18.6
        return value
    with pytest.raises(ValueError, match="promoted to experimental"):
        CLAIMS.audit_decision_followup(changed)
