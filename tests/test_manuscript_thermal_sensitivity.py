"""Frozen thermal-setting scenarios cannot become prediction intervals or tuned primary values."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cutoff_claims", ROOT / "manuscript/check_claims.py")
CLAIMS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLAIMS)
EPOCH = ROOT / "manuscript/evidence/2026-10-23_thermal_cutoff_revision"


def load(name):
    return json.loads((EPOCH / name).read_text())


def test_all_cutoff_endpoint_and_reaction_aggregates_reconstruct():
    result = CLAIMS.audit_thermal_settings(load)
    assert result["benchmark_endpoints"] == 39 and result["reaction_models"] == 87
    assert result["maximum_endpoint_span_mV"] == pytest.approx(25.4275196959)
    assert result["maximum_reaction_span_kJ_mol"] == pytest.approx(4.45953084353)
    assert result["reaction_sign_crossings"] == 0 and not result["prediction_intervals_established"]


@pytest.mark.parametrize("change,error", [
    (lambda x: x["benchmark"]["metrics_by_cutoff"]["25"].update(MAE_V=0), "statistics"),
    (lambda x: x["benchmark"]["endpoints"].pop(), "coverage"),
    (lambda x: x["benchmark"]["endpoints"][0]["E_by_cutoff_V"].update({"50": 0}), "Primary benchmark"),
    (lambda x: x["grafting"]["reactions"][0].update(range_kJ_mol=100), "Reaction cutoff span"),
    (lambda x: x.update(candidate_prediction_intervals_established=True), "promoted beyond"),
])
def test_relabelled_or_changed_cutoff_results_are_rejected(change, error):
    def altered(name):
        value = load(name)
        if name == "thermal_cutoff_sensitivity.json":
            change(value)
        return value
    with pytest.raises(ValueError, match=error):
        CLAIMS.audit_thermal_settings(altered)
