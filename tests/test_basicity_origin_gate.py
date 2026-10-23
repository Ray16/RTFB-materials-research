import importlib.util
import json
from pathlib import Path
import pytest

SPEC=importlib.util.spec_from_file_location('basicity_origin',Path(__file__).resolve().parents[1]/'scripts/pipeline/review_basicity_origins.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def test_computed_value_can_support_only_a_model_comparison():
    r=dict(origin_class='computed_PCM_estimate',origin_conflict=True,independent_experimental_anchor=False)
    assert not MODULE.origin_role(r,'model_sensitivity')['independent_experimental_origin_eligible']
    with pytest.raises(ValueError,match='measured'):MODULE.origin_role(r,'experimental_validation')


def test_relabeling_computed_value_as_measured_is_rejected():
    r=dict(origin_class='computed_PCM_estimate',independent_experimental_anchor=True)
    with pytest.raises(ValueError,match='promoted'):MODULE.origin_role(r,'model_sensitivity')


def test_measured_origin_alone_does_not_validate_full_route():
    r=dict(origin_class='experimental_measurement',origin_conflict=False,independent_experimental_anchor=True)
    output=MODULE.origin_role(r,'experimental_validation')
    assert output['independent_experimental_origin_eligible'] and not output['full_route_accuracy_validated']


@pytest.mark.parametrize('field,value',[('solvent','DMSO'),('reported_value',20.0)])
def test_source_review_does_not_admit_a_different_value_or_medium(field,value):
    reference=json.loads((Path(__file__).resolve().parents[1]/'config/basicity_reference_origins.json').read_text())['references'][0]
    MODULE.reviewed_DIPEA_context(reference);reference[field]=value
    with pytest.raises(ValueError,match='value/context'):MODULE.reviewed_DIPEA_context(reference)
