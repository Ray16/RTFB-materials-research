import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('counterpoise_point_coordinator',ROOT/'scripts/pipeline/coordinate_counterpoise_point.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_memory_probe_requires_one_finite_numeric_capacity():
    assert module.parse_memory('32768\n')==32768
    for value in ('','15360\n32768','nan','-1','unavailable'):
        with pytest.raises(ValueError):module.parse_memory(value)


def test_only_explicit_no_start_is_safe_to_retry():
    for status in ('atomic_placement_unavailable','resource_capacity_unavailable'):
        assert module.classify_launch_response({'status':status})=='not_started'
    assert module.classify_launch_response({'status':'gate_started_scientific_start_pending','owner':{'pid':5}})=='started'
    for outcome in ({'status':'gate_started_scientific_start_pending'}, {'status':'unexpected'}, {}, None):
        assert module.classify_launch_response(outcome)=='preserve'


def test_prepared_geometry_must_be_reviewed_for_the_exact_case(tmp_path):
    case=tmp_path/'case'
    checked={'status':'prepared_motif_counterpoise_geometry_independently_checked_no_scientific_start',
             'campaign':str(case),'primary_candidate_admission':False,'solvent_counterpoise_computed':False}
    module.validate_geometry_review(checked,case)
    with pytest.raises(ValueError):module.validate_geometry_review({**checked,'campaign':str(tmp_path/'other')},case)
    with pytest.raises(ValueError):module.validate_geometry_review({**checked,'status':'running'},case)
