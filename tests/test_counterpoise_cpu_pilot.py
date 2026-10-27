import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('counterpoise_cpu_pilot',ROOT/'scripts/pipeline/counterpoise_cpu_pilot.py')
pilot=importlib.util.module_from_spec(spec);spec.loader.exec_module(pilot)


def test_one_spacer_case_is_bound_to_independently_reviewed_geometry():
    case=(ROOT/'results/validation/local_counterpoise_2026-10-26_cpu_one_spacer/occ_2_4__pose1__cap109').resolve()
    review=(ROOT/'results/reviews/2026-10-26/counterpoise_control/prepared_one_spacer_cpu_geometry.json').resolve()
    plan,sources=pilot.source_hashes(case,review)
    assert plan['source_construction_key']=='occ_2_4__pose1__cap109'
    assert sources[str(case/'plan.json')]
    assert not plan['solvent_counterpoise_computed'] and not plan['primary_candidate_admission']


def test_launch_requires_shared_durable_storage(tmp_path):
    with pytest.raises(ValueError,match='durable'):
        pilot.launch(tmp_path,tmp_path,tmp_path/'raw',8)
