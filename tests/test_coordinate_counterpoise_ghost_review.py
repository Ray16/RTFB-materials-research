import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ghost_review_coordinator',ROOT/'scripts/pipeline/coordinate_counterpoise_ghost_review.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_upstream_and_ghost_states_do_not_promote_live_or_failed_work():
    assert module.classify_upstream({'status':'waiting'},True)=='waiting'
    assert module.classify_upstream({'status':'waiting'},False)=='upstream_owner_stopped_without_result'
    assert module.classify_upstream({'status':'ghost_pilot_started_unqualified','ghost_started':True},False)=='ghost_started'
    assert module.classify_ghost({'status':'CPU_ghost_pilot_SCF_running_unqualified'},True)=='waiting_for_ghost'
    assert module.classify_ghost({'status':'CPU_ghost_pilot_converged_independent_review_pending'},True)=='waiting_for_owner_exit'
    assert module.classify_ghost({'status':'CPU_ghost_pilot_converged_independent_review_pending'},False)=='review_ready'
    assert module.classify_ghost({'status':'CPU_ghost_pilot_interrupted_or_failed_preserved'},False)=='ghost_failed_preserved'
    with pytest.raises(ValueError,match='Unknown'):
        module.classify_upstream({'status':'ghost_pilot_started_unqualified','ghost_started':False},False)


def test_water_evaluator_is_bound_to_exact_code():
    water=module.water_control()
    assert water['energy_error_Ha']<=1e-8
    assert water['orbital_gradient_norm_au']<=water['orbital_gradient_threshold_au']
