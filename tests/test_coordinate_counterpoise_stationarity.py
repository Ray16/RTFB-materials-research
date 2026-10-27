import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('stationarity_coordinator',ROOT/'scripts/pipeline/coordinate_counterpoise_stationarity.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_only_successful_checkpoint_review_releases_stationarity():
    assert module.classify({'status':'waiting_for_pilot'},True)=='waiting'
    assert module.classify({'status':'waiting_for_pilot'},False)=='review_watcher_stopped_without_result'
    assert module.classify({'status':'pilot_checkpoint_review_complete_gradient_pending','review_complete':True},False)=='ready'
    assert module.classify({'status':'review_failed_preserved'},False)=='prior_stage_failed_preserved'
    with pytest.raises(ValueError,match='Unknown'):
        module.classify({'status':'pilot_checkpoint_review_complete_gradient_pending','review_complete':False},True)


def test_recomputed_water_control_is_source_bound():
    result=module.water_control()
    assert result['passed'] and result['orbital_gradient_norm_au']<=result['orbital_gradient_threshold_au']
    assert result['reference_energy_error_Ha']<=1e-8
