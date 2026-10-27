import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ghost_coordinator',ROOT/'scripts/pipeline/coordinate_counterpoise_cpu_ghost.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_only_passed_own_stationarity_releases_ghost():
    assert module.classify({'status':'waiting'},True)=='waiting'
    assert module.classify({'status':'waiting'},False)=='upstream_owner_stopped_without_result'
    assert module.classify({'status':'stationarity_failed_preserved'},False)=='upstream_failed_preserved'
    assert module.classify({'status':'pilot_stationarity_complete_no_full_interaction','stationarity_complete':True},False)=='ready'
    with pytest.raises(ValueError,match='Unknown'):
        module.classify({'status':'pilot_stationarity_complete_no_full_interaction','stationarity_complete':False},True)


def test_water_ghost_control_is_source_bound():
    control=module.water_control()
    assert control['passed'] and control['role']=='A_ghost'
    assert control['orbital_gradient_norm_au']<=control['orbital_gradient_threshold_au']
