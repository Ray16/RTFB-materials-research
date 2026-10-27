import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('checked_cp_collection',ROOT/'scripts/pipeline/collect_counterpoise_when_checked.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def test_live_or_stopped_observer_does_not_promote_an_unqualified_result():
    assert m.classify({'status':'waiting_for_refined_AB_checkpoint'},True)=='waiting_for_independent_AB_check'
    assert m.classify({'status':'waiting_for_refined_AB_checkpoint'},False)=='observer_stopped_without_qualified_result'
    assert m.classify({'status':'refined_AB_independent_stationarity_failed_preserved'},True)=='upstream_failed_preserved'


def test_independently_passed_refinement_can_enter_full_collector():
    assert m.classify({'status':'refined_AB_independent_stationarity_passed'},False)=='ready'
