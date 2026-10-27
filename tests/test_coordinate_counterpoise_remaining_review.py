import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('remaining_review_coordinator',ROOT/'scripts/pipeline/coordinate_counterpoise_remaining_review.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_review_waits_for_all_terms_and_owner_exit():
    assert module.classify_upstream({'status':'remaining_terms_in_progress_unqualified'},True)=='waiting'
    assert module.classify_upstream({'status':'remaining_terms_in_progress_unqualified'},False)=='upstream_owner_stopped_without_result'
    assert module.classify_upstream({'status':'remaining_terms_completed_checkpoint_review_pending'},True)=='waiting_for_owner_exit'
    assert module.classify_upstream({'status':'remaining_terms_completed_checkpoint_review_pending'},False)=='ready'
    assert module.classify_upstream({'status':'remaining_term_failure_or_ownership_held'},False)=='upstream_failed_preserved'
    with pytest.raises(ValueError,match='Unknown'):
        module.classify_upstream({'status':'invented'},False)


def test_stationarity_workers_reserve_unconsumed_memory(monkeypatch):
    monkeypatch.setattr(module.psutil,'Process',lambda pid:SimpleNamespace(memory_info=lambda:SimpleNamespace(rss={1:0,2:8*module.GIB}[pid])))
    assert module.reserved_future_gib({'AB':SimpleNamespace(pid=1),'B_own':SimpleNamespace(pid=2)})==56


def test_water_control_covers_exact_evaluator_and_all_roles():
    control=module.water_control()
    assert control['passed'] and {r['role'] for r in control['terms']}==set(module.ROLES)
