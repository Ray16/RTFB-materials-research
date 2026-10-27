import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('remaining_coordinator',ROOT/'scripts/pipeline/coordinate_counterpoise_remaining.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_upstream_requires_successful_ghost_stationarity():
    assert module.classify_upstream({'status':'waiting'},True)=='waiting'
    assert module.classify_upstream({'status':'waiting'},False)=='upstream_owner_stopped_without_result'
    assert module.classify_upstream({'status':'ghost_stationarity_failed_preserved'},False)=='upstream_failed_preserved'
    assert module.classify_upstream({'status':'ghost_stationarity_complete_no_full_interaction','ghost_stationarity_complete':True},False)=='ready'
    with pytest.raises(ValueError,match='Unknown'):
        module.classify_upstream({'status':'ghost_stationarity_complete_no_full_interaction','ghost_stationarity_complete':False},True)


def test_memory_reserves_each_live_workers_unconsumed_ceiling(monkeypatch):
    def process(pid):return SimpleNamespace(memory_info=lambda:SimpleNamespace(rss={1:0,2:8*module.GIB}[pid]))
    monkeypatch.setattr(module.psutil,'Process',process)
    states={'a':{'status':'starting','pid':1},'b':{'status':'running','owner':{'pid':2}},'c':{'status':'completed_review_pending'}}
    assert module.future_reservation_gib(states)==56


def test_failed_or_missing_term_never_counts_as_complete(tmp_path):
    assert module.classify_term(tmp_path/'absent')['status']=='pending'
    raw=tmp_path/'failed';raw.mkdir();(raw/'status.json').write_text(json.dumps({'status':'CPU_remaining_term_interrupted_or_failed_preserved'}))
    assert module.classify_term(raw)['status']=='failed_preserved'


def test_water_controls_are_source_bound():
    control=module.water_control()
    assert control['passed'] and {x['role'] for x in control['terms']}==set(module.ROLES)
