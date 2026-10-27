import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('cp_manifest',ROOT/'scripts/validation/collect_counterpoise_manifest.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def fixture(tmp_path,monkeypatch):
    monkeypatch.setattr(m,'ROOT',tmp_path)
    return dict(schema=1,remaining_terms={r:dict(raw=str(tmp_path/'calcs'/r),review=str(tmp_path/'calcs'/r/'review.json'),stationarity=str(tmp_path/'calcs'/r/'stationarity.json')) for r in m.ROLES})


def test_explicit_term_locations_are_distinct_and_project_bound(tmp_path,monkeypatch):
    data=fixture(tmp_path,monkeypatch);assert set(m.locations(data))==m.ROLES
    data['remaining_terms']['AB']['raw']=data['remaining_terms']['B_own']['raw']
    with pytest.raises(ValueError):m.locations(data)


def test_missing_role_and_outside_paths_are_rejected(tmp_path,monkeypatch):
    data=fixture(tmp_path,monkeypatch);data['remaining_terms']['AB']['raw']='/outside'
    with pytest.raises(ValueError):m.locations(data)
    del data['remaining_terms']['AB']
    with pytest.raises(ValueError):m.locations(data)


def test_collector_uses_the_existing_physical_gate_template():
    module=m.api.load_module('actual_counterpoise_collector',m.ORIGINAL)
    fn=m.adapted_collect(module,{},{})
    assert callable(fn)
