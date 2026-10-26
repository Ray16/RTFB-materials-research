import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('energy_capacity_test',ROOT/'scripts/pipeline/dispatch_smooth_energy_capacity.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_energy_preflight_uses_declared_bound_plus_headroom(tmp_path,monkeypatch):
    helper=module.api.load_module('energy_capacity_preflight_template_test',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    (tmp_path/'batch.json').write_text(json.dumps(dict(maximum_RSS_GiB=256)))
    calls=[]
    def run(argv,**kwargs):
        calls.append(argv)
        return SimpleNamespace(returncode=0,stdout=json.dumps(dict(ok=True,python='profile-python')))
    monkeypatch.setattr(helper.subprocess,'run',run)
    probe=module.capacity_preflight(helper)
    assert probe(tmp_path,'lambda-test',{'python':'profile-python','library_path':'profile-libs'})
    assert calls[0][calls[0].index('--min-available-ram-mib')+1]=='294912'
    assert 'preflight_cpu_dependencies.py' in calls[0][1]


def test_failed_capacity_probe_is_rechecked_without_waiting_ten_minutes(tmp_path,monkeypatch):
    helper=module.api.load_module('energy_capacity_negative_cache_test',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    (tmp_path/'batch.json').write_text(json.dumps(dict(maximum_RSS_GiB=256)))
    (tmp_path/'lambda-test_preflight.json').write_text(json.dumps(dict(ok=False,python='profile-python')))
    calls=[]
    def run(argv,**kwargs):
        calls.append(argv);return SimpleNamespace(returncode=0,stdout=json.dumps(dict(ok=True,python='profile-python')))
    monkeypatch.setattr(helper.subprocess,'run',run)
    assert module.capacity_preflight(helper)(tmp_path,'lambda-test',{'python':'profile-python','library_path':'profile-libs'})
    assert len(calls)==1


def test_capacity_placement_template_compiles_before_any_GPU_action():
    assert callable(module.capacity_placement({'library_path':'profile-libs'}).launch_local)


def test_native_preflight_requires_its_own_bound_and_headroom(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('native_capacity_preflight_test',ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py')
    coordinator=importlib.util.module_from_spec(spec);spec.loader.exec_module(coordinator)
    helper=module.api.load_module('native_capacity_preflight_template_test',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    (tmp_path/'batch.json').write_text(json.dumps(dict(maximum_RSS_GiB=80)))
    calls=[]
    def run(argv,**kwargs):
        calls.append(argv);return SimpleNamespace(returncode=0,stdout=json.dumps(dict(ok=True,python='profile-python')))
    monkeypatch.setattr(helper.subprocess,'run',run)
    assert coordinator.cpu_preflight(helper,tmp_path,'lambda-test',{'python':'profile-python','library_path':'profile-libs'})
    assert calls[0][calls[0].index('--min-available-ram-mib')+1]=='114688'
