"""Geometry refinement keeps exact state/method and cannot inherit T4 admission."""
import json
import shlex
from types import SimpleNamespace
import pytest

from redox.core.common import ROOT
from redox.workflows.route_geometry import validate_plan
from redox.workflows.molecular_native_namespace import load


@pytest.fixture
def plan():
    return json.loads((ROOT/'validation/phenothiazine_geometry_2026-11-02/minus_plan.json').read_text())


def test_signed_geometry_contract_keeps_full_workflow_resource_and_physics_gates(plan):
    assert validate_plan(plan)['minimum_GPU_memory_MiB']==31744
    assert plan['charge']==0 and plan['mult']==1 and not plan['native_G_qualified']


@pytest.mark.parametrize('field,value',[('charge',1),('mult',3),('nao',450),('atoms',42),
    ('minimum_GPU_memory_MiB',15360),('native_G_qualified',True),('geometry_memory_pilot',True),
    ('optimizer_convergence',{}),('method',{})])
def test_changed_geometry_state_method_or_resources_need_separate_admission(plan,field,value):
    plan[field]=value
    with pytest.raises(ValueError):validate_plan(plan)


@pytest.mark.parametrize('launcher',[None,'redox.workflows.route_geometry'])
def test_shared_remote_launcher_preserves_legacy_default_and_explicit_scientific_adapter(tmp_path,monkeypatch,launcher):
    module=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','tested_scientific_launcher_adapter')
    row=dict(namespace=str(tmp_path/'namespace.json'))
    if launcher is not None:row['launcher_module']=launcher
    calls=[]
    def run(argv,**kwargs):
        calls.append(argv);return SimpleNamespace(stdout='{"status":"recorded"}')
    monkeypatch.setattr(module.subprocess,'run',run)
    assert module.launch(row,'host',0,'/python',tmp_path/'claim',1)['status']=='recorded'
    remote=shlex.split(calls[0][-1]);assert remote[-1]==(launcher or 'redox.workflows.screening_control')


def test_unmaintained_launcher_cannot_be_executed(tmp_path):
    module=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','tested_scientific_launcher_adapter')
    with pytest.raises(ValueError):module.launch(dict(namespace=str(tmp_path/'namespace.json'),launcher_module='other.command'),'host',0,'/python',tmp_path/'claim',1)
