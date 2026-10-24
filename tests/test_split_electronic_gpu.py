"""GPU energy handoff keeps real backend identity and permits native reuse."""
from copy import deepcopy
import json

import pytest

from redox.core.protocol import ACTIVE_SP,geom_sha1,load_record
from redox.workflows.split_electronic import composite_result


def parents():
    geometry=dict(charge=0,mult=1,backend='gpu',optimized=True,geometry_only_result=True,
        geometry_optimization_converged=True,optimization_SCF_converged=True,opt_basis='def2-svp',opt_xc='r2scan')
    active=dict(protocol=deepcopy(ACTIVE_SP),geom_sha1='coordinates',charge=0,mult=1,
        source='computed:gpu',converged_smd=True,converged_gas=True,
        e_smd_eV=-10.,e_gas_eV=-9.9,density_fit=True,conv_tol=1e-9)
    return geometry,active


def test_GPU_parent_is_not_relabelled_and_inputs_stay_unchanged():
    geometry,active=parents();before=deepcopy((geometry,active))
    result=composite_result(geometry,active,'coordinates','/gpu/energy.json',backend='gpu')
    assert (geometry,active)==before
    assert result['electronic_backend']=='gpu' and result['geometry_backend']=='gpu'
    assert result['electronic_record_source']=='/gpu/energy.json' and not result['geometry_only_result']


@pytest.mark.parametrize('change',[dict(source='computed:cpu'),dict(geom_sha1='other'),dict(charge=False),dict(converged_gas=False),dict(e_gas_eV=float('nan'))])
def test_wrong_backend_state_geometry_or_incomplete_energy_cannot_be_joined(change):
    geometry,active=parents();active.update(change)
    with pytest.raises(ValueError):composite_result(geometry,active,'coordinates','/gpu/energy.json',backend='gpu')


def test_native_handoff_does_not_repeat_the_checked_GPU_energy(tmp_path,monkeypatch):
    from redox.qm import sp,dft
    (tmp_path/'opt.xyz').write_text('3\ncontract control\nO 0 0 0\nH .96 0 0\nH -.24 .93 0\n')
    geometry,active=parents();identity=geom_sha1(tmp_path/'opt.xyz');active['geom_sha1']=identity
    result=composite_result(geometry,active,identity,'/gpu/original.json',backend='gpu')
    (tmp_path/'result.json').write_text(json.dumps(result))
    def forbidden(*args,**kwargs):raise AssertionError('Do not recompute a checked GPU electronic parent')
    monkeypatch.setattr(dft,'dft_smd',forbidden)
    status,record=sp.ensure('control','rest',tmp_path,0,1,backend='gpu')
    assert status=='adopted' and record['e_smd_eV']==active['e_smd_eV']
    assert load_record(tmp_path,0,1)['e_gas_eV']==active['e_gas_eV']


def test_general_CPU_join_matches_the_retained_CPU_contract():
    from redox.workflows.split_polymer import composite_result as original
    geometry,active=parents();active['source']='computed:cpu'
    assert composite_result(geometry,active,'coordinates','/cpu/original.json',backend='cpu')==original(geometry,active,'coordinates','/cpu/original.json')


def test_duplicate_source_hashes_merge_only_when_the_epochs_agree():
    import importlib.util
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('GPU_handoff_hash_merge',root/'scripts/pipeline/finish_polymer_electronic_gpu.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    assert module.merge_hashes({'same':'sha'},{'same':'sha','new':'other'})=={'same':'sha','new':'other'}
    with pytest.raises(ValueError,match='Conflicting source epochs'):
        module.merge_hashes({'same':'old'},{'same':'new'})


def test_preparation_keeps_native_input_handoff_and_exact_geometry(tmp_path,monkeypatch):
    import importlib.util
    from pathlib import Path
    from types import SimpleNamespace
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('GPU_prepare_complete',root/'scripts/pipeline/finish_polymer_electronic_gpu.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    geometry,pilot,review,out=(tmp_path/p for p in ('geometry','pilot','review','new'))
    for p in (geometry,pilot,review):p.mkdir()
    task=dict(key='state',charge=0,mult=1)
    plan=dict(tasks=[task],quantum_settings={'same':'settings'},thermochemistry={'temperature_K':298.15},
        solvent='acetonitrile',producer_sha256={},input_sha256={})
    for p in (geometry,pilot):(p/'plan.json').write_text(json.dumps(dict(plan,maximum_memory_pilot_MiB=31744,execution_workspace_GiB=28)))
    (geometry/'controls').mkdir();live=geometry/'live_inputs/state';live.mkdir(parents=True)
    xyz=live/'input.xyz';xyz.write_text('1\nsource\nHe 0 0 0\n')
    (live/'task.json').write_text(json.dumps(dict(input_xyz=str(xyz),input_sha256={str(xyz):module.api.digest(xyz)})))
    (geometry/'geometry').mkdir();(geometry/'geometry/status.json').write_text(json.dumps(dict(
        status='stationary_geometry_checked_active_parent_and_native_minimum_pending',owner={},output_sha256={},state_dir='source-state')))
    (review/'summary.json').write_text(json.dumps(dict(status='electronic_resource_point_independently_reviewed',
        source=str(pilot),source_sha256={},host_RSS_bound_GiB=283)))
    monkeypatch.setattr(module.api,'load_module',lambda name,path:SimpleNamespace(owner_alive=lambda owner:False,checked_controls=lambda controls:{}))
    module.prepare(geometry,pilot,review,tmp_path/'controls',out)
    new=json.loads((out/'plan.json').read_text());ready=json.loads((out/'live_inputs/state/task.json').read_text())
    assert (out/'live_inputs/state/input.xyz').read_bytes()==xyz.read_bytes()
    assert ready['input_xyz']==str(out/'live_inputs/state/input.xyz')
    assert new['tasks']==plan['tasks'] and new['quantum_settings']==plan['quantum_settings']
    module.api.verify_sources(new['input_sha256']);module.api.verify_sources(new['producer_sha256'])
