"""Checkpoint durability and protection against mixing oligomer identities."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    loaded=importlib.util.module_from_spec(spec);spec.loader.exec_module(loaded)
    return loaded


def test_checkpoint_retains_each_evaluation_without_claiming_convergence(tmp_path):
    recovery=module('recovery','scripts/pipeline/recover_pi_stacking.py')
    mol=SimpleNamespace(natm=2,atom_coords=lambda unit:np.array([[0.,0.,0.],[0.,0.,.75]]),
                        atom_pure_symbol=lambda i:'H')
    callback=recovery.checkpoint_callback(tmp_path)
    for energy in (-1.,-1.01):
        callback(dict(mol=mol,energy=energy,gradients=np.zeros((2,3)),
                      g_scanner=SimpleNamespace(converged=True)))
    assert len((tmp_path/'evaluated_trajectory.xyz').read_text().splitlines())==8
    records=[json.loads(line) for line in (tmp_path/'evaluations.jsonl').read_text().splitlines()]
    assert [r['energy_Ha'] for r in records]==[-1.,-1.01]
    assert not json.loads((tmp_path/'latest.json').read_text())['solution_geometry_qualified']
    assert 'Unqualified SCF evaluation 2' in (tmp_path/'last_evaluated.xyz').read_text()


@pytest.mark.parametrize('change',[{'stereoisomer':1},{'solvent':'water'},{'overwrite':True}])
def test_recovery_cannot_replace_a_different_or_already_verified_slot(tmp_path,monkeypatch,change):
    collect=module('collector','scripts/pipeline/collect_pi_stacking.py')
    monkeypatch.setattr(collect,'ROOT',tmp_path);parent=tmp_path/'parent';parent.mkdir()
    monkeypatch.setattr(collect,'BASE',parent);monkeypatch.setattr(collect,'OUT',tmp_path)
    campaign=tmp_path/'recovery';(campaign/'done').mkdir(parents=True);(campaign/'claims').mkdir()
    task=dict(key='case',candidate_id='ndi',model_key='n3',model_type='covalent_oligomer',
              repeats=3,stereoisomer=0,target='open',smiles='CC',fragment='[*:1]N',atoms=8,charge=0,mult=1)
    plan=dict(parent_campaign='parent',electronic_protocol='same',solvent=change.get('solvent','acetonitrile'),raw_root='raw')
    plan_path=campaign/'plan.json';plan_path.write_text(json.dumps(plan))
    (tmp_path/'recovery_campaigns.json').write_text(json.dumps(dict(parent_campaign='parent',campaigns=[dict(path='recovery',plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest())])))
    (campaign/'done/case.json').write_text('{}')
    record=dict(task=dict(task,stereoisomer=change.get('stereoisomer',0)),E_smd_eV=-10)
    monkeypatch.setattr(collect,'verify_state_checkpoint',lambda path:record)
    parent_plan=dict(tasks=[task],electronic_protocol='same',solvent='acetonitrile')
    if 'solvent' in change:
        with pytest.raises(ValueError,match='protocol'):collect.recovery_records(parent_plan,{})
    else:
        records,invalid,_,_=collect.recovery_records(parent_plan,{'case':record} if change.get('overwrite') else {})
        assert not records and len(invalid)==1
