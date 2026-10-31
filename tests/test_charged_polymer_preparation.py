"""Preparation failures stay visible and cannot qualify quantum/material results."""
import importlib.util
import json
from pathlib import Path

import pytest

from redox.build.ionic_polymer import capacity_bound
from redox.build.polymer_cases import build_case

ROOT=Path(__file__).parents[1]
SPEC=importlib.util.spec_from_file_location('charged_preparation',ROOT/'scripts/pipeline/prepare_charged_polymer.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)
CONFIG=json.loads((ROOT/'config/viologen_polymer_preparation.json').read_text())
FRAGMENT='[*:1][n+]1ccc(-c2cc[n+](C)cc2)cc1'


def test_bounded_random_fallback_keeps_real_geometry_and_stereo_audits(tmp_path,monkeypatch):
    original=MODULE.AllChem.EmbedMultipleConfs;calls=[]
    def embed(mol,**kwargs):
        calls.append(kwargs['params'].useRandomCoords)
        return [] if len(calls)==1 else original(mol,**kwargs)
    monkeypatch.setattr(MODULE.AllChem,'EmbedMultipleConfs',embed)
    case=build_case(FRAGMENT,2,[1],[],['CW'])
    case['composition_bounds']={'2':capacity_bound(FRAGMENT,CONFIG['ionic_model'],2,1,0,2)}
    settings=dict(CONFIG['initial_geometry'],embeds_per_graph=1,embedding_max_iterations=100)
    result=MODULE.geometry(case,FRAGMENT,settings,tmp_path,20261031)
    assert calls==[False,True]
    assert result['selected']['connectivity']==result['selected']['stereochemistry']=='ok'
    assert result['XYZ_backbone_carbon_count']==4 and result['XYZ_backbone_CC_bonds']==3
    assert not result['native_G_qualified'] and not result['solvent_optimized']
    assert result['embedding_attempts'][0]['conformers']==0
    assert result['embedding_attempts'][1]['conformers']==1
    assert (tmp_path/'embedding.json').is_file()


def test_total_embedding_failure_records_both_attempts(tmp_path,monkeypatch):
    monkeypatch.setattr(MODULE.AllChem,'EmbedMultipleConfs',lambda *a,**k: [])
    case=build_case(FRAGMENT,2,[1],[],['CW'])
    with pytest.raises(ValueError,match='bounded attempts'):
        MODULE.geometry(case,FRAGMENT,CONFIG['initial_geometry'],tmp_path,20261031)
    report=json.loads((tmp_path/'embedding.json').read_text())
    assert len(report['attempts'])==2
    assert not list(tmp_path.glob('trial_*.xyz'))


def test_failed_cases_leave_exact_graphs_sources_and_terminal_partial_result(tmp_path,monkeypatch):
    def failure(*args):raise ValueError('Deliberate embedding failure')
    monkeypatch.setattr(MODULE,'geometry',failure)
    source=ROOT/'results/candidates/pool/2026-10-30_viologen_checked'
    output=tmp_path/'prepared'
    result=MODULE.prepare(ROOT/'config/viologen_polymer_preparation.json',source/'registry.json',source/'inputs/route_evidence.json',source/'inputs/assessment.json',output)
    assert result['status']=='initial_geometry_partial' and len(result['failures'])==12
    assert not result['cases'] and result['quantum_jobs_submitted']==0
    assert not result['material_admission'] and not result['quantum_calculation_admitted']
    assert len(list(output.glob('*/case.json')))==12
    assert len(list(output.glob('*/failure.json')))==12
    assert (output/'plan.json').is_file() and (output/'preparation.json').is_file()
    assert all((output/name).is_file() for name in result['snapshot_files'].values())
    with pytest.raises(FileExistsError):
        MODULE.prepare(ROOT/'config/viologen_polymer_preparation.json',source/'registry.json',source/'inputs/route_evidence.json',source/'inputs/assessment.json',output)


@pytest.mark.parametrize('mutation',['none','capacity','site_binding','cycle','qualification','inventory'])
def test_independent_actual_preparation_replay_and_negative_controls(tmp_path,mutation):
    spec=importlib.util.spec_from_file_location('check_charged_preparation',ROOT/'scripts/validation/check_charged_polymer_preparation.py')
    checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)
    source=ROOT/'results/polymer_models/viologen_preparation_2026-10-31_v2/preparation.json'
    data=json.loads(source.read_text())
    if mutation=='none':
        result=checker.check(source)
        assert result['graphs']==12 and result['balanced_cycles']==4
        assert all(r['backbone_edges']==9 and r['chain_ends']==2 for r in result['records'])
        assert not result['native_G_qualified'] and not result['material_admission']
        return
    # Keep the real immutable input files in place; alter only the result claim.
    data['files_sha256']={str((source.parent/p).resolve()):h for p,h in data['files_sha256'].items()}
    data['snapshot_files']={p:str((source.parent/f).resolve()) for p,f in data['snapshot_files'].items()}
    (tmp_path/'plan.json').write_bytes((source.parent/'plan.json').read_bytes())
    if mutation=='capacity':data['cases'][1]['case']['composition_bounds']['2']['capacity_upper_bound_mAh_g']+=1
    if mutation=='site_binding':data['cases'][4]['case']['occupied_sites']=[2,4]
    if mutation=='cycle':data['cycles'][0]['coefficients'][next(iter(data['cycles'][0]['coefficients']))]+=1
    if mutation=='qualification':data['native_G_eV']=0.
    if mutation=='inventory':data['files_sha256'][str(source.resolve())]='0'*64
    changed=tmp_path/'changed.json';changed.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker.check(changed)
