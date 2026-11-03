from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from redox.screening.qualified_cases import attach,read_case,validate_case

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'results/candidates/verified_cases/phenothiazine/2026-11-03_primary_redox_context/candidate.json'


def test_actual_case_replays_and_preserves_unresolved_axes():
    case=read_case(CASE)
    assert case['properties']['first_oxidation']['value']==pytest.approx(.6260868628669414)
    assert case['all_axis_decision']['status']=='pending'
    assert case['properties']['grafting_free_energy'] is None
    assert case['properties']['electronic_coupling'] is None
    assert not case['material_admission']


@pytest.mark.parametrize('change',['identity','value','solvent','thermal','reference','policy','missing','invented','material','spin','model_scope'])
def test_changed_case_is_not_a_qualified_property(change):
    case=read_case(CASE);prop=case['properties']['first_oxidation']
    if change=='identity':case['product_smiles']='CCO'
    if change=='value':prop['value']+=.1
    if change=='solvent':prop['context']['solvent']='water'
    if change=='thermal':prop['context']['thermal_definition']['selected_rotor_cutoff_cm']=100
    if change=='reference':prop['context']['reference']['electrode']='SHE'
    if change=='policy':prop['component_numerical_policies']['graft_neu']['Hessian_gradient_step_Bohr']=1e-4
    if change=='missing':case['properties'].pop('electronic_coupling')
    if change=='invented':case['properties']['grafting_free_energy']=-100.
    if change=='material':case['material_admission']=True
    if change=='spin':prop['component_states']['minus']['graft_ox']['multiplicity']=1
    if change=='model_scope':prop['context']['model_scope']='full polymer'
    with pytest.raises(ValueError):validate_case(case)


def test_registry_join_retains_prior_evaluations_and_current_route_hold():
    case=read_case(CASE)
    registry=json.loads((CASE.parent/'inputs/registry.json').read_text())
    target=next(r for r in registry['candidates'] if r['graph_key']==case['registry_graph_key'])
    target['route_status']='hold';target['evaluations']=[dict(source='historical',reported_values={'other':1})]
    before=deepcopy(registry['candidates'])
    attach(registry['candidates'],case,'qualified','a'*64)
    assert target['route_status']=='hold' and target['material_decision']=='pending'
    assert target['evaluations'][0]==before[registry['candidates'].index(target)]['evaluations'][0]
    assert target['evaluations'][-1]['all_axis_status']=='pending'
    for old,new in zip(before,registry['candidates']):
        if old['graph_key']!=case['registry_graph_key']:assert old==new
    with pytest.raises(ValueError,match='already attached'):attach(registry['candidates'],case,'qualified','a'*64)


def test_frozen_case_inputs_cannot_be_changed(tmp_path):
    dest=tmp_path/'case';shutil.copytree(CASE.parent,dest)
    p=dest/'inputs/redox.json';d=json.loads(p.read_text());d['redox_potential_qualified']=False;p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='snapshot changed'):read_case(dest/'candidate.json')
