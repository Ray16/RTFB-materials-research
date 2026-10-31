import hashlib
import json
from pathlib import Path

import pytest

from redox.workflows.charged_polymer_pilot import resting_state
from redox.workflows.matched_polymer_reference import select_seed


def test_neutral_reference_requires_explicit_scope_and_zero_grafts():
    case=dict(smiles='CC',charge=0,occupied_sites=[],proposed_global_states=[dict(charge=0,mult=1)])
    with pytest.raises(ValueError):resting_state(case)
    state,empty=resting_state(case,True);assert empty and state['charge']==0
    case['occupied_sites']=[2]
    with pytest.raises(ValueError):resting_state(case,True)


@pytest.fixture
def seeds(tmp_path):
    model=dict(name='same model');case=dict(configuration='all_same',charge=0,occupied_sites=[],graph_sha256='exact',smiles='C[C@H](F)Cl',repeat_units=5,reactive_sites=[2,3,4])
    rows={}
    for i,value in enumerate([-1.,-2.]):
        xyz=tmp_path/f'{i}.xyz';xyz.write_text('source geometry '+str(i));review=tmp_path/f'{i}.json'
        task=dict(case,mult=1);review.write_text(json.dumps(dict(status='alternate_native_state_independently_reconstructed',task=task,experimental_model=model,G_by_cutoff_eV={'50':value})))
        rows[f'all_same__occ_none__conf{i+1:03d}']=dict(review=str(review),review_sha256=hashlib.sha256(review.read_bytes()).hexdigest(),geometry=str(xyz),geometry_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest())
    analysis=tmp_path/'analysis.json';analysis.write_text(json.dumps(dict(inputs=rows)));return analysis,case,model,tmp_path


def test_seed_selection_is_exact_bound_and_never_imports_old_thermochemistry(seeds):
    analysis,case,model,root=seeds;selected=select_seed(analysis,case,model)
    assert selected['key'].endswith('002') and selected['candidates_checked']==2
    assert selected['old_G_used_for_seed_selection_only'] and not selected['old_E_or_G_admitted_to_new_cycle']


@pytest.mark.parametrize('change',['bytes','graph','model','occupied'])
def test_incompatible_or_changed_reference_is_rejected(seeds,change):
    analysis,case,model,root=seeds
    if change=='bytes':(root/'0.xyz').write_text('changed')
    if change=='graph':case['graph_sha256']='different'
    if change=='model':model={'name':'different model'}
    if change=='occupied':case['occupied_sites']=[2]
    with pytest.raises(ValueError):select_seed(analysis,case,model)
