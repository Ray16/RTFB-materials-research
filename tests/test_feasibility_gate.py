"""Admission is an intersection of complete, identity-bound axis constraints."""
from copy import deepcopy
import pytest
from redox.screening.feasibility import assess_candidate,feasible_ids

SOURCES=[dict(path='checked_record.json',sha256='a'*64)]


@pytest.fixture
def policy():
    return dict(profile_id='test_v1',decision_scope='conditional model',required_axes=[
        dict(id='capacity',kind='numeric',scope='polymer',unit='mAh/g',minimum=100,maximum=None,
             accepted_interval_kinds=['validated_prediction_interval'],criterion_source='Test acceptance specification'),
        dict(id='interaction',kind='numeric',scope='motif_pair',unit='kJ/mol',minimum=-10,maximum=10,
             accepted_interval_kinds=['qualified_scenario_envelope'],criterion_source='Test acceptance specification'),
        dict(id='route',kind='categorical',scope='exact_route',accepted_values=['proceed'],rejected_values=['exclude'],criterion_source='Test exact-route review')])


@pytest.fixture
def candidate():
    return dict(id='a',product_key='graph_a',axes={
        'capacity':dict(product_key='graph_a',scope='polymer',status='qualified',unit='mAh/g',interval=[140,160],value=150,interval_kind='validated_prediction_interval',sources=SOURCES),
        'interaction':dict(product_key='graph_a',scope='motif_pair',status='qualified',unit='kJ/mol',interval=[-5,5],interval_kind='qualified_scenario_envelope',sources=SOURCES),
        'route':dict(product_key='graph_a',scope='exact_route',status='qualified',value='proceed',sources=SOURCES)})


def test_all_axes_pass_before_front(candidate,policy):
    result=assess_candidate(candidate,policy)
    assert result['status']=='eligible'
    assert feasible_ids([result],policy=policy)=={'a'}


def test_great_capacity_cannot_rescue_strong_interaction(candidate,policy):
    candidate['axes']['capacity']['interval']=[1000,1200];candidate['axes']['capacity']['value']=1100
    candidate['axes']['interaction']['interval']=[-40,-30]
    result=assess_candidate(candidate,policy)
    assert result['status']=='excluded' and result['failed_axes']==['interaction']
    assert feasible_ids([result],policy=policy)==set()


def test_missing_axis_is_pending_not_good_or_bad(candidate,policy):
    del candidate['axes']['interaction']
    result=assess_candidate(candidate,policy)
    assert result['status']=='pending' and result['pending_axes']==['interaction'] and not result['failed_axes']


def test_boundary_overlap_pending(candidate,policy):
    candidate['axes']['interaction']['interval']=[-12,-8]
    assert assess_candidate(candidate,policy)['status']=='pending'


def test_unvalidated_point_is_not_zero_uncertainty(candidate,policy):
    candidate['axes']['interaction'].pop('interval')
    candidate['axes']['interaction']['value']=0
    assert assess_candidate(candidate,policy)['status']=='pending'


def test_sigma_prior_cannot_impersonate_validated_interval(candidate,policy):
    candidate['axes']['capacity']['interval_kind']='uncalibrated_sigma_prior'
    assert assess_candidate(candidate,policy)['status']=='pending'


def test_unset_target_is_explicitly_pending(candidate,policy):
    policy['required_axes'][0]['minimum']=None
    assert assess_candidate(candidate,policy)['status']=='pending'


@pytest.mark.parametrize('change',[('product_key','another_graph'),('scope','whole_segment'),('unit','eV')])
def test_incompatible_identity_scope_or_units_refused(candidate,policy,change):
    candidate['axes']['interaction'][change[0]]=change[1]
    with pytest.raises(ValueError):assess_candidate(candidate,policy)


def test_route_hold_cannot_be_overridden_by_thermodynamics(candidate,policy):
    candidate['axes']['route']['value']='hold'
    assert assess_candidate(candidate,policy)['status']=='pending'


def test_control_never_becomes_candidate(candidate,policy):
    candidate['role']='comparison_control'
    r=assess_candidate(candidate,policy)
    assert r['all_required_axes_pass'] and not r['eligible_for_Pareto']


def test_changed_cutoff_invalidates_cached_decision(candidate,policy):
    decision=assess_candidate(candidate,policy);changed=deepcopy(policy);changed['required_axes'][0]['minimum']=200
    with pytest.raises(ValueError,match='profile'):feasible_ids([decision],policy=changed)


def test_dropping_axis_cannot_manufacture_admission(candidate,policy):
    decision=assess_candidate(candidate,policy);decision['axes'].pop()
    with pytest.raises(ValueError,match='profile'):feasible_ids([decision],policy=policy)


def test_decision_front_uses_qualified_scoped_axes(candidate,policy):
    from redox.screening.pareto import decision_front
    other=deepcopy(candidate);other['id']='b';other['product_key']='graph_b'
    for record in other['axes'].values():record['product_key']='graph_b'
    other['axes']['capacity']['interval']=[105,110];other['axes']['capacity']['value']=108
    candidate['axes']['interaction']['interval']=[0,1]
    other['axes']['interaction']['interval']=[3,4]
    result=decision_front([candidate,other],policy,[dict(axis='capacity',direction='maximize'),dict(axis='interaction',direction='minimize')])
    assert result['eligible_ids']==['a','b'] and result['front_ids']==['a']
    # A hero descriptor cannot survive a bad required route.
    candidate['axes']['route']['value']='exclude'
    result=decision_front([candidate,other],policy,[dict(axis='capacity',direction='maximize'),dict(axis='interaction',direction='minimize')])
    assert result['front_ids']==['b']


def test_overlapping_intervals_do_not_manufacture_dominance(candidate,policy):
    from redox.screening.feasibility import pareto_after_feasibility
    other=deepcopy(candidate);other['id']='b'
    result=pareto_after_feasibility([candidate,other],policy,[dict(axis='capacity',direction='maximize'),dict(axis='interaction',direction='minimize')])
    assert result['front_ids']==['a','b']
