from copy import deepcopy
import pytest
from redox.screening.feasibility import pareto_after_feasibility


def fixture():
    axes=[dict(id=k,kind='numeric',scope=k,unit='test',minimum=-100,maximum=100,
        accepted_interval_kinds=['qualified_scenario_envelope'],criterion_source='Test criterion') for k in ('capacity','coupling')]
    policy=dict(profile_id='scope_test',decision_scope='Test model',required_axes=axes,context=dict(pools_evaluated_separately=True))
    candidate=dict(id='a',product_key='graph_a',electrode_pools=['anolyte'],axes={
        axis['id']:dict(product_key='graph_a',scope=axis['id'],status='qualified',unit='test',interval=[1,1],interval_kind='qualified_scenario_envelope',sources=[dict(path='record.json',sha256='a'*64)]) for axis in axes})
    return policy,candidate


def test_electrode_roles_are_not_mixed():
    policy,a=fixture();b=deepcopy(a);b['id']='b';b['electrode_pools']=['catholyte']
    objectives=[dict(axis='capacity',direction='maximize'),dict(axis='coupling',direction='minimize_magnitude')]
    with pytest.raises(ValueError,match='electrode pool'):pareto_after_feasibility([a,b],policy,objectives)
    result=pareto_after_feasibility([a,b],policy,objectives,pool='anolyte')
    assert result['front_ids']==['a']


def test_signed_coupling_uses_magnitude():
    policy,a=fixture();b=deepcopy(a);b['id']='b'
    a['axes']['coupling']['interval']=[-.05,-.03];b['axes']['coupling']['interval']=[-.1,-.08]
    result=pareto_after_feasibility([a,b],policy,[dict(axis='capacity',direction='maximize'),dict(axis='coupling',direction='minimize_magnitude')],pool='anolyte')
    assert result['front_ids']==['a']


def test_best_axes_cannot_be_combined_from_different_electrolytes():
    from redox.screening.feasibility import assess_candidate
    policy,a=fixture();policy['shared_context_groups']=[dict(id='operating',axes=['capacity','coupling'],required_context=dict(medium='MeCN'))]
    a['axes']['capacity']['context']=dict(medium='MeCN',salt='LiPF6')
    a['axes']['coupling']['context']=dict(medium='MeCN',salt='TBAPF6')
    with pytest.raises(ValueError,match='incompatible contexts'):assess_candidate(a,policy)


def test_missing_operating_context_cannot_pass():
    from redox.screening.feasibility import assess_candidate
    policy,a=fixture();policy['shared_context_groups']=[dict(id='operating',axes=['capacity','coupling'],required_context=dict(medium='MeCN'))]
    assert assess_candidate(a,policy)['status']=='pending'
