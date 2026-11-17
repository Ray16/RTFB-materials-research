from copy import deepcopy
import pytest
from redox.properties.context import combine_properties


def records():
    ctx=dict(quantity='G',unit='eV',electronic_protocol={'xc':'example','basis':'example'},
        thermal_definition='example_qrrho',solvent='MeCN',solvent_model='SMD',temperature_K=298.15,
        standard_state='1M',reference='absolute_component_free_energy')
    return {k:dict(value=v,context=deepcopy(ctx),state=dict(identity=k,charge=0,multiplicity=1),
                   qc='independently_checked',sources={'source':'test-only'}) for k,v in [('a',2.),('b',1.)]}


def test_combination_retains_context_states_and_sources():
    result=combine_properties(records(),{'a':1,'b':-1})
    assert result['value']==1 and result['context']['unit']=='eV'
    assert result['components']['a']['state']['multiplicity']==1
    assert result['components']['b']['sources']=={'source':'test-only'}
    assert result['empirical_accuracy_validated'] is False


@pytest.mark.parametrize('field',['solvent','solvent_model','reference','thermal_definition','electronic_protocol','standard_state','unit'])
def test_incompatible_physical_context_cannot_be_combined(field):
    x=records();x['b']['context'][field]='different'
    with pytest.raises(ValueError,match='contexts differ'):combine_properties(x,{'a':1,'b':-1})


def test_missing_spin_qc_or_nonfinite_value_is_rejected():
    for field in ['spin','qc','energy']:
        x=records()
        if field=='spin':x['a']['state'].pop('multiplicity')
        if field=='qc':x['a']['qc']='pending'
        if field=='energy':x['a']['value']=float('nan')
        with pytest.raises(ValueError):combine_properties(x,{'a':1,'b':-1})
