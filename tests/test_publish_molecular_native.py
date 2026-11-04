import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('tested_exact_molecular_redox',ROOT/'scripts/pipeline/publish_molecular_native.py')
publisher=importlib.util.module_from_spec(definition);definition.loader.exec_module(publisher)


@pytest.fixture
def values():
    energy=lambda x:{str(c):x for c in (25,50,100)}
    states=[dict(state=s,charge=q,mult=m,product_inchikey='same',G_by_cutoff_eV=energy(g))
        for s,q,m,g in [('ox2',2,1,-10.),('ox1',1,2,-13.4),('neu',0,1,-16.1)]]
    fc=dict(neu=dict(G_by_cutoff_eV=energy(-10.)),ox=dict(G_by_cutoff_eV=energy(-6.)))
    return states,fc


def test_adjacent_native_redox_steps_use_the_same_unfitted_Fc_reference(values):
    result=publisher.ladder(*values)
    assert result['events'][0]['potential_V_vs_Fc_by_cutoff']['50']==pytest.approx(-.6)
    assert result['events'][1]['potential_V_vs_Fc_by_cutoff']['50']==pytest.approx(-1.3)
    assert result['disproportionation_kJ_mol_by_cutoff']['50']==pytest.approx(.7*96.48533212331)
    assert result['fitted_offset_V']==0.
    assert not result['experimental_exact_graft_accuracy_validated']


@pytest.mark.parametrize('change',['nonfinite','identity','duplicate_charge','charge_gap','missing_cutoff'])
def test_missing_or_mismatched_states_cannot_enter_the_redox_ladder(values,change):
    states,fc=values
    if change=='nonfinite':states[0]['G_by_cutoff_eV']['50']=float('nan')
    if change=='identity':states[1]['product_inchikey']='another'
    if change=='duplicate_charge':states[1]['charge']=2
    if change=='charge_gap':states[0]['charge']=3
    if change=='missing_cutoff':fc['ox']['G_by_cutoff_eV'].pop('25')
    with pytest.raises(ValueError):publisher.ladder(states,fc)
