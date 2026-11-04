from copy import deepcopy
import math

import pytest

from redox.core.common import HARTREE_EV,EV_KJ
from redox.workflows.ionic_native_cycles import native_cycle,validate_component


@pytest.fixture
def cycle():
    definition=dict(coefficients=dict(both=1,empty=1,first=-1,second=-1))
    cases={k:dict(smiles=s,charge=0) for k,s in [('both','[H][H]'),('empty','CC'),('first','C'),('second','C')]}
    context=dict(quantity='native_component_G',unit='eV',electronic_protocol={'basis':'test'},
        thermal_definition={'protocol':'test'},solvent='acetonitrile',solvent_model='test',
        temperature_K=298.15,standard_state='test convention',reference='same absolute-component model')
    rows={}
    for index,(key,energy) in enumerate([('both',-5.),('empty',-1.),('first',-2.),('second',-3.)]):
        nuclear={c:(index+1)*.1+int(c)*.0001 for c in ('25','50','100')}
        rows[key]=dict(charge=0,mult=1,native_G_qualified=True,primary_model_admission=False,
            material_admission=False,conformer_ensemble_qualified=False,E_geometry_level_Ha=energy,
            nuclear_by_cutoff_eV=nuclear,G_by_cutoff_eV={c:energy*HARTREE_EV+v for c,v in nuclear.items()},
            context=deepcopy(context),sources={'fixture':'arithmetic-only test'})
    return definition,rows,cases


def test_complete_matched_cycle_closes_and_preserves_context(cycle):
    definition,rows,cases=cycle;out=native_cycle(*cycle)
    assert out['native_cycle_qualified'] and not out['conformer_ensemble_qualified']
    assert not out['isolated_motif_binding_energy'] and not out['electronic_coupling_computed']
    assert out['delta_delta_E_kJ_mol']==pytest.approx(-HARTREE_EV*EV_KJ)
    for c in ('25','50','100'):
        expected=math.fsum(coefficient*rows[k]['G_by_cutoff_eV'][c] for k,coefficient in definition['coefficients'].items())*EV_KJ
        assert out['delta_delta_G_by_cutoff_kJ_mol'][c]==pytest.approx(expected)
        assert out['property_assemblies'][c]['context_verified']
        assert out['property_assemblies'][c]['context']['thermal_definition']['selected_rotor_cutoff_cm']==int(c)


def test_missing_native_components_stay_blank_even_with_other_G_values(cycle):
    definition,rows,cases=cycle;rows.pop('empty');out=native_cycle(*cycle)
    assert out['missing']==['empty'] and out['delta_delta_G_by_cutoff_kJ_mol'] is None
    assert out['delta_delta_E_kJ_mol'] is None and not out['native_cycle_qualified']


@pytest.mark.parametrize('field,value',[('solvent','water'),('standard_state','different'),('electronic_protocol',{'basis':'other'}),('reference','other')])
def test_mixed_physical_context_cannot_form_cycle(cycle,field,value):
    _,rows,_=cycle;rows['first']['context'][field]=value
    with pytest.raises(ValueError,match='contexts differ'):native_cycle(*cycle)


@pytest.mark.parametrize('change',['pending','charge','spin','missing_cutoff','nonfinite','sum','ensemble','material','atoms'])
def test_unqualified_or_mismatched_components_cannot_be_combined(cycle,change):
    _,rows,cases=cycle;r=rows['first']
    if change=='pending':r['native_G_qualified']=False
    if change=='charge':r['charge']=1
    if change=='spin':r['mult']=3
    if change=='missing_cutoff':r['G_by_cutoff_eV'].pop('100')
    if change=='nonfinite':r['nuclear_by_cutoff_eV']['25']=float('nan')
    if change=='sum':r['G_by_cutoff_eV']['50']+=.01
    if change=='ensemble':r['conformer_ensemble_qualified']=True
    if change=='material':r['material_admission']=True
    if change=='atoms':cases['empty']['smiles']='CCC'
    with pytest.raises(ValueError):native_cycle(*cycle)


def test_single_checked_component_still_requires_finite_values(cycle):
    _,rows,_=cycle;r=rows['first'];r['E_geometry_level_Ha']=float('nan')
    with pytest.raises(ValueError):validate_component(r,0)
