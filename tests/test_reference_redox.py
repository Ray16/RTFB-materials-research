from copy import deepcopy
import pytest

from redox.properties.reference_redox import neutral_oxidation_vs_fc


def components():
    context = dict(quantity='composite_solution_G', unit='eV', electronic_protocol={'model':'test'},
        thermal_definition={'native_protocol':'test'}, solvent='acetonitrile', solvent_model='SMD',
        temperature_K=298.15, standard_state='same source convention', reference='absolute energy')
    values = dict(graft_ox=-95., graft_neu=-100., fc_ox=-195.6, fc_neu=-200.)
    rows = {}
    for key, value in values.items():
        oxidized = key.endswith('_ox'); fc = key.startswith('fc_')
        rows[key] = dict(state=dict(identity=key, graph_identity='Fc' if fc else 'graft',
            charge=int(oxidized), multiplicity=2 if oxidized else 1),
            symbols=['Fe']+['C']*10+['H']*10 if fc else ['C','S','H','H'],
            context=deepcopy(context), G_by_cutoff_eV={c:value for c in ('25','50','100')},
            native_numerical_policy=dict(component='primary legacy SMD CDS',
                Hessian_gradient_step_Bohr=1e-4 if fc else 1e-6,
                energy_function_changed=False, gradient_function_changed=False, global_defaults_changed=False),
            native_G_qualified=True, thermodynamic_entities=1, qc='independently_checked',
            sources={'synthetic':'unit-test'})
    return rows


def test_oxidation_sign_reference_cancellation_and_policy_retention():
    data=components(); result=neutral_oxidation_vs_fc(data)
    assert result['potential_by_cutoff_V']['50']==pytest.approx(.6)
    assert result['one_atm_to_one_M_correction_V']==0
    assert result['components']['fc_neu']['native_numerical_policy']['Hessian_gradient_step_Bohr']==1e-4
    assert result['components']['graft_neu']['native_numerical_policy']['Hessian_gradient_step_Bohr']==1e-6
    for row in data.values():
        row['G_by_cutoff_eV']={c:v+123. for c,v in row['G_by_cutoff_eV'].items()}
    assert neutral_oxidation_vs_fc(data)['potential_by_cutoff_V']['50']==pytest.approx(.6)


@pytest.mark.parametrize('change', ['solvent','spin','atoms','reference','graph','policy','pending','nan'])
def test_invalid_components_do_not_produce_a_potential(change):
    data=components()
    if change=='solvent':data['graft_ox']['context']['solvent']='water'
    if change=='spin':data['graft_ox']['state']['multiplicity']=1
    if change=='atoms':data['graft_ox']['symbols'].append('H')
    if change=='reference':
        for key in ('fc_neu','fc_ox'):data[key]['symbols'][0]='Ru'
    if change=='graph':data['graft_ox']['state']['graph_identity']='different graft'
    if change=='policy':data['fc_neu']['native_numerical_policy']['energy_function_changed']=True
    if change=='pending':data['graft_ox']['native_G_qualified']=False
    if change=='nan':data['graft_neu']['G_by_cutoff_eV']['50']=float('nan')
    with pytest.raises(ValueError):neutral_oxidation_vs_fc(data)
