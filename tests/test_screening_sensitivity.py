"""Reference cancellation, signs and unsupported chemical-state rejection."""
import copy
import json

import pytest

from redox.core.common import ROOT
from redox.validation.screening_sensitivity import reconstruct,molecular_derivatives


@pytest.fixture
def states():
    rows=[]
    for key,alias,q in [('fo','benchmark:ferrocene:ox',1),('fr','benchmark:ferrocene:neu',0),
                         ('ox','benchmark:toy:neu',0),('red','benchmark:toy:red1',-1)]:
        rows.append(dict(key=key,aliases=[alias],charge=q,stored_Ha=-1.,installed_Ha=-1.,corrected_Ha=-1.))
    return rows


def test_shifted_reference_and_endpoint_contributions_have_the_correct_sign(states):
    endpoint=[dict(id='toy',event='neu->red1')]
    states[0]['corrected_Ha']+=.125
    result=reconstruct(states,endpoint,32.)
    assert result['endpoint_shifts'][0]['electronic_contribution_to_potential_shift_V']['corrected_minus_installed_V']==-4.
    states[2]['corrected_Ha']+=.25
    result=reconstruct(states,endpoint,32.)
    assert result['endpoint_shifts'][0]['electronic_contribution_to_potential_shift_V']['corrected_minus_installed_V']==4.


def test_common_energy_offset_cancels_exactly(states):
    for row in states:row['corrected_Ha']+=.5
    result=reconstruct(states,[dict(id='toy',event='neu->red1')],32.)
    assert all(v==0 for v in result['maximum_electronic_potential_shifts_V'].values())
    assert not result['new_redox_potentials_qualified']


def test_same_state_pair_from_two_studies_keeps_both_endpoints(states):
    rows=[dict(id='toy',event='neu->red1',doi=doi) for doi in ['study-A','study-B']]
    result=reconstruct(states,rows,32.)
    assert result['endpoints']==2 and [r['doi'] for r in result['endpoint_shifts']]==['study-A','study-B']


@pytest.mark.parametrize('change',['missing','duplicate','nonfinite','two_electrons'])
def test_missing_or_changed_electronic_state_cannot_enter_replay(states,change):
    rows=copy.deepcopy(states)
    if change=='missing':rows.pop()
    if change=='duplicate':rows[1]['aliases'].append(rows[0]['aliases'][0])
    if change=='nonfinite':rows[0]['corrected_Ha']=float('nan')
    if change=='two_electrons':rows[-1]['charge']=-2
    with pytest.raises(ValueError):reconstruct(rows,[dict(id='toy',event='neu->red1')],32.)


def test_actual_four_molecular_controls_replay_all_saved_columns():
    data=json.loads((ROOT/'validation/screening_chemical_impact_2026-10-31/molecular_inputs.json').read_text())
    result=molecular_derivatives(data['cases'])
    assert result['cases']==4 and result['checks']==144 and result['maximum_relative_Hessian_error']<.001


@pytest.mark.parametrize('change',['force','Hessian','displacement','missing_axis','unconverged'])
def test_altered_saved_quantum_derivative_evidence_is_rejected(change):
    data=json.loads((ROOT/'validation/screening_chemical_impact_2026-10-31/molecular_inputs.json').read_text());data['cases']=data['cases'][:1]
    variant=data['cases'][0]['variants']['corrected'];point=variant['points'][1]
    if change=='force':point['gradient_Ha_Bohr'][0][0]+=.01
    if change=='Hessian':variant['Hessian'][0][0]+=.01
    if change=='displacement':point['coordinates_Bohr'][0][0]+=.1
    if change=='missing_axis':point['gradient_Ha_Bohr'][0].pop()
    if change=='unconverged':point['SCF_converged']=False
    with pytest.raises(ValueError):molecular_derivatives(data['cases'])
