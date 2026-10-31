import copy
import json
from itertools import combinations
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from redox.build.ionic_polymer import ionic_definition,capacity_bound,charged_panel
from redox.build.polymer_cases import build_case,panel
from redox.build.oligomer_topology import chain_topology

FRAGMENT='[*:1][n+]1ccc(-c2cc[n+](C)cc2)cc1'
CONFIG=json.loads((Path(__file__).parents[1]/'config/viologen_polymer_preparation.json').read_text())


def spec():return dict(copy.deepcopy(CONFIG),backbone_stereo_tags=CONFIG['configurations']['all_same'])


def test_charged_graph_panel_is_balanced_and_does_not_enter_neutral_QM_path():
    cases,cycles,definition=charged_panel(spec(),FRAGMENT)
    assert len(cases)==6 and len(cycles)==2 and definition['charge_per_graft']==2
    for case in cases:
        assert case['charge']==2*len(case['occupied_sites'])
        assert case['site_states'][0]==case['site_states'][-1]=='styrene'
        assert not case['quantum_calculation_admitted'] and not case['native_G_qualified']
        assert not case['charge_or_spin_localization_established']
        topology=chain_topology(case['smiles'],5,[2,3,4])
        assert topology['backbone_carbon_count']==10 and topology['backbone_CC_bonds']==9 and topology['chain_ends']==2
        assert topology['para_graft_sites'][0]['benzyl_atom'] is None
    with pytest.raises(ValueError,match='counterion'):panel(spec(),FRAGMENT)


def test_quaternization_mass_and_counterions_match_every_actual_occupancy_graph():
    model=CONFIG['ionic_model'];pf6=Descriptors.MolWt(Chem.MolFromSmiles('F[P-](F)(F)(F)(F)F'))
    for g in range(4):
        masses=[]
        for occupied in combinations([2,3,4],g):
            case=build_case(FRAGMENT,5,[2,3,4],list(occupied),['CW']*4)
            result=capacity_bound(FRAGMENT,model,5,3,g,2)
            assert result['polymer_cation_mass_g_mol']==pytest.approx(Descriptors.MolWt(Chem.MolFromSmiles(case['smiles'])),abs=1e-8)
            assert result['counterion_mass_g_mol']==pytest.approx(2*g*pf6)
            assert result['mass_g_mol']==pytest.approx(result['polymer_cation_mass_g_mol']+2*g*pf6)
            assert result['bare_polymer_charge']==2*g and result['reference_salt_charge']==0
            assert not result['accessible_capacity_established'];masses.append(result['mass_g_mol'])
        assert max(masses)-min(masses)<1e-8


def test_inactive_mass_and_electron_budget_are_explicit_scenarios():
    one=capacity_bound(FRAGMENT,CONFIG['ionic_model'],5,3,2,1)
    two=capacity_bound(FRAGMENT,CONFIG['ionic_model'],5,3,2,2)
    diluted=capacity_bound(FRAGMENT,CONFIG['ionic_model'],5,3,2,2,100.)
    assert two['capacity_upper_bound_mAh_g']==pytest.approx(2*one['capacity_upper_bound_mAh_g'])
    assert diluted['capacity_upper_bound_mAh_g']<two['capacity_upper_bound_mAh_g']
    assert not two['mobile_ion_mass_during_redox_established']


@pytest.mark.parametrize('change',['no_convention','missing_ions','wrong_charge','neutral_ion','wrong_parity','missing_state','duplicate_handles'])
def test_undefined_ion_or_spin_models_cannot_silently_pass(change):
    config=spec()
    if change=='no_convention':config['ionic_model']['calculation_convention']='implicit'
    if change=='missing_ions':config['ionic_model']['composition_reference']['counterions']=[]
    if change=='wrong_charge':config['ionic_model']['composition_reference']['counterions'][0]['per_graft']=1
    if change=='neutral_ion':config['ionic_model']['composition_reference']['counterions'][0]['smiles']='Cl'
    if change=='wrong_parity':config['ionic_model']['proposed_global_states_by_graft_count']['2'][0]['mult']=2
    if change=='missing_state':config['ionic_model']['proposed_global_states_by_graft_count']['1']=[]
    if change=='duplicate_handles':config['reactive_sites'].append(2)
    with pytest.raises(ValueError):charged_panel(config,FRAGMENT)
