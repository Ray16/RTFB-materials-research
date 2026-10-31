import copy

import pytest

from redox.build.grafting import grafting_reaction
from redox.core.common import EV_KJ
from redox.properties.grafting import quaternization_delta_g


@pytest.fixture
def cycle():
    reaction=grafting_reaction('[*:1][n+]1ccc(-c2cc[n+](C)cc2)cc1',
        'Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1')
    from rdkit import Chem
    values={'product':-100.125,'byproduct':-20.,'precursor':-70.,'site':-50.}
    rows={role:dict(smiles=reaction[role+'_smiles'],charge=Chem.GetFormalCharge(Chem.MolFromSmiles(reaction[role+'_smiles'])),
        G_by_cutoff_eV={c:value for c in ['25','50','100']}) for role,value in values.items()}
    return reaction,rows


def test_chloride_quaternization_has_no_hcl_or_base_shift(cycle):
    result=quaternization_delta_g(*cycle)
    assert result['deltaG_kJ_mol_by_cutoff']['50']==pytest.approx(-.125*EV_KJ)
    assert result['delta_n_solutes']==0
    assert result['standard_state_correction_kJ_mol']==0
    assert result['base_term_kJ_mol']==0
    assert result['material_admission'] is False


@pytest.mark.parametrize('change',['HCl','wrong_charge','wrong_core','missing_G','nan','base_route','forged_balance'])
def test_incompatible_species_and_missing_thermochemistry_cannot_enter_cycle(cycle,change):
    reaction,rows=copy.deepcopy(cycle)
    if change=='HCl':rows['byproduct']['smiles']='Cl'
    if change=='wrong_charge':rows['product']['charge']=1
    if change=='wrong_core':rows['precursor']['smiles']='c1ccncc1'
    if change=='missing_G':del rows['site']['G_by_cutoff_eV']['25']
    if change=='nan':rows['product']['G_by_cutoff_eV']['100']=float('nan')
    if change=='base_route':reaction['base_correction_required']=True
    if change=='forged_balance':reaction['site_smiles']=rows['site']['smiles']='CCl'
    with pytest.raises(ValueError):quaternization_delta_g(reaction,rows)
