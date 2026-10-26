import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('local_counterpoise_test',ROOT/'scripts/validation/local_counterpoise.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_ghosts_add_basis_without_nuclei_or_electrons():
    coords=[[0,0,0],[0,.76,.58],[0,-.76,.58],[0,0,3],[0,.76,3.58],[0,-.76,3.58]]
    mols=module.molecules(['O','H','H']*2,coords,[0,1,2])
    assert mols['AB'].nelectron==20
    assert mols['A_ghost'].nelectron==10 and list(mols['A_ghost'].atom_charges())==[8,1,1,0,0,0]
    assert mols['A_ghost'].nao_nr()==mols['AB'].nao_nr()
    assert mols['A_own'].nao_nr()*2==mols['AB'].nao_nr()


def test_counterpoise_shift_is_added_once_and_gas_only():
    row=module.sums(dict(AB=-200.02,A_own=-100.,B_own=-100.,A_ghost=-100.001,B_ghost=-100.002))
    assert row['counterpoise_gas_kJ_mol']-row['uncorrected_common_grid_kJ_mol']==pytest.approx(row['AO_basis_expansion_shift_kJ_mol'])
    assert row['AO_basis_expansion_shift_kJ_mol']>0 and not row['solvent_counterpoise_computed']


@pytest.mark.parametrize('change',['missing','nan'])
def test_missing_or_nonfinite_counterpoise_terms_are_refused(change):
    row=dict(AB=-200.02,A_own=-100.,B_own=-100.,A_ghost=-100.001,B_ghost=-100.002)
    if change=='missing':row.pop('A_ghost')
    else:row['AB']=float('nan')
    with pytest.raises(ValueError):module.sums(row)
