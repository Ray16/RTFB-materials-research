"""The reported NaH bead route is exact-graph scoped and supplies no calculated G."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

from rdkit import Chem
import pytest

from redox.workflows.ionic_stationary_cycles import validate_cycle_balance

ROOT=Path(__file__).resolve().parents[1]


def test_primary_bead_route_uses_N_attachment_and_a_balanced_hydride_reaction():
    review=json.loads((ROOT/'validation/phenothiazine_postfunctionalization_2026-11-02/reviewed_route.json').read_text())
    binding=review['route_binding'];product=Chem.MolFromSmiles(binding['mapped_product_smiles'])
    mapped={a.GetAtomMapNum():a for a in product.GetAtoms()}
    assert mapped[binding['new_bond_atom_maps'][1]].GetAtomicNum()==7
    route=review['actual_condition_stoichiometry']
    validate_cycle_balance(route,route['species'])
    assert route['species']['hydride']['smiles']=='[Na+].[H-]'
    assert route['species']['hydrogen']['smiles']=='[H][H]'
    assert route['reaction_deltaG_kJ_mol'] is None and not route['actual_condition_thermochemistry_computed']
    assert not review['full_material_admission']
    bad=deepcopy(route);bad['coefficients'].pop('hydrogen')
    with pytest.raises(ValueError):validate_cycle_balance(bad,bad['species'])


def test_unsubstituted_route_does_not_admit_ring_substitution_or_ether_attachment():
    spec=importlib.util.spec_from_file_location('primary_ptz_gate',ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
    review=json.loads((ROOT/'validation/phenothiazine_postfunctionalization_2026-11-02/reviewed_route.json').read_text())
    panel=json.loads((ROOT/'results/candidates/pool/proposals/ptz_direct_ring_2026-11-02.json').read_text())
    assessment={'candidates':{}}
    baseline=next(r for r in panel if r['id']=='phenothiazine')
    assert gate.screen(baseline,{'phenothiazine':review},assessment)['proceed_to_full_redox']
    derivative=next(r for r in panel if r['id']!='phenothiazine')
    assert not gate.screen(derivative,{derivative['id']:review},assessment)['proceed_to_full_redox']
    ether=dict(id='phenothiazine',frag='c1ccc2c(c1)Sc1ccccc1N2CCO[*:1]')
    assert not gate.screen(ether,{'phenothiazine':review},assessment)['proceed_to_full_redox']
