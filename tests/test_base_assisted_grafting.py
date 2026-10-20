"""Thermodynamic-cycle invariants for the base sensitivity pilot."""
import importlib.util
import math
from pathlib import Path

from rdkit import Chem
from redox.build.grafting import atom_inventory
from collections import Counter

spec=importlib.util.spec_from_file_location('base_pilot',Path(__file__).resolve().parents[1]/'scripts/pipeline/base_assisted_grafting.py')
P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)


def test_protonation_and_pair_reactions_are_atom_charge_balanced():
    hcl=Chem.MolFromSmiles('Cl');chloride=Chem.MolFromSmiles('[Cl-]')
    for identity in P.BASES.values():
        neutral=Chem.MolFromSmiles(identity['neutral']);protonated=Chem.MolFromSmiles(identity['protonated'])
        pair=Chem.MolFromSmiles(identity['protonated']+'.[Cl-]')
        assert atom_inventory(neutral)+atom_inventory(hcl)==atom_inventory(pair)
        assert atom_inventory(protonated)+atom_inventory(chloride)==atom_inventory(pair)
        assert Chem.GetFormalCharge(neutral)+Chem.GetFormalCharge(hcl)==Chem.GetFormalCharge(pair)==0
        assert Chem.GetFormalCharge(protonated)+Chem.GetFormalCharge(chloride)==0


def test_concentration_convention_counts_independent_solutes():
    assert P.concentration_correction(0)==0
    assert math.isclose(P.concentration_correction(-1),-7.9259,abs_tol=.002)
    # Dissociation and association corrections cancel in a closed cycle.
    assert P.concentration_correction(1)+P.concentration_correction(-1)==0


def test_contact_pair_and_free_ion_cycles_close_and_missing_is_not_zero():
    e={'dbu_neutral':-100.,'dbu_protonated':-101.,'dbu_chloride_pair':-110.5}
    result=P.corrections(e,'dbu',-9.,-9.2)
    assert math.isclose(result['contact_pair_kJ_mol'],result['free_ions_kJ_mol']+result['ion_pair_association_kJ_mol'],abs_tol=1e-9)
    del e['dbu_protonated']
    result=P.corrections(e,'dbu',-9.,-9.2)
    assert result['free_ions_kJ_mol'] is None and result['ion_pair_association_kJ_mol'] is None
    assert result['contact_pair_kJ_mol'] is not None
