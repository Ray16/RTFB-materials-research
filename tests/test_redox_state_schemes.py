"""Chemical state drawings must conserve the molecule and represent the added electrons."""
import importlib.util
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def module():
    path=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates/redox_state_schemes.py'
    spec=importlib.util.spec_from_file_location('redox_state_drawings',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def test_quinone_reduction_is_not_accidentally_hydrogenation():
    states=module().reduction_forms('O=C1C=CC(=O)C=C1')
    assert [rdMolDescriptors.CalcMolFormula(m) for m in states]==['C6H4O2','C6H4O2-','C6H4O2-2']
    assert [sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()) for m in states]==[0,1,0]
    assert all(a.GetIsAromatic() for a in states[2].GetAtoms() if a.GetAtomicNum()==6)


def test_viologen_state_charges_are_consistent_with_quinoid_neutral_form():
    states=module().reduction_forms('C[n+]1ccc(-c2cc[n+](C)cc2)cc1')
    assert [Chem.GetFormalCharge(m) for m in states]==[2,1,0]
    assert [rdMolDescriptors.CalcMolFormula(m) for m in states]==['C12H14N2+2','C12H14N2+','C12H14N2']
    expected=Chem.MolFromSmiles('CN1C=CC(=C2C=CN(C)C=C2)C=C1')
    assert Chem.MolToSmiles(states[2])==Chem.MolToSmiles(expected)


def test_diimide_forms_preserve_alkyl_caps_and_all_carbonyl_oxygens():
    states=module().reduction_forms('Cn1c(=O)c2cc3c(=O)n(C)c(=O)c3cc2c1=O')
    assert [rdMolDescriptors.CalcMolFormula(m) for m in states]==['C12H8N2O4','C12H8N2O4-','C12H8N2O4-2']
    assert sum(a.GetAtomicNum()==8 for a in states[2].GetAtoms())==4


def test_unrecognized_chemistry_cannot_receive_a_fabricated_redox_diagram():
    with pytest.raises(ValueError,match='No supported valid'):
        module().reduction_forms('CC')
