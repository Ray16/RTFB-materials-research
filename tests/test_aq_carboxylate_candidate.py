import importlib.util
from pathlib import Path

import pytest
from rdkit import Chem


path = Path(__file__).resolve().parents[1] / 'scripts/pipeline/prepare_aq_carboxylate_candidate.py'
spec = importlib.util.spec_from_file_location('aq_carboxylate_candidate', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_exact_aq_acid_is_reconstructed_as_balanced_ester():
    result = module.build('C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=C(C=C3)C(=O)O')
    assert result['precursor_inchikey'] == 'ASDLSKCKYGVMAI-UHFFFAOYSA-N'
    assert result['reaction']['atom_balanced'] and result['reaction']['charge_balanced']
    assert result['reaction']['base_correction_required']
    assert result['product_formula'] == 'C23H16O4'
    assert Chem.MolFromSmiles(result['frag']) is not None
    assert len(result['route_binding']['new_bond_atom_maps']) == 2


@pytest.mark.parametrize('smiles', ['not a smiles', 'CCO', 'OC(=O)CC(=O)O'])
def test_missing_ambiguous_or_invalid_acid_is_rejected(smiles):
    with pytest.raises(ValueError):
        module.build(smiles)
