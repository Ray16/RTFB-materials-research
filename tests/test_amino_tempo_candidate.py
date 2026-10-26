import importlib.util
from pathlib import Path

import pytest
from rdkit import Chem


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/pipeline/prepare_amino_tempo_candidate.py'
spec = importlib.util.spec_from_file_location('amino_tempo_candidate', SCRIPT)
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)

IDENTITY = dict(CID=550942, InChIKey='XUXUHDYTLNCYQQ-UHFFFAOYSA-N',
                ConnectivitySMILES='CC1(CC(CC(N1[O])(C)C)N)C')


def test_exact_mapped_graft_preserves_radical_and_mass_balance():
    row = candidate.build(IDENTITY)
    assert row['reaction']['atom_balanced'] and row['reaction']['charge_balanced']
    assert row['route_binding']['heavy_atom_mapping']
    product = Chem.MolFromSmiles(row['route_binding']['product_smiles'])
    assert sum(a.GetNumRadicalElectrons() for a in product.GetAtoms()) == 1
    assert row['precursor_inchikey'] == IDENTITY['InChIKey']
    assert row['product_formula'] == 'C17H27N2O'


@pytest.mark.parametrize('field,value', [
    ('CID', 137994), ('InChIKey', 'bad'), ('ConnectivitySMILES', 'CCO'),
])
def test_rejects_wrong_source_identity(field, value):
    identity = {**IDENTITY, field: value}
    with pytest.raises(ValueError):
        candidate.build(identity)


def test_unreviewed_route_cannot_enter_full_redox():
    gate_path = SCRIPT.parent / 'screen_grafting_routes.py'
    gate_spec = importlib.util.spec_from_file_location('route_gate', gate_path)
    gate = importlib.util.module_from_spec(gate_spec)
    gate_spec.loader.exec_module(gate)
    result = gate.screen(candidate.build(IDENTITY), {}, {})
    assert result['route_status'] == 'hold'
    assert result['proceed_to_full_redox'] is False
    assert 'conditions_defined' in result['missing_route_evidence']
