import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('pmdi_hydroxyethyl',ROOT/'scripts/pipeline/prepare_pmdi_hydroxyethyl_control.py')
candidate=importlib.util.module_from_spec(spec);spec.loader.exec_module(candidate)


def test_exact_o_attachment_preserves_competing_imide_handle():
    row=candidate.build()
    assert row['precursor_formula']=='C12H8N2O5'
    assert row['precursor_inchikey']==candidate.PRECURSOR_KEY
    assert row['reaction']['atom_balanced'] and row['reaction']['charge_balanced']
    for smiles in (row['route_binding']['precursor_smiles'],row['product_smiles']):
        assert len(Chem.MolFromSmiles(smiles).GetSubstructMatches(Chem.MolFromSmarts('[nH]')))==1
    product=Chem.MolFromSmiles(row['route_binding']['mapped_product_smiles'])
    maps={a.GetAtomMapNum():a.GetIdx() for a in product.GetAtoms()}
    a,b=row['route_binding']['new_bond_atom_maps']
    assert {product.GetAtomWithIdx(maps[a]).GetSymbol(),product.GetAtomWithIdx(maps[b]).GetSymbol()}=={'C','O'}
    assert product.GetBondBetweenAtoms(maps[a],maps[b]) is not None


def test_access_report_does_not_bypass_missing_selectivity():
    gate_spec=importlib.util.spec_from_file_location('route_gate',ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate=importlib.util.module_from_spec(gate_spec);gate_spec.loader.exec_module(gate)
    row=candidate.build()
    review={row['id']:{'exact_precursor_access_reviewed':True}}
    result=gate.screen(row,review,{})
    assert result['route_status']=='hold' and not result['proceed_to_full_redox']
    assert 'attachment_selectivity_reviewed' in result['missing_route_evidence']


def test_all_axis_assessment_refuses_foreign_route_identity(tmp_path):
    spec=importlib.util.spec_from_file_location('source_control_assess',ROOT/'scripts/pipeline/assess_pmdi_hydroxyethyl_control.py')
    assessor=importlib.util.module_from_spec(spec);spec.loader.exec_module(assessor)
    source=ROOT/'results/discovery/pmdi_hydroxyethyl_source_control_2026-10-26'
    gate=json.loads((source/'route_screen.json').read_text())
    gate['candidates'][0]['evidence']['exact_design_key']='different-product'
    altered=tmp_path/'gate.json';altered.write_text(json.dumps(gate))
    with pytest.raises(ValueError,match='exact candidate'):
        assessor.assess(source/'candidates.json',altered,ROOT/'config/candidate_feasibility.json')
