import importlib.util
import json
from pathlib import Path
from rdkit import Chem


def test_declared_ring_panel_preserves_linker_and_unique_graphs():
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('ring_panel',root/'scripts/candidates/expand_ring_panel.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=module.prepare(json.loads((root/'config/candidate_ring_expansion.json').read_text()))
    assert len(rows)==9 and sum(len(r['recipes']) for r in rows)==17
    assert len({r['product_smiles'] for r in rows})==9
    anchor=Chem.MolFromSmarts('[#8]-[#6]-[#6]-[#7]')
    assert all(Chem.MolFromSmiles(r['product_smiles']).HasSubstructMatch(anchor) for r in rows)
    assert all(not r['properties_computed'] for r in rows)
    assert sum(r['recipes'][0]['substituent']=='H' for r in rows)==1
