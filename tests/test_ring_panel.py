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


def test_direct_nitrogen_panel_has_exact_mapped_pairs_without_inheriting_ether_route():
    root=Path(__file__).resolve().parents[1]
    def load(name,path):
        spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
    panel=load('direct_nitrogen_panel',root/'scripts/candidates/expand_ring_panel.py')
    gate=load('independent_direct_nitrogen_gate',root/'scripts/pipeline/screen_grafting_routes.py')
    rows=panel.prepare(json.loads((root/'config/candidate_ring_expansion_direct_ptz.json').read_text()))
    ether=panel.prepare(json.loads((root/'config/candidate_ring_expansion.json').read_text()))
    assert len(rows)==9 and sum(len(r['recipes']) for r in rows)==17
    assert not {r['product_smiles'] for r in rows}&{r['product_smiles'] for r in ether}
    assert next(r for r in rows if r['id']=='phenothiazine')['recipes'][0]['substituent']=='H'
    for row in rows:
        assert row['states']==[['neu',0,1,0],['ox',1,2,1]] and not row['properties_computed']
        binding=row['route_binding'];assert binding==gate.exact_route_binding(row)
        mapped=Chem.MolFromSmiles(binding['mapped_product_smiles'])
        atoms={a.GetAtomMapNum():a for a in mapped.GetAtoms()}
        site,anchor=binding['new_bond_atom_maps']
        assert atoms[site].GetAtomicNum()==6 and atoms[anchor].GetAtomicNum()==7
        assert mapped.GetBondBetweenAtoms(atoms[site].GetIdx(),atoms[anchor].GetIdx()) is not None
        assert row['reference_reaction']['atom_balanced'] and row['reference_reaction']['charge_balanced']
