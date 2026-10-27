import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pmdi_matrix', ROOT/'scripts/pipeline/prepare_pmdi_linker_matrix.py')
matrix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matrix)


def detached_imide_core(fragment):
    molecule = Chem.MolFromSmiles(fragment)
    dummy = next(a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==0)
    nitrogens = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==7 and a.GetIsAromatic()]
    path = min((Chem.GetShortestPath(molecule,dummy,n) for n in nitrogens),key=len)
    editable = Chem.RWMol(molecule);editable.RemoveBond(path[-2],path[-1])
    nitrogen = editable.GetAtomWithIdx(path[-1]);nitrogen.SetNumExplicitHs(1);nitrogen.SetNoImplicit(True)
    detached=editable.GetMol();Chem.SanitizeMol(detached)
    cores=[part for part in Chem.GetMolFrags(detached,asMols=True,sanitizeFrags=True)
           if sum(a.GetAtomicNum()==7 for a in part.GetAtoms())==2]
    assert len(cores)==1
    return Chem.MolToSmiles(cores[0])


def test_matrix_binds_exact_products_and_reuses_existing_flexible_graphs():
    library = json.loads((ROOT/'results/discovery/systematic_families_2026-10-18/structures.json').read_text())
    rows = matrix.build(library)
    assert len(rows) == len({r['product_inchikey'] for r in rows}) == 4
    assert sum(r['origin']=='reused_existing_graph' for r in rows) == 2
    assert {(r['design_factors']['linker'],r['design_factors']['edge_substituent']) for r in rows} == {
        ('flexible','H'),('flexible','Me'),('cyclic','H'),('cyclic','Me')}
    cores = {(r['design_factors']['linker'],r['design_factors']['edge_substituent']):detached_imide_core(r['frag']) for r in rows}
    assert cores['flexible','H'] == cores['cyclic','H']
    assert cores['flexible','Me'] == cores['cyclic','Me']
    assert cores['flexible','H'] != cores['flexible','Me']
    for row in rows:
        assert row['reaction']['atom_balanced'] and row['reaction']['charge_balanced']
        assert row['route_binding']['heavy_atom_mapping'] and not row['properties_computed']
        mol = Chem.MolFromSmiles(row['route_binding']['precursor_smiles'])
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if row['design_factors']['linker']=='cyclic':
            assert len(centers)==2 and all(value!='?' for _,value in centers)
        else:
            assert centers == []


def test_unreviewed_matrix_routes_are_held():
    library = json.loads((ROOT/'results/discovery/systematic_families_2026-10-18/structures.json').read_text())
    gate_spec = importlib.util.spec_from_file_location('matrix_route_gate', ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate = importlib.util.module_from_spec(gate_spec);gate_spec.loader.exec_module(gate)
    assert all(not gate.screen(row, {}, {})['proceed_to_full_redox'] for row in matrix.build(library))


def test_matrix_all_axis_decisions_reject_mismatched_route(tmp_path):
    source = ROOT/'results/discovery/pmdi_linker_matrix_2026-10-26'
    assess_spec = importlib.util.spec_from_file_location('matrix_assess', ROOT/'scripts/pipeline/assess_pmdi_linker_matrix.py')
    module = importlib.util.module_from_spec(assess_spec);assess_spec.loader.exec_module(module)
    gate = json.loads((source/'route_screen.json').read_text())
    gate['candidates'][0]['evidence']['exact_design_key'] = 'wrong-product'
    altered = tmp_path/'gate.json';altered.write_text(json.dumps(gate))
    with pytest.raises(ValueError, match='exact design'):
        module.assess(source/'matrix.json',altered,ROOT/'config/candidate_feasibility.json')
