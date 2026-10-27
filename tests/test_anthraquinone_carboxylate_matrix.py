import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('aq_matrix',ROOT/'scripts/pipeline/prepare_anthraquinone_carboxylate_matrix.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def inputs():
    design=json.loads((ROOT/'config/anthraquinone_carboxylate_matrix_2026-10-26.json').read_text())
    baseline=json.loads((ROOT/'results/discovery/aq_carboxylate_2026-10-26/candidates_mapped.json').read_text())
    ledger=json.loads((ROOT/'results/selection/campaigns/20261024_UMA_default/decisions.json').read_text())['candidates']
    return design,baseline,ledger


def detached_core(fragment):
    molecule=Chem.MolFromSmiles(fragment);dummy=next(a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==0)
    aromatic=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==6 and a.GetIsAromatic()]
    path=min((Chem.GetShortestPath(molecule,dummy,n) for n in aromatic),key=len)
    edit=Chem.RWMol(molecule);edit.RemoveBond(path[-2],path[-1]);detached=edit.GetMol();Chem.SanitizeMol(detached)
    carbonyl=Chem.MolFromSmarts('[#6]=[O]')
    cores=[p for p in Chem.GetMolFrags(detached,asMols=True,sanitizeFrags=True) if len(p.GetSubstructMatches(carbonyl))==2]
    assert len(cores)==1
    return Chem.MolToSmiles(cores[0])


def linker_stereochemistry(smiles):
    molecule=Chem.MolFromSmiles(smiles);Chem.AssignStereochemistry(molecule,cleanIt=True,force=True);result={}
    for atom in molecule.GetAtoms():
        if not atom.HasProp('_CIPCode'):continue
        key='aryl' if any(n.GetIsAromatic() for n in atom.GetNeighbors()) else 'carboxyl'
        result[key]=atom.GetProp('_CIPCode')
    return result


def test_matched_products_reuse_baseline_and_preserve_quinone_core():
    rows=module.build(*inputs());assert len(rows)==len({r['product_inchikey'] for r in rows})==4
    assert sum(r['origin']=='reused_source_proposal' for r in rows)==1 and sum(r['origin']=='new_exact_graph' for r in rows)==3
    cores={(r['design_factors']['linker'],r['design_factors']['substituent']):detached_core(r['frag']) for r in rows}
    assert cores['direct','H']==cores['cyclic','H'] and cores['direct','Me']==cores['cyclic','Me']
    assert cores['direct','H']!=cores['direct','Me']
    stereos={r['design_factors']['substituent']:linker_stereochemistry(r['route_binding']['precursor_smiles'])
             for r in rows if r['design_factors']['linker']=='cyclic'}
    assert stereos['H']==stereos['Me'] and set(stereos['H'])=={'aryl','carboxyl'}
    carbonyl=Chem.MolFromSmarts('[#6]=[O]')
    for row in rows:
        assert row['reaction']['atom_balanced'] and row['reaction']['charge_balanced'] and row['route_binding']['heavy_atom_mapping']
        assert len(row['states'])==2 and not row['properties_computed']
        assert len(Chem.MolFromSmiles(row['product_smiles']).GetSubstructMatches(carbonyl))==3


def test_missing_design_factor_is_rejected():
    design,baseline,ledger=inputs();design['anchors'].pop('cyclic')
    with pytest.raises(ValueError,match='design factors'):module.build(design,baseline,ledger)


def test_shared_panel_gate_keeps_aq_pending_and_rejects_foreign_route(tmp_path):
    source=ROOT/'results/discovery/anthraquinone_carboxylate_matrix_2026-10-26'
    gate_spec=importlib.util.spec_from_file_location('aq_panel_gate',ROOT/'scripts/pipeline/assess_matched_route_panel.py')
    gate=importlib.util.module_from_spec(gate_spec);gate_spec.loader.exec_module(gate)
    policy=ROOT/'config/candidate_feasibility.json';result=gate.assess(source/'matrix.json',source/'route_screen.json',policy)
    assert len(result['decisions'])==4 and all(x['status']=='pending' for x in result['decisions'])
    rows=json.loads((source/'route_screen.json').read_text());rows['candidates'][0]['evidence']['exact_design_key']='wrong'
    altered=tmp_path/'route.json';altered.write_text(json.dumps(rows))
    with pytest.raises(ValueError,match='exact candidate'):gate.assess(source/'matrix.json',altered,policy)
