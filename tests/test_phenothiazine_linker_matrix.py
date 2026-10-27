import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ptz_matrix',ROOT/'scripts/pipeline/prepare_phenothiazine_linker_matrix.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def inputs():
    design=json.loads((ROOT/'config/phenothiazine_linker_matrix_2026-10-26.json').read_text())
    baseline=json.loads((ROOT/'results/discovery/ptz_hydroxyethyl_source_2026-10-26/candidates.json').read_text())
    ledger=json.loads((ROOT/'results/selection/campaigns/20261024_UMA_default/decisions.json').read_text())['candidates']
    return design,baseline,ledger


def detached_core(fragment):
    molecule=Chem.MolFromSmiles(fragment);dummy=next(a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==0)
    nitrogens=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==7 and a.IsInRing()]
    path=min((Chem.GetShortestPath(molecule,dummy,n) for n in nitrogens),key=len)
    edit=Chem.RWMol(molecule);edit.RemoveBond(path[-2],path[-1]);atom=edit.GetAtomWithIdx(path[-1]);atom.SetNumExplicitHs(1);atom.SetNoImplicit(True)
    detached=edit.GetMol();Chem.SanitizeMol(detached)
    cores=[p for p in Chem.GetMolFrags(detached,asMols=True,sanitizeFrags=True) if any(a.GetAtomicNum()==16 for a in p.GetAtoms())]
    assert len(cores)==1
    return Chem.MolToSmiles(cores[0])


def linker_stereochemistry(smiles):
    molecule=Chem.MolFromSmiles(smiles);Chem.AssignStereochemistry(molecule,cleanIt=True,force=True)
    result={}
    for atom in molecule.GetAtoms():
        if not atom.HasProp('_CIPCode'):continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() in ('O','N'):
                result[neighbor.GetSymbol()]=atom.GetProp('_CIPCode')
    return result


def test_matched_products_reuse_baseline_and_preserve_core_identity():
    rows=module.build(*inputs());assert len(rows)==len({r['product_inchikey'] for r in rows})==4
    assert sum(r['origin']=='reused_source_proposal' for r in rows)==1 and sum(r['origin']=='new_exact_graph' for r in rows)==3
    cores={(r['design_factors']['linker'],r['design_factors']['substituent']):detached_core(r['frag']) for r in rows}
    assert cores['flexible','H']==cores['cyclic','H']
    assert cores['flexible','3_7_dimethoxy']==cores['cyclic','3_7_dimethoxy']
    assert cores['flexible','H']!=cores['flexible','3_7_dimethoxy']
    feed=json.loads((ROOT/'results/candidates/current.json').read_text())
    prior=next(r for r in feed['original39']+feed['additional_candidates'] if r['id']=='exp_ptz_37dimethoxy')
    assert detached_core(prior['frag'])==cores['flexible','3_7_dimethoxy']
    stereos={r['design_factors']['substituent']:linker_stereochemistry(r['route_binding']['precursor_smiles'])
             for r in rows if r['design_factors']['linker']=='cyclic'}
    assert stereos['H']==stereos['3_7_dimethoxy'] and set(stereos['H'])=={'O','N'}
    for row in rows:
        assert row['reaction']['atom_balanced'] and row['reaction']['charge_balanced'] and row['route_binding']['heavy_atom_mapping']
        assert len(row['states'])==2 and not row['properties_computed']
        centers=Chem.FindMolChiralCenters(Chem.MolFromSmiles(row['route_binding']['precursor_smiles']),includeUnassigned=True)
        assert len(centers)==(2 if row['design_factors']['linker']=='cyclic' else 0)
        assert all(label!='?' for _,label in centers)


def test_missing_design_factor_is_rejected():
    design,baseline,ledger=inputs();design['anchors'].pop('cyclic')
    with pytest.raises(ValueError,match='design factors'):module.build(design,baseline,ledger)
