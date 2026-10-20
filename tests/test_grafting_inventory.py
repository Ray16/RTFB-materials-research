"""Independent stoichiometric and site-identity checks for the new model inventories."""
from collections import Counter
import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem

from redox.build.grafting import grafting_reaction, precursor
from redox.build.oligomers import decorate_sites
from redox.build.structures import decorate, unassigned_stereo

ROOT = Path(__file__).resolve().parents[1]


def counts(smi):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    return Counter(a.GetAtomicNum() for a in mol.GetAtoms())


@pytest.mark.parametrize('fragment,kind,precursor_smi', [
    ('[*:1]Oc1ccccc1', 'NH_OH_substitution', 'Oc1ccccc1'),
    ('[*:1]N1C(=O)CCC1=O', 'NH_OH_substitution', 'O=C1CCC(=O)N1'),
    ('[*:1][n+]1ccccc1', 'N_quaternization', 'c1ccncc1'),
])
def test_precursors_and_balanced_reactions(fragment, kind, precursor_smi):
    product = Chem.MolToSmiles(decorate('Cc1ccc(C[*:1])cc1', fragment))
    r = grafting_reaction(fragment, product)
    assert r['reaction_class'] == kind
    assert r['precursor_smiles'] == Chem.MolToSmiles(Chem.MolFromSmiles(precursor_smi))
    assert counts(r['precursor_smiles'])+counts(r['site_smiles']) == counts(r['product_smiles'])+counts(r['byproduct_smiles'])
    q = lambda s: Chem.GetFormalCharge(Chem.MolFromSmiles(s))
    assert q(r['precursor_smiles'])+q(r['site_smiles']) == q(r['product_smiles'])+q(r['byproduct_smiles'])


def test_isomeric_wrong_attachment_is_not_accepted_by_formula():
    with pytest.raises(ValueError, match='attachment site'):
        grafting_reaction('[*:1]Oc1ccccc1', 'Cc1cccc(COc2ccccc2)c1')


def test_all_current_reactions_preserve_declared_graphs():
    feed = json.loads((ROOT/'results/candidates/current.json').read_text())
    rows = feed['original39']+feed['additional_candidates']
    assert len(rows) == 75
    for row in rows:
        r = grafting_reaction(row['frag'], row['exact_grafted_smiles'])
        assert r['precursor_smiles'] == Chem.MolToSmiles(Chem.MolFromSmiles(row['precursor_smiles']))


@pytest.mark.parametrize('fragment', ['[*:1]Oc1ccccc1', '[*:1][n+]1ccccc1'])
def test_two_site_cooperativity_cycle_is_balanced(fragment):
    AA = decorate_sites(fragment, 2, (1, 2))
    ClCl = decorate_sites(fragment, 2, ())
    ACl = decorate_sites(fragment, 2, (1,))
    ClA = decorate_sites(fragment, 2, (2,))
    assert counts(AA[0])+counts(ClCl[0]) == counts(ACl[0])+counts(ClA[0])
    charge = lambda smi: Chem.GetFormalCharge(Chem.MolFromSmiles(smi))
    assert charge(AA[0])+charge(ClCl[0]) == charge(ACl[0])+charge(ClA[0])
    for smi in AA+ClCl+ACl+ClA:
        mol = Chem.MolFromSmiles(smi)
        assert len(Chem.GetMolFrags(mol)) == 1
        assert unassigned_stereo(mol) == 0
        assert all(a.GetAtomicNum() != 0 for a in mol.GetAtoms())
    # Hydrogen-capped covalent dimer cannot be treated as two monomer fragments.
    monomer = decorate_sites(fragment, 1, (1,))[0]
    assert counts(AA[0]) != counts(monomer)+counts(monomer)


def test_reference_energy_sign_and_missing_components():
    spec = importlib.util.spec_from_file_location('collector', ROOT/'scripts/pipeline/collect_grafting_results.py')
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    r = dict(id='example', precursor_id='precursor', reaction_class='NH_OH_substitution')
    components = dict(precursor=-7., merrifield_local_chloride=-10., hydrogen_chloride=-2., chloride=-1.)
    assert mod.reaction_delta_g(r, -16., components) == (-1., [])
    r['reaction_class']='N_quaternization'
    assert mod.reaction_delta_g(r, -16., components) == (0., [])
    assert mod.reaction_delta_g(r, None, components) == (None, ['product_example'])
    del components['chloride']
    assert mod.reaction_delta_g(r, -16., components) == (None, ['chloride'])
