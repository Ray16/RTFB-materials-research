import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('local_motif_fragments_test',ROOT/'scripts/pipeline/local_motif_fragments.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
TEMPLATE='[*:1]N1C(=O)c2ccccc2C1=O'


@pytest.fixture
def parent():
    data=json.loads((ROOT/'results/figures/polymer_models/spacing_conformers.json').read_text())['states'][0]
    mol=Chem.AddHs(Chem.MolFromSmiles(data['task']['smiles']));conf=Chem.Conformer(mol.GetNumAtoms())
    for i,r in enumerate(data['coordinates_A']):conf.SetAtomPosition(i,r)
    mol.AddConformer(conf);return mol


def test_exact_parent_cores_and_hydrogens_are_retained_without_relaxation(parent):
    before=np.array(parent.GetConformer().GetPositions());parts=module.fragments(parent,TEMPLATE)
    assert len(parts)==2
    for p in parts:
        mol=p['mol'];coords=mol.GetConformer().GetPositions()
        assert mol.GetNumAtoms()==19 and Chem.GetFormalCharge(mol)==0
        assert Chem.MolToSmiles(Chem.RemoveHs(mol))==Chem.MolToSmiles(Chem.MolFromSmiles('CN1C(=O)c2ccccc2C1=O'))
        assert np.array_equal(coords[:-1],before[p['parent_indices']])
        assert Chem.GetMolFrags(mol)==(tuple(range(19)),)
        anchor=p['parent_to_fragment'][p['anchor_parent_index']]
        assert np.linalg.norm(coords[p['cap_index']]-coords[anchor])==pytest.approx(1.09)
        direction=before[p['cut_parent_bond'][1]]-before[p['cut_parent_bond'][0]]
        assert np.dot(coords[p['cap_index']]-coords[anchor],direction)>0
    assert np.array_equal(parent.GetConformer().GetPositions(),before)


def test_parent_motif_pair_keeps_relative_spatial_contacts(parent):
    from redox.build.polymer_contacts import contact_metrics,pi_core_groups
    parts=module.fragments(parent,TEMPLATE);pair=Chem.CombineMols(*(p['mol'] for p in parts))
    old=contact_metrics(parent.GetConformer().GetPositions(),pi_core_groups(parent,TEMPLATE))[0]
    new=contact_metrics(pair.GetConformer().GetPositions(),pi_core_groups(pair,TEMPLATE))[0]
    assert len(Chem.GetMolFrags(pair))==2
    for key in ['centroid_distance_A','interplane_angle_deg','plane_separation_A','nearest_core_atom_A']:
        assert new[key]==pytest.approx(old[key],abs=1e-10)


@pytest.mark.parametrize('length',[float('nan'),.8,1.4])
def test_invalid_caps_are_not_silently_used(parent,length):
    with pytest.raises(ValueError):module.fragments(parent,TEMPLATE,length)


def test_multiple_attachment_or_missing_motif_templates_fail_closed(parent):
    for template in ['[*:1]NC[*:2]','CCO']:
        with pytest.raises(ValueError):module.fragments(parent,template)


def test_an_ethyl_or_methyl_precursor_is_not_misidentified_as_a_polymer_benzyl_site():
    for smiles in ['CN1C(=O)c2ccccc2C1=O','CCN1C(=O)c2ccccc2C1=O']:
        mol=Chem.AddHs(Chem.MolFromSmiles(smiles));assert AllChem.EmbedMolecule(mol,randomSeed=32)==0
        with pytest.raises(ValueError,match='primary benzyl'):module.fragments(mol,TEMPLATE)


def test_cap_length_sensitivity_does_not_reposition_the_core(parent):
    short=module.fragments(parent,TEMPLATE,1.07);long=module.fragments(parent,TEMPLATE,1.11)
    for a,b in zip(short,long):
        x=a['mol'].GetConformer().GetPositions();y=b['mol'].GetConformer().GetPositions()
        assert np.array_equal(x[:-1],y[:-1])
        assert np.linalg.norm(x[-1]-y[-1])==pytest.approx(.04)
