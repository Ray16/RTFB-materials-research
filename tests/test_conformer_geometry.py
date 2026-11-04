from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from redox.properties.conformer_geometry import alignment_molecule, symmetry_rmsd_matrix


def structure(smiles):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=163) == 0
    AllChem.MMFFOptimizeMolecule(molecule)
    return molecule


def save(molecule, path):
    path.write_text(Chem.MolToXYZBlock(molecule))
    return path


def test_arbitrary_xyz_order_and_proper_rotation_preserve_geometry(tmp_path):
    smiles = 'C[C@H](O)CC'
    first = structure(smiles)
    second = Chem.RenumberAtoms(first, list(reversed(range(first.GetNumAtoms()))))
    coords = second.GetConformer().GetPositions()
    for i, value in enumerate(coords):
        second.GetConformer().SetAtomPosition(i, [value[1]+3, -value[0]-2, value[2]+1])
    molecules = [alignment_molecule(smiles, save(m, tmp_path/f'{i}.xyz')) for i, m in enumerate((first, second))]
    before = molecules[0].GetConformer().GetPositions().copy()
    matrix, count = symmetry_rmsd_matrix(molecules)
    assert matrix[0, 1] < 1e-5
    assert np.array_equal(before, molecules[0].GetConformer().GetPositions())
    assert count >= 1


def test_reflected_configured_stereocenter_is_rejected(tmp_path):
    smiles = 'C[C@H](O)CC'
    molecule = structure(smiles)
    for i, xyz in enumerate(molecule.GetConformer().GetPositions()):
        molecule.GetConformer().SetAtomPosition(i, [-xyz[0], xyz[1], xyz[2]])
    with pytest.raises(ValueError, match='stereochemistry'):
        alignment_molecule(smiles, save(molecule, tmp_path/'mirror.xyz'))


def test_unconfigured_gauche_mirrors_are_not_merged_by_reflection(tmp_path):
    first = structure('CCCC')
    rdMolTransforms.SetDihedralDeg(first.GetConformer(), 0, 1, 2, 3, 60)
    second = Chem.Mol(first)
    for i, xyz in enumerate(second.GetConformer().GetPositions()):
        second.GetConformer().SetAtomPosition(i, [-xyz[0], xyz[1], xyz[2]])
    molecules = [alignment_molecule('CCCC', save(m, tmp_path/f'{i}.xyz')) for i, m in enumerate((first, second))]
    matrix, _ = symmetry_rmsd_matrix(molecules)
    assert matrix[0, 1] > .05


def test_symmetry_limit_cannot_silently_truncate_matches(tmp_path):
    molecule = alignment_molecule('CC', save(structure('CC'), tmp_path/'ethane.xyz'))
    with pytest.raises(ValueError, match='truncated'):
        symmetry_rmsd_matrix([molecule], maximum_matches=1)


def test_pairwise_atom_mapping_handles_different_molecule_index_orders():
    first = Chem.RemoveHs(structure('CCO'))
    second = Chem.RenumberAtoms(first, [2, 1, 0])
    matrix, _ = symmetry_rmsd_matrix([first, second])
    assert matrix[0, 1] < 1e-6
