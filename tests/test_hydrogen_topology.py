"""Proton transfer must not masquerade as an unchanged exact molecular graph."""
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from redox.qm.uma import geometry_checks
from redox.validation.protonation import hydrogen_topology_status


def write_geometry(tmp_path, smiles, *, reverse=False, force_columns=False):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=14) == 0
    if reverse:
        molecule = Chem.RenumberAtoms(molecule, list(reversed(range(molecule.GetNumAtoms()))))
    xyz = Chem.MolToXYZBlock(molecule)
    if force_columns:
        lines = xyz.splitlines()
        xyz = '\n'.join(lines[:2] + [line + ' 0 0 0' for line in lines[2:]]) + '\n'
    path = tmp_path / 'geometry.xyz'
    path.write_text(xyz)
    return path


@pytest.mark.parametrize(('declared', 'actual'), [
    ('NCCC(=O)O', '[NH3+]CCC(=O)[O-]'),
    ('CC=O', 'C=CO'),
])
def test_proton_transfer_with_unchanged_composition_and_heavy_graph_is_detected(tmp_path, declared, actual):
    path = write_geometry(tmp_path, actual)
    assert not geometry_checks(declared, path, 0, 1)  # Existing heavy-atom gate misses this.
    assert hydrogen_topology_status(declared, path) == 'hydrogen_connectivity_mismatch'


@pytest.mark.parametrize('smiles', ['NCCC(=O)O', '[NH3+]CCC(=O)[O-]', 'CCO'])
def test_unchanged_proton_topology_allows_atom_permutations_and_ase_columns(tmp_path, smiles):
    path = write_geometry(tmp_path, smiles, reverse=True, force_columns=True)
    assert hydrogen_topology_status(smiles, path) == 'ok'


def test_missing_hydrogen_and_missing_file_do_not_pass(tmp_path):
    path = write_geometry(tmp_path, 'CCO')
    lines = path.read_text().splitlines()
    path.write_text(str(int(lines[0]) - 1) + '\n' + '\n'.join(lines[1:-1]) + '\n')
    assert hydrogen_topology_status('CCO', path) == 'composition_mismatch'
    assert hydrogen_topology_status('CCO', tmp_path / 'absent.xyz') == 'unreadable'
