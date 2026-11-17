"""Chirality-preserving heavy-atom comparisons of declared XYZ structures."""
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign

from redox.build.geometry import _skeleton, geometry_stereo_status
from redox.validation.protonation import hydrogen_topology_status


def alignment_molecule(smiles, xyz):
    if (hydrogen_topology_status(smiles, xyz) != 'ok'
            or geometry_stereo_status(smiles, xyz) not in ('ok', 'not_applicable')):
        raise ValueError('Conformer graph, hydrogen attachment or stereochemistry differs')
    expected = Chem.AddHs(Chem.MolFromSmiles(smiles))
    lines = Path(xyz).read_text().splitlines()
    plain = '\n'.join(lines[:2]+[' '.join(line.split()[:4]) for line in lines[2:]])+'\n'
    observed = Chem.MolFromXYZBlock(plain)
    rdDetermineBonds.DetermineConnectivity(observed, charge=Chem.GetFormalCharge(expected))
    mapping = _skeleton(observed).GetSubstructMatch(_skeleton(expected))
    if len(mapping) != expected.GetNumAtoms():
        raise ValueError('Incomplete full-atom geometry mapping')
    conformer = Chem.Conformer(expected.GetNumAtoms())
    conformer.Set3D(True)
    for i, j in enumerate(mapping):
        conformer.SetAtomPosition(i, observed.GetConformer().GetAtomPosition(j))
    expected.AddConformer(conformer)
    return Chem.RemoveHs(expected)


def symmetry_rmsd_matrix(molecules, maximum_matches=10000):
    if not molecules or type(maximum_matches) is not int or maximum_matches < 1:
        raise ValueError('Nonempty geometry set and positive symmetry-enumeration limit required')
    if len({Chem.MolToSmiles(m, isomericSmiles=True) for m in molecules}) != 1:
        raise ValueError('One exact stereochemical graph per geometry comparison required')
    for molecule in molecules:
        if molecule.GetNumConformers() != 1 or not np.isfinite(molecule.GetConformer().GetPositions()).all():
            raise ValueError('One finite geometry per molecule required')
        if any(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms()):
            raise ValueError('Heavy-atom geometry comparisons require explicit hydrogen removal')
    matches = molecules[0].GetSubstructMatches(molecules[0], uniquify=False, useChirality=True,
                                               maxMatches=maximum_matches+1)
    if not matches or len(matches) > maximum_matches:
        raise ValueError('Symmetry enumeration exceeds its declared limit; no truncated comparison allowed')
    matrix = np.zeros((len(molecules), len(molecules)))
    for i, first in enumerate(molecules):
        for j, second in enumerate(molecules[:i]):
            pairs = second.GetSubstructMatches(first, uniquify=False, useChirality=True,
                                               maxMatches=maximum_matches+1)
            if not pairs or len(pairs) > maximum_matches:
                raise ValueError('Symmetry enumeration exceeds its limit; no truncated comparison allowed')
            maps = [list(enumerate(match)) for match in pairs]
            # Explicit graph automorphisms preserve configured stereocenters;
            # proper rotations do not identify reflected conformational wells.
            value = rdMolAlign.GetBestAlignmentTransform(first, second, map=maps,
                                                         reflect=False, numThreads=1)[0]
            if not np.isfinite(value):
                raise ValueError('Nonfinite symmetry-aware RMSD')
            matrix[i, j] = matrix[j, i] = value
    return matrix, len(matches)
