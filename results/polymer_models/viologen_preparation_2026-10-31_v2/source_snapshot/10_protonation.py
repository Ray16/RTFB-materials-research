"""Check explicit hydrogen attachments without guessing charges or bond orders.

An element count and heavy-atom graph alone cannot distinguish a declared acid
from a proton-transferred zwitterion. This check compares the full connectivity
graphs, allowing arbitrary atom order and symmetry-equivalent proton sites.
"""
from pathlib import Path


def hydrogen_topology_status(smiles: str, xyz_path) -> str:
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from redox.build.geometry import _skeleton

    try:
        lines = Path(xyz_path).read_text().splitlines()
        # ASE output may append force columns after the first four XYZ columns.
        xyz = '\n'.join(lines[:2] + [' '.join(line.split()[:4]) for line in lines[2:]]) + '\n'
        observed = Chem.MolFromXYZBlock(xyz)
        declared = Chem.MolFromSmiles(smiles)
        if observed is None or declared is None:
            return 'unreadable'
        declared = Chem.AddHs(declared)
        if sorted(a.GetAtomicNum() for a in observed.GetAtoms()) != sorted(a.GetAtomicNum() for a in declared.GetAtoms()):
            return 'composition_mismatch'
        rdDetermineBonds.DetermineConnectivity(observed, charge=0)
        actual, expected = _skeleton(observed), _skeleton(declared)
        if actual.HasSubstructMatch(expected) and expected.HasSubstructMatch(actual):
            return 'ok'
        heavy_actual = _skeleton(Chem.RemoveHs(observed, sanitize=False))
        heavy_expected = _skeleton(Chem.RemoveHs(declared))
        if heavy_actual.HasSubstructMatch(heavy_expected) and heavy_expected.HasSubstructMatch(heavy_actual):
            return 'hydrogen_connectivity_mismatch'
        return 'heavy_connectivity_mismatch'
    except (OSError, ValueError, RuntimeError, IndexError):
        return 'unreadable'
