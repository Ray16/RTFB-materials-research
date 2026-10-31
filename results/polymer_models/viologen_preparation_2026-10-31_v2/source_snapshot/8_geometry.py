"""Map declared molecular graphs to stored XYZ atom order and verify stereochemistry.

Heavy cheminformatics dependencies are imported only when a check is requested.
"""
from pathlib import Path

# --- SMILES <-> stored-geometry atom mapping -------------------------------------------------
# WHY THIS EXISTS: build_candidates.py embeds the mol returned by decorate() but records
# Chem.MolToSmiles(mol) -- the CANONICAL SMILES -- in the manifest. Canonicalization reorders
# atoms, so for every GRAFTED candidate the manifest-SMILES atom indices do NOT correspond to
# the atom order in library/*.xyz or calcs/dft/*/opt.xyz. Standalone entries (built from a
# direct SMILES) happen to round-trip to identity, so anything that indexes a geometry by
# SMILES index silently works on the standalones and is scrambled on the grafted molecules --
# which are the actual screening candidates. Always go through this mapper.
#
# Bond ORDERS cannot be perceived reliably from geometry for charged aromatics (RDKit's
# DetermineBonds fails on viologen dications and radical anions), but CONNECTIVITY can, so the
# match is done on the element+connectivity skeleton with charges and bond orders stripped.

def _skeleton(mol):
    """Element + connectivity only: no formal charges, no explicit-H counts, no bond orders."""
    from rdkit import Chem
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0); a.SetNumExplicitHs(0); a.SetNoImplicit(True); a.SetIsAromatic(False)
        a.SetNumRadicalElectrons(0)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    return rw.GetMol()


def smiles_to_xyz_map(smiles: str, xyz_path):
    """Heavy-atom index in `smiles` -> heavy-atom index in the geometry at `xyz_path`.

    Returns a tuple m where m[i] is the xyz index of SMILES heavy atom i, or None if the
    geometry cannot be read / matched. Verified to map all 12 current candidate molecules.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    try:
        ref = Chem.MolFromXYZFile(str(xyz_path))
        if ref is None:
            return None
        rdDetermineBonds.DetermineConnectivity(ref, charge=0)
        ref = _skeleton(Chem.RemoveHs(ref, sanitize=False))
        q = Chem.MolFromSmiles(smiles)
        if q is None:
            return None
        return ref.GetSubstructMatch(_skeleton(Chem.RemoveHs(q)), useChirality=False) or None
    except Exception:
        return None


def geometry_stereo_status(smiles: str, xyz_path):
    """Check declared tetrahedral/E-Z stereochemistry separately from graph mapping.

    Reconstruct the DECLARED bond-order graph using coordinates matched by element
    and connectivity. This avoids attempting to infer charged redox bond orders.
    Extra XYZ force columns are ignored. Unknown/flattened stereochemistry fails
    explicitly; achiral inputs have no stereochemical requirement.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    try:
        expected = Chem.AddHs(Chem.MolFromSmiles(smiles))
        specified = [element for element in Chem.FindPotentialStereo(expected)
                     if element.specified == Chem.StereoSpecified.Specified]
        if not specified:
            return 'not_applicable'
        supported = (Chem.StereoType.Atom_Tetrahedral, Chem.StereoType.Bond_Double)
        if any(element.type not in supported for element in specified):
            return 'unsupported_stereo_type'
        lines = Path(xyz_path).read_text().splitlines()
        plain = '\n'.join(lines[:2] + [' '.join(line.split()[:4]) for line in lines[2:]]) + '\n'
        observed = Chem.MolFromXYZBlock(plain)
        if observed is None or observed.GetNumAtoms() != expected.GetNumAtoms():
            return 'unreadable_or_composition_mismatch'
        rdDetermineBonds.DetermineConnectivity(observed, charge=0)
        mapping = _skeleton(observed).GetSubstructMatch(_skeleton(expected), useChirality=False)
        if len(mapping) != expected.GetNumAtoms():
            return 'connectivity_mismatch'
        rebuilt = Chem.Mol(expected)
        rebuilt.RemoveAllConformers()
        conformer = Chem.Conformer(expected.GetNumAtoms())
        conformer.Set3D(True)
        for atom_index, xyz_index in enumerate(mapping):
            conformer.SetAtomPosition(atom_index, observed.GetConformer().GetAtomPosition(xyz_index))
        rebuilt.AddConformer(conformer)
        Chem.RemoveStereochemistry(rebuilt)
        Chem.AssignStereochemistryFrom3D(rebuilt, replaceExistingTags=True)
        # Compare canonical stereochemical representations, allowing symmetric atom
        # permutations without accepting an inverted stereocenter or alkene.
        expected_smiles = Chem.MolToSmiles(Chem.RemoveHs(expected), isomericSmiles=True)
        actual_smiles = Chem.MolToSmiles(Chem.RemoveHs(rebuilt), isomericSmiles=True)
        return 'ok' if expected_smiles == actual_smiles else 'stereochemistry_mismatch'
    except (ValueError, RuntimeError, OSError, AttributeError):
        return 'stereochemistry_unverifiable'


def redox_core_xyz_indices(smiles: str, xyz_path, tether_smarts="[CH3]c1ccc([CH2])cc1"):
    """xyz heavy-atom indices of the REDOX CORE, i.e. ring atoms outside the Merrifield
    benzyl tether plus the exocyclic atoms double-bonded to them (quinone / imide C=O).

    Acyclic substituents (N-methyl, N-ethyl, ammoniopropyl, methoxy) are EXCLUDED: they are
    conformational degrees of freedom, not redox-active, and letting them into an RMSD is
    what makes a whole-molecule test too blunt to separate a tether rotation from a genuine
    inner-sphere distortion. Returns (core_idx, tether_idx) as sorted lists of xyz indices.
    """
    from rdkit import Chem
    m = smiles_to_xyz_map(smiles, xyz_path)
    if m is None:
        return None, None
    q = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
    teth = set(q.GetSubstructMatch(Chem.MolFromSmarts(tether_smarts)))
    ri = q.GetRingInfo()
    ring = {i for i in range(q.GetNumAtoms()) if ri.NumAtomRings(i) > 0}
    core = ring - teth
    for a in list(core):
        for b in q.GetAtomWithIdx(a).GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE:
                o = b.GetOtherAtomIdx(a)
                if o not in ring:
                    core.add(o)
    return sorted(m[i] for i in core), sorted(m[i] for i in teth)
