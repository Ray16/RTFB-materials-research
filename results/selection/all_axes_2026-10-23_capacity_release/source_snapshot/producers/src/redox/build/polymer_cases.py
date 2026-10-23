"""Configurable, site-labelled H-capped polymer scenarios and matched controls.

Use a common tagged backbone before decorating it. CIP R/S labels may change
when substituents change; mapped-neighbor tetrahedral parity must not change.
No occupancy probability, physical loading or whole-network claim is implied.
"""

import hashlib
from collections import Counter
from rdkit import Chem

BACKBONE_MAP = 1000


def resolve_motif(spec, feed):
    """Resolve a registered motif or an explicit config fragment without inheritance."""
    if "motif" in spec:
        definition = spec["motif"]
        fragment = definition["fragment"]
        if Chem.MolFromSmiles(fragment) is None:
            raise ValueError("Invalid inline motif fragment")
        from redox.build.grafting import precursor

        molecule, reaction_class = precursor(fragment)
        if reaction_class != "NH_OH_substitution":
            raise ValueError("Inline charged/quaternizing motifs need an ionic model")
        declared = definition.get("precursor_smiles")
        if declared:
            declared_mol = Chem.MolFromSmiles(declared)
            if declared_mol is None or Chem.MolToSmiles(
                declared_mol
            ) != Chem.MolToSmiles(molecule):
                raise ValueError(
                    "Inline precursor does not match the attachment fragment"
                )
        return dict(id=definition["id"], frag=fragment), "inline configuration"
    identifier = spec["motif_id"]
    rows = feed["original39"] + feed["additional_candidates"]
    matches = [row for row in rows if row["id"] == identifier]
    if len(matches) != 1:
        raise ValueError(f"Unknown registered motif: {identifier}")
    return matches[0], feed["snapshot_path"]


def _path(mol, n):
    atoms = {
        a.GetAtomMapNum() - BACKBONE_MAP: a.GetIdx()
        for a in mol.GetAtoms()
        if BACKBONE_MAP < a.GetAtomMapNum() <= BACKBONE_MAP + 2 * n
    }
    if set(atoms) != set(range(1, 2 * n + 1)):
        raise ValueError("Backbone atom labels are incomplete")
    for position, index in atoms.items():
        atom = mol.GetAtomWithIdx(index)
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            raise ValueError("Backbone must be an open carbon chain")
        neighbors = {
            a.GetAtomMapNum() - BACKBONE_MAP
            for a in atom.GetNeighbors()
            if a.GetIdx() in atoms.values()
        }
        if neighbors != {i for i in (position - 1, position + 1) if 1 <= i <= 2 * n}:
            raise ValueError("Unexpected branch or closure in mapped backbone")
    return atoms


def stereo_fingerprint(mapped_smiles, n):
    mol = Chem.MolFromSmiles(mapped_smiles)
    atoms = _path(mol, n)
    result = []
    for position in range(2, 2 * n, 2):
        atom = mol.GetAtomWithIdx(atoms[position])
        tag = atom.GetChiralTag()
        if tag not in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        ):
            raise ValueError("Unspecified backbone stereochemistry")
        order = []
        for neighbor in atom.GetNeighbors():
            mapping = neighbor.GetAtomMapNum()
            order.append(
                0
                if mapping == BACKBONE_MAP + position - 1
                else 1
                if mapping == BACKBONE_MAP + position + 1
                else 2
            )
        if sorted(order) != [0, 1, 2] or atom.GetTotalNumHs() != 1:
            raise ValueError("Invalid tetrahedral backbone center")
        order.append(3)  # Implicit H follows the three stored heavy neighbors.
        inversions = sum(order[i] > order[j] for i in range(4) for j in range(i + 1, 4))
        sign = (1 if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW else -1) * (
            -1 if inversions % 2 else 1
        )
        result.append([position, sign])
    return result


def build_case(fragment, n, reactive_sites, occupied_sites, stereo):
    if type(n) is not int or not 1 <= n <= 100:
        raise ValueError(
            "Specify 1–100 repeat units; computational admission is separate"
        )
    reactive = set(reactive_sites)
    occupied = set(occupied_sites)
    if not reactive <= set(range(1, n + 1)) or not occupied <= reactive:
        raise ValueError("Grafts must occupy declared reactive sites")
    if len(stereo) != n - 1 or any(t not in ("CW", "CCW") for t in stereo):
        raise ValueError("Specify every backbone tetrahedral tag")

    def aryl(i):
        return f"c1ccc(C[*:{i}])cc1"

    scaffold = "".join("CC(" + aryl(i) + ")" for i in range(1, n)) + "CC" + aryl(n)
    mol = Chem.MolFromSmiles(scaffold)
    # The longest acyclic carbon component is the backbone (the benzyl linkers
    # are isolated carbon components before attachment).
    nodes = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()
    }
    start = next(i for i in nodes if mol.GetAtomWithIdx(i).GetTotalNumHs() == 3)
    chain = [start]
    while True:
        next_atoms = [
            a.GetIdx()
            for a in mol.GetAtomWithIdx(chain[-1]).GetNeighbors()
            if a.GetIdx() in nodes and a.GetIdx() not in chain
        ]
        if not next_atoms:
            break
        if len(next_atoms) != 1:
            raise ValueError("Ambiguous initial backbone")
        chain.append(next_atoms[0])
    if len(chain) != 2 * n:
        raise ValueError("Initial backbone length mismatch")
    for i, index in enumerate(chain, 1):
        mol.GetAtomWithIdx(index).SetAtomMapNum(BACKBONE_MAP + i)
    for position, tag in zip(range(2, 2 * n, 2), stereo):
        mol.GetAtomWithIdx(chain[position - 1]).SetChiralTag(
            Chem.ChiralType.CHI_TETRAHEDRAL_CW
            if tag == "CW"
            else Chem.ChiralType.CHI_TETRAHEDRAL_CCW
        )
    rw = Chem.RWMol(mol)
    remove = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        site = atom.GetAtomMapNum()
        if site not in reactive:
            remove.extend([atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()])
        elif site not in occupied:
            target = rw.GetAtomWithIdx(atom.GetIdx())
            target.SetAtomicNum(17)
            target.SetAtomMapNum(0)
            target.SetNoImplicit(True)
    for index in sorted(remove, reverse=True):
        rw.RemoveAtom(index)
    combined = rw.GetMol()
    Chem.SanitizeMol(combined)
    for site in sorted(occupied):
        motif = Chem.MolFromSmiles(fragment)
        dummies = [a for a in motif.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) != 1:
            raise ValueError("One terminal motif attachment required")
        dummies[0].SetAtomMapNum(site)
        combined = Chem.CombineMols(combined, motif)
    product = Chem.molzip(combined) if occupied else combined
    Chem.SanitizeMol(product)
    if any(a.GetAtomicNum() == 0 for a in product.GetAtoms()):
        raise ValueError("Unresolved attachment")
    if len(Chem.GetMolFrags(product)) != 1:
        raise ValueError("Expected one covalent polymer chain")
    backbone = _path(product, n)
    states = []
    for site in range(1, n + 1):
        center = product.GetAtomWithIdx(backbone[2 * site])
        ring_atom = next(a for a in center.GetNeighbors() if a.GetIsAromatic())
        ring = next(
            set(r)
            for r in product.GetRingInfo().AtomRings()
            if ring_atom.GetIdx() in r and len(r) == 6
        )
        links = [
            (i, a.GetIdx())
            for i in ring
            for a in product.GetAtomWithIdx(i).GetNeighbors()
            if a.GetAtomicNum() == 6
            and not a.IsInRing()
            and a.GetIdx() not in backbone.values()
        ]
        if site in reactive:
            if (
                len(links) != 1
                or len(Chem.GetShortestPath(product, ring_atom.GetIdx(), links[0][0]))
                - 1
                != 3
            ):
                raise ValueError("Incorrect para benzyl handle")
        elif links:
            raise ValueError(
                "Plain styrene padding unexpectedly contains a benzyl handle"
            )
        states.append(
            "grafted"
            if site in occupied
            else "benzyl_chloride"
            if site in reactive
            else "styrene"
        )
    mapped = Chem.MolToSmiles(product, isomericSmiles=True)
    fingerprint = stereo_fingerprint(mapped, n)
    plain = Chem.Mol(product)
    for atom in plain.GetAtoms():
        atom.SetAtomMapNum(0)
    smiles = Chem.MolToSmiles(plain, isomericSmiles=True)
    return dict(
        smiles=smiles,
        mapped_smiles=mapped,
        repeat_units=n,
        reactive_sites=sorted(reactive),
        occupied_sites=sorted(occupied),
        site_states=states,
        backbone_stereo_tags=list(stereo),
        backbone_stereo_fingerprint=fingerprint,
        backbone_carbon_count=2 * n,
        backbone_CC_bonds=2 * n - 1,
        chain_ends=2,
        atoms=Chem.AddHs(plain).GetNumAtoms(),
        charge=Chem.GetFormalCharge(plain),
        model_loading_per_repeat=len(occupied) / n,
        model_conversion_per_reactive_site=len(occupied) / len(reactive)
        if reactive
        else None,
        graph_sha256=hashlib.sha256(smiles.encode()).hexdigest(),
    )


def panel(spec, fragment):
    if len({scenario["name"] for scenario in spec["scenarios"]}) != len(
        spec["scenarios"]
    ):
        raise ValueError("Scenario names must be unique")
    n = spec["repeat_units"]
    reactive = spec["reactive_sites"]
    stereo = spec["backbone_stereo_tags"]
    occupied_sets = {()}
    pairs = []
    for scenario in spec["scenarios"]:
        sites = tuple(sorted(scenario["occupied_sites"]))
        if len(sites) != 2 or len(set(sites)) != 2:
            raise ValueError(
                "The first comparison panel requires two distinct occupied sites"
            )
        if any(i in (1, n) for i in sites):
            raise ValueError("Pair sites require at least one terminal padding unit")
        occupied_sets.update([sites, (sites[0],), (sites[1],)])
        pairs.append(
            dict(
                name=scenario["name"],
                sites=list(sites),
                intervening_repeat_units=sites[1] - sites[0] - 1,
                coefficients={
                    "occ_" + "_".join(map(str, sites)): 1,
                    "occ_none": 1,
                    "occ_" + str(sites[0]): -1,
                    "occ_" + str(sites[1]): -1,
                },
            )
        )
    cases = []
    for sites in sorted(occupied_sets, key=lambda s: (len(s), s)):
        case = build_case(fragment, n, reactive, sites, stereo)
        case["key"] = "occ_" + "_".join(map(str, sites)) if sites else "occ_none"
        if case["charge"] != 0:
            raise ValueError(
                "Only neutral controls are admitted; charged panels require "
                "explicit spin/counterion models"
            )
        cases.append(case)
    fingerprints = {str(c["backbone_stereo_fingerprint"]) for c in cases}
    if len(fingerprints) != 1:
        raise ValueError(
            "Controls do not preserve the same mapped backbone stereochemistry"
        )
    lookup = {c["key"]: c for c in cases}
    for cycle in pairs:
        balance = Counter()
        charge = 0
        for key, coefficient in cycle["coefficients"].items():
            mol = Chem.AddHs(Chem.MolFromSmiles(lookup[key]["smiles"]))
            for atom in mol.GetAtoms():
                balance[atom.GetAtomicNum()] += coefficient
            charge += coefficient * Chem.GetFormalCharge(mol)
        if any(balance.values()) or charge:
            raise ValueError("Site-nonadditivity cycle is not atom/charge balanced")
    return cases, pairs
