"""Exact precursor/product and heavy-atom attachment mapping for local grafts."""

from rdkit import Chem

from redox.build.contracts import validate_group
from redox.build.grafting import grafting_reaction, precursor
from redox.workflows.chunked import content_hash


def canonical(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("One valid connected molecular graph required")
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def exact_attachment(task):
    """Derive the new bond from the same two dummies used by the real builder.

    Maps identify heavy atoms, including the chloride leaving atom. Proton/base
    transfer is represented by the separate atom/charge-balanced reaction; this
    is not a fully proton-mapped reaction or a whole-polymer environment.
    """
    identity, payload = task["identity"], task["payload"]
    fragment = payload["attachment_fragment"]
    if fragment != identity["attachment_fragment"]:
        raise ValueError("Payload fragment differs from the exact design identity")
    group = dict(id="checked_attachment", name="local_graft", family=payload["family"],
                 frag=fragment, states=payload["states"])
    expected = validate_group(identity["scaffold_smiles"], group)
    product = canonical(payload["product_smiles"])
    if product != canonical(identity["product_smiles"]) or product != canonical(Chem.MolToSmiles(expected)):
        raise ValueError("Declared product differs from the constructed graft")
    reaction = grafting_reaction(fragment, product)
    precursor_smiles = canonical(reaction["precursor_smiles"])
    if any(canonical(value) != precursor_smiles for value in
           (payload["precursor_smiles"], identity["precursor_smiles"])):
        raise ValueError("Declared precursor differs from the attachment route")
    # Canonicalize traversal before mapping so alternate SMILES of one graph
    # produce the same heavy-atom attachment identity.
    scaffold = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(identity["scaffold_smiles"])))
    motif = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(fragment)))
    maps = []
    offset = 0
    for part in (scaffold, motif):
        dummies = [a for a in part.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) != 1 or dummies[0].GetAtomMapNum() != 1 or dummies[0].GetDegree() != 1:
            raise ValueError("One terminal [*:1] required per attachment partner")
        for atom in part.GetAtoms():
            if atom.GetAtomicNum():
                atom.SetAtomMapNum(1001 + offset + atom.GetIdx())
        maps.append(dummies[0].GetNeighbors()[0].GetAtomMapNum())
        offset += part.GetNumAtoms()
    mapped = Chem.molzip(Chem.CombineMols(scaffold, motif))
    Chem.SanitizeMol(mapped)
    indices = {atom.GetAtomMapNum(): atom.GetIdx() for atom in mapped.GetAtoms()}
    bond = mapped.GetBondBetweenAtoms(indices[maps[0]], indices[maps[1]])
    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
        raise ValueError("Attachment is not the intended single bond")
    mapped_precursor, kind = precursor(Chem.MolToSmiles(motif))
    site = Chem.RWMol(scaffold)
    leaving = next(a for a in site.GetAtoms() if a.GetAtomicNum() == 0)
    leaving.SetAtomicNum(17)
    leaving.SetAtomMapNum(9001)
    Chem.SanitizeMol(site)
    if canonical(Chem.MolToSmiles(mapped)) != product:
        raise ValueError("Heavy-atom mapping changed the product")
    if canonical(Chem.MolToSmiles(mapped_precursor)) != precursor_smiles:
        raise ValueError("Heavy-atom mapping changed the precursor")
    if canonical(Chem.MolToSmiles(site)) != canonical(reaction["site_smiles"]):
        raise ValueError("Support-handle model differs from the balanced reaction")
    binding = dict(product_smiles=product, precursor_smiles=precursor_smiles,
                   site_smiles=canonical(reaction["site_smiles"]), reaction_class=kind,
                   mapped_product_smiles=Chem.MolToSmiles(mapped, isomericSmiles=True),
                   mapped_precursor_smiles=Chem.MolToSmiles(mapped_precursor, isomericSmiles=True),
                   mapped_site_smiles=Chem.MolToSmiles(site, isomericSmiles=True),
                   new_bond_atom_maps=maps, leaving_atom_map=9001,
                   heavy_atom_mapping=True, model_scope="4-methylbenzyl local site on a preformed polymer")
    return dict(exact_design_key=content_hash(binding), route_binding=binding, balanced_reaction=reaction)
