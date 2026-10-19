"""Validate motif definitions before creating any files or scheduling calculations."""
import re

from rdkit import Chem

from redox.build.structures import decorate, unassigned_stereo


def validate_group(scaffold, group):
    for field in ("id", "name", "family", "frag", "states"):
        if not group.get(field):
            raise ValueError(f"missing candidate field: {field}")
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", group["id"]):
        raise ValueError("candidate id must be a safe directory name")
    mol = decorate(scaffold, group["frag"])
    if unassigned_stereo(mol):
        raise ValueError(f"{group['id']}: unspecified stereochemistry")
    states = group["states"]
    labels, charges, resting = set(), set(), []
    nuclei = sum(a.GetAtomicNum() for a in Chem.AddHs(mol).GetAtoms())
    for state in states:
        if len(state) != 4:
            raise ValueError("states require (label, charge, multiplicity, n_e)")
        label, charge, mult, ne = state
        if not isinstance(label, str) or not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", label):
            raise ValueError("state label must be a safe directory name")
        if label in labels or charge in charges:
            raise ValueError("duplicate state label or charge")
        if any(type(x) is not int for x in (charge, mult, ne)):
            raise ValueError("charge, multiplicity and n_e must be integers")
        electrons = nuclei - charge
        if mult < 1 or electrons < mult - 1 or (electrons - (mult - 1)) % 2:
            raise ValueError(f"{label}: multiplicity incompatible with electron count")
        labels.add(label); charges.add(charge)
        if ne == 0:
            resting.append(charge)
    if len(resting) != 1 or resting[0] != Chem.GetFormalCharge(mol):
        raise ValueError("exactly one n_e=0 state must match the grafted SMILES charge")
    if any(q - resting[0] != ne for _, q, _, ne in states):
        raise ValueError("n_e must equal charge minus resting charge")
    if len(states) < 2 or max(charges) - min(charges) + 1 != len(charges):
        raise ValueError("declare a contiguous ladder of at least two charge states")
    return mol
