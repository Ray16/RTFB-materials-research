"""Ferrocene is written without metal bonds in SMILES; its eta-5 Fe-C contacts must not fail the identity gate."""
from redox.structures.validation import ferrocene_xyz
from redox.backends.uma import geometry_checks

FERROCENE = "[Fe].c1ccc[cH-]1.c1ccc[cH-]1"


def test_ferrocene_sandwich_passes_identity_gate(tmp_path):
    path = tmp_path / "fc.xyz"
    path.write_text(ferrocene_xyz())
    assert geometry_checks(FERROCENE, path, 0, 1) == {}


def test_dissociated_metal_is_still_rejected(tmp_path):
    lines = ferrocene_xyz().splitlines()
    moved = [" ".join([t[0], t[1], t[2], str(float(t[3]) + 6.0)]) if t and t[0] == "Fe" else line
             for line, t in ((line, line.split()) for line in lines[2:])]
    path = tmp_path / "fc_far.xyz"
    path.write_text("\n".join(lines[:2] + moved) + "\n")
    assert geometry_checks(FERROCENE, path, 0, 1).get("uncoordinated_metal") is True
