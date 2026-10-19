"""Two imide units do not establish the measured PMDI/NDI redox-core scope."""
import pytest
from rdkit import Chem

from redox.core.common import load_config
from redox.screening.discovery import families, identify


@pytest.mark.parametrize("core", ["pmdi", "ndi"])
@pytest.mark.parametrize("anchor,cap", [("C", "CC"), ("[H]", "C"), ("[H]", "[H]")])
def test_measured_imide_core_recognition_accepts_substituted_and_nh_cores(core, anchor, cap):
    template = load_config("family_search").CORES[core]["template"]
    molecule = Chem.MolFromSmiles(template.format(anchor=anchor, cap=cap))
    assert molecule is not None
    assert "bis-imide" in families(molecule)
    assert "multiple-imide-unvalidated" not in families(molecule)


@pytest.mark.parametrize("smiles", [
    "O=C1c2ccccc2C(=O)N1N1C(=O)c2ccccc2C1=O",
    "O=C1c2ccccc2C(=O)N1CCN1C(=O)c2ccccc2C1=O",
])
def test_linked_phthalimides_retain_unknown_two_wave_behavior(smiles):
    found = families(Chem.MolFromSmiles(smiles))
    assert "multiple-imide-unvalidated" in found
    assert "bis-imide" not in found


def test_unvalidated_multiple_imide_graft_keeps_electron_count_and_capacity_blank():
    smiles = "O=C(O)c1cccc2c1C(=O)N(N1C(=O)c3ccccc3C1=O)C2=O"
    rows = identify([dict(source_id="external_NN_phthalimide", smiles=smiles)], 96485)
    assert rows
    for row in rows:
        assert row["status"] == "multiple imide groups; redox core requires validation"
        assert row["multi_e"] is False
        assert row["n_e_nominal"] is None
        assert row["capacity_nominal_mAh_g"] is None
        assert "unestablished" in row["evidence"]
        assert "0.61-0.64" not in row["evidence"]


def test_registered_nominal_capacity_uses_declared_ladder_not_family_evidence():
    import importlib.util
    from redox.core.common import ROOT
    spec = importlib.util.spec_from_file_location("identify_nominal_test", ROOT / "scripts/mining/identify_candidates.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    group = dict(family="phenothiazine", states=[("neu", 0, 1, 0), ("ox1", 1, 2, 1), ("ox2", 2, 1, 2)])
    assert module.declared_electrons(group) == 2
    assert module.declared_electrons(dict(states=[("neu", 0, 1, 0), ("ox2", 2, 1, 2)])) == 0
