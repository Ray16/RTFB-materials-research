"""Partial SCF recovery requires exact state and HDF continuity, not an exit flag."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from pyscf import gto
from pyscf.lib import chkfile


@pytest.fixture
def recovery(tmp_path):
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("gradient_recovery", root / "scripts/validation/qualify_polymer_gradient_workspace.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="def2-svp", verbose=0)
    raw = tmp_path / "raw"
    (tmp_path / "plan.json").write_text(json.dumps(dict(raw_root=str(raw))))
    reference = dict(result=dict(e_smd_Ha=-1.1, e_gas_Ha=-1.0))
    for index, energy, solvent in [(0, -1.1, "acetonitrile"), (1, -1.0, None)]:
        folder = raw / "pilot_execution/SCF_checkpoints" / f"{index:03d}"
        folder.mkdir(parents=True)
        hdf = folder / "SCF.h5"
        chkfile.dump(str(hdf), "mol", mol.dumps())
        chkfile.dump(str(hdf), "scf/e_tot", energy)
        identity = dict(charge=0, mult=1, atoms=list(mol.elements), basis=repr(mol.basis),
                        coordinates_Bohr=mol.atom_coords().tolist(), solvent=solvent)
        data = dict(SCF_converged=True, qualified_result=False, identity=identity,
                    checkpoint_sha256=module.api.digest(hdf), E_Ha=energy)
        (folder / "progress.json").write_text(json.dumps(data))
    return module, tmp_path, raw, reference


def test_saved_points_remain_unqualified_after_continuity_check(recovery):
    module, parent, _, reference = recovery
    rows = module.checkpoint_continuity(parent, reference)
    assert len(rows) == 2
    assert all(r["difference_from_previous_eV"] == 0 for r in rows)
    assert not any(r["optimized_state_qualified"] for r in rows)


@pytest.mark.parametrize("field,value", [("mult", 3), ("solvent", "water"),
                                         ("coordinates_Bohr", [[0, 0, 0], [0, 0, 1]])])
def test_same_energy_cannot_hide_changed_state_or_geometry(recovery, field, value):
    module, parent, raw, reference = recovery
    p = raw / "pilot_execution/SCF_checkpoints/000/progress.json"
    d = json.loads(p.read_text())
    d["identity"][field] = value
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError):
        module.checkpoint_continuity(parent, reference)


def test_changed_hdf_is_rejected_before_using_energy(recovery):
    module, parent, raw, reference = recovery
    hdf = raw / "pilot_execution/SCF_checkpoints/000/SCF.h5"
    chkfile.dump(str(hdf), "scf/e_tot", -1.11)
    with pytest.raises(ValueError, match="checksum"):
        module.checkpoint_continuity(parent, reference)


def test_cross_hardware_energy_discrepancy_is_not_silently_accepted(recovery):
    module, parent, _, reference = recovery
    reference["result"]["e_smd_Ha"] += 1e-5
    with pytest.raises(ValueError, match="continuity"):
        module.checkpoint_continuity(parent, reference)
