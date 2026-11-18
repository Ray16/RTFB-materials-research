"""The gas-phase single point runs once, on the winning seed only, with unchanged results.

Seeds are ranked on E_smd alone, so computing the gas energy for losing seeds was pure
repeated work. These tests pin the contract of the optimized path with a fake DFT backend:
  * the gas single point is requested exactly once, for the lowest-E_smd seed, on that seed's
    in-memory optimized coordinates (not a re-read of the 8-decimal .xyz);
  * result.json carries the same gas fields and dG_solv_eV the combined call produced;
  * a seed whose gas SCF raises is discarded and the next-lowest seed wins, exactly as when
    the gas point ran inside each seed's call.
"""
import json

import numpy as np
import pytest
from ase import Atoms
from ase.io import write as ase_write

from redox.backends import dft
from redox.core.common import HARTREE_EV

E_SMD_EV = [-1.0, -1.2, -1.1]          # seed 1 wins, seed 2 is the runner-up
E_GAS_HA = {0: -0.010, 1: -0.020, 2: -0.030}


@pytest.fixture
def harness(tmp_path, monkeypatch):
    seeds = []
    for i in range(3):
        p = tmp_path / f"conf_{i:02d}.xyz"
        ase_write(str(p), Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74 + 0.01 * i]]))
        seeds.append(p)
    monkeypatch.setattr(dft, "DFT", tmp_path / "dft")
    monkeypatch.setattr(dft, "conformer_seeds", lambda gid, st, n: seeds[:n])
    monkeypatch.setattr(dft, "_uma_mult", lambda gid, st, fallback: fallback)
    calls = []

    def fake(xyz_path, q, mult, do_opt=True, do_gas=True, do_smd=True, do_freq=True,
             opt_out=None, backend="cpu", atoms_override=None, return_atoms=False, **kw):
        if do_smd:                                   # per-seed SMD optimization + single point
            si = int(str(xyz_path).split("_seed_")[1][:2])
            assert do_gas is False, "seeds must not compute the gas single point"
            atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7412345678901 + si]])
            if opt_out is not None:
                ase_write(str(opt_out), atoms)
            e = E_SMD_EV[si]
            out = dict(e_smd_eV=e, e_smd_Ha=e / HARTREE_EV, converged_smd=True,
                       sp_xc="wb97m-v", sp_basis="def2-tzvpd", sp_disp=None, sp_nlc="vv10",
                       solvent="acetonitrile", scf_guesses=["minao"],
                       timing_s=dict(opt=1.0, sp_smd=2.0))
            if return_atoms:
                out["_atoms"] = atoms
            return out
        # gas-only call: must be on the winner's exact in-memory coordinates, at its sp level
        assert do_gas and not do_opt and atoms_override is not None
        assert kw["sp_xc"] == "wb97m-v" and kw["sp_basis"] == kw["sp_basis_anion"] == "def2-tzvpd"
        assert kw["sp_nlc"] == "vv10" and kw["solvent"] == "acetonitrile"
        si = int(round(atoms_override.positions[1, 2] - 0.7412345678901))
        assert atoms_override.positions[1, 2] == pytest.approx(0.7412345678901 + si, abs=0)
        calls.append(si)
        if harness_fail.get(si):
            raise RuntimeError("gas SCF failed")
        return dict(e_gas_Ha=E_GAS_HA[si], e_gas_eV=E_GAS_HA[si] * HARTREE_EV,
                    converged_gas=True, scf_guess_energies_gas=[dict(guess="minao")],
                    gas_occupation_guesses=[], gas_homo_eV=-5.0, gas_s_squared=0.0,
                    timing_s=dict(sp_gas=3.0))

    harness_fail = {}
    monkeypatch.setattr(dft, "dft_smd", fake)
    row = dict(id="h2", state="neu", charge=0, n_e=0, mult_hint=1)
    return tmp_path, row, calls, harness_fail


def _run(tmp_path, row):
    failures = dft.run_batch([row], force=False, do_opt=True, do_freq=False, nconf=3, preopt="none")
    return failures, json.loads((tmp_path / "dft" / "h2" / "neu" / "result.json").read_text())


def test_gas_single_point_runs_once_on_the_winner(harness):
    tmp_path, row, calls, _ = harness
    failures, rec = _run(tmp_path, row)
    assert failures == 0
    assert calls == [1]
    assert rec["winning_seed"] == 1
    assert rec["e_gas_Ha"] == E_GAS_HA[1]
    assert rec["dG_solv_eV"] == (E_SMD_EV[1] / HARTREE_EV - E_GAS_HA[1]) * HARTREE_EV
    assert rec["gas_homo_eV"] == -5.0 and rec["converged_gas"] is True
    assert rec["seed_e_smd_eV"] == E_SMD_EV
    assert rec["timing_s"]["sp_gas_winner"] == 3.0
    assert "_atoms" not in rec


def test_gas_failure_falls_back_to_the_next_lowest_seed(harness):
    tmp_path, row, calls, fail = harness
    fail[1] = True
    failures, rec = _run(tmp_path, row)
    assert failures == 0
    assert calls == [1, 2]
    assert rec["winning_seed"] == 2
    assert rec["e_gas_Ha"] == E_GAS_HA[2]
    assert rec["seed_e_smd_eV"] == [E_SMD_EV[0], None, E_SMD_EV[2]]
    assert np.isclose(rec["seed_spread_eV"], abs(E_SMD_EV[0] - E_SMD_EV[2]))
