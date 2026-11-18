"""Open-shell SCF guess policy: minao + huckel for organics, escalation on competing solutions.

Evidence for the policy is in the comment above redox.backends.dft.OPEN_SHELL_GUESSES. The
positive control is the OROP-170 radical anion in SMD, where minao and atom converge to the
same solution 0.112 eV above the ground state that only huckel finds.
"""
import pytest

from redox.backends import dft, sp
from redox.core.common import HARTREE_EV
from redox.core import protocol

OROP170_SPLIT_HA = 0.112 / HARTREE_EV


def test_guess_policy_by_state_kind():
    assert dft._guess_policy(0, ["C", "H"]) == (("minao",), ())
    assert dft._guess_policy(1, ["C", "H", "N"]) == (("minao", "huckel"), ("atom",))
    assert dft._guess_policy(1, ["Fe", "C", "H"]) == (("minao", "atom", "huckel"), ())
    assert dft._guess_policy(2, ["Cu", "N"]) == (("minao", "atom", "huckel"), ())


def _run_policy(monkeypatch, energies, converged):
    monkeypatch.setattr(dft, "_free_gpu_pool", lambda: None)

    class FakeMF:
        def __init__(self):
            self.init_guess = None

    calls = []

    def fake_kernel(mf, dm0=None):
        calls.append(mf.init_guess)
        return energies[mf.init_guess], converged[mf.init_guess]

    monkeypatch.setattr(dft, "_kernel_robust", fake_kernel)
    out = dft._lowest_scf_policy(FakeMF, ("minao", "huckel"), ("atom",), post=lambda m: m.init_guess)
    return calls, out


def test_agreeing_guesses_do_not_escalate(monkeypatch):
    calls, (e, conv, post, per, spread, escalated) = _run_policy(
        monkeypatch, {"minao": -1.0, "huckel": -1.0 - 1e-8, "atom": -9.0},
        {"minao": True, "huckel": True, "atom": True})
    assert calls == ["minao", "huckel"]            # atom never run
    assert escalated is False and spread < dft.SCF_SPLIT_EV
    assert (e, conv, post) == (-1.0 - 1e-8, True, "huckel")


def test_orop170_like_split_escalates_and_keeps_the_ground_state(monkeypatch):
    # minao and atom share the upper solution; huckel finds the lower one.
    energies = {"minao": -1.0, "atom": -1.0, "huckel": -1.0 - OROP170_SPLIT_HA}
    calls, (e, conv, post, per, spread, escalated) = _run_policy(
        monkeypatch, energies, dict.fromkeys(energies, True))
    assert calls == ["minao", "huckel", "atom"]
    assert escalated is True and spread == pytest.approx(0.112)
    assert (e, conv, post) == (energies["huckel"], True, "huckel")
    assert [x["guess"] for x in per] == ["minao", "huckel", "atom"]


def test_nonconverged_guess_escalates(monkeypatch):
    calls, (e, conv, post, per, spread, escalated) = _run_policy(
        monkeypatch, {"minao": -1.0, "huckel": -5.0, "atom": -1.2},
        {"minao": True, "huckel": False, "atom": True})
    assert calls == ["minao", "huckel", "atom"]
    assert escalated is True and spread is None
    assert (e, conv, post) == (-1.2, True, "atom")  # lowest CONVERGED solution


def test_adoption_requires_huckel_or_metal_occupation_seeds():
    P = protocol.ACTIVE_SP
    base = dict(optimized=True, sp_xc=P["xc"], sp_basis=P["basis"], sp_nlc=P["nlc"],
                sp_disp=P["disp"], solvent=P["solvent"], converged_smd=True,
                converged_gas=True, e_smd_eV=-1.0, e_gas_eV=-1.0, mult=2, charge=0,
                energy_protocol=P, energy_geom_sha1="g", density_fit=True,
                conv_tol=P["conv_tol"])
    ident = dict(charge=0, mult=2, geom_hash="g")
    assert sp._adoptable(dict(base, scf_guesses=["minao", "huckel"]), **ident)
    assert sp._adoptable(dict(base, scf_guesses=["minao", "atom", "huckel"]), **ident)
    assert sp._adoptable(dict(base, scf_guesses=["rm0", "rm1"]), **ident)
    assert not sp._adoptable(dict(base, scf_guesses=["minao", "atom"]), **ident)
