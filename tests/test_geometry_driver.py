"""The optimizer must preserve gradient settings and report actual convergence."""
from types import SimpleNamespace

import pytest

from redox.qm import dft


def test_geometry_driver_exposes_kernel_convergence():
    from pyscf.geomopt.geometric_solver import kernel
    assert dft._geomopt_fn("cpu") is kernel
    assert dft._geomopt_fn("gpu") is kernel


@pytest.mark.parametrize("geometry_ok,scf_ok,strict", [
    (True, True, True), (False, True, True),
    (False, True, False), (True, False, True), (True, False, False),
])
def test_geometry_convergence_and_full_grid_response(
    tmp_path, monkeypatch, geometry_ok, scf_ok, strict
):
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("2\n\nH 0 0 0\nH 0 0 .74\n")
    target = tmp_path / "opt.xyz"
    target.write_text("preserved previous geometry\n")
    original = target.read_bytes()
    scanner = SimpleNamespace(converged=scf_ok)
    gradient = SimpleNamespace(grid_response=False, as_scanner=lambda: scanner)
    seen = {}

    def build(mol, *args):
        seen["mol"] = mol
        return SimpleNamespace(nuc_grad_method=lambda: gradient)

    def optimize(actual_scanner, **kwargs):
        assert actual_scanner is scanner
        assert gradient.grid_response is True
        assert kwargs["assert_convergence"] is strict
        assert kwargs["maxsteps"] == 3
        return geometry_ok, seen["mol"]

    monkeypatch.setattr(dft, "_dft_module", lambda backend: None)
    monkeypatch.setattr(dft, "_build_smd_mf", build)
    monkeypatch.setattr(dft, "_geomopt_fn", lambda backend: optimize)
    from redox.core import protocol
    monkeypatch.setattr(protocol, "software_versions", lambda: {})
    kwargs = dict(do_gas=False, do_smd=False, do_freq=False,
                  opt_out=target, max_opt_steps=3, assert_convergence=strict)
    if strict and not (geometry_ok and scf_ok):
        with pytest.raises(RuntimeError, match="Geometry optimization did not converge"):
            dft.dft_smd(xyz, 0, 1, **kwargs)
        assert target.read_bytes() == original
    else:
        result = dft.dft_smd(xyz, 0, 1, **kwargs)
        assert result["optimized"] is (geometry_ok and scf_ok)
        assert result["geometry_optimization_converged"] is geometry_ok
        assert result["optimization_SCF_converged"] is scf_ok
        assert result["opt_gradient_grid_response"] is True
        assert target.read_bytes() != original


def test_single_point_does_not_claim_optimization(tmp_path, monkeypatch):
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("2\n\nH 0 0 0\nH 0 0 .74\n")
    monkeypatch.setattr(dft, "_dft_module", lambda backend: None)
    monkeypatch.setattr(dft, "_geomopt_fn", lambda backend: pytest.fail("single point optimized"))
    from redox.core import protocol
    monkeypatch.setattr(protocol, "software_versions", lambda: {})
    result = dft.dft_smd(xyz, 0, 1, do_opt=False, do_freq=False, do_smd=False, do_gas=False)
    assert result["optimized"] is False
    assert "geometry_optimization_converged" not in result
