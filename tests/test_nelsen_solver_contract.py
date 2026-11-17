"""The established fresh-gas API retains eV units and refuses stopped solvers."""
from types import SimpleNamespace

import numpy as np
import pytest

from redox.properties import nelsen


@pytest.mark.parametrize('converged,energy',[(False,-1.),(True,float('nan'))])
def test_nelsen_single_points_reject_invalid_solver_values(monkeypatch,converged,energy):
    monkeypatch.setattr(nelsen,'mol',lambda *a:None)
    monkeypatch.setattr(nelsen,'mean_field',lambda *a:SimpleNamespace(kernel=lambda:energy,converged=converged))
    with pytest.raises(RuntimeError):nelsen.energy(None,0,0,'pbe','sto-3g',backend='cpu')


def test_nelsen_energy_keeps_its_public_eV_units(monkeypatch):
    monkeypatch.setattr(nelsen,'mol',lambda *a:None)
    monkeypatch.setattr(nelsen,'mean_field',lambda *a:SimpleNamespace(kernel=lambda:-2.,converged=True))
    assert nelsen.energy(None,0,0,'pbe','sto-3g',backend='cpu')==pytest.approx(-2*27.211386245988)


@pytest.mark.parametrize('converged',[False,True])
def test_nelsen_gas_geometry_uses_strict_driver_and_does_not_rescore_failure(monkeypatch,converged):
    from pyscf.geomopt import geometric_solver
    gradient=SimpleNamespace(grid_response=False,as_scanner=lambda:SimpleNamespace(converged=True))
    monkeypatch.setattr(nelsen,'mol',lambda *a:None)
    monkeypatch.setattr(nelsen,'mean_field',lambda *a:SimpleNamespace(nuc_grad_method=lambda:gradient))
    molecule=SimpleNamespace(atom_coords=lambda:np.zeros((2,3)))
    def kernel(scanner,**kwargs):
        assert gradient.grid_response and kwargs['maxsteps']==4
        return converged,molecule
    monkeypatch.setattr(geometric_solver,'kernel',kernel)
    atoms=object();monkeypatch.setattr(nelsen,'atoms_of',lambda m:atoms)
    monkeypatch.setattr(nelsen,'energy',lambda *a,**k: -12. if converged else pytest.fail('stopped geometry rescored'))
    if converged:
        assert nelsen.gas_opt(None,0,0,'pbe','sto-3g',backend='cpu',maxsteps=4)==(atoms,-12.)
    else:
        with pytest.raises(RuntimeError):nelsen.gas_opt(None,0,0,'pbe','sto-3g',backend='cpu',maxsteps=4)
