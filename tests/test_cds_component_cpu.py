"""The real primary legacy CDS function can be evaluated without a GPU import."""
from types import SimpleNamespace
import sys

import numpy as np
from pyscf import gto

from redox.validation.cds_hessian_component import cpu_legacy
from redox.qm.cds_finite import energy_gradient


def test_installed_sulfur_cds_energy_and_gradient_are_cpu_only_and_consistent():
    legacy,sources,descriptors=cpu_legacy()
    assert len(descriptors)==8 and len(sources)==2 and 'cupy' not in sys.modules
    mol=gto.M(atom='S 0 0 0; H 1.3 0 0; H -0.4 1.2 0',basis='def2-svp',verbose=0)
    solvent=SimpleNamespace(mol=mol,solvent='acetonitrile',solvent_descriptors=None)
    coordinates=mol.atom_coords(unit='Bohr');energy,gradient=legacy(solvent)
    def evaluate(point):
        mol.set_geom_(point,unit='Bohr')
        return legacy(solvent)[0]
    numerical=energy_gradient(coordinates,evaluate,step=1e-5)
    assert np.isfinite(energy) and gradient.shape==(3,3)
    assert abs(numerical-gradient).max()<1e-6
    assert 'cupy' not in sys.modules
