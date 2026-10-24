"""A surface diagnostic must distinguish a CDS kink from other force errors."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cds_surface_diagnostic",
    ROOT / "scripts/validation/diagnose_cds_surface.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.mark.parametrize("unexplained_error", [0., 2e-5])
def test_component_isolation_does_not_hide_a_non_CDS_error(unexplained_error):
    # A kink at +0.00075 is crossed only by the coarsest symmetric probe.
    cds = lambda x: 0.003*x + 0.00015*max(x-0.00075, 0.)
    quantum = lambda x: 0.12*x + 0.002*x**3
    derivatives = {str(h): (cds(h)-cds(-h))/(2*h) for h in m.STEPS}
    h = 0.001
    numeric = (quantum(h)+cds(h)-quantum(-h)-cds(-h))/(2*h)
    check = dict(atom=0, axis=0, step_Bohr=h, numerical=numeric,
                 analytic=0.123-unexplained_error)
    result = m.component_diagnostic(check, derivatives)
    assert result["CDS_coarse_minus_fine_Ha_Bohr"] == pytest.approx(1.875e-5)
    assert result["remaining_difference_Ha_Bohr"] == pytest.approx(unexplained_error+2e-9)


@pytest.mark.parametrize("converged,extra", [(False, 0.), (True, 1e-5)])
def test_displaced_records_require_convergence_and_exact_other_coordinates(converged, extra):
    original = np.zeros((2, 3))
    point = original.copy()
    point[0, 1] = 0.0005
    point[1, 2] = extra
    with pytest.raises(ValueError, match="named coordinate"):
        m.checked_displacement(original,
            dict(coordinates_Bohr=point.tolist(), SCF_converged=converged), 0, 1, 0.0005, 1)


def test_CDS_wrapper_preserves_nuclear_coordinates_and_units():
    xyz = np.array([[1., 2., 3.], [4., 5., 6.]])
    mol = m.CoordinateMolecule([6, 1], xyz)
    assert mol.natm == 2
    np.testing.assert_array_equal(mol.atom_charges(), [6, 1])
    np.testing.assert_array_equal(mol.atom_coords("B"), xyz)
    with pytest.raises(ValueError, match="unit"):
        mol.atom_coords("Angstrom")
