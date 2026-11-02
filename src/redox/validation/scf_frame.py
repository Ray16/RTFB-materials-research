"""Validate the molecular frame actually serialized in an SCF checkpoint."""
import numpy as np
from pyscf.lib import chkfile


def checkpoint_frame(path, *, charge, mult, basis, coordinates_Bohr, symbols):
    molecule=chkfile.load_mol(str(path))
    actual=np.asarray(molecule.atom_coords(),dtype=float);expected=np.asarray(coordinates_Bohr,dtype=float)
    elements=[molecule.atom_pure_symbol(i) for i in range(molecule.natm)]
    if (actual.shape!=expected.shape or actual.shape!=(len(symbols),3)
            or not np.isfinite(actual).all() or not np.isfinite(expected).all()
            or elements!=list(symbols) or (molecule.charge,molecule.spin)!=(charge,mult-1)
            or not isinstance(molecule.basis,str) or molecule.basis.lower()!=basis.lower()):
        raise ValueError('Actual checkpoint molecular geometry/state/basis differs')
    error=float(np.max(abs(actual-expected)))
    if error>1e-10:raise ValueError('Actual checkpoint coordinate frame differs')
    return dict(atoms=molecule.natm,charge=molecule.charge,multiplicity=molecule.spin+1,
                basis=molecule.basis,maximum_coordinate_error_Bohr=error)
