"""Validate the molecular frame actually serialized in an SCF checkpoint."""
import numpy as np
from pyscf.lib import chkfile
import math
from pathlib import Path
import hashlib

SOURCE_BYTES=Path(__file__).read_bytes()


def freeze_reviewer(case):
    """Keep the exact CPU reviewer available beside its immutable receipt."""
    from redox.core.hashing import sha256_file
    digest=hashlib.sha256(SOURCE_BYTES).hexdigest()
    frozen=Path(case)/'review_sources'/digest/'scf_frame.py'
    if not frozen.exists():
        frozen.parent.mkdir(parents=True,exist_ok=True);frozen.write_bytes(SOURCE_BYTES)
    if sha256_file(frozen)!=digest:raise ValueError('Frozen checkpoint reviewer differs')
    return frozen


def check_spin_contamination(value, limit=.1):
    if (type(value) not in (int,float) or not math.isfinite(value)
            or not 0<=abs(value)<=limit):
        raise ValueError('Finite spin contamination within the declared limit required')
    return value


def checkpoint_frame(path, *, charge, mult, basis, coordinates_Bohr, symbols, maximum_error_Bohr=1e-10):
    if (type(charge) is not int or type(mult) is not int or mult<1
            or not isinstance(basis,str) or not basis or not symbols):
        raise ValueError('Explicit checkpoint charge, multiplicity, basis and elements required')
    if type(maximum_error_Bohr) not in (int,float) or not math.isfinite(maximum_error_Bohr) or not 0<maximum_error_Bohr<=1e-10:
        raise ValueError('An explicit checkpoint coordinate tolerance no looser than 1e-10 Bohr required')
    molecule=chkfile.load_mol(str(path))
    actual=np.asarray(molecule.atom_coords(),dtype=float);expected=np.asarray(coordinates_Bohr,dtype=float)
    elements=[molecule.atom_pure_symbol(i) for i in range(molecule.natm)]
    if (actual.shape!=expected.shape or actual.shape!=(len(symbols),3)
            or not np.isfinite(actual).all() or not np.isfinite(expected).all()
            or elements!=list(symbols) or (molecule.charge,molecule.spin)!=(charge,mult-1)
            or not isinstance(molecule.basis,str) or molecule.basis.lower()!=basis.lower()):
        raise ValueError('Actual checkpoint molecular geometry/state/basis differs')
    error=float(np.max(abs(actual-expected)))
    if error>maximum_error_Bohr:raise ValueError('Actual checkpoint coordinate frame differs')
    return dict(atoms=molecule.natm,charge=molecule.charge,multiplicity=molecule.spin+1,
                basis=molecule.basis,maximum_coordinate_error_Bohr=error)


def checkpoint_spin(path):
    """Reconstruct global S² from saved occupied orbitals and the AO metric.

    No SCF or GPU operation is performed. Integer RKS/UKS occupations, electron
    counts and occupied-subspace orthonormality are checked independently.
    """
    import h5py
    molecule=chkfile.load_mol(str(path))
    with h5py.File(path,'r') as stream:
        coefficients=np.asarray(stream['scf/mo_coeff'][()]);occupations=np.asarray(stream['scf/mo_occ'][()])
    if not np.isfinite(coefficients).all() or not np.isfinite(occupations).all():
        raise ValueError('Finite actual SCF orbitals/occupations required')
    overlap=molecule.intor_symmetric('int1e_ovlp');na,nb=molecule.nelec
    if coefficients.ndim==2 and occupations.ndim==1:
        if molecule.spin!=0 or coefficients.shape!=(molecule.nao_nr(),len(occupations)):
            raise ValueError('Closed-shell restricted checkpoint required')
        if not np.all((abs(occupations)<1e-8)|(abs(occupations-2)<1e-8)) or abs(sum(occupations)-molecule.nelectron)>1e-7:
            raise ValueError('Restricted occupations/electron count differ')
        ca=cb=coefficients[:,occupations>1]
    elif coefficients.ndim==3 and occupations.ndim==2:
        if coefficients.shape!=(2,molecule.nao_nr(),occupations.shape[1]) or occupations.shape[0]!=2:
            raise ValueError('Explicit unrestricted orbital dimensions required')
        if (not np.all((abs(occupations)<1e-8)|(abs(occupations-1)<1e-8))
                or max(abs(occupations.sum(axis=1)-np.asarray([na,nb])))>1e-7):
            raise ValueError('Unrestricted occupations/electron counts differ')
        ca=coefficients[0][:,occupations[0]>.5];cb=coefficients[1][:,occupations[1]>.5]
    else:raise ValueError('Supported restricted/unrestricted checkpoint required')
    errors=[float(np.max(abs(c.conj().T@overlap@c-np.eye(c.shape[1])))) if c.shape[1] else 0. for c in (ca,cb)]
    if max(errors)>1e-7:raise ValueError('Actual occupied SCF orbitals are not orthonormal')
    target=abs(na-nb)/2;expected=target*(target+1)
    cross=ca.conj().T@overlap@cb
    squared=float(expected+min(na,nb)-np.sum(abs(cross)**2))
    contamination=check_spin_contamination(squared-expected)
    return dict(s_squared=squared,expected_s_squared=expected,spin_contamination=contamination,
        alpha_electrons=na,beta_electrons=nb,maximum_occupied_orthonormality_error=max(errors),
        wavefunction_spin_reconstructed=True)
