"""Vibrational file contracts shared by independent Hessian validation tools."""
from pathlib import Path
import numpy as np
import re
from redox.core.common import HARTREE_EV


def project_rigid_modes(matrix, coordinates, masses):
    """Remove rigid translation/rotation in mass-weighted Cartesian space.

    Coordinates may use any one consistent length unit. The returned Hessian
    has the input's Cartesian units; the second value is the rigid-space rank
    (six for a nonlinear molecule, five for a linear molecule). Internal
    negative curvature is retained. This operation does not establish that a
    geometry is stationary or that the input Hessian is physically qualified.
    """
    hessian = np.asarray(matrix, dtype=float)
    xyz = np.asarray(coordinates, dtype=float)
    mass = np.asarray(masses, dtype=float)
    if (xyz.ndim != 2 or xyz.shape[1:] != (3,) or len(xyz) == 0
            or mass.shape != (len(xyz),) or hessian.shape != (3 * len(xyz),) * 2
            or not all(np.isfinite(value).all() for value in (hessian, xyz, mass))
            or np.any(mass <= 0)):
        raise ValueError('finite Cartesian Hessian/coordinates and positive atomic masses required')
    centered = xyz - np.average(xyz, axis=0, weights=mass)
    roots = np.sqrt(mass)
    rigid = []
    for axis in np.eye(3):
        rigid.append((roots[:, None] * axis).ravel())
        rigid.append((roots[:, None] * np.cross(axis, centered)).ravel())
    left, singular, _ = np.linalg.svd(np.column_stack(rigid), full_matrices=False)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-10))
    basis = left[:, :rank]
    projector = np.eye(3 * len(xyz)) - basis @ basis.T
    repeated = np.repeat(roots, 3)
    weighted = (hessian + hessian.T) / (2 * np.outer(repeated, repeated))
    projected = (projector @ weighted @ projector) * np.outer(repeated, repeated)
    return (projected + projected.T) / 2, rank


def write_xtb_hessian(path,matrix):
    """xTB6.7.1 Turbomole reader: columns, chunks of5, fixed(5x,5f15.10).

    The numeric matrix remains in atomic units (Hartree/Bohr²), not mass weighted.
    Source: grimme-lab/xtb v6.7.1 src/io/reader/turbomole.f90.
    """
    array=np.asarray(matrix)
    if array.ndim!=2 or array.shape[0]!=array.shape[1] or not np.isfinite(array).all():
        raise ValueError('a finite square Cartesian Hessian is required')
    with Path(path).open('w') as file:
        file.write('$hessian\n')
        for column in array.T:
            for begin in range(0,len(column),5):
                fields=[f'{value:15.10f}' for value in column[begin:begin+5]]
                if any(len(field)!=15 for field in fields):raise ValueError('Hessian value exceeds format width')
                file.write('     '+''.join(fields)+'\n')


def parse_xtb_thermal_eV(text):
    """Use the precise RRHO contribution, avoiding the rounded summary table."""
    matches=re.findall(r'G\(RRHO\) contrib\.\s+([-+0-9.Ee]+)\s+Eh',text)
    if not matches:raise ValueError('precise G(RRHO) contribution absent')
    result=float(matches[-1])*HARTREE_EV
    if not np.isfinite(result):raise ValueError('nonfinite RRHO contribution')
    return result
