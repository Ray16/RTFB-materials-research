"""Shared actual-frame, occupied-spin and energy checks for molecular SCFs."""
import math
from pathlib import Path

import h5py
import numpy as np
from ase.io import read as read_atoms
from pyscf.lib import param

from redox.validation.scf_frame import checkpoint_frame, checkpoint_spin

SOURCE_BYTES = Path(__file__).read_bytes()


def checkpoint_series(directory, geometry, expected_settings, read, used, *,
                      expected_count=None, fixed_geometry=False, native_diagnostics=None,
                      require_all_converged=True):
    atoms = read_atoms(geometry)
    points = sorted(Path(directory).glob('*/*/settings.json'))
    if not points or (expected_count is not None and len(points) != expected_count):
        raise ValueError('The complete declared SCF checkpoint series is required')
    frame_errors, spin_rows, energies, converged_energies = [], [], [], []
    for point in points:
        settings = read(point)
        solved = read(point.with_name('SCF_result.json'))
        checkpoint = point.with_name('SCF.h5')
        used.append(checkpoint)
        if (any(settings.get(k) != v for k, v in expected_settings.items())
                or (require_all_converged and not solved['converged']) or not math.isfinite(solved['energy_Ha'])):
            raise ValueError('Actual molecular SCF settings or convergence differ')
        coordinates = atoms.positions / param.BOHR if fixed_geometry else settings['coordinates_Bohr']
        frame = checkpoint_frame(checkpoint, charge=expected_settings['charge'], mult=expected_settings['mult'],
            basis=expected_settings['basis'], coordinates_Bohr=coordinates,
            symbols=atoms.get_chemical_symbols(), maximum_error_Bohr=1e-12)
        frame_errors.append(frame['maximum_coordinate_error_Bohr'])
        spin_row = checkpoint_spin(checkpoint)
        spin_rows.append(spin_row)
        if len(spin_rows) == 1 and native_diagnostics is not None:
            if (round(spin_row['s_squared'], 4) != native_diagnostics['s_squared']
                    or round(spin_row['spin_contamination'], 4) != native_diagnostics['spin_contam']):
                raise ValueError('Occupied-orbital spin differs from the rounded producer diagnostic')
        with h5py.File(checkpoint, 'r') as stream:
            actual = float(stream['scf/e_tot'][()])
        if not math.isfinite(actual) or abs(actual - solved['energy_Ha']) > 1e-10:
            raise ValueError('Actual molecular checkpoint energy differs')
        energies.append(actual)
        if solved['converged']:
            converged_energies.append(actual)
    return dict(frame_errors=frame_errors, spins=spin_rows, energies_Ha=energies,
                converged_energies_Ha=converged_energies,
                actual_SCF_frames_checked=len(points))
