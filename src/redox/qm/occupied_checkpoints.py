"""Durable, occupied-orbital checkpoints for restricted single-point SCF.

All occupied coefficients, occupations and orbital energies retain their native
precision. Virtual coefficients are explicitly absent. This represents the same
Slater determinant and density, and supports a density restart; it is not a
full-orbital Hessian/response restart. No solver or physical setting is changed.
"""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import time

import h5py
import numpy as np

FORMAT = 'redox.restricted_occupied_SCF.v1'


def host(value):
    return np.asarray(value.get() if hasattr(value, 'get') else value)


def save(path, molecule, coefficients, occupations, energies, energy, *, converged, cycle):
    """Atomically replace one checkpoint, retaining the preceding complete file."""
    path = Path(path)
    occ = host(occupations)
    if (molecule.spin != 0 or occ.ndim != 1 or not np.isfinite(occ).all()
            or not np.all((occ == 0) | (occ == 2))
            or int(occ.sum()) != molecule.nelectron):
        raise ValueError('Integer closed-shell restricted occupations required')
    indices = np.flatnonzero(occ > 0)
    # Slice on the device before transfer: never serialize a full virtual space.
    selected = host(coefficients[:, indices])
    levels = host(energies)
    if (selected.shape != (molecule.nao_nr(), len(indices))
            or levels.shape != occ.shape or not np.isfinite(selected).all()
            or not np.isfinite(levels).all() or not np.isfinite(energy)):
        raise ValueError('Finite atom/basis-matched restricted orbitals required')
    temporary = path.with_name('.' + path.name + '.' + str(os.getpid()) + '.tmp')
    with h5py.File(temporary, 'w') as stream:
        stream.attrs.update(format=FORMAT, virtual_coefficients_retained=False,
                            density_restart_supported=True, response_restart_supported=False)
        stream.create_dataset('mol', data=molecule.dumps())
        for name, value in dict(mo_coeff_occupied=selected, occupied_indices=indices,
                                mo_occ=occ, mo_energy=levels, e_tot=float(energy),
                                converged=bool(converged), cycle=int(cycle)).items():
            stream.create_dataset('scf/' + name, data=value)
        stream.flush()
    temporary.replace(path)


def read(path):
    """Read a determinant without silently treating it as a full SCF checkpoint."""
    from pyscf.lib import chkfile
    with h5py.File(path, 'r') as stream:
        if (stream.attrs.get('format') != FORMAT
                or bool(stream.attrs.get('virtual_coefficients_retained', True))
                or not bool(stream.attrs.get('density_restart_supported', False))
                or bool(stream.attrs.get('response_restart_supported', True))
                or 'scf/mo_coeff' in stream):
            raise ValueError('Explicit occupied-only density-restart format required')
        data = {name: stream['scf/' + name][()] for name in
                ('mo_coeff_occupied', 'occupied_indices', 'mo_occ', 'mo_energy',
                 'e_tot', 'converged', 'cycle')}
    molecule = chkfile.load_mol(str(path))
    occ = data['mo_occ']; indices = data['occupied_indices']; coeff = data['mo_coeff_occupied']
    if (molecule.spin != 0 or occ.ndim != 1 or not np.all((occ == 0) | (occ == 2))
            or occ.sum() != molecule.nelectron or not np.array_equal(indices, np.flatnonzero(occ > 0))
            or coeff.shape != (molecule.nao_nr(), len(indices))
            or data['mo_energy'].shape != occ.shape
            or any(not np.isfinite(data[k]).all() for k in
                   ('mo_coeff_occupied', 'mo_occ', 'mo_energy', 'e_tot'))):
        raise ValueError('Occupied-only determinant dimensions/electron counts differ')
    return molecule, data


def density(path):
    """Return the exact saved density suitable for kernel(dm0=...)."""
    _, data = read(path)
    coeff = data['mo_coeff_occupied']
    return (coeff * data['mo_occ'][data['occupied_indices']]) @ coeff.conj().T


def review(path):
    """Independently check occupied orthonormality and restricted global spin."""
    molecule, data = read(path)
    coeff = data['mo_coeff_occupied']
    metric = molecule.intor_symmetric('int1e_ovlp')
    error = float(np.max(abs(coeff.conj().T @ metric @ coeff - np.eye(coeff.shape[1]))))
    if not error <= 1e-7:
        raise ValueError('Occupied orbitals fail the actual AO-metric check')
    return dict(energy_Ha=float(data['e_tot']), converged=bool(data['converged']),
                electrons=molecule.nelectron, occupied_orbitals=coeff.shape[1],
                basis_functions=molecule.nao_nr(), maximum_orthonormality_error=error,
                restricted_s_squared=0., virtual_coefficients_retained=False,
                density_restart_supported=True, response_restart_supported=False)


def _json(path, data):
    path = Path(path); temporary = path.with_name('.' + path.name + '.' + str(os.getpid()) + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n'); temporary.replace(path)


class Context:
    def __init__(self, directory):
        self.directory = Path(directory); self.directory.mkdir(parents=True, exist_ok=False)
        self.current = None; self.counter = 0

    def begin(self, mf):
        folder = self.directory / f'{self.counter:04d}'; self.counter += 1; folder.mkdir()
        self.current = folder
        _json(folder / 'settings.json', dict(format=FORMAT, molecule=mf.mol.dumps(),
            charge=mf.mol.charge, mult=mf.mol.spin + 1, xc=getattr(mf, 'xc', None),
            disp=getattr(mf, 'disp', None), nlc=getattr(mf, 'nlc', None),
            conv_tol=mf.conv_tol, conv_tol_grad=getattr(mf, 'conv_tol_grad', None),
            max_cycle=mf.max_cycle, init_guess=mf.init_guess,
            grid_level=getattr(getattr(mf, 'grids', None), 'level', None),
            solvent=getattr(getattr(mf, 'with_solvent', None), 'solvent', None),
            level_shift=mf.level_shift, damp=mf.damp, started=time.time(), qualified=False))
        return folder


def attach(mf, directory):
    """Observe every kernel and iteration without retaining mean-field objects."""
    if mf.mol.spin != 0:
        raise ValueError('This storage profile supports restricted singlets only')
    mf._redox_occupied_context = Context(directory)
    mf._keys = set(mf._keys) | {'_redox_occupied_context'}
    # PySCF's normal dump_chk would retain the full virtual coefficient array.
    # This format intentionally supplies a density restart through density().
    mf.chkfile = None
    previous = getattr(mf, 'callback', None)
    context = mf._redox_occupied_context

    def callback(env):
        if callable(previous):
            previous(env)
        if context.current is None:
            raise RuntimeError('SCF iteration arrived before its kernel checkpoint')
        save(context.current / 'occupied.h5', env['mol'], env['mo_coeff'], env['mo_occ'],
             env['mo_energy'], float(env['e_tot']), converged=False, cycle=int(env['cycle']))
        _json(context.current / 'progress.json', dict(cycle=int(env['cycle']),
              energy_Ha=float(env['e_tot']), updated=time.time(), qualified=False))

    mf.callback = callback
    cls = type(mf); original = cls.kernel
    if not getattr(original, '_redox_occupied_wrapper', False):
        def kernel(self, *args, **kwargs):
            active = getattr(self, '_redox_occupied_context', None)
            if active is None:
                return original(self, *args, **kwargs)
            folder = active.begin(self)
            result = original(self, *args, **kwargs)
            save(folder / 'occupied.h5', self.mol, self.mo_coeff, self.mo_occ, self.mo_energy,
                 float(self.e_tot), converged=bool(self.converged), cycle=-1)
            _json(folder / 'result.json', dict(converged=bool(self.converged),
                  energy_Ha=float(self.e_tot), finished=time.time(), qualified=False))
            return result
        kernel._redox_occupied_wrapper = True; cls.kernel = kernel
    return mf


@contextmanager
def scf_checkpoints(directory):
    """Keep the existing multiple-guess solver and record each factory/solve."""
    from redox.qm import dft
    original = dft._lowest_scf; robust = dft._kernel_robust; counter = 0
    def lowest(factory, *args, **kwargs):
        def build():
            nonlocal counter
            folder = Path(directory) / f'{counter:04d}'; counter += 1
            return attach(factory(), folder)
        return original(build, *args, **kwargs)
    def solve(mf, *args, **kwargs):
        result = robust(mf, *args, **kwargs)
        context = getattr(mf, '_redox_occupied_context', None)
        if context is not None:
            # The existing solver copies a Newton fallback's final orbitals
            # back onto mf. Record that actual result even if the Newton class
            # overrides kernel and bypasses the ordinary iteration wrapper.
            folder = context.directory / 'final'; folder.mkdir()
            save(folder / 'occupied.h5', mf.mol, mf.mo_coeff, mf.mo_occ, mf.mo_energy,
                 float(mf.e_tot), converged=bool(mf.converged), cycle=-1)
            if abs(float(mf.e_tot) - result[0]) > 1e-12:
                raise ValueError('Saved final determinant differs from the returned SCF energy')
        return result
    dft._lowest_scf = lowest; dft._kernel_robust = solve
    try:
        yield
    finally:
        dft._lowest_scf = original; dft._kernel_robust = robust
