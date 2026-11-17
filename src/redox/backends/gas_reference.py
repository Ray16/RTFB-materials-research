"""Strict solver admission for fresh gas-phase reference calculations.

Convergence is separate from minimum, spin, binding and experimental validation.
These helpers do not change functionals, bases, SCF tolerances or fit any values.
"""
from contextlib import contextmanager
import fcntl
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path

import numpy as np

HARTREE_EV = 27.211386245988
SOURCE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
SOLVER_CONTRACT = dict(schema=1, finite_SCF_energy=True, SCF_converged=True,
                       optimizer_converged=True, optimization_SCF_converged=True,
                       full_grid_response=True)


def energy_Ha(mean_field):
    """Return a finite, converged SCF energy in Hartree, or refuse it."""
    value = mean_field.kernel()
    if isinstance(value, (bool, str, bytes)) or getattr(getattr(value, 'dtype', None), 'kind', None) == 'b':
        raise RuntimeError('SCF returned a nonnumeric energy')
    value = float(value)
    if not mean_field.converged or not math.isfinite(value):
        raise RuntimeError('SCF energy is unconverged or nonfinite')
    return value


def optimized_molecule(mean_field, *, maxsteps=100):
    """Retain geomeTRIC and final-SCF convergence; do not return a stopped geometry."""
    if getattr(mean_field, 'with_solvent', None) is not None:
        raise ValueError('This helper is for gas reference calculations')
    from pyscf.geomopt.geometric_solver import kernel
    gradient = mean_field.nuc_grad_method()
    # DFT quadrature response belongs in the conservative gradient, also in gas.
    if hasattr(gradient, 'grid_response'):
        gradient.grid_response = True
    scanner = gradient.as_scanner()
    converged, molecule = kernel(scanner, maxsteps=maxsteps, assert_convergence=True)
    if not converged or not scanner.converged or molecule is None:
        raise RuntimeError('Gas geometry optimizer or final SCF did not converge')
    if not np.isfinite(molecule.atom_coords()).all():
        raise RuntimeError('Gas optimizer returned nonfinite coordinates')
    return molecule


def homo_eV(mean_field):
    """Occupied HOMO in eV; explicitly transfer device arrays to host."""
    def host(value):
        if hasattr(value, 'get'):
            value = value.get()
        elif isinstance(value, (list, tuple)):
            value = [v.get() if hasattr(v, 'get') else v for v in value]
        return np.asarray(value)
    try:
        energies, occupancies = host(mean_field.mo_energy), host(mean_field.mo_occ)
        if (energies.shape != occupancies.shape or energies.ndim not in (1, 2)
                or not np.isfinite(occupancies).all()):
            return None
        occupied = energies[occupancies > 0]
        if not occupied.size or not np.isfinite(occupied).all():
            return None
        return float(occupied.max()) * HARTREE_EV
    except (AttributeError, TypeError, ValueError, RuntimeError):
        return None


def verification(identity, producer):
    software = {}
    for name in ('numpy', 'pyscf', 'geometric', 'gpu4pyscf', 'gpu4pyscf-cuda12x', 'gpu4pyscf-cuda11x'):
        try:
            software[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    return dict(identity=dict(identity), contract=dict(SOLVER_CONTRACT),
                solver_sha256=SOURCE_SHA256,
                software=software,
                producer_sha256=hashlib.sha256(Path(producer).read_bytes()).hexdigest())


def existing_result(path, identity, producer, *, required_fields=()):
    """Reuse only exactly qualified inputs; preserve every incompatible/legacy file."""
    path = Path(path)
    if not path.exists():
        return None
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise RuntimeError(f'Unreadable reference record preserved: {path}; choose a fresh output directory') from exc
    if (not isinstance(record, dict) or record.get('status') != 'ok'
            or record.get('solver_verification') != verification(identity, producer)):
        raise RuntimeError(f'Unqualified or incompatible reference record preserved: {path}; choose a fresh output directory')
    for name in required_fields:
        value = record.get(name)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise RuntimeError(f'Incomplete reference output preserved: {path} ({name}); choose a fresh output directory')
    return record


@contextmanager
def record_guard(path):
    """Exclude two writers of the same reference case across shared-node workers."""
    path = Path(path)
    directory = path.parent / '.locks'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / (path.name + '.lock')).open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Reference case has a live writer: {path}') from exc
        yield
