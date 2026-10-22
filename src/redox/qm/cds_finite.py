"""Conservative finite derivatives of the unchanged SMD CDS energy function."""

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import time

import numpy as np


def energy_gradient(coordinates, energy, step=1e-4):
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Positive finite CDS gradient step required")
    coordinates = np.asarray(coordinates, dtype=float)
    result = np.empty_like(coordinates)
    for atom in range(len(coordinates)):
        for axis in range(3):
            plus = coordinates.copy()
            plus[atom, axis] += step
            minus = coordinates.copy()
            minus[atom, axis] -= step
            result[atom, axis] = (energy(plus) - energy(minus)) / (2 * step)
    if not np.isfinite(result).all():
        raise ValueError("CDS energy gradient is nonfinite")
    return result


def gradient_hessian(coordinates, gradient, step=1e-5):
    coordinates = np.asarray(coordinates, dtype=float)
    matrix = np.empty((len(coordinates), len(coordinates), 3, 3))
    for atom in range(len(coordinates)):
        for axis in range(3):
            plus = coordinates.copy()
            plus[atom, axis] += step
            minus = coordinates.copy()
            minus[atom, axis] -= step
            matrix[atom, :, axis] = (gradient(plus) - gradient(minus)) / (2 * step)
    if not np.isfinite(matrix).all():
        raise ValueError("CDS gradient Hessian is nonfinite")
    return matrix


@contextmanager
def conservative_cds_gradient(log, step=1e-4):
    """Run-local geometry fallback; energies and solvent descriptors stay unchanged.

    Use only after an independently checked derivative diagnostic. This context
    is for geometry forces, not for nesting numerical Hessian evaluations.
    """
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    from gpu4pyscf.solvent.grad import smd as gradient_module
    from gpu4pyscf.solvent import smd as solvent_module

    original = gradient_module.get_cds
    log = Path(log)

    def evaluate(solvent):
        started = time.time()
        mol = solvent.mol.copy()
        temporary = solvent.copy()
        positions = mol.atom_coords(unit="Bohr")

        def energy(coordinates):
            mol.set_geom_(coordinates, unit="Bohr")
            temporary.reset(mol)
            return solvent_module.get_cds_legacy(temporary)[0]

        legacy = np.asarray(original(solvent))
        numeric = energy_gradient(positions, energy, step)
        with log.open("a") as stream:
            stream.write(
                json.dumps(
                    dict(
                        step_Bohr=step,
                        source_coordinates_sha256=hashlib.sha256(
                            positions.tobytes()
                        ).hexdigest(),
                        maximum_legacy_difference=float(abs(numeric - legacy).max()),
                        elapsed_seconds=time.time() - started,
                        CDS_energy_function_unchanged=True,
                    )
                )
                + "\n"
            )
        return numeric

    gradient_module.get_cds = evaluate
    try:
        yield
    finally:
        gradient_module.get_cds = original


@contextmanager
def finer_cds_hessian(step=1e-5):
    """Use an already qualified finite-gradient step; admission needs its evidence."""
    if step not in (1e-5, 1e-6):
        raise ValueError("Only independently qualified CDS Hessian steps are allowed")
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    from gpu4pyscf.solvent.hessian import smd as hessian_module
    from gpu4pyscf.solvent.grad import smd as gradient_module

    original = hessian_module.get_cds

    def evaluate(solvent):
        mol = solvent.mol.copy()
        temporary = solvent.copy()

        def gradient(coordinates):
            mol.set_geom_(coordinates, unit="Bohr")
            temporary.reset(mol)
            return np.asarray(gradient_module.get_cds(temporary)).copy()

        return gradient_hessian(mol.atom_coords(unit="Bohr"), gradient, step)

    hessian_module.get_cds = evaluate
    try:
        yield
    finally:
        hessian_module.get_cds = original
