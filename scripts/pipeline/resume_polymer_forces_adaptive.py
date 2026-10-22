#!/usr/bin/env python
"""Source-preserving continuation with adaptive conservative CDS derivatives."""

from contextlib import contextmanager
import json
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
BASE_SOURCE = ROOT / "scripts/pipeline/resume_polymer_forces.py"
base = api.load_module("fixed_step_force_parent", BASE_SOURCE)
base.SCRIPT = SCRIPT


def adaptive_gradient(positions, energy, tolerance=1e-6):
    """Refine only unresolved components; never relax the acceptance tolerance."""
    positions = np.asarray(positions, dtype=float)
    result = np.empty_like(positions)
    evidence = []
    for atom in range(len(positions)):
        for axis in range(3):
            previous = None
            attempts = []
            for h in (1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8):
                minus, plus = positions.copy(), positions.copy()
                minus[atom, axis] -= h
                plus[atom, axis] += h
                value = (energy(plus) - energy(minus)) / (2 * h)
                if not np.isfinite(value):
                    raise ValueError("Nonfinite conservative CDS derivative")
                attempts.append(dict(step_Bohr=h, derivative_Ha_Bohr=float(value)))
                if previous is not None and abs(value - previous) <= tolerance:
                    result[atom, axis] = value
                    evidence.append(dict(atom=atom, axis=axis, attempts=attempts,
                                         final_step_difference=abs(value - previous)))
                    break
                previous = value
            else:
                raise ValueError(f"CDS component {atom}/{axis} remains step inconsistent")
    return result, evidence


@contextmanager
def adaptive_forces(log):
    from redox.workflows.resources import verify_gate_owner
    from redox.workflows.geometry_checkpoints import write_geometry
    from gpu4pyscf.solvent.grad import smd as grad
    from gpu4pyscf.solvent import smd

    verify_gate_owner()
    original = grad.get_cds
    folder = Path(str(log) + ".evaluations")
    folder.mkdir()
    index = 0

    def evaluate(solvent):
        nonlocal index
        directory = folder / f"{index:06d}"
        index += 1
        directory.mkdir()
        write_geometry(solvent.mol, directory / "input.xyz", "Conservative force input; unqualified until all gates pass")
        positions = solvent.mol.atom_coords(unit="Bohr")
        mol, temporary = solvent.mol.copy(), solvent.copy()
        started = time.time()
        status = dict(status="computing", geometry_sha256=api.digest(directory / "input.xyz"),
                      source_sha256=api.digest(SCRIPT), tolerance_Ha_Bohr=1e-6, started=started)
        atomic_write_json(directory / "status.json", status)

        def energy(coords):
            mol.set_geom_(coords, unit="Bohr")
            temporary.reset(mol)
            return smd.get_cds_legacy(temporary)[0]

        try:
            result, evidence = adaptive_gradient(positions, energy)
            np.save(directory / "gradient_Ha_Bohr.npy", result)
            atomic_write_json(directory / "component_checks.json", evidence)
            status.update(status="derivative_step_checks_passed", elapsed=time.time() - started,
                          refined_components=sum(len(c["attempts"]) > 2 for c in evidence),
                          maximum_step_difference=max(c["final_step_difference"] for c in evidence),
                          CDS_energy_function_unchanged=True)
            return result
        except BaseException as error:
            status.update(status="failed_preserved", error=repr(error), elapsed=time.time() - started)
            raise
        finally:
            atomic_write_json(directory / "status.json", status)

    grad.get_cds = evaluate
    try:
        yield
    finally:
        grad.get_cds = original


previous_prepare = base.prepare


def prepare(spec_path):
    previous_prepare(spec_path)
    campaign = Path(json.loads(spec_path.read_text())["campaign"])
    plan = api.read_plan(campaign)
    plan["producer_sha256"][str(BASE_SOURCE)] = api.digest(BASE_SOURCE)
    plan["derivative_settings"].update(adaptive_steps_Bohr=[1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8],
                                       componentwise_refinement=True, checkpoint_before_each_derivative=True)
    atomic_write_json(campaign / "plan.json", plan)
    pointer = ROOT / "results/polymer_models/recovery_campaigns.json"
    registry = json.loads(pointer.read_text())
    for entry in registry["campaigns"]:
        if entry["path"] == str(campaign.relative_to(ROOT)):
            entry["plan_sha256"] = api.digest(campaign / "plan.json")
    atomic_write_json(pointer, registry)


base.prepare = prepare
base.checked_forces = adaptive_forces


if __name__ == "__main__":
    import faulthandler

    faulthandler.enable()
    faulthandler.dump_traceback_later(600, repeat=True)
    base.main()
