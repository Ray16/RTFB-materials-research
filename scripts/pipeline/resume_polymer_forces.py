#!/usr/bin/env python
"""Restart a topology-checked oligomer with conservative, step-checked CDS forces.

Each spec names an immutable prior failed task, seed and derivative diagnostic.
The ordinary optimizer, active electronic SP, topology and gradient gates remain.
No minimum or free-energy claim follows from this electronic-only continuation.
"""

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()


def checked_diagnostic(path, xyz):
    report = json.loads(Path(path).read_text())
    if report["geometry_sha256"] != api.digest(xyz):
        raise ValueError("CDS force diagnostic has a different geometry")
    if report["steps_Bohr"] != [1e-4, 5e-5]:
        raise ValueError("Diagnostic must qualify the two execution steps")
    if report["maximum_step_disagreement"] > 1e-6:
        raise ValueError("CDS energy derivative is not step converged")
    checks = report.get("full_energy_difference_checks", [])
    if any(c["absolute_error"] > 1e-6 for c in checks):
        raise ValueError("Corrected total force fails full energy differences")
    return report


def prepare(spec_path):
    spec = json.loads(spec_path.read_text())
    campaign, parent = Path(spec["campaign"]), Path(spec["parent_campaign"])
    prior = api.read_plan(parent)
    api.verify_sources(prior["producer_sha256"])
    key, xyz = spec["key"], Path(spec["xyz"])
    failure = parent / "failures" / (key + ".json")
    failed = json.loads(failure.read_text())
    owner_check = api.load_module("force_owner", ROOT / "scripts/pipeline/run_polymer_spacing.py")
    if owner_check.owner_alive(failed["owner"]) or (parent / "claims" / key).exists():
        raise RuntimeError("Prior worker still owns this task")
    if "did not converge" not in failed["error"]:
        raise ValueError("Expected a preserved optimizer-convergence failure")
    task_path = parent / "live_inputs" / key / "task.json"
    task = json.loads(task_path.read_text())
    diagnostic = Path(spec["diagnostic"])
    checked_diagnostic(diagnostic, xyz)
    audit = api.load_module("force_audit", ROOT / "scripts/validation/audit_polymer_geometries.py")
    sanity = audit.audit(task, xyz)
    if sanity["issues"]:
        raise ValueError("Seed connectivity or geometry needs correction")
    campaign.mkdir()
    for name in ("claims", "tasks", "done", "failures", "workers", "launches", "pilot", "live_inputs"):
        (campaign / name).mkdir()
    directory = campaign / "live_inputs" / key
    directory.mkdir()
    shutil.copyfile(xyz, directory / "input.xyz")
    task.update(input_xyz=str(directory / "input.xyz"), input_sha256={str(xyz): api.digest(xyz),
                str(task_path): api.digest(task_path), str(directory / "input.xyz"): api.digest(directory / "input.xyz")})
    atomic_write_json(directory / "task.json", task)
    atomic_write_json(campaign / "checkpoint_geometry_review.json", sanity)
    shutil.copyfile(parent / "pilot/status.json", campaign / "pilot/status.json")
    shutil.copyfile(spec_path, campaign / "continuation_spec.json")
    inputs = [failure, task_path, xyz, diagnostic, directory / "task.json", campaign / "pilot/status.json",
              campaign / "continuation_spec.json", campaign / "checkpoint_geometry_review.json", parent / "plan.json"]
    plan = dict(prior, created=time.time(), tasks=[next(t for t in prior["tasks"] if t["key"] == key)],
                raw_root=str(ROOT / "calcs/polymer_models" / campaign.name),
                input_sha256={str(p): api.digest(p) for p in inputs},
                recovery_reason="Same SMD energy; conservative two-step checked CDS geometry derivative",
                derivative_settings=dict(steps_Bohr=[1e-4, 5e-5], maximum_step_disagreement_Ha_Bohr=1e-6,
                                         energy_function_unchanged=True, full_final_energy_gradient_check=True))
    plan["producer_sha256"].update({str(SCRIPT): api.digest(SCRIPT),
                                  str(ROOT / "src/redox/qm/cds_finite.py"): api.digest(ROOT / "src/redox/qm/cds_finite.py")})
    atomic_write_json(campaign / "plan.json", plan)
    # Branches may be studied independently; only a chosen qualified branch is published.
    if spec.get("register", True):
        pointer = ROOT / "results/polymer_models/recovery_campaigns.json"
        registry = json.loads(pointer.read_text())
        registry["campaigns"].append(dict(path=str(campaign.relative_to(ROOT)), plan_sha256=api.digest(campaign / "plan.json")))
        atomic_write_json(pointer, registry)
    print(json.dumps(dict(campaign=str(campaign), seed_geometry_QC=sanity["status"])))


@contextmanager
def checked_forces(log):
    from redox.qm.cds_finite import conservative_cds_gradient, energy_gradient
    from gpu4pyscf.solvent.grad import smd as grad
    from gpu4pyscf.solvent import smd

    with conservative_cds_gradient(log, 1e-4):
        coarse = grad.get_cds

        def evaluate(solvent):
            first = coarse(solvent)
            mol, temporary = solvent.mol.copy(), solvent.copy()

            def energy(position):
                mol.set_geom_(position, unit="Bohr")
                temporary.reset(mol)
                return smd.get_cds_legacy(temporary)[0]

            finer = energy_gradient(mol.atom_coords(unit="Bohr"), energy, 5e-5)
            difference = float(abs(first - finer).max())
            with Path(str(log) + ".step_checks.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(maximum_step_difference_Ha_Bohr=difference, passed=difference <= 1e-6)) + "\n")
            if difference > 1e-6:
                raise ValueError("Conservative CDS force step convergence failed; checkpoint preserved")
            return first

        grad.get_cds = evaluate
        try:
            yield
        finally:
            grad.get_cds = coarse


def final_energy_check(xyz, out):
    from ase.io import read
    from pyscf import gto
    from redox.qm.dft import _build_smd_mf, _free_gpu_pool
    from redox.qm.stable_geometry_scf import stabilize

    atoms = read(xyz)
    positions = atoms.positions / 0.52917721067
    out.mkdir()
    baseline = None

    def solve(coords, name, gradient=False):
        mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(), coords.tolist())), unit="Bohr",
                    basis="def2-svp", charge=0, spin=0, max_memory=6000, verbose=0)
        mf = _build_smd_mf(mol, "r2scan", "acetonitrile", "gpu", disp="d4")
        mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-12, 1e-8, 300
        stabilize(mf, out / "scf.jsonl")
        mf.kernel(dm0=baseline)
        if not mf.converged:
            raise ValueError("Final energy-difference SCF failed")
        result = dict(E_Ha=float(mf.e_tot), coordinates_Bohr=coords.tolist(), SCF_converged=True)
        density = mf.make_rdm1().get()
        force = None
        if gradient:
            method = mf.nuc_grad_method()
            method.grid_response = True
            raw = method.kernel()
            force = np.asarray(raw.get() if hasattr(raw, "get") else raw)
            result["gradient_Ha_Bohr"] = force.tolist()
        atomic_write_json(out / (name + ".json"), result)
        del mf, mol
        _free_gpu_pool()
        return result["E_Ha"], density, force

    _, baseline, force = solve(positions, "baseline", True)
    indices = np.argsort(abs(force).ravel())[-2:]
    checks = []
    for index in indices:
        atom, axis = map(int, np.unravel_index(index, force.shape))
        for step in (.0005, .001):
            energies = []
            for sign in (-1, 1):
                shifted = positions.copy()
                shifted[atom, axis] += sign * step
                energies.append(solve(shifted, f"atom{atom}_axis{axis}_h{step}_sign{sign}")[0])
            numerical = (energies[1] - energies[0]) / (2 * step)
            checks.append(dict(atom=atom, axis=axis, step_Bohr=step, analytic=float(force[atom, axis]),
                               numerical=numerical, absolute_error=abs(numerical - force[atom, axis])))
    report = dict(geometry_sha256=api.digest(xyz), checks=checks,
                  passed=all(c["absolute_error"] <= 1e-6 for c in checks), acceptance_Ha_Bohr=1e-6)
    atomic_write_json(out / "summary.json", report)
    if not report["passed"]:
        raise ValueError("Final conservative total force fails energy differences")
    return report


def runner(campaign):
    import redox.qm.dft as dft
    from redox.qm.stable_geometry_scf import stabilize

    continuation = api.load_module("force_continuation", ROOT / "scripts/pipeline/continue_polymer_models.py")
    continuation.RECOVERY, continuation.SCRIPT = campaign, SCRIPT
    base = continuation.recovery_runner()
    previous, builder = base.calculate, dft._build_smd_mf

    def build(*args, **kwargs):
        mf = builder(*args, **kwargs)
        if (args[1] if len(args) > 1 else kwargs.get("xc")) == "r2scan":
            stabilize(mf, Path(api.read_plan(campaign)["raw_root"]) / "scf_recovery.jsonl")
        return mf

    def calculate(task, plan, pilot=False):
        raw = Path(plan["raw_root"])
        raw.mkdir(parents=True, exist_ok=True)
        state = raw / task["key"] / "rest"
        try:
            with checked_forces(raw / (task["key"] + "_cds.jsonl")):
                result = previous(task, plan, pilot)
                report = final_energy_check(state / "opt.xyz", state / "energy_gradient_check")
            status = campaign / "tasks" / task["key"] / "status.json"
            record = json.loads(status.read_text())
            record["conservative_force_energy_check"] = report
            record["output_sha256"].update({str(p): api.digest(p) for p in state.rglob("*.json")})
            atomic_write_json(status, record)
            return result
        except BaseException as error:
            status = campaign / "tasks" / task["key"] / "status.json"
            if status.exists():
                record = json.loads(status.read_text())
                record.update(status="failed_preserved", error=repr(error), finished=time.time())
                atomic_write_json(status, record)
            raise

    dft._build_smd_mf = build
    base.calculate = calculate
    return base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker"))
    parser.add_argument("--spec", type=Path)
    parser.add_argument("--campaign", type=Path)
    # The inherited worker command passes role; campaign comes from its environment.
    parser.add_argument("--role", choices=("queue_large",), default="queue_large")
    args = parser.parse_args()
    if args.action == "prepare":
        if not args.spec:
            parser.error("--spec is required")
        prepare(args.spec.resolve())
    else:
        campaign = args.campaign or Path(os.environ["REDOX_FORCE_CAMPAIGN"])
        runner(campaign.resolve()).worker(args.role)


if __name__ == "__main__":
    main()
