#!/usr/bin/env python
"""Recover a fixed-geometry solvent-gradient check without repeating hybrid SCFs.

The earlier point calculation remains failed. Its HDF checkpoints are inspected
for identity/energy continuity, never promoted to an optimized state or G.
"""
import argparse
from contextlib import redirect_stderr, redirect_stdout
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/smd_gradient_workspace.py"
THREAD_KEYS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
               "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")


def checkpoint_continuity(parent, reference):
    import numpy as np
    from pyscf.lib import chkfile

    plan = api.read_plan(parent)
    folder = Path(plan["raw_root"]) / "pilot_execution/SCF_checkpoints"
    records = []
    for index, solvent, field in ((0, "acetonitrile", "e_smd_Ha"), (1, None, "e_gas_Ha")):
        progress = folder / f"{index:03d}" / "progress.json"
        hdf = progress.with_name("SCF.h5")
        data = json.loads(progress.read_text())
        if (api.digest(hdf) != data["checkpoint_sha256"] or not data["SCF_converged"]
                or data["qualified_result"] or data["identity"]["solvent"] != solvent):
            raise ValueError("Saved SCF identity/convergence/checksum differs")
        mol = chkfile.load_mol(str(hdf))
        identity = data["identity"]
        if (mol.charge != identity["charge"] or mol.spin + 1 != identity["mult"]
                or list(mol.elements) != identity["atoms"] or repr(mol.basis) != identity["basis"]
                or not np.allclose(mol.atom_coords(), identity["coordinates_Bohr"], rtol=0, atol=1e-12)):
            raise ValueError("HDF molecule differs from the saved point identity")
        energy = float(chkfile.load(str(hdf), "scf/e_tot"))
        difference = abs(energy - reference["result"][field]) * 27.211386245988
        if abs(energy - data["E_Ha"]) > 1e-10 or difference > 1e-6:
            raise ValueError("HDF energy fails continuity with preserved point calculation")
        records.append(dict(progress=str(progress), progress_sha256=api.digest(progress),
            HDF=str(hdf), HDF_sha256=api.digest(hdf), E_Ha=energy,
            difference_from_previous_eV=difference, optimized_state_qualified=False))
    return records


def checked_controls(folder, budget):
    helper = api.load_module("workspace_control_policy", ROOT / "scripts/validation/qualify_polymer_workspace.py")
    paths = [str(folder / "controls" / str(spin) / "status.json") for spin in (1, 2)]
    proof = helper.checked_controls(paths, budget)
    for spin in (1, 2):
        path = folder / "controls" / str(spin) / "compact/execution/status.json"
        data = json.loads(path.read_text())
        if (data.get("gradient_solver_calls", 0) < 1
                or data.get("direct_LU_solver_retained") is not True
                or data["source_sha256"].get(str(ADAPTER)) != api.digest(ADAPTER)):
            raise ValueError("Bounded-gradient implementation was not exercised")
        proof[str(path)] = api.digest(path)
    return proof


def prepare(args):
    parent = args.parent.resolve()
    old = api.read_plan(parent)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    if len(old["tasks"]) != 1 or old["tasks"][0]["charge"] != 0 or old["tasks"][0]["mult"] != 1:
        raise ValueError("This recovery is scoped to one existing neutral singlet")
    failure = parent / "footprint_qualification.json"
    if json.loads(failure.read_text())["status"] != "footprint_attempt_failed_preserved":
        raise ValueError("An observed failed point attempt is required")
    reference = Path(old["source_failed_resource_pilot"]).parent / "energy_and_gradient.json"
    data = json.loads(reference.read_text())
    task = old["tasks"][0]
    ready = json.loads((parent / "live_inputs" / task["key"] / "task.json").read_text())
    if api.digest(data["result"]["xyz"]) != api.digest(ready["input_xyz"]):
        raise ValueError("Previous gradient and recovery warm start differ")
    continuity = checkpoint_continuity(parent, data)
    controls = checked_controls(args.controls.resolve(), old["execution_workspace_GiB"])
    out = args.campaign.resolve()
    out.mkdir(parents=True, exist_ok=False)
    xyz = out / "input.xyz"
    shutil.copyfile(ready["input_xyz"], xyz)
    sources = dict(old["producer_sha256"], **controls)
    for p in (Path(__file__).resolve(), ADAPTER, ROOT / "scripts/validation/qualify_polymer_workspace.py"):
        sources[str(p)] = api.digest(p)
    inputs = {str(p): api.digest(p) for p in (parent / "plan.json", failure, reference, xyz)}
    for r in continuity:
        for key in ("progress", "HDF"):
            inputs[r[key]] = r[key + "_sha256"]
    plan = dict(old, raw_root=str(ROOT / "calcs/validation" / out.parent.name / out.name),
        producer_sha256=sources, input_sha256=inputs, input_xyz=str(xyz),
        parent_failed_point=str(parent), previous_reference=str(reference),
        gradient_controls=str(args.controls.resolve()), SCF_continuity=continuity,
        purpose="Fixed-geometry gradient and sampled total-energy derivative resource qualification only",
        repeated_hybrid_SCFs=False, full_large_geometry_workflow_qualified=False,
        native_Hessian_qualified=False, free_energy_computed=False)
    atomic_write_json(out / "plan.json", plan)
    print(json.dumps(dict(status="gradient_recovery_prepared", SCF_checkpoints_verified=2,
        hybrid_SCFs_to_repeat=0, exact_task=task["key"])))


def run(args):
    from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability
    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != "2" for k in THREAD_KEYS):
        raise ValueError("Two disjoint CPU cores and native caps=2 required")
    campaign = args.campaign.resolve()
    plan = api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    checked_controls(Path(plan["gradient_controls"]), plan["execution_workspace_GiB"])
    host_capacity(2, minimum_ram_gib=32)
    if gpu_suitability(int(os.environ["CUDA_VISIBLE_DEVICES"]))["memory_MiB"] < plan["minimum_large_GPU_MiB"]:
        raise ValueError("Device below the qualified test profile")
    path = campaign / "gradient_qualification.json"
    if path.exists():
        raise FileExistsError("Preserve the previous gradient attempt")
    claim = campaign / "gradient_claim"
    claim.mkdir()
    owner = dict(host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    atomic_write_json(claim / "owner.json", owner)
    report = dict(owner, status="gradient_check_starting", started=time.time(),
        CPU_affinity=sorted(os.sched_getaffinity(0)), GPU=os.environ["CUDA_VISIBLE_DEVICES"],
        full_large_geometry_workflow_qualified=False, native_Hessian_qualified=False,
        free_energy_computed=False, synthesis_admission=False)
    memory = None
    try:
        atomic_write_json(path, report)
        raw = Path(plan["raw_root"])
        raw.mkdir(parents=True, exist_ok=False)
        driver = api.load_module("isolated_gradient_geometry", ROOT / "scripts/pipeline/refine_polymer_conformers.py")
        engine = api.load_module("isolated_gradient_monitor", ROOT / "scripts/pipeline/run_pi_stacking.py")
        adapter = api.load_module("qualified_gradient_workspace", ADAPTER)
        with engine.MemoryMonitor() as memory, (raw / "calculation.log").open("x") as log, redirect_stdout(log), redirect_stderr(log):
            with driver.geometry_protocol(campaign, plan["tasks"][0]) as forces:
                with adapter.compact_smd_execution(raw / "execution", budget_GiB=plan["execution_workspace_GiB"]):
                    checks = forces.base.final_energy_check(Path(plan["input_xyz"]), raw / "energy_force_check")
        baseline = json.loads((raw / "energy_force_check/baseline.json").read_text())
        maximum = max(abs(v) for row in baseline["gradient_Ha_Bohr"] for v in row)
        previous = json.loads(Path(plan["previous_reference"]).read_text())["gradient_max_Ha_Bohr"]
        difference = abs(maximum - previous)
        if difference > 1e-6 or memory.peak > plan["maximum_memory_pilot_MiB"]:
            raise ValueError("Starting-gradient continuity or device-workspace bound failed")
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        report.update(status="fixed_geometry_gradient_and_sampled_derivatives_passed", checks=checks,
            gradient_max_Ha_Bohr=maximum, previous_gradient_max_difference=difference,
            sampled_components=2, displacement_steps_Bohr=[.0005, .001],
            output_sha256={str(p): api.digest(p) for p in raw.rglob("*") if p.is_file()})
    except BaseException as error:
        report.update(status="gradient_attempt_failed_preserved", error=repr(error))
        raise
    finally:
        report.update(finished=time.time(), observed_peak_GPU_MiB=memory.peak if memory else None)
        atomic_write_json(path, report)
        (claim / "owner.json").unlink()
        claim.rmdir()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--controls", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        if not args.parent or not args.controls:
            parser.error("Preparation requires the failed parent and passed controls")
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
