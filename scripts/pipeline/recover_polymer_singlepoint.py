#!/usr/bin/env python
"""Replay-checked SMD allocation repair for the one missing spacing conformer.

Reuse the saved converged geometry. A same-geometry control must reproduce its
published energies and independent gradient before the new target is admitted.
"""

import argparse
from contextlib import redirect_stdout, redirect_stderr
import fcntl
import json
import os
from pathlib import Path
import shutil
import socket
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

SCRIPT = Path(__file__).resolve()
MODULE = ROOT / "src/redox/qm/compact_smd.py"
PARENT = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_continuation"
CAMPAIGN = (
    ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_singlepoint_repair"
)
KEY = "occ_none__pose1"


def prepare():
    old = api.read_plan(PARENT)
    api.verify_sources(old["producer_sha256"])
    control_receipt = PARENT / "energy/done/occ_none__pose0.json"
    baseline = verify_state_checkpoint(control_receipt)
    failure = PARENT / "energy/failures" / f"{KEY}.json"
    failed = json.loads(failure.read_text())
    owner = api.load_module(
        "sp_owner", ROOT / "scripts/pipeline/run_polymer_spacing.py"
    )
    if owner.owner_alive(failed["owner"]) or (PARENT / "energy/claims" / KEY).exists():
        raise RuntimeError("The parent failure still has a live owner")
    if "out of memory" not in failed["error"]:
        raise ValueError("This repair addresses the recorded allocation failure only")
    CAMPAIGN.mkdir()
    (CAMPAIGN / "energy/done").mkdir(parents=True)
    cases = {}
    inputs = [control_receipt, failure, PARENT / "plan.json"]
    for role, key in (("control", "occ_none__pose0"), ("target", KEY)):
        ready_path = PARENT / "live_inputs" / key / "task.json"
        task = json.loads(ready_path.read_text())
        source = Path(old["raw_root"]) / "energy" / key / "rest"
        xyz = source / "opt.xyz"
        returned = source / "geometry_checkpoints/optimizer_returned.xyz"
        if "convergence=True" not in returned.read_text().splitlines()[1]:
            raise ValueError("The saved optimizer did not report convergence")
        geometry = api.load_module(
            "sp_audit", ROOT / "scripts/validation/audit_polymer_geometries.py"
        ).audit(task, xyz)
        if geometry["issues"]:
            raise ValueError("Saved geometry fails independent sanity checks")
        state = ROOT / "calcs/polymer_models" / CAMPAIGN.name / role / "rest"
        state.mkdir(parents=True)
        shutil.copyfile(xyz, state / "opt.xyz")
        cases[role] = dict(
            task=task,
            state=str(state),
            source_geometry=str(xyz),
            optimizer_returned=str(returned),
            geometry_review=geometry,
        )
        inputs += [ready_path, xyz, returned, state / "opt.xyz"]
    plan = dict(
        tasks=[next(t for t in old["tasks"] if t["key"] == KEY)],
        cases=cases,
        baseline=baseline,
        baseline_result=json.loads(
            (Path(baseline["state_dir"]) / "result.json").read_text()
        ),
        producer_sha256=dict(
            old["producer_sha256"],
            **{str(SCRIPT): api.digest(SCRIPT), str(MODULE): api.digest(MODULE)},
        ),
        input_sha256={str(path): api.digest(path) for path in inputs},
        parent=str(PARENT),
        protocol=(
            "Unchanged active electronic protocol, same geometry; "
            "independent gradient required"
        ),
        tolerances=dict(energy_eV=1e-5, gradient_Ha_Bohr=2e-7),
    )
    atomic_write_json(CAMPAIGN / "plan.json", plan)


def finalize():
    with (CAMPAIGN / "admission.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        paths = [CAMPAIGN / f"{role}.json" for role in ("control", "target")]
        if not all(path.exists() for path in paths):
            return
        control, target = [json.loads(path.read_text()) for path in paths]
        if (
            control["status"] != "control_replay_passed"
            or target["status"] != "computed_pending_control"
        ):
            return
        plan = api.read_plan(CAMPAIGN)
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(control["output_sha256"])
        api.verify_sources(target["output_sha256"])
        report = dict(
            target, status="state_records_verified", control_sha256=api.digest(paths[0])
        )
        directory = CAMPAIGN / "energy/tasks" / KEY
        directory.mkdir(parents=True, exist_ok=True)
        atomic_write_json(directory / "status.json", report)
        receipt = dict(
            key=KEY,
            status_path=str(directory / "status.json"),
            status_sha256=api.digest(directory / "status.json"),
        )
        atomic_write_json(CAMPAIGN / "energy/done" / f"{KEY}.json", receipt)
        verify_state_checkpoint(CAMPAIGN / "energy/done" / f"{KEY}.json")
        atomic_write_json(
            CAMPAIGN / "summary.json",
            dict(
                status="replay_checked_singlepoint_complete",
                control=control["differences"],
                target_E_smd_eV=target["E_smd_eV"],
                receipt=str(CAMPAIGN / "energy/done" / f"{KEY}.json"),
            ),
        )


def worker(role):
    from redox.workflows.resources import verify_gate_owner
    from redox.qm.compact_smd import compact_solvent_build
    from redox.qm.dft import dft_smd
    from redox.qm.sp import ensure
    from redox.core.protocol import load_record
    import pyscf.lib

    verify_gate_owner()
    pyscf.lib.num_threads(2)
    if len(os.sched_getaffinity(0)) != 2:
        raise RuntimeError("Two disjoint CPU slots required")
    plan = api.read_plan(CAMPAIGN)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    case = plan["cases"][role]
    state = Path(case["state"])
    status_path = CAMPAIGN / f"{role}.json"
    if status_path.exists():
        raise RuntimeError("Preserve the previous attempt")
    report = dict(
        status="computing",
        role=role,
        task=case["task"],
        host=socket.gethostname(),
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        cpus=sorted(os.sched_getaffinity(0)),
        gpu=os.environ["CUDA_VISIBLE_DEVICES"],
        started=time.time(),
    )
    atomic_write_json(status_path, report)
    print(f"Actual compact SMD {role} start", flush=True)
    backend = api.load_module(
        "sp_memory_backend", ROOT / "scripts/pipeline/run_pi_stacking.py"
    )
    try:
        with (
            backend.MemoryMonitor() as memory,
            (state / "calculation.log").open("w") as stream,
            redirect_stdout(stream),
            redirect_stderr(stream),
        ):
            with compact_solvent_build() as source:
                (state / "executed_smd_build.py").write_text(source)
                raw = dft_smd(
                    state / "opt.xyz", 0, 1, do_opt=False, do_freq=False, backend="gpu"
                )
                gradient = backend.gradient(state / "opt.xyz")
            if (
                not raw.get("converged_smd")
                or not raw.get("converged_gas")
                or gradient > 0.00045
            ):
                raise ValueError("SP convergence or independent stationarity failed")
            raw.update(
                optimized=True,
                geometry_optimization_converged=True,
                optimization_SCF_converged=True,
                opt_out=str(state / "opt.xyz"),
                geometry_source_review=case["geometry_review"],
                geometry_source_evidence=plan["input_sha256"],
                execution_compact_SMD=True,
            )
            atomic_write_json(state / "result.json", raw)
            ensure(case["task"]["key"], "rest", state, 0, 1, backend="gpu")
            active = load_record(state, 0, 1)
            if active is None:
                raise ValueError("Active SP record not admitted")
            report.update(
                E_smd_eV=active["e_smd_eV"],
                gradient_max_Ha_Bohr=gradient,
                state_dir=str(state),
                physical_QC="passed",
            )
            if role == "control":
                differences = dict(
                    smd_eV=abs(raw["e_smd_eV"] - plan["baseline"]["E_smd_eV"]),
                    gas_eV=abs(raw["e_gas_eV"] - plan["baseline_result"]["e_gas_eV"]),
                    gradient_Ha_Bohr=abs(
                        gradient - plan["baseline"]["gradient_max_Ha_Bohr"]
                    ),
                )
                if (
                    max(differences["smd_eV"], differences["gas_eV"])
                    > plan["tolerances"]["energy_eV"]
                    or differences["gradient_Ha_Bohr"]
                    > plan["tolerances"]["gradient_Ha_Bohr"]
                ):
                    raise ValueError(f"Same-geometry replay failed: {differences}")
                report.update(status="control_replay_passed", differences=differences)
            else:
                report["status"] = "computed_pending_control"
            report["output_sha256"] = {
                str(p): api.digest(p)
                for p in state.rglob("*")
                if p.is_file() and p.suffix in (".json", ".xyz", ".py")
            }
        report["observed_peak_GPU_MiB"] = memory.peak
        if memory.peak <= 0:
            raise RuntimeError("GPU memory observation unavailable")
        api.verify_sources(plan["producer_sha256"])
    except BaseException as error:
        report.update(status="failed_preserved", error=repr(error))
        raise
    finally:
        report["finished"] = time.time()
        atomic_write_json(status_path, report)
    finalize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker", "finalize"))
    parser.add_argument("--role", choices=("control", "target"))
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "worker":
        if args.role is None:
            parser.error("worker requires a role")
        worker(args.role)
    else:
        finalize()
