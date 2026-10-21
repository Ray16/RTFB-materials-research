#!/usr/bin/env python
"""Resume the frozen spacing cohort with checked solvent-memory execution only.

Original energy receipts and failures remain immutable. Full singlet/radical
Hessian replays gate a fresh large-model thermal pilot before queue expansion.
"""

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import socket
import sys
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
PARENT = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_continuation"
CAMPAIGN = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_memory_resume"
ADAPTER = ROOT / (
    "results/validation/resume_campaign_2026-10-12/physical_campaign/"
    "pcm_memory/host_pcm_adapter.py"
)


def area_contraction(tensor, left, right, xp=np):
    """The installed d2A @ (left*right), without a full strided-tensor copy."""
    if tensor.ndim != 5 or tensor.shape[2:4] != (3, 3):
        raise ValueError("Unexpected analytic area tensor")
    vector = left * right
    result = xp.empty(tensor.shape[:-1], dtype=tensor.dtype)
    for atom in range(tensor.shape[0]):
        block = xp.ascontiguousarray(tensor[atom])
        result[atom] = block @ vector
        del block
    return result


@contextmanager
def solvent_memory(metrics):
    """Use unchanged analytic PCM equations, with independently checked sources."""
    import cupy as cp
    from gpu4pyscf.solvent.hessian import pcm
    from redox.qm.hessian_memory import pcm_certificate
    from redox.qm.native_thermal import NATIVE_THERMAL

    _, certificate = pcm_certificate(NATIVE_THERMAL)
    adapter = api.load_module("polymer_host_pcm_memory", ADAPTER)
    if api.digest(ADAPTER) != certificate["adapter_sha256"]:
        raise RuntimeError("Host PCM adapter differs from its qualification")
    original = pcm.get_v_dot_d2A_dot_q
    metrics.update(area_calls=0, maximum_area_block_bytes=0, host_curvature={})

    def blocked(tensor, left, right):
        metrics["area_calls"] += 1
        metrics["maximum_area_block_bytes"] = max(
            metrics["maximum_area_block_bytes"], tensor[0].nbytes
        )
        return area_contraction(tensor, left, right, cp)

    pcm.get_v_dot_d2A_dot_q = blocked
    try:
        with adapter.host_curvature(metrics["host_curvature"]):
            for field in ("cpu_pcm_library_sha256", "gpu_pcm_library_sha256"):
                if metrics["host_curvature"][field] != certificate[field]:
                    raise RuntimeError("Installed PCM source differs from qualification")
            yield
    finally:
        pcm.get_v_dot_d2A_dot_q = original


def prepare():
    from redox.workflows.checkpoints import verify_state_checkpoint

    plan = api.read_plan(PARENT)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    for stage in ("energy", "thermal"):
        if any((PARENT / stage / "claims").iterdir()):
            raise RuntimeError("Preserve the live parent campaign")
    ownership = api.load_module("thermal_parent_ownership", ROOT / "scripts/pipeline/run_polymer_spacing.py")
    previous = [PARENT / "coordinator.json"] + list(PARENT.glob("*/workers/*.json"))
    if any(ownership.owner_alive(json.loads(path.read_text())) for path in previous):
        raise RuntimeError("A parent coordinator or worker is still live")
    if json.loads((PARENT / "thermal/pilot/status.json").read_text()).get(
        "error"
    ) != "CUDARuntimeError('cudaErrorMemoryAllocation: out of memory')":
        raise ValueError("This continuation requires the reviewed memory failure")
    CAMPAIGN.mkdir()
    for stage in ("energy", "thermal"):
        for folder in ("done", "claims", "tasks", "failures", "workers", "pilot", "launches"):
            (CAMPAIGN / stage / folder).mkdir(parents=True)
    for name in ("cases.json", "reaction_terms.json"):
        shutil.copyfile(PARENT / name, CAMPAIGN / name)
    shutil.copytree(PARENT / "live_inputs", CAMPAIGN / "live_inputs")
    borrowed = {}
    for receipt in (PARENT / "energy/done").glob("*.json"):
        verify_state_checkpoint(receipt)
        shutil.copyfile(receipt, CAMPAIGN / "energy/done" / receipt.name)
        borrowed[str(receipt)] = api.digest(receipt)
    # Keep the exhausted alternate conformer explicit. It is not silently retried
    # or counted as a completed conformer; all six graphs have qualified samples.
    for failure in (PARENT / "energy/failures").glob("*.json"):
        shutil.copyfile(failure, CAMPAIGN / "energy/failures" / failure.name)
        borrowed[str(failure)] = api.digest(failure)
    shutil.copyfile(PARENT / "energy/pilot/status.json", CAMPAIGN / "energy/pilot/status.json")
    plan.update(
        created=time.time(),
        parent_campaign=str(PARENT),
        raw_root=str(ROOT / "calcs/polymer_models" / CAMPAIGN.name),
        execution_change="Atom-block area contraction and certified host PCM curvature; unchanged equations, surface and physical gates",
        borrowed_energy_receipts=borrowed,
    )
    plan["producer_sha256"].update({str(SCRIPT): api.digest(SCRIPT), str(ADAPTER): api.digest(ADAPTER)})
    plan["input_sha256"].update(borrowed)
    for path in (PARENT / "plan.json", PARENT / "thermal/pilot/status.json"):
        plan["input_sha256"][str(path)] = api.digest(path)
    for path in CAMPAIGN.rglob("*.json"):
        plan["input_sha256"][str(path)] = api.digest(path)
    atomic_write_json(CAMPAIGN / "plan.json", plan)
    api.collect(CAMPAIGN)


def replay(mult):
    from redox.qm import native_thermal as native
    from redox.qm.hessian_memory import pcm_certificate
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    import cupy as cp
    import pyscf.lib

    pyscf.lib.num_threads(2)
    if len(os.sched_getaffinity(0)) != 2:
        raise RuntimeError("Replay requires two disjoint CPU slots")
    plan = api.read_plan(CAMPAIGN)
    api.verify_sources(plan["producer_sha256"])
    directory = CAMPAIGN / "controls" / str(mult)
    directory.mkdir(parents=True)
    status = dict(status="computing", mult=mult, host=socket.gethostname(), pid=os.getpid(),
                  process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                  cpus=sorted(os.sched_getaffinity(0)), gpu=os.environ["CUDA_VISIBLE_DEVICES"], started=time.time())
    atomic_write_json(directory / "status.json", status)
    print(f"Actual full native Hessian replay start: multiplicity {mult}", flush=True)
    try:
        rng = np.random.default_rng(20261021)
        tensor = cp.asarray(rng.normal(size=(129, 11, 11, 3, 3))).transpose(1, 2, 3, 4, 0)
        left, right = cp.asarray(rng.normal(size=(2, 129)))
        maximum = float(cp.max(cp.abs(area_contraction(tensor, left, right, cp) - tensor @ (left * right))).get())
        if maximum > 1e-10:
            raise ValueError("GPU strided area contraction replay failed")
        status["GPU_primitive_maximum_difference"] = maximum
        del tensor, left, right
        _, certificate = pcm_certificate(native.NATIVE_THERMAL)
        control = next(row for row in certificate["controls"] if row["mult"] == mult)
        reference = next(ROOT / name for name in control["evidence_sha256"]
                         if "/native_thermal/" in name and name.endswith(".json"))
        source = reference.parent.parent
        work = ROOT / "calcs/polymer_models" / CAMPAIGN.name / "controls" / str(mult)
        work.mkdir(parents=True)
        for name in ("opt.xyz", "result.json"):
            shutil.copyfile(source / name, work / name)
        baseline = json.loads(reference.read_text())
        metrics = {}
        with solvent_memory(metrics):
            _, record = native.ensure_native_thermal(
                baseline["id"], baseline["state"], state_dir=work,
                memory_cap_gib=24, qualification=True,
            )
        path, _ = native.cache_path(work, baseline["identity"]["charge"], mult)
        matrix = path.with_suffix("") / "hessian_cartesian.npy"
        expected = np.load(reference.with_suffix("") / "hessian_cartesian.npy")
        actual = np.load(matrix)
        difference = actual - expected
        errors = dict(relative_matrix_error=float(np.linalg.norm(difference) / np.linalg.norm(expected)),
                      maximum_matrix_difference=float(abs(difference).max()),
                      maximum_thermal_difference_eV=max(abs(record["nuclear_by_cutoff_eV"][cutoff] - baseline["nuclear_by_cutoff_eV"][cutoff]) for cutoff in ("25", "50", "100")))
        if any(not np.isfinite(value) or value > 1e-6 for value in errors.values()):
            raise ValueError(f"Full native replay differs: {errors}")
        if not metrics["area_calls"] or not metrics["host_curvature"]["uploaded_blocks"]:
            raise ValueError("Memory adapters were not exercised")
        api.verify_sources(plan["producer_sha256"])
        status.update(status="full_native_replay_passed", native_gates_passed=True,
                      errors=errors, metrics=metrics,
                      evidence_sha256={str(p): api.digest(p) for p in (path, matrix, reference)})
    except BaseException as error:
        status.update(status="failed_preserved", error=repr(error))
        raise
    finally:
        status["finished"] = time.time()
        atomic_write_json(directory / "status.json", status)


def checked_replays():
    paths = [CAMPAIGN / "controls" / str(mult) / "status.json" for mult in (1, 2)]
    rows = [json.loads(path.read_text()) for path in paths]
    if any(row["status"] != "full_native_replay_passed" or row.get("native_gates_passed") is not True for row in rows):
        raise RuntimeError("Both spin controls must pass full native replay")
    for row in rows:
        if set(row["errors"]) != {"relative_matrix_error", "maximum_matrix_difference", "maximum_thermal_difference_eV"} or not row["evidence_sha256"]:
            raise ValueError("Incomplete replay evidence")
        api.verify_sources(row["evidence_sha256"])
        if any(not np.isfinite(value) or value > 1e-6 for value in row["errors"].values()):
            raise ValueError("Replay tolerance failed")
    return {str(path): api.digest(path) for path in paths}


def runner():
    import redox.qm.native_thermal as native
    from redox.qm.hessian_memory import pcm_certificate
    from redox.qm.projected_thermal import ensure_projected_thermal
    import redox.qm.projected_thermal as projected

    replays = checked_replays()
    base = api.load_module("polymer_memory_runner", ROOT / "scripts/pipeline/run_polymer_spacing.py")
    base.SCRIPT = SCRIPT
    original = native.ensure_native_thermal

    def evaluate(*args, **kwargs):
        metrics = {}
        state = Path(kwargs["state_dir"])
        proof = state / "solvent_memory_execution.json"
        prior = projected.ensure_projected_thermal

        def seal(path, geometry):
            record = json.loads(Path(path).read_text())
            atomic_write_json(proof, dict(metrics=metrics, replay_evidence_sha256=replays,
                                         adapter_sha256=api.digest(SCRIPT), scientific_model_changed=False))
            seal_hash, certificate = pcm_certificate(native.NATIVE_THERMAL)
            evidence = {str(p.relative_to(ROOT)): api.digest(p) for p in (proof, Path(path).with_suffix("") / "hessian_cartesian.npy")}
            evidence.update({str(Path(p).relative_to(ROOT)): value for p, value in replays.items()})
            record["execution_pcm_curvature"] = dict(
                mode=certificate["mode"], qualification_only=False, production_admission=True,
                scientific_model_changed=False, surface_grid_changed=False,
                qualification_sha256=seal_hash, source_geometry_sha256=api.digest(geometry),
                target_evidence_sha256=evidence,
                **{field: certificate[field] for field in ("adapter_sha256", "cpu_pcm_library_sha256", "gpu_pcm_library_sha256", "isolated_proof_sha256")},
            )
            record["execution_area_contraction"] = dict(
                adapter_sha256=api.digest(SCRIPT), replay_evidence_sha256=replays,
                atom_block=True, scientific_model_changed=False,
            )
            atomic_write_json(path, record)
            return ensure_projected_thermal(path, geometry)

        projected.ensure_projected_thermal = seal
        try:
            with solvent_memory(metrics):
                return original(*args, **kwargs)
        finally:
            projected.ensure_projected_thermal = prior
            if not proof.exists():
                atomic_write_json(proof, dict(metrics=metrics, status="incomplete_unqualified"))

    native.ensure_native_thermal = evaluate
    # Host curvature arrays can be larger than GPU memory. Bound parallelism by
    # the conservative measured pilot RAM projection, never a fixed worker cap.
    from redox.workflows import resources
    original_free = resources.free_slots

    def ram_suitable(hosts=None):
        slots = original_free(hosts)
        proofs = list((Path(api.read_plan(CAMPAIGN)["raw_root"]) / "thermal").glob("*/rest/solvent_memory_execution.json"))
        projected_bytes = [event["projected_conservative_host_peak_bytes"]
                           for path in proofs for event in json.loads(path.read_text()).get("metrics", {}).get("host_curvature", {}).get("calls", [])]
        if hosts == [socket.gethostname().split(".")[0]] and projected_bytes:
            import psutil
            count = int(psutil.virtual_memory().available * .75 // max(projected_bytes))
            return slots[:count]
        return slots

    resources.free_slots = ram_suitable
    return base


def coordinate(hours):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    import fcntl

    allowed_cpus = os.sched_getaffinity(0)
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    with (CAMPAIGN / "replay_coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.time() + 3600 * hours
        while time.time() < deadline:
            rows = [json.loads(path.read_text()) for path in (CAMPAIGN / "controls").glob("*/status.json")]
            if any(row["status"] == "failed_preserved" for row in rows):
                atomic_write_json(CAMPAIGN / "replay_coordinator.json", dict(status="replay_requires_review", controls=rows))
                return
            if len(rows) == 2 and all(row["status"] == "full_native_replay_passed" for row in rows):
                atomic_write_json(CAMPAIGN / "replay_qualification.json", dict(status="full_native_replays_passed", evidence_sha256=checked_replays()))
                # The reused coordinator allocates its own CPU from full allowed
                # affinity; restore it after our single-core waiting allocation.
                os.sched_setaffinity(0, allowed_cpus)
                runner().coordinate(CAMPAIGN, max(.01, (deadline-time.time())/3600))
                return
            atomic_write_json(CAMPAIGN / "replay_coordinator.json", dict(status="awaiting_full_native_replays", controls=rows, deadline=deadline))
            time.sleep(20)
        atomic_write_json(CAMPAIGN / "replay_coordinator.json", dict(status="finite_resource_window_ended"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "replay", "coordinate", "worker", "launch", "collect"))
    parser.add_argument("--mult", type=int, choices=(1, 2))
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN)
    parser.add_argument("--stage", choices=("energy", "thermal"), default="thermal")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if args.campaign.resolve() != CAMPAIGN or args.max_hours <= 0:
        parser.error("Use the frozen campaign and a positive resource window")
    if args.action == "prepare":
        prepare()
    elif args.action == "replay":
        if args.mult is None:
            parser.error("replay requires --mult")
        replay(args.mult)
    elif args.action == "coordinate":
        coordinate(args.max_hours)
    elif args.action == "collect":
        api.collect(CAMPAIGN)
    else:
        getattr(runner(), args.action)(CAMPAIGN, args.stage, args.pilot)
