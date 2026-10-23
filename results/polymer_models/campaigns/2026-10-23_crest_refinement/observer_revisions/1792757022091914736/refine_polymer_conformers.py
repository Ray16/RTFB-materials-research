#!/usr/bin/env python
"""DFT re-ranking and native finite partitions for a frozen CREST-derived case.

Per-case energy and thermal pilots gate fan-out. Corrected conservative geometry
forces, unchanged SMD energies, independent physical checks and native Hessian
admission remain separate requirements. Missing states never receive a proxy G.
"""
from contextlib import contextmanager
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import socket
import subprocess
import time
from types import SimpleNamespace

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()


def finite_batches(records, distances, temperature=298.15):
    """Compare fixed nested prefixes, deduplicating within each prefix itself."""
    from redox.properties.conformer_ensemble import cluster_minima, summarize, batch_change
    import numpy as np

    result = {}
    for cutoff in (0.05, 0.1, 0.2):
        batches = {}
        for n in (10, 20):
            if len(records) < n:
                batches[str(n)] = dict(status="insufficient_declared_starts", G_eV=None)
                continue
            prefix = records[:n]
            matrix = np.asarray(distances)[:n, :n]
            groups = cluster_minima(prefix, matrix, rmsd_cutoff_A=cutoff)
            batches[str(n)] = summarize(prefix, groups, [row["key"] for row in prefix], temperature)
        batches["convergence"] = batch_change(batches["10"], batches["20"])
        result[str(cutoff)] = batches
    return result


def collect(campaign):
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from redox.qm.native_thermal import load_native_thermal, NATIVE_THERMAL
    from redox.core.protocol import load_record

    plan = api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    energy, invalid = api.admitted_records(campaign, "energy")
    thermal, thermal_invalid = api.admitted_records(campaign, "thermal")
    invalid += thermal_invalid
    geometry = api.load_module("refined_conformer_mapping", ROOT / "scripts/pipeline/polymer_conformer_thermals.py")
    records, molecules = [], []
    partition_protocol = None
    for task in plan["tasks"]:
        key = task["key"]
        row = dict(key=key, search_rank=task["search_rank"], E_eV=None, G_eV=None,
            G_by_rotor_cutoff_eV={str(cutoff): None for cutoff in (25, 50, 100)})
        if key in energy:
            parent = energy[key]
            if parent.get("conservative_total_force_check", {}).get("passed") is not True:
                raise ValueError("Energy parent lacks its independent total-energy force check")
            xyz = Path(parent["state_dir"]) / "opt.xyz"
            row["E_eV"] = parent["E_smd_eV"]
            molecules.append(geometry.alignment_molecule(task["smiles"], xyz))
            if key in thermal:
                report = thermal[key]
                state = Path(report["state_dir"])
                native = load_native_thermal(state, 0, 1)
                active = load_record(state, 0, 1)
                if native is None or active is None or api.digest(state / "opt.xyz") != api.digest(xyz):
                    raise ValueError("Native parent is absent or belongs to another geometry")
                protocol = native["thermal_partition_protocol"]
                if partition_protocol is None:
                    partition_protocol = protocol
                if protocol != partition_protocol or native["identity"]["protocol"] != NATIVE_THERMAL:
                    raise ValueError("Conformer thermochemistry definitions differ")
                value = active["e_smd_eV"] + native["g_thermal_total_eV"]
                if not math.isfinite(value) or abs(value - report["G_eV"]) > 1e-8 or abs(active["e_smd_eV"] - row["E_eV"]) > 1e-8:
                    raise ValueError("Native G does not reconstruct from its electronic/thermal parents")
                row["G_eV"] = value
                row["G_by_rotor_cutoff_eV"] = {str(cutoff): active["e_smd_eV"] + native["g_spin_partition_eV"] + native["nuclear_by_cutoff_eV"][str(cutoff)] for cutoff in (25, 50, 100)}
        records.append(row)
    complete_energy = len(energy) == len(plan["tasks"]) and not invalid
    batches = {}
    if complete_energy:
        matrix = np.zeros((len(records), len(records)))
        for i in range(len(records)):
            for j in range(i):
                matrix[i, j] = matrix[j, i] = rdMolAlign.GetBestRMS(Chem.Mol(molecules[i]), molecules[j], maxMatches=10000)
        for rotor in (25, 50, 100):
            rows = [dict(row, G_eV=row["G_by_rotor_cutoff_eV"][str(rotor)]) for row in records]
            batches[str(rotor)] = finite_batches(rows, matrix)
        rmsd = matrix.tolist()
    else:
        rmsd = None
    primary = batches.get("50", {}).get("0.1", {}).get("20", {})
    result = dict(status="finite_20_start_partition_complete" if primary.get("status") == "finite_set_partition_computed" else "pending_checked_conformer_components",
        case=plan["graph_cases"][0], expected_starts=len(plan["tasks"]), checked_energy=len(energy), checked_native_G=len(thermal),
        invalid_receipts=invalid, records=records, heavy_atom_RMSD_A=rmsd,
        batches_by_rotor_cutoff_cm=batches, primary_G_eV=primary.get("G_eV"),
        temperature_K=298.15, solvent="acetonitrile", partition_protocol=partition_protocol,
        primary_duplicate_RMSD_A=0.1, primary_rotor_cutoff_cm=50,
        global_ensemble_convergence_proven=False, experimental_population_weights=None,
        source_plan_sha256=api.digest(campaign / "plan.json"), updated=time.time(),
        scope="Finite conformer set for one exact occupancy/backbone configuration; not pure motif binding or an ET rate")
    atomic_write_json(campaign / "current.json", result)
    return result


@contextmanager
def geometry_protocol(campaign, task):
    """Run-local force refinement; no live shared producer file is changed."""
    from redox.qm import dft
    from redox.qm.stable_geometry_scf import stabilize

    module = api.load_module("refinement_conservative_CDS", ROOT / "scripts/pipeline/resume_polymer_forces_adaptive.py")
    original = dft._build_smd_mf
    raw = Path(api.read_plan(campaign)["raw_root"]) / "energy"
    raw.mkdir(parents=True, exist_ok=True)

    def build(*args, **kwargs):
        mf = original(*args, **kwargs)
        if (args[1] if len(args) > 1 else kwargs.get("xc")) == "r2scan":
            stabilize(mf, raw / (task["key"] + "_SCF.jsonl"))
        return mf

    dft._build_smd_mf = build
    try:
        with module.adaptive_forces(raw / (task["key"] + "_CDS")):
            yield module
    finally:
        dft._build_smd_mf = original


def calculate_energy(campaign, task, pilot=False):
    with geometry_protocol(campaign, task) as forces:
        api.calculate_energy(campaign, task, pilot=pilot)
        if not pilot:
            state = Path(api.read_plan(campaign)["raw_root"]) / "energy" / task["key"] / "rest"
            proof = forces.base.final_energy_check(state / "opt.xyz", state / "total_energy_force_check")
            path = campaign / "energy/tasks" / task["key"] / "status.json"
            record = json.loads(path.read_text())
            record["conservative_total_force_check"] = proof
            record["output_sha256"].update({str(p): api.digest(p) for p in state.rglob("*") if p.is_file() and p.suffix in (".json", ".npy", ".xyz")})
            directory = Path(str(state.parents[1] / (task["key"] + "_CDS")) + ".evaluations")
            record["output_sha256"].update({str(p): api.digest(p) for p in directory.rglob("*") if p.is_file()})
            atomic_write_json(path, record)


def backend(campaign, stage, quantum=False):
    if stage == "thermal" and quantum:
        helper = api.load_module("refinement_native_adapter", ROOT / "scripts/pipeline/polymer_conformer_thermals.py")
        base = helper.backend(campaign, quantum=True)
    else:
        base = api.load_module("refinement_phase_queue", ROOT / "scripts/pipeline/run_polymer_spacing.py")
    base.SCRIPT = SCRIPT

    def phase_tasks(folder, plan, phase):
        if phase == "energy":
            return [base.task_input(folder, task) for task in plan["tasks"]]
        records, invalid = api.admitted_records(folder, "energy")
        if invalid or len(records) != len(plan["tasks"]):
            return []
        completed, invalid = api.admitted_records(folder, "thermal")
        if invalid:
            raise ValueError("Invalid completed native receipt")
        return [dict(row["task"], model_key=key) for key, row in records.items() if key not in completed]

    base.phase_tasks = phase_tasks
    base.campaign_api = SimpleNamespace(read_plan=api.read_plan, verify_sources=api.verify_sources,
        selected_graph_records=lambda folder: api.admitted_records(folder, "energy"),
        calculate_energy=calculate_energy, calculate_native_G=api.calculate_native_G, digest=api.digest, collect=collect)
    return base


def launch(campaign, stage, pilot):
    from redox.workflows import resources
    import psutil

    original = resources.free_slots
    original_suitability = resources.gpu_suitability
    observations = []

    def verified_device(gpu):
        try:
            return original_suitability(gpu)
        except (OSError, subprocess.SubprocessError) as error:
            observations.append(dict(gpu=gpu, status="identity_temporarily_unverified", error=repr(error)))
            raise RuntimeError("GPU identity could not be verified; no task attempt made") from error

    def suitable(hosts=None):
        try:
            slots = original(hosts)
        except (OSError, subprocess.SubprocessError) as error:
            observations.append(dict(status="availability_temporarily_unverified", error=repr(error)))
            return []
        # Account for promised growth of already launched workers. Counting
        # available RAM alone can admit too many Hessians during their small-SCF
        # startup phase, before their bounded curvature cache is allocated.
        promises = host_promises(campaign.parent)
        available = psutil.virtual_memory().available
        budget = 256 if stage == "thermal" else 32
        count = additional_host_slots(available, promises, budget)
        atomic_write_json(campaign / stage / (socket.gethostname().split(".")[0] + "_RAM_admission.json"),
            dict(available_RAM_bytes=available, worker_budget_GiB=budget, promises=promises,
                maximum_additional_workers=count, updated=time.time()))
        return slots[:count]

    resources.free_slots = suitable
    resources.gpu_suitability = verified_device
    try:
        backend(campaign, stage).launch(campaign, stage, pilot)
        return observations
    finally:
        resources.free_slots = original
        resources.gpu_suitability = original_suitability


def additional_host_slots(available_bytes, promises, budget_GiB):
    remaining = sum(max(0, row["budget_bytes"] - row["observed_unique_bytes"]) for row in promises)
    return max(0, int((available_bytes * 0.75 - remaining) // (budget_GiB * 1024**3)))


def host_promises(root):
    """Keep unknown startup/ownership reservations; release only finished work."""
    import psutil

    host = socket.gethostname().split(".")[0]
    paths = list(root.glob("*/energy/launches/*.json")) + list(root.glob("*/thermal/launches/*.json"))
    baseline = ROOT / "results/polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2/thermal"
    paths += list((baseline / "launches").glob("*.json"))
    rows = []
    for path in paths:
        launch = json.loads(path.read_text())
        if launch["host"] != host:
            continue
        stage = launch["stage"]
        for worker in launch["workers"]:
            status_path = path.parents[1] / "workers" / f"{host}_gpu{worker['gpu']}_{launch['role']}.json"
            status = json.loads(status_path.read_text()) if status_path.exists() else {}
            if "finished" in status:
                continue
            unique = 0
            try:
                proc = Path("/proc") / str(status["pid"])
                if (proc / "stat").read_text().rsplit(")", 1)[1].split()[19] == status["process_start"]:
                    unique = psutil.Process(status["pid"]).memory_full_info().uss
            except (OSError, KeyError, psutil.Error):
                pass  # Missing/unreachable/reused handles do not release a promise.
            rows.append(dict(owner_record=str(status_path), budget_bytes=(256 if stage == "thermal" else 32) * 1024**3,
                observed_unique_bytes=unique))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("worker", "launch", "collect"))
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--stage", choices=("energy", "thermal"), default="energy")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    if args.action == "collect":
        from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        result = collect(campaign)
        print(json.dumps({key: result[key] for key in ("status", "checked_energy", "checked_native_G")}))
    elif args.action == "launch":
        host = socket.gethostname().split(".")[0]
        role = "pilot" if args.pilot else "queue"
        path = campaign / args.stage / f"launch_observer_{host}_{role}.json"
        with path.with_suffix(".lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            report = dict(host=host, pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                status="checking_resources_and_submitting", stage=args.stage, role=role, started=time.time())
            atomic_write_json(path, report)
            try:
                report["resource_observations"] = launch(campaign, args.stage, args.pilot)
                report["status"] = "launch_observation_completed_outputs_preserved"
            except BaseException as error:
                report.update(status="launch_observation_failed_requires_review", error=repr(error))
                raise
            finally:
                report["finished"] = time.time()
                atomic_write_json(path, report)
    elif args.stage == "thermal":
        from redox.qm.cds_finite import finer_cds_hessian
        with finer_cds_hessian(1e-5):
            backend(campaign, args.stage, quantum=True).worker(campaign, args.stage, args.pilot)
    else:
        backend(campaign, args.stage).worker(campaign, args.stage, args.pilot)


if __name__ == "__main__":
    main()
