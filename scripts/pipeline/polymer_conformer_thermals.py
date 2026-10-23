#!/usr/bin/env python
"""Native thermochemistry and finite partitions for every declared DFT pose.

This adapter fills six missing thermochemical records in the existing phthalimide
two-pose spacing panel. It does not replace the wider CREST searches. The finite
partition functions are independent of the motif; this campaign is a control.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import time
from types import SimpleNamespace

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

SCRIPT = Path(__file__).resolve()
PARENT = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_cds_recovery"
HOST_BUDGET_GIB = 256  # Bounded 128-GiB curvature cache plus 128-GiB headroom.


def alignment_molecule(smiles, xyz):
    """Map declared connectivity to XYZ before symmetry-aware alignment."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from redox.build.geometry import _skeleton, geometry_stereo_status
    from redox.validation.protonation import hydrogen_topology_status

    if hydrogen_topology_status(smiles, xyz) != "ok" or geometry_stereo_status(smiles, xyz) not in ("ok", "not_applicable"):
        raise ValueError("Conformer graph/stereochemistry differs")
    expected = Chem.AddHs(Chem.MolFromSmiles(smiles))
    observed = Chem.MolFromXYZFile(str(xyz))
    rdDetermineBonds.DetermineConnectivity(observed, charge=0)
    mapping = _skeleton(observed).GetSubstructMatch(_skeleton(expected))
    if len(mapping) != expected.GetNumAtoms():
        raise ValueError("Incomplete full-atom geometry mapping")
    conformer = Chem.Conformer(expected.GetNumAtoms())
    conformer.Set3D(True)
    for i, j in enumerate(mapping):
        conformer.SetAtomPosition(i, observed.GetConformer().GetAtomPosition(j))
    expected.AddConformer(conformer)
    return Chem.RemoveHs(expected)


def collect(campaign):
    """Recheck native parents; keep absent G and finite-sampling limits explicit."""
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from redox.core.protocol import load_record
    from redox.qm.native_thermal import load_native_thermal, NATIVE_THERMAL
    from redox.properties.conformer_ensemble import cluster_minima, neighbor_cycle, summarize

    plan = api.read_plan(campaign)
    energy, invalid = api.admitted_records(campaign, "energy")
    thermal, thermal_invalid = api.admitted_records(campaign, "thermal")
    invalid += thermal_invalid
    graphs = []
    partition_protocol = None
    for case in plan["graph_cases"]:
        tasks = [task for task in plan["tasks"] if task["model_key"] == case["key"]]
        rows, molecules = [], []
        if any(task["key"] not in energy for task in tasks):
            raise ValueError("The frozen energy panel is incomplete")
        for task in tasks:
            key = task["key"]
            parent = energy[key]
            row = dict(key=key, E_eV=parent["E_smd_eV"], G_eV=None)
            molecules.append(alignment_molecule(task["smiles"], Path(parent["state_dir"]) / "opt.xyz"))
            if key in thermal:
                report = thermal[key]
                state = Path(report["state_dir"])
                native = load_native_thermal(state, 0, 1)
                active = load_record(state, 0, 1)
                if native is None or active is None:
                    raise ValueError("A thermal receipt lacks a currently admitted native parent")
                protocol = native["thermal_partition_protocol"]
                if partition_protocol is None:
                    partition_protocol = protocol
                if protocol != partition_protocol or native["identity"]["protocol"] != NATIVE_THERMAL:
                    raise ValueError("Conformer thermochemical protocols differ")
                if api.digest(state / "opt.xyz") != api.digest(Path(parent["state_dir"]) / "opt.xyz"):
                    raise ValueError("Thermal geometry differs from the frozen energy parent")
                value = active["e_smd_eV"] + native["g_thermal_total_eV"]
                if not math.isfinite(value) or abs(value - report["G_eV"]) > 1e-8 or abs(active["e_smd_eV"] - parent["E_smd_eV"]) > 1e-8:
                    raise ValueError("Conformer free energy does not reconstruct")
                row.update(G_eV=value, thermal_receipt=str(campaign / "thermal/done" / (key + ".json")))
            rows.append(row)
        matrix = np.zeros((len(rows), len(rows)))
        for i in range(len(rows)):
            for j in range(i):
                matrix[i, j] = matrix[j, i] = rdMolAlign.GetBestRMS(Chem.Mol(molecules[i]), molecules[j], maxMatches=10000)
        by_cutoff = {}
        for cutoff in (0.05, 0.1, 0.2):
            groups = cluster_minima(rows, matrix, rmsd_cutoff_A=cutoff)
            by_cutoff[str(cutoff)] = summarize(rows, groups, [task["key"] for task in tasks])
        graphs.append(dict(case=case, conformers=rows, heavy_atom_RMSD_A=matrix.tolist(), partitions_by_duplicate_RMSD_A=by_cutoff))
    lookup = {row["case"]["key"]: row for row in graphs}
    cycles = []
    for i, j in ((2, 3), (2, 4)):
        keys = ["occ_none", f"occ_{i}", f"occ_{j}", f"occ_{i}_{j}"]
        sensitivity = {}
        for cutoff in ("0.05", "0.1", "0.2"):
            values = [lookup[key]["partitions_by_duplicate_RMSD_A"][cutoff]["G_eV"] for key in keys]
            sensitivity[cutoff] = neighbor_cycle(*values) if all(value is not None for value in values) else None
        cycles.append(dict(sites=[i, j], deltaDeltaG_kJ_mol_by_duplicate_RMSD_A=sensitivity))
    complete = len(thermal) == len(plan["tasks"]) and not invalid
    current = dict(
        status="finite_two_pose_partitions_complete" if complete else "pending_native_conformer_G",
        expected_DFT_poses=len(plan["tasks"]), checked_native_G=len(thermal),
        invalid_receipts=invalid, cases=graphs, cycles=cycles,
        temperature_K=298.15, solvent="acetonitrile", backbone_configuration="all_same",
        partition_protocol=partition_protocol, sampling_convergence_proven=False,
        experimental_occupancy_weights=None, new_synthesis_candidates_admitted=0,
        scope="Two-pose phthalimide control only; no pure motif binding, atactic population or intersite ET rate",
        source_plan_sha256=api.digest(campaign / "plan.json"), updated=time.time(),
    )
    atomic_write_json(campaign / "current.json", current)
    return current


def prepare(campaign):
    """Keep all existing records; calculate only previously absent pose G values."""
    energy, invalid = api.admitted_records(PARENT, "energy")
    thermal, thermal_invalid = api.admitted_records(PARENT, "thermal")
    if invalid or thermal_invalid or len(energy) != 12 or len(thermal) != 6:
        raise ValueError("Expected checked 12-energy/6-native baseline")
    campaign.mkdir(parents=True, exist_ok=False)
    for stage in ("energy", "thermal"):
        for directory in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (campaign / stage / directory).mkdir(parents=True)
        for receipt in (PARENT / stage / "done").glob("*.json"):
            verify_state_checkpoint(receipt)
            shutil.copyfile(receipt, campaign / stage / "done" / receipt.name)
    shutil.copytree(PARENT / "controls", campaign / "controls")
    parent = api.read_plan(PARENT)
    sources = list((ROOT / "src/redox").rglob("*.py")) + [SCRIPT]
    sources += [ROOT / "scripts/pipeline" / name for name in (
        "run_polymer_spacing.py", "resume_bounded_polymer_thermal.py", "resume_polymer_thermal.py",
        "resume_polymer_cds.py", "run_pi_stacking.py",
    )]
    sources += [ROOT / "scripts/validation/diagnose_saved_hessian.py"]
    inputs = [PARENT / "plan.json"] + list((campaign / "energy/done").glob("*.json"))
    inputs += list((campaign / "thermal/done").glob("*.json")) + list((campaign / "controls").glob("*/status.json"))
    plan = dict(
        parent, raw_root=str(ROOT / "calcs/polymer_models" / campaign.name),
        producer_sha256={str(path): api.digest(path) for path in sources},
        input_sha256={str(path): api.digest(path) for path in inputs},
        parent_campaign=str(PARENT), previously_checked_native_G=6,
        requested_new_native_G=6, all_declared_poses_require_G=True,
        Hessians=True,
        host_budget_GiB_per_thermal_worker=HOST_BUDGET_GIB,
        scope="Finite two-pose baseline for comparison with broader conformer searches; sampling adequacy unresolved",
    )
    atomic_write_json(campaign / "plan.json", plan)
    collect(campaign)


def backend(campaign, quantum=False):
    """Use the checked native worker for all poses, with private run adapters."""
    if quantum:
        from redox.qm import projected_thermal as projected, sp
        from redox.qm.hessian_resume import annotate_finer_hessian

        bounded = api.load_module("ensemble_bounded_native", ROOT / "scripts/pipeline/resume_bounded_polymer_thermal.py")
        bounded.CAMPAIGN = campaign
        previous = projected.ensure_projected_thermal

        def seal(path, geometry):
            annotate_finer_hessian(path, geometry)
            return previous(path, geometry)

        projected.ensure_projected_thermal = seal
        try:
            base = bounded.configure().runner()
        finally:
            projected.ensure_projected_thermal = previous
        original_ensure = sp.ensure

        def seed_sp(gid, state_name, state, q, m, **kwargs):
            energy, invalid = api.admitted_records(campaign, "energy")
            if invalid or gid not in energy:
                raise ValueError("Missing checked electronic conformer")
            source = Path(energy[gid]["state_dir"])
            if api.digest(source / "opt.xyz") != api.digest(Path(state) / "opt.xyz"):
                raise ValueError("SP seed belongs to a different geometry")
            (Path(state) / "sp").mkdir(exist_ok=True)
            for path in (source / "sp").glob("*.json"):
                target = Path(state) / "sp" / path.name
                if not target.exists():
                    shutil.copyfile(path, target)
            return original_ensure(gid, state_name, state, q, m, **kwargs)

        sp.ensure = seed_sp
    else:
        base = api.load_module("ensemble_native_queue", ROOT / "scripts/pipeline/run_polymer_spacing.py")
    base.SCRIPT = SCRIPT

    def phase_tasks(folder, plan, stage):
        if stage != "thermal":
            raise ValueError("This campaign reuses its checked electronic poses")
        records, invalid = api.admitted_records(folder, "energy")
        if invalid or len(records) != len(plan["tasks"]):
            raise ValueError("All declared electronic parents must be checked")
        completed, thermal_invalid = api.admitted_records(folder, "thermal")
        if thermal_invalid:
            raise ValueError("Invalid completed thermochemical receipt")
        # The original backend chooses one energy leader per graph. This private
        # adapter indexes every pose separately while retaining its original task
        # identity in the actual energy/thermal records.
        return [dict(record["task"], model_key=key) for key, record in records.items() if key not in completed]

    base.phase_tasks = phase_tasks
    base.campaign_api = SimpleNamespace(
        read_plan=api.read_plan, verify_sources=api.verify_sources,
        selected_graph_records=lambda folder: api.admitted_records(folder, "energy"),
        calculate_native_G=api.calculate_native_G, digest=api.digest, collect=collect,
    )
    return base


def launch(campaign, pilot):
    import psutil
    from redox.workflows import resources

    original = resources.free_slots

    def memory_suitable(hosts=None):
        # One finite campaign owns this launch. Reserve capacity conservatively
        # for the whole submitted batch before any Hessian reaches its RAM peak.
        available = int(psutil.virtual_memory().available * 0.75 // (HOST_BUDGET_GIB * 1024**3))
        return original(hosts)[:available]

    resources.free_slots = memory_suitable
    try:
        backend(campaign).launch(campaign, "thermal", pilot)
    finally:
        resources.free_slots = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "collect", "worker", "launch"))
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--stage", choices=("thermal",), default="thermal")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    if args.action in ("prepare", "collect"):
        from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        if args.action == "prepare":
            prepare(campaign)
        else:
            result = collect(campaign)
            print(json.dumps({key: result[key] for key in ("status", "expected_DFT_poses", "checked_native_G")}))
    elif args.action == "launch":
        launch(campaign, args.pilot)
    else:
        from redox.qm.cds_finite import finer_cds_hessian
        with finer_cds_hessian(1e-5):
            backend(campaign, quantum=True).worker(campaign, "thermal", args.pilot)


if __name__ == "__main__":
    main()
