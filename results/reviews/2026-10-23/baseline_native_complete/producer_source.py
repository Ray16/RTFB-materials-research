#!/usr/bin/env python
"""Independently reconstruct available baseline native partitions and pair sums.

The thermodynamic equations and Hessian projection come from the independent
audit reader, not production thermochemistry. Geometry mapping is reused. This
does not repeat the QM Hessian or establish conformer/experimental convergence.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from redox.core.common import ROOT
from redox.core.protocol import load_record, record_path
from redox.qm.native_thermal import cache_path, load_native_thermal
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    args.out.mkdir(parents=True, exist_ok=False)
    plan = api.read_plan(args.campaign)
    snapshot = args.out / "current_snapshot.json"
    snapshot.write_bytes((args.campaign / "current.json").read_bytes())
    current = json.loads(snapshot.read_text())
    audit_source = ROOT / "scripts/validation/audit/native_partition.py"
    geometry_source = ROOT / "scripts/pipeline/polymer_conformer_thermals.py"
    audit = api.load_module("independent_baseline_nuclear_partition", audit_source)
    geometry = api.load_module("baseline_geometry_mapping", geometry_source)
    sources = {str(p): api.digest(p) for p in (Path(__file__).resolve(), args.campaign / "plan.json", snapshot, audit_source, geometry_source)}
    published = {row["key"]: row for case in current["cases"] for row in case["conformers"] if row["G_eV"] is not None}
    checked = {}
    for task in plan["tasks"]:
        key = task["key"]
        if key not in published:
            continue
        receipt = args.campaign / "thermal/done" / (key + ".json")
        report = verify_state_checkpoint(receipt)
        state = Path(report["state_dir"])
        native = load_native_thermal(state, 0, 1)
        energy = load_record(state, 0, 1)
        surface, _ = cache_path(state, 0, 1)
        partition = audit.read_partition(surface, state / "opt.xyz")
        if native is None or energy is None or partition is None:
            raise ValueError("A published G lacks an admitted/independently reconstructed parent")
        value = energy["e_smd_eV"] + partition["g_thermal_total_eV"]
        if abs(value - published[key]["G_eV"]) > 1e-8 or abs(value - report["G_eV"]) > 1e-8:
            raise ValueError("Published G differs from independent reconstruction")
        sources.update(partition["evidence_sha256"])
        for path in (receipt, Path(json.loads(receipt.read_text())["status_path"]), record_path(state, 0, 1)):
            sources[str(path)] = api.digest(path)
        checked[key] = dict(key=key, G_eV=value, E_eV=energy["e_smd_eV"], state_dir=str(state),
            G_by_rotor_cutoff_eV={str(cutoff): energy["e_smd_eV"] + partition["nuclear_by_cutoff_eV"][str(cutoff)] + native["g_spin_partition_eV"] for cutoff in (25, 50, 100)},
            min_frequency_cm=min(native["frequencies_cm"]), gradient_max_Ha_Bohr=native["gradient_max"],
            maximum_Hessian_asymmetry=native["maximum_asymmetry"],
            maximum_sampled_column_error=max(row["relative_error"] for row in native["derivative_checks"]))
    finite_pairs = []
    for case in current["cases"]:
        keys = [row["key"] for row in case["conformers"]]
        if len(keys) != 2:
            raise ValueError("This baseline audit requires exactly two declared poses per graph")
        if any(key not in checked for key in keys):
            continue
        rows = [checked[key] for key in keys]
        molecules = [geometry.alignment_molecule(case["case"]["smiles"], Path(row["state_dir"]) / "opt.xyz") for row in rows]
        rmsd = rdMolAlign.GetBestRMS(Chem.Mol(molecules[0]), molecules[1], maxMatches=10000)
        values = [row["G_eV"] for row in rows]
        kt = 1.380649e-23 / 1.602176634e-19 * 298.15
        by_cutoff = {}
        for cutoff in (0.05, 0.1, 0.2):
            duplicate = rmsd <= cutoff and abs(rows[0]["E_eV"] - rows[1]["E_eV"]) * 96.48533212331001 <= 0.2
            minimum = min(values)
            weights = [math.exp(-(value - minimum) / kt) for value in values]
            G = minimum if duplicate else minimum - kt * math.log(math.fsum(weights))
            expected = case["partitions_by_duplicate_RMSD_A"][str(cutoff)]
            if expected["G_eV"] is None or abs(G - expected["G_eV"]) > 1e-8:
                raise ValueError("Finite pair partition differs from independent log-sum")
            by_cutoff[str(cutoff)] = G
        finite_pairs.append(dict(graph_key=case["case"]["key"], conformer_keys=keys, heavy_atom_RMSD_A=rmsd,
            independent_G_eV_by_duplicate_RMSD_A=by_cutoff,
            configurational_lowering_kJ_mol=(min(values) - by_cutoff["0.1"]) * 96.48533212331001,
            conditional_two_pose_populations=[value / math.fsum(weights) for value in weights] if not duplicate else None))
        rotor_values = {}
        for rotor in (25, 50, 100):
            free = [row["G_by_rotor_cutoff_eV"][str(rotor)] for row in rows]
            minimum = min(free)
            rotor_values[str(rotor)] = minimum - kt * math.log(math.fsum(math.exp(-(value - minimum) / kt) for value in free))
        # All six baseline pairs are separated by >3 A, independently replayed
        # above. Do not silently generalize this all-distinct shortcut.
        if rmsd <= 0.2:
            raise ValueError("Rotor audit requires explicit duplicate handling for this pair")
        finite_pairs[-1]["independent_G_eV_by_rotor_cutoff_cm"] = rotor_values
    lookup = {row["graph_key"]: row for row in finite_pairs}
    cycles = []
    terms_path = Path(plan["parent_campaign"]) / "reaction_terms.json"
    terms = json.loads(terms_path.read_text())
    api.verify_sources(terms["source_sha256"])
    sources.update(terms["source_sha256"])
    sources[str(terms_path)] = api.digest(terms_path)
    for i, j in ((2, 3), (2, 4)):
        keys = ["occ_none", f"occ_{i}", f"occ_{j}", f"occ_{i}_{j}"]
        if any(key not in lookup for key in keys):
            continue
        by_rotor = {}
        for rotor in ("25", "50", "100"):
            values = [lookup[key]["independent_G_eV_by_rotor_cutoff_cm"][rotor] for key in keys]
            by_rotor[rotor] = math.fsum((values[3], values[0], -values[1], -values[2])) * 96.48533212331001
        published_cycle = next(row for row in current["cycles"] if row["sites"] == [i, j])
        if abs(by_rotor["50"] - published_cycle["deltaDeltaG_kJ_mol_by_duplicate_RMSD_A"]["0.1"]) > 1e-6:
            raise ValueError("Conditional cycle differs from independently reconstructed partitions")
        empty, first, second, both = [lookup[key]["independent_G_eV_by_rotor_cutoff_cm"]["50"] for key in keys]
        routes = {}
        for model, constant in terms["constant_eV"].items():
            routes[model] = dict(
                graft_i_on_empty_kJ_mol=(first - empty + constant) * 96.48533212331001,
                graft_j_on_empty_kJ_mol=(second - empty + constant) * 96.48533212331001,
                graft_j_given_i_kJ_mol=(both - first + constant) * 96.48533212331001,
                graft_i_given_j_kJ_mol=(both - second + constant) * 96.48533212331001)
        cycles.append(dict(sites=[i, j], deltaDeltaG_kJ_mol_by_rotor_cutoff_cm=by_rotor,
            conditional_DBU_grafting=routes, scope="Two-pose neutral phthalimide control; frozen reagent models, no kinetics or material admission"))
    if any(api.digest(path) != digest for path, digest in sources.items()):
        raise ValueError("Independent review input changed")
    report = dict(status="available_native_partitions_and_complete_pairs_independently_reconstructed", checked_at=time.time(),
        checked_native_G=len(checked), expected_native_G=len(plan["tasks"]),
        pending_keys=[task["key"] for task in plan["tasks"] if task["key"] not in checked],
        states=list(checked.values()), finite_pairs=finite_pairs, cycles=cycles, source_sha256=sources,
        scope="Numerical native partition/energy reconstruction and finite-pair sum; no new QM, global sampling or experimental validation")
    (args.out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.out / "producer_source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({key: report[key] for key in ("status", "checked_native_G", "expected_native_G", "cycles")}))


if __name__ == "__main__":
    main()
