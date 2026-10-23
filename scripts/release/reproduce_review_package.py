#!/usr/bin/env python3
"""Offline primary-result replay using only Python's standard library.

All original absolute paths are descriptive provenance. This program opens only
files inside the supplied package and does not invoke the cluster or QM tools.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path


def local(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Package path escapes root")
    return path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def replay(root):
    root = Path(root).resolve()
    package = json.loads((root / "package_manifest.json").read_text())
    for relative, expected in package["files_sha256"].items():
        if digest(local(root, relative)) != expected:
            raise ValueError("Package artifact differs: " + relative)
    spec = importlib.util.spec_from_file_location("portable_cached_reader", local(root, "code/examples/cached_thermochemistry/reproduce.py"))
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    # Reuse the existing independently tested geometry/state/balance reader.
    cached = reader.reconstruct(root / "data")
    data = json.loads((root / "data/manifest.json").read_text())
    state = {r["id"]: r["checked_G_eV"] for r in data["states"]}
    aliases = data["aliases"]
    g = lambda label: state[aliases[label]]
    epoch = local(root, package["evidence_epoch"])
    load = lambda name: json.loads((epoch / name).read_text())
    reference = g("benchmark:ferrocene:ox") - g("benchmark:ferrocene:neu")
    errors = []
    for endpoint in load("benchmark.json")["outliers"]:
        left, right = endpoint["event"].split("->")
        value = g(f"benchmark:{endpoint['id']}:{left}") - g(f"benchmark:{endpoint['id']}:{right}") - reference
        if abs(value-endpoint["calc_V"]) > 1e-8:
            raise ValueError("Raw-state benchmark endpoint does not reconstruct")
        errors.append(value-endpoint["exp_V"])
    metrics = dict(MAE_V=math.fsum(abs(v) for v in errors)/len(errors),
        RMSE_V=math.sqrt(math.fsum(v*v for v in errors)/len(errors)), bias_V=math.fsum(errors)/len(errors))
    if any(abs(value-load("benchmark.json")["metrics"][key]) > 1e-8 for key, value in metrics.items()):
        raise ValueError("Raw-state benchmark metrics differ")
    finite = load("finite_ensemble.json")
    by_key = {r["id"]: r for r in data["states"]}
    conversion = 96.48533212331002
    rt = 8.31446261815324e-3*298.15/conversion
    partitions = {}
    for pair in finite["finite_pairs"]:
        cutoffs = {}
        for cutoff in ("25", "50", "100"):
            values = []
            for key in pair["conformer_keys"]:
                row = by_key[aliases["polymer:"+key]]
                part = json.loads(local(root / "data", row["files"]["partition"]).read_text())
                energy = json.loads(local(root / "data", row["files"]["energy"]).read_text())
                values.append(energy["e_smd_eV"] + part["nuclear_by_cutoff_eV"][cutoff] + part["g_spin_partition_eV"])
            minimum = min(values)
            value = minimum - rt*math.log(math.fsum(math.exp(-(v-minimum)/rt) for v in values))
            if abs(value-pair["independent_G_eV_by_rotor_cutoff_cm"][cutoff]) > 1e-8:
                raise ValueError("Portable finite partition differs")
            cutoffs[cutoff] = value
        partitions[pair["graph_key"]] = cutoffs
    cycles = []
    for cycle in finite["cycles"]:
        i, j = cycle["sites"]
        values = {}
        for cutoff in ("25", "50", "100"):
            value = math.fsum((partitions[f"occ_{i}_{j}"][cutoff], partitions["occ_none"][cutoff],
                -partitions[f"occ_{i}"][cutoff], -partitions[f"occ_{j}"][cutoff]))*conversion
            if abs(value-cycle["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"][cutoff]) > 1e-6:
                raise ValueError("Portable polymer cycle differs")
            values[cutoff] = value
        cycles.append(dict(sites=cycle["sites"], deltaDeltaG_kJ_mol=values))
    relative = load("proton_transfer.json")
    for row in relative["comparisons"]:
        base = row["base"]
        value = (g("relative_basicity:" + base + "_protonated") + g("base:dbu_neutral")
            - g("relative_basicity:" + base + "_neutral") - g("base:dbu_protonated"))*conversion
        if abs(value-row["DFT_deltaG_kJ_mol"]) > 1e-6:
            raise ValueError("Portable relative-basicity control differs")
    spec = importlib.util.spec_from_file_location("frozen_paper_claims", local(root, "code/manuscript/check_claims.py"))
    claims = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(claims)
    paper = claims.audit()
    return dict(status="primary_cached_results_reproduced_offline", raw_states=len(state), benchmark_endpoints=len(errors),
        benchmark_metrics=metrics, reference_reactions=sum(r["name"].startswith("reference:") for r in data["reactions"]),
        base_assisted_values=sum(not r["name"].startswith("reference:") for r in data["reactions"]),
        finite_graph_partitions=len(partitions), finite_cycles=cycles,
        relative_basicity_controls=len(relative["comparisons"]), frozen_manuscript_claims_replayed=True,
        all_axis_eligible=paper["decision_followup"]["all_axis_eligible"],
        all_reaction_sums_atom_charge_balanced=True, PMDI_first_reduction_V_vs_Fc=cached["PMDI_first_reduction_V_vs_Fc"],
        original_cluster_paths_opened=False, quantum_calculations_run=False, nuclear_partition_reintegrated=False,
        independent_physical_QC_reestablished=False, experimental_accuracy_inferred=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = replay(args.package)
    if args.out:
        with args.out.open("x") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
