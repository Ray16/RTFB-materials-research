#!/usr/bin/env python
"""Linked cutoff sensitivity from checked archived nuclear partitions.

Primary values stay at 50 cm^-1. Electronic energies, geometries, Hessians,
spin factors, standard states and all source results remain unchanged.
"""
import argparse
from collections import defaultdict
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

CUTOFFS = (25, 50, 100)
CONVERSION = 1.602176634e-19 * 6.02214076e23 / 1000


def metrics(errors):
    return dict(MAE_V=math.fsum(abs(x) for x in errors) / len(errors),
                RMSE_V=math.sqrt(math.fsum(x*x for x in errors) / len(errors)),
                bias_V=math.fsum(errors) / len(errors))


def state_shifts(root, data, local):
    """Return only nuclear deltas; unchanged spin/translation conventions cancel."""
    shifts = {}
    for row in data["states"]:
        if row["partition_kind"] == "monatomic":
            shifts[row["id"]] = {str(k): 0. for k in CUTOFFS}
            continue
        part = json.loads(local(root / "data", row["files"]["partition"]).read_text())
        values = part["nuclear_by_cutoff_eV"]
        if (not all(math.isfinite(values[str(k)]) for k in CUTOFFS)
                or abs(values["50"] - part["g_nuclear_qrrho_eV"]) > 1e-8
                or abs(part["g_nuclear_qrrho_eV"] + part["g_spin_partition_eV"] - part["g_thermal_total_eV"]) > 1e-8):
            raise ValueError("Primary nuclear/spin partition does not reconstruct")
        shifts[row["id"]] = {str(k): values[str(k)] - values["50"] for k in CUTOFFS}
    return shifts


def analyze(root):
    replay_path = root / "code/scripts/release/reproduce_review_package.py"
    spec = importlib.util.spec_from_file_location("checked_cutoff_source_package", replay_path)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    replay = reader.replay(root)  # Checks archived bytes, state identities and primary arithmetic.
    package = json.loads((root / "package_manifest.json").read_text())
    data = json.loads((root / "data/manifest.json").read_text())
    epoch = reader.local(root, package["evidence_epoch"])
    benchmark = json.loads((epoch / "benchmark.json").read_text())
    shifts = state_shifts(root, data, reader.local)
    states = {r["id"]: r["checked_G_eV"] for r in data["states"]}
    aliases = data["aliases"]
    delta = lambda label, k: shifts[aliases[label]][str(k)]
    fc = {str(k): states[aliases["benchmark:ferrocene:ox"]] - states[aliases["benchmark:ferrocene:neu"]]
          + delta("benchmark:ferrocene:ox", k) - delta("benchmark:ferrocene:neu", k) for k in CUTOFFS}
    endpoints, by_molecule = [], defaultdict(dict)
    for row in benchmark["outliers"]:
        left, right = row["event"].split("->")
        wave = 1 if "neu" in (left, right) else 2
        primary = states[aliases[f"benchmark:{row['id']}:{left}"]] - states[aliases[f"benchmark:{row['id']}:{right}"]] - fc["50"]
        if abs(primary - row["calc_V"]) > 1e-8:
            raise ValueError("Primary benchmark potential differs from the checked state sums")
        values = {str(k): row["calc_V"] + delta(f"benchmark:{row['id']}:{left}", k)
            - delta(f"benchmark:{row['id']}:{right}", k)
            - delta("benchmark:ferrocene:ox", k) + delta("benchmark:ferrocene:neu", k) for k in CUTOFFS}
        entry = dict(id=row["id"], event=row["event"], wave=wave, exp_V=row["exp_V"],
            E_by_cutoff_V=values, range_mV=1000*(max(values.values())-min(values.values())),
            shift_from_50_mV={k: 1000*(v-values["50"]) for k, v in values.items()})
        endpoints.append(entry)
        by_molecule[row["id"]][wave] = values
    stats = {str(k): metrics([r["E_by_cutoff_V"][str(k)]-r["exp_V"] for r in endpoints]) for k in CUTOFFS}
    if any(abs(stats["50"][k] - benchmark["metrics"][k]) > 1e-9 for k in stats["50"]):
        raise ValueError("Primary benchmark statistics changed")
    spacings = [dict(id=key, spacing_by_cutoff_V={str(k): waves[1][str(k)]-waves[2][str(k)] for k in CUTOFFS})
                for key, waves in by_molecule.items() if set(waves) == {1, 2}]
    reactions = []
    for row in data["reactions"]:
        coefficients = row["coefficients"]
        primary_eV = math.fsum(states[key]*nu for key, nu in coefficients.items())
        standard = row["delta_n_solutes"]*.00831446261815324*data["temperature_K"]*math.log(.08205736608096*data["temperature_K"])
        source = primary_eV*data["source_energy_conversion_kJ_mol_per_eV"] + standard
        if abs(source-row["source_kJ_mol"]) > 1e-6:
            raise ValueError("Primary reaction sum differs from the checked source")
        values = {str(k): (primary_eV+math.fsum(shifts[key][str(k)]*nu for key, nu in coefficients.items()))*CONVERSION+standard for k in CUTOFFS}
        reactions.append(dict(name=row["name"], source_primary_kJ_mol=row["source_kJ_mol"],
            SI_deltaG_by_cutoff_kJ_mol=values, range_kJ_mol=max(values.values())-min(values.values()),
            signs_cross_zero=min(values.values()) < 0 < max(values.values())))
    references = {r["name"].split(":", 1)[1]: r for r in reactions if r["name"].startswith("reference:")}
    base_offsets = {}
    for base in ("dbu", "dipea"):
        for model in ("free_ions_kJ_mol", "contact_pair_kJ_mol"):
            rows = [r for r in reactions if r["name"].startswith(base+":"+model+":")]
            offsets = {}
            for k in CUTOFFS:
                differences = [r["SI_deltaG_by_cutoff_kJ_mol"][str(k)]
                    - references[r["name"].split(":", 2)[2]]["SI_deltaG_by_cutoff_kJ_mol"][str(k)] for r in rows]
                if max(differences)-min(differences) > 1e-6:
                    raise ValueError("Shared base offset was not constant at fixed cutoff")
                offsets[str(k)] = math.fsum(differences)/len(differences)
            base_offsets[base+":"+model] = offsets
    return dict(status="checked_nuclear_cutoff_sensitivity_reconstructed", checked=time.time(),
        settings=dict(cutoffs_cm=list(CUTOFFS), primary_cutoff_cm=50, temperature_K=data["temperature_K"],
            SI_kJ_mol_per_eV=CONVERSION, electronic_energies_changed=False, geometries_changed=False,
            Hessians_changed=False, spin_factors_changed=False),
        benchmark=dict(endpoints=endpoints, metrics_by_cutoff=stats, Fc_absolute_V_by_cutoff=fc,
            max_endpoint_range_mV=max(r["range_mV"] for r in endpoints),
            max_abs_shift_from_50_mV=max(abs(v) for r in endpoints for v in r["shift_from_50_mV"].values()),
            paired_spacings=spacings),
        grafting=dict(reactions=reactions, fixed_base_offsets_by_cutoff_kJ_mol=base_offsets,
            max_range_kJ_mol=max(r["range_kJ_mol"] for r in reactions),
            zero_crossings=[r["name"] for r in reactions if r["signs_cross_zero"]]),
        source_package=str(root), source_package_manifest_sha256=api.digest(root/"package_manifest.json"),
        source_replay_status=replay["status"], source_sha256={str(p): api.digest(p) for p in
            (Path(__file__).resolve(), replay_path, root/"package_manifest.json", root/"data/manifest.json", epoch/"benchmark.json")},
        quantum_calculations_run=False, primary_results_modified=False,
        experimental_accuracy_inferred=False, candidate_prediction_intervals_established=False,
        interpretation="Deterministic thermal-setting sensitivity only; not a sampling, electronic-structure or experimental uncertainty interval")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--package", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise FileExistsError("Preserve existing sensitivity output")
    result = analyze(args.package.resolve())
    args.out.mkdir(parents=True)
    atomic_write_json(args.out/"summary.json", result)
    print(json.dumps(dict(status=result["status"], endpoints=len(result["benchmark"]["endpoints"]),
        max_endpoint_range_mV=result["benchmark"]["max_endpoint_range_mV"],
        reactions=len(result["grafting"]["reactions"]), max_reaction_range_kJ_mol=result["grafting"]["max_range_kJ_mol"],
        sign_crossings=result["grafting"]["zero_crossings"])))


if __name__ == "__main__":
    main()
