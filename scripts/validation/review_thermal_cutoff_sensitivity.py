#!/usr/bin/env python
"""Independently check cutoff-sensitivity outputs from full archived G sums."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def review(path):
    result = json.loads(path.read_text())
    api.verify_sources(result["source_sha256"])
    if (result["quantum_calculations_run"] or result["primary_results_modified"]
            or result["experimental_accuracy_inferred"] or result["candidate_prediction_intervals_established"]
            or result["settings"]["primary_cutoff_cm"] != 50
            or result["settings"]["cutoffs_cm"] != [25, 50, 100]):
        raise ValueError("Thermal-setting sensitivity acquired an unsupported scope")
    root = Path(result["source_package"])
    package = json.loads((root / "package_manifest.json").read_text())
    if api.digest(root / "package_manifest.json") != result["source_package_manifest_sha256"]:
        raise ValueError("Source package identity changed")
    for name, h in package["files_sha256"].items():
        p = (root / name).resolve()
        if not p.is_relative_to(root.resolve()) or api.digest(p) != h:
            raise ValueError("Source package artifact changed or escapes its root")
    data = json.loads((root / "data/manifest.json").read_text())
    values = {}
    for row in data["states"]:
        if row["partition_kind"] == "monatomic":
            values[row["id"]] = {k: row["checked_G_eV"] for k in ("25", "50", "100")}
        else:
            part = json.loads((root / "data" / row["files"]["partition"]).read_text())
            values[row["id"]] = {k: row["checked_G_eV"] + part["nuclear_by_cutoff_eV"][k]
                - part["g_nuclear_qrrho_eV"] for k in ("25", "50", "100")}
    aliases = data["aliases"]
    endpoint_checks = []
    for row in result["benchmark"]["endpoints"]:
        left, right = row["event"].split("->")
        for k in ("25", "50", "100"):
            get = lambda label: values[aliases[label]][k]
            e = get(f"benchmark:{row['id']}:{left}") - get(f"benchmark:{row['id']}:{right}") \
                - get("benchmark:ferrocene:ox") + get("benchmark:ferrocene:neu")
            if abs(e-row["E_by_cutoff_V"][k]) > 1e-8:
                raise ValueError("Full-state potential differs from cutoff result")
            endpoint_checks.append(abs(e-row["E_by_cutoff_V"][k]))
        span = 1000*(max(row["E_by_cutoff_V"].values())-min(row["E_by_cutoff_V"].values()))
        if abs(span-row["range_mV"]) > 1e-8:
            raise ValueError("Endpoint cutoff range differs")
    if abs(max(r["range_mV"] for r in result["benchmark"]["endpoints"])-result["benchmark"]["max_endpoint_range_mV"]) > 1e-9:
        raise ValueError("Maximum benchmark cutoff range differs")
    recipes = {r["name"]: r for r in data["reactions"]}
    reaction_checks = []
    conversion = 1.602176634e-19*6.02214076e23/1000
    for row in result["grafting"]["reactions"]:
        source = recipes.pop(row["name"])
        correction = source["delta_n_solutes"]*.00831446261815324*data["temperature_K"]*math.log(.08205736608096*data["temperature_K"])
        for k in ("25", "50", "100"):
            dg = math.fsum(nu*values[state][k] for state, nu in source["coefficients"].items())*conversion+correction
            error = abs(dg-row["SI_deltaG_by_cutoff_kJ_mol"][k])
            if error > 1e-6:
                raise ValueError("Full-state reaction differs from cutoff result")
            reaction_checks.append(error)
        span = max(row["SI_deltaG_by_cutoff_kJ_mol"].values())-min(row["SI_deltaG_by_cutoff_kJ_mol"].values())
        if abs(span-row["range_kJ_mol"]) > 1e-8:
            raise ValueError("Reaction cutoff range differs")
    if recipes:
        raise ValueError("Cutoff reaction coverage is incomplete")
    crossings = [r["name"] for r in result["grafting"]["reactions"] if
        min(r["SI_deltaG_by_cutoff_kJ_mol"].values()) < 0 < max(r["SI_deltaG_by_cutoff_kJ_mol"].values())]
    if crossings != result["grafting"]["zero_crossings"]:
        raise ValueError("Reaction sign-crossing inventory differs")
    return dict(status="thermal_cutoff_full_state_reconstruction_checked", checked=time.time(),
        result_sha256=api.digest(path), potential_sums=len(endpoint_checks), reaction_sums=len(reaction_checks),
        maximum_potential_reconstruction_error_V=max(endpoint_checks),
        maximum_reaction_reconstruction_error_kJ_mol=max(reaction_checks),
        source_sha256={str(path): api.digest(path), str(Path(__file__).resolve()): api.digest(__file__)},
        scope="Independent sums and scenario ranges; no fresh physical QC or experimental prediction interval")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--result", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise FileExistsError("Preserve prior independent review")
    value = review(args.result.resolve())
    atomic_write_json(args.out, value)
    print(json.dumps(value))


if __name__ == "__main__":
    main()
