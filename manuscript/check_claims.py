"""Independently reconstruct the numerical statements in the frozen paper snapshot."""

import csv
import json
import math
import statistics as descriptive
from pathlib import Path


def audit():
    root = Path(__file__).resolve().parents[1]
    pointer = json.loads((root / "manuscript/evidence/current.json").read_text())
    folder = (root / pointer["manifest"]).parent
    load = lambda name: json.loads((folder / name).read_text())
    benchmark = load("benchmark.json")
    panel = load("benchmark_panel.json")
    rows = list(csv.DictReader((folder / "benchmark_points.csv").open()))
    errors = [float(r["calc_V"]) - float(r["exp_V"]) for r in rows]
    keys = [(r["id"], r["event"], r["doi"]) for r in rows]
    targets = {(r["id"], r["event"], r["doi"]): r for r in panel["targets"]}
    if len(set(keys)) != len(keys) or not set(keys) <= set(targets):
        raise ValueError("Duplicated or non-panel manuscript benchmark point")
    if any(abs(float(r["exp_V"]) - targets[key]["exp_V"]) > 1e-12 for r, key in zip(rows, keys)):
        raise ValueError("Manuscript experimental label differs from frozen panel")
    statistics = dict(MAE_V=math.fsum(abs(e) for e in errors) / len(errors),
                      RMSE_V=math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
                      bias_V=math.fsum(errors) / len(errors), maximum_absolute_error_V=max(abs(e) for e in errors))
    if any(abs(statistics[key] - benchmark["metrics"][key]) > 1e-12 for key in statistics):
        raise ValueError("Manuscript benchmark arithmetic mismatch")
    if (len(rows) != benchmark["admitted_events"] or len(targets) != benchmark["target_events"]
            or len({r["id"] for r in rows}) != benchmark["admitted_molecules"]
            or len({r["doi"] for r in rows}) != benchmark["original_studies"]):
        raise ValueError("Manuscript coverage mismatch")
    bias = load("reference_bias.json")
    if bias["primary_potentials_modified"] or bias["external_validation_claimed"]:
        raise ValueError("Diagnostic offsets replaced the primary benchmark")
    cv_checks = 0
    for model in bias["cross_validation"]:
        residuals = []
        for index, row in enumerate(rows):
            label = model["holdout"]
            members = [j for j, other in enumerate(rows)
                       if (j != index if label == "endpoint" else other[label] != row[label])
                       and (model["model"] == "constant_offset"
                            or targets[keys[j]]["wave"] == targets[keys[index]]["wave"])]
            if not members:
                continue
            offset = math.fsum(errors[j] for j in members) / len(members)
            residuals.append(errors[index] - offset)
        checked = dict(MAE_V=math.fsum(abs(e) for e in residuals)/len(residuals),
                       RMSE_V=math.sqrt(math.fsum(e*e for e in residuals)/len(residuals)),
                       bias_V=math.fsum(residuals)/len(residuals))
        if any(abs(checked[k]-model["metrics"][k]) > 1e-12 for k in checked):
            raise ValueError("Grouped offset diagnostic does not reconstruct")
        cv_checks += 1
    for interval in bias["uncertainty"]:
        primary = next(r for r in benchmark["uncertainty"] if r["cluster"] == interval["cluster"])
        if interval["intervals"] != primary["intervals"]:
            raise ValueError("Printed bootstrap interval differs from primary evidence")
    base = load("base_assisted.json")
    base_values = 0
    for candidate in base["candidates"]:
        for name, modes in candidate["base_models"].items():
            for mode, value in modes.items():
                if value is None:
                    continue
                if abs(value - candidate["reference_deltaG_kJ_mol"] - base["corrections"][name][mode]) > 1e-7:
                    raise ValueError("Manuscript base-assisted cycle mismatch")
                base_values += 1
    polymer = load("polymer_validation.json")
    cases = {r["key"]: r for r in polymer["cases"]}
    conversion_differences = []
    for cycle in polymer["cycles"]:
        difference_eV = math.fsum(cases[key]["E_eV"] * weight for key, weight in cycle["coefficients"].items())
        # Existing source-frozen records use this documented five-digit factor.
        # Also record its difference from the SI conversion rather than hiding it.
        value = difference_eV * 96.485
        conversion_differences.append(abs(difference_eV * 96.48533212 - value))
        if abs(value - cycle["electronic_deltaDeltaE_kJ_mol"]) > 1e-6:
            raise ValueError("Manuscript electronic occupancy cycle mismatch")
        if any(cases[key]["G_eV"] is None for key in cycle["coefficients"]):
            if cycle["free_energy_deltaDeltaG_kJ_mol"] is not None:
                raise ValueError("Missing native thermochemistry was filled in")
        else:
            dg = math.fsum(cases[key]["G_eV"] * weight for key, weight in cycle["coefficients"].items()) * 96.485
            if abs(dg - cycle["free_energy_deltaDeltaG_kJ_mol"]) > 1e-6:
                raise ValueError("Manuscript native occupancy cycle mismatch")
            for cutoff, expected in cycle["rotor_cutoff_sensitivity_kJ_mol"].items():
                value = math.fsum(cases[key]["cutoff_G_eV"][cutoff] * weight for key, weight in cycle["coefficients"].items()) * 96.485
                if abs(value - expected) > 1e-6:
                    raise ValueError("Manuscript rotor-cutoff sensitivity mismatch")
            for terms in cycle["conditional_grafting_kJ_mol"].values():
                for second, first in (("second_j_given_i", "first_j"), ("second_i_given_j", "first_i")):
                    if abs(terms[second] - terms[first] - dg) > 1e-6:
                        raise ValueError("Conditional grafting difference does not close")
    routes = load("exact_routes.json")
    if routes["completed"] != routes["expected"] or sum(routes["counts"].values()) != routes["completed"]:
        raise ValueError("Exact route handoff coverage mismatch")
    acid = load("acid_base.json")
    for row in acid["rows"]:
        if row["reference"] is None:
            if row["external_capture_kJ_mol"] is not None:
                raise ValueError("Missing exact-base anchor filled in")
            continue
        expected = -.00831446261815324 * acid["temperature_K"] * math.log(10) * (row["reference"]["pKaH_MeCN"] - acid["acid_anchor"]["pKa_MeCN"])
        if abs(expected - row["external_capture_kJ_mol"]) > 1e-10:
            raise ValueError("External acid-capture conversion mismatch")
    concurrency = load("cpu_concurrency.json")
    if concurrency["run_count"] != len(concurrency["runs"]) or concurrency["serial_discrepancies"]:
        raise ValueError("Concurrency experiment count or independent comparison differs")
    for condition in concurrency["conditions"]:
        runs = [r for r in concurrency["runs"] if r["condition"] == condition["condition"]]
        if len(runs) != condition["repeats"]:
            raise ValueError("Concurrency repeat count differs")
        for key in ("launch_to_observed_completion_seconds", "workers_executing_chunks"):
            values = [r[key] for r in runs]
            expected = dict(median=descriptive.median(values), minimum=min(values), maximum=max(values))
            if expected != condition[key]:
                raise ValueError("Concurrency summary does not reconstruct from actual runs")
    fresh = load("fresh_checkout.json")
    if not fresh["figure_CSV_exact_match"] or not fresh["figure_RGB_pixels_exact_match"] or fresh["serial_discrepancies"]:
        raise ValueError("Fresh-checkout example did not reproduce the checked results")
    transfer = load("proton_transfer.json")
    components = {r["id"]: r["G_eV"] for r in transfer["states"]}
    references = load("proton_transfer_sources.json")
    for row, anchor in zip(transfer["comparisons"], references["comparisons"]):
        if row["base"] != anchor["base"]:
            raise ValueError("Relative-basicity source identity mismatch")
        key = row["base"]
        dg = (components[key+"_protonated"] + components["dbu_neutral"]
              - components[key+"_neutral"] - components["dbu_protonated"]) * 96.48533212331002
        external = .00831446261815324 * 298.15 * math.log(10) * (
            anchor["reference"]["pKaH_MeCN"]-anchor["source"]["pKaH_MeCN"])
        if abs(dg-row["DFT_deltaG_kJ_mol"]) > 1e-6 or abs(external-row["external_deltaG_kJ_mol"]) > 1e-9:
            raise ValueError("Relative-basicity comparison does not reconstruct")
    restart = load("quantum_restart.json")
    for control in restart["controls"]:
        if any(control["differences"][key] > tolerance for key, tolerance in restart["tolerances"].items()):
            raise ValueError("Quantum restart exceeds its predeclared tolerance")
    return dict(status="numerical_snapshot_reconstructed", benchmark=dict(admitted=len(rows), target=len(targets),
                molecules=len({r["id"] for r in rows}), studies=len({r["doi"] for r in rows}), **statistics),
                base_assisted_values=base_values, qualified_spacing_G=sum(r["G_eV"] is not None for r in cases.values()),
                spacing_graphs=len(cases), exact_route_entries=len(load("route_evidence.json")["candidates"]),
                exact_routes_checked=routes["completed"], route_decisions=routes["counts"],
                acid_base_scope=acid["status"], conditional_native_cycles=len(polymer["cycles"]),
                checked_concurrency_runs=len(concurrency["runs"]), fresh_checkout_designs=fresh["exact_designs"],
                source_energy_conversion_kJ_mol_per_eV=96.485,
                maximum_SI_conversion_difference_kJ_mol=max(conversion_differences),
                checked_offset_models=cv_checks, relative_basicity_controls=len(transfer["comparisons"]),
                quantum_restart_controls=len(restart["controls"]),
                experimental_accuracy_inferred=False, submission_ready=False)


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
