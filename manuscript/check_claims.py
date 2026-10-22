"""Independently reconstruct the numerical statements in the frozen paper snapshot."""

import csv
import json
import math
from pathlib import Path


def audit():
    folder = Path(__file__).resolve().parent / "evidence/2026-10-21"
    load = lambda name: json.loads((folder / name).read_text())
    benchmark = load("benchmark.json")
    panel = load("benchmark_panel.json")
    rows = list(csv.DictReader((folder / "benchmark_points.csv").open()))
    errors = [float(r["calc_V"]) - float(r["exp_V"]) for r in rows]
    keys = [(r["id"], r["event"], r["doi"]) for r in rows]
    expected = {(r["id"], r["event"], r["doi"]): r for r in panel["targets"]}
    if len(set(keys)) != len(keys) or not set(keys) <= set(expected):
        raise ValueError("Duplicated or non-panel manuscript benchmark point")
    if any(abs(float(r["exp_V"]) - expected[key]["exp_V"]) > 1e-12 for r, key in zip(rows, keys)):
        raise ValueError("Manuscript experimental label differs from frozen panel")
    statistics = dict(MAE_V=math.fsum(abs(e) for e in errors) / len(errors),
                      RMSE_V=math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
                      bias_V=math.fsum(errors) / len(errors), maximum_absolute_error_V=max(abs(e) for e in errors))
    if any(abs(statistics[key] - benchmark["metrics"][key]) > 1e-12 for key in statistics):
        raise ValueError("Manuscript benchmark arithmetic mismatch")
    if (len(rows) != benchmark["admitted_events"] or len(expected) != benchmark["target_events"]
            or len({r["id"] for r in rows}) != benchmark["admitted_molecules"]
            or len({r["doi"] for r in rows}) != benchmark["original_studies"]):
        raise ValueError("Manuscript coverage mismatch")
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
    return dict(status="numerical_snapshot_reconstructed", benchmark=dict(admitted=len(rows), target=len(expected),
                molecules=len({r["id"] for r in rows}), studies=len({r["doi"] for r in rows}), **statistics),
                base_assisted_values=base_values, qualified_spacing_G=sum(r["G_eV"] is not None for r in cases.values()),
                spacing_graphs=len(cases), exact_route_entries=len(load("route_evidence.json")["candidates"]),
                source_energy_conversion_kJ_mol_per_eV=96.485,
                maximum_SI_conversion_difference_kJ_mol=max(conversion_differences),
                experimental_accuracy_inferred=False, submission_ready=False)


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
