#!/usr/bin/env python
"""Show independently checked occupancy cycles and conditional DBU grafting G."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, default=ROOT / "results/polymer_models/workflow_validation.json")
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/polymer_models")
    args = parser.parse_args()
    source = json.loads(args.review.read_text())
    if source["qualified_native_G"] != source["expected_graphs"] or source["qualified_spacing_conformers"] != source["expected_spacing_conformers"]:
        raise ValueError("All selected native partitions and electronic conformers must qualify")
    if any(row["native_QC"] != "passed" or row["geometry_QC"] != "passed" for row in source["cases"]):
        raise ValueError("Unqualified molecular field in plotted cohort")
    cases = {row["key"]: row for row in source["cases"]}
    rows = source["cycles"]
    for row in rows:
        value = sum(cases[key]["G_eV"] * n for key, n in row["coefficients"].items()) * 96.485
        if abs(value - row["free_energy_deltaDeltaG_kJ_mol"]) > 1e-6 or not row["atoms_and_charge_balanced"]:
            raise ValueError("Plotted free-energy cycle does not close")
        for mode in ("free_ions", "contact_pair"):
            values = row["conditional_grafting_kJ_mol"][mode]
            if abs(values["second_j_given_i"] - values["first_j"] - value) > 1e-6:
                raise ValueError("Conditional route differs from the neighboring-site cycle")
    sources = {str(path): sha(path) for path in (args.review, Path(__file__))}
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.2, 4.15))
    fig.subplots_adjust(left=.09, right=.98, top=.83, bottom=.27, wspace=.43)
    positions = np.arange(len(rows))
    labels = ["Adjacent (2,3)", "One spacer (2,4)"]
    energies = [row["electronic_deltaDeltaE_kJ_mol"] for row in rows]
    free = np.array([row["free_energy_deltaDeltaG_kJ_mol"] for row in rows])
    low = np.array([min(row["rotor_cutoff_sensitivity_kJ_mol"].values()) for row in rows])
    high = np.array([max(row["rotor_cutoff_sensitivity_kJ_mol"].values()) for row in rows])
    left.axhline(0, color=".6", lw=.7)
    left.bar(positions - .18, energies, width=.31, color=".65", label="Electronic ΔΔE")
    left.bar(positions + .18, free, width=.31, color="#009E73", label="Native ΔΔG",
             yerr=np.vstack([free - low, high - free]), capsize=3, error_kw=dict(lw=.8))
    for index, (energy, value) in enumerate(zip(energies, free)):
        left.text(index - .18, energy - .6, f"{energy:.2f}", ha="center", va="top", fontsize=7)
        left.text(index + .18, low[index] - .6, f"{value:.2f}", ha="center", va="top", fontsize=7, weight="bold")
    left.set(xticks=positions, xticklabels=labels, ylabel="Neighboring-site nonadditivity (kJ mol⁻¹)",
             ylim=(-18, 2), title="Matched occupancy cycle")
    left.legend(frameon=False, fontsize=7, loc="lower left")
    handles = []
    for mode, color, shift, label in [("free_ions", "#0072B2", -.17, "Separated ions"),
                                       ("contact_pair", "#D55E00", .17, "Contact pair")]:
        for index, row in enumerate(rows):
            values = row["conditional_grafting_kJ_mol"][mode]
            first, second = values["first_j"], values["second_j_given_i"]
            right.plot([index + shift] * 2, [first, second], color=color, lw=1)
            right.scatter(index + shift, first, marker="o", facecolor="white", edgecolor=color, s=38, lw=1, zorder=3)
            right.scatter(index + shift, second, marker="o", color=color, s=38, zorder=3)
        handles.append(Line2D([], [], marker="o", color=color, linestyle="", label=label, markersize=5))
    right.set(xticks=positions, xticklabels=labels, ylabel="DBU-assisted grafting ΔG (kJ mol⁻¹)",
              title="Grafting at site j", ylim=(-124, -98), xlim=(-.55, len(rows) - .45))
    right.legend(handles=handles, frameon=False, loc="upper left", fontsize=7)
    for axis, letter in zip((left, right), ("a", "b")):
        axis.text(-.15, 1.07, letter, transform=axis.transAxes, weight="bold", fontsize=10)
    fig.suptitle("Conditional grafting on a matched five-repeat backbone", y=.98, fontsize=11)
    fig.text(.09, .90, "Neutral phthalimide method control · SMD MeCN · 298.15 K", fontsize=8)
    fig.text(.09, .16, "a: whiskers vary the rotor cutoff (25/50/100 cm⁻¹); they are not confidence intervals.", fontsize=7)
    fig.text(.09, .115, "b: open markers graft j onto an empty backbone; filled markers graft j when i is occupied.", fontsize=7)
    fig.text(.09, .065, "Nonadditivity includes conformation and solvent response; no occupancy distribution or electron-transfer rate.", fontsize=7)
    args.out.mkdir(parents=True, exist_ok=True)
    stem = args.out / "conditional_grafting"
    files = []
    for suffix in ("pdf", "png", "svg"):
        path = stem.with_suffix("." + suffix)
        fig.savefig(path, dpi=600)
        files.append(path)
    plt.close(fig)
    exported = []
    for row in rows:
        for mode, values in row["conditional_grafting_kJ_mol"].items():
            exported.append(dict(scenario=row["name"], mode=mode,
                            electronic_deltaDeltaE_kJ_mol=row["electronic_deltaDeltaE_kJ_mol"],
                            free_energy_deltaDeltaG_kJ_mol=row["free_energy_deltaDeltaG_kJ_mol"],
                            **values, **{"cutoff_" + k + "_kJ_mol": v for k, v in row["rotor_cutoff_sensitivity_kJ_mol"].items()}))
    with stem.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(exported[0]))
        writer.writeheader()
        writer.writerows(exported)
    files.append(stem.with_suffix(".csv"))
    if any(sha(p) != h for p, h in sources.items()):
        raise ValueError("Plot inputs changed during rendering")
    atomic_write_json(stem.with_suffix(".json"), dict(status="checked_selected_conformer_native_cycles",
                      source_sha256=sources, output_sha256={str(p.relative_to(ROOT)): sha(p) for p in files},
                      original_oligomer_panel_complete=source["original_checked"] == source["original_expected"],
                      pure_core_interaction=False, quantum_accuracy_externally_validated=False))


if __name__ == "__main__":
    main()
