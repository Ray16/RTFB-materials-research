#!/usr/bin/env python
"""Plot independently reconstructed finite-baseline cycles and well entropy."""
import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/polymer_models/conditional_ensemble_baseline")
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    data = json.loads(args.review.read_text())
    if data["checked_native_G"] != 12 or len(data["cycles"]) != 2:
        raise ValueError("The figure requires all twelve independently checked baseline G values")
    api.verify_sources(data["source_sha256"])
    states = {row["key"]: row for row in data["states"]}
    pairs = {row["graph_key"]: row for row in data["finite_pairs"]}
    minimum_E = {key: min(states[pose]["E_eV"] for pose in row["conformer_keys"]) for key, row in pairs.items()}
    electronic, free, lower, upper = [], [], [], []
    for cycle in data["cycles"]:
        i, j = cycle["sites"]
        keys = ["occ_none", f"occ_{i}", f"occ_{j}", f"occ_{i}_{j}"]
        values = [minimum_E[key] for key in keys]
        electronic.append((values[3] + values[0] - values[1] - values[2]) * 96.48533212331001)
        sensitivity = cycle["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"]
        primary = sensitivity["50"]
        free.append(primary)
        lower.append(primary - min(sensitivity.values()))
        upper.append(max(sensitivity.values()) - primary)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "pdf.fonttype": 42, "svg.fonttype": "none"})
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3), layout="constrained", gridspec_kw={"width_ratios": [1.05, 1]})
    x = np.arange(2)
    axes[0].bar(x - 0.17, electronic, width=0.3, color="#56B4E9", edgecolor="black", linewidth=0.5, label="Lowest sampled E")
    axes[0].bar(x + 0.17, free, width=0.3, color="#E69F00", hatch="///", edgecolor="black", linewidth=0.5, label="Two-pose ensemble G")
    axes[0].errorbar(x + 0.17, free, yerr=[lower, upper], fmt="none", ecolor="black", capsize=3, linewidth=0.8)
    axes[0].axhline(0, color="black", linewidth=0.7)
    axes[0].set_xticks(x, ["Adjacent\n(2,3)", "One spacer\n(2,4)"])
    axes[0].set_ylabel("Site nonadditivity (kJ mol$^{-1}$)")
    axes[0].set_ylim(-16, 3.5)
    axes[0].legend(frameon=False, fontsize=7, loc="upper right")
    order = ["occ_none", "occ_2", "occ_3", "occ_4", "occ_2_3", "occ_2_4"]
    lowering = [pairs[key]["configurational_lowering_kJ_mol"] for key in order]
    axes[1].bar(np.arange(6), lowering, color="#777777", edgecolor="black", linewidth=0.5)
    axes[1].set_xticks(np.arange(6), ["Empty", "2", "3", "4", "2,3", "2,4"])
    axes[1].set_xlabel("Occupied sites")
    axes[1].set_ylabel("Finite-set lowering (kJ mol$^{-1}$)")
    axes[1].set_ylim(0, max(lowering) * 1.3)
    for label, axis in zip(("a", "b"), axes):
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(-0.1, 1.02, label, transform=axis.transAxes, fontweight="bold", fontsize=11)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "svg", "png"):
        figure.savefig(args.out.with_suffix("." + extension), dpi=600)
    plt.close(figure)
    metadata = dict(source_sha256={str(path): api.digest(path) for path in (Path(__file__).resolve(), args.review)},
        output_sha256={str(args.out.with_suffix("." + extension)): api.digest(args.out.with_suffix("." + extension)) for extension in ("pdf", "svg", "png")},
        electronic_kJ_mol=electronic, ensemble_G_kJ_mol=free, configurational_lowering_kJ_mol=dict(zip(order, lowering)),
        caption="Neutral phthalimide five-repeat baseline, one backbone configuration and two poses per graph. (a) Electronic site nonadditivity and two-pose ensemble free-energy nonadditivity. Whiskers span the prescribed 25, 50 and 100 cm^-1 rotor cutoffs, not statistical confidence intervals. (b) min(G_i) minus the finite-set G at 298.15 K and the primary 50 cm^-1 cutoff. These results do not establish broader conformer convergence, experimental spacing statistics, pure motif binding or electron-transfer rates.")
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
