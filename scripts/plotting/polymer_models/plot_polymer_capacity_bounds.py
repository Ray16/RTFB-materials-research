#!/usr/bin/env python
"""Plot composition-defined capacity bounds without an assumed loading prior."""
import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/polymer_models/pmdi_capacity_bounds")
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    report = json.loads(args.source.read_text())
    api.verify_sources(report["source_sha256"])
    if report["accessible_capacity_established"] or report["selected_application_target_mAh_g"] is not None:
        raise ValueError("Expected conditional upper bounds with no fitted application target")
    conversion = [row for row in report["composition_grid"] if row["reactive_handles"] == 3]
    density = [row for row in report["composition_grid"] if row["grafts"] == row["reactive_handles"]]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "pdf.fonttype": 42, "svg.fonttype": "none"})
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3), layout="constrained", sharey=True)
    for index, (axis, rows) in enumerate(zip(axes, (conversion, density))):
        x = [row["grafts"] / 3 if index == 0 else row["reactive_handles"] / row["repeat_units"] for row in rows]
        for electrons, color, marker in ((1, "#0072B2", "s"), (2, "#D55E00", "o")):
            y = [row[str(electrons) + "e"]["capacity_upper_bound_mAh_g"] for row in rows]
            axis.plot(x, y, marker=marker, color=color, linewidth=1, linestyle="--", markersize=4,
                label=f"{electrons} e per graft assumed")
        axis.axhline(report["reference_capped_two_e_nominal_mAh_g"], color="#777777", linestyle=":", linewidth=1,
            label="Capped molecule, 2 e nominal")
        axis.set_xlim(-0.04, 1.04)
        axis.set_ylim(0, 170)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(-0.1, 1.02, "ab"[index], transform=axis.transAxes, fontweight="bold", fontsize=11)
    axes[0].set_xlabel("Conversion of reactive handles")
    axes[0].set_xticks([0, 1/3, 2/3, 1], ["0", "1/3", "2/3", "1"])
    axes[0].set_ylabel("Capacity upper bound (mAh g$^{-1}$)")
    axes[1].set_xlabel("Reactive handles per repeat unit")
    axes[0].legend(frameon=True, facecolor="white", edgecolor="none", framealpha=1, fontsize=7, loc="upper left")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "svg", "png"):
        figure.savefig(args.out.with_suffix("." + extension), dpi=600)
    plt.close(figure)
    metadata = dict(source_sha256={str(path): api.digest(path) for path in (Path(__file__).resolve(), args.source)},
        output_sha256={str(args.out.with_suffix("." + extension)): api.digest(args.out.with_suffix("." + extension)) for extension in ("pdf", "svg", "png")},
        caption="Conditional PMDI stoichiometric capacity bounds. (a) Three reactive handles on five repeat units, with varied conversion. (b) Alternative handle-density scenarios at complete conversion. Points are discrete composition scenarios; connecting lines are guides. One- and two-electron budgets assume full utilization. The dry neutral H-capped chain mass retains residual chloride and excludes unlisted crosslinker/electrolyte mass. Accessible polymer capacity and actual graft loading are not established; no application cutoff or probability distribution is fitted.")
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
