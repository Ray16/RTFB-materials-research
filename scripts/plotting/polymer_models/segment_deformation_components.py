#!/usr/bin/env python
"""Plot signed individual fragment changes and their small sum."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from redox.core.common import ROOT
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def draw(rows, colors):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 10, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(7.2, 3.1), layout="constrained")
    x = np.arange(len(rows))
    for j, offset in enumerate((-0.16, 0.16)):
        values = [row["fragments"][j]["deformation_kJ_mol"] for row in rows]
        bars = ax.bar(x + offset, values, width=0.27, color=colors[j], edgecolor="black", linewidth=0.6,
            hatch=("///", "...")[j], label=f"Fragment {j + 1}")
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, value + (0.04 if value >= 0 else -0.04),
                f"{value:+.2f}", ha="center", va="bottom" if value >= 0 else "top", fontsize=8)
    totals = [row["deformation_kJ_mol_CODATA"] for row in rows]
    ax.plot(x, totals, "D", color="black", markersize=5, label="Sum", zorder=4)
    for i, value in enumerate(totals):
        side = -1 if rows[i]["fragments"][0]["deformation_kJ_mol"] <= rows[i]["fragments"][1]["deformation_kJ_mol"] else 1
        ax.text(i + side * 0.07, value + 0.075, f"{value:+.2f}", ha="right" if side < 0 else "left", va="bottom", fontsize=8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x, ["Phthalimide\nParallel", "Phthalimide\nAntiparallel", "PMDI\nAntiparallel"])
    ax.set_ylabel("Change from sampled monomer\n(kJ mol$^{-1}$)")
    ax.set_ylim(-0.95, 1.65)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", ncols=3, frameon=False, fontsize=9)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/polymer_models/segment_deformation_components")
    parser.add_argument("--grayscale-preview", type=Path)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    data = json.loads(args.source.read_text())
    if data["status"] != "segment_deformation_and_saved_relaxation_audited" or len(data["pairs"]) != 3:
        raise ValueError("Expected the checked three-pair decomposition audit")
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    for path, digest in data["source_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Audited deformation input changed")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure = draw(data["pairs"], ("#E69F00", "#56B4E9"))
    paths = []
    for extension in ("pdf", "svg", "png"):
        path = args.out.with_suffix("." + extension)
        figure.savefig(path, dpi=600)
        paths.append(path)
    plt.close(figure)
    if args.grayscale_preview:
        preview = draw(data["pairs"], ("#aaaaaa", "#eeeeee"))
        args.grayscale_preview.parent.mkdir(parents=True, exist_ok=True)
        preview.savefig(args.grayscale_preview, dpi=200)
        plt.close(preview)
    metadata = dict(source=str(args.source), source_sha256=sha(args.source), producer_sha256=sha(__file__),
        output_sha256={str(path): sha(path) for path in paths},
        caption="Signed frozen-fragment energy changes from one sampled independently optimized monomer. Opposing terms partially cancel in two pairs. Total deformation is not a positive-definite strain penalty; geometry and single-point levels differ. This operational SMD decomposition is uncorrected for BSSE and does not isolate bare motif binding.")
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
