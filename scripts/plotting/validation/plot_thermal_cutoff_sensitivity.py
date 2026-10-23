#!/usr/bin/env python
"""Plot deterministic nuclear-cutoff scenarios without treating them as CIs."""
import argparse
import csv
import json
import os
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--result", type=Path, required=True)
    p.add_argument("--review", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    data, review = (json.loads(p.read_text()) for p in (args.result, args.review))
    if (review["status"] != "thermal_cutoff_full_state_reconstruction_checked"
            or review["result_sha256"] != api.digest(args.result)
            or review["potential_sums"] != 117 or review["reaction_sums"] != 261):
        raise ValueError("Expected independently checked 39-endpoint/87-reaction sensitivity")
    api.verify_sources(data["source_sha256"])
    api.verify_sources(review["source_sha256"])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
        "axes.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1.25, 1]})
    ordered = sorted(data["benchmark"]["endpoints"], key=lambda r: (-r["range_mV"], r["id"], r["event"]))
    x = np.arange(1, len(ordered)+1)
    ax = axes[0]
    lower = [min(r["shift_from_50_mV"].values()) for r in ordered]
    upper = [max(r["shift_from_50_mV"].values()) for r in ordered]
    ax.vlines(x, lower, upper, color=".75", linewidth=1.1, zorder=1)
    ax.scatter(x, [r["shift_from_50_mV"]["25"] for r in ordered], s=14, marker="o", color="#0072b2", label="25 cm$^{-1}$", zorder=3)
    ax.scatter(x, [r["shift_from_50_mV"]["100"] for r in ordered], s=17, marker="^", color="#d55e00", label="100 cm$^{-1}$", zorder=3)
    ax.axhline(0, color=".35", linestyle="--", linewidth=.8)
    ax.set(xlabel="Benchmark endpoint (ordered by range)", ylabel="ΔE versus 50 cm$^{-1}$ (mV)",
           xlim=(0, 40), ylim=(-27, 13), xticks=[1, 10, 20, 30, 39])
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.text(-.15, 1.04, "a", transform=ax.transAxes, weight="bold", fontsize=11)
    ax = axes[1]
    jitter = np.random.default_rng(20261023)
    groups = (("reference:", "HCl reference", "#777777", "o"),
              ("dbu:", "DBU", "#0072b2", "s"), ("dipea:", "DIPEA", "#d55e00", "^"))
    for index, (prefix, label, color, marker) in enumerate(groups):
        rows = sorted((r for r in data["grafting"]["reactions"] if r["name"].startswith(prefix)), key=lambda r: r["range_kJ_mol"])
        values = [r["range_kJ_mol"] for r in rows]
        ax.scatter(index+jitter.uniform(-.16, .16, len(rows)), values, color=color, marker=marker, s=17, alpha=.8, linewidth=.3)
        ax.plot([index-.24, index+.24], [statistics.median(values)]*2, color="black", lw=1.2)
    ax.set(xticks=[0, 1, 2], xticklabels=[g[1] for g in groups], ylabel="ΔG range (kJ mol$^{-1}$)",
           xlabel="Reaction model", xlim=(-.5, 2.5), ylim=(0, 5))
    ax.text(-.18, 1.04, "b", transform=ax.transAxes, weight="bold", fontsize=11)
    fig.subplots_adjust(left=.085, right=.99, bottom=.23, top=.92, wspace=.42)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if any(args.out.with_suffix("."+ext).exists() for ext in ("pdf", "png", "svg", "json")):
        raise FileExistsError("Preserve existing figure outputs")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(args.out.with_suffix("."+ext), dpi=600, facecolor="white")
    plt.close(fig)
    table_paths = []
    for suffix, fields, rows in (
        ("benchmark", ["order", "id", "event", "wave", "range_mV", "shift25_mV", "shift100_mV"],
         [dict(order=i+1, **{k:r[k] for k in ("id", "event", "wave", "range_mV")},
               shift25_mV=r["shift_from_50_mV"]["25"], shift100_mV=r["shift_from_50_mV"]["100"]) for i,r in enumerate(ordered)]),
        ("metrics", ["cutoff_cm", "MAE_V", "RMSE_V", "bias_V"],
         [dict(cutoff_cm=k, **data["benchmark"]["metrics_by_cutoff"][k]) for k in ("25", "50", "100")]),
        ("reactions", ["name", "range_kJ_mol", "G25_SI_kJ_mol", "G50_SI_kJ_mol", "G100_SI_kJ_mol"],
         [dict(name=r["name"], range_kJ_mol=r["range_kJ_mol"], **{f"G{k}_SI_kJ_mol":r["SI_deltaG_by_cutoff_kJ_mol"][k] for k in ("25", "50", "100")}) for r in data["grafting"]["reactions"]])):
        path = args.out.with_name(args.out.name+"_"+suffix).with_suffix(".csv")
        with path.open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        table_paths.append(path)
    output = [args.out.with_suffix("."+ext) for ext in ("pdf", "png", "svg")]+table_paths
    meta = dict(status="checked_cutoff_scenarios_drawn", source_sha256={str(p.resolve()):api.digest(p) for p in
        (Path(__file__), args.result, args.review)}, output_sha256={str(p.resolve()):api.digest(p) for p in output},
        benchmark_endpoints=39, reaction_models=87, display_jitter_seed=20261023,
        caption="Nuclear rotor-cutoff sensitivity at fixed geometries, electronic energies and spin factors. Panel a: potentials relative to the primary 50 cm^-1 treatment; endpoints are ordered by decreasing range over 25/50/100 cm^-1. Panel b: 67 HCl-reference and 10 DBU/10 DIPEA reaction-model ranges; horizontal marks are medians. Shared chemical components make reaction points correlated. These ranges are deterministic settings scenarios, not confidence or candidate prediction intervals.")
    args.out.with_suffix(".json").write_text(json.dumps(meta, indent=2)+"\n")
    print(json.dumps(dict(status=meta["status"], files=len(output))))


if __name__ == "__main__":
    main()
