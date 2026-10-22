#!/usr/bin/env python
"""Plot a checked retrospective benchmark without dropping pending or outlier rows."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.validation.benchmark_panel import review_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--panel", type=Path, default=ROOT / "results/validation/publication_benchmark/panel_v1.json")
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/validation")
    args = parser.parse_args()
    panel = json.loads(args.panel.read_text())
    report = review_report(panel, args.review)
    targets = {(t["id"], t["event"], t["doi"]): t for t in panel["targets"]}
    rows = []
    with (args.review / "comparisons.csv").open() as stream:
        for row in csv.DictReader(stream):
            key = row["id"], row["event"], row["doi"]
            if key not in targets or row["kind"] != "E" or row["status"] != "checked_native_parent":
                continue
            target = targets[key]
            rows.append(dict(id=row["id"], event=row["event"], doi=row["doi"], family=target["family"], wave=target["wave"],
                             exp_V=float(row["exp_V"]), calc_V=float(row["calc_unrounded_V"]),
                             error_V=float(row["calc_unrounded_V"]) - float(row["exp_V"])))
    errors = np.array([row["error_V"] for row in rows])
    if len(rows) != report["admitted_events"] or not np.isfinite(errors).all():
        raise ValueError("Figure rows differ from independently admitted coverage")
    if abs(np.mean(abs(errors)) - report["metrics"]["MAE_V"]) > 1e-12:
        raise ValueError("Figure arithmetic differs from unrounded benchmark")
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    sources = {str(p): sha(p) for p in [args.panel, Path(__file__), *args.review.glob("*.json"), args.review / "comparisons.csv"]}
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.2, 3.8), gridspec_kw={"width_ratios": [1, 1.15]})
    fig.subplots_adjust(left=.09, right=.98, bottom=.23, top=.82, wspace=.48)
    families = sorted(report["by_family"])
    colors = dict(zip(families, ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00"]))
    markers = {1: "o", 2: "^"}
    limits = (min(min(r["exp_V"], r["calc_V"]) for r in rows) - .15,
              max(max(r["exp_V"], r["calc_V"]) for r in rows) + .15)
    left.plot(limits, limits, color=".5", lw=.8, zorder=1)
    right.axvline(0, color=".5", lw=.8)
    family_counts = {name: 0 for name in families}
    for row in rows:
        name = row["family"]
        i = family_counts[name]
        family_counts[name] += 1
        jitter = ((i % 5) - 2) * .055
        style = dict(color=colors[name], marker=markers[row["wave"]], s=25,
                     edgecolor="black", linewidth=.35, alpha=.9, zorder=3)
        left.scatter(row["exp_V"], row["calc_V"], **style)
        right.scatter(row["error_V"], families.index(name) + jitter, **style)
    left.set(xlim=limits, ylim=limits, xlabel="Measured potential (V vs Fc/Fc⁺)",
             ylabel="Calculated potential (V vs Fc/Fc⁺)")
    left.set_aspect("equal", adjustable="box")
    right.set(yticks=range(len(families)), yticklabels=[f"{f}\n(n={family_counts[f]})" for f in families],
              xlabel="Calculated − measured (V)", ylim=(-.5, len(families) - .5))
    right.invert_yaxis()
    right.tick_params(axis="y", length=0, labelsize=7)
    for ax, letter in zip((left, right), ("a", "b")):
        ax.grid(axis="x", alpha=.12)
        ax.text(-.15, 1.05, letter, transform=ax.transAxes, weight="bold", fontsize=10)
    left.text(.03, .97, f'MAE {report["metrics"]["MAE_V"]:.3f} V\nRMSE {report["metrics"]["RMSE_V"]:.3f} V\nBias {report["metrics"]["bias_V"]:+.3f} V',
              transform=left.transAxes, va="top", fontsize=7,
              bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
    fig.suptitle(f'Exact free molecules: {len(rows)}/{report["target_events"]} potentials · {report["admitted_molecules"]} molecules · {report["original_studies"]} studies', y=.98, fontsize=10)
    handles = [Line2D([], [], marker=markers[w], color="black", linestyle="", label=f"Redox step {w}", markersize=5) for w in (1, 2)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .92), ncol=2, frameon=False, fontsize=7)
    fig.text(.09, .10, f'Retrospective MeCN panel; {report["pending_events"]} pending. All admitted outliers retained; no fitted correction.', fontsize=7)
    fig.text(.09, .045, 'Matched free-parent accuracy only; grafted structures, polymers and kinetics remain unvalidated.', fontsize=7)
    args.out.mkdir(parents=True, exist_ok=True)
    stem = args.out / "retrospective_redox_benchmark"
    outputs = []
    for suffix in ("pdf", "png", "svg"):
        path = stem.with_suffix("." + suffix)
        fig.savefig(path, dpi=600)
        outputs.append(path)
    plt.close(fig)
    path = stem.with_suffix(".csv")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    outputs.append(path)
    tex = ROOT / "slides/generated/retrospective_redox_benchmark.tex"
    tex.write_text("\n".join([r"\begin{frame}{Retrospective validation: expanded exact-parent panel}", r"\centering",
                   r"\includegraphics[width=.98\linewidth,height=.76\textheight,keepaspectratio]{../results/figures/validation/retrospective_redox_benchmark.pdf}",
                   r"\par\smallskip\footnotesize " + f'{len(rows)}/{report["target_events"]} admitted potentials. MAE {report["metrics"]["MAE_V"]:.3f} V; RMSE {report["metrics"]["RMSE_V"]:.3f} V. ',
                   r"No fitted correction; all outliers retained. Exact-graft accuracy remains unvalidated.", r"\end{frame}"]) + "\n")
    outputs.append(tex)
    if any(sha(p) != h for p, h in sources.items()):
        raise ValueError("Benchmark snapshot changed while drawing")
    atomic_write_json(stem.with_suffix(".json"), dict(status="checked_retrospective_snapshot", report=report,
                      source_sha256=sources, output_sha256={str(p.relative_to(ROOT)): sha(p) for p in outputs},
                      style_scope="Color identifies family; shape identifies wave; labels retain grayscale interpretation"))


if __name__ == "__main__":
    main()
