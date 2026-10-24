#!/usr/bin/env python
"""Draw the independently replayed CDS step-isolation diagnostic."""
import argparse
import hashlib
import json
import os
from pathlib import Path

from redox.core.io import atomic_write_json


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--review", type=Path, required=True)
    p.add_argument("--out-prefix", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    sources = [args.input.resolve(), args.review.resolve(), Path(__file__).resolve()]
    hashes = {str(f): sha(f) for f in sources}
    data, review = (json.loads(f.read_text()) for f in sources[:2])
    if review["source_sha256"].get(str(sources[0])) != hashes[str(sources[0])]:
        raise ValueError("Independent review belongs to different diagnostic bytes")
    outputs = [args.out_prefix.with_suffix(s) for s in (".pdf", ".png", ".svg", ".json")]
    if any(p.exists() for p in outputs):
        raise FileExistsError("Preserve existing figure outputs")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), layout="constrained")
    colors = ("#0072B2", "#D55E00")
    bars, labels = [], []
    for record, color in zip(data["records"], colors):
        pairs = sorted((float(k), v*1e6) for k, v in record["CDS_derivative_by_step_Ha_Bohr"].items())
        axes[0].semilogx([h for h, _ in pairs], [v for _, v in pairs], "o-", ms=4,
            color=color, label=f"Atom {record['atom']}, y")
        for r in record["total_energy_comparisons"]:
            bars.append((abs(r["total_numeric_minus_analytic_Ha_Bohr"]), abs(r["remaining_difference_Ha_Bohr"])))
            labels.append(f"{record['atom']}/y\n{r['step_Bohr']*1000:g}")
    axes[0].set(xlabel="Displacement (bohr)", ylabel=r"CDS derivative ($\mu E_\mathrm{h}/a_0$)", title="CDS term versus displacement")
    axes[0].legend(frameon=False, fontsize=8)
    positions = np.arange(len(bars))
    for offset, column, label, color in ((-.18, 0, "Total-force discrepancy", colors[0]),
            (.18, 1, "After CDS step isolation", colors[1])):
        axes[1].bar(positions+offset, [v[column] for v in bars], width=.36, label=label, color=color)
    axes[1].axhline(1e-6, color=".35", linestyle="--", lw=.8, label="Original tolerance")
    axes[1].set(yscale="log", ylim=(1e-9, 1e-3), xticks=positions, xticklabels=labels,
        xlabel="Atom/axis; displacement (millibohr)", ylabel=r"Absolute difference ($E_\mathrm{h}/a_0$)",
        title="Isolating the sampled discrepancy")
    axes[1].legend(frameon=False, fontsize=7, loc="upper left")
    for ax, letter in zip(axes, "ab"):
        ax.text(-.15, 1.03, letter, weight="bold", transform=ax.transAxes)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for path in outputs[:3]:
        fig.savefig(path, dpi=300)
    plt.close(fig)
    if any(sha(p) != h for p, h in hashes.items()):
        raise ValueError("Figure inputs changed")
    caption = ("Fixed-geometry site-2 PMDI CDS diagnostic. Atom indices are zero-based. "
        "Panel a uses seven finite steps; panel b isolates the CDS contribution to four stored total-force checks. "
        "The remaining differences are not newly computed total forces. Original failures and the 1e-6 Ha/Bohr gate remain unchanged; "
        "the separate fine-step quantum diagnostic, native Hessian and G are not qualified by this plot.")
    atomic_write_json(outputs[-1], dict(status="checked_component_diagnostic_drawn",
        source_sha256=hashes, output_sha256={str(f.resolve()): sha(f) for f in outputs[:3]}, caption=caption))
    print(json.dumps(dict(status="component_diagnostic_drawn", outputs=list(map(str, outputs)))))


if __name__ == "__main__":
    main()
