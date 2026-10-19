#!/usr/bin/env python
"""E° accuracy figure: computed vs experimental redox potential (V vs Fc/Fc+, MeCN).

Two panels split by COUPLE charge class (the honest ranking metric is WITHIN class, FINDINGS #2):
  left : 0/-1 couples  (neutral <-> radical anion)
  right: +1/0 couples  (cation <-> neutral; includes reductions of cations, e.g. pyridinium)
Grounded anchors of other classes (+2/+1, -1/-2) are listed in the footnote.
Grey points = OROP external benchmark at the production protocol (uniform wB97M-V/def2-TZVPD,
SMD and explicitly matched live Fc thermochemistry). Historical/mismatched rows are excluded. Coloured diamonds = our in-house experimental anchors
(the sourced MeCN benchmark and config/validation.py), same protocol. Each panel states n, MAE, signed bias and Spearman for
both sets, so the offset disagreement between them is visible, not averaged away. Nothing is
fitted or rescaled.

  results/orop_benchmark.csv, results/redox_potentials.csv, config/validation.py
  python scripts/plotting/validation/plot_orop_benchmark.py -> results/figures/validation/orop_benchmark_parity.png
"""
from __future__ import annotations
import csv
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import C, apply_style, parity_panel  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "results"
FIGDIR = ROOT / "results" / "figures" / "validation"


def _validation():
    s = importlib.util.spec_from_file_location("validation", ROOT / "config" / "validation.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m.VALIDATION


def _inhouse(results=None):
    from redox.screening.scorecard import experimental_e_points
    results = Path(results) if results is not None else RESULTS
    rows = list(csv.DictReader(open(results / "redox_potentials.csv")))
    rows = [row for row in rows if row.get("status") == "ok" and row.get("thermal_qc") == "ok"]
    charge = {(row["id"], row["event"]): int(row["q_ox"]) for row in rows}
    # The same explicitly grounded MeCN anchors used by uncertainty estimation:
    # sourced tier-A benchmark events plus independently grounded validation events.
    return [dict(exp=point["exp"], calc=point["calc"],
                 q_ox=charge[point["id"], point["event"]], label=point["id"])
            for point in experimental_e_points(rows)]


def _stats(x, y):
    d = np.asarray(y) - np.asarray(x)
    rho = spearmanr(x, y).correlation if len(x) > 2 else float("nan")
    return len(x), float(np.mean(np.abs(d))), float(np.mean(d)), rho


def main(results=None, out_dir=None):
    results = Path(results) if results is not None else RESULTS
    out_dir = Path(out_dir) if out_dir is not None else FIGDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    apply_style()
    from redox.validation.orop import calibration_row_ok
    from redox.properties.potentials import fc_live, thermal_definition
    from redox.core.common import read_result
    reference = fc_live()
    definition = thermal_definition(read_result("ferrocene", "ox"))
    all_orop = list(csv.DictReader(open(results / "orop_benchmark.csv")))
    orop = [row for row in all_orop if calibration_row_ok(row, reference, definition)]
    inh = _inhouse(results)
    panels = [(0, "0 / −1 couples", (-3.4, 0.6)),
              (1, "+1 / 0 couples", (-2.1, 2.6))]
    fig, axs = plt.subplots(1, 2, figsize=(17, 8.6))
    for ax, (q, title, lim) in zip(axs, panels):
        o = [r for r in orop if int(r["charge_ox"]) == q]
        xo = [float(r["exp"]) for r in o]; yo = [float(r["calc"]) for r in o]
        if o:
            parity_panel(ax, xo, yo, ["#8C8C8C"] * len(xo), lim=lim, band=0.2,
                         point_size=46, marginals=False)
        else:
            ax.plot(lim, lim, color="#333333", lw=1.6)
            ax.fill_between(lim, [v-.2 for v in lim], [v+.2 for v in lim], color="#B8C0C8", alpha=.28)
            ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
        # in-house anchors in this charge class (q_ox >= 1 -> oxidation-type panel)
        ih = [r for r in inh if r["q_ox"] == q]
        if ih:
            ax.scatter([r["exp"] for r in ih], [r["calc"] for r in ih], marker="D", s=150,
                       color=C["dft"], edgecolor="white", lw=1.2, zorder=5)
        if o:
            n, mae, bias, rho = _stats(xo, yo)
            txt = f"OROP  n={n}  MAE {mae:.2f} V  bias {bias:+.2f} V  ρ {rho:.2f}"
        else:
            txt = "OROP: matching thermochemistry pending"
        if ih:
            n2, mae2, bias2, _ = _stats([r["exp"] for r in ih], [r["calc"] for r in ih])
            txt += f"\nours  n={n2}  MAE {mae2:.2f} V  bias {bias2:+.2f} V"
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=16,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#BBBBBB", alpha=0.95))
        ax.set_title(title)
        ax.set_xlabel("experimental E° (V vs Fc/Fc⁺)")
        ax.set_ylabel("computed E° (V vs Fc/Fc⁺)")
    h = []
    if orop:
        h.append(plt.Line2D([0], [0], marker="o", ls="", color="#8C8C8C", ms=10, label="matched OROP benchmark"))
    if inh:
        h.append(plt.Line2D([0], [0], marker="D", ls="", color=C["dft"], ms=11, label="grounded in-house anchors"))
    h.append(plt.Rectangle((0, 0), 1, 1, fc="#B8C0C8", alpha=0.5, label="±0.2 V guide"))
    fig.legend(handles=h, loc="lower center", ncol=len(h), frameon=False, bbox_to_anchor=(0.5, .03))
    other = [r for r in inh if r["q_ox"] not in (0, 1)]
    if other:
        fig.text(0.5, -.025, "other charge classes (not shown): " + ", ".join(
            f"{r['label']} q_ox={r['q_ox']:+d} err {r['calc'] - r['exp']:+.2f} V" for r in other),
            ha="center", fontsize=10, color="#444")
    fig.suptitle("E° at the active protocol: SMD MeCN and matched live Fc thermochemistry", fontsize=18)
    fig.text(.5, .005, f"OROP historical/pending rows excluded: {len(all_orop)-len(orop)}; grounded in-house anchors only",
             ha="center", fontsize=11, color="#444")
    fig.tight_layout(rect=(0, 0.13, 1, 0.95))
    out = out_dir / "orop_benchmark_parity.png"
    fig.savefig(out); print(f"  -> {out}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--out-dir", type=Path, default=FIGDIR)
    args = parser.parse_args()
    main(args.results, args.out_dir)
