#!/usr/bin/env python
"""E° accuracy figure: computed vs experimental redox potential (V vs Fc/Fc+, MeCN).

Two panels split by COUPLE charge class (the honest ranking metric is WITHIN class, FINDINGS #2):
  left : 0/-1 couples  (neutral <-> radical anion)
  right: +1/0 couples  (cation <-> neutral; includes reductions of cations, e.g. pyridinium)
In-house anchors of other classes (+2/+1, -1/-2) have no OROP counterpart and are counted in
the footnote, not plotted.
Grey points = OROP external benchmark at the production protocol (uniform wB97M-V/def2-TZVPD,
SMD, xTB RRHO, live level-matched Fc). Coloured diamonds = our in-house experimental anchors
(config/validation.py), same protocol. Each panel states n, MAE, signed bias and Spearman for
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
FIGDIR = ROOT / "results" / "figures" / "validation"
FIGDIR.mkdir(parents=True, exist_ok=True)


def _validation():
    s = importlib.util.spec_from_file_location("validation", ROOT / "config" / "validation.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m.VALIDATION


def _inhouse():
    comp = {}
    for r in csv.DictReader(open(ROOT / "results" / "redox_potentials.csv")):
        if (r.get("E_vs_Fc_V") or "").strip():
            comp[(r["id"], r["event"])] = (float(r["E_vs_Fc_V"]), int(r["q_ox"]))
    out = []
    for e in _validation():
        if e["id"] == "ferrocene":          # the reference itself (0 by construction)
            continue
        for ev in e.get("events", []):
            c = comp.get((e["id"], ev["event"]))
            if c:
                out.append(dict(exp=ev["exp_V_vs_Fc"], calc=c[0], q_ox=c[1], label=e["id"]))
    return out


def _stats(x, y):
    d = np.asarray(y) - np.asarray(x)
    rho = spearmanr(x, y).correlation if len(x) > 2 else float("nan")
    return len(x), float(np.mean(np.abs(d))), float(np.mean(d)), rho


def main():
    apply_style()
    orop = [r for r in csv.DictReader(open(ROOT / "results" / "orop_benchmark.csv"))]
    inh = _inhouse()
    panels = [(0, "0 / −1 couples", (-3.4, 0.6)),
              (1, "+1 / 0 couples", (-2.1, 2.6))]
    fig, axs = plt.subplots(1, 2, figsize=(17, 8.6))
    for ax, (q, title, lim) in zip(axs, panels):
        o = [r for r in orop if int(r["charge_ox"]) == q]
        xo = [float(r["exp"]) for r in o]; yo = [float(r["calc"]) for r in o]
        parity_panel(ax, xo, yo, ["#8C8C8C"] * len(xo), lim=lim, band=0.2, point_size=46,
                     marginals=False)
        # in-house anchors in this charge class (q_ox >= 1 -> oxidation-type panel)
        ih = [r for r in inh if r["q_ox"] == q]
        if ih:
            ax.scatter([r["exp"] for r in ih], [r["calc"] for r in ih], marker="D", s=150,
                       color=C["dft"], edgecolor="white", lw=1.2, zorder=5)
        n, mae, bias, rho = _stats(xo, yo)
        txt = f"OROP  n={n}  MAE {mae:.2f} V  bias {bias:+.2f} V  ρ {rho:.2f}"
        if ih:
            n2, mae2, bias2, _ = _stats([r["exp"] for r in ih], [r["calc"] for r in ih])
            txt += f"\nours  n={n2}  MAE {mae2:.2f} V  bias {bias2:+.2f} V"
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=16,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#BBBBBB", alpha=0.95))
        ax.set_title(title)
        ax.set_xlabel("experimental E° (V vs Fc/Fc⁺)")
        ax.set_ylabel("computed E° (V vs Fc/Fc⁺)")
    h = [plt.Line2D([0], [0], marker="o", ls="", color="#8C8C8C", ms=10, label="OROP benchmark"),
         plt.Line2D([0], [0], marker="D", ls="", color=C["dft"], ms=11, label="in-house anchors"),
         plt.Rectangle((0, 0), 1, 1, fc="#B8C0C8", alpha=0.5, label="±0.2 V")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.04))
    other = [r for r in inh if r["q_ox"] not in (0, 1)]
    if other:
        fig.text(0.5, -0.085, "not plotted (no OROP class): " + ", ".join(
            f"{r['label']} q_ox={r['q_ox']:+d} err {r['calc'] - r['exp']:+.2f} V" for r in other),
            ha="center", fontsize=15, color="#444")
    fig.suptitle("E° at the production protocol (ωB97M-V/def2-TZVPD, SMD MeCN, xTB RRHO, "
                 "level-matched Fc) — nothing fitted", fontsize=18)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = FIGDIR / "orop_benchmark_parity.png"
    fig.savefig(out); print(f"  -> {out}")


if __name__ == "__main__":
    main()
