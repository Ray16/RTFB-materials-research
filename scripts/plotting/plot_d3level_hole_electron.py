#!/usr/bin/env python
"""D3TaLES-level match, per couple: do we reproduce D3TaLES's hole & electron reorg when we run
their EXACT level (lc_wpbe + their tuned omega + def2-svp, gas)?

Two parity panels (hole | electron), reading results/d3tales_reorg_validation_d3level/calc/*.json.
The story: we reproduce their CLEAN entries closely; the scatter is where D3TaLES's own value is
an outlier (their diffuse-free def2-SVP makes the electron/anion column noisy).

  PYTHONPATH=src python scripts/plotting/plot_d3level_hole_electron.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, parity_panel, FAM_COLOR  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CALC = ROOT / "results" / "d3tales_reorg_validation_d3level" / "calc"
OUT = ROOT / "results" / "figures" / "reorg" / "d3level_hole_electron_parity.png"
LO, HI = 0.0, 1.4


def _f(x):
    try:
        v = float(x); return v if v == v else None
    except (TypeError, ValueError):
        return None


def _fam(fam):
    for k, c in FAM_COLOR.items():
        if k in (fam or ""):
            return c
    return "#888888"


def _family_map():
    import csv
    m = {}
    src = ROOT / "results" / "d3tales_reorg_validation_d3level" / "molecules.csv"
    if src.exists():
        for r in csv.DictReader(src.open()):
            m[r["id"]] = r.get("family", "")
    return m


def main():
    fam = _family_map()
    recs = []
    for p in CALC.glob("*.json"):
        try:
            d = json.loads(p.read_text())
        except (OSError, ValueError):
            continue
        if d.get("status") in ("ok", "partial"):
            d["family"] = fam.get(d.get("id"), "")
            recs.append(d)
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 8.0))
    fig.subplots_adjust(top=0.74, bottom=0.12, wspace=0.42, left=0.07, right=0.965)
    for ax, (ours_k, d3_k, title, tag) in zip(axes, [
            ("our_hole", "d3_hole", "HOLE couple (oxidation)", "cation ↔ neutral"),
            ("our_electron", "d3_electron", "ELECTRON couple (reduction)", "neutral ↔ anion")]):
        pts = [(_f(r.get(d3_k)), _f(r.get(ours_k)), _fam(r.get("family"))) for r in recs]
        pts = [(x, y, c) for x, y, c in pts if x is not None and y is not None and LO < x < HI and LO < y < HI]
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; cs = [p[2] for p in pts]
        s = parity_panel(ax, xs, ys, cs, lim=(LO, HI), band=0.2)
        n_tail = round(s["n"] * (1 - s["frac_in_band"]))
        ax.set_xlabel("D3TaLES reported $\\lambda$ (eV)")
        ax.set_ylabel("our LC-$\\omega$HPBE $\\lambda$ (eV)")
        ax.set_title(f"{title}\n" + r"$\it{" + tag.replace(" ", r"\ ") + "}$",
                     fontsize=17, fontweight="bold", pad=40, color="#222222")
        box = (f"n = {s['n']}\nmedian |$\\Delta$| = {s['median']:.3f} eV\n"
               f"MAD = {s['mad']:.3f} eV\nr = {s['r']:.2f}\n"
               f"{100*s['frac_in_band']:.0f}% within 0.2 eV\n({n_tail} in outlier tail)")
        ax.text(0.045, 0.955, box, transform=ax.transAxes, va="top", ha="left", fontsize=13,
                bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#cccccc", alpha=0.94), zorder=5)
    handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=12, label=k)
               for k, c in FAM_COLOR.items() if k != "viologen"]
    handles.append(plt.Line2D([0], [0], color="#B8C0C8", lw=9, alpha=0.6, label="$\\pm$0.2 eV band"))
    axes[1].legend(handles=handles, loc="lower right", fontsize=12, framealpha=0.95)
    fig.suptitle("D3TaLES level reproduced (LC-$\\omega$HPBE + their tuned $\\omega$ + def2-SVP, gas): "
                 "the core agrees; scatter is D3TaLES's own outlier tail",
                 fontsize=16.5, fontweight="bold", y=1.02)
    fig.text(0.5, 0.02, "Shaded = $\\pm$0.2 eV agreement band; top/right histograms are the marginal "
             "$\\lambda$ distributions. D3TaLES's diffuse-free def2-SVP makes its electron column noisiest.",
             ha="center", fontsize=12.5, color="#555555")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  -> {OUT}  ({len(recs)} molecules)")


if __name__ == "__main__":
    main()
