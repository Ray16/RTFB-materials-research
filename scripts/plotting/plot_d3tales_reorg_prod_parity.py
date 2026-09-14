#!/usr/bin/env python
"""Track-2 D3TaLES reorganization sweep: 3-panel parity showing WHY the gap vs D3TaLES is not a
solvent effect.

  (a) our LC-wHPBE (D3TaLES's level)  vs  D3TaLES(reported)   <- match test
  (b) our production (SMD)            vs  D3TaLES(reported)    <- the observed gap
  (c) our production                 vs  our LC-wHPBE         <- level change only: tiny

Molecules whose inner-sphere lambda is untrustworthy (unbound vertical anion; see
finalize_reorg_qc.py) are excluded and counted. Reads comparison.csv.
  PYTHONPATH=src python scripts/plotting/plot_d3tales_reorg_prod_parity.py
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, parity_panel, FAM_COLOR  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "results" / "d3tales_reorg_validation_prod" / "comparison.csv"
OUT = ROOT / "results" / "figures" / "reorg" / "reorg_prod_parity_3panel.png"
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


def main():
    rows = list(csv.DictReader(CSV.open()))
    apply_style()
    panels = [
        ("our_d3level_eV", "d3tales_eV", "our LC-$\\omega$HPBE  vs  D3TaLES",
         "D3TaLES reported $\\lambda$ (eV)", "our LC-$\\omega$HPBE $\\lambda$ (eV)", "match test"),
        ("our_smd_eV", "d3tales_eV", "our production  vs  D3TaLES",
         "D3TaLES reported $\\lambda$ (eV)", "our production $\\lambda$ (eV)", "observed gap"),
        ("our_smd_eV", "our_d3level_eV", "our production  vs  our LC-$\\omega$HPBE",
         "our LC-$\\omega$HPBE $\\lambda$ (eV)", "our production $\\lambda$ (eV)", "level change only"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(20.5, 7.8))
    fig.subplots_adjust(top=0.74, bottom=0.12, wspace=0.42, left=0.05, right=0.985)
    for ax, (xk, yk, title, xl, yl, tag) in zip(axes, panels):
        rel = [r for r in rows if str(r.get("reliable", "True")) == "True"]
        n_excl = len(rows) - len(rel)
        pts = [(_f(r[xk]), _f(r[yk]), _fam(r["family"])) for r in rel]
        pts = [(x, y, c) for x, y, c in pts if x is not None and y is not None and LO < x < HI and LO < y < HI]
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; cs = [p[2] for p in pts]
        s = parity_panel(ax, xs, ys, cs, lim=(LO, HI), band=0.2)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        # bold title above the marginal histogram; the italic tag sits just under it
        ax.set_title(f"{title}\n" + r"$\it{" + tag.replace(" ", r"\ ") + "}$",
                     fontsize=17, fontweight="bold", pad=40, color="#222222")
        box = (f"n = {s['n']}\nMAD = {s['mad']:.3f} eV\nmedian = {s['median']:.3f}\n"
               f"bias = {s['bias']:+.3f}\nr = {s['r']:.2f}\n"
               f"{100*s['frac_in_band']:.0f}% within 0.2 eV")
        ax.text(0.045, 0.955, box, transform=ax.transAxes, va="top", ha="left", fontsize=13,
                bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#cccccc", alpha=0.94), zorder=5)
    handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=12, label=k)
               for k, c in FAM_COLOR.items() if k != "viologen"]
    handles.append(plt.Line2D([0], [0], color="#B8C0C8", lw=9, alpha=0.6, label="$\\pm$0.2 eV band"))
    axes[1].legend(handles=handles, loc="lower right", fontsize=12, title="", framealpha=0.95)
    fig.suptitle("Inner-sphere $\\lambda$ (electron): matching D3TaLES's IP-tuned LC-$\\omega$HPBE reproduces their "
                 "$\\lambda$; the production gap is functional+basis, not solvent",
                 fontsize=16.5, fontweight="bold", y=1.02)
    fig.text(0.5, 0.02, "426 reliable molecules shown; 16 unbound-anion outliers excluded (see FINDINGS #10). "
             "Shaded = $\\pm$0.2 eV agreement band; top/right histograms are the marginal $\\lambda$ distributions.",
             ha="center", fontsize=12.5, color="#555555")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  -> {OUT}")


if __name__ == "__main__":
    main()
