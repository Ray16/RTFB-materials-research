#!/usr/bin/env python
"""Manuscript Figure 4 — validation parity plot (JCIM).

Computed vs experimental redox potential (V vs Fc/Fc+, MeCN), one point per scored
endpoint (redox couple / wave), COLOURED BY FAMILY. Single panel per house rule
(scripts/plotting/plot_style.py): 18 pt text, 300 dpi, CVD-safe Okabe-Ito palette.

Includes the y=x identity line, a shaded +/-MAE agreement band, and an annotation box
with MAE / RMSE / bias / n (and the molecule + study counts). All numbers are READ from
the frozen benchmark.json; none are hardcoded or rescaled. A companion residual histogram
is written as a SEPARATE file (fig4b_residuals) so each figure stays one panel.
PNG only, 300 dpi, per manuscript/figures/STYLE_CONTRACT.md (no .pdf).

  python scripts/plotting/manuscript/fig4_validation_parity.py

Data (frozen benchmark record, READ ONLY):
  manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json
Outputs:
  manuscript/figures/fig4_validation_parity.png
  manuscript/figures/fig4b_residuals.png
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, C, OKABE_ITO, grid_x  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json"
FIGDIR = ROOT / "manuscript" / "figures"

# Family -> colour per STYLE_CONTRACT.md section 2 (identical hues across all figures).
# nitroxide uses yellow #F0E442 to stay distinct from any orange; phenothiazine is purple.
FAM_COLOR = {
    "quinone":        OKABE_ITO["blue"],          # #0072B2
    "aromatic imide": OKABE_ITO["green"],         # #009E73
    "viologen":       OKABE_ITO["reddishpurple"], # #CC79A7
    "phenothiazine":  "#762A83",                  # purple
    "nitroxide":      OKABE_ITO["yellow"],        # #F0E442
}
FAM_LABEL = {
    "quinone": "quinone (n-type)",
    "aromatic imide": "aromatic imide (n-type)",
    "viologen": "viologen (n-type)",
    "phenothiazine": "phenothiazine (p-type)",
    "nitroxide": "nitroxide (p-type)",
}


def load():
    d = json.loads(DATA.read_text())
    rows = d["outliers"]  # despite the name, this is the full scored-endpoint list
    m = d["metrics"]
    return d, rows, m


def main():
    d, rows, m = load()
    n = m["n"]
    mae, rmse, bias = m["MAE_V"], m["RMSE_V"], m["bias_V"]
    maxerr = m["maximum_absolute_error_V"]
    n_mol = d["admitted_molecules"]
    n_std = d["original_studies"]

    # integrity check: the data must agree with the numbers quoted in the manuscript.
    assert n == 39, f"expected n=39, got {n}"
    assert len(rows) == 39, f"expected 39 endpoints, got {len(rows)}"
    assert n_mol == 30 and n_std == 7, (n_mol, n_std)
    assert abs(mae - 0.247) < 0.001 and abs(rmse - 0.279) < 0.001, (mae, rmse)
    assert round(bias, 2) == 0.19, bias

    x = np.array([r["exp_V"] for r in rows], float)
    y = np.array([r["calc_V"] for r in rows], float)
    fams = [r["family"] for r in rows]

    apply_style()

    # ---- Figure 4: parity --------------------------------------------------
    fig, ax = plt.subplots(figsize=(10.5, 10.5))
    lo, hi = -2.3, 1.45  # keeps every point incl. the worst phenothiazine outlier visible

    # +/-MAE agreement band + identity line
    band = mae
    ax.fill_between([lo, hi], [lo - band, hi - band], [lo + band, hi + band],
                    color="#B8C0C8", alpha=0.30, lw=0, zorder=1,
                    label=rf"$\pm$0.25 V (MAE {mae:.3f})")
    ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=2, label="1:1")

    # one scatter call per family (for a clean legend)
    for fam in ["quinone", "aromatic imide", "phenothiazine", "viologen", "nitroxide"]:
        sel = [i for i, f in enumerate(fams) if f == fam]
        if not sel:
            continue
        ax.scatter(x[sel], y[sel], s=95, zorder=4,
                   facecolors=FAM_COLOR[fam], edgecolors="#222222", linewidths=0.6,
                   alpha=0.9, label=FAM_LABEL[fam])

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(r"experimental $E$ (V vs Fc/Fc$^+$, MeCN)")
    ax.set_ylabel(r"computed $E$ (V vs Fc/Fc$^+$)")
    ax.set_axisbelow(True)
    ax.grid(True, color=C["grid"], lw=0.6, alpha=0.3)

    # stats annotation (read, not hardcoded). 18 pt floor per style contract.
    ax.text(0.025, 0.975,
            f"MAE = {mae:.3f} V\nRMSE = {rmse:.3f} V\n"
            f"bias = {bias:+.2f} V\n"
            f"n = {n} ({n_mol} molecules, {n_std} studies)",
            transform=ax.transAxes, va="top", ha="left", fontsize=18,
            bbox=dict(boxstyle="round,pad=0.5", fc="#FBFBF7", ec="#B9B9B9", lw=1.2))

    ax.legend(loc="lower right", framealpha=0.95, fontsize=18)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / "fig4_validation_parity.png", bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 4b: signed residual histogram (separate, one panel) --------
    err = np.array([r["error_V"] for r in rows], float)  # calc - exp
    fig2, ax2 = plt.subplots(figsize=(10.0, 7.0))
    bins = np.arange(-0.6, 0.65, 0.1)
    counts, _, _ = ax2.hist(err, bins=bins, color=OKABE_ITO["skyblue"],
                            edgecolor="#222222", lw=0.8, zorder=3)
    ax2.set_ylim(0, counts.max() * 1.30)  # headroom so the 18 pt box clears the tall bars
    ax2.axvline(0.0, color="#333333", lw=1.4, zorder=4, label="zero error")
    ax2.axvline(bias, color=C["warn"], lw=2.0, ls="--", zorder=5,
                label=f"mean bias {bias:+.2f} V")
    ax2.set_xlabel(r"signed residual  $E_{\mathrm{calc}} - E_{\mathrm{exp}}$  (V)")
    ax2.set_ylabel("endpoints")
    ax2.set_xlim(-0.6, 0.6)
    grid_x(ax2)
    ax2.set_axisbelow(True)
    ax2.grid(True, axis="y", color=C["grid"], lw=0.7, alpha=0.4)
    ax2.legend(loc="upper left", framealpha=0.95, fontsize=18)
    ax2.text(0.975, 0.975,
             f"n = {n}\nmax|err| = {maxerr:.2f} V",
             transform=ax2.transAxes, va="top", ha="right", fontsize=18,
             bbox=dict(boxstyle="round,pad=0.5", fc="#FBFBF7", ec="#B9B9B9", lw=1.2))
    fig2.savefig(FIGDIR / "fig4b_residuals.png", bbox_inches="tight")
    plt.close(fig2)

    print(f"n={n}  MAE={mae:.3f}  RMSE={rmse:.3f}  bias={bias:+.3f}  "
          f"max|err|={maxerr:.3f}  mols={n_mol}  studies={n_std}")
    print(f"  -> {FIGDIR/'fig4_validation_parity.png'}  (PNG only)")
    print(f"  -> {FIGDIR/'fig4b_residuals.png'}  (PNG only)")


if __name__ == "__main__":
    main()
