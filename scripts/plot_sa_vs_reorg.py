#!/usr/bin/env python
"""2x2 landscape: synthetic accessibility (SA) vs reorganization energy (lambda), split by
   BEFORE vs AFTER functionalization  x  GAS-phase vs SOLVATED reorganization.

Axes:
  - x = SA score (RDKit sascorer; lower = more synthesizable). Same in gas & solvent.
  - y = reorganization energy (eV):
        GAS      = lambda_i (inner-sphere, Nelsen 4-point; results/reorganization.csv)
        SOLVATED = lambda_i + lambda_o, lambda_o = Marcus/Born outer-sphere (Sharma 2021 Eq.13,
                   src/redox/solvated_reorg.py) with a SASA hydrodynamic-radius proxy in MeCN.

Rows = before (parent core) / after (Merrifield-functionalized monomer).
Cols = gas / solvated.

CAVEAT baked into the figure: lambda_o (hence the SOLVATED column) is a semi-quantitative
ESTIMATE (crude single-sphere Born model + structural r_hyd proxy, ~factor-2 uncertainty).
The RELATIVE story (low-SA + low-lambda corner, before->after shift) is the reliable readout;
lambda_i (gas) is the validated, trustworthy reorganization descriptor.

  PYTHONPATH=src python scripts/plot_sa_vs_reorg.py
"""
from __future__ import annotations
import csv
import importlib.util
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox.solvated_reorg import hydrodynamic_radius_A, lambda_outer_eV  # noqa: E402

RESULTS = ROOT / "results"


def _cfg(mod):
    s = importlib.util.spec_from_file_location(mod, ROOT / "config" / f"{mod}.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


# family + before/after (parent core) classification for the validation/original set
FAMILY = {
    "anthraquinone": "quinone", "anthraquinone_parent": "quinone",
    "methyl_viologen": "viologen", "viologen": "viologen",
    "methylpyridinium": "pyridinium", "pyridinium": "pyridinium", "cyanopyridinium": "pyridinium",
    "phenothiazine": "phenothiazine", "phenothiazine_parent": "phenothiazine",
    "tempo": "nitroxide", "tempo_parent": "nitroxide",
}
PARENTS = {"anthraquinone_parent", "methyl_viologen", "methylpyridinium",
           "phenothiazine_parent", "tempo_parent"}
# Okabe-Ito colorblind-safe, print-friendly palette
FAM_COLOR = {"quinone": "#0072B2", "viologen": "#CC79A7", "pyridinium": "#009E73",
             "phenothiazine": "#E69F00", "nitroxide": "#D55E00"}
FAM_ORDER = ["quinone", "viologen", "pyridinium", "phenothiazine", "nitroxide"]
NICE = {"anthraquinone": "anthraquinone", "anthraquinone_parent": "anthraquinone",
        "methyl_viologen": "methyl viologen", "viologen": "benzyl viologen",
        "methylpyridinium": "N-methylpyridinium", "pyridinium": "N-benzylpyridinium",
        "cyanopyridinium": "cyanopyridinium", "phenothiazine": "phenothiazine",
        "phenothiazine_parent": "phenothiazine", "tempo": "TEMPO", "tempo_parent": "TEMPO"}


def load():
    ele = _cfg("electrolyte")
    eop, er = ele.SOLVENT["eps_optical"], ele.SOLVENT["eps_r"]
    probe = ele.SOLVENT["solvent_probe_radius_A"]

    smiles = {r["id"]: r["smiles"] for r in csv.DictReader((ROOT / "library" / "manifest.csv").open())}
    sa = {r["id"]: float(r["SA_score"]) for r in csv.DictReader((RESULTS / "capacity_and_proxies.csv").open())
          if r.get("SA_score") not in (None, "")}
    lam_i = {}
    for r in csv.DictReader((RESULTS / "reorganization.csv").open()):
        v = r.get("lambda_i_eV")
        if v not in (None, ""):
            lam_i.setdefault(r["id"], []).append(float(v))

    rows = []
    for gid, lis in lam_i.items():
        if gid == "ferrocene" or gid not in FAMILY or gid not in sa:
            continue
        li = statistics.mean(lis)
        a = hydrodynamic_radius_A(smiles[gid], probe, "sasa")
        lo = lambda_outer_eV(a, eop, er, z=1) if a else None
        rows.append(dict(id=gid, family=FAMILY[gid], stage=("before" if gid in PARENTS else "after"),
                         SA=sa[gid], lam_gas=li,
                         lam_solv=(li + lo) if lo is not None else None,
                         r_hyd=a, lam_o=lo))
    return rows


def _journal_style(plt):
    plt.rcParams.update({
        "font.size": 18, "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.linewidth": 1.4, "axes.edgecolor": "#2b2b2b",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 6, "ytick.major.size": 6,
        "xtick.major.width": 1.3, "ytick.major.width": 1.3,
        "savefig.dpi": 300, "figure.dpi": 120,
    })


def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from adjustText import adjust_text
    _journal_style(plt)

    xs = [r["SA"] for r in rows]
    ys = [r["lam_gas"] for r in rows] + [r["lam_solv"] for r in rows if r["lam_solv"] is not None]
    xpad = (max(xs) - min(xs)) * 0.30 + 0.1
    xlim = (min(xs) - xpad, max(xs) + xpad)
    ylim = (0, max(ys) * 1.18)

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 12.6), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.155, right=0.985, top=0.9, bottom=0.085, hspace=0.1, wspace=0.07)
    col_key = ["lam_gas", "lam_solv"]
    col_hdr = ["Gas phase   (λ$_i$)", "In acetonitrile   (λ$_i$ + λ$_o$)"]
    row_hdr = ["Before functionalization\n(parent core)", "After functionalization\n(Merrifield monomer)"]
    panel = [["a", "b"], ["c", "d"]]

    for i, stage in enumerate(["before", "after"]):
        for j, key in enumerate(col_key):
            ax = axes[i][j]
            pts = [r for r in rows if r["stage"] == stage and r[key] is not None]
            for r in pts:
                ax.scatter(r["SA"], r[key], s=150, color=FAM_COLOR[r["family"]],
                           edgecolor="white", linewidth=1.6, zorder=3, alpha=0.98)
            texts = [ax.text(r["SA"], r[key], NICE.get(r["id"], r["id"]),
                             fontsize=13, color="#1c1c1c", zorder=5) for r in pts]
            adjust_text(texts, x=[r["SA"] for r in pts], y=[r[key] for r in pts], ax=ax,
                        force_static=(1.3, 1.7), force_text=(0.7, 1.0),
                        force_explode=(0.7, 1.1), expand=(2.5, 2.9),
                        min_arrow_len=8, iter_lim=400,
                        arrowprops=dict(arrowstyle="-", color="#9a9a9a", lw=0.9))
            ax.set_xlim(*xlim); ax.set_ylim(*ylim)
            ax.grid(True, axis="y", alpha=0.25, linewidth=0.7)
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(labelsize=17)
            ax.text(0.035, 0.955, panel[i][j], transform=ax.transAxes, fontsize=23,
                    fontweight="bold", va="top", ha="left")
            if i == 0:
                ax.set_title(col_hdr[j], fontsize=20, fontweight="bold", pad=14)
            if j == 0:
                ax.set_ylabel("λ  (eV)", fontsize=19)
            if i == 1:
                ax.set_xlabel("Synthetic accessibility score", fontsize=19)

    # row descriptors on the far left
    for yc, lab in zip([0.70, 0.29], row_hdr):
        fig.text(0.028, yc, lab, rotation=90, va="center", ha="center",
                 fontsize=17, fontweight="bold", color="#2b2b2b")

    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=FAM_COLOR[f],
                      markersize=15, label=f) for f in FAM_ORDER]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False,
               bbox_to_anchor=(0.57, 1.0), fontsize=17, handletextpad=0.3, columnspacing=1.3)
    out = RESULTS / "figures" / "sa_vs_reorg_2x2.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")
    return out


def main():
    rows = load()
    print(f"{'id':22s}{'stage':8s}{'family':13s}{'SA':>5s}{'r_hyd':>7s}{'lam_i':>7s}{'lam_o':>7s}{'lam_tot':>8s}")
    for r in sorted(rows, key=lambda r: (r["stage"], r["family"])):
        print(f"{r['id']:22s}{r['stage']:8s}{r['family']:13s}{r['SA']:>5.2f}"
              f"{r['r_hyd']:>7.2f}{r['lam_gas']:>7.3f}{r['lam_o']:>7.3f}{r['lam_solv']:>8.3f}")
    plot(rows)


if __name__ == "__main__":
    main()
