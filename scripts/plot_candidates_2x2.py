#!/usr/bin/env python
"""2x2 for the six starting candidates: STRUCTURE (standalone vs grafted) x ENVIRONMENT
(gas vs acetonitrile). SA vs reorganization energy, one point per molecule (mean over its
charge-state couples). Rows = standalone (un-tethered core) / grafted (Merrifield monomer);
cols = gas (lambda_i) / MeCN (lambda_i + lambda_o).

  PYTHONPATH=src python scripts/plot_candidates_2x2.py
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

# (nice name, family, standalone id, grafted id)
MOL = [
    ("methyl viologen", "viologen", "methyl_viologen", "viologen"),
    ("ethyl viologen", "viologen", "ethylviologen_sa", "ethylviologen"),
    ("PMDI", "imide", "pmdi_sa", "pmdi"),
    ("ammonium-NDI", "imide", "ndi_ammonium_sa", "ndi_ammonium"),
    ("MeO-phenyl quinone", "quinone", "mophquinone_sa", "mophquinone"),
    ("(MeO)$_2$-phenyl quinone", "quinone", "dmophquinone_sa", "dmophquinone"),
]
FAM_COLOR = {"viologen": "#CC79A7", "imide": "#009E73", "quinone": "#0072B2"}


def _cfg(mod):
    s = importlib.util.spec_from_file_location(mod, ROOT / "config" / f"{mod}.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


def load():
    from rdkit.Chem import RDConfig
    sys.path.append(f"{RDConfig.RDContribDir}/SA_Score"); import sascorer
    from rdkit import Chem
    from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")
    ele = _cfg("electrolyte")
    eop, er, probe = ele.SOLVENT["eps_optical"], ele.SOLVENT["eps_r"], ele.SOLVENT["solvent_probe_radius_A"]
    smiles = {r["id"]: r["smiles"] for r in csv.DictReader((ROOT / "library" / "manifest.csv").open())}
    lam = {}
    for r in csv.DictReader((RESULTS / "reorganization.csv").open()):
        if r.get("lambda_i_eV") not in (None, ""):
            lam.setdefault(r["id"], []).append(float(r["lambda_i_eV"]))

    def one(gid):
        if gid not in lam or gid not in smiles:
            return None
        li = statistics.mean(lam[gid])
        a = hydrodynamic_radius_A(smiles[gid], probe, "sasa")
        lo = lambda_outer_eV(a, eop, er, z=1)
        sa = sascorer.calculateScore(Chem.MolFromSmiles(smiles[gid]))
        return dict(SA=sa, lam_i=li, lam_solv=li + lo, lam_o=lo)

    rows = []
    for nice, fam, sa_id, gr_id in MOL:
        st, gr = one(sa_id), one(gr_id)
        rows.append(dict(nice=nice, family=fam, standalone=st, grafted=gr,
                         sa_id=sa_id, gr_id=gr_id))
    return rows


def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from adjustText import adjust_text
    plt.rcParams.update({"font.size": 18, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "axes.linewidth": 1.4, "axes.edgecolor": "#2b2b2b",
                         "xtick.major.width": 1.3, "ytick.major.width": 1.3})
    # sharey per column (gas vs MeCN differ in range; standalone/grafted rows share a column)
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 12.6), sharex=True, sharey="col")
    fig.subplots_adjust(left=0.10, right=0.985, top=0.90, bottom=0.075, hspace=0.10, wspace=0.10)
    allsa = [d[k]["SA"] for r in rows for d, k in [(r, "standalone"), (r, "grafted")] if d[k]]
    xpad = (max(allsa) - min(allsa)) * 0.22 + 0.15
    xlim = (min(allsa) - xpad, max(allsa) + xpad)
    col_hdr = ["Gas phase   (λ$_i$)", "In acetonitrile   (λ$_i$ + λ$_o$)"]
    row_hdr = ["Standalone\n(un-tethered core)", "Grafted\n(Merrifield monomer)"]
    panel = [["a", "b"], ["c", "d"]]

    for i, rowkey in enumerate(["standalone", "grafted"]):
        for j, ykey in enumerate(["lam_i", "lam_solv"]):
            ax = axes[i][j]
            pts = [(r, r[rowkey]) for r in rows if r[rowkey]]
            for r, d in pts:
                ax.scatter(r["SA"] if False else d["SA"], d[ykey], s=220,
                           color=FAM_COLOR[r["family"]], edgecolor="white", linewidth=2.0, zorder=3)
            texts = [ax.text(d["SA"], d[ykey], r["nice"], fontsize=13, color="#1c1c1c", zorder=5)
                     for r, d in pts]
            adjust_text(texts, x=[d["SA"] for _, d in pts], y=[d[ykey] for _, d in pts], ax=ax,
                        force_static=(1.1, 1.5), force_text=(0.6, 0.9), force_explode=(0.6, 1.0),
                        expand=(2.2, 2.6), min_arrow_len=8, iter_lim=400,
                        arrowprops=dict(arrowstyle="-", color="#9a9a9a", lw=0.9))
            ax.set_xlim(*xlim)
            ax.grid(True, axis="y", alpha=0.25, linewidth=0.7); ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=16)
            ax.text(0.035, 0.955, panel[i][j], transform=ax.transAxes, fontsize=22,
                    fontweight="bold", va="top")
            if i == 0:
                ax.set_title(col_hdr[j], fontsize=20, fontweight="bold", pad=12)
            if i == 1:
                ax.set_xlabel("Synthetic accessibility score", fontsize=18)
            if j == 0:
                ax.set_ylabel("λ  (eV)", fontsize=18)
    for yc, lab in zip([0.70, 0.29], row_hdr):
        fig.text(0.022, yc, lab, rotation=90, va="center", ha="center",
                 fontsize=16, fontweight="bold", color="#2b2b2b")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=FAM_COLOR[f],
                      markersize=15, label=f) for f in ["viologen", "imide", "quinone"]]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.55, 1.0), fontsize=17)
    out = RESULTS / "figures" / "candidates_2x2_standalone_grafted.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print("wrote", out)


def main():
    rows = load()
    ready = sum(1 for r in rows if r["standalone"] and r["grafted"])
    for r in rows:
        s = f"{r['standalone']['lam_i']:.3f}" if r["standalone"] else "—"
        g = f"{r['grafted']['lam_i']:.3f}" if r["grafted"] else "—"
        print(f"  {r['nice']:24s} standalone λi={s}  grafted λi={g}")
    print(f"{ready}/6 candidates have both standalone+grafted")
    if ready >= 5:
        plot(rows)
    else:
        print("waiting for more standalone reorg")


if __name__ == "__main__":
    main()
