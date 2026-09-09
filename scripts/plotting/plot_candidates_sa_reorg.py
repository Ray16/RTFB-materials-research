#!/usr/bin/env python
"""SA vs reorganization energy for the SIX starting candidates (Candidates.xlsx), gas vs
solvated. Only the grafted (functionalized) monomers exist for these, so this is a 1x2
gas | solvated panel (no before/after). Publication style.

  PYTHONPATH=src python scripts/plotting/plot_candidates_sa_reorg.py
"""
from __future__ import annotations
import csv
import importlib.util
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.solvated_reorg import hydrodynamic_radius_A, lambda_outer_eV  # noqa: E402
RESULTS = ROOT / "results"

ORDER = ["viologen", "ethylviologen", "pmdi", "ndi_ammonium", "mophquinone", "dmophquinone"]
FAMILY = {"viologen": "viologen", "ethylviologen": "viologen",
          "pmdi": "imide", "ndi_ammonium": "imide",
          "mophquinone": "quinone", "dmophquinone": "quinone"}
NICE = {"viologen": "methyl viologen", "ethylviologen": "ethyl viologen",
        "pmdi": "PMDI", "ndi_ammonium": "ammonium-NDI",
        "mophquinone": "MeO-phenyl quinone", "dmophquinone": "(MeO)$_2$-phenyl quinone"}
# per-molecule colours: two shades per family (viologen=pink, imide=green, quinone=blue)
MOL_COLOR = {"viologen": "#B23A78", "ethylviologen": "#E68FBE",
             "pmdi": "#127A52", "ndi_ammonium": "#6BC79E",
             "mophquinone": "#00629A", "dmophquinone": "#6FA8D8"}


def _cfg(mod):
    s = importlib.util.spec_from_file_location(mod, ROOT / "config" / f"{mod}.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


def _qlab(q):
    q = int(q)
    return "0" if q == 0 else f"{q:+d}"


def load():
    from rdkit.Chem import RDConfig
    sys.path.append(f"{RDConfig.RDContribDir}/SA_Score"); import sascorer
    from rdkit import Chem
    from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

    ele = _cfg("electrolyte")
    eop, er, probe = ele.SOLVENT["eps_optical"], ele.SOLVENT["eps_r"], ele.SOLVENT["solvent_probe_radius_A"]
    smiles = {r["id"]: r["smiles"] for r in csv.DictReader((ROOT / "library" / "manifest.csv").open())}
    # per-molecule list of couples (one per charge-state transition)
    couples = {}
    for r in csv.DictReader((RESULTS / "reorganization.csv").open()):
        if r.get("lambda_i_eV") in (None, ""):
            continue
        couples.setdefault(r["id"], []).append(
            dict(label=f"{_qlab(r['q_ox'])}/{_qlab(r['q_red'])}",
                 q_ox=int(r["q_ox"]), lam_i=float(r["lambda_i_eV"])))
    rows = []
    for gid in ORDER:
        if gid not in couples:
            print(f"  [skip] {gid}: no reorg yet"); continue
        sa = sascorer.calculateScore(Chem.MolFromSmiles(smiles[gid]))
        a = hydrodynamic_radius_A(smiles[gid], probe, "sasa")
        lo = lambda_outer_eV(a, eop, er, z=1)   # per-molecule (r_hyd is charge-state independent)
        cs = sorted(couples[gid], key=lambda c: -c["q_ox"])   # high charge -> low
        for c in cs:
            c["lam_solv"] = c["lam_i"] + lo
        rows.append(dict(id=gid, family=FAMILY[gid], SA=sa, lam_o=lo, r_hyd=a, couples=cs))
    return rows


FAM_COLOR = {"viologen": "#CC79A7", "imide": "#009E73", "quinone": "#0072B2"}


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
    # one point per molecule: mean lambda over its charge-state couples (clean demo-2x2 style)
    for r in rows:
        r["lam_i"] = sum(c["lam_i"] for c in r["couples"]) / len(r["couples"])
        r["lam_solv"] = sum(c["lam_solv"] for c in r["couples"]) / len(r["couples"])
    allsa = [r["SA"] for r in rows]
    xpad = (max(allsa) - min(allsa)) * 0.22 + 0.15
    xlim = (min(allsa) - xpad, max(allsa) + xpad)

    # per-panel y-zoom: candidates cluster tightly in lambda, so each panel scales to its own
    # range (gas lambda_i ~0.45-0.55; solvated ~1.0-1.1) to give the labels room to spread.
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.6), sharex=True, sharey=False)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.135, wspace=0.16)
    for j, (key, hdr, lab) in enumerate([("lam_i", "Gas phase   (λ$_i$)", "a"),
                                         ("lam_solv", "In acetonitrile   (λ$_i$ + λ$_o$)", "b")]):
        ax = axes[j]
        yy = [r[key] for r in rows]
        ypad = (max(yy) - min(yy)) * 0.55 + 0.02
        ylim = (min(yy) - ypad, max(yy) + ypad)
        for r in rows:
            ax.scatter(r["SA"], r[key], s=230, color=FAM_COLOR[r["family"]],
                       edgecolor="white", linewidth=2.0, zorder=3)
        texts = [ax.text(r["SA"], r[key], NICE[r["id"]], fontsize=14, color="#1c1c1c", zorder=5)
                 for r in rows]
        adjust_text(texts, x=[r["SA"] for r in rows], y=[r[key] for r in rows], ax=ax,
                    force_static=(1.1, 1.5), force_text=(0.6, 0.9), force_explode=(0.6, 1.0),
                    expand=(2.2, 2.6), min_arrow_len=8, iter_lim=400,
                    arrowprops=dict(arrowstyle="-", color="#9a9a9a", lw=0.9))
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.7); ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=17)
        ax.set_title(hdr, fontsize=20, fontweight="bold", pad=12)
        ax.set_xlabel("Synthetic accessibility score", fontsize=19)
        ax.text(0.03, 1.05, lab, transform=ax.transAxes, fontsize=23, fontweight="bold")
        ax.set_ylabel("λ  (eV)", fontsize=19)   # both panels labelled (independent y-zoom)
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=FAM_COLOR[f],
                      markersize=15, label=f) for f in ["viologen", "imide", "quinone"]]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.53, 1.0), fontsize=17)
    fig.text(0.5, 0.028, "one point per molecule (mean over its charge states); note the "
             "independent y-zoom per panel — solvation adds $\\sim$0.5--0.6 eV outer-sphere.",
             ha="center", fontsize=13, color="#666", style="italic")
    out = RESULTS / "figures" / "candidates" / "candidates_sa_reorg.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print("wrote", out)


def main():
    rows = load()
    for r in rows:
        print(f"  {r['id']:14s} {r['family']:9s} SA={r['SA']:.2f} r_hyd={r['r_hyd']:.2f} "
              f"lam_o={r['lam_o']:.3f}")
        for c in r["couples"]:
            print(f"       {c['label']:8s} lam_i={c['lam_i']:.3f} lam_solv={c['lam_solv']:.3f}")
    if len(rows) >= 6:
        plot(rows)
    else:
        print(f"only {len(rows)}/6 candidates have reorg — waiting for the rest.")


if __name__ == "__main__":
    main()
