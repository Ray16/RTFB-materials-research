#!/usr/bin/env python
"""Computed discovered leads (config/discovered_candidates.py) vs the existing ranked candidates,
from the SAME production tables (results/redox_potentials.csv, stability_disproportionation.csv,
reorganization.csv, scorecard.csv). Only computed values are plotted: a lead missing a quantity
appears only in the figures whose quantities it has, and is listed as incomplete. PNG only.

  leads_vs_existing_E1_disp.png      first reduction E1 vs disproportionation free energy
  leads_vs_existing_lambda_cap.png   inner-sphere lambda (QC-clean) vs specific capacity

Within-family comparison only: quinone E1 is offset +0.31 V vs experiment (FINDINGS #24), the
same offset for every quinone, so the ORDER among quinones is meaningful, absolute values are not.

  python scripts/plotting/candidates/plot_leads_vs_existing.py
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, grid_xy  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import load_config  # noqa: E402

R = ROOT / "results"
OUT = ROOT / "results" / "figures" / "candidates"
FAM = {"quinone (n-type)": "#0072B2", "imide (n-type)": "#009E73", "pyridine-multi-e": "#CC79A7"}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _rows(name):
    p = R / f"{name}.csv"
    return list(csv.DictReader(open(p))) if p.exists() else []


def gather():
    existing = {}
    for cfg in ("starting_candidates", "merrifield_multielectron", "redox_groups"):
        for g in load_config(cfg).GROUPS:
            existing[g["id"]] = g
    ranked = {r["id"] for r in _rows("scorecard")}
    existing = {k: v for k, v in existing.items() if k in ranked and v.get("rankable", True)}
    leads = {g["id"]: g for g in load_config("discovered_candidates").GROUPS}
    E = {}
    for r in _rows("redox_potentials"):
        E.setdefault(r["id"], {})[r["event"]] = _f(r["E_vs_Fc_V"]) if r.get("status", "ok") == "ok" else None
    D = {r["id"]: _f(r["dG_disp_kJmol"]) for r in _rows("stability_disproportionation")}
    L = {}
    for r in _rows("reorganization"):
        if not (r.get("flag") or "").strip() and _f(r.get("lambda_i_eV")) is not None:
            L.setdefault(r["id"], []).append(_f(r["lambda_i_eV"]))
    cap = {r["id"]: _f(r.get("specific_capacity_mAh_g")) for r in _rows("scorecard")
           if r.get("pool") in ("anolyte", "")}
    # states whose multiplicity was chosen by the UMA scan against the hint and has no DFT
    # spin check (audit plausibility flag): anything derived from them is PROVISIONAL
    spin_unverified = {r["id"] for r in _rows("validation/audit_sanity")
                       if r.get("check") == "mult_differs_from_hint"}
    out = []
    for kind, groups in (("existing", existing), ("lead", leads)):
        for gid, g in groups.items():
            first = "neu->red1" if "neu->red1" in E.get(gid, {}) else "ox2->ox1" if "ox2->ox1" in E.get(gid, {}) else "ox->red1"
            out.append(dict(id=gid, name=g["name"], family=g.get("family"), kind=kind,
                            E1=E.get(gid, {}).get(first), disp=D.get(gid),
                            lam=(sum(L[gid]) / len(L[gid])) if gid in L else None,
                            cap=cap.get(gid), flag=g.get("flag", ""),
                            provisional=gid in spin_unverified))
    return out


def _scatter(rows, x, y, xlabel, ylabel, fname, note):
    apply_style()
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.09, 0.40, 0.88, 0.56])
    from adjustText import adjust_text
    texts, key = [], []
    num = 0
    for r in rows:
        if r[x] is None or r[y] is None:
            continue
        num += 1
        c = FAM.get(r["family"], "#888888")
        lead = r["kind"] == "lead"
        prov = r["provisional"] and y in ("disp", "lam")
        ax.scatter(r[x], r[y], s=520 if lead else 420, marker="D" if lead else "*",
                   facecolor="white" if lead else c, edgecolor="#999999" if prov else (c if lead else "black"),
                   lw=2.6 if lead else 1.2, ls="--" if prov else "-", zorder=6)
        texts.append(ax.text(r[x], r[y], f"{num}{'?' if prov else ''}", fontsize=18,
                             fontweight="bold", color="#999999" if prov else c, zorder=7))
        key.append(f"{num}{'?' if prov else ''}  {'[new] ' if lead else ''}{r['name']}"
                   + ("  (provisional)" if prov else ""))
    adjust_text(texts, ax=ax, expand=(1.6, 1.8), arrowprops=dict(arrowstyle="-", color="#555", lw=1.2))
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); grid_xy(ax)
    hs = [plt.Line2D([], [], ls="", marker="*", ms=20, mfc="#888", mec="black", label="existing candidate"),
          plt.Line2D([], [], ls="", marker="D", ms=14, mfc="white", mec="#0072B2", mew=2.4,
                     label="new lead (discovered)")]
    hs += [plt.Line2D([], [], ls="", marker="s", ms=13, color=c, label=f) for f, c in FAM.items()]
    ax.legend(handles=hs, loc="best", fontsize=18)
    half = (len(key) + 1) // 2
    for col, chunk in enumerate((key[:half], key[half:])):
        fig.text(0.02 + 0.5 * col, 0.31, "\n".join(chunk), va="top", ha="left", fontsize=18)
    fig.text(0.02, 0.355, note, fontsize=18, va="top")
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / fname; fig.savefig(p); plt.close(fig); print(" ->", p)


def main():
    rows = gather()
    inc = [r for r in rows if r["kind"] == "lead" and None in (r["E1"], r["disp"], r["lam"])]
    print("leads with a missing quantity:", [(r["id"], {k: r[k] is not None for k in ("E1", "disp", "lam")}) for r in inc])
    _scatter(rows, "E1", "disp", "computed first reduction E1 (V vs Fc/Fc$^+$; quinones +0.31 V vs exp.)",
             "computed ΔG$_{disp}$ of the radical (kJ mol$^{-1}$; over-estimated ~24)",
             "leads_vs_existing_E1_disp.png",
             "Both E1 and ΔG_disp computed; compare within a family. '?' = provisional: dianion\n"
             "multiplicity chosen by the UMA scan (triplet), no DFT spin check yet.")
    _scatter(rows, "cap", "lam", "specific capacity (mAh g$^{-1}$, computed accessible electrons)",
             "inner-sphere λ (eV, QC-clean couples, mean)", "leads_vs_existing_lambda_cap.png",
             "Only molecules with a QC-clean λ and a scorecard capacity are shown.")


if __name__ == "__main__":
    main()
