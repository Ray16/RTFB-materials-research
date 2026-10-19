#!/usr/bin/env python
"""Computed discovered leads (config/discovered_candidates.py) vs the existing ranked candidates,
from the SAME production tables (results/redox_potentials.csv, stability_disproportionation.csv,
reorganization.csv, scorecard.csv). Only computed values are plotted: a lead missing a quantity
appears only in the figures whose quantities it has, and is listed as incomplete. PNG only.

  leads_vs_existing_E1_disp.png      first reduction E1 vs disproportionation free energy
  leads_vs_existing_lambda_cap.png   inner-sphere lambda (QC-clean) vs specific capacity

Historical points retain their original reference and thermal treatment. Native points
use the independently audited matched Fc reference; experimental error of this revised
model is uncalibrated. Historical empirical offsets do not qualify the native model.

  python scripts/plotting/candidates/plot_leads_vs_existing.py
"""
import csv
import json
import math
import textwrap
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, grid_xy  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import load_config  # noqa: E402
from verified_overlay import load as load_overlay,short,family,snapshot_inputs_current

R = ROOT / "results"
OUT = ROOT / "results" / "figures" / "candidates"
FAM = {"quinone (n-type)": "#0072B2", "imide (n-type)": "#009E73", "pyridine-multi-e": "#CC79A7",
       'phenothiazine (p-type)':'#762A83'}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _rows(name):
    p = R / f"{name}.csv"
    return list(csv.DictReader(open(p))) if p.exists() else []


def gather(feed=None):
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
                            provisional=gid in spin_unverified,
                            withdrawn_spin=(gid=='d3_80izcq'),gas_events=len(L.get(gid,[])),
                            declared_events=len(g['states'])-1))
    # Preserve earlier plotted coordinates. Their thermal/reference definition is
    # historical; new points come only from the audited native feed.
    payload,allrows,admitted=load_overlay(feed)
    if payload:
        couples=list(csv.DictReader(open(Path(payload['snapshot_path'])/'couples.csv')))
        for r in admitted:
            gas=[]
            for c in couples:
                if c['proposal_id']!=r['id']:continue
                value=_f(c.get('gas_Nelsen_sum_eV'))
                if (value is not None and math.isfinite(value)
                        and c.get('gas_Nelsen_status')=='QC-clean gas diagnostic'
                        and not c.get('gas_Nelsen_flags')):gas.append(value)
            out.append(dict(id=r['id'],name=short(r['id']),family=family(r['family']),kind='expansion',
                            E1=r.get('E1_vs_Fc_V'),disp=r.get('dG_disproportionation_kJ_mol'),
                            lam=sum(gas)/len(gas) if gas else None,cap=r.get('accessible_capacity_mAh_g'),
                            provisional=False,flag='',withdrawn_spin=False,
                            gas_events=len(gas),declared_events=r['declared_states']-1,
                            snapshot=payload['as_of_Chicago']))
    # The assessment independently reconstructed a small prior-cohort subset
    # under the SAME native definition. Add these without removing old points.
    prior=ROOT/'results/validation/resume_campaign_2026-10-12/candidate_assessment/existing_candidate_native_reconstruction.json'
    if prior.exists() and payload and snapshot_inputs_current(payload['summary']):
        data=json.loads(prior.read_text())
        old_reference=data.get('reference_native_Fc_absolute_V')
        current_reference=payload['summary'].get('reference_native_Fc_absolute_V')
        matched=(isinstance(old_reference,(int,float)) and isinstance(current_reference,(int,float))
                 and math.isfinite(old_reference) and math.isfinite(current_reference)
                 and abs(old_reference-current_reference)<1e-8)
        for r in data['candidates'] if matched else []:
            group=existing.get(r['id']) or leads.get(r['id'])
            if not group:continue
            states=r.get('states',[]);events=r.get('couples',[])
            if len(states)!=3 or len(events)!=2 or any(s.get('G_eV') is None for s in states):continue
            hashes={p:h for s in states for p,h in s.get('source_sha256',{}).items()}
            import hashlib
            if any(not Path(p).exists() or hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in hashes.items()):continue
            first,second=[e.get('E_native_vs_Fc_V') for e in events]
            if first is None or second is None:continue
            out.append(dict(id=r['id'],name=r['id'],family=group.get('family'),kind='prior_native',
                            E1=first,disp=(first-second)*load_config('electrolyte').FARADAY/1000,
                            lam=None,cap=None,provisional=False,withdrawn_spin=False))
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
        expansion=r['kind']=='expansion'
        matched=r['kind']=='prior_native'
        prov = r["provisional"] and y in ("disp", "lam")
        ax.scatter(r[x], r[y], s=520 if lead else 420, marker='s' if matched else 'o' if expansion else "D" if lead else "*",
                   facecolor="white" if lead else c, edgecolor='#D55E00' if expansion else "#999999" if prov else (c if lead else "black"),
                   lw=2.6 if lead else 1.2, ls="--" if prov else "-", zorder=6)
        if r.get('withdrawn_spin'):
            ax.scatter(r[x],r[y],s=220,marker='x',color='#666666',linewidth=2.2,zorder=8)
        texts.append(ax.text(r[x], r[y], f"{num}{'?' if prov else ''}", fontsize=18,
                             fontweight="bold", color="#999999" if prov else c, zorder=7))
        suffix=(' (withdrawn spin)' if r.get('withdrawn_spin') else ' (provisional)' if prov else '')
        if y=='lam':suffix+=f" [gas {r.get('gas_events','?')}/{r.get('declared_events','?')}]"
        key.append(textwrap.shorten(f"{num}{'?' if prov else ''}  {'[prior native] ' if matched else '[native] ' if expansion else '[earlier lead] ' if lead else ''}{r['name']}"+suffix,
                                   width=66,placeholder='…'))
    adjust_text(texts, ax=ax, expand=(1.6, 1.8), arrowprops=dict(arrowstyle="-", color="#555", lw=1.2))
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); grid_xy(ax)
    hs = [plt.Line2D([], [], ls="", marker="*", ms=20, mfc="#888", mec="black", label="previous candidates (historical table)"),
          plt.Line2D([], [], ls="", marker="D", ms=14, mfc="white", mec="#0072B2", mew=2.4,
                     label="earlier discovered leads (historical table)"),
          plt.Line2D([],[],ls='',marker='o',ms=13,mfc='#888',mec='#D55E00',mew=2.5,label='new expansion: admitted native results')]
    if any(r['kind']=='prior_native' and r[x] is not None and r[y] is not None for r in rows):
        hs.append(plt.Line2D([],[],ls='',marker='s',ms=12,mfc='#0072B2',mec='black',label='prior candidates: current matched native model'))
    hs += [plt.Line2D([], [], ls="", marker="s", ms=13, color=c, label=f) for f, c in FAM.items()]
    if any(r.get('withdrawn_spin') for r in rows):
        hs.append(plt.Line2D([],[],ls='',marker='x',ms=12,color='#666666',label='historical incorrect spin; excluded from ranking'))
    ax.legend(handles=hs, loc="best", fontsize=18)
    half = (len(key) + 1) // 2
    for col, chunk in enumerate((key[:half], key[half:])):
        fig.text(0.02 + 0.5 * col, 0.31, "\n".join(chunk), va="top", ha="left", fontsize=16)
    fig.text(0.02, 0.355, note, fontsize=18, va="top")
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / fname; fig.savefig(p);fig.savefig(p.with_suffix('.pdf'))
    renderer=fig.canvas.get_renderer()
    box=ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.03,1.06)
    fig.savefig(OUT/(p.stem+'_slide.png'),bbox_inches=box,dpi=220)
    plt.close(fig); print(" ->", p)


def main(feed=None):
    rows = gather(feed)
    inc = [r for r in rows if r["kind"] == "lead" and None in (r["E1"], r["disp"], r["lam"])]
    print("leads with a missing quantity:", [(r["id"], {k: r[k] is not None for k in ("E1", "disp", "lam")}) for r in inc])
    _scatter(rows, "E1", "disp", "first redox potential E1 (V vs Fc/Fc$^+$)",
             "computed ΔG$_{disp}$ of the intermediate (kJ mol$^{-1}$)",
             "leads_vs_existing_E1_disp.png",
             "Earlier points retain their historical reference/thermal treatment; native additions use matched solution Fc.\n"
             "Compare current native points within a family. Experimental uncertainty of the revised model is uncalibrated.")
    _scatter(rows, "cap", "lam", "specific capacity (mAh g$^{-1}$, computed accessible electrons)",
             "gas Nelsen sum (eV, QC-clean available couples, mean)", "leads_vs_existing_lambda_cap.png",
             "Gas diagnostics only: every new two-electron point currently covers one of its two events; no solution-rate claim.\n"
             "Native capacities use window/integrity bookkeeping; earlier capacities retain their historical table definition.")


if __name__ == "__main__":
    main()
