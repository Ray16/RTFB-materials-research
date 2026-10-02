#!/usr/bin/env python
"""One-slide, presentation-honest summary of the merrifield multi-electron candidates.
Panel A: multi-electron redox map (E1,E2 vs Fc inside the electrolyte window) -> all genuine 2e.
Panel B: inner-sphere lambda, grafted vs standalone; CLEAN values solid, conformer-flagged open
         (preliminary) so nothing is over-claimed.
Panel C: 2e specific capacity (bars) with SA score annotated.
"""
from __future__ import annotations
import csv
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, C  # noqa
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
R = ROOT / "results"
FIG = R / "figures" / "candidates"; FIG.mkdir(parents=True, exist_ok=True)

NAME = {
    "aq_benzyloxy":  "2-(benzyloxy)-\nanthraquinone",
    "aq_benzylamino":"2-(benzylamino)-\nanthraquinone",
    "nq_benzyloxy":  "2-(benzyloxy)-1,4-\nnaphthoquinone",
    "bisviologen":   "bis-benzyl\nviologen",
    "dtbc_phenol":   "mono-O-benzyl\nDTB-catechol (1e)",
}
# verified SA (RDKit sascorer) + 2e specific capacity (mAh/g) from the grafted monomer model
SA  = {"aq_benzyloxy":1.81,"aq_benzylamino":1.86,"nq_benzyloxy":1.99,"bisviologen":2.34,"dtbc_phenol":2.18}
CAP = {"aq_benzyloxy":163.2,"aq_benzylamino":163.7,"nq_benzyloxy":192.6,"bisviologen":146.3,"dtbc_phenol":82.1}
NE  = {"aq_benzyloxy":2,"aq_benzylamino":2,"nq_benzyloxy":2,"bisviologen":2,"dtbc_phenol":1}

def load_me():
    d={}
    for r in csv.DictReader(open(R/"multielectron_stability.csv")):
        d[r["id"]]=(float(r["E1_vs_Fc_V"]), float(r["E2_vs_Fc_V"]))
    return d

def load_lambda():
    """primary-couple lambda (grafted) + standalone, with clean/flagged status."""
    prim={"aq_benzyloxy":"neu->red1","aq_benzylamino":"neu->red1","nq_benzyloxy":"neu->red1",
          "bisviologen":"ox2->ox1","dtbc_phenol":"ox->neu"}
    graft={}; sa={}
    for r in csv.DictReader(open(R/"reorganization.csv")):
        i=r["id"]; c=r["couple"]
        try: lam=float(r["lambda_i_eV"])
        except: continue
        flag=bool(r.get("flag"))
        if i in prim and c==prim[i]: graft[i]=(lam,flag)
        base=i[:-3] if i.endswith("_sa") else None
        if base in prim and c==prim.get(base): sa[base]=(lam,flag)
    return graft, sa

apply_style()
me=load_me(); glam,slam=load_lambda()
order=["aq_benzyloxy","aq_benzylamino","nq_benzyloxy","bisviologen","dtbc_phenol"]
labels=[NAME[i] for i in order]
y=np.arange(len(order))[::-1]

fig,axes=plt.subplots(1,3,figsize=(25,8.2),gridspec_kw=dict(wspace=0.55))
QC="#0072B2"; VIO="#CC79A7"; GREY="#8A8A8A"
famc={"aq_benzyloxy":QC,"aq_benzylamino":QC,"nq_benzyloxy":QC,"bisviologen":VIO,"dtbc_phenol":"#762A83"}

# --- Panel A: multi-electron redox map ---
axA=axes[0]
axA.axvspan(-3.0,2.0,color="#1B7837",alpha=0.06,zorder=0)
axA.axvline(0,color=GREY,ls=":",lw=1.2)
for yi,i in zip(y,order):
    col=famc[i]
    if i in me:
        e1,e2=me[i]
        axA.plot([e2,e1],[yi,yi],"-",color=col,lw=2.5,zorder=2)
        axA.scatter([e1],[yi],s=130,color=col,zorder=3,edgecolor="white",linewidth=1.2)
        axA.scatter([e2],[yi],s=130,color=col,zorder=3,edgecolor="white",linewidth=1.2)
        axA.text(e1+0.05,yi+0.14,"$E_1$",fontsize=15,color=col)
        axA.text(e2-0.05,yi+0.14,"$E_2$",fontsize=15,color=col,ha="right")
        axA.text(min(e1,e2)-0.15,yi,"2e",fontsize=15,fontweight="bold",va="center",ha="right",color=col)
    else:  # dtbc = 1e
        axA.scatter([-0.0],[yi],s=0)
        axA.text(0.05,yi,"1e (p-type,\nout of 2e scope)",fontsize=13,va="center",color="#762A83")
axA.set_yticks(y); axA.set_yticklabels(labels,fontsize=15)
axA.set_xlabel("Redox potential (V vs Fc/Fc$^+$)")
axA.set_xlim(-2.6,0.6)
axA.set_title("Genuine 2e (both couples in-window)")

# --- Panel B: lambda grafted vs standalone ---
axB=axes[1]
for yi,i in zip(y,order):
    col=famc[i]
    if i in glam:
        lam,flag=glam[i]
        axB.scatter([lam],[yi+0.12],s=150,color=col,zorder=3,
                    edgecolor="white",linewidth=1.2,
                    marker="o" if not flag else "o",
                    facecolor=col if not flag else "none")
        if flag:
            axB.scatter([lam],[yi+0.12],s=150,facecolor="none",edgecolor=col,linewidth=2.2,zorder=3)
        axB.text(lam,yi+0.30,f"{lam:.2f}",fontsize=15,ha="center",color=col)
    if i in slam:
        lam,flag=slam[i]
        if flag:   # flagged standalone -> open diamond (preliminary), don't imply clean
            axB.scatter([lam],[yi-0.14],s=120,facecolor="none",edgecolor=GREY,marker="D",zorder=3,linewidth=2.0)
        else:
            axB.scatter([lam],[yi-0.14],s=110,color=GREY,marker="D",zorder=3,edgecolor="white",linewidth=1)
        axB.text(lam,yi-0.36,f"{lam:.2f}",fontsize=15,ha="center",color=GREY)
axB.set_yticks(y); axB.set_yticklabels(labels,fontsize=15)
axB.set_xlabel("Inner-sphere $\\lambda_i$ (eV)  — lower = faster kinetics")
axB.set_xlim(0.2,1.25)
axB.set_title("Reorganization energy $\\lambda_i$ (primary couple)",fontsize=19)
h=[plt.Line2D([],[],marker="o",color="#0072B2",ls="",ms=11,label="grafted (filled = clean)"),
   plt.Line2D([],[],marker="o",color="#0072B2",ls="",ms=11,mfc="none",mew=2,label="grafted (open = prelim, conformer-flagged)"),
   plt.Line2D([],[],marker="D",color=GREY,ls="",ms=10,label="standalone core (_sa)")]
axB.legend(handles=h,loc="upper center",bbox_to_anchor=(0.5,-0.17),fontsize=14,frameon=False)
axB.set_ylim(min(y)-0.75, max(y)+0.75)

# --- Panel C: capacity + SA ---
axC=axes[2]
caps=[CAP[i] for i in order]
bars=axC.barh(y,caps,color=[famc[i] for i in order],alpha=0.85,height=0.6,edgecolor="white")
for yi,i in zip(y,order):
    axC.text(CAP[i]+3,yi,f"{CAP[i]:.0f}",va="center",fontsize=15,fontweight="bold",color="#333")
    axC.text(4,yi,f"SA {SA[i]:.2f}  ({NE[i]}e)",va="center",fontsize=15,color="white",fontweight="bold")
axC.set_yticks(y); axC.set_yticklabels(labels,fontsize=15)
axC.set_xlabel("Specific capacity (mAh g$^{-1}$, grafted monomer)")
axC.set_xlim(0,215)
axC.set_title("Capacity & synthetic accessibility",fontsize=19)

fig.suptitle("Merrifield-grafted multi-electron redox candidates — DFT+SMD(MeCN), wB97M-V/def2-TZVPD (uniform)",
             fontsize=21,fontweight="bold",y=1.03)
fig.text(0.5,-0.2,
         "Grafted $\\lambda_i$ marked 'prelim' are inflated by flexible-tether rotation in the isolated model "
         "(conformer-matching in progress;\nthe real polymer tether is backbone-anchored). "
         "Standalone (_sa) cores are the clean reference.",
         ha="center",fontsize=15,style="italic",color="#555")
fig.tight_layout()
out=FIG/"merrifield_summary.png"
fig.savefig(out,bbox_inches="tight")
print("wrote",out)
