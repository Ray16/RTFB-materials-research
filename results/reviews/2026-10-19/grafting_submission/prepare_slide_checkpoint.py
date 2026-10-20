"""Freeze a consistent initial grafting-results checkpoint for the Backup slides."""
from pathlib import Path
import json
import shutil
import fcntl
from datetime import datetime
from zoneinfo import ZoneInfo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[4]
out=Path(__file__).resolve().parent
snap=out/'first_results'
if not snap.exists():
    snap.mkdir()
    with (root/'results/grafting/.collection.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        for name in ['current.json','current.csv','status.json']:
            shutil.copyfile(root/'results/grafting'/name,snap/name)
        for name in ['deltaG_ranking.pdf','coverage.pdf','reaction_schemes.pdf']:
            shutil.copyfile(root/'results/grafting/figures'/name,snap/name)
feed=json.loads((snap/'current.json').read_text())
ready=[r for r in feed['candidates'] if r['deltaG_reference_eV'] is not None]
preferred=['exp_phthalimide_control','hyp_cata_ndi_ethyl','design_pmdi_ethyl','design_ndi_ethyl','exp_bq_oxy_trimethyl','exp_aq_1oxy']
selection=[next(r for r in ready if r['id']==key) for key in preferred if any(r['id']==key for r in ready)]
aliases={'design_pmdi_ethyl':'PMDI / ethyl','exp_aq_1oxy':'AQ / 1-O tether','exp_bq_oxy_trimethyl':'BQ / 3Me','exp_bq_oxy_methoxy':'BQ / OMe','exp_phthalimide_control':'Phthalimide control','design_ndi_ethyl':'NDI / ethyl','hyp_cata_ndi_ethyl':'Cata-NDI / ethyl'}
plt.rcParams.update({'font.size':13,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(7.7,3.25),layout='constrained')
vals=[r['deltaG_reference_kJ_mol'] for r in selection]
ax.scatter(vals,range(len(vals)),s=65,color='#0072B2',zorder=3)
ax.set_yticks(range(len(vals)),[aliases.get(r['id'],r['id']) for r in selection])
ax.invert_yaxis();ax.axvline(0,color='.6',lw=1)
for i,value in enumerate(vals):
    ax.annotate(f'{value:+.1f}',(value,i),xytext=(7,3),textcoords='offset points',fontsize=12)
span=max(vals)-min(vals)
ax.set_xlim(min(vals)-max(5,span*.08),max(vals)+max(8,span*.2))
ax.set_xlabel('Reference grafting ΔG (kJ mol⁻¹)');ax.tick_params(axis='y',length=0)
fig.savefig(snap/'grafting_slide.pdf');fig.savefig(snap/'grafting_slide.png',dpi=250)
plt.close(fig)
when=datetime.fromtimestamp(feed['status']['updated_unix'],ZoneInfo('America/Chicago')).strftime('%Y-%m-%d %H:%M Chicago')
frame=r'''\begin{frame}{Grafting: define the reaction before ranking candidates}
\large
\begin{block}{Local attachment models in acetonitrile, 298.15 K}
N--H / O--H:\quad $\mathrm{site{-}Cl + Nu{-}H \rightarrow site{-}Nu + HCl}$\\[0.8em]
N quaternization:\quad $\mathrm{site{-}Cl + Nu \rightarrow site{-}Nu^{+} + Cl^{-}}$
\end{block}
\normalsize Compare within each reaction class. An explicit base contribution is still pending.
\par\medskip\footnotesize These are capped-site thermodynamic models. Activation barriers, yields, site accessibility and polymer interactions are separate questions.
\end{frame}
\begin{frame}{Grafting thermodynamics: the first computed results}
\centering
\includegraphics[width=.97\linewidth,height=.68\textheight,keepaspectratio]{../results/reviews/2026-10-19/grafting_submission/first_results/grafting_slide.pdf}
\par\smallskip\footnotesize N--H/O--H reference: $\mathrm{site{-}Cl + Nu{-}H \rightarrow site{-}Nu + HCl}$. Base capture is uncomputed; values do not predict grafting rates or yields.
\par\smallskip\scriptsize CHECKPOINT
\end{frame}
'''.replace('CHECKPOINT',f'{len(ready)}/75 reactions qualified at {when}; continued calculations update the visual report.')
(root/'slides/generated/grafting_reference_2026-10-19.tex').write_text(frame)
print('Frozen first-result checkpoint:',len(ready),'qualified reactions; slide subset',len(selection))
