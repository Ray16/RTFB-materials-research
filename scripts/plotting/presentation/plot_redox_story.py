#!/usr/bin/env python3
"""Presentation-sized charts from the frozen, checked discussion evidence."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
EVIDENCE=ROOT/'results/reviews/2026-10-26/colleague_workflow/evidence.json'
OUT=ROOT/'results/figures/presentation/story'
D=json.loads(EVIDENCE.read_text())
INK='#15343B';TEAL='#007F86';GOLD='#BD6928'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':18,'axes.labelsize':20,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
outputs={}
def save(fig,name):
 for ext in ['pdf','png']:
  p=OUT/(name+'.'+ext);fig.savefig(p,bbox_inches='tight',dpi=180,facecolor='white');outputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
 plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.2),layout='constrained')
rows=[r for r in D['systematic_matches'] if r['core']=='pmdi']
contrasts=['ring_F','ring_CN','cap_length','linker_C2'];labels=['Ring F','Ring CN','Ethyl → propyl cap','Direct → ether linker']
values=[next(r['delta_E1_V'] for r in rows if r['contrast']==k) for k in contrasts]
ax.barh(range(4),values,color=[TEAL,TEAL,GOLD,GOLD],height=.52)
for y,v in enumerate(values):ax.text(v+(.008 if v>0 else -.008),y,f'{v:+.3f}',va='center',ha='left' if v>0 else 'right',fontsize=22,weight='bold')
ax.set_yticks(range(4),labels);ax.invert_yaxis();ax.axvline(0,color=INK,lw=1);ax.set_xlim(-.14,.25);ax.set_xlabel('First-potential shift vs matched PMDI baseline (V)');ax.grid(axis='x',alpha=.12);ax.set_axisbelow(True)
save(fig,'matched_pmdi')

fig,ax=plt.subplots(figsize=(6.5,5.5),layout='constrained')
for wave,color,marker,label in [(1,TEAL,'o','First electron'),(2,GOLD,'s','Second electron')]:
 rows=[r for r in D['benchmark']['points'] if r['wave']==wave]
 ax.scatter([r['exp_V'] for r in rows],[r['calc_V'] for r in rows],color=color,marker=marker,s=82,label=label,edgecolor='white',linewidth=.7,zorder=3)
ax.plot([-2.4,1.6],[-2.4,1.6],ls='--',color=INK,lw=1)
ax.set(xlim=(-2.4,1.6),ylim=(-2.4,1.6),xlabel='Measured vs Fc/Fc⁺ (V)',ylabel='Computed vs Fc/Fc⁺ (V)');ax.set_aspect('equal');ax.legend(loc='upper left',fontsize=16,frameon=False);ax.grid(alpha=.12)
save(fig,'redox_validation')

fig,ax=plt.subplots(figsize=(10,4.6),layout='constrained')
for y,gid,label in [(0,'exp_bq_oxy_methoxy','BQ–OMe'),(1,'exp_bq_oxy_trimethyl','BQ–Me₃')]:
 row=next(r for r in D['base_assisted']['candidates'] if r['id']==gid)
 a=row['reference_deltaG_kJ_mol'];b=row['base_models']['dbu']['free_ions_kJ_mol']
 ax.plot([b,a],[y,y],color='#B9C4C7',lw=3,zorder=1)
 for value,color,marker in [(a,GOLD,'o'),(b,TEAL,'s')]:
  ax.scatter(value,y,c=color,marker=marker,s=140,zorder=3);ax.text(value,y-.18,f'{value:+.1f}',ha='center',fontsize=23,weight='bold',color=color)
ax.axvline(0,color=INK,ls=':',lw=1.4);ax.set_yticks([0,1],['BQ–OMe','BQ–Me₃']);ax.set_ylim(1.55,-.6);ax.set_xlim(-70,82);ax.set_xlabel('Reaction ΔG in MeCN (kJ/mol)');ax.text(-49,-.48,'With DBU: separated ions',ha='center',fontsize=17,color=TEAL);ax.text(52,-.48,'HCl reference',ha='center',fontsize=17,color=GOLD);ax.grid(axis='x',alpha=.12)
save(fig,'base_bq')
(OUT/'manifest.json').write_text(json.dumps({'evidence':str(EVIDENCE.relative_to(ROOT)),'evidence_sha256':hashlib.sha256(EVIDENCE.read_bytes()).hexdigest(),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'outputs':outputs,'quantum_jobs_submitted':0},indent=2)+'\n')
