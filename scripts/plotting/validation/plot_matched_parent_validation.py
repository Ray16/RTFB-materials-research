#!/usr/bin/env python
"""Maintain an actual-potential parity/error figure from the audited parent snapshot."""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
LABELS={
'bm_1_benzyl_1_methyl_4_4_bipyridinium':'Benzyl–methyl viologen',
'bm_n_n_diisopropylpyromellitic_diimide':'PMDI / isopropyl',
'bm_n_n_bis_2_methoxyethyl_pyromellitic_diimide':'PMDI / methoxyethyl',
'bm_2_6_dimethyl_1_4_benzoquinone':'Dimethyl benzoquinone',
'bm_n_ethylphenothiazine':'N-ethylphenothiazine',
'bm_n_2_2_methoxyethoxy_ethyl_phenothiazine':'MEEPT',
'bm_n_ethyl_3_7_dimethoxyphenothiazine':'DMeOEPT',
'bm_3_7_dimethoxy_n_2_2_methoxyethoxy_ethyl_phenothiazine':'DMeOMEEPT'}

def step(row):
    first='ox1->neu' if 'phenothiazine' in row['family'] else ('ox2->ox1' if 'viologen' in row['family'] else 'neu->red1')
    return 1 if row['event']==first else 2


def load_snapshot(folder):
    summary=json.loads((folder/'summary.json').read_text())
    assert summary['status']=='checked' and summary['arithmetic_discrepancies']==0 and not summary['fitted_correction']
    with (folder/'comparisons.csv').open() as stream:all_rows=list(csv.DictReader(stream))
    rows=[]
    for source in all_rows:
        if source['tier']!='A' or source['kind']!='E':continue
        row=dict(source,step=step(source),label=LABELS[source['id']])
        for key in ('exp_V','calc_V','error_V'):row[key]=float(row[key]);assert math.isfinite(row[key])
        assert abs(row['calc_V']-row['exp_V']-row['error_V'])<1e-12
        rows.append(row)
    errors=np.array([r['error_V'] for r in rows]);stats=dict(n=len(rows),MAE=float(np.mean(abs(errors))),RMSE=float(np.sqrt(np.mean(errors**2))),bias=float(np.mean(errors)))
    for key,value in stats.items():assert abs(value-summary['tier_A_potential_stats'][key])<1e-12,key
    assert len(rows)==summary['tier_A_potential_events']==12 and len({r['id'] for r in rows})==8
    assert len({r['doi'] for r in rows})==summary['tier_A_original_studies']==4
    return rows,stats,summary


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot',type=Path,default=ROOT/'results/reviews/2026-10-19/full_matched_parent_final');ap.add_argument('--out',type=Path,default=ROOT/'results/figures/validation');args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    rows,stats,summary=load_snapshot(args.snapshot)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    sources={str(p):sha(p) for p in (args.snapshot/'comparisons.csv',args.snapshot/'summary.json',Path(__file__))}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcdefaults();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(7.2,4.2));left=fig.add_axes([.095,.22,.315,.59]);right=fig.add_axes([.69,.22,.285,.59])
    colors={1:'#0072B2',2:'#D55E00'};markers={1:'o',2:'D'}
    order=[gid for gid in LABELS if any(r['id']==gid for r in rows)]
    left.plot([-2.3,1.5],[-2.3,1.5],color='#777777',lw=.9,zorder=1)
    right.axvline(0,color='#777777',lw=.9,zorder=1)
    for row in rows:
        s=row['step'];left.scatter(row['exp_V'],row['calc_V'],c=colors[s],marker=markers[s],s=33,edgecolors='white',linewidths=.45,zorder=3)
        y=order.index(row['id'])+(-.11 if s==1 else .11)
        right.scatter(row['error_V'],y,c=colors[s],marker=markers[s],s=33,edgecolors='white',linewidths=.45,zorder=3)
        if row['id']=='bm_n_ethyl_3_7_dimethoxyphenothiazine' and s==2:
            right.annotate(f'{row["error_V"]:+.3f}',xy=(row['error_V'],y),xytext=(4,1),textcoords='offset points',fontsize=7,color=colors[s],va='center')
    left.set(xlim=(-2.3,1.5),ylim=(-2.3,1.5),xlabel='Measured potential (V vs Fc/Fc⁺)',ylabel='Calculated potential (V vs Fc/Fc⁺)',title='Measured versus calculated')
    left.set_aspect('equal',adjustable='box');left.grid(alpha=.12)
    left.text(.04,.96,f'12 potentials · 8 parents\nMAE  {stats["MAE"]:.3f} V\nRMSE  {stats["RMSE"]:.3f} V',transform=left.transAxes,va='top',fontsize=7.5,bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=3))
    right.set(xlim=(-.27,.72),ylim=(len(order)-.55,-.55),yticks=list(range(len(order))),yticklabels=[LABELS[g] for g in order],xlabel='Calculated − measured (V)',title='Errors by parent')
    right.set_xticks([-.2,0,.2,.4,.6]);right.tick_params(axis='y',length=0);right.grid(axis='x',alpha=.12)
    for axis,letter in [(left,'a'),(right,'b')]:axis.text(-.15,1.06,letter,transform=axis.transAxes,fontweight='bold',fontsize=11)
    handles=[Line2D([],[],marker=markers[s],color=colors[s],linestyle='',label=f'{s}'+('st' if s==1 else 'nd')+' redox step',markersize=6) for s in (1,2)]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.97),ncol=2,frameon=False)
    fig.text(.095,.075,'No fitted correction. DMeOEPT second step: measured 0.655 V; calculated 1.202 V.',fontsize=7)
    fig.text(.095,.035,'Selected free-parent benchmark in MeCN; exact grafts, polymer behavior and kinetics are outside this validation.',fontsize=6.8)
    stem=args.out/'matched_parent_validation'
    for suffix in ('png','pdf','svg'):fig.savefig(stem.with_suffix('.'+suffix),dpi=600)
    plt.close(fig)
    columns=['id','label','step','event','exp_V','calc_V','error_V','doi','spin_qc']
    with stem.with_suffix('.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    tex=[r'% Actual potentials from the audited 12-event parent snapshot.',r'\begin{frame}{Validation: actual potentials and model errors}',r'\centering',r'\includegraphics[width=.98\linewidth,height=.71\textheight,keepaspectratio]{../results/figures/validation/matched_parent_validation.pdf}',r'\par\smallskip\textcolor{teal}{\textbf{Observed MAE '+f'{stats["MAE"]:.2f}'+r' V; RMSE '+f'{stats["RMSE"]:.2f}'+r' V. Largest error: '+f'{max(r["error_V"] for r in rows):+.2f}'+r' V.}}',r'\par\smallskip\footnotesize Twelve free-parent potentials; exact-graft accuracy and kinetics remain unvalidated.',r'\end{frame}',r'\begin{frame}{Measured and calculated parent potentials}',r'\small\centering\renewcommand{\arraystretch}{1.2}',r'\begin{tabular}{@{}lrrrr@{}}\toprule',r'Parent & \multicolumn{2}{c}{First redox step} & \multicolumn{2}{c}{Second redox step} \\',r'& Measured & Calculated & Measured & Calculated \\ \midrule']
    for gid in order:
        values={r['step']:r for r in rows if r['id']==gid};numbers=[]
        for s in (1,2):
            for field in ('exp_V','calc_V'):
                value=values.get(s,{}).get(field);numbers.append('---' if value is None else (f'{value:.4f}' if gid=='bm_2_6_dimethyl_1_4_benzoquinone' and field=='exp_V' else f'{value:.3f}'))
        tex.append(LABELS[gid].replace('–','--')+' & '+' & '.join(numbers)+r' \\')
    tex.extend([r'\bottomrule\end{tabular}',r'\par\medskip\raggedright\footnotesize All values V vs Fc/Fc$^+$ in MeCN. First/second refer to successive reductions for imides, quinone and viologen; successive oxidations for phenothiazines.',r'\par\smallskip Dashes: no admitted measurement for that step; excluded waves are not plotted. MAE 0.209 V; RMSE 0.246 V. No correction fitted.',r'\vfill\scriptsize Twelve potentials, eight exact free parents, four original studies. Source snapshot: 19 October 2026.',r'\end{frame}'])
    content='\n'.join(tex)+'\n'
    include=ROOT/'slides/generated/matched_parent_validation_figure.tex';include.write_text(content)
    split=content.index(r'\begin{frame}{Measured and calculated parent potentials}')
    (include.parent/'matched_parent_validation_plot.tex').write_text(content[:split])
    (include.parent/'matched_parent_actual_values.tex').write_text(content[split:])
    assert all(sha(p)==h for p,h in sources.items())
    files=[p for p in args.out.glob('matched_parent_validation.*') if p.suffix!='.json']
    files.extend(include.parent/name for name in ('matched_parent_validation_figure.tex','matched_parent_validation_plot.tex','matched_parent_actual_values.tex'))
    record=dict(status='checked_12_event_parent_visualization',n=12,parents=8,studies=4,stats=stats,source_snapshot=str(args.snapshot),source_sha256=sources,files={str(p.relative_to(ROOT)):sha(p) for p in files},scope=summary['scope'],fitted_correction=False)
    stem.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ('status','n','parents','stats')}))


if __name__=='__main__':main()
