#!/usr/bin/env python
"""Show raw step-dependent errors and keep internal offset diagnostics separate."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.chunked import file_hash


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis',type=Path,required=True)
    parser.add_argument('--out',type=Path,default=ROOT/'results/figures/validation')
    args=parser.parse_args();source=json.loads((args.analysis/'summary.json').read_text())
    rows=list(csv.DictReader((args.analysis/'points.csv').open()))
    if source['primary_potentials_modified'] or source['external_validation_claimed']:
        raise ValueError('This figure is a post hoc diagnostic, not a replacement benchmark')
    sources={str(p):file_hash(p) for p in [Path(__file__),args.analysis/'summary.json',args.analysis/'points.csv']}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,
        'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,(left,right)=plt.subplots(1,2,figsize=(7.2,3.2));fig.subplots_adjust(left=.09,right=.98,bottom=.25,top=.87,wspace=.45)
    for step,color in [(1,'#0072B2'),(2,'#D55E00')]:
        selected=sorted([r for r in rows if int(r['wave'])==step],key=lambda r:r['id'])
        values=np.array([float(r['error_V']) for r in selected]);x=step+np.linspace(-.16,.16,len(values))
        left.scatter(x,values,color=color,s=16,edgecolor='white',linewidth=.3)
        left.plot([step-.24,step+.24],[values.mean()]*2,color='black',lw=1)
        left.text(step,.61,f'n={len(values)}',ha='center',fontsize=7)
    left.axhline(0,color='.6',lw=.7);left.set(xticks=[1,2],xticklabels=['First step','Second step'],
        ylabel='Calculated − measured potential (V)',ylim=(-.6,.7),xlim=(.5,2.5),title='Unfitted residuals')
    models=[source['primary_unfitted']]+[next(r['metrics'] for r in source['cross_validation']
        if r['holdout']=='id' and r['model']==model) for model in ('constant_offset','separate_step_offsets')]
    x=np.arange(3)
    for key,color,shift,label in [('MAE_V','#0072B2',-.18,'MAE'),('RMSE_V','#E69F00',.18,'RMSE')]:
        values=[r[key] for r in models];right.bar(x+shift,values,width=.32,color=color,label=label)
        for i,v in enumerate(values):right.text(i+shift,v+.008,f'{v:.3f}',ha='center',fontsize=6.5)
    right.set(xticks=x,xticklabels=['Unfitted\n(primary)','Molecule CV\none offset','Molecule CV\ntwo step offsets'],
        ylabel='Potential error (V)',ylim=(0,.36),title='Internal calibration diagnostics')
    right.legend(frameon=False,fontsize=7,loc='upper right')
    for axis,letter in zip((left,right),('a','b')):axis.text(-.14,1.06,letter,transform=axis.transAxes,fontweight='bold')
    args.out.mkdir(parents=True,exist_ok=True);stem=args.out/'reference_bias';outputs=[]
    for suffix in ('pdf','png','svg'):
        p=stem.with_suffix('.'+suffix);fig.savefig(p,dpi=600);outputs.append(p)
    plt.close(fig)
    if any(file_hash(p)!=v for p,v in sources.items()):raise ValueError('Figure source changed')
    atomic_write_json(stem.with_suffix('.json'),dict(status='reference_bias_diagnostic_figure',source_sha256=sources,
        output_sha256={str(p):file_hash(p) for p in outputs},primary_model_changed=False,external_validation=False))


if __name__=='__main__':main()
