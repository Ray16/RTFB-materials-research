#!/usr/bin/env python
"""Maintained base-assisted sensitivity figure; missing results stay visibly pending."""
import hashlib
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from redox.core.common import ROOT
from redox.core.io import atomic_write_json

OUT=ROOT/'results/grafting/base_assisted'
LABELS={'exp_bq_oxy_methoxy':'BQ · OMe / O tether','exp_bq_oxy_trimethyl':'BQ · 3Me / O tether',
        'exp_aq_1oxy':'AQ · 1-O tether','exp_aq_2oxy_6methoxy':'AQ · 2-O / 6-OMe',
        'design_ndi_ethyl':'NDI · ethyl','design_pmdi_ethyl':'PMDI · ethyl'}


def main():
    payload=json.loads((OUT/'current.json').read_text());rows=payload['candidates'];status=payload['status']
    figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(8.2,4.5))
    for i,row in enumerate(rows):
        ref=row['reference_deltaG_kJ_mol']
        if ref is not None:ax.scatter(ref,i,s=48,facecolors='white',edgecolors='#777777',zorder=3)
        else:ax.text(.98,i,'Reference held',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=8,color='#777777')
        for b,offset,color in [('dbu',-.16,'#0072B2'),('dipea',.16,'#D55E00')]:
            model=row['base_models'][b];a=model['free_ions_kJ_mol'];c=model['contact_pair_kJ_mol']
            if a is not None and c is not None:ax.plot([a,c],[i+offset]*2,color=color,lw=1.5,alpha=.65)
            for value,marker in [(a,'o'),(c,'s')]:
                if value is not None:ax.scatter(value,i+offset,s=42,marker=marker,color=color,zorder=4)
    ax.axvline(0,color='#444444',lw=.9)
    ax.set_yticks(range(len(rows)),[LABELS.get(r['id'],r['name']) for r in rows]);ax.invert_yaxis()
    ax.set_xlabel('Grafting reaction ΔG° (kJ mol⁻¹) · MeCN, 298 K')
    ax.set_title('Base-assisted grafting: thermodynamic sensitivity',loc='left',fontsize=12,pad=17)
    for spine in ('top','right','left'):ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y',length=0);ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
    legend=[Line2D([],[],marker='o',markerfacecolor='white',markeredgecolor='#777777',color='none',label='HCl reference'),
            Line2D([],[],color='#0072B2',lw=2,label='DBU'),Line2D([],[],color='#D55E00',lw=2,label='DIPEA'),
            Line2D([],[],marker='o',color='#444444',ls='none',label='Separated ions'),
            Line2D([],[],marker='s',color='#444444',ls='none',label='Contact pair')]
    fig.legend(handles=legend,loc='lower center',ncol=3,frameon=False,fontsize=8,bbox_to_anchor=(.5,.01))
    fig.tight_layout(rect=(0,.15,1,1))
    for ext in ('pdf','svg','png'):fig.savefig(figdir/f'base_comparison.{ext}',dpi=300,bbox_inches='tight',metadata={'Creator':'redox pipeline'} if ext=='pdf' else None)
    plt.close(fig)
    def value(v):return 'Pending' if v is None else f'{v:+.1f}'
    table=''
    for r in rows:
        cells=[LABELS.get(r['id'],r['name']),value(r['reference_deltaG_kJ_mol'])]
        for b in ('dbu','dipea'):
            cells += [value(r['base_models'][b]['free_ions_kJ_mol']),value(r['base_models'][b]['contact_pair_kJ_mol'])]
        table+='<tr>'+''.join('<td>'+html.escape(s)+'</td>' for s in cells)+'</tr>'
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Base-assisted grafting pilot</title><style>body{{font:16px system-ui;max-width:1050px;margin:36px auto;padding:0 20px;color:#222}}img{{max-width:100%}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:10px;border-bottom:1px solid #ddd;text-align:left}}.table{{overflow-x:auto}}.note{{background:#fff4df;padding:16px}}a{{color:#0072B2}}</style>
<p><a href="../index.html">Grafting reference results</a> · <a href="../../../docs/PRACTICAL_GRAFTING.md">Practical chemistry</a></p>
<h1>Base-assisted grafting pilot</h1><p><b>{status['components_qualified']}/6 components qualified</b>; {status['failures']} preserved failures. Reload for updates.</p>
<p class="note">Exploratory MeCN thermodynamics. No candidate is cleared for synthesis by this pilot. Reaction barriers, competing base alkylation, salt precipitation and experimental yields remain unestablished.</p>
<p><b>Base required; BQ candidates retained for evaluation.</b> Both BQ designs and AQ-1-O have favorable modeled DBU/DIPEA reactions. Positive HCl-reference values alone do not exclude them. AQ-2-O/6-OMe remains unresolved because of precursor QC.</p>
<p>Each base is modeled as the neutral reagent, protonated ion and chloride contact pair. Connected endpoints show two ionic models, not statistical error bars. The HCl-reference results remain separate.</p>
<img src="figures/base_comparison.svg" alt="Reference and base-assisted grafting energies, with uncomputed values omitted">
<p><a href="figures/base_comparison.pdf">Vector PDF</a> · <a href="current.json">Data and source receipts</a> · <a href="status.json">Status</a></p>
<div class="table"><table><thead><tr><th>Candidate</th><th>HCl reference</th><th>DBU ions</th><th>DBU pair</th><th>DIPEA ions</th><th>DIPEA pair</th></tr></thead><tbody>{table}</tbody></table></div>
<p>All values: kJ/mol, 298.15 K. Explicit 1-atm → 1-M correction is applied when the number of independently translating solutes changes. Contact-pair entropy and conformer sampling are approximate.</p>
<p>Precedents: <a href="https://doi.org/10.1002/pat.1968">AQ carboxylate/DBU in NMP</a>; <a href="https://doi.org/10.1002/chem.202200149">ferrocene carboxylate/DIPEA in DMF</a>. Applying these bases to our exact phenolic candidates in MeCN is a proposed test.</p></html>'''
    (OUT/'index.html').write_text(page)
    files=list(figdir.glob('base_comparison.*'))+[OUT/'index.html']
    atomic_write_json(figdir/'manifest.json',dict(source_sha256=hashlib.sha256((OUT/'current.json').read_bytes()).hexdigest(),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}))


if __name__=='__main__':main()
