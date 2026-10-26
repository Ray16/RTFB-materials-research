#!/usr/bin/env python
"""Render discussion figures from one checked, frozen evidence JSON.

Presentation figures are kept separate from manuscript and scientific outputs.
The script does not change any calculation, candidate, policy or result pointer.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

INK='#15304A';TEAL='#087F8C';AMBER='#D58A17';BLUE='#4875B5'
FAMILIES={'quinone':'#C57328','aromatic imide':TEAL,'viologen':BLUE,
          'phenothiazine':'#A26396','nitroxide':'#778638'}


def style():
    plt.rcParams.update({'font.size':13,'axes.titlesize':16,'axes.labelsize':14,
        'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,
        'text.color':INK,'xtick.color':INK,'ytick.color':INK,
        'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none',
        'legend.frameon':False,'axes.titleweight':'bold','axes.edgecolor':'#B0BAC2'})


def library(d):
    fig,ax=plt.subplots(figsize=(10.7,3.6),layout='constrained')
    names=['Viologen','NDI','PMDI'];counts=[d['enumeration']['families'][k] for k in ['viologen','ndi','pmdi']]
    bars=ax.barh(names[::-1],counts[::-1],height=.6,color=[TEAL,BLUE,AMBER][::-1])
    for b in bars:ax.text(b.get_width()+9,b.get_y()+b.get_height()/2,f'{b.get_width():.0f}',va='center',fontsize=20,weight='bold')
    ax.set_xlim(0,540);ax.set_xlabel('Distinct exact local graft-site designs');ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    return fig


def structures(d):
    fig,axes=plt.subplots(1,3,figsize=(11.5,3.45),layout='constrained')
    for ax,row,label in zip(axes,d['exact_examples'],['Baseline PMDI','Ring F','Ring CN']):
        mol=Chem.MolFromSmiles(row['exact_grafted_smiles'])
        draw=rdMolDraw2D.MolDraw2DCairo(900,650)
        opts=draw.drawOptions();opts.bondLineWidth=3;opts.baseFontSize=.9;opts.padding=.08
        highlight=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==9]
        if label=='Ring CN':
            highlight=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==7 and any(b.GetBondType()==Chem.BondType.TRIPLE for b in a.GetBonds())]
        draw.DrawMolecule(mol,highlightAtoms=highlight);draw.FinishDrawing()
        ax.imshow(Image.open(io.BytesIO(draw.GetDrawingText())));ax.axis('off');ax.set_title(label)
    return fig


def pmdi_model(d):
    row=d['exact_examples'][0]
    fig,ax=plt.subplots(figsize=(5,4.4),layout='constrained')
    mol=Chem.MolFromSmiles(row['exact_grafted_smiles'])
    draw=rdMolDraw2D.MolDraw2DCairo(1000,850)
    draw.drawOptions().bondLineWidth=3;draw.drawOptions().padding=.05
    draw.DrawMolecule(mol);draw.FinishDrawing()
    ax.imshow(Image.open(io.BytesIO(draw.GetDrawingText())));ax.axis('off')
    return fig


def matched(d):
    contrasts=['cap_length','ring_F','ring_CN','ring_OMe','ring_Me','linker_C2']
    labels=['Propyl cap','Ring F','Ring CN','Ring OMe','Ring Me','Ether linker']
    rows={(r['core'],r['contrast']):r for r in d['systematic_matches']}
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),sharex=True,sharey=True,layout='constrained')
    values=[r[f'delta_E{w}_V'] for r in rows.values() for w in [1,2]]
    bound=max(abs(x) for x in values)*1.15
    for ax,core,title in zip(axes,['viologen','pmdi','ndi'],['Viologen','PMDI','NDI']):
        for y,c in enumerate(contrasts):
            r=rows[core,c];a=r['delta_E1_V'];b=r['delta_E2_V']
            ax.plot([a,b],[y,y],color='#AEB9C3',lw=2)
            ax.scatter(a,y,s=70,c=TEAL,marker='o',zorder=3,label='$E_1$' if y==0 else None)
            ax.scatter(b,y,s=65,c=AMBER,marker='s',zorder=3,label='$E_2$' if y==0 else None)
        ax.axvline(0,c=INK,lw=.9,ls=':');ax.set_title(title);ax.set_xlim(-bound,bound)
        ax.set_yticks(range(6),labels);ax.grid(axis='x',alpha=.12);ax.set_xlabel('Shift vs matched baseline (V)')
    axes[0].invert_yaxis();axes[-1].legend(loc='lower right',ncol=2)
    return fig


def benchmark(d):
    b=d['benchmark'];fig,(ax,bar)=plt.subplots(1,2,figsize=(11,4.3),gridspec_kw={'width_ratios':[1.15,1]},layout='constrained')
    vals=[r[k] for r in b['points'] for k in ['exp_V','calc_V']];lo=min(vals)-.12;hi=max(vals)+.12
    ax.plot([lo,hi],[lo,hi],c=INK,ls='--',lw=1)
    for family,color in FAMILIES.items():
        for wave,marker in [(1,'o'),(2,'s')]:
            rows=[r for r in b['points'] if r['family']==family and r['wave']==wave]
            ax.scatter([r['exp_V'] for r in rows],[r['calc_V'] for r in rows],c=color,marker=marker,
                       s=65,edgecolors='white',linewidths=.7)
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Measured potential vs Fc (V)',ylabel='Computed potential vs Fc (V)')
    ax.set_aspect('equal');ax.grid(alpha=.12);ax.set_title('39 exact-parent endpoints')
    names=['quinone','aromatic imide','phenothiazine','viologen','nitroxide']
    ys=np.arange(len(names));mae=[b['by_family'][k]['MAE_V'] for k in names]
    bar.barh(ys,mae,color=[FAMILIES[k] for k in names],height=.6)
    bar.set_yticks(ys,[f'{k.title()} (n={b["by_family"][k]["n"]})' for k in names]);bar.invert_yaxis()
    bar.axvline(b['metrics']['MAE_V'],ls=':',c=INK);bar.set_xlim(0,.36);bar.set_xlabel('Unfitted MAE (V)')
    bar.set_title('Family coverage is uneven');bar.grid(axis='x',alpha=.12);bar.set_axisbelow(True)
    return fig


def bias(d):
    b=d['benchmark'];fig,(ax,bar)=plt.subplots(1,2,figsize=(11.2,4.1),layout='constrained')
    for wave,color in [(1,TEAL),(2,AMBER)]:
        rows=[r for r in b['points'] if r['wave']==wave]
        jitter=np.linspace(-.15,.15,len(rows))
        errors=[r['calc_V']-r['exp_V'] for r in rows]
        ax.scatter(wave+jitter,errors,c=color,s=45,alpha=.85)
        mean=b['by_wave'][str(wave)]['bias_V'];ax.plot([wave-.28,wave+.28],[mean,mean],c=INK,lw=2)
        ax.text(wave+.30,mean,f'{mean:+.3f} V',va='center',fontsize=11)
    ax.axhline(0,c=INK,lw=.8,ls=':');ax.set_xticks([1,2],['First step (30)','Second step (9)'])
    ax.set_xlim(.55,2.95);ax.set_ylabel('Computed − measured potential (V)');ax.set_title('The error depends on the step')
    heights=[b['metrics']['MAE_V']]+[r['metrics']['MAE_V'] for r in b['molecule_held_out']]
    bar.bar([0,1,2],heights,color=[INK,TEAL,AMBER],width=.65)
    for i,v in enumerate(heights):bar.text(i,v+.01,f'{v:.3f}',ha='center',weight='bold',fontsize=18)
    bar.set_xticks([0,1,2],['Unfitted\n(primary)','One offset\n(molecule-held-out)','Step offsets\n(molecule-held-out)'],fontsize=11)
    bar.set_ylim(0,.31);bar.set_ylabel('MAE (V)');bar.set_title('Offsets diagnose structure in the error')
    return fig


def sensitivity(d):
    rows=d['sensitivity']['basis_rows'];fig,(a,b)=plt.subplots(1,2,figsize=(11.3,3.8),layout='constrained')
    labels=['BQ','Naphthoquinone','Phenanthrenequinone','N-methylphthalimide','PMDI (iPr)']
    nuclear=np.array([r['delta_thermal_eV'] for r in rows])*1000
    geom=np.array([r['delta_geometry_SP_eV'] for r in rows])*1000
    a.barh(range(5),nuclear,color=TEAL,label='Nuclear term');a.barh(range(5),geom,left=nuclear,color=AMBER,label='Geometry / SP')
    a.set_yticks(range(5),labels);a.invert_yaxis();a.set_xlabel('Change in first potential (mV)');a.set_title('Five neutral SVPD controls')
    a.legend(fontsize=10,loc='lower left',bbox_to_anchor=(0,1.02),ncol=2)
    a.set_title('Five neutral SVPD controls',pad=38);a.set_xlim(-45,0)
    values=sorted([float(r['range_mV']) for r in d['sensitivity']['rotor_cutoff_rows']],reverse=True)
    b.plot(range(1,40),values,c=TEAL,lw=2);b.scatter(range(1,40),values,c=TEAL,s=18)
    b.set_xlabel('Benchmark endpoint, ordered by span');b.set_ylabel('Potential span (mV)')
    b.set_title('25 / 50 / 100 cm$^{-1}$ rotor settings');b.set_ylim(0,30);b.grid(alpha=.12)
    b.text(18,24,f'Maximum: {values[0]:.2f} mV',fontsize=16,weight='bold')
    return fig


def pareto(d):
    selection=d['selection'];front=set(selection['descriptor_front']['descriptor_front_ids'])
    # Same domain as the source figure: only reduction-active two-potential models.
    root=Path(__file__).resolve().parents[3]
    import sys
    sys.path.insert(0,str(root/'scripts/plotting/candidates'))
    from candidate_highlights import eligible_anolytes
    rows=selection['descriptor_candidates']
    helper=root/'scripts/plotting/candidates/candidate_highlights.py'
    assert hashlib.sha256(helper.read_bytes()).hexdigest()==selection['descriptor_front']['source_sha256'][str(helper)]
    filtered=eligible_anolytes(rows)
    assert len(filtered)==selection['descriptor_front']['plotted']==56
    # The selected IDs are supplied by the source-frozen three-objective comparison.
    fig,ax=plt.subplots(figsize=(8.4,4.2),layout='constrained')
    sc=ax.scatter([r['E1_vs_Fc_V'] for r in filtered],[r['capacity_nominal_mAh_g'] for r in filtered],
        c=[r['E2_vs_Fc_V'] for r in filtered],cmap='viridis',s=60,alpha=.85)
    fr=[r for r in filtered if r['id'] in front]
    ax.scatter([r['E1_vs_Fc_V'] for r in fr],[r['capacity_nominal_mAh_g'] for r in fr],
               facecolors='none',edgecolors=AMBER,s=175,lw=2,label='Three-objective descriptor front')
    ax.set_xlabel('First reduction potential vs Fc (V)');ax.set_ylabel('Nominal molecular capacity (mAh/g)')
    fig.colorbar(sc,ax=ax,label='Second potential vs Fc (V)');ax.legend(fontsize=10,loc='lower left')
    ax.grid(alpha=.12)
    assert len(fr)==len(front)==8
    return fig


def base_shift(d):
    data=d['base_assisted'];rows=data['candidates'];fig,ax=plt.subplots(figsize=(10.8,3.7),layout='constrained')
    labels={'design_pmdi_ethyl':'PMDI','design_ndi_ethyl':'NDI','exp_bq_oxy_methoxy':'BQ OMe',
            'exp_bq_oxy_trimethyl':'BQ Me$_3$','exp_aq_1oxy':'AQ-1-O'}
    good=[r for r in rows if isinstance(r.get('reference_deltaG_kJ_mol'),(int,float))]
    x=np.arange(len(good))
    ref=np.array([r['reference_deltaG_kJ_mol'] for r in good])
    for label,delta,color,marker in [('HCl reference',0,INK,'o'),('DIPEA (separated ions)',data['corrections']['dipea']['free_ions_kJ_mol'],AMBER,'s'),
                                    ('DBU (separated ions)',data['corrections']['dbu']['free_ions_kJ_mol'],TEAL,'^')]:
        ax.plot(x,ref+delta,marker=marker,ms=8,lw=2,color=color,label=label)
    ax.axhline(0,c='#7B8791',ls=':',lw=1);ax.set_xticks(x,[labels.get(r['id'],r['id']) for r in good])
    ax.set_ylabel('Defined reaction ΔG (kJ/mol)');ax.legend(ncol=3,fontsize=11,loc='lower center',bbox_to_anchor=(.5,1.01))
    ax.grid(axis='y',alpha=.12)
    return fig


PLOTS={'library':library,'exact_variants':structures,'pmdi_model':pmdi_model,'matched_redox':matched,'benchmark':benchmark,
       'bias':bias,'sensitivity':sensitivity,'descriptor_projection':pareto,'base_shift':base_shift}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--only',choices=list(PLOTS));a=p.parse_args();d=json.loads(a.evidence.read_text());style();a.out.mkdir(parents=True,exist_ok=True)
    output={}
    for name,fn in PLOTS.items():
        if a.only and a.only!=name:continue
        fig=fn(d)
        for ext in ['pdf','png']:
            dest=a.out/f'{name}.{ext}';fig.savefig(dest,dpi=190,bbox_inches='tight',facecolor='white')
            output[str(dest)]=hashlib.sha256(dest.read_bytes()).hexdigest()
        plt.close(fig)
    (a.out/('figures.json' if not a.only else a.only+'_figure.json')).write_text(json.dumps(dict(
        status='figures_rendered_visual_review_pending',evidence=str(a.evidence),
        evidence_sha256=hashlib.sha256(a.evidence.read_bytes()).hexdigest(),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),outputs=output),indent=2)+'\n')


if __name__=='__main__':main()
