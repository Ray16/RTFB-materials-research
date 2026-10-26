#!/usr/bin/env python
"""Plot sampled PMDI polymer-context penalties from checked cached records."""
import argparse
import hashlib
import json
from pathlib import Path


def series(record):
    rows={row['site']:row for row in record['singles']}
    if set(rows)!={2,3} or record['conformer_convergence_established'] or record['synthesis_admission']:
        raise ValueError('Expected the explicitly finite, non-admitting two-site comparison')
    cutoffs=[25,50,100]
    penalties={site:[rows[site]['by_cutoff'][str(k)]['polymer_minus_capped_kJ_mol'] for k in cutoffs] for site in (2,3)}
    placement=[rows[2]['by_cutoff'][str(k)]['first_graft_kJ_mol']['free_ions']-
        rows[3]['by_cutoff'][str(k)]['first_graft_kJ_mol']['free_ions'] for k in cutoffs]
    return dict(cutoffs_cm= cutoffs,polymer_minus_local_kJ_mol=penalties,site2_minus_site3_kJ_mol=placement)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':9,'legend.fontsize':8,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    data=series(json.loads(args.source.read_text()))
    figure,axes=plt.subplots(1,2,figsize=(6.5,2.65),layout='constrained')
    for site,color,marker in [(2,'#0072B2','o'),(3,'#D55E00','s')]:
        axes[0].plot(data['cutoffs_cm'],data['polymer_minus_local_kJ_mol'][site],
            marker=marker,color=color,lw=1.2,ms=4,label=f'Site {site}')
    axes[0].set(ylabel=r'$\Delta G_{\mathrm{polymer}}-\Delta G_{\mathrm{local}}$ (kJ mol$^{-1}$)',ylim=(-.2,7.2),title='Polymer-context penalty')
    axes[0].legend(frameon=False,loc='lower right')
    axes[1].plot(data['cutoffs_cm'],data['site2_minus_site3_kJ_mol'],color='#665191',marker='D',ms=4,lw=1.2)
    axes[1].set(ylabel=r'$\Delta G_{\mathrm{site\ 2}}-\Delta G_{\mathrm{site\ 3}}$ (kJ mol$^{-1}$)',ylim=(-1,.4),title='Site preference changes sign')
    for label,axis in zip('ab',axes):
        axis.axhline(0,color='.45',ls='--',lw=.7,zorder=0)
        axis.set(xlabel=r'Rotor cutoff (cm$^{-1}$)',xticks=data['cutoffs_cm'])
        axis.text(-.16,1.08,label,transform=axis.transAxes,fontweight='bold',fontsize=11)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    paths=[args.out.with_suffix('.'+extension) for extension in ('pdf','png','svg')]
    for path in paths:figure.savefig(path,dpi=300)
    plt.close(figure)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    proof=dict(status='cached_first_site_figure_generated',source=str(args.source),source_sha256=sha(args.source),
        producer=str(Path(__file__).resolve()),producer_sha256=sha(__file__),data=data,
        outputs_sha256={str(p):sha(p) for p in paths},scope='One sampled minimum per occupied graph; rotor sensitivity, no ensemble or statistical interval',
        base_terms_cancel_from_penalty=True,isolated_motif_interaction_inferred=False)
    args.out.with_suffix('.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(dict(status=proof['status'],outputs=[str(p) for p in paths])))


if __name__=='__main__':
    main()
