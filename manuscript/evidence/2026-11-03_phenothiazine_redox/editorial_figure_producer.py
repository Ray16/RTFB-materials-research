"""Publication layout of the frozen four-frame electronic comparison."""
import json
import os
from pathlib import Path

def main():
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    source=ROOT/'manuscript/evidence/2026-11-02_occupied_sampling_checked/motif_definition_comparison.json'
    data=json.loads(source.read_text());rows=data['states']
    assert len(rows)==data['compared_central_cap_frames']==4
    keys=['all_same__occ_2_3','all_same__occ_2_4','alternating__occ_2_3','alternating__occ_2_4']
    rows=[next(r for r in rows if r['parent_key']==key) for key in keys]
    values=[[r['whole_chain_nonadditivity_E_kJ_mol'] for r in rows],
            [r['frozen_methyl_motif_interaction_E_kJ_mol'] for r in rows]]
    assert np.isfinite(values).all()
    target=ROOT/'results/figures/editorial/2026-11-03';target.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
                         'pdf.fonttype':42,'ps.fonttype':42})
    outputs={}
    for gray in (False,True):
        fig,ax=plt.subplots(figsize=(7.1,3.45),layout='constrained');y=np.arange(4)
        colors=['#333333','#bbbbbb'] if gray else ['#007F86','#BD6928']
        for offset,series,color,hatch,label in zip((-.18,.18),values,colors,('', '///'),
                ('Whole-chain occupancy cycle','Frozen methyl-capped motif pair')):
            ax.barh(y+offset,series,height=.32,color=color,edgecolor='#25343A',linewidth=.45,hatch=hatch,label=label)
            for pos,value in zip(y+offset,series):
                ax.text(value+(.35 if value>=0 else -.35),pos,f'{value:+.2f}',
                        ha='left' if value>=0 else 'right',va='center',fontsize=9)
        ax.axvline(0,color='#333333',linewidth=.9)
        ax.set(yticks=y,yticklabels=['All-same · adjacent','All-same · one spacer',
            'Alternating · adjacent','Alternating · one spacer'],xlim=(-17.5,23.5),
            xlabel=r'Electronic response / kJ mol$^{-1}$')
        ax.invert_yaxis();ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
        ax.legend(loc='lower center',bbox_to_anchor=(.48,1.02),ncols=2,frameon=False,fontsize=9)
        stem='viologen_interaction_definitions_final'+('_grayscale' if gray else '')
        for ext in ('.pdf','.png'):
            p=target/(stem+ext);fig.savefig(p,dpi=300);outputs[str(p.relative_to(ROOT))]=sha256_file(p)
        plt.close(fig)
    proof=dict(status='editorial_layout_only_quantities_unchanged',source=str(source),source_sha256=sha256_file(source),
        keys=keys,whole_chain_E_kJ_mol=values[0],motif_E_kJ_mol=values[1],cap_length_A=1.09,
        error_bars_or_population_weights_added=False,grayscale_distinguishable_by_hatching=True,
        output_sha256=outputs,producer_sha256=sha256_file(__file__))
    (target/'figure_review.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(dict(status=proof['status'],outputs=len(outputs))))

if __name__=='__main__':main()
