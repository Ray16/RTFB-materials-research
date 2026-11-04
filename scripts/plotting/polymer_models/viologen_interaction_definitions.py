#!/usr/bin/env python3
"""Show two scoped E definitions at matched ionic parent frames."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from redox.core.io import atomic_write_json
from redox.core.hashing import sha256_file

SOURCE_BYTES=Path(__file__).read_bytes()


def plot(source,output):
    data=json.loads(source.read_text());rows=data['states']
    if data['compared_central_cap_frames']!=4 or any(r['frozen_methyl_motif_interaction_E_kJ_mol'] is None for r in rows):
        raise ValueError('All four checked definition comparisons required')
    output.parent.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.labelsize':10,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42})
    fig,ax=plt.subplots(figsize=(7.1,3.3),layout='constrained');x=np.arange(4);width=.34
    whole=[r['whole_chain_nonadditivity_E_kJ_mol'] for r in rows];motif=[r['frozen_methyl_motif_interaction_E_kJ_mol'] for r in rows]
    for offset,values,color,label in [(-width/2,whole,'#0072B2','Whole-chain occupancy cycle'),(width/2,motif,'#D55E00','Frozen methyl-capped motif pair')]:
        bars=ax.bar(x+offset,values,width,color=color,label=label)
        for bar,value in zip(bars,values):ax.annotate(f'{value:.1f}',(bar.get_x()+bar.get_width()/2,value),
            xytext=(0,3 if value>=0 else -4),textcoords='offset points',ha='center',va='bottom' if value>=0 else 'top',fontsize=8)
    ax.axhline(0,color='#444444',linewidth=.7);ax.set_ylim(-18,24)
    ax.set_xticks(x,['All-same\nadjacent','All-same\none spacer','Alternating\nadjacent','Alternating\none spacer'])
    ax.set_ylabel('Electronic response (kJ mol$^{-1}$)');ax.spines[['top','right']].set_visible(False)
    ax.legend(loc='upper left',ncols=2,frameon=False)
    files=[]
    for suffix in ('.pdf','.png'):
        path=output.with_suffix(suffix);fig.savefig(path,dpi=300);files.append(path)
    plt.close(fig)
    from redox.validation.execution_review import freeze_source
    frozen=freeze_source(output.parent,Path(__file__),SOURCE_BYTES)
    result=dict(status='four_frame_interaction_definition_figure_generated',source=str(source),source_sha256=sha256_file(source),
        figure_sha256={str(p):sha256_file(p) for p in files},producer=str(frozen),producer_sha256=hashlib.sha256(SOURCE_BYTES).hexdigest(),
        quantities='electronic energies only; whole-chain cycle and frozen capped-pair interaction are distinct',
        cap_length_A=1.09,cap_variation_is_not_a_confidence_interval=True)
    atomic_write_json(output.with_suffix('.json'),result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(plot(a.source.resolve(),a.out.resolve())))
