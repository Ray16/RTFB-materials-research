#!/usr/bin/env python3
"""Plot two checked exact-candidate potentials without inventing confidence bars."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(source,output):
    result=json.loads(source.read_text())
    if (result['status']!='exact_molecular_native_redox_published_with_field_QC'
            or len(result['redox']['events'])!=2 or result['material_admission']):
        raise ValueError('This plot requires the checked three-state molecular publication')
    events=result['redox']['events'];y=np.array([r['potential_V_vs_Fc_by_cutoff']['50'] for r in events])
    all_y=np.array([list(r['potential_V_vs_Fc_by_cutoff'].values()) for r in events])
    colors=['#007F86','#BD6928']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(5.6,3.7),layout='constrained')
    for i,value in enumerate(y):
        ax.errorbar(i,value,yerr=[[value-all_y[i].min()],[all_y[i].max()-value]],fmt='o',
            markersize=8,color=colors[i],capsize=5,elinewidth=1.4)
        ax.annotate(f'{value:+.3f} V',(i,value),xytext=(0,12),textcoords='offset points',ha='center',color=colors[i],fontweight='bold')
    charge=lambda q:f'{q:+d}' if q else '0'
    labels=[charge(e['initial_charge'])+' → '+charge(e['final_charge'])+'\n'
        +f"multiplicity {e['initial_multiplicity']} → {e['final_multiplicity']}" for e in events]
    ax.set_xticks([0,1],labels)
    ax.set_ylabel('Computed reduction potential vs Fc/Fc⁺ (V)')
    ax.set_xlim(-.45,1.45);ax.set_ylim(float(all_y.min())-.16,float(all_y.max())+.18)
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    report=dict(source=str(source.resolve()),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        figure_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),primary_potentials_V=y.tolist(),
        rotor_sensitivity_V=[[float(row.min()),float(row.max())] for row in all_y],
        statistical_intervals_plotted=False,experimental_exact_graft_accuracy_claimed=False,
        scope='Declared-spin capped graft; bars span 25–100 cm-1 rotor sensitivity, not a prediction interval')
    output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(plot(a.source,a.out)))
