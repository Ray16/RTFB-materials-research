#!/usr/bin/env python3
"""Plot both checked electronic solution-relaxation directions for one exact graft."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from redox.core.hashing import sha256_file


def plot(source,output):
    data=json.loads(source.read_text())
    if data['status']!='exact_candidate_solution_descriptor_published_with_field_QC':
        raise ValueError('Checked candidate descriptor required')
    rows=data['components'];geometry=np.array([r['equilibrium_geometry_component_eV'] for r in rows])
    slow=np.array([r['frozen_slow_increment_eV'] for r in rows]);totals=geometry+slow
    if not np.isfinite(totals).all() or np.min(geometry)<0 or np.min(slow)<0:
        raise ValueError('This positive-component plot needs an explicitly reviewed alternate layout')
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError('Preserve earlier checked figures')
    with plt.rc_context({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.spines.left':False,'pdf.fonttype':42}):
        fig,ax=plt.subplots(figsize=(6.5,2.75))
        y=np.array([1,0]);ax.barh(y,geometry,color='#007F86',height=.45,label='Equilibrium geometry component')
        ax.barh(y,slow,left=geometry,color='#D48A48',height=.45,label='Frozen-slow increment')
        for position,g,s,total in zip(y,geometry,slow,totals):
            ax.text(g/2,position,f'{g:.3f}',ha='center',va='center',color='white',weight='bold')
            ax.text(g+s/2,position,f'{s:.3f}',ha='center',va='center',color='#15343B',weight='bold')
            ax.text(total+.035,position,f'{total:.3f}',ha='left',va='center',weight='bold',color='#15343B')
        ax.set_yticks(y,['Reduction\n(+1 → 0)','Oxidation\n(0 → +1)'])
        ax.set_xlim(0,max(totals)*1.17);ax.set_ylim(-.5,1.65)
        ax.set_xlabel('Electronic solution-relaxation energy (eV)');ax.tick_params(axis='y',length=0,pad=10)
        ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.03,1.18),fontsize=10,ncol=1)
        fig.subplots_adjust(left=.22,right=.98,bottom=.26,top=.79)
        fig.savefig(output,bbox_inches='tight');plt.close(fig)
    report=dict(status='checked_directional_solution_pair_figure',source_sha256=sha256_file(source),figure_sha256=sha256_file(output),
                labels=['reduction','oxidation'],totals_eV=totals.tolist(),statistical_intervals_plotted=False,
                symmetric_scalar_qualified=data['qualified_symmetric_single_body_eV'] is not None,
                quantum_calculations_run=False)
    output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();print(json.dumps(plot(a.source,a.out)))
