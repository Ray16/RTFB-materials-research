#!/usr/bin/env python3
"""Plot replayed input-state coverage; counts do not imply ensemble convergence."""
import argparse
import json
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def plot(census_path,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    source_hash=sha256_file(census_path);data=json.loads(census_path.read_text())
    if data['status']!='distinct_source_replayed_native_state_census' or data['global_ensemble_convergence_proven']:
        raise ValueError('An explicitly scoped native-state census is required')
    configurations=['all_same','alternating'];sites=[[],[2],[3],[4],[2,3],[2,4]]
    rows={(r['configuration'],tuple(r['occupied_sites'])):r for r in data['groups']}
    if len(rows)!=12 or any(r['expected']!=20 for r in rows.values()):raise ValueError('Expected two configurations, six graphs and twenty declared starts')
    matrix=np.array([[rows[c,tuple(s)]['checked'] for s in sites] for c in configurations])
    if int(matrix.sum())!=data['checked_states'] or np.any(matrix<0) or np.any(matrix>20):raise ValueError('Census counts differ')
    if any(output.with_suffix('.'+ext).exists() for ext in ['pdf','png','json']):raise FileExistsError('Preserve prior coverage plots')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(10.4,3.8));fig.subplots_adjust(left=.16,right=.97,top=.77,bottom=.23)
    ax.imshow(matrix,cmap=LinearSegmentedColormap.from_list('coverage',['#f0f6f5','#007f86']),vmin=0,vmax=20,aspect='auto')
    ax.set_xticks(range(6),['Empty','Site 2','Site 3','Site 4','Sites 2+3','Sites 2+4']);ax.set_yticks(range(2),['All same','Alternating'])
    ax.tick_params(length=0,pad=10)
    for spine in ax.spines.values():spine.set_visible(False)
    for i in range(2):
        for j in range(6):ax.text(j,i,f'{matrix[i,j]}/20',ha='center',va='center',fontsize=17,color='white' if matrix[i,j]>=12 else '#15343b')
    fig.text(.16,.94,'Native state coverage of the phthalimide spacing control',fontsize=16,weight='bold',color='#15343b')
    fig.text(.16,.85,f"{data['checked_states']} independently replayed input states · two declared backbone configurations",fontsize=11,color='#4b6266')
    fig.text(.16,.055,'Counts include duplicate minima as separate input states. Clustering and complete balanced\nensemble comparisons are still required; no experimental occupancy distribution is inferred.',fontsize=10,color='#4b6266')
    output.parent.mkdir(parents=True,exist_ok=True);artifacts={}
    for ext in ['pdf','png']:
        target=output.with_suffix('.'+ext);fig.savefig(target,dpi=190,facecolor='white');artifacts[str(target)]=sha256_file(target)
    plt.close(fig)
    if sha256_file(census_path)!=source_hash:raise ValueError('Census changed during rendering')
    atomic_write_json(output.with_suffix('.json'),dict(status='replayed_input_state_coverage_figure',census=str(census_path),census_sha256=source_hash,
        producer_sha256=sha256_file(__file__),checked_states=data['checked_states'],outputs_sha256=artifacts,ensemble_convergence_claimed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--census',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    plot(a.census.resolve(),a.out.resolve())
