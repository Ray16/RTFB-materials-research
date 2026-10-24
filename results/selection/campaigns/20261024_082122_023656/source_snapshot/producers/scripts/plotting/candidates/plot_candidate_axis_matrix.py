#!/usr/bin/env python
"""Show exact-case decision status without turning supporting proxies into passes."""
import argparse,hashlib,json,os
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.screening.feasibility import feasible_ids
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

SELECTION=[('exp_bq_oxy_methoxy','Methoxy BQ'),('exp_bq_oxy_trimethyl','Trimethyl BQ'),('exp_aq_1oxy','1-O AQ'),
 ('design_pmdi_ethyl','PMDI ethyl'),('design_ndi_ethyl','NDI ethyl'),('exp_ptz_37dimethoxy','Dimethoxy phenothiazine'),
 ('hyp_ptz_37difluoro','Difluoro phenothiazine'),('exp_viologen_methoxyethyl','Ether-cap viologen')]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--decisions',type=Path,required=True)
    parser.add_argument('--out',type=Path,default=ROOT/'results/figures/selection');args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    policy_path=args.decisions/'source_snapshot/policy.json';data_path=args.decisions/'decisions.json'
    policy=json.loads(policy_path.read_text());rows=json.loads(data_path.read_text())['candidates'];feasible_ids(rows,policy=policy)
    by_id={r['id']:r for r in rows};axes=policy['required_axes'];selected=[(by_id[key],name) for key,name in SELECTION if key in by_id]
    if len(selected)!=len(SELECTION):raise ValueError('A declared displayed case is absent')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    sources={str(p):sha(p) for p in (Path(__file__),policy_path,data_path)}
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    code={'pending':0,'pass':1,'fail':2};matrix=np.array([[code[next(a['status'] for a in r['axes'] if a['axis']==rule['id'])] for rule in axes] for r,_ in selected])
    fig,ax=plt.subplots(figsize=(8.2,4.4));fig.subplots_adjust(left=.24,right=.98,bottom=.30,top=.9)
    ax.imshow(matrix,cmap=ListedColormap(['#E5E7EA','#009E73','#D55E00']),vmin=0,vmax=2,aspect='auto',interpolation='nearest')
    labels=['Redox\nwindow','Material\ncapacity','State\nstability','Intended\nET','Graft\nroute','Graft\nΔG','Motif\ninteraction','Intersite\ncoupling','Loading/\nspacing','Access\nmodel']
    ax.set(xticks=range(len(axes)),xticklabels=labels,yticks=range(len(selected)),yticklabels=[name for _,name in selected])
    ax.tick_params(length=0);ax.set_xticks(np.arange(-.5,len(axes),1),minor=True);ax.set_yticks(np.arange(-.5,len(selected),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=2);ax.tick_params(which='minor',length=0)
    support=[]
    for i,(row,name) in enumerate(selected):
        for j,rule in enumerate(axes):
            record=row['evidence'][rule['id']]
            if record.get('supporting_data_available'):
                ax.scatter(j,i,s=16,color='#0072B2');support.append(dict(id=row['id'],axis=rule['id'],meaning='Supporting model/proxy/review only; not an axis pass'))
    for spine in ax.spines.values():spine.set_visible(False)
    fig.suptitle('All required axes govern candidate admission',fontsize=11,y=.99)
    fig.legend(handles=[Line2D([],[],marker='s',color='#E5E7EA',linestyle='',label='Pending'),Line2D([],[],marker='s',color='#009E73',linestyle='',label='Pass'),
        Line2D([],[],marker='s',color='#D55E00',linestyle='',label='Fails declared criterion'),Line2D([],[],marker='o',color='#0072B2',linestyle='',label='Supporting model/proxy/review')],
        loc='lower center',bbox_to_anchor=(.56,.09),ncol=2,frameon=False,fontsize=7)
    fig.text(.24,.025,'Pending includes missing evidence or an undeclared quantitative target; it does not mean poor chemistry.',fontsize=7)
    args.out.mkdir(parents=True,exist_ok=True);stem=args.out/'candidate_axis_matrix';outputs=[]
    for ext in ('pdf','png','svg'):
        p=stem.with_suffix('.'+ext);fig.savefig(p,dpi=600);outputs.append(p)
    plt.close(fig)
    if any(sha(p)!=h for p,h in sources.items()):raise ValueError('Matrix input changed')
    atomic_write_json(stem.with_suffix('.json'),dict(status='all_axes_readiness_matrix',displayed_ids=[r['id'] for r,_ in selected],
        raw_statuses=matrix.tolist(),code=code,supporting_information_markers=support,source_sha256=sources,output_sha256={str(p):sha(p) for p in outputs}))


if __name__=='__main__':main()
