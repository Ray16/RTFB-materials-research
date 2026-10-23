#!/usr/bin/env python
"""Show the molecular descriptor front alongside full-axis eligibility status."""
import argparse,hashlib,importlib.util,json,os
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.screening.feasibility import feasible_ids
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def nondominated(rows):
    values={r['id']:(-r['E1_vs_Fc_V'],-r['E2_vs_Fc_V'],r['capacity_nominal_mAh_g']) for r in rows}
    if len(values)!=len(rows):raise ValueError('Duplicate candidate identity in descriptor comparison')
    return {key for key,v in values.items() if not any(other!=key and all(a>=b for a,b in zip(w,v))
        and any(a>b for a,b in zip(w,v)) for other,w in values.items())}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--decisions',type=Path,required=True)
    parser.add_argument('--out',type=Path,default=ROOT/'results/figures/selection');args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    folder=args.decisions/'source_snapshot';feed_path=folder/'feed.json';policy_path=folder/'policy.json';decision_path=args.decisions/'decisions.json'
    feed=json.loads(feed_path.read_text());policy=json.loads(policy_path.read_text());decisions=json.loads(decision_path.read_text())['candidates']
    eligible=feasible_ids(decisions,policy=policy);by_key={r['product_key']:r for r in decisions}
    helper=ROOT/'scripts/plotting/candidates/candidate_highlights.py'
    spec=importlib.util.spec_from_file_location('descriptor_plot_eligibility',helper);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=module.eligible_anolytes(feed['original39']+feed['additional_candidates']);front=nondominated(rows)
    if any(r['exact_grafted_inchikey'] not in by_key for r in rows):raise ValueError('A plotted graph is absent from the decision ledger')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    sources={str(p):sha(p) for p in [Path(__file__),helper,feed_path,policy_path,decision_path]}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4),sharey=True);fig.subplots_adjust(left=.1,right=.98,bottom=.25,top=.86,wspace=.15)
    for i,(axis,field,title) in enumerate(zip(axes,['E1_vs_Fc_V','E2_vs_Fc_V'],['First reduction','Second reduction'])):
        for row in rows:
            decision=by_key[row['exact_grafted_inchikey']];qualified=decision['id'] in eligible;is_front=row['id'] in front
            axis.scatter(row[field],row['capacity_nominal_mAh_g'],s=38 if is_front else 20,
                facecolors='#009E73' if qualified else 'white',edgecolors='#7B3294' if is_front else '#AEB5BD',linewidths=1.2 if is_front else .6,zorder=3 if is_front else 2)
        axis.grid(axis='y',alpha=.15);axis.set(xlabel='Calculated potential (V vs Fc/Fc⁺)',title=title)
        axis.text(-.12,1.08,chr(97+i),transform=axis.transAxes,fontweight='bold')
    axes[0].set_ylabel('Nominal molecular capacity (mAh g⁻¹)')
    handles=[Line2D([],[],marker='o',markerfacecolor='white',markeredgecolor='#AEB5BD',linestyle='',label='Other molecular models'),
        Line2D([],[],marker='o',markerfacecolor='white',markeredgecolor='#7B3294',linestyle='',label='Descriptor Pareto set')]
    qualified_plotted=[r['id'] for r in rows if by_key[r['exact_grafted_inchikey']]['id'] in eligible]
    if qualified_plotted:handles.append(Line2D([],[],marker='o',color='#009E73',linestyle='',label='All-axis eligible in declared model'))
    fig.legend(handles=handles,
        loc='upper center',bbox_to_anchor=(.54,.985),ncol=2,frameon=False,fontsize=7)
    fig.text(.1,.09,'Open markers: at least one required feasibility axis is unresolved.',fontsize=7)
    fig.text(.1,.035,f'Three descriptor objectives: E₁, E₂ and nominal capacity. All-axis eligible among plotted models: {len(qualified_plotted)}.',fontsize=7)
    args.out.mkdir(parents=True,exist_ok=True);stem=args.out/'descriptor_pareto';outputs=[]
    for ext in ('pdf','png','svg'):
        p=stem.with_suffix('.'+ext);fig.savefig(p,dpi=600);outputs.append(p)
    plt.close(fig)
    if any(sha(p)!=h for p,h in sources.items()):raise ValueError('Plot inputs changed')
    atomic_write_json(stem.with_suffix('.json'),dict(status='partial_molecular_descriptor_front_with_full_axis_overlay',
        input_candidate_as_of=feed.get('as_of_Chicago'),plotted=len(rows),descriptor_front_ids=sorted(front),
        all_axes_eligible_plotted_ids=qualified_plotted,polymer_performance_predicted=False,
        axes=dict(E1_vs_Fc_V='minimize',E2_vs_Fc_V='minimize',capacity_nominal_mAh_g='maximize'),
        no_candidate_uncertainty_intervals_imputed=True,source_sha256=sources,output_sha256={str(p):sha(p) for p in outputs}))
    print(json.dumps(dict(plotted=len(rows),descriptor_front=len(front),all_axes_eligible=len(eligible))))


if __name__=='__main__':main()
