"""Descriptive tradeoff highlights; never replace the qualified primary Pareto policy."""
import json
import math


def finite(value):
    return type(value) in (int,float) and math.isfinite(value)


def eligible_anolytes(rows):
    eligible=[]
    for row in rows:
        states=json.loads(row['states']) if isinstance(row['states'],str) else row['states']
        if any(state[3]>0 for state in states):continue
        if row.get('model_primary_admissible') is not True:continue
        if row.get('solution_potentials_status')!='verified native calculation':continue
        if row.get('potential_display_status','Admitted calculation')!='Admitted calculation':continue
        if row.get('n_e_nominal')!=2:continue
        if not all(finite(row.get(k)) for k in ('E1_vs_Fc_V','E2_vs_Fc_V','capacity_nominal_mAh_g','accessible_capacity_mAh_g')):continue
        # Published nominal capacities are rounded to 0.1 mAh/g; a one-electron
        # integrity/window path cannot be promoted to a two-electron comparison.
        if row['accessible_capacity_mAh_g']<row['capacity_nominal_mAh_g']-.1:continue
        eligible.append(row)
    return eligible


def highlights(rows):
    from redox.screening.pareto import _dominates
    eligible=eligible_anolytes(rows)
    values={r['id']:dict(first=(-r['E1_vs_Fc_V'],0.),second=(-r['E2_vs_Fc_V'],0.),
                        capacity=(r['capacity_nominal_mAh_g'],0.)) for r in eligible}
    front=[r for r in eligible if not any(_dominates(values[o['id']],values[r['id']],k=0)
                                          for o in eligible if o['id']!=r['id'])]
    selected=[]
    for row in sorted(front,key=lambda r:(-r['capacity_nominal_mAh_g'],r['id']))[:2]:selected.append(row['id'])
    if front:
        best=min(front,key=lambda r:(r['E1_vs_Fc_V'],r['id']))['id']
        if best not in selected:selected.append(best)
    return dict(status='partial_information_descriptor_comparison',eligible_ids=[r['id'] for r in eligible],
        partial_front_ids=[r['id'] for r in front],highlighted_ids=selected,
        axes=dict(E1_vs_Fc_V='minimize',E2_vs_Fc_V='minimize',capacity_nominal_mAh_g='maximize'),
        selection_rule='Highlight the two highest nominal capacities on the displayed front and its most negative E1.',
        SA_is_annotation_not_objective=True,
        qualified_solution_lambda_count=sum(finite(r.get('lambda_i_eV')) for r in eligible),
        qualified_primary_Pareto_front_published=False,
        uncertainty_treatment='Exact descriptor dominance only; no validated candidate-specific error bounds available.',
        polymer_synthesis_ranking=False,
        molecular_scope='p-methylbenzyl-capped local graft models; vinyl-bearing monomers and polymers require separate qualification')


def matched_highlights(payload, admitted):
    by_id={r['id']:r for r in admitted};groups={}
    for comparison in payload.get('matched_comparisons',[]):
        row=by_id.get(comparison['id']);control=by_id.get(comparison['baseline'])
        if row is None or control is None:continue
        if len(eligible_anolytes([row,control]))!=2:continue
        first=row['E1_vs_Fc_V']-control['E1_vs_Fc_V'];second=row['E2_vs_Fc_V']-control['E2_vs_Fc_V']
        if not first<0 or not second<0:continue
        core=row.get('search_core') or ('ndi' if 'ndi' in row['id'] else 'pmdi' if 'pmdi' in row['id'] else None)
        if core is None:continue
        record=dict(id=row['id'],baseline=control['id'],delta_E1_V=first,delta_E2_V=second,
                    delta_nominal_capacity_mAh_g=row['capacity_nominal_mAh_g']-control['capacity_nominal_mAh_g'],core=core)
        if core not in groups or min(-first,-second)>min(-groups[core]['delta_E1_V'],-groups[core]['delta_E2_V']):groups[core]=record
    return [groups[core] for core in sorted(groups)]


def focused_figures(rows, report, matches, out, names, colors):
    import matplotlib.pyplot as plt
    by_id={r['id']:r for r in rows}
    fig,ax=plt.subplots(figsize=(7.2,4.5))
    for row in rows:
        if not all(finite(row.get(k)) for k in ('SA_grafted','capacity_nominal_mAh_g')):continue
        ax.scatter(row['SA_grafted'],row['capacity_nominal_mAh_g'],s=20,color='#bdc5cc',alpha=.6)
    for gid in report['partial_front_ids']:
        row=by_id[gid]
        ax.scatter(row['SA_grafted'],row['capacity_nominal_mAh_g'],s=40,color=colors[row['family']],edgecolor='white',linewidth=.5)
    for gid in report['highlighted_ids']:
        row=by_id[gid]
        ax.scatter(row['SA_grafted'],row['capacity_nominal_mAh_g'],s=100,marker='*',color='#111',zorder=5)
        ax.annotate(names[gid],(row['SA_grafted'],row['capacity_nominal_mAh_g']),xytext=(7,7),textcoords='offset points',fontsize=7)
    ax.set(xlabel='Synthetic accessibility score (lower is easier; heuristic)',ylabel='Nominal capacity (mAh g⁻¹)')
    ax.grid(alpha=.15)
    fig.subplots_adjust(left=.12,right=.97,top=.94,bottom=.2)
    fig.text(.12,.06,'Colored: potential–capacity partial front. Stars: selected descriptor leaders.\nGray: other or pending designs. SA does not establish a polymer synthesis route.',fontsize=7)
    for suffix in ('pdf','svg','png'):fig.savefig(out/f'synthesis_tradeoff.{suffix}',dpi=600)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,3.8))
    labels=[]
    for i,match in enumerate(matches):
        ax.scatter(match['delta_E1_V'],i+.1,color='#0072B2',marker='o',s=40,label='First event' if i==0 else None)
        ax.scatter(match['delta_E2_V'],i-.1,color='#D55E00',marker='s',s=40,label='Second event' if i==0 else None)
        label=names.get(match['id'],match['id'])+'\n− '+names.get(match['baseline'],match['baseline'])
        labels.append(label)
    ax.set_yticks(range(len(labels)),labels,fontsize=8)
    ax.axvline(0,color='#777',linewidth=.8);ax.grid(axis='x',alpha=.15)
    ax.set_xlabel('Potential shift from matched control (V)')
    if matches:ax.legend(frameon=False,fontsize=8,loc='best')
    else:ax.text(.5,.5,'No eligible matched shifts',ha='center',transform=ax.transAxes)
    fig.subplots_adjust(left=.42,right=.97,top=.95,bottom=.27)
    fig.text(.08,.05,'One contrast per core with the strongest decrease shared by both events.\nComputed local graft models; a potential shift does not prove polymer performance.',fontsize=7)
    for suffix in ('pdf','svg','png'):fig.savefig(out/f'new_candidate_shifts.{suffix}',dpi=600)
    plt.close(fig)
