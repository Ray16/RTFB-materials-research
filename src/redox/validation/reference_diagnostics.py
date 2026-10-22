"""Post hoc reference diagnostics; never replace the primary unfitted benchmark."""
from collections import defaultdict
import math


def metrics(values):
    values=list(values)
    if not values or not all(math.isfinite(v) for v in values):
        return dict(n=len(values),MAE_V=None,RMSE_V=None,bias_V=None)
    return dict(n=len(values),MAE_V=math.fsum(abs(v) for v in values)/len(values),
        RMSE_V=math.sqrt(math.fsum(v*v for v in values)/len(values)),bias_V=math.fsum(values)/len(values))


def cross_validated_offsets(rows,group,step_specific=False):
    """Fit only on training groups; leave absent step models explicitly missing."""
    if group not in ('endpoint','id','doi'):raise ValueError('Unknown holdout group')
    totals=defaultdict(lambda:[0.,0]);excluded=defaultdict(lambda:[0.,0])
    for index,row in enumerate(rows):
        step=row['wave'] if step_specific else 'all'
        unit=index if group=='endpoint' else row[group]
        totals[step][0]+=row['error_V'];totals[step][1]+=1
        excluded[unit,step][0]+=row['error_V'];excluded[unit,step][1]+=1
    predictions=[]
    for index,row in enumerate(rows):
        step=row['wave'] if step_specific else 'all'
        unit=index if group=='endpoint' else row[group]
        count=totals[step][1]-excluded[unit,step][1]
        offset=(totals[step][0]-excluded[unit,step][0])/count if count else None
        predictions.append(dict(id=row['id'],event=row['event'],doi=row['doi'],wave=row['wave'],
            held_group=unit,training_count=count,fitted_offset_V=offset,
            corrected_error_V=row['error_V']-offset if offset is not None else None))
    available=[r['corrected_error_V'] for r in predictions if r['corrected_error_V'] is not None]
    return dict(holdout=group,model='separate_step_offsets' if step_specific else 'constant_offset',
        coverage=len(available),target=len(rows),metrics=metrics(available),predictions=predictions,
        interpretation='Post hoc internal cross-validation; not an untouched external test')


def independent_cluster_intervals(rows,field,replicates=2000,seed=20261021):
    """Use cluster sufficient statistics, independent of the production resampler."""
    import numpy as np
    groups=defaultdict(list)
    for row in rows:groups[row[field]].append(row['error_V'])
    if len(groups)<5:return None
    values=[groups[key] for key in sorted(groups)]
    sufficient=np.array([[len(v),math.fsum(v),math.fsum(abs(x) for x in v),math.fsum(x*x for x in v)] for v in values])
    indices=np.random.default_rng(seed).integers(0,len(values),size=(replicates,len(values)))
    sums=sufficient[indices].sum(axis=1)
    estimates=dict(MAE_V=sums[:,2]/sums[:,0],RMSE_V=np.sqrt(sums[:,3]/sums[:,0]),bias_V=sums[:,1]/sums[:,0])
    return {key:np.percentile(value,[2.5,97.5]).tolist() for key,value in estimates.items()}


def organic_reference(rows,identifier):
    anchors=[r for r in rows if r['id']==identifier and r['wave']==1]
    if len(anchors)!=1:raise ValueError('An organic reference needs one exact first-step endpoint')
    anchor=anchors[0]
    targets=[r for r in rows if r['id']!=identifier]
    residuals=[]
    for row in targets:
        relative_calc=row['calc_V']-anchor['calc_V'];relative_exp=row['exp_V']-anchor['exp_V']
        residuals.append(relative_calc-relative_exp)
    return dict(anchor=anchor,targets_exclude_entire_anchor_molecule=True,
        uncorrected_same_targets=metrics(r['error_V'] for r in targets),
        Fc_independent_relative_metrics=metrics(residuals),
        interpretation='Reference cancellation diagnostic; the organic anchor has its own model and experimental errors')
