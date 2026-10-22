#!/usr/bin/env python
"""Diagnose fixed-panel bias without changing its primary uncorrected potentials."""
import argparse
import csv
import json
import math
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.validation.benchmark_panel import review_report
from redox.validation.reference_diagnostics import (metrics,cross_validated_offsets,
    independent_cluster_intervals,organic_reference)
from redox.workflows.chunked import file_hash
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review',type=Path,required=True)
    parser.add_argument('--panel',type=Path,default=ROOT/'results/validation/publication_benchmark/panel_v1.json')
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    panel=json.loads(args.panel.read_text());report=review_report(panel,args.review)
    if report['admitted_events']!=report['target_events']:raise ValueError('Use the complete fixed panel for this analysis')
    targets={(r['id'],r['event'],r['doi']):r for r in panel['targets']};rows=[]
    with (args.review/'comparisons.csv').open() as stream:
        for row in csv.DictReader(stream):
            key=row['id'],row['event'],row['doi']
            if row['kind']!='E' or key not in targets:continue
            if row['status']!='checked_native_parent':raise ValueError('An endpoint is not independently checked')
            target=targets[key];calc=float(row['calc_unrounded_V']);exp=float(row['exp_V'])
            rows.append(dict(id=row['id'],event=row['event'],doi=row['doi'],family=target['family'],wave=target['wave'],
                             calc_V=calc,exp_V=exp,error_V=calc-exp))
    primary=metrics(r['error_V'] for r in rows)
    if any(abs(primary[k]-report['metrics'][k])>1e-12 for k in ('MAE_V','RMSE_V','bias_V')):
        raise ValueError('Primary statistics changed')
    intervals=[]
    for old in report['uncertainty']:
        checked=independent_cluster_intervals(rows,old['cluster'],old['replicates'],old['seed'])
        if any(abs(a-b)>1e-12 for key in checked for a,b in zip(checked[key],old['intervals'][key])):
            raise ValueError('Independent cluster resampling differs')
        intervals.append(dict(old,independently_reconstructed=True,
            caution='Seven unequal study clusters; percentile intervals are weakly calibrated' if old['cluster']=='doi'
                    else 'Retrospective molecule-cluster resampling; studies and families remain shared'))
    by_step={str(step):metrics(r['error_V'] for r in rows if r['wave']==step) for step in (1,2)}
    cv=[cross_validated_offsets(rows,group,step) for group in ('endpoint','id','doi') for step in (False,True)]
    config=json.loads((ROOT/'config/project.json').read_text())['reference'];fc=config['fc_abs_computed_V']
    orbital=-8.617333262145e-5*298.15*math.log(2)
    anchors=['tempo_parent','bm_n_ethylphenothiazine','bm_n_methylphthalimide','bm_1_4_benzoquinone']
    sources=[Path(__file__),ROOT/'src/redox/validation/reference_diagnostics.py',args.panel,
        ROOT/'config/project.json',*args.review.glob('*.json'),args.review/'comparisons.csv']
    result=dict(status='post_hoc_reference_and_bias_diagnostics_checked',primary_unfitted=primary,
        by_step=by_step,first_step_all_errors_positive=all(r['error_V']>0 for r in rows if r['wave']==1),
        mean_squared_error_V2=primary['RMSE_V']**2,squared_bias_V2=primary['bias_V']**2,
        step_offset_difference_V=by_step['1']['bias_V']-by_step['2']['bias_V'],
        cross_validation=cv,uncertainty=intervals,signed_wave_spacings=report['signed_wave_spacings'],
        all_spacing_magnitudes_too_large=all(abs(r['calc_V'])>abs(r['exp_V']) for r in report['signed_wave_spacings']),
        organic_references=[organic_reference(rows,key) for key in anchors],
        reference=dict(computed_Fc_absolute_V=fc,stored_other_method_absolute_V=config['fc_abs_ref_MeCN_V'],
            difference_from_computed_V={k:v-fc for k,v in config['fc_abs_ref_MeCN_V'].items()},
            stored_values_status='Cross-method OROP metadata; not level-matched experimental replacements',
            unverified_SHE_conversion_constants_not_used=True,
            limiting_extra_twofold_orbital_partition=dict(Fc_free_energy_change_eV=orbital,
                Fc_absolute_potential_change_V=orbital,target_vs_Fc_potential_change_V=-orbital,
                status='Sensitivity only; orbital/vibronic/solvent splitting is not modeled; this shift would increase positive target bias')),
        common_reference_shift_changes_step_spacing=False,
        interpretation='Fc error is testable but not identified uniquely; state/solvation/ion-pairing errors can coexist',
        primary_potentials_modified=False,external_validation_claimed=False,
        source_sha256={str(p.resolve()):file_hash(p) for p in sources})
    args.out.mkdir(parents=True,exist_ok=False);atomic_write_json(args.out/'summary.json',result)
    with (args.out/'points.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(primary=primary,by_step=by_step,cv=[dict(holdout=r['holdout'],model=r['model'],**r['metrics']) for r in cv]),indent=2))


if __name__=='__main__':main()
