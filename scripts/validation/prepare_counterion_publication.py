#!/usr/bin/env python3
"""Freeze checked primary chloride components and their matched scalar evidence."""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.counterion_thermochemistry import checked_component
from redox.workflows.molecular_native_namespace import load,verify


def prepare(summary_path,primary_inputs,duplicate_review,output):
    if output.exists():raise FileExistsError('Preserve all previous evidence')
    summary=json.loads(summary_path.read_text())
    if (summary['checked_components'] not in (4,6) or summary['distinct_original_starts_with_checked_G']!=summary['checked_components']-1
            or summary['declared_original_starts']!=7 or summary['conformer_ensemble_qualified']
            or summary['populations_inferred'] or summary['global_thermodynamic_pose_selected']):
        raise ValueError('Four checked finite outcomes and three original starts required')
    parents=json.loads(primary_inputs.read_text());parent=parents['components']['graft_ox2']
    parent_paths={Path(p) for p in parent['sources']}
    for path,digest in parent['sources'].items():
        if sha256_file(path)!=digest:raise ValueError('Primary parent source changed')
    reviews=[json.loads(p.read_text()) for p in parent_paths]
    parent_review=next(r for r in reviews if r.get('status')=='primary_molecular_native_control_independently_reconstructed')
    spec=verify(Path(parent_review['namespace']));case=Path(spec['case'])
    result=json.loads((case/'result.json').read_text());native=Path(result['native_record_path'])
    partitions=[p for p in (case/'state/projected_thermal').glob('*.json')
        if json.loads(p.read_text()).get('identity',{}).get('raw_native_sha256')==sha256_file(native)]
    if len(partitions)!=1:raise ValueError('One exact primary projected thermal definition required')
    thermal=json.loads(partitions[0].read_text())['identity']['protocol']
    output.mkdir(parents=True)
    shutil.copyfile(primary_inputs,output/'counterion_primary_inputs.json')
    shutil.copyfile(duplicate_review,output/'counterion_duplicate_review.json')
    atomic_write_json(output/'counterion_parent_thermal_protocol.json',thermal)
    identities=[('N0_minus',r'N0 above, $-$ mode'),('N0_plus',r'N0 above, $+$ mode'),
                ('N1_above','N1 above'),('N1_below','N1 below')]
    if summary['checked_components']==6:identities.extend([('N0_below_minus',r'N0 below, $-$ mode'),('N0_outer_minus',r'N0 outer, $-$ mode')])
    if len(summary['records'])!=len(identities):raise ValueError('Complete declared publication outcomes required')
    outcomes=[];sources={str(p):sha256_file(p) for p in [summary_path,primary_inputs,duplicate_review,partitions[0],Path(__file__)]}
    for item,(key,label) in zip(summary['records'],identities):
        review=Path(item['review']);ledger=Path(item['ledger'])
        if sha256_file(review)!=item['review_sha256'] or sha256_file(ledger)!=item['ledger_sha256']:
            raise ValueError('Summary component or matched ledger changed')
        report,spec=checked_component(review)
        if report['thermal_protocol']!=thermal:raise ValueError('Pair and actual primary parent thermal definitions differ')
        geometry=Path(spec['case'])/'state/opt.xyz'
        if sha256_file(geometry)!=report['geometry_sha256']:raise ValueError('Actual checked contact geometry differs')
        row=dict(key=key,label=label,review='counterion_review_'+key+'.json',
                 geometry='counterion_geometry_'+key+'.xyz',ledger='counterion_ledger_'+key+'.json')
        for source,name in [(review,row['review']),(geometry,row['geometry']),(ledger,row['ledger'])]:
            shutil.copyfile(source,output/name);sources[str(source)]=sha256_file(source)
        outcomes.append(row)
    checker=load(ROOT/'manuscript/check_counterion_association.py','frozen_counterion_publication_checker')
    catalog=dict(checked_outcomes=len(outcomes),declared_original_starts=7,outcomes=outcomes,
        conformer_ensemble_qualified=False,experimental_accuracy_validated=False,material_admission=False,
        reaction_rate_or_yield_inferred=False,source_sha256=sources)
    atomic_write_json(output/'counterion_inputs.json',catalog)
    rows,_=checker.reconstruct(output)
    (output/'counterion_values.tex').write_text(checker.values(rows))
    (output/'counterion_table.tex').write_text(checker.table(rows))
    report=checker.audit(output);atomic_write_json(output/'counterion_replay.json',report)
    shutil.copyfile(ROOT/'manuscript/check_counterion_association.py',output/'counterion_checker_source.py')
    shutil.copyfile(Path(__file__),output/'counterion_publication_source.py')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('summary','primary-inputs','duplicate-review','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(args.summary.resolve(),args.primary_inputs.resolve(),args.duplicate_review.resolve(),args.out.resolve())))
