#!/usr/bin/env python
"""Independent record/source/count checks for the current all-axis publication.

Does not repeat electronic calculations or turn the absence of a shortlisted
candidate into a chemical failure. Positive/negative rule behavior is covered
by separate adversarial tests; this replay checks the actual published records.
"""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
from rdkit import Chem
from redox.core.common import ROOT
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False);pointer=json.loads((ROOT/'results/selection/current.json').read_text())
    folder=ROOT/pointer['report'];summary=json.loads((folder/'summary.json').read_text());manifest=json.loads((folder/'publication_manifest.json').read_text())
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    for path,h in manifest['source_sha256'].items():
        if sha(path)!=h:raise ValueError('Current source changed')
    for path,h in manifest['output_sha256'].items():
        if sha(ROOT/path)!=h:raise ValueError('Published output changed')
    data=json.loads((folder/'decisions.json').read_text())['candidates'];policy=json.loads((folder/'source_snapshot/policy.json').read_text())
    feed=json.loads((folder/'source_snapshot/feed.json').read_text());pool=json.loads((folder/'source_snapshot/pool.json').read_text())
    keys=set()
    for row in pool:
        keys.add(Chem.MolToInchiKey(Chem.MolFromSmiles(row['grafted_smiles'])))
    for row in feed['original39']+feed['additional_candidates']+feed.get('comparison_controls',[]):
        keys.add(Chem.MolToInchiKey(Chem.MolFromSmiles(row['exact_grafted_smiles'])))
    if keys!={r['product_key'] for r in data} or len(keys)!=len(data):raise ValueError('Exact graph union or deduplication differs')
    failed=pending=controls=eligible=0
    for row in data:
        if Chem.MolToInchiKey(Chem.MolFromSmiles(row['grafted_smiles']))!=row['product_key']:raise ValueError('Decision graph key differs')
        f=[r['axis'] for r in row['axes'] if r['status']=='fail'];p=[r['axis'] for r in row['axes'] if r['status']=='pending']
        if f!=row['failed_axes'] or p!=row['pending_axes'] or row['all_required_axes_pass']!=(not f and not p):raise ValueError('Decision conjunction differs')
        is_eligible=not f and not p and row['role']=='candidate'
        if row['eligible_for_Pareto']!=is_eligible:raise ValueError('A nonpassing/control row was admitted')
        eligible+=is_eligible;controls+=row['role']=='comparison_control';failed+=bool(f);pending+=bool(p) and row['role']=='candidate'
    unset=[r['id'] for r in policy['required_axes'] if r['kind']=='numeric' and r.get('minimum') is None and r.get('maximum') is None]
    if unset!=summary['unset_numeric_cutoffs']:raise ValueError('Unconfigured targets were hidden')
    if eligible!=len(summary['eligible_ids']) or eligible!=pointer['eligible']:raise ValueError('Eligibility count differs')
    old=json.loads(subprocess.check_output(['git','show','9ba043c:results/reviews/2026-10-18/verified_candidates/latest.json'],cwd=ROOT,text=True))
    before={r['exact_grafted_inchikey']:r for r in old['original39']+old['additional_candidates']}
    changes=[]
    for row in feed['original39']+feed['additional_candidates']:
        previous=before[row['exact_grafted_inchikey']]
        for field in ('E1_vs_Fc_V','E2_vs_Fc_V','capacity_nominal_mAh_g','accessible_capacity_mAh_g'):
            a,b=previous.get(field),row.get(field)
            if a!=b and (a is None or b is None or abs(a-b)>1e-9):changes.append(dict(id=row['id'],field=field,old=a,new=b))
    report=dict(status='all_axes_published_records_and_sources_checked',checked=time.time(),exact_products=len(keys),
        pending_candidate_designs=pending,comparison_controls=controls,failed_designs=failed,eligible=eligible,
        unset_numeric_targets=unset,candidate_feed_as_of=feed['as_of_Chicago'],candidate_arithmetic_checks=feed['summary']['arithmetic_checks'],
        candidate_arithmetic_discrepancies=feed['summary']['arithmetic_discrepancies'],changed_capped_property_values_since_previous_commit=changes,
        source_current=True,record_logic_checked=True,missing_scientific_properties_validated=False,
        source_sha256={str(p):sha(p) for p in (Path(__file__),ROOT/'results/selection/current.json',folder/'summary.json',folder/'publication_manifest.json')})
    (args.out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'}))


if __name__=='__main__':main()
