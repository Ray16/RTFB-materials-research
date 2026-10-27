#!/usr/bin/env python3
"""Apply the shared material gate to the source-reported precursor control."""
import argparse
import hashlib
import json
from pathlib import Path

from redox.screening.feasibility import assess_candidate


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assess(candidate_path,gate_path,policy_path):
    rows=json.loads(candidate_path.read_text());screen=json.loads(gate_path.read_text())['candidates']
    policy=json.loads(policy_path.read_text())
    if len(rows)!=1 or len(screen)!=1:
        raise ValueError('Exactly one candidate and route decision required')
    row=rows[0];route=screen[0]
    if (row['id']!=route['id'] or route['route_status']!='hold' or route['proceed_to_full_redox']
            or route['evidence'].get('exact_design_key')!=row['exact_design_key']
            or route['evidence'].get('route_binding')!=row['route_binding']
            or not route['evidence'].get('exact_precursor_access_reviewed')):
        raise ValueError('Held route and source-identified precursor do not bind the exact candidate')
    rule=next(a for a in policy['required_axes'] if a['id']=='attachment_route')
    axis=dict(product_key=row['product_inchikey'],scope=rule['scope'],status='partial',
              reason='Precursor isolation reported; attachment selectivity and conditions unresolved')
    decision=assess_candidate(dict(id=row['id'],product_key=row['product_inchikey'],role='candidate',
                                   axes={'attachment_route':axis}),policy)
    if decision['status']!='pending' or decision['all_required_axes_pass']:
        raise ValueError('Partial precursor evidence cannot establish material eligibility')
    return dict(status='source_precursor_control_all_axis_pending',decision=decision,
                source_sha256={str(p):sha(p) for p in (candidate_path,gate_path,policy_path,Path(__file__).resolve())},
                quantum_jobs_submitted=0,primary_decision_ledger_modified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('candidate','gate','policy','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=assess(a.candidate.resolve(),a.gate.resolve(),a.policy.resolve())
    with a.out.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(status=result['decision']['status'],pending_axes=len(result['decision']['pending_axes']))))
