#!/usr/bin/env python3
"""Apply the shared all-axis policy to an exact route-held design panel."""
import argparse
import hashlib
import json
from pathlib import Path

from redox.screening.feasibility import assess_candidate


def assess(matrix_path,gate_path,policy_path):
    matrix=json.loads(matrix_path.read_text());gate=json.loads(gate_path.read_text());policy=json.loads(policy_path.read_text())
    rows=matrix['candidates'];screened={row['id']:row for row in gate['candidates']}
    if not rows or len(screened)!=len(gate['candidates']) or len({r['id'] for r in rows})!=len(rows):
        raise ValueError('Unique candidate and route identities required')
    if len({r['product_inchikey'] for r in rows})!=len(rows):raise ValueError('Duplicate exact products')
    rule=next(a for a in policy['required_axes'] if a['id']=='attachment_route');decisions=[]
    for row in rows:
        if row['id'] not in screened:raise ValueError('A candidate lacks its route decision')
        route=screened.pop(row['id'])
        if (route['route_status']!='hold' or route['proceed_to_full_redox']
                or route['evidence'].get('exact_design_key')!=row['exact_design_key']
                or route['evidence'].get('route_binding')!=row['route_binding']):
            raise ValueError('Route evidence does not bind the held exact candidate')
        key=row['product_inchikey'];axis=dict(product_key=key,scope=rule['scope'],status='partial',reason='Candidate route has unresolved evidence requirements')
        decision=assess_candidate(dict(id=row['id'],product_key=key,role='candidate',axes={'attachment_route':axis}),policy)
        if decision['status']!='pending' or decision['all_required_axes_pass']:
            raise ValueError('Partial route evidence cannot establish material eligibility')
        decisions.append(decision)
    if screened:raise ValueError('Unexpected route-screen candidate')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    return dict(status='matched_panel_all_axis_decisions_pending',decisions=decisions,
                source_sha256={str(p):sha(p) for p in (matrix_path,gate_path,policy_path,Path(__file__).resolve())},
                quantum_jobs_submitted=0,primary_decision_ledger_modified=False,
                scope='Conjunctive material policy; no proposed graph or route is promoted to a qualified material axis')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('matrix','gate','policy','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=assess(a.matrix.resolve(),a.gate.resolve(),a.policy.resolve())
    with a.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(status=result['status'],pending=len(result['decisions']))))
