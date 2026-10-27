#!/usr/bin/env python3
"""Apply the shared all-axis policy to route-screened exact PMDI proposals."""
import argparse
import hashlib
import json
from pathlib import Path

from redox.screening.feasibility import assess_candidate


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assess(matrix_path, gate_path, policy_path):
    matrix = json.loads(matrix_path.read_text())
    gate = json.loads(gate_path.read_text())
    policy = json.loads(policy_path.read_text())
    screened = {row['id']:row for row in gate['candidates']}
    if len(screened) != len(gate['candidates']) or len(matrix['candidates']) != 4:
        raise ValueError('Expected four unique route-screened designs')
    route_rule = next(a for a in policy['required_axes'] if a['id']=='attachment_route')
    decisions = []
    for row in matrix['candidates']:
        route = screened.pop(row['id'])
        if (route['route_status'] != 'hold' or route['proceed_to_full_redox']
                or route['evidence'].get('exact_design_key') != row['exact_design_key']
                or route['evidence'].get('route_binding') != row['route_binding']):
            raise ValueError('Route evidence does not bind the held exact design')
        key = row['product_inchikey']
        axes = {'attachment_route':dict(product_key=key, scope=route_rule['scope'], status='partial',
                    reason='Exact precursor access, reaction conditions, selectivity and redox retention remain unreviewed')}
        decision = assess_candidate(dict(id=row['id'],product_key=key,role='candidate',axes=axes),policy)
        if decision['status'] != 'pending' or decision['all_required_axes_pass']:
            raise ValueError('Insufficient material evidence cannot be eligible')
        decisions.append(decision)
    if screened:
        raise ValueError('Unexpected route-screen candidate')
    return dict(status='matrix_all_axis_decisions_pending',decisions=decisions,
                source_sha256={str(p):sha(p) for p in (matrix_path,gate_path,policy_path,Path(__file__).resolve())},
                primary_decision_ledger_modified=False, quantum_jobs_submitted=0,
                scope='Conjunctive material policy; no descriptor or route proposal is substituted for a measured or computed material axis')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('matrix','gate','policy','out'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=assess(a.matrix.resolve(),a.gate.resolve(),a.policy.resolve())
    with a.out.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(status=result['status'],pending=len(result['decisions']))))
