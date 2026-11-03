#!/usr/bin/env python3
"""Wait for capacity for one source-bound saved ionic Hessian validation.

Scientific admission is specific to this recovery. The shared prepared-namespace
coordinator preserves ownership, source, resource and prelaunch-refusal rules.
"""
import argparse
import json
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.workflows.prelaunch_accounting import nonresource_prelaunch_failures as noncontention_prelaunch_failures
from redox.workflows.prepared_namespace_queue import run


def validate(row):
    namespace=Path(row['namespace']);spec=json.loads(namespace.read_text())
    path=Path(spec['case'])/'plan.json';plan=json.loads(path.read_text())
    if (sha256_file(path)!=spec['plan_sha256'] or not plan.get('charged_saved_Hessian_validation')
            or not plan.get('saved_full_Hessian')):
        raise ValueError('Only a prepared source-bound ionic saved-matrix validation is admitted')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('spec','out','storage-peers'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--python',required=True);parser.add_argument('--hours',type=float,default=12)
    parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    namespace=args.spec.resolve();spec=json.loads(namespace.read_text())
    plan=json.loads((Path(spec['case'])/'plan.json').read_text())
    row=dict(key=plan['prepared_case']['key'],namespace=str(namespace),namespace_sha256=sha256_file(namespace))
    result=run([row],args.out.resolve(),args.storage_peers.resolve(),args.python,args.hours,
        finished_status='charged_polymer_native_component_finished_review_pending',validate=validate,
        sources=[Path(__file__).resolve()],resume=args.resume,label='saved_ionic_validation')
    print(json.dumps(dict(status=result['status'],scientific_retries_performed=0)))


if __name__=='__main__':main()
