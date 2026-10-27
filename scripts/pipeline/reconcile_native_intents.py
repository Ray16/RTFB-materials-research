#!/usr/bin/env python
"""Reconcile only source-bound unstarted native-stage placement intents."""
import argparse
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

FORBIDDEN=('launch.json','dispatch_attempts','experimental_native','residency_guard','host_guard','tensor_runtime')


def no_attempt(case,raw):
    return not raw.exists() and not any((case/name).exists() for name in FORBIDDEN)


def collect(audit):
    cases={}
    for row in audit['batches']:
        batch=Path(row['batch']);manifest=json.loads((batch/'batch.json').read_text())
        api.verify_sources(manifest['source_sha256'])
        if manifest['stage']!='native' or api.digest(batch/'batch.json')!=row['batch_sha256']:
            raise ValueError('Audited native batch changed')
        state=json.loads((batch/'dispatch.json').read_text())
        held={e['key'] for e in state['entries'] if e['status']=='submission_ownership_unverified_preserved'}
        for entry in manifest['entries']:
            if entry['key'] not in held:continue
            case=Path(entry['campaign']);plan=api.read_plan(case);raw=Path(plan['raw_root'])
            intent_bytes=(case/'dispatch_intent.json').read_bytes();intent=json.loads(intent_bytes)
            if (entry['existing_checked_reference'] or intent.get('outcome') is not None or not no_attempt(case,raw)
                    or plan.get('primary_candidate_admission') is not False
                    or plan.get('experimental_method_id')!=manifest['experimental_method_id']):
                raise ValueError('Held native case is not an exact unstarted input')
            api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
            if entry['key'] in cases:raise ValueError('Duplicate native case')
            cases[entry['key']]=dict(case=case,raw=raw,intent_bytes=intent_bytes,intent=intent,plan_sha256=api.digest(case/'plan.json'))
    return cases


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    recovery=api.load_module('native_intent_all_host_reconciliation',ROOT/'scripts/pipeline/reconcile_geometry_intents.py')
    recovery.collect=collect;recovery.no_attempt=no_attempt
    result=recovery.reconcile(a.audit.resolve(),a.out.resolve())
    print(json.dumps(dict(status=result['status'],reopened=result['reopened'],scope='Unstarted native intents; no computed G supplied')))
