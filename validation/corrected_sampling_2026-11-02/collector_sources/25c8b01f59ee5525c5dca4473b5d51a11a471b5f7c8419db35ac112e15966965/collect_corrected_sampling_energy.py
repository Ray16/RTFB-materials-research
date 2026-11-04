#!/usr/bin/env python3
"""Review finished corrected inputs and freeze full-roster, field-specific status.

No job is restarted. Missing E/G remains null; a completed producer is distinct
from an independently reviewed electronic state and from native thermochemistry.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import time
import hashlib

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.sampling_receipts import checked_energy as checked_report
from redox.validation.sampling_receipts import SOURCE_BYTES as RECEIPT_SOURCE_BYTES

SOURCE_BYTES=Path(__file__).read_bytes()


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def collect(manifest_path,reviews,output):
    if output.exists():raise FileExistsError('Preserve earlier complete-roster snapshots')
    manifest=json.loads(manifest_path.read_text());rows=manifest['states']
    if manifest['target_states']!=240 or len(rows)!=240 or len({r['key'] for r in rows})!=240:
        raise ValueError('All 240 distinct declared input slots must be retained')
    reviews.mkdir(parents=True,exist_ok=True)
    reviewer=load(ROOT/'scripts/validation/review_corrected_sampling_energy.py','cohort_actual_electronic_review')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','cohort_electronic_owners')
    records=[];errors=[];sources={str(manifest_path):sha256_file(manifest_path)}
    for name,content in [(Path(__file__).name,SOURCE_BYTES),('sampling_receipts.py',RECEIPT_SOURCE_BYTES)]:
        digest=hashlib.sha256(content).hexdigest();frozen=output.parent/'collector_sources'/digest/name
        if not frozen.exists():frozen.parent.mkdir(parents=True,exist_ok=True);frozen.write_bytes(content)
        if sha256_file(frozen)!=digest:raise ValueError('As-used collection source differs')
        sources[str(frozen)]=digest
    for row in rows:
        namespace=Path(row['namespace']);spec=json.loads(namespace.read_text());case=Path(spec['case']);key=row['key']
        if sha256_file(namespace)!=row['namespace_sha256']:raise ValueError('Prepared input namespace changed')
        state=dict(key=key,namespace=str(namespace),E_smd_eV=None,E_gas_eV=None,G_eV=None,native_G_qualified=False)
        receipt=reviews/(key+'.json');status_path=case/'status.json'
        try:
            if receipt.exists():
                report=checked_report(receipt,namespace,key,manifest['experimental_method_id'])
            elif status_path.exists():
                stage=json.loads(status_path.read_text())
                if stage['status']=='corrected_sampling_energy_finished_independent_review_pending':
                    live=owners.alive(stage['owner'])
                    if live is not False or not list(case.parent.glob('placement*/completion.json')):
                        state['status']='execution_closure_pending';records.append(state);continue
                    report=reviewer.review(namespace);atomic_write_json(receipt,report)
                    report=checked_report(receipt,namespace,key,manifest['experimental_method_id'])
                else:
                    state.update(status=stage['status'],error=stage.get('error'));records.append(state);continue
            else:
                state['status']='not_scientifically_started';records.append(state);continue
            state.update(status='electronic_independently_checked',review=str(receipt),review_sha256=sha256_file(receipt),
                E_smd_eV=report['E_smd_eV'],E_gas_eV=report['E_gas_eV'])
            sources[str(receipt)]=sha256_file(receipt)
        except Exception as error:
            state.update(status='independent_review_failed_preserved',error=repr(error));errors.append(dict(key=key,error=repr(error)))
        records.append(state)
    result=dict(status='complete_roster_corrected_electronic_review_snapshot',checked=time.time(),declared_inputs=240,
        checked_electronic_inputs=sum(r['status']=='electronic_independently_checked' for r in records),
        states=records,independent_review_errors=errors,experimental_method_id=manifest['experimental_method_id'],
        source_sha256=sources,native_G_qualified=False,occupied_ensemble_cycles_qualified=False,
        pending_values_are_null=True,jobs_restarted=False)
    atomic_write_json(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('manifest','reviews','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=collect(a.manifest.resolve(),a.reviews.resolve(),a.out.resolve())
    print(json.dumps(dict(status=result['status'],checked_inputs=result['checked_electronic_inputs'],errors=result['independent_review_errors'])))
