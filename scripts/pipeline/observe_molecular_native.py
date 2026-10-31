#!/usr/bin/env python3
"""Finite, source-bound collection of exact molecular native/resource outcomes."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.molecular_native import review,ReviewPending
from redox.workflows.molecular_native_namespace import verify


def execution_closed(case):
    """Wait for all closing records, including delayed shared-file visibility."""
    paths=[case/'residency_guard/owner.json',case/'host_guard/status.json',case/'tensor_runtime/status.json']
    if any(not p.exists() for p in paths):return False
    guard,host,runtime=[json.loads(p.read_text()) for p in paths]
    return guard.get('guard_stopped') is True and 'finished' in host and 'finished' in runtime


def checked_cached(output,namespace):
    result=json.loads(output.read_text())
    if result['namespace_sha256']!=sha256_file(namespace) or result['namespace']!=str(namespace):
        raise ValueError('Cached molecular review belongs to another exact namespace')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--hours',type=float,default=24);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    a.out.mkdir(parents=True,exist_ok=False)
    paths=[Path(v) for v in json.loads(a.prepared.read_text())['namespaces']]
    if not paths or len(paths)!=len(set(paths)):raise ValueError('Distinct exact-state namespaces required')
    definition=importlib.util.spec_from_file_location('molecular_native_owner_check',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    source=[a.prepared.resolve(),*paths,Path(__file__).resolve(),ROOT/'src/redox/validation/molecular_native.py',
        ROOT/'src/redox/workflows/molecular_native_namespace.py',ROOT/'src/redox/validation/solution_pair_review.py']
    proof={str(p):sha256_file(p) for p in source}
    record=dict(status='molecular_native_outcomes_pending',owner=dict(host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),source_sha256=proof,
        started=time.time(),deadline=time.time()+a.hours*3600,native_free_energies_inferred_from_process_exit=False)
    done={}
    while time.time()<record['deadline']:
        if any(sha256_file(p)!=h for p,h in proof.items()):raise ValueError('Molecular observer source/inputs changed')
        rows=[]
        for path in paths:
            if str(path) in done:rows.append(done[str(path)]);continue
            spec=verify(path);case=Path(spec['case']);stage=case/'status.json'
            row=dict(namespace=str(path),state=spec['state'],status='prepared_not_started')
            if stage.exists():
                state=json.loads(stage.read_text());row['status']='running_or_ownership_unverified'
                if 'finished' in state:
                    if state['status']!='molecular_native_state_finished_independent_review_pending':
                        row.update(status='failed_native_attempt_preserved',error=state.get('error'));done[str(path)]=row
                    else:
                        try:
                            if owners.alive(state['owner']):row['status']='finished_owner_closing'
                            elif not execution_closed(case):row['status']='finished_execution_records_closing'
                            else:
                                output=a.out/(sha256_file(path)[:16]+'__'+spec['source_id']+'__'+spec['state']+'_independent_review.json')
                                result=checked_cached(output,path) if output.exists() else review(path,output)
                                row.update(status=result['status'],review=str(output),review_sha256=sha256_file(output),
                                    G_by_cutoff_eV=result['G_by_cutoff_eV']);done[str(path)]=row
                        except (OSError,subprocess.SubprocessError,ReviewPending) as error:
                            row.update(status='review_observation_unverified_preserved',error=repr(error))
                        except ValueError as error:
                            row.update(status='independent_review_failed_preserved',error=repr(error));done[str(path)]=row
            rows.append(row)
        record.update(updated=time.time(),entries=rows,terminal=len(done),expected=len(paths),
            checked=sum(r['status']=='exact_molecular_native_G_independently_checked' for r in rows))
        atomic_write_json(a.out/'status.json',record)
        if len(done)==len(paths):
            record.update(status='molecular_native_terminal_outcomes_separate',finished=time.time())
            atomic_write_json(a.out/'status.json',record);return
        time.sleep(20)
    record.update(status='molecular_native_observation_window_ended_workers_untouched',finished=time.time())
    atomic_write_json(a.out/'status.json',record)


if __name__=='__main__':main()
