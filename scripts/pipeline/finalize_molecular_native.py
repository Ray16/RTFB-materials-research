#!/usr/bin/env python3
"""Finish source-bound molecular analysis once every selected state is checked."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--hours',type=float,default=24);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    a.out.mkdir(parents=True,exist_ok=False);plan=json.loads(a.plan.read_text())
    prepared=Path(plan['prepared']);wanted=json.loads(prepared.read_text())['namespaces'];collections=[Path(p) for p in plan['collections']]
    paths=[a.plan.resolve(),prepared,Path(__file__).resolve(),ROOT/'scripts/pipeline/collect_molecular_native.py',
           ROOT/'scripts/pipeline/publish_molecular_native.py',ROOT/'scripts/plotting/candidates/plot_molecular_redox.py']
    proof={str(p):sha256_file(p) for p in paths}
    record=dict(status='waiting_for_selected_checked_molecular_states',owner=dict(host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),started=time.time(),
        deadline=time.time()+a.hours*3600,source_sha256=proof,material_admission=False)
    try:
        while time.time()<record['deadline']:
            if any(sha256_file(p)!=h for p,h in proof.items()):raise ValueError('Molecular finalization source/plan changed')
            observed={}
            for path in collections:
                if path.exists():
                    for row in json.loads(path.read_text())['entries']:
                        if row['namespace'] in wanted:observed[row['namespace']]=row['status']
            if any(status in ('failed_native_attempt_preserved','independent_review_failed_preserved') for status in observed.values()):
                record.update(status='selected_molecular_failure_preserved',observed=observed);return
            record.update(updated=time.time(),observed=observed)
            atomic_write_json(a.out/'status.json',record)
            if len(observed)==len(wanted) and all(s=='exact_molecular_native_G_independently_checked' for s in observed.values()):
                collector=load(ROOT/'scripts/pipeline/collect_molecular_native.py','final_exact_molecular_collection')
                checked=a.out/'checked_state_collection.json';collector.collect(prepared,collections,checked)
                publisher=load(ROOT/'scripts/pipeline/publish_molecular_native.py','final_exact_molecular_publication')
                result=publisher.publish(checked,Path(plan['reference']),Path(plan['archive']),Path(plan['registry']),
                    Path(plan['policy']),Path(plan['publication']),plan['electrode_pool'])
                plotting=load(ROOT/'scripts/plotting/candidates/plot_molecular_redox.py','final_exact_molecular_figure')
                source=Path(plan['publication'])/'candidate.json';figure=Path(plan['figure']);plot=plotting.plot(source,figure)
                pointer=Path(plan['current_pointer'])
                if not pointer.resolve().is_relative_to(ROOT/'results/candidates/verified_cases'):
                    raise ValueError('Current pointer must stay in the exact-candidate result collection')
                atomic_write_json(pointer,dict(candidate=str(source.relative_to(ROOT)),
                    numeric_inputs=str((source.parent/'numeric_inputs.json').relative_to(ROOT)),
                    figure=str(figure.relative_to(ROOT)),source_snapshot=str((a.out/'status.json').relative_to(ROOT)),
                    all_axis_decision=result['all_axis_decision']['status'],figure_visual_review='pending'))
                record.update(status='molecular_redox_outputs_numerically_checked_visual_review_pending',publication=plan['publication'],
                    candidate_sha256=sha256_file(source),figure=plan['figure'],figure_sha256=plot['figure_sha256'],
                    current_pointer=str(pointer),all_axis_decision=result['all_axis_decision']['status'])
                return
            time.sleep(20)
        record['status']='molecular_finalization_window_ended_source_workers_untouched'
    except BaseException as error:
        record.update(status='molecular_finalization_failed_preserved',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(a.out/'status.json',record)


if __name__=='__main__':main()
