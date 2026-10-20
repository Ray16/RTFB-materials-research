#!/usr/bin/env python
"""Publish immutable checked covalent warm starts as individual ensembles finish."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import time

from prepare_pi_stacking import BASE,choose
from redox.core.io import atomic_write_json
from redox.build.geometry import geometry_stereo_status
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    owner=BASE/'conformer_stager.json'
    assert not owner.exists(),'Inspect original stager ownership before restarting'
    plan=json.loads((BASE/'plan.json').read_text());plan_hash=sha(BASE/'plan.json')
    report=dict(host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],started=time.time(),status='watching_completed_ensembles')
    (BASE/'input_holds').mkdir(exist_ok=True)
    tasks=[t for t in plan['tasks'] if t['model_type']=='covalent_oligomer' and t['repeats']>1]
    deadline=time.time()+4*3600
    try:
        while time.time()<deadline:
            assert sha(BASE/'plan.json')==plan_hash
            staged=0;held=0
            for task in tasks:
                directory=BASE/'live_inputs'/task['key'];receipt=directory/'task.json'
                hold=BASE/'input_holds'/(task['key']+'.json')
                if receipt.exists():staged+=1;continue
                if hold.exists():held+=1;continue
                source=BASE/'seeds'/task['model_key']/'ensemble.json'
                if not source.exists():continue
                try:
                    ensemble=json.loads(source.read_text())
                    record,fallback=choose(ensemble,task['target'])
                    path=Path(record['xyz'])
                    stereo=geometry_stereo_status(task['smiles'],path)
                    topology=hydrogen_topology_status(task['smiles'],path)
                    if stereo not in ('ok','not_applicable') or topology!='ok':
                        raise ValueError('Warm-start structural gate: '+stereo+'; '+topology)
                    directory.mkdir(exist_ok=True)
                    destination=directory/'input.xyz';shutil.copyfile(path,destination)
                    ready=dict(task,input_xyz=str(destination),pi_groups=record['pi_groups'],
                        seed_stage='biased_FF_start; DFT_required' if fallback else 'converged_topology_checked_xTB_ALPB',
                        input_stereochemistry_QC=stereo,
                        input_sha256={str(destination):sha(destination),str(source):sha(source),str(path):sha(path)})
                    atomic_write_json(receipt,ready);staged+=1
                except (ValueError,AssertionError) as exc:
                    atomic_write_json(hold,dict(key=task['key'],stage='conformer_preparation',error=repr(exc),actual_DFT_attempts=0));held+=1
            report.update(staged_covalent_cases=staged,preparation_holds=held,expected_covalent_cases=len(tasks),updated=time.time())
            if staged+held==len(tasks):report['status']='finite_conformer_staging_terminal';atomic_write_json(owner,report);return
            atomic_write_json(owner,report);time.sleep(10)
        report['status']='bounded_wait_ended_remaining_preparation_preserved'
    except BaseException as exc:report.update(status='stopped_preserving_inputs',error=repr(exc));raise
    finally:report['finished']=time.time();atomic_write_json(owner,report)


if __name__=='__main__':main()
