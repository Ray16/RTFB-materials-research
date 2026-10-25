#!/usr/bin/env python
"""Restart interrupted smooth-model geometries from durable unqualified frames."""
import argparse
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(sources, out):
    owner = api.load_module('resident_recovery_stopped_owner', ROOT/'scripts/pipeline/run_polymer_spacing.py')
    audit = api.load_module('resident_recovery_exact_graph', ROOT/'scripts/validation/audit_polymer_geometries.py')
    worker = ROOT/'scripts/pipeline/run_smooth_geometry.py'
    coordinator = ROOT/'scripts/pipeline/coordinate_smooth_geometry.py'
    declarations = []
    for source in sources:
        plan = api.read_plan(source)
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
        status_path = source/'geometry/status.json'; state = json.loads(status_path.read_text())
        interruption = source/'residency_guard/interruption.json'; collision = json.loads(interruption.read_text())
        if (state['status']!='geometry_stage_failed_preserved' or owner.owner_alive(state['owner'])
                or collision['pid']!=state['owner']['pid'] or not collision['foreign_resident_PIDs']
                or collision['foreign_processes_signaled'] or len(plan['tasks'])!=1):
            raise ValueError('A stopped, preserved exclusive-residency interruption is required')
        sd = Path(state['state_dir']); xyz = sd/'geometry_checkpoints/last_evaluated.xyz'
        latest_path=sd/'geometry_checkpoints/latest.json';latest=json.loads(latest_path.read_text())
        if not latest['SCF_converged'] or latest['solution_geometry_qualified']:
            raise ValueError('The warm-start frame must be explicitly unqualified and SCF converged')
        task=plan['tasks'][0];graph=audit.audit(task,xyz)
        if graph['issues']:raise ValueError('Interrupted checkpoint graph/stereochemistry failed')
        declarations.append((source,plan,state,status_path,interruption,xyz,latest_path,latest,graph))
    out.mkdir(parents=True,exist_ok=False);entries=[]
    for source,old,state,status_path,interruption,xyz,latest_path,latest,graph in declarations:
        key=old['tasks'][0]['key'];case=out/key;case.mkdir()
        shutil.copytree(source/'live_inputs',case/'live_inputs')
        ready=case/'live_inputs'/key;target=ready/'input.xyz';shutil.copyfile(xyz,target)
        row=json.loads((ready/'task.json').read_text());row['input_xyz']=str(target)
        row['input_sha256'][str(target)]=api.digest(target)
        row['source_restart_checkpoint']=dict(path=str(xyz),evaluation=latest['evaluation'],qualified=False)
        atomic_write_json(ready/'task.json',row);atomic_write_json(ready/'initial_graph_review.json',graph)
        inputs=dict(old['input_sha256'])
        for p in (source/'plan.json',status_path,interruption,xyz,latest_path,target,
                  ready/'task.json',ready/'initial_graph_review.json'):inputs[str(p)]=api.digest(p)
        producers=dict(old['producer_sha256']);producers[str(Path(__file__).resolve())]=api.digest(__file__)
        for stage in ('energy','thermal'):
            for name in ('tasks','done','claims','failures','workers','pilot','launches'):
                (case/stage/name).mkdir(parents=True)
        plan=dict(old,input_sha256=inputs,producer_sha256=producers,
            raw_root=str(ROOT/'calcs/polymer_models'/out.name/key),
            source_interrupted_campaign=str(source),restart_checkpoint=str(xyz),
            restart_evaluation=latest['evaluation'],source_checkpoint_is_qualified=False,
            purpose='Restart same model/graph from a residency-interrupted unqualified frame; repeat full geometry QC',
            created=time.time())
        atomic_write_json(case/'plan.json',plan)
        entries.append(dict(key=key,campaign=str(case),plan_sha256=api.digest(case/'plan.json')))
    files=[Path(__file__).resolve(),worker,coordinator,ROOT/'scripts/pipeline/coordinate_geometry_batch.py',
           ROOT/'scripts/fleet/preflight_cpu_dependencies.py',ROOT/'scripts/pipeline/coordinate_with_cpu_preflight.py',
           ROOT/'environments/cutensor22/runtime.json']
    atomic_write_json(out/'batch.json',dict(status='residency_interrupted_geometry_restarts_prepared',
        stage='geometry',entries=entries,worker=str(worker),maximum_RSS_GiB=32,minimum_free_RAM_GiB=32,
        minimum_GPU_MiB=32768,experimental_model=declarations[0][1]['experimental_model'],
        primary_candidate_admission=False,quantum_results_transferred=False,
        source_sha256={str(p):api.digest(p) for p in files}))
    print(json.dumps(dict(prepared=len(entries),old_attempts_preserved=True,qualified_geometries_reused=False)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,action='append',required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    prepare([p.resolve() for p in args.source],args.out.resolve())


if __name__=='__main__':main()
