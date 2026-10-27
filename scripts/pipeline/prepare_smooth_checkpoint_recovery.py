#!/usr/bin/env python
"""Prepare isolated restarts of stopped geometry I/O/residency failures."""
import argparse
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(sources,out):
    if out.exists():raise FileExistsError('Preserve previous attempts')
    owner=api.load_module('geometry_checkpoint_recovery_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    auditor=api.load_module('geometry_checkpoint_recovery_graph',ROOT/'scripts/validation/audit_polymer_geometries.py')
    ready=[]
    for source in sources:
        plan=api.read_plan(source);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        status=source/'geometry/status.json';state=json.loads(status.read_text())
        if state['status']!='geometry_stage_failed_preserved' or owner.owner_alive(state['owner']) or len(plan['tasks'])!=1:
            raise ValueError('Require a stopped, preserved single-geometry failure')
        launch_path=source/'launch.json';launch=json.loads(launch_path.read_text())
        if owner.owner_alive(launch['owner']):raise ValueError('Original gate owner is still live')
        error=state.get('error','');lock_error='BlockingIOError' in error and 'lock file' in error
        interruption=source/'residency_guard/interruption.json'
        if not lock_error:
            collision=json.loads(interruption.read_text())
            if (not collision['foreign_resident_PIDs'] or collision['foreign_processes_signaled']
                    or collision['pid']!=state['owner']['pid']):raise ValueError('Failure is not a declared I/O or residency interruption')
        task=plan['tasks'][0];directory=Path(state['state_dir'])/'geometry_checkpoints';latest=directory/'latest.json'
        if latest.exists():
            point=json.loads(latest.read_text());xyz=directory/'last_evaluated.xyz'
            if point['SCF_converged'] is not True or point['solution_geometry_qualified'] is not False:
                raise ValueError('Saved warm start must remain SCF-converged and explicitly unqualified')
            evaluation=point['evaluation']
        else:xyz=source/'live_inputs'/task['key']/'input.xyz';evaluation=None
        graph=auditor.audit(task,xyz)
        if graph['issues']:raise ValueError('Restart graph or stereochemistry changed')
        proof=[source/'plan.json',status,launch_path,xyz]
        if latest.exists():proof.append(latest)
        if interruption.exists():proof.append(interruption)
        ready.append((source,plan,task,xyz,evaluation,graph,proof))
    if len({r[2]['key'] for r in ready})!=len(ready):raise ValueError('Distinct failed inputs required')
    out.mkdir(parents=True);entries=[]
    producer=ROOT/'scripts/pipeline/run_smooth_checkpoint_recovery.py'
    for source,old,task,xyz,evaluation,graph,proof in ready:
        key=task['key'];case=out/key;case.mkdir();shutil.copytree(source/'live_inputs',case/'live_inputs')
        directory=case/'live_inputs'/key;target=directory/'input.xyz';shutil.copyfile(xyz,target)
        metadata=json.loads((directory/'task.json').read_text());metadata['input_xyz']=str(target);metadata['input_sha256'][str(target)]=api.digest(target)
        metadata['source_restart_checkpoint']=dict(path=str(xyz),evaluation=evaluation,qualified=False)
        atomic_write_json(directory/'task.json',metadata);atomic_write_json(directory/'initial_graph_review.json',graph)
        for stage in ('energy','thermal'):
            for name in ('tasks','done','claims','failures','workers','pilot','launches'):(case/stage/name).mkdir(parents=True)
        inputs={**old['input_sha256'],**{str(p):api.digest(p) for p in [*proof,target,directory/'task.json',directory/'initial_graph_review.json']}}
        producers={**old['producer_sha256'],**{str(p):api.digest(p) for p in (Path(__file__).resolve(),producer,ROOT/'src/redox/workflows/checkpoint_io.py')}}
        fresh=dict(old,input_sha256=inputs,producer_sha256=producers,raw_root=str(ROOT/'calcs/polymer_models'/out.name/key),
                   execution_driver=str(producer),source_interrupted_campaign=str(source),restart_checkpoint=str(xyz),
                   source_checkpoint_is_qualified=False,checkpoint_IO_policy='bounded EAGAIN HDF5-lock retries; locking retained',created=time.time())
        atomic_write_json(case/'plan.json',fresh);entries.append(dict(key=key,campaign=str(case),plan_sha256=api.digest(case/'plan.json')))
    watched=(Path(__file__).resolve(),producer,ROOT/'src/redox/workflows/checkpoint_io.py',ROOT/'scripts/pipeline/run_smooth_geometry.py',ROOT/'environments/cutensor22/runtime.json')
    manifest=dict(status='stopped_geometry_checkpoint_restarts_prepared',stage='geometry',entries=entries,worker=str(producer),
                  maximum_RSS_GiB=32,minimum_free_RAM_GiB=32,minimum_GPU_MiB=32768,
                  experimental_model=ready[0][1]['experimental_model'],primary_candidate_admission=False,
                  quantum_results_transferred=False,source_sha256={str(p):api.digest(p) for p in watched})
    atomic_write_json(out/'batch.json',manifest);return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    m=prepare([s.resolve() for s in a.source],a.out.resolve());print(json.dumps(dict(prepared=len(m['entries']),qualified_coordinates_reused=False)))
