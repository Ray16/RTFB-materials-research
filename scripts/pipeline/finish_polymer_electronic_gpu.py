#!/usr/bin/env python
"""Canonical GPU energies after separate geometry and resource qualification."""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import json
import os
from pathlib import Path
import shutil
import signal
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability

GC=ROOT/'scripts/pipeline/smd_gc_workspace.py'
CHECKER=ROOT/'scripts/validation/qualify_gc_workspace.py'


def merge_hashes(*maps):
    result={}
    for mapping in maps:
        for key,value in mapping.items():
            if key in result and result[key]!=value:raise ValueError('Conflicting source epochs: '+key)
            result[key]=value
    return result


def prepare(geometry,pilot,review,controls,out):
    original,point=api.read_plan(geometry),api.read_plan(pilot)
    for p in (original,point):api.verify_sources(p['producer_sha256']);api.verify_sources(p['input_sha256'])
    parent_path=geometry/'geometry/status.json';parent=json.loads(parent_path.read_text())
    proof_path=review/'summary.json';proof=json.loads(proof_path.read_text())
    if parent['status']!='stationary_geometry_checked_active_parent_and_native_minimum_pending':
        raise ValueError('The optimizer, graph, stationarity and finite differences must pass first')
    if proof['status']!='electronic_resource_point_independently_reviewed' or proof['source']!=str(pilot):
        raise ValueError('The exact completed GPU resource point requires independent review')
    api.verify_sources(parent['output_sha256']);api.verify_sources(proof['source_sha256'])
    for key in ('tasks','quantum_settings','thermochemistry','solvent'):
        if original[key]!=point[key]:raise ValueError('Geometry and resource point belong to different tasks/protocols')
    if len(original['tasks'])!=1:raise ValueError('One exact stationary state is required')
    alive=api.load_module('GPU_energy_predecessor_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py').owner_alive
    if alive(parent['owner']):raise ValueError('The geometry producer must finish before handoff')
    witness=api.load_module('GPU_energy_GC_controls',CHECKER).checked_controls(controls)
    inputs=merge_hashes(original['input_sha256'],point['input_sha256'],parent['output_sha256'],proof['source_sha256'],witness)
    for p in (geometry/'plan.json',pilot/'plan.json',parent_path,proof_path):inputs[str(p)]=api.digest(p)
    sources=merge_hashes(original['producer_sha256'],point['producer_sha256'],witness)
    for p in (Path(__file__).resolve(),GC,CHECKER,ROOT/'src/redox/workflows/split_electronic.py'):sources[str(p)]=api.digest(p)
    plan=dict(original,raw_root=str(ROOT/'calcs/polymer_models'/out.parent.name/out.name),
        producer_sha256=sources,input_sha256=inputs,source_geometry_report=parent,
        GPU_resource_qualification=str(pilot),GPU_resource_review=str(proof_path),
        maximum_RSS_GiB=proof['host_RSS_bound_GiB'],minimum_free_RAM_GiB=32,
        maximum_peak_GPU_MiB=point['maximum_memory_pilot_MiB'],execution_workspace_GiB=point['execution_workspace_GiB'],
        energy_only_execution=True,geometry_only_execution=False,native_pipeline_enabled=False,
        externally_managed_stages=True,created=time.time(),
        purpose='Canonical GPU active SP on exact checked coordinates; native thermochemistry is separate')
    out.mkdir(parents=True,exist_ok=False);shutil.copytree(geometry/'controls',out/'controls')
    shutil.copytree(geometry/'live_inputs',out/'live_inputs')
    key=original['tasks'][0]['key'];ready_path=out/'live_inputs'/key/'task.json'
    ready=json.loads(ready_path.read_text());xyz=ready_path.parent/'input.xyz'
    ready['input_xyz']=str(xyz);ready['input_sha256'][str(xyz)]=api.digest(xyz)
    atomic_write_json(ready_path,ready)
    inputs[str(ready_path)]=api.digest(ready_path);inputs[str(xyz)]=api.digest(xyz)
    for stage in ('energy','thermal'):
        for name in ('tasks','done','claims','failures','workers','pilot','launches'):(out/stage/name).mkdir(parents=True)
    atomic_write_json(out/'plan.json',plan)


def run(campaign):
    verify_gate_owner()
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):raise ValueError('Two disjoint CPUs and matching caps required')
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    host_capacity(2,minimum_ram_gib=plan['maximum_RSS_GiB']+32)
    if gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))['memory_MiB']<32768:raise ValueError('A 32 GiB GPU is required')
    guard=json.loads((campaign/'host_guard/status.json').read_text())
    if guard['pid']!=os.getpid() or guard['maximum_RSS_GiB']!=plan['maximum_RSS_GiB'] or guard['minimum_free_GiB']!=32:
        raise ValueError('The exact owned host-resource guard is required')
    task=plan['tasks'][0];key=task['key'];claim=campaign/'energy/claims'/key;claim.mkdir()
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    atomic_write_json(claim/'owner.json',owner)
    directory=campaign/'energy/tasks'/key;directory.mkdir(exist_ok=False)
    state=Path(plan['raw_root'])/'energy'/key/'rest';state.mkdir(parents=True,exist_ok=False)
    source=Path(plan['source_geometry_report']['state_dir'])
    for name in ('opt.xyz','result.json'):shutil.copyfile(source/name,state/name)
    shutil.copyfile(source/'result.json',state/'geometry_parent.json')
    record=dict(status='computing',owner=owner,task=task,state_dir=str(state),started=time.time(),
        geometry_parent=str(source),GPU_resource_qualification=plan['GPU_resource_qualification'],native_G_computed=False)
    atomic_write_json(directory/'status.json',record)
    def stop(signum,frame):raise KeyboardInterrupt('Owned GPU energy stop; checkpoints preserved')
    previous=signal.signal(signal.SIGTERM,stop)
    try:
        from redox.core import protocol
        from redox.qm import sp
        from redox.workflows.split_electronic import composite_result
        memory=api.load_module('split_GPU_energy_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        adapter=api.load_module('split_GPU_energy_workspace',GC)
        with memory.MemoryMonitor() as used,adapter.compact_smd_execution(state/'execution',plan['execution_workspace_GiB']), (directory/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
            sp.ensure(key,'rest',state,task['charge'],task['mult'],backend='gpu')
            active=protocol.load_record(state,task['charge'],task['mult'])
            if active is None:raise ValueError('The canonical GPU energy record failed normal admission')
            combined=composite_result(json.loads((state/'geometry_parent.json').read_text()),active,
                protocol.geom_sha1(state/'opt.xyz'),protocol.record_path(state,task['charge'],task['mult']),backend='gpu')
            atomic_write_json(state/'result.json',combined)
        if not 0<used.peak<=plan['maximum_peak_GPU_MiB']:raise ValueError('GPU active-energy resource bound failed')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        record.update(status='state_records_verified',E_smd_eV=active['e_smd_eV'],observed_peak_GPU_MiB=used.peak,
            physical_QC='passed',physical_QC_source='Exact-byte inherited stationary geometry and independent total-energy derivatives',
            topology_QC='ok',gradient_max_Ha_Bohr=plan['source_geometry_report']['gradient_max_Ha_Bohr'],
            conservative_total_force_check=plan['source_geometry_report']['conservative_total_force_check'],
            output_sha256={str(p):api.digest(p) for p in state.rglob('*') if p.is_file() and p.suffix in ('.json','.xyz','.npy')},finished=time.time())
        atomic_write_json(directory/'status.json',record)
        receipt=campaign/'energy/done'/(key+'.json')
        atomic_write_json(receipt,dict(key=key,status_path=str(directory/'status.json'),status_sha256=api.digest(directory/'status.json'),owner=owner))
        from redox.workflows.checkpoints import verify_state_checkpoint
        verify_state_checkpoint(receipt)
        atomic_write_json(campaign/'workspace_energy_status.json',dict(status='stationary_electronic_case_checked_native_G_pending',
            source_receipt=str(receipt),source_receipt_sha256=api.digest(receipt),native_G_computed=False))
    except BaseException as error:
        record.update(status='failed_preserved',error=repr(error),finished=time.time());atomic_write_json(directory/'status.json',record);raise
    finally:
        signal.signal(signal.SIGTERM,previous);shutil.rmtree(claim)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run'))
    for name in ('geometry','pilot','review','controls','campaign'):p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.campaign is None:p.error('--campaign required')
    if a.action=='prepare':
        if any(getattr(a,k) is None for k in ('geometry','pilot','review','controls')):p.error('prepare requires geometry, pilot, review and controls')
        prepare(a.geometry.resolve(),a.pilot.resolve(),a.review.resolve(),a.controls.resolve(),a.campaign.resolve())
    else:run(a.campaign.resolve())


if __name__=='__main__':main()
