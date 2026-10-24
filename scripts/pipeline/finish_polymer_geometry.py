#!/usr/bin/env python
"""Finish force validation after a converged optimizer's resource interruption."""
import argparse
import json
import math
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
CONTROL=ROOT/'scripts/validation/qualify_gc_workspace.py'


def validate_predecessor(status,result):
    if (status.get('status')!='geometry_stage_failed_preserved'
            or 'out of memory' not in status.get('error','').lower()):
        raise ValueError('Only the preserved resource-interrupted geometry stage is eligible')
    if any(result.get(k) is not True for k in ('optimized','geometry_optimization_converged','optimization_SCF_converged','geometry_only_result')):
        raise ValueError('The completed optimizer and SCF flags must be explicit')
    if (result.get('charge'),result.get('mult'),result.get('solvent'))!=(0,1,'acetonitrile'):
        raise ValueError('This force-recovery profile covers a neutral-singlet MeCN geometry')


def prepare(source,controls,out):
    old=api.read_plan(source);api.verify_sources(old['producer_sha256']);api.verify_sources(old['input_sha256'])
    report_path=source/'geometry/status.json';report=json.loads(report_path.read_text())
    state=Path(report['state_dir']);result_path=state/'result.json';result=json.loads(result_path.read_text())
    validate_predecessor(report,result)
    review_path=source/'geometry/geometry_review.json'
    if json.loads(review_path.read_text())['issues']:raise ValueError('The returned geometry failed structural checks')
    if api.load_module('force_recovery_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py').owner_alive(report['owner']):
        raise ValueError('The predecessor must finish before recovery')
    witness=api.load_module('force_recovery_GC_controls',CONTROL).checked_controls(controls)
    if len(old['tasks'])!=1:raise ValueError('One exact returned geometry is required')
    out.mkdir(parents=True,exist_ok=False);shutil.copytree(source/'controls',out/'controls')
    key=old['tasks'][0]['key'];target=out/'live_inputs'/key;target.mkdir(parents=True)
    shutil.copyfile(state/'opt.xyz',target/'input.xyz')
    ready=json.loads((source/'live_inputs'/key/'task.json').read_text())
    ready['input_xyz']=str(target/'input.xyz')
    ready['input_sha256'][ready['input_xyz']]=api.digest(target/'input.xyz')
    ready['source_optimizer_result']=str(result_path)
    atomic_write_json(target/'task.json',ready)
    inputs=dict(old['input_sha256'],**witness)
    for path in (source/'plan.json',report_path,result_path,state/'opt.xyz',review_path,target/'input.xyz',target/'task.json'):
        inputs[str(path)]=api.digest(path)
    sources=dict(old['producer_sha256'],**witness)
    for path in (Path(__file__).resolve(),GC,CONTROL):sources[str(path)]=api.digest(path)
    for stage in ('energy','thermal'):
        for name in ('tasks','done','claims','failures','workers','pilot','launches'):(out/stage/name).mkdir(parents=True)
    policy=dict(old['geometry_resource_policy'],minimum_host_RAM_GiB=64,maximum_RSS_GiB=64)
    plan=dict(old,input_sha256=inputs,producer_sha256=sources,
        raw_root=str(ROOT/'calcs/polymer_models'/out.parent.name/out.name),
        geometry_resource_policy=policy,GC_controls=str(controls),source_optimizer_result=str(result_path),
        source_geometry_report=str(report_path),source_optimization_reused=True,
        purpose='Independent final force/energy differences on exact converged optimizer coordinates; no repeated optimization or active SP',
        created=time.time(),externally_managed_stages=True)
    atomic_write_json(out/'plan.json',plan)


def run(campaign):
    verify_gate_owner()
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):raise ValueError('Two disjoint CPUs and matching caps required')
    host_capacity(2,minimum_ram_gib=96)
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))['memory_MiB']<32768:raise ValueError('32 GiB GPU required')
    task=plan['tasks'][0];key=task['key'];claim=campaign/'energy/claims'/key;claim.mkdir()
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    atomic_write_json(claim/'owner.json',owner)
    directory=campaign/'geometry';directory.mkdir(exist_ok=False)
    raw=Path(plan['raw_root']);state=raw/'energy'/key/'rest';state.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(campaign/'live_inputs'/key/'input.xyz',state/'opt.xyz')
    result=json.loads(Path(plan['source_optimizer_result']).read_text())
    result.update(opt_out=str(state/'opt.xyz'),source_optimizer_result=plan['source_optimizer_result'])
    atomic_write_json(state/'result.json',result)
    record=dict(status='returned_geometry_force_checks_running',owner=owner,task=task,state_dir=str(state),
        started=time.time(),source_optimization_reused=True,optimization_recomputed=False,
        active_electronic_energy_computed=False,native_Hessian_qualified=False,G_computed=False)
    atomic_write_json(directory/'status.json',record)
    def stop(signum,frame):raise KeyboardInterrupt('Owned final-force recovery stop; checkpoints preserved')
    previous=signal.signal(signal.SIGTERM,stop)
    try:
        engine=api.load_module('returned_geometry_checkpoint_engine',ROOT/'scripts/pipeline/optimize_polymer_geometry.py')
        forces=api.load_module('unchanged_final_energy_checks',ROOT/'scripts/pipeline/resume_polymer_forces.py')
        gc=api.load_module('returned_geometry_GC',GC)
        audit=api.load_module('returned_geometry_audit',ROOT/'scripts/validation/audit_polymer_geometries.py')
        memory=api.load_module('returned_geometry_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        statistics=dict(collections=0,unreachable_objects_collected=0)
        with memory.MemoryMonitor() as used,engine.geometry_engine(raw),gc.collect_before_build(statistics):
            check=forces.final_energy_check(state/'opt.xyz',state/'total_energy_force_check')
        baseline=json.loads((state/'total_energy_force_check/baseline.json').read_text())
        maximum=max(abs(v) for row in baseline['gradient_Ha_Bohr'] for v in row)
        if not math.isfinite(maximum) or maximum>4.5e-4 or not check['passed']:
            raise ValueError('Original stationarity/finite-difference thresholds failed')
        if not 0<used.peak<=plan['geometry_resource_policy']['maximum_peak_GPU_MiB']:
            raise ValueError('Final-force GPU resource bound failed')
        review=audit.audit(task,state/'opt.xyz')
        if review['issues']:raise ValueError('Returned geometry identity/topology failed')
        atomic_write_json(directory/'geometry_review.json',review)
        atomic_write_json(directory/'GC_execution.json',statistics)
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        record.update(status='stationary_geometry_checked_active_parent_and_native_minimum_pending',
            gradient_max_Ha_Bohr=maximum,conservative_total_force_check=check,observed_peak_GPU_MiB=used.peak,
            output_sha256={str(p):api.digest(p) for p in [state/'opt.xyz',state/'result.json',directory/'geometry_review.json',
                directory/'GC_execution.json',*list((state/'total_energy_force_check').glob('*.json'))]})
    except BaseException as error:
        record.update(status='returned_geometry_force_recovery_failed_preserved',error=repr(error));raise
    finally:
        signal.signal(signal.SIGTERM,previous);record['finished']=time.time();atomic_write_json(directory/'status.json',record);shutil.rmtree(claim)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run'))
    for name in ('source','controls','campaign'):p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.campaign is None:p.error('--campaign required')
    if a.action=='prepare':
        if a.source is None or a.controls is None:p.error('prepare requires source and controls')
        prepare(a.source.resolve(),a.controls.resolve(),a.campaign.resolve())
    else:run(a.campaign.resolve())


if __name__=='__main__':main()
