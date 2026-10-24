#!/usr/bin/env python
"""Replay-gated native polymer G with streamed PCM and component checkpoints."""
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import signal
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability

ADAPTER=ROOT/'scripts/pipeline/streamed_pcm_execution.py'
CHECKER=ROOT/'scripts/validation/qualify_streamed_pcm.py'


def prepare(source, controls, response, out):
    checker=api.load_module('streamed_native_full_replay_check',CHECKER)
    evidence=checker.checked_controls(controls)
    original=sys.argv
    try:
        sys.argv=['switch_curvature_polymer_campaign.py','prepare-native','--source',str(source),
            '--out',str(out),'--controls',str(ROOT/'results/validation/pcm_switch_curvature_guarded_2026-10-23'),
            '--response-proof',str(response),'--workspace-GiB','28',
            '--minimum-GPU-MiB','32768','--maximum-native-GPU-MiB','31744']
        module=api.load_module('existing_exact_native_handoff',ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py')
        module.main()
    finally:sys.argv=original
    plan=api.read_plan(out)
    if len(plan['tasks'])!=1 or (plan['tasks'][0]['charge'],plan['tasks'][0]['mult'])!=(0,1):
        raise ValueError('This polymer native driver currently covers one neutral singlet rest state')
    plan.update(execution_driver=str(Path(__file__).resolve()),
        execution_profile='streamed_FP64_PCM_curvature',streamed_controls=str(controls),
        streamed_control_evidence=evidence,host_required_bytes=128*2**30,
        maximum_host_RSS_GiB=80,externally_managed_stages=True,
        native_host_budget_GiB=128,
        purpose='Independent large streamed-native resource/minimum/G qualification; exact electronic parent retained')
    plan['producer_sha256'].update(evidence)
    for p in (Path(__file__).resolve(),ADAPTER,CHECKER,ROOT/'scripts/fleet/checkpoint_scf_worker.py'):
        plan['producer_sha256'][str(p)]=api.digest(p)
    atomic_write_json(out/'plan.json',plan)
    print(json.dumps(dict(status='streamed_native_case_prepared',host_RSS_limit_GiB=80,
        native_G_qualified=False,source=str(source))))


@contextmanager
def checkpoint_components(directory):
    import numpy as np
    from gpu4pyscf.solvent import smd
    from gpu4pyscf.solvent.hessian import smd as hsmd
    directory.mkdir(parents=True,exist_ok=False)
    previous_hess=smd.SMD.hess
    previous_cds=hsmd.get_cds
    def save(kind,solvent,value):
        array=np.asarray(value.get() if hasattr(value,'get') else value)
        np.save(directory/(kind+'.npy'),array)
        atomic_write_json(directory/(kind+'.json'),dict(charge=solvent.mol.charge,
            mult=solvent.mol.spin+1,coordinates_Bohr=solvent.mol.atom_coords().tolist(),
            phase=kind,qualified_native_G=False,updated=time.time()))
    def solvent_hess(self,*args,**kwargs):
        value=previous_hess(self,*args,**kwargs)
        save('explicit_PCM',self,value)
        return value
    def cds(self,*args,**kwargs):
        value=previous_cds(self,*args,**kwargs)
        save('CDS',self,value)
        return value
    smd.SMD.hess,hsmd.get_cds=solvent_hess,cds
    try:yield
    finally:smd.SMD.hess,hsmd.get_cds=previous_hess,previous_cds


def run(campaign):
    verify_gate_owner()
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):
        raise ValueError('Two disjoint CPUs and matching native caps required')
    plan=api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    host_capacity(2,minimum_ram_gib=128)
    if gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))['memory_MiB']<32768:
        raise ValueError('The full native target still requires a 32-GiB-class GPU')
    checker=api.load_module('large_streamed_native_controls',CHECKER)
    witnesses=checker.checked_controls(Path(plan['streamed_controls']))
    energy,invalid=api.admitted_records(campaign,'energy')
    if invalid or len(energy)!=1:raise ValueError('The exact stationary electronic parent is required')
    key,parent=next(iter(energy.items()))
    if (campaign/'thermal/tasks'/key).exists():
        raise ValueError('Preserve an existing attempt; prepare fresh recovery storage')
    claim=campaign/'thermal/claims'/key;claim.mkdir()
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    atomic_write_json(claim/'owner.json',owner)
    runtime=campaign/'streamed_native';runtime.mkdir(exist_ok=False)
    report=dict(status='streamed_native_resource_and_physical_pilot_running',owner=owner,started=time.time(),
        source_energy_state=parent['state_dir'],native_G_qualified=False,
        control_evidence=witnesses,host_RSS_limit_GiB=plan['maximum_host_RSS_GiB'])
    atomic_write_json(runtime/'status.json',report)
    def stop(signum,frame):raise KeyboardInterrupt('Owned streamed-native stop; components preserved')
    previous_handler=signal.signal(signal.SIGTERM,stop)
    import redox.qm.projected_thermal as projected
    previous_projection=projected.ensure_projected_thermal
    execution=Path(plan['raw_root'])/'streamed_execution'
    def seal(path,geometry):
        from redox.qm.hessian_resume import annotate_finer_hessian
        record=json.loads(Path(path).read_text())
        if 'execution_derivatives' not in record:
            annotate_finer_hessian(path,geometry)
        record=json.loads(Path(path).read_text())
        record['execution_streamed_PCM']=dict(mode=plan['execution_profile'],
            adapter_sha256=api.digest(ADAPTER),full_spin_replay_evidence=witnesses,
            source_geometry_sha256=api.digest(geometry),scientific_model_changed=False,
            precision='float64',surface_grid_changed=False)
        atomic_write_json(path,record)
        return previous_projection(path,geometry)
    projected.ensure_projected_thermal=seal
    try:
        from redox.qm.cds_finite import finer_cds_hessian
        memory_module=api.load_module('streamed_native_GPU_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        adapter=api.load_module('large_streamed_PCM_execution',ADAPTER)
        with memory_module.MemoryMonitor() as memory,finer_cds_hessian(1e-5),adapter.compact_smd_execution(execution,28),checkpoint_components(execution/'components'):
            api.calculate_native_G(campaign,parent)
        stats=json.loads((execution/'streamed_PCM.json').read_text())
        if (not 0<memory.peak<=31744 or stats['status']!='streamed_PCM_finished_requires_full_replay_and_physical_QC'
                or min(stats[k] for k in ('switch_contractions','diagonal_contractions','offdiagonal_primitive_calls'))<1
                or stats['peak_GPU_block_cache_bytes']>128*1024**2):
            raise ValueError('The full native resource/operator qualification failed')
        for name in ('explicit_PCM','CDS'):
            if not (execution/'components'/(name+'.npy')).is_file():
                raise ValueError('An additive Hessian component checkpoint is missing')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        path=campaign/'thermal/tasks'/key/'status.json'
        state=json.loads(path.read_text())
        state.update(observed_peak_GPU_MiB=memory.peak,execution_profile=plan['execution_profile'],
            full_spin_replay_evidence=witnesses)
        state['output_sha256'].update({str(p):api.digest(p) for p in execution.rglob('*') if p.is_file() and p.suffix in ('.json','.npy','.py')})
        atomic_write_json(path,state)
        from redox.workflows.checkpoints import verify_state_checkpoint
        receipt=campaign/'thermal/done'/(key+'.json')
        atomic_write_json(receipt,dict(key=key,owner=owner,status_path=str(path),status_sha256=api.digest(path),finished=time.time()))
        verify_state_checkpoint(receipt)
        report.update(status='streamed_native_case_passed_independent_review_pending',
            observed_peak_GPU_MiB=memory.peak,native_G_qualified=True,G_eV=state['G_eV'])
    except BaseException as error:
        report.update(status='streamed_native_pilot_failed_preserved',error=repr(error))
        raise
    finally:
        projected.ensure_projected_thermal=previous_projection
        signal.signal(signal.SIGTERM,previous_handler)
        report['finished']=time.time();atomic_write_json(runtime/'status.json',report)
        shutil.rmtree(claim)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run'))
    for name in ('source','controls','response','campaign'):parser.add_argument('--'+name,type=Path)
    args=parser.parse_args()
    if args.campaign is None:parser.error('--campaign required')
    if args.action=='prepare':
        if any(getattr(args,k) is None for k in ('source','controls','response')):
            parser.error('prepare requires --source, --controls and --response')
        prepare(args.source.resolve(),args.controls.resolve(),args.response.resolve(),args.campaign.resolve())
    else:run(args.campaign.resolve())


if __name__=='__main__':main()
