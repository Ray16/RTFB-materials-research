#!/usr/bin/env python
"""Reserved gas counterpoise controls; no solvent or thermochemical correction."""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

CORE=ROOT/'scripts/validation/local_counterpoise.py'
PROFILE=ROOT/'environments/cutensor22/runtime.json'


def prepare_control(out):
    core=api.load_module('counterpoise_declared_protocol',CORE)
    out.mkdir(parents=True,exist_ok=False);case=out/'water_dimer';case.mkdir()
    coords=[[0,0,0],[0,.76,.58],[0,-.76,.58],[0,0,2.9],[0,.76,3.48],[0,-.76,3.48]]
    geometry=dict(symbols=['O','H','H']*2,coordinates_A=coords,A_indices=[0,1,2])
    atomic_write_json(case/'geometry.json',geometry)
    files=[Path(__file__).resolve(),CORE,ROOT/'scripts/pipeline/local_motif_energy.py',
        ROOT/'scripts/pipeline/run_pi_stacking.py',ROOT/'scripts/fleet/exclusive_python_worker.py',
        ROOT/'scripts/fleet/tensor_runtime_worker.py',ROOT/'scripts/fleet/bounded_host_python_worker.py']
    plan=dict(task=dict(key='water_dimer',atoms=6,charge=0,mult=1,role='gas_counterpoise_CPU_GPU_control'),
        settings=core.SETTINGS,input_geometry=str(case/'geometry.json'),
        input_sha256={str(p):api.digest(p) for p in (case/'geometry.json',PROFILE)},
        producer_sha256={str(p):api.digest(p) for p in files},runtime_profile=str(PROFILE),
        execution_driver=str(Path(__file__).resolve()),maximum_host_RSS_GiB=32,
        raw_root=str(ROOT/'calcs/validation'/out.name/'water_dimer'),
        primary_candidate_admission=False,solvent_counterpoise_computed=False,maximum_GPU_MiB=14336)
    atomic_write_json(case/'plan.json',plan)
    atomic_write_json(out/'campaign.json',dict(status='gas_counterpoise_control_prepared',entries=[dict(key='water_dimer',
        campaign=str(case),plan_sha256=api.digest(case/'plan.json'))],CPU_GPU_comparison_required=True,
        source_sha256={str(p):api.digest(p) for p in files},primary_candidate_admission=False))


def run(case):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();plan=api.read_plan(case);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or set(os.environ.get(k) for k in caps)!={'2'}:raise ValueError('Two assigned CPUs and capped native threads required')
    geometry=json.loads(Path(plan['input_geometry']).read_text());raw=Path(plan['raw_root']);raw.mkdir(parents=True,exist_ok=False)
    stage=case/'counterpoise';stage.mkdir(exist_ok=False);path=stage/'status.json'
    report=dict(status='gas_counterpoise_control_running',started=time.time(),raw_root=str(raw),
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        solvent_counterpoise_computed=False,primary_candidate_admission=False)
    atomic_write_json(path,report)
    try:
        core=api.load_module('counterpoise_actual_five_term_evaluation',CORE)
        if plan['settings']!=core.SETTINGS:raise ValueError('Counterpoise settings changed')
        monitor=api.load_module('counterpoise_actual_GPU_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        with monitor.MemoryMonitor() as observed,(raw/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
            result=core.compute(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices'],'gpu',raw/'GPU')
        if not 0<observed.peak<=plan['maximum_GPU_MiB']:raise ValueError('Counterpoise resource ceiling exceeded')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        report.update(status='gas_counterpoise_GPU_complete_CPU_comparison_and_review_pending',observed_peak_GPU_MiB=observed.peak,
            sums=result['sums'],output_sha256={str(p):api.digest(p) for p in raw.rglob('*') if p.is_file() and p.suffix in ('.json','.npy','.h5')})
    except BaseException as error:report.update(status='gas_counterpoise_control_failed_preserved',error=repr(error));raise
    finally:report['finished']=time.time();atomic_write_json(path,report)


def launch(case,gpu,preflight):
    # Reuse the already checked one-process reservation/host/exclusivity path.
    # Its __file__ identifies only the child program; numerical code stays here.
    launcher=api.load_module('counterpoise_reserved_resource_placement',ROOT/'scripts/pipeline/local_motif_energy.py')
    launcher.__file__=str(Path(__file__).resolve())
    return launcher.launch(case,gpu,preflight)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run','launch'));p.add_argument('--campaign',type=Path,required=True);p.add_argument('--gpu',type=int);p.add_argument('--preflight',type=Path);a=p.parse_args()
    if a.action=='prepare':
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
        prepare_control(a.campaign.resolve())
    elif a.action=='run':run(a.campaign.resolve())
    else:print(json.dumps(launch(a.campaign.resolve(),a.gpu,a.preflight.resolve())))


if __name__=='__main__':main()
