#!/usr/bin/env python3
"""Reconstruct a closed three-SCF gradient diagnostic without admitting native G."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows.molecular_native_namespace import verify

PRODUCER_BYTES=Path(__file__).read_bytes()

def analyze(namespace,output,*,recover_post_scf_context=False):
    if output.exists():raise FileExistsError('Preserve earlier diagnostic analyses')
    spec=verify(namespace);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text());stage=json.loads((case/'status.json').read_text())
    guard=json.loads((case/'residency_guard/owner.json').read_text());residency=json.loads((case/'residency_guard/state.json').read_text())
    host=json.loads((case/'host_guard/status.json').read_text());runtime=json.loads((case/'tensor_runtime/status.json').read_text())
    post_error="ValueError('No response solve exercised')"
    recovered=(recover_post_scf_context and stage['status']=='polymer_gradient_diagnostic_failed_preserved'
        and stage.get('error')==post_error and host.get('error')==post_error and runtime.get('error')==post_error)
    if recovered:
        response=json.loads((case/'execution/response_solver.json').read_text())
        if (response.get('calls')!=[] or response.get('error')!=post_error
                or 'graceful_stop_requested' in host or host['free_RAM_bytes']<32*2**30):
            raise ValueError('Only the closed unused-response-context exception is recoverable')
    if ((not recovered and stage['status']!='polymer_gradient_components_written_independent_analysis_pending') or len(stage['SCF_points'])!=3
            or plan['role']!='polymer_gradient_diagnostic' or not spec.get('scientific_sources_read_only')
            or not guard.get('guard_stopped') or residency['foreign_resident_PIDs']
            or not residency['own_residency_observed'] or (case/'residency_guard/interruption.json').exists()
            or not stage['owner']['pid']==guard['pid']==runtime['pid']==host['pid']
            or not stage['owner']['process_start']==guard['process_start']==host['process_start']
            or host['status']!=('effective_memory_guard_failed_preserved' if recovered else 'effective_memory_guard_finished_profile_review_pending')
            or host['maximum_RSS_GiB']!=80 or host['minimum_free_GiB']!=32
            or not 0<max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])<=80*2**30
            or not host['CPU_affinity']==guard['CPU_affinity']==runtime['CPU_affinity'] or len(host['CPU_affinity'])!=2
            or runtime['status']!=('runtime_or_task_failed_preserved' if recovered else 'task_returned_under_checked_cutensor_runtime')
            or runtime['program']!=str(Path(spec['sources'])/'worker.py')
            or (not recovered and not 0<stage['observed_peak_GPU_MiB']<=31744)):
        raise ValueError('Closed exact diagnostic/resource chain required')
    gradients=[];parts=[];records=[];used=[]
    for index,point in enumerate(plan['points']):
        folder=case/f'point{index}';original=json.loads(Path(point['evaluation']).read_text());row=json.loads((folder/'result.json').read_text())
        gradient=np.load(folder/'gradient.npy',allow_pickle=False)
        component={name:np.load(folder/(name+'.npy'),allow_pickle=False) for name in ('solute','solvent','cds')}
        if (gradient.shape!=(len(plan['symbols']),3) or not np.isfinite(gradient).all()
                or any(v.shape!=gradient.shape or not np.isfinite(v).all() for v in component.values())
                or np.max(abs(gradient-sum(component.values())))>1e-10 or row['SCF_converged'] is not True):
            raise ValueError('Converged component gradients must sum to the full force')
        checkpoints=sorted((case/'SCF_evaluations'/f'{index:04d}').glob('*/SCF_result.json'))
        if len(checkpoints)!=1:raise ValueError('One observed SCF per diagnostic point required')
        solved=json.loads(checkpoints[0].read_text());settings=json.loads(checkpoints[0].with_name('settings.json').read_text())
        if (solved['converged'] is not True or abs(solved['energy_Ha']-row['E_Ha'])>1e-10
                or any(settings[k]!=v for k,v in dict(charge=0,mult=1,basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,conv_tol_grad=1e-8,grid_level=3).items())
                or np.max(abs(np.asarray(settings['coordinates_Bohr'])-np.asarray(original['coordinates_Bohr'])))>1e-10):
            raise ValueError('Actual diagnostic SCF state/settings/geometry differ from the source point')
        difference=float(np.max(abs(gradient-np.asarray(original['gradient_Ha_Bohr']))))
        if abs(difference-row['maximum_gradient_difference_Ha_Bohr'])>1e-12 or abs(row['energy_difference_Ha']-(row['E_Ha']-original['E_Ha']))>1e-12:
            raise ValueError('Reported repeated energy/force differences do not reconstruct')
        records.append(dict(row,maximum_gradient_difference_Ha_Bohr=difference));gradients.append(gradient);parts.append(component)
        used.extend([folder/'result.json',folder/'gradient.npy',*[folder/(name+'.npy') for name in component],checkpoints[0],checkpoints[0].with_name('settings.json'),Path(point['evaluation'])])
    original_spec=json.loads(Path(plan['original_namespace']).read_text());reader=NamespaceSources(original_spec)
    hessian_path=reader.path(original_spec['original_raw_root'])/'hessian.npy';hessian=np.load(hessian_path,allow_pickle=False)
    selected=plan['failed_column'];column=3*selected['atom']+selected['axis'];step=selected['step_Bohr']
    finite=(gradients[2]-gradients[1]).ravel()/(2*step);error=float(np.linalg.norm(finite-hessian[:,column])/np.linalg.norm(finite))
    component_norms={name:float(np.linalg.norm((parts[2][name]-parts[1][name]).ravel()/(2*step))) for name in parts[0]}
    output.parent.mkdir(parents=True,exist_ok=True);snapshot=output.with_name(output.stem+'_producer.py');snapshot.write_bytes(PRODUCER_BYTES)
    used.extend([namespace,case/'plan.json',case/'status.json',case/'residency_guard/owner.json',case/'residency_guard/state.json',case/'host_guard/status.json',case/'tensor_runtime/status.json',hessian_path,snapshot])
    if recovered:used.append(case/'execution/response_solver.json')
    sources={str(p):sha256_file(p) for p in used}
    result=dict(status='three_SCF_gradient_diagnostic_reconstructed_no_native_admission',points=records,
        diagnostic_namespace=str(namespace),original_namespace=plan['original_namespace'],
        original_failed_relative_error=selected['relative_error'],repeated_relative_error_against_source_Hessian=error,
        repeated_column_meets_one_percent=error<=.01,component_finite_difference_column_norms=component_norms,
        component_column_norms_are_not_component_errors=True,full_native_G_qualified=False,source_sha256=sources,
        recovered_after_unused_response_context_error=recovered,full_execution_qualified=not recovered,
        GPU_peak_qualified=not recovered,
        interpretation='Repeated full gradients and their components only. A source Hessian comparison does not supply a fully qualified native state or identify the causal component without further controls.')
    atomic_write_json(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--wait-hours',type=float,default=0)
    p.add_argument('--recover-post-scf-context',action='store_true');a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if not 0<=a.wait_hours<=48:p.error('Use a finite 0–48 hour observation window')
    spec=verify(a.spec);case=Path(spec['case']);deadline=time.time()+a.wait_hours*3600
    helper=ROOT/'scripts/pipeline/observe_archived_native_namespace.py';definition=importlib.util.spec_from_file_location('diagnostic_owner_check',helper);owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    while True:
        state=case/'status.json';ready=False
        if state.exists():
            status=json.loads(state.read_text())
            if 'finished' in status:
                if status['status']=='polymer_gradient_diagnostic_failed_preserved' and not a.recover_post_scf_context:raise RuntimeError(status['error'])
                try:ready=owners.alive(status['owner']) is False and json.loads((case/'residency_guard/owner.json').read_text()).get('guard_stopped') is True
                except (OSError,subprocess.SubprocessError):ready=False
        if ready:break
        if time.time()>=deadline:raise TimeoutError('Diagnostic still pending; existing worker and checkpoints are untouched')
        time.sleep(10)
    result=analyze(a.spec.resolve(),a.out.resolve(),recover_post_scf_context=a.recover_post_scf_context);print(json.dumps({k:result[k] for k in ['status','original_failed_relative_error','repeated_relative_error_against_source_Hessian']}))
