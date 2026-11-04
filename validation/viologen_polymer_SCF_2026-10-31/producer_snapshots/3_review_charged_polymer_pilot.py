#!/usr/bin/env python3
"""Check charged initial-state SCF/forces and compare execution profiles."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path

import h5py
import numpy as np
from ase.io import read as read_atoms
from pyscf.lib import param

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);used=[namespace]
    def read(path):used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json');stage=read(case/'status.json');placement=case.parent/'placement'
    completion=read(placement/'completion.json');launch=read(placement/'launch.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    owner=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    definition=importlib.util.spec_from_file_location('ionic_pilot_closed_owner',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    if (stage['status']!='charged_polymer_SCF_gradient_resource_pilot_finished_review_pending'
            or stage['native_G_qualified'] or stage['material_admission'] or owners.alive(stage['owner']) is not False
            or completion['returncode']!=0 or completion['budget_stop'] is not None
            or not owner['guard_stopped'] or not residency['own_residency_observed'] or residency['foreign_resident_PIDs']
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or host['maximum_RSS_GiB']!=40 or host['free_RAM_bytes']<32*2**30
            or max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])>40*2**30
            or not 0<stage['observed_peak_GPU_MiB']<=plan['maximum_peak_GPU_MiB']
            or launch['device']['memory_MiB']<plan['minimum_GPU_memory_MiB']
            or not stage['owner']['pid']==owner['pid']==runtime['pid']==host['pid']
            or not stage['owner']['process_start']==owner['process_start']==host['process_start']
            or not owner['CPU_affinity']==runtime['CPU_affinity']==host['CPU_affinity']==launch['CPU_affinity']):
        raise ValueError('Closed charged-state resource and ownership checks required')
    route=read(Path(spec['sources'])/'prospective_route.json')
    if not route['proceed_to_full_redox'] or not route['exact_route_binding_verified']:
        raise ValueError('The exact graft route must pass prospectively')
    xyz=Path(spec['sources'])/'input.xyz';atoms=read_atoms(xyz);used.append(xyz)
    result=read(case/'evaluation.json');gradient_path=case/'gradient.npy';gradient=np.load(gradient_path,allow_pickle=False);used.append(gradient_path)
    settings_path=case/'SCF_evaluations/0000/0000/settings.json';settings=read(settings_path);solved=read(settings_path.with_name('SCF_result.json'));checkpoint=settings_path.with_name('SCF.h5');used.append(checkpoint)
    if (result['charge']!=plan['charge'] or result['mult']!=plan['mult'] or not result['SCF_converged']
            or abs(result['spin'].get('spin_contam',math.inf))>.1 or result['geometry_optimized']
            or result['native_G_qualified'] or gradient.shape!=(plan['atoms'],3) or not np.isfinite(gradient).all()
            or not np.allclose(result['coordinates_Bohr'],atoms.positions/param.BOHR,rtol=0,atol=1e-10)
            or not np.array_equal(gradient,result['gradient_Ha_Bohr'])
            or np.max(abs(gradient-sum(np.asarray(v) for v in result['gradient_components_Ha_Bohr'].values())))>1e-10
            or any(settings[k]!=v for k,v in dict(charge=plan['charge'],mult=plan['mult'],basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,conv_tol_grad=1e-8,max_cycle=200,grid_level=3).items())
            or not solved['converged'] or abs(solved['energy_Ha']-result['E_Ha'])>1e-10):
        raise ValueError('Actual fixed charged geometry, SCF settings or force components differ')
    with h5py.File(checkpoint,'r') as handle:
        if abs(float(handle['scf/e_tot'][()])-result['E_Ha'])>1e-10:raise ValueError('SCF checkpoint energy differs')
    return dict(status='initial_charged_polymer_SCF_gradient_and_resources_independently_checked',
        namespace=str(namespace),candidate_alias=plan['exact_candidate_alias'],graph_key=plan['prepared_case']['key'],
        geometry_sha256=sha256_file(xyz),charge=plan['charge'],mult=plan['mult'],atoms=plan['atoms'],
        execution_profile=plan['execution_profile'],E_Ha=result['E_Ha'],gradient_Ha_Bohr=gradient.tolist(),
        initial_gradient_max_Ha_Bohr=float(abs(gradient).max()),peak_GPU_MiB=stage['observed_peak_GPU_MiB'],
        peak_host_RSS_bytes=max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes']),
        solvent_model=plan['solvent_model'],calculation_ion_convention=plan['calculation_ion_convention'],
        solvent_optimization_complete=False,conformer_sampling_qualified=False,interaction_strength_computed=False,
        native_G_qualified=False,material_admission=False,source_sha256={str(p):sha256_file(p) for p in used},
        producer_sha256=sha256_file(__file__))


def compare(paths):
    rows=[json.loads(p.read_text()) for p in paths]
    if len(rows)!=2 or {r['execution_profile'] for r in rows}!={'compact_28GiB','direct_large_GPU_reference'}:
        raise ValueError('One direct and one compact independently checked pilot required')
    for row in rows:
        if row['status']!='initial_charged_polymer_SCF_gradient_and_resources_independently_checked' or row['native_G_qualified']:
            raise ValueError('Only checked initial-state pilots can compare')
        if any(sha256_file(p)!=h for p,h in row['source_sha256'].items()):raise ValueError('Pilot backing changed')
    fields=('candidate_alias','graph_key','geometry_sha256','charge','mult','solvent_model','calculation_ion_convention')
    if any(rows[0][k]!=rows[1][k] for k in fields):raise ValueError('Different chemical or solution scenarios cannot compare')
    delta=rows[1]['E_Ha']-rows[0]['E_Ha'];force=float(abs(np.asarray(rows[1]['gradient_Ha_Bohr'])-rows[0]['gradient_Ha_Bohr']).max())
    return dict(status='charged_polymer_execution_profiles_compared',energy_difference_Ha=delta,
        maximum_gradient_difference_Ha_Bohr=force,energy_limit_Ha=1e-7,gradient_limit_Ha_Bohr=1e-6,
        initial_state_execution_agreement_passed=abs(delta)<=1e-7 and force<=1e-6,
        native_G_qualified=False,polymer_occupancy_result_qualified=False,
        source_sha256={str(p):sha256_file(p) for p in paths},producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('review','compare'))
    parser.add_argument('--spec',type=Path);parser.add_argument('--review',type=Path,action='append');parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.out.exists():raise FileExistsError('Preserve earlier charged-pilot reviews')
    result=review(args.spec.resolve()) if args.action=='review' else compare([p.resolve() for p in args.review])
    atomic_write_json(args.out,result);print(json.dumps({k:v for k,v in result.items() if k in ['status','charge','atoms','initial_gradient_max_Ha_Bohr','peak_GPU_MiB','energy_difference_Ha','maximum_gradient_difference_Ha_Bohr','initial_state_execution_agreement_passed']}))
