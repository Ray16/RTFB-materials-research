#!/usr/bin/env python3
"""Review closed benchmark energy controls and independently qualify resources."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path

import h5py
import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def review(namespace,placement):
    from ase.io import read as read_atoms
    from pyscf.lib import param
    spec=verify(namespace);case=Path(spec['case']);used=[namespace]
    def read(path):used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json');stage=read(case/'status.json')
    launch=read(placement/'launch.json');completion=read(placement/'completion.json')
    guard=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    definition=importlib.util.spec_from_file_location('energy_control_owner',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    if (stage['status']!='screening_energy_control_finished_review_pending' or stage['native_G_qualified']
            or completion['returncode']!=0 or completion['budget_stop'] is not None
            or owners.alive(stage['owner']) is not False or not guard['guard_stopped']
            or not residency['own_residency_observed'] or residency['foreign_resident_PIDs']
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or host['maximum_RSS_GiB']!=40 or host['minimum_free_GiB']!=32
            or max(host['self_VmHWM_bytes'],host['peak_tree_RSS_bytes'])>40*2**30
            or host['free_RAM_bytes']<32*2**30 or not 0<stage['observed_peak_GPU_MiB']<=31744
            or not stage['owner']['pid']==guard['pid']==runtime['pid']==host['pid']
            or not stage['owner']['process_start']==guard['process_start']==host['process_start']
            or not guard['CPU_affinity']==runtime['CPU_affinity']==host['CPU_affinity']==launch['CPU_affinity']
            or len(host['CPU_affinity'])!=2 or runtime['program']!=str(Path(spec['sources'])/'worker.py')):
        raise ValueError('Closed source-exact energy control and normal resource/owner checks required')
    expected_coordinates=read_atoms(Path(spec['sources'])/'input.xyz').positions/param.BOHR
    energy_results={};checkpoint_counts={}
    for variant in ('installed','corrected'):
        root=case/variant;result=read(root/'result.json')
        if (result['energy_protocol']!=plan['electronic_protocol'] or result['optimized']
                or result['charge']!=plan['source_state']['charge'] or result['mult']!=plan['source_state']['mult']
                or result['converged_smd'] is not True or not math.isfinite(result['e_smd_Ha'])
                or abs(result.get('smd_spin_contam',math.inf))>.1):
            raise ValueError('Matched converged unoptimized electronic parent required')
        energies=[];count=0
        for path in sorted((root/'SCF_evaluations').glob('*/*/SCF_result.json')):
            solved=read(path);settings=read(path.with_name('settings.json'));checkpoint=path.with_name('SCF.h5')
            used.append(checkpoint)
            if (settings['charge'],settings['mult'])!=(plan['source_state']['charge'],plan['source_state']['mult']):continue
            if (settings['basis']!=plan['electronic_protocol']['basis'] or settings['xc']!=plan['electronic_protocol']['xc']
                    or settings['disp']!=plan['electronic_protocol']['disp'] or settings['conv_tol']!=1e-9
                    or settings['max_cycle']!=200 or settings['solvent_class']!='SMD'
                    or np.max(abs(np.asarray(settings['coordinates_Bohr'])-expected_coordinates))>1e-10):
                raise ValueError('Recorded actual SCF settings differ')
            with h5py.File(checkpoint,'r') as handle:
                actual=float(handle['scf/e_tot'][()])
                if not math.isfinite(actual) or abs(actual-solved['energy_Ha'])>1e-9:
                    raise ValueError('Checkpoint and recorded solver energy differ')
            count+=1
            if solved['converged']:energies.append(actual)
        if not energies or abs(min(energies)-result['e_smd_Ha'])>1e-9:
            raise ValueError('The selected lowest converged state is not backed by actual SCF checkpoints')
        checkpoint_counts[variant]=count;energy_results[variant]=result['e_smd_Ha']
        if variant=='corrected':
            implementation=read(root/'implementation/implementation.json')
            if implementation['installed_library_files_modified'] or implementation['physics_model_changed']:
                raise ValueError('Only the isolated screening implementation may change')
            for name in ('screening_kernel.cu','corrected_host.py','original_host.py'):used.append(root/'implementation'/name)
    old=plan['reference_energy_Ha']
    differences=dict(installed_minus_stored_Ha=energy_results['installed']-old,
        corrected_minus_stored_Ha=energy_results['corrected']-old,
        corrected_minus_installed_Ha=energy_results['corrected']-energy_results['installed'])
    if any(abs(v-stage['differences'][k])>1e-14 for k,v in differences.items()):
        raise ValueError('Implementation differences do not reconstruct')
    verify(namespace)
    return dict(status='benchmark_screening_energy_control_independently_checked',namespace=str(namespace),
        key=plan['source_id'],source_state=plan['source_state'],nao=plan['nao'],atoms=plan['atoms'],
        differences=differences,agreement_tolerance_Ha=plan['comparison_absolute_tolerance_Ha'],
        declared_agreement_passed=all(abs(v)<=plan['comparison_absolute_tolerance_Ha'] for v in differences.values()),
        actual_SCF_checkpoint_counts=checkpoint_counts,peak_GPU_MiB=stage['observed_peak_GPU_MiB'],
        peak_host_RSS_bytes=max(host['self_VmHWM_bytes'],host['peak_tree_RSS_bytes']),
        resources_checked=True,native_G_qualified=False,new_potential_qualified=False,
        source_sha256={str(p):sha256_file(p) for p in used},producer_sha256=sha256_file(__file__))


def pilots(paths,manifest_path):
    manifest=json.loads(manifest_path.read_text());states=manifest['states']
    expected={max([r for r in states if (r['source_state']['mult']>1)==open_shell],
        key=lambda r:(r['nao'],-abs(r['source_state']['charge'])))['key'] for open_shell in (False,True)}
    rows=[json.loads(p.read_text()) for p in paths]
    if len(rows)!=2 or {r['key'] for r in rows}!=expected or any(
            r['status']!='benchmark_screening_energy_control_independently_checked' or not r['resources_checked']
            or r['native_G_qualified'] or r['peak_host_RSS_bytes']>40*2**30 or not 0<r['peak_GPU_MiB']<=31744
            for r in rows):raise ValueError('Both largest independently checked RKS/UKS resource points required')
    sources={str(p):sha256_file(p) for p in [manifest_path,*paths]}
    for row in rows:
        if any(sha256_file(p)!=h for p,h in row['source_sha256'].items()):raise ValueError('Pilot backing changed')
        sources.update(row['source_sha256'])
    return dict(status='largest_RKS_UKS_energy_control_resources_independently_checked',keys=sorted(expected),
        peak_GPU_MiB=max(r['peak_GPU_MiB'] for r in rows),peak_host_RSS_bytes=max(r['peak_host_RSS_bytes'] for r in rows),
        chemical_implementation_agreement_required_for_resource_release=False,
        native_G_qualified=False,source_sha256=sources,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('review','pilots'))
    parser.add_argument('--spec',type=Path);parser.add_argument('--placement',type=Path)
    parser.add_argument('--review',type=Path,action='append');parser.add_argument('--manifest',type=Path)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.out.exists():raise FileExistsError('Preserve previous independent qualification reports')
    result=(review(args.spec.resolve(),args.placement.resolve()) if args.action=='review'
            else pilots([p.resolve() for p in args.review],args.manifest.resolve()))
    atomic_write_json(args.out,result)
    print(json.dumps({k:v for k,v in result.items() if k in ('status','key','keys','differences','declared_agreement_passed','peak_GPU_MiB','peak_host_RSS_bytes')}))
