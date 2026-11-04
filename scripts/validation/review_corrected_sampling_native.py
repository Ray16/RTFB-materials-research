#!/usr/bin/env python3
"""Reuse original native scientific reconstruction with current closure/frame QC."""
import argparse
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
from types import SimpleNamespace
from ase.io import read as read_atoms

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify

SOURCE_BYTES=Path(__file__).read_bytes()


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def adapt_resource_checks(body):
    """Only map old guard fields; retain every original scientific expression."""
    changes={"host['status']!='program_returned_review_required'":"host['status']!='effective_memory_guard_finished_profile_review_pending'",
        "host['peak_RSS_bytes']":"max(host['self_VmHWM_bytes'],host['peak_tree_RSS_bytes'])"}
    if body.count(next(iter(changes)))!=1 or body.count("host['peak_RSS_bytes']")!=2:
        raise ValueError('The original independently reviewed native resource template changed')
    for old,new in changes.items():body=body.replace(old,new)
    return body


def review(namespace):
    from redox.workflows.sampling_native import profile
    from redox.validation.execution_review import review_closure,freeze_source,SOURCE_BYTES as closure_bytes
    from redox.validation.scf_frame import checkpoint_frame,checkpoint_spin,freeze_reviewer
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);archive=Path(spec['archive']);used=[namespace]
    def read(path):used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json')
    if plan['role']!='corrected_sampling_native' or any(plan.get(k)!=v for k,v in profile().items()):
        raise ValueError('Exact corrected neutral native profile required')
    thermal_path=archive/'src/redox/qm/qrrho.py';thermal=load(thermal_path,'sampling_native_checked_thermal_definition')
    if plan['thermal_protocol']!=thermal.PROTOCOL:
        raise ValueError('Native thermal definition differs from its actual pinned kernel protocol')
    used.append(thermal_path)
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','sampling_native_review_owner')
    def stopped(owner):
        value=owners.alive(owner)
        if value is None:raise ValueError('Source ownership is unverified, not stopped')
        return value
    stage=review_closure(case,read,stopped,
        expected_status='corrected_sampling_native_finished_independent_review_pending',host_GiB=80,
        GPU_MiB=31744,worker=sources/'worker.py')
    producer=read(case/'experimental_native/status.json')
    if producer['owner']!=stage['owner']:raise ValueError('Original and corrected native owners differ')
    used += [freeze_source(case,Path(__file__),SOURCE_BYTES),freeze_reviewer(case),
        freeze_source(case,ROOT/'src/redox/validation/execution_review.py',closure_bytes)]
    legacy_path=archive/'scripts/validation/review_smooth_polymer_native.py'
    legacy=load(legacy_path,'sampling_original_native_reconstruction');used.append(legacy_path)
    body=adapt_resource_checks(inspect.getsource(legacy.review))
    frozen=case/'review_sources'/hashlib.sha256(body.encode()).hexdigest()/'resource_adapted_native_review.py'
    if not frozen.exists():frozen.parent.mkdir(parents=True,exist_ok=True);frozen.write_text(body)
    if frozen.read_text()!=body:raise ValueError('Frozen scientific review adapter differs')
    used.append(frozen)
    def scientific_module(name,path):
        path=Path(path)
        if path==ROOT/'scripts/pipeline/run_polymer_spacing.py':return SimpleNamespace(owner_alive=stopped)
        actual=archive/path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
        used.append(actual);return load(actual,name)
    api=SimpleNamespace(read_plan=legacy.api.read_plan,verify_sources=legacy.api.verify_sources,
        digest=legacy.api.digest,load_module=scientific_module)
    scope=dict(vars(legacy),api=api);local={};exec(compile(body,str(frozen),'exec'),scope,local)
    numerical_path=frozen.with_name('numerical_review.json')
    if numerical_path.exists():raise FileExistsError('Preserve earlier independent native reconstructions')
    local['review'](case,numerical_path);numerical=read(numerical_path)
    raw=Path(plan['raw_root']);atoms=read_atoms(raw/'geometry.xyz');used.append(raw/'geometry.xyz')
    points=sorted((raw/'SCF_evaluations').glob('*/evaluation.json'))
    frames=sorted((raw/'SCF_frames').glob('*/*/SCF.h5'))
    if len(points)!=len(frames):raise ValueError('Every native gradient needs its actual SCF checkpoint')
    spins=[];frame_errors=[]
    for point,checkpoint in zip(points,frames):
        evaluation=read(point);settings=read(checkpoint.with_name('settings.json'))
        if (settings['xc']!='r2scan' or settings['disp']!='d4' or settings['conv_tol']!=1e-12
                or settings['conv_tol_grad']!=1e-8 or settings['max_cycle']!=200
                or settings['grid_level']!=3 or settings['solvent_class']!='SMD'):
            raise ValueError('Actual native solver settings differ')
        frame=checkpoint_frame(checkpoint,charge=0,mult=1,basis='def2-svp',
            coordinates_Bohr=evaluation['coordinates_Bohr'],symbols=atoms.get_chemical_symbols(),maximum_error_Bohr=1e-12)
        import h5py
        with h5py.File(checkpoint,'r') as stream:energy=float(stream['scf/e_tot'][()])
        if abs(energy-evaluation['E_Ha'])>1e-10:raise ValueError('Actual native checkpoint energy differs')
        frame_errors.append(frame['maximum_coordinate_error_Bohr']);spins.append(checkpoint_spin(checkpoint));used.append(checkpoint)
    solver=read(raw/'execution/response_solver.json')
    if (solver['status']!='response_equations_solved_full_physical_QC_separate' or solver['equations_changed']
            or solver['level_shift_added'] or not solver['calls']
            or any(r['solver_info']!=0 or r['true_residual_norm']>r['solver_absolute_target'] for r in solver['calls'])):
        raise ValueError('Original-equation native response residual checks failed')
    implementation=read(case/'screening_implementation/implementation.json')
    if implementation['installed_library_files_modified'] or implementation['physics_model_changed']:
        raise ValueError('Only the isolated screening correction is permitted')
    restart=read(raw/'restart_frame.json')
    if not restart['CPHF_grids_matched'] or restart['native_G_qualified']:
        raise ValueError('Exact saved center/grid frame required without implied qualification')
    for key,hashkey in [('SCF_checkpoint','checkpoint_sha256'),('grids','grids_sha256')]:
        path=Path(restart[key]);used.append(path)
        if sha256_file(path)!=restart[hashkey]:raise ValueError('Saved center/grid frame changed')
    verify(namespace)
    partition=read(raw/'experimental_partition.json');electronic=read(Path(plan['electronic_review']))
    return dict(status='corrected_sampling_native_independently_reconstructed',namespace=str(namespace),key=plan['source_id'],
        task=plan['tasks'][0],experimental_method_id=plan['experimental_method_id'],experimental_model=plan['experimental_model'],
        electronic_protocol=plan['electronic_protocol'],native_settings=plan['experimental_native_settings'],
        electronic_E_smd_eV=electronic['E_smd_eV'],nuclear_by_cutoff_eV=partition['nuclear']['nuclear_by_cutoff_eV'],
        spin_partition_eV=partition['spin_partition_eV'],temperature_K=298.15,
        thermal_protocol=plan['thermal_protocol'],standard_state=plan['thermal_protocol']['pressure'],
        thermal_source=partition['nuclear'],geometry_sha256=plan['scored_geometry_sha256'],charge=0,multiplicity=1,
        G_by_cutoff_eV=numerical['G_by_cutoff_eV'],minimum_frequency_cm=numerical['minimum_frequency_cm'],
        maximum_column_error=numerical['maximum_sampled_column_error'],actual_SCF_frames_checked=len(frames),
        observed_peak_GPU_MiB=stage['observed_peak_GPU_MiB'],peak_host_RSS_bytes=numerical['peak_host_RSS_bytes'],
        maximum_frame_error_Bohr=max(frame_errors),occupied_wavefunction_spins=spins,
        response_blocks_checked=len(solver['calls']),native_G_qualified=True,
        primary_candidate_admission=False,conformer_ensemble_qualified=False,
        source_sha256={**numerical['source_sha256'],**{str(p):sha256_file(p) for p in used}})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--namespace',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier native reviews')
    result=review(a.namespace.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','key','actual_SCF_frames_checked','minimum_frequency_cm')}))
