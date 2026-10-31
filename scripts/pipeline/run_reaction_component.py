#!/usr/bin/env python3
"""UMA warm start, matched DFT/SMD and native G for one exact route reagent."""
import argparse
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    from redox.validation.protonation import hydrogen_topology_status
    from ase.io import read,write
    verify_gate_owner();case=a.plan.parent;plan=json.loads(a.plan.read_text());sources=case.parent/'sources'
    if plan['role']!='reaction_component' or any(sha(k)!=v for k,v in plan['source_sha256'].items()):
        raise ValueError('Prepared reagent inputs or source snapshot changed')
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or {os.environ.get(k) for k in caps}!={'2'}:
        raise ValueError('Two disjoint CPUs and six native thread caps required')
    host_capacity(2,minimum_ram_gib=plan['maximum_host_RSS_GiB']+32)
    device=gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    if device['compute_capability']!='7.0' and not plan.get('UMA_cross_device_reference'):
        raise ValueError('A non-V100 standard execution requires a fixed-geometry V100 comparison')
    guard=json.loads((case/'host_guard/status.json').read_text())
    if guard['pid']!=os.getpid() or guard['maximum_RSS_GiB']!=plan['maximum_host_RSS_GiB']:
        raise ValueError('Prospective same-PID host guard required')
    if not plan['route_gate']['proceed_to_full_redox'] or not plan['route_gate']['exact_route_binding_verified']:
        raise ValueError('Exact linked route admission required')
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status=dict(status='reaction_component_running_unqualified',role='reaction_component',owner=owner,
        started=time.time(),plan_sha256=sha(a.plan),source_id=plan['source_id'],native_G_qualified=False,material_admission=False)
    atomic_write_json(case/'status.json',status)
    def stop(signum,frame):raise KeyboardInterrupt('Owned reagent job stopped; retain checkpoints')
    signal.signal(signal.SIGTERM,stop)
    try:
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','reagent_memory_monitor')
        with monitor.MemoryMonitor() as memory:
            from redox.qm.uma_model_identity import model_identity
            from redox.qm.uma_current import calculator_for_profile
            from ase.optimize import BFGS
            policy=json.loads((ROOT/'config/uma_execution.json').read_text())
            os.environ['FAIRCHEM_CACHE_DIR']=policy['fairchem_cache']
            identity=model_identity(policy['model'])
            if (identity['checkpoint_sha256']!=policy['checkpoint_sha256']
                    or identity['reference_sha256']!=policy['reference_sha256']):raise ValueError('UMA weights/references differ')
            if policy['model'] not in json.loads(Path(identity['registry_path']).read_text()):
                raise ValueError('Reagent jobs may not modify the installed UMA registry')
            atoms=read(case/'input.xyz');atoms.info.update(charge=plan['charge'],spin=plan['mult'])
            calc=calculator_for_profile(policy['model'],'cuda','legacy',checkpoint_root=case/'UMA_evaluations')
            if plan.get('UMA_cross_device_reference'):
                import numpy as np
                from ase import Atoms
                reference=json.loads(Path(plan['UMA_cross_device_reference']).read_text())
                if reference['model_identity']!=identity:raise ValueError('Cross-device UMA model identity differs')
                fixture=Atoms(numbers=reference['numbers'],positions=reference['coordinates_A'])
                fixture.info.update(charge=plan['charge'],spin=plan['mult']);fixture.calc=calc
                repeats=[]
                for _ in range(3):
                    calc.reset();force=fixture.get_forces();energy=fixture.get_potential_energy()
                    repeats.append(dict(energy_eV=float(energy),forces_eV_A=force.tolist(),
                        energy_error_eV=abs(energy-reference['energy_eV']),
                        force_max_error_eV_A=float(np.abs(force-np.asarray(reference['forces_eV_A'])).max())))
                atomic_write_json(case/'UMA_cross_device_check.json',dict(reference_sha256=sha(plan['UMA_cross_device_reference']),
                    actual_profile='legacy_standard',compute_capability=device['compute_capability'],repeats=repeats,
                    energy_tolerance_eV=1e-4,force_tolerance_eV_A=1e-4))
                if any(not np.isfinite([r['energy_error_eV'],r['force_max_error_eV_A']]).all()
                       or r['energy_error_eV']>1e-4 or r['force_max_error_eV_A']>1e-4 for r in repeats):
                    raise ValueError('Fixed-geometry standard UMA differs from the V100 energy/forces')
                fixture.calc=None;del fixture
            atoms.calc=calc;recipe=plan['recipe'];warm=case/'uma';warm.mkdir()
            status['stage']='UMA_standard_warm_start';atomic_write_json(case/'status.json',status)
            optimizer=BFGS(atoms,logfile=str(warm/'optimization.log'),trajectory=str(warm/'optimization.traj'))
            converged=optimizer.run(fmax=recipe['UMA_fmax_eV_A'],steps=recipe['UMA_max_steps'])
            forces=atoms.get_forces();maximum=float((forces**2).sum(axis=1).max()**.5)
            write(warm/'relaxed.xyz',atoms)
            warm_record=dict(status='UMA_warm_start_only',charge=plan['charge'],mult=plan['mult'],
                converged=bool(converged),fmax_eV_A=maximum,energy_eV=float(atoms.get_potential_energy()),
                model_identity=identity,actual_mode='legacy_standard',compute_capability=device['compute_capability'],steps=optimizer.nsteps,
                geometry_sha256=sha(warm/'relaxed.xyz'),solution_energy=False)
            atomic_write_json(warm/'result.json',warm_record)
            if not converged or maximum>recipe['UMA_fmax_eV_A'] or hydrogen_topology_status(plan['smiles'],warm/'relaxed.xyz')!='ok':
                raise ValueError('UMA geometry or convergence gate failed; warm-start output preserved')
            atoms.calc=None;del calc,optimizer,atoms;gc.collect()
            import torch
            torch.cuda.empty_cache()
            # Recheck the numerical engine after importing Torch/FAIRChem. Its
            # outputs remain warm starts; DFT keeps the declared FP64 runtime.
            tensor=load(ROOT/'scripts/fleet/tensor_runtime_worker.py','reagent_post_UMA_tensor')
            atomic_write_json(case/'post_UMA_runtime.json',tensor.validate(json.loads((ROOT/'environments/cutensor22/runtime.json').read_text())))
            checkpoint=load(sources/'checkpoints.py','reagent_SCF_checkpoints')
            from redox.qm import dft
            dft._build_smd_mf=checkpoint.recording_builder(dft._build_smd_mf,case/'SCF_evaluations',owner)
            from pyscf.geomopt import geometric_solver
            geometric_solver.kernel=checkpoint.recording_optimizer(geometric_solver.kernel,case/'optimizer_returns',owner)
            directory=Path(plan['target_directory']);directory.mkdir(parents=True)
            status['stage']='DFT_SMD_geometry_and_uniform_diffuse_SP';atomic_write_json(case/'status.json',status)
            result=dft.dft_smd(warm/'relaxed.xyz',plan['charge'],plan['mult'],
                sp_basis='def2-tzvpd',sp_basis_anion='def2-tzvpd',solvent='acetonitrile',
                do_gas=True,do_smd=True,do_opt=True,do_freq=False,opt_out=directory/'opt.xyz',backend='gpu')
            result.update(id=plan['source_id'],state=plan['state'],role='reaction_component',smiles=plan['smiles'])
            atomic_write_json(directory/'result.json',result)
            if not result['optimized'] or hydrogen_topology_status(plan['smiles'],directory/'opt.xyz')!='ok':
                raise ValueError('Optimized DFT reagent connectivity/convergence failed')
            initial_optimized_sha=sha(directory/'opt.xyz')
            from redox.qm.sp import ensure
            ensure(plan['source_id'],plan['state'],directory,plan['charge'],plan['mult'])
            from redox.qm.native_thermal import select_policy,ensure_native_thermal,load_native_thermal,free_energy_native,cache_path
            from redox.qm import native_thermal_repair
            installer=load(sources/'native_worker.py','reagent_native_physical_installer')
            native_thermal_repair.install_repaired_state=installer.physical_installer(native_thermal_repair.install_repaired_state,ROOT,case)
            select_policy(directory);status['stage']='native_DFT_Hessian_and_partition';atomic_write_json(case/'status.json',status)
            try:ensure_native_thermal(plan['source_id'],plan['state'])
            except RuntimeError as error:
                if not any(s in str(error) for s in ('not stationary','not a minimum')):raise
                native_thermal_repair.repair_minimum(plan['source_id'],plan['state'],str(error))
            from redox.core.protocol import load_record,record_path
            energy=load_record(directory,plan['charge'],plan['mult']);thermal=load_native_thermal(directory,plan['charge'],plan['mult'])
            total,quality=free_energy_native(plan['source_id'],plan['state'])
        if energy is None or thermal is None or total is None or quality!='ok' or not 0<memory.peak<=plan['maximum_peak_GPU_MiB']:
            raise ValueError('Reagent native/energy/resource QC failed')
        surface,_=cache_path(directory,plan['charge'],plan['mult'])
        produced=dict(status='reaction_component_native_finished_independent_review_pending',role='reaction_component',
            source_id=plan['source_id'],state=plan['state'],charge=plan['charge'],mult=plan['mult'],inchikey=plan['inchikey'],
            route_candidate=plan['route_candidate'],G_eV=total,active_E_smd_eV=energy['e_smd_eV'],
            G_by_cutoff_eV={str(c):energy['e_smd_eV']+thermal['nuclear_by_cutoff_eV'][str(c)]+thermal['g_spin_partition_eV'] for c in (25,50,100)},
            native_record_path=str(surface),active_record_path=str(record_path(directory,plan['charge'],plan['mult'])),
            final_geometry_sha256=sha(directory/'opt.xyz'),initial_optimized_geometry_sha256=initial_optimized_sha,
            observed_peak_GPU_MiB=memory.peak,native_G_qualified=False,material_admission=False,
            global_conformer_convergence_proven=False)
        atomic_write_json(case/'result.json',produced);status.update(status=produced['status'],observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:
        status.update(status='reaction_component_failed_preserved',error=repr(error));raise
    finally:
        status['finished']=time.time();atomic_write_json(case/'status.json',status)


if __name__=='__main__':main()
