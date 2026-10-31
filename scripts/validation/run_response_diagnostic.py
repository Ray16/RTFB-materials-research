#!/usr/bin/env python3
"""Trace a legacy frozen-solvent/solute-H1 execution diagnostic; no native QC."""
import argparse
import faulthandler
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import time


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    verify_gate_owner();host_capacity(2,minimum_ram_gib=112);gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two disjoint CPUs and all six thread caps=2 are required')
    case=args.plan.parent;plan=json.loads(args.plan.read_text())
    if (plan['role']!='orbital_response_diagnostic' or plan.get('available_response_memory_cap_GiB')!=3
            or plan.get('original_native_memory_policy_context_required') is not True
            or any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items())):raise ValueError('Diagnostic source identity or prospective batching profile changed')
    helper=load(case.parent/'sources/helper.py','frozen_response_diagnostic_helpers')
    checkpoint,arrays,h1=helper.read_checkpoint(Path(plan['raw_source']))
    if checkpoint!=plan['checkpoint']:raise ValueError('Read-only response checkpoint changed')
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    record=dict(status='orbital_response_diagnostic_running',started=time.time(),owner=owner,native_G_qualified=False,phase='initializing')
    def phase(value):record.update(phase=value,updated=time.time());atomic_write_json(case/'status.json',record)
    if plan.get('debugger') is not None:
        phase('waiting_for_explicit_owned_debugger')
        debugger=load(case.parent/'sources/debugger.py','frozen_response_owned_debugger')
        debugger.permit(case,plan['debugger']['sha256'])
    phase('imports')
    import cupy as cp
    import numpy as np
    from pyscf import gto
    from redox.qm import dft
    from gpu4pyscf.dft.rks import initialize_grids
    from gpu4pyscf.hessian import rhf as response_rhf
    adapter=load(ROOT/'scripts/pipeline/residual_response_solver.py','diagnostic_original_response_solver')
    model=load(ROOT/'scripts/pipeline/smooth_cds_control.py','diagnostic_original_smooth_model')
    monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','diagnostic_original_GPU_monitor')
    native=load(ROOT/'scripts/pipeline/smooth_polymer_native.py','diagnostic_original_native_memory_policy')
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    fatal=(case/'fatal_stack.log').open('x');faulthandler.enable(file=fatal,all_threads=True)
    checkpoint_identity=dict(source_namespace=plan['source_namespace'],checkpoint=plan['checkpoint'])
    try:
        with helper.trace_solver(adapter,case/'trace',checkpoint_root=case/'response_blocks' if plan.get('checkpoint_response') else None,checkpoint_identity=checkpoint_identity),monitor.MemoryMonitor() as memory,\
                native.execution_budget(),helper.response_memory_cap(response_rhf,3*2**30),\
                adapter.compact_smd_execution(case/'execution',28),\
                model.smooth_cds(case/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'):
            phase('recreating_saved_MF_and_initial_grid_recipe')
            symbols=plan['symbols'];coords=checkpoint['identity']['coordinates_Bohr']
            mol=gto.M(atom=list(zip(symbols,coords)),unit='Bohr',basis='def2-svp',charge=0,spin=0,verbose=0,max_memory=6000)
            mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3
            # This auxiliary probe reproduces the historical execution model.
            # Full native Hessians instead use Hessian.kernel(), which enables
            # equilibrium response and adds PCM to the solute H1 contribution.
            mf.with_solvent.equilibrium_solvation=False
            mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
            # Original grids are not present in the historical checkpoint.
            # Recreate the original default initial-guess/pruning recipe and
            # retain that limitation; this is not an exact Hessian resume.
            initial=mf.get_init_guess(mol,mf.init_guess);initialize_grids(mf,mol,initial)
            for name,value in arrays.items():setattr(mf,name,cp.asarray(value))
            mf.e_tot=checkpoint['energy_Ha'];mf.converged=True;density=mf.make_rdm1()
            phase('checking_saved_density_energy_on_recreated_grid')
            potential=mf.get_veff(mol,density);energy=float(mf.energy_tot(density,mf.get_hcore(),potential))
            difference=energy-checkpoint['energy_Ha']
            atomic_write_json(case/'energy_reproduction.json',dict(saved_energy_Ha=checkpoint['energy_Ha'],recreated_energy_Ha=energy,difference_Ha=difference,original_grid_arrays_available=False,native_G_qualified=False))
            if not np.isfinite(energy) or abs(difference)>1e-8:raise ValueError('Recreated fixed-density energy does not reproduce the source within 1e-8 Hartree')
            mf.cphf_grids=mf.grids.copy();hessian=mf.Hessian();hessian.grid_response=True;hessian.max_cycle=100
            if plan.get('checkpoint_response'):
                np.savez(case/'recreated_grids.npz',coords=cp.asnumpy(mf.grids.coords),weights=cp.asnumpy(mf.grids.weights))
                checkpoint_identity['recreated_grids']=dict(path=str(case/'recreated_grids.npz'),sha256=helper.digest(case/'recreated_grids.npz'))
            phase('building_orbital_response_operator')
            operator=hessian.gen_vind(mf.mo_coeff,mf.mo_occ)
            phase('solving_orbital_response_from_saved_first_derivatives')
            mo1,e1=hessian.solve_mo1(mf.mo_energy,mf.mo_coeff,mf.mo_occ,h1,operator)
            hashes={}
            for name,value in [('mo1',mo1),('e1',e1)]:
                array=cp.asnumpy(value) if isinstance(value,cp.ndarray) else np.asarray(value)
                if not np.isfinite(array).all():raise ValueError('Nonfinite diagnostic orbital response')
                hashes[name]=dict(shape=list(array.shape),sha256=hashlib.sha256(array.tobytes()).hexdigest())
                if plan.get('checkpoint_response'):hashes[name]['checkpoint']=helper.save_array(case/(name+'.npy'),array)
            atomic_write_json(case/'response_summary.json',dict(arrays=hashes,original_grid_arrays_available=False,native_G_qualified=False))
        if not 0<memory.peak<=31744:raise ValueError('Diagnostic GPU memory ceiling failed')
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):raise ValueError('Diagnostic source changed during execution')
        record.update(status='orbital_response_diagnostic_finished_no_native_admission',observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:record.update(status='orbital_response_diagnostic_failed_preserved',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(case/'status.json',record)
        faulthandler.disable();fatal.close()


if __name__=='__main__':main()
