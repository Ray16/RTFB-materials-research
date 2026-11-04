#!/usr/bin/env python3
"""Qualify an isolated screening-kernel correction against independent masks."""
import argparse
import ctypes
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def synthetic_case(constants):
    import numpy as np
    rng=np.random.default_rng(20261031);coords=[]
    for repeat in range(2):
        for tile in range(16):coords.append(np.column_stack([np.full(256,30.*tile),np.zeros(256),np.zeros(256)])+rng.normal(0,.05,(256,3)))
    coords=np.asfortranarray(np.vstack(coords));atm=np.zeros((17,constants['ATM_SLOTS']),dtype=np.int32);env=[]
    for i in range(17):
        atm[i,constants['PTR_COORD']]=len(env);env.extend([30.*i if i<16 else 1000.,0.,0.])
    bas=np.zeros((512,constants['BAS_SLOTS']),dtype=np.int32)
    for shell in range(len(bas)):
        bas[shell,constants['ATOM_OF']]=shell//32 if shell<480 else 16;bas[shell,constants['ANG_OF']]=shell%3
        nprim=1+shell%4;bas[shell,constants['NPRIM_OF']]=nprim
        bas[shell,constants['PTR_EXP']]=len(env);env.extend(np.linspace(.5,1.5,nprim));bas[shell,constants['PTR_COEFF']]=len(env);env.extend(np.full(nprim,1./nprim))
    return coords,atm,bas,np.asarray(env,dtype=np.float64)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);args=p.parse_args()
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();plan=json.loads(args.plan.read_text());case=args.plan.parent;sources=case.parent/'sources'
    if plan['role']!='screening_barrier_control' or any(hashlib.sha256(Path(path).read_bytes()).hexdigest()!=value for path,value in plan['source_sha256'].items()):raise ValueError('Control source changed')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):raise ValueError('Two CPUs and six matching native caps required')
    record=dict(status='screening_barrier_control_running',started=time.time(),owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),native_G_qualified=False,physics_model_changed=False)
    atomic_write_json(case/'status.json',record)
    try:
        import cupy as cp
        import numpy as np
        cp.get_default_memory_pool().set_limit(size=2*2**30 if plan.get('water_reference') else 512*2**20)
        from pyscf.gto import mole
        from gpu4pyscf.dft import gen_grid
        helper=load(sources/'helper.py','frozen_screening_barrier_helper');constants={name:int(getattr(mole,name)) for name in ('ATM_SLOTS','PTR_COORD','BAS_SLOTS','ATOM_OF','ANG_OF','NPRIM_OF','PTR_EXP','PTR_COEFF')}
        coords,atm,bas,env=synthetic_case(constants);log_cutoff=float(np.log(1e-12));block_size=4096
        reference=helper.reference_mask(coords,atm,bas,env,block_size,log_cutoff,constants);np.save(case/'independent_mask.npy',reference,allow_pickle=False)
        cuda=helper.standalone_source((sources/'kernel.cu').read_text(),constants);kernel=cp.RawKernel(cuda,'redox_screen_index',options=('--std=c++11',));kernel.compile()
        arrays=[cp.asarray(value,order='F' if i==0 else 'C') for i,value in enumerate([coords,atm,bas,env])];gpu_coords,gpu_atm,gpu_bas,gpu_env=arrays;stream=cp.cuda.get_current_stream()
        rows=[]
        for variant in ('installed','corrected'):
            for index in range(30):
                mask=cp.zeros(reference.shape,dtype=np.int8)
                if variant=='installed':
                    code=gen_grid.libgdft.GDFTscreen_index(ctypes.c_void_p(stream.ptr),ctypes.c_void_p(mask.data.ptr),ctypes.c_double(log_cutoff),ctypes.c_void_p(gpu_coords.data.ptr),ctypes.c_int(len(coords)),ctypes.c_int(block_size),ctypes.c_void_p(gpu_atm.data.ptr),ctypes.c_int(len(atm)),ctypes.c_void_p(gpu_bas.data.ptr),ctypes.c_int(len(bas)),ctypes.c_void_p(gpu_env.data.ptr))
                    if code:raise RuntimeError('Installed screening launch failed')
                else:helper.launch(mask,log_cutoff,gpu_coords,gpu_atm,gpu_bas,gpu_env,block_size,kernel,stream)
                actual=cp.asnumpy(mask);mismatch=int(np.count_nonzero(actual!=reference));row=dict(variant=variant,repeat=index,mismatched_shell_blocks=mismatch,mask_sha256=hashlib.sha256(actual.tobytes()).hexdigest());rows.append(row)
                atomic_write_json(case/'progress.json',dict(rows=rows,native_G_qualified=False))
                if mismatch:np.save(case/f'{variant}_{index:03d}_mask.npy',actual,allow_pickle=False)
        corrected=[row for row in rows if row['variant']=='corrected']
        if any(row['mismatched_shell_blocks'] for row in corrected):raise ValueError('Corrected screening mask differs from independent CPU bound')
        water_results=[]
        if plan.get('water_reference'):
            from contextlib import nullcontext
            from pyscf import gto
            from gpu4pyscf.dft import numint
            reference_case=json.loads((sources/'water.json').read_text())
            molecule=gto.M(atom=';'.join(reference_case['atoms']),basis=reference_case['basis'],unit=reference_case['unit'],verbose=0)
            density=cp.ones((molecule.nao,molecule.nao))
            for variant in ('installed','corrected'):
                context=nullcontext() if variant=='installed' else helper.corrected_screening(sources/'kernel.cu',case/'water_implementation')
                with context:
                    values=[]
                    for sorted_grids in (True,False):
                        grids=gen_grid.Grids(molecule);grids.atom_grid=tuple(reference_case['atom_grid']);grids.build(sort_grids=sorted_grids)
                        number,xc_energy,_=numint.NumInt().nr_rks(molecule,grids,reference_case['xc'],density)
                        values.append(dict(sorted=sorted_grids,numeric_density_integral=float(number),XC_energy_Ha=float(xc_energy),grid_points=grids.size))
                        del grids;gc.collect();cp.get_default_memory_pool().free_all_blocks()
                    row=dict(variant=variant,values=values,density_integral_difference=values[0]['numeric_density_integral']-values[1]['numeric_density_integral'],XC_energy_difference_Ha=values[0]['XC_energy_Ha']-values[1]['XC_energy_Ha'])
                    water_results.append(row);atomic_write_json(case/'water_comparison.json',dict(rows=water_results,interpretation=reference_case['fixed_density'],native_G_qualified=False))
            fixed=water_results[-1];limit=reference_case['sorted_unsorted_absolute_limit']
            if max(abs(fixed['density_integral_difference']),abs(fixed['XC_energy_difference_Ha']))>limit:raise ValueError('Corrected official integration-order control exceeds its original absolute limit')
        if any(hashlib.sha256(Path(path).read_bytes()).hexdigest()!=value for path,value in plan['source_sha256'].items()):raise ValueError('Screening-control source changed during evaluation')
        record.update(status='screening_barrier_masks_and_requested_integration_controls_checked',comparisons=rows,shape=list(reference.shape),reference_sha256=hashlib.sha256(reference.tobytes()).hexdigest(),installed_mismatches=sum(row['mismatched_shell_blocks'] for row in rows if row['variant']=='installed'),corrected_mismatches=0,water_integration_control=water_results,real_SCF_gradient_Hessian_qualification_complete=False,implementation_qualified=False)
    except BaseException as error:record.update(status='screening_barrier_control_failed_preserved',error=repr(error));raise
    finally:record['finished']=time.time();atomic_write_json(case/'status.json',record)


if __name__=='__main__':main()
