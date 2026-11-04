#!/usr/bin/env python3
"""Compute a first primary native component with explicit sulfur FD coverage."""
import argparse
from contextlib import contextmanager
import gc
import importlib.util
import json
import math
import os
from pathlib import Path
import time


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def recorded_center(case,charge,mult,count):
    """Use the pinned engine's actual center, preserving its unit convention."""
    import numpy as np
    paths=sorted((case/'SCF_evaluations').glob('*/*/settings.json'))
    if not paths:raise ValueError('Actual native center evidence required')
    center=json.loads(paths[0].read_text());coordinates=np.asarray(center['coordinates_Bohr'])
    if (coordinates.shape!=(count,3) or not np.isfinite(coordinates).all()
            or (center['charge'],center['mult'])!=(charge,mult)):
        raise ValueError('The actual native center must define every added displacement')
    return paths[0],coordinates


@contextmanager
def numerical_cds_context(plan,sources,case,sha):
    if not plan.get('CDS_Hessian_numerical_control'):
        yield
        return
    import numpy as np
    from gpu4pyscf.solvent.hessian import smd as hessian_module
    module=load(sources/'cds_finite.py','private_qualified_CDS_derivative_step')
    policy=plan['native_numerical_policy']
    with module.finer_cds_hessian(policy['Hessian_gradient_step_Bohr']):
        finer=hessian_module.get_cds
        def evaluate(solvent):
            matrix=np.asarray(finer(solvent));coordinates=solvent.mol.atom_coords(unit='Bohr')
            event=dict(policy=policy,coordinates_Bohr=coordinates.tolist(),all_finite=bool(np.isfinite(matrix).all()),
                helper_sha256=sha(sources/'cds_finite.py'),shape=list(matrix.shape))
            with (case/'CDS_numerical_events.jsonl').open('a') as stream:stream.write(json.dumps(event)+'\n')
            return matrix
        hessian_module.get_cds=evaluate
        try:yield
        finally:hessian_module.get_cds=finer


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=a.plan.parent;sources=case.parent/'sources';plan=json.loads(a.plan.read_text())
    sha=load(sources/'hashing.py','route_native_hashing').sha256_file
    load(sources/'admission.py','route_native_admission').validate_plan(plan)
    def verify():
        if any(sha(p)!=h for p,h in plan['source_sha256'].items()):raise ValueError('Frozen first-native source changed')
    verify()
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two disjoint CPUs and matching thread caps required')
    owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    stage=dict(status='route_primary_native_component_running',owner=owner,started=time.time(),native_G_qualified=False)
    atomic_write_json(case/'status.json',stage)
    try:
        import numpy as np
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft
        from redox.qm.native_thermal import ensure_native_thermal,cache_path,NATIVE_THERMAL
        from redox.validation.protonation import hydrogen_topology_status
        if NATIVE_THERMAL!=plan['native_protocol'] or hydrogen_topology_status(plan['smiles'],case/'state/opt.xyz')!='ok':
            raise ValueError('Exact primary native protocol and graph required')
        helper=load(sources/'helper.py','route_native_screening')
        checkpoints=load(sources/'checkpoints.py','route_native_checkpoints')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','route_native_memory')
        builder=dft._build_smd_mf;dft._build_smd_mf=checkpoints.recording_builder(builder,case/'SCF_evaluations',owner,save_density=False)
        with helper.corrected_screening(sources/'kernel.cu',case/'screening_implementation'),numerical_cds_context(plan,sources,case,sha),monitor.MemoryMonitor() as memory:
            kind,thermal=ensure_native_thermal(plan['native_id'],plan['native_state'],state_dir=case/'state',memory_cap_gib=6)
            if kind!='computed':raise ValueError('Fresh own native thermochemistry required')
            native_path,_=cache_path(case/'state',plan['charge'],plan['mult'])
            matrix=np.load(native_path.with_suffix('')/'hessian_cartesian.npy',allow_pickle=False)
            atoms=read(case/'state/opt.xyz');symbols=atoms.get_chemical_symbols()
            center_path,coordinates=recorded_center(case,plan['charge'],plan['mult'],39)
            extra=[]
            for atom,axis in plan['extra_columns']:
                if symbols[atom]!='S':raise ValueError('Actual sulfur-centered coverage required')
                for step in (.0005,.001):
                    gradients=[]
                    for sign in (-1,1):
                        point=coordinates.copy();point[atom,axis]+=sign*step
                        mol=gto.M(atom=list(zip(symbols,point.tolist())),unit='Bohr',basis='def2-svp',charge=plan['charge'],spin=plan['mult']-1,verbose=0,max_memory=6000)
                        mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
                        mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.kernel()
                        if not mf.converged:raise ValueError('Sulfur displaced SCF failed')
                        obj=mf.nuc_grad_method();obj.grid_response=True;value=obj.kernel()
                        gradients.append(np.asarray(value.get() if hasattr(value,'get') else value))
                        diag=dft._orbital_diag(mf,plan['mult']-1)
                        if abs(diag.get('spin_contam',math.inf))>.1:raise ValueError('Sulfur displaced spin check failed')
                        del mf,obj;gc.collect();dft._free_gpu_pool()
                    numerical=((gradients[1]-gradients[0])/(2*step)).ravel()
                    error=float(np.linalg.norm(numerical-matrix[3*atom+axis])/np.linalg.norm(numerical))
                    if not math.isfinite(error) or error>.01:raise ValueError('Extra sulfur finite-gradient Hessian column failed')
                    extra.append(dict(atom_index=atom,axis=axis,step_Bohr=step,relative_error=error))
                    atomic_write_json(case/'extra_derivatives.json',extra)
        dft._build_smd_mf=builder
        if not 0<memory.peak<=31744:raise ValueError('Primary native GPU envelope exceeded')
        energy=json.loads((sources/'corrected_electronic.json').read_text())
        values={c:energy['e_smd_eV']+thermal['nuclear_by_cutoff_eV'][c]+thermal['g_spin_partition_eV'] for c in ('25','50','100')}
        if not all(math.isfinite(v) for v in values.values()):raise ValueError('Finite own component sums required')
        result=dict(key=plan['source_id'],charge=plan['charge'],mult=plan['mult'],
            corrected_electronic_eV=energy['e_smd_eV'],native_record_path=str(native_path),
            nuclear_by_cutoff_eV=thermal['nuclear_by_cutoff_eV'],spin_partition_eV=thermal['g_spin_partition_eV'],
            G_by_cutoff_eV=values,extra_derivative_checks=extra,extra_FD_center_settings=str(center_path),
            native_G_qualified=False,redox_potential_qualified=False,material_admission=False)
        if plan.get('CDS_Hessian_numerical_control'):
            result.update(native_numerical_policy=plan['native_numerical_policy'],
                CDS_numerical_events=str(case/'CDS_numerical_events.jsonl'),
                full_Hessian_recomputed=True,old_saved_matrix_reused=False)
        atomic_write_json(case/'result.json',result)
        verify();stage.update(status='route_primary_native_component_finished_review_pending',observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:
        stage.update(status='route_primary_native_component_failed_preserved',error=repr(error));raise
    finally:stage['finished']=time.time();atomic_write_json(case/'status.json',stage)


if __name__=='__main__':main()
