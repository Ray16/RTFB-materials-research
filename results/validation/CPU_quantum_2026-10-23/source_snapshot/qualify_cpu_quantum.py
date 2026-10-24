#!/usr/bin/env python
"""Compare CPU execution with stored full GPU spin controls before CPU fallback."""
import argparse
import json
import os
from pathlib import Path
import socket
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def run(controls, out, threads):
    with allocation_guard():
        host_capacity(threads, minimum_ram_gib=32)
        os.sched_setaffinity(0, cpu_allocation(threads))
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if any(os.environ.get(k) != str(threads) for k in caps):
        raise ValueError('CPU affinity and all native thread caps must match')
    import numpy as np
    from ase.io import read
    import pyscf.lib
    from pyscf import gto
    from redox.qm import dft
    pyscf.lib.num_threads(threads)
    execution_path=ROOT/'scripts/pipeline/cpu_quantum_execution.py'
    execution=api.load_module('qualified_CPU_execution',execution_path)
    checker=api.load_module('CPU_GPU_reference_checker',ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py')
    witnesses=checker.checked_controls(controls,28)
    out.mkdir(parents=True,exist_ok=False)
    report=dict(status='CPU_GPU_comparison_running',host=socket.gethostname(),pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),cases=[],
        tolerances=dict(active_energy_eV=1e-5,geometry_energy_Ha=1e-7,gradient_Ha_Bohr=1e-6),
        source_sha256={**witnesses,**{str(p):api.digest(p) for p in
            (Path(__file__).resolve(),execution_path,Path(dft.__file__))}},
        full_polymer_resource_qualified=False,native_Hessian_qualified=False,G_computed=False)
    atomic_write_json(out/'status.json',report)
    try:
        with execution.cpu_execution(out/'execution',8192):
            for mult in (1,2):
                control=controls/'controls'/str(mult)
                reference=json.loads((control/'status.json').read_text())['arms']['original']
                xyz=control/'original/opt.xyz'
                active=dft.dft_smd(xyz,0,mult,backend='cpu',do_opt=False,do_freq=False)
                if not active['converged_smd'] or not active['converged_gas']:
                    raise ValueError('CPU active SCF did not converge')
                atoms=read(xyz)
                mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),
                    basis='def2-svp',charge=0,spin=mult-1,verbose=0)
                mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','cpu',disp='d4')
                mf.grids.level=3
                mf.conv_tol,mf.conv_tol_grad,mf.max_cycle=1e-12,1e-8,200
                mf.chkfile=str(out/f'geometry_mult{mult}.h5')
                mf.kernel()
                if not mf.converged:
                    raise ValueError('CPU geometry-level SCF did not converge')
                gradient=mf.nuc_grad_method();gradient.grid_response=True
                values=np.asarray(gradient.kernel())
                expected=np.load(control/'original/gradient.npy',allow_pickle=False)
                differences=dict(active_energy_eV=max(abs(active['e_'+phase+'_eV']-reference['E_'+phase+'_eV'])
                    for phase in ('smd','gas')),
                    geometry_energy_Ha=abs(float(mf.e_tot)-reference['geometry_level_E_Ha']),
                    gradient_Ha_Bohr=float(np.max(abs(values-expected))))
                if any(not np.isfinite(v) or v>report['tolerances'][k] for k,v in differences.items()):
                    raise ValueError('CPU/GPU equivalence failed: '+repr(differences))
                np.save(out/f'gradient_mult{mult}.npy',values)
                atomic_write_json(out/f'active_mult{mult}.json',active)
                report['cases'].append(dict(mult=mult,differences=differences))
                atomic_write_json(out/'status.json',report)
        api.verify_sources(report['source_sha256'])
        report['status']='two_spin_CPU_GPU_energy_and_gradient_equivalence_passed'
    except BaseException as error:
        report.update(status='CPU_GPU_comparison_failed_preserved',error=repr(error))
        raise
    finally:
        report['finished']=time.time()
        atomic_write_json(out/'status.json',report)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threads',type=int,required=True)
    a=p.parse_args();run(a.controls.resolve(),a.out.resolve(),a.threads)


if __name__=='__main__':
    main()
