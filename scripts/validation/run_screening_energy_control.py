#!/usr/bin/env python3
"""Repeat one frozen benchmark electronic state under both screening variants.

Uses the original SMD single-point solver, including multiple initial guesses
and ferrocene occupation seeds. It never combines new energies with old thermal
terms into a qualified potential or reaction free energy.
"""
import argparse
from contextlib import nullcontext
import gc
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time


def load(path, name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    case=args.plan.parent;sources=case.parent/'sources'
    plan=json.loads(args.plan.read_text())
    if plan['role']!='screening_energy_control' or plan['native_G_qualified']:
        raise ValueError('A separately identified electronic implementation control is required')
    def check_sources():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):
            raise ValueError('Frozen energy-control source changed')
    check_sources()
    if len(os.sched_getaffinity(0))!=2 or os.environ.get('NUMPY_MADVISE_HUGEPAGE')!='0' or any(
            os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS',
                'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two CPUs and matching native runtime caps required')
    owner=dict(host=os.uname().nodename,pid=os.getpid(),
               process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    record=dict(status='screening_energy_control_running',owner=owner,started=time.time(),
                source_state=plan['source_state'],native_G_qualified=False,production_default_changed=False)
    atomic_write_json(case/'status.json',record)
    try:
        import cupy as cp
        from ase.io import read
        from redox.qm import dft
        from redox.core.protocol import ACTIVE_SP
        helper=load(sources/'helper.py','frozen_energy_screening_helper')
        checkpoints=load(sources/'checkpoints.py','frozen_energy_checkpoint_helper')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','frozen_energy_memory_monitor')
        state=plan['source_state'];protocol=plan['electronic_protocol'];results={}
        if protocol!=ACTIVE_SP:raise ValueError('Benchmark active electronic protocol changed')
        atoms=read(sources/'input.xyz')
        if len(atoms)!=plan['atoms']:raise ValueError('Frozen molecular size differs')
        original_builder=dft._build_smd_mf
        with monitor.MemoryMonitor() as memory:
            for variant in ('installed','corrected'):
                folder=case/variant;folder.mkdir()
                dft._build_smd_mf=checkpoints.recording_builder(original_builder,folder/'SCF_evaluations',owner,save_density=False)
                context=nullcontext() if variant=='installed' else helper.corrected_screening(
                    sources/'kernel.cu',folder/'implementation')
                try:
                    with context:
                        result=dft.dft_smd(sources/'input.xyz',state['charge'],state['mult'],
                            sp_xc=protocol['xc'],sp_basis=protocol['basis'],sp_basis_anion=protocol['basis'],
                            sp_disp=protocol['disp'],sp_nlc=protocol['nlc'],solvent=protocol['solvent'],
                            do_opt=False,do_freq=False,do_gas=False,do_smd=True,backend='gpu')
                finally:
                    dft._build_smd_mf=original_builder
                if (result['converged_smd'] is not True or not math.isfinite(result['e_smd_Ha'])
                        or result['energy_protocol']!=protocol
                        or abs(result.get('smd_spin_contam',math.inf))>.1):
                    raise ValueError('Matched SCF/protocol/spin gate failed')
                atomic_write_json(folder/'result.json',result)
                results[variant]=result
                record['finished_variants']=list(results)
                atomic_write_json(case/'status.json',record)
                gc.collect();cp.get_default_memory_pool().free_all_blocks()
        if not 0<memory.peak<=31744:
            raise ValueError('Electronic control exceeds the declared 32-GB GPU profile')
        old=plan['reference_energy_Ha'];installed=results['installed']['e_smd_Ha'];corrected=results['corrected']['e_smd_Ha']
        differences=dict(installed_minus_stored_Ha=installed-old,corrected_minus_stored_Ha=corrected-old,
                         corrected_minus_installed_Ha=corrected-installed)
        limit=plan['comparison_absolute_tolerance_Ha']
        check_sources()
        record.update(status='screening_energy_control_finished_review_pending',differences=differences,
            declared_agreement_passed=all(abs(v)<=limit for v in differences.values()),
            observed_peak_GPU_MiB=memory.peak,chemical_free_energy_assembled=False)
    except BaseException as error:
        record.update(status='screening_energy_control_failed_preserved',error=repr(error))
        raise
    finally:
        record['finished']=time.time();atomic_write_json(case/'status.json',record)


if __name__=='__main__':main()
