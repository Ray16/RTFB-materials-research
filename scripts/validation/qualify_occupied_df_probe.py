#!/usr/bin/env python
"""Isolate the large fixed-orbital DF derivative before another native retry.

The saved solvent-SCF orbitals are inputs. This runs no SCF, optimization, XC,
core or PCM derivative and establishes no complete Hessian or free energy.
"""
import argparse
from contextlib import redirect_stdout, redirect_stderr
import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / 'scripts/pipeline/df_occupied_workspace.py'
CHECKER = ROOT / 'scripts/validation/qualify_df_occupied_workspace.py'


def prepare(source, controls, cap, out):
    plan = api.read_plan(source)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    status_path = source / 'streamed_native/status.json'; status = json.loads(status_path.read_text())
    if status['status'] != 'streamed_native_pilot_failed_preserved' or not status.get('finished'):
        raise ValueError('This probe is scoped to a preserved failed native attempt')
    owner = api.load_module('DF_probe_stopped_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    if owner.owner_alive(status['owner']):
        raise ValueError('The source scientific owner must have ended')
    checkpoint = source / 'SCF_checkpoints/SCF_0000/SCF.h5'
    settings = source / 'SCF_checkpoints/SCF_0000/settings.json'
    inputs = {str(p): api.digest(p) for p in (source / 'plan.json', status_path, checkpoint, settings)}
    identity = json.loads(settings.read_text())
    if (identity['charge'], identity['mult'], identity['xc'], identity['solvent']) != (0, 1, 'r2scan', 'acetonitrile'):
        raise ValueError('The saved solvent-SCF state differs from the declared probe')
    proof = api.load_module('DF_probe_full_spin_caps', CHECKER).checked_controls(controls)
    out.mkdir(parents=True, exist_ok=False)
    record = dict(source=str(source), controls=str(controls), occupied_cap=cap,
        checkpoint=str(checkpoint), settings=str(settings), input_sha256=inputs,
        producer_sha256={**proof, **{str(p): api.digest(p) for p in (Path(__file__).resolve(), ADAPTER, CHECKER)}},
        raw_root=str(ROOT / 'calcs/validation' / out.parent.name / out.name),
        maximum_GPU_MiB=14336, minimum_GPU_MiB=15360, maximum_RSS_GiB=64,
        scope='Restored fixed-solvated-orbital DF potential derivative only; no SCF or full Hessian',
        native_G_qualified=False)
    atomic_write_json(out / 'plan.json', record)


def run(out):
    from redox.workflows.resources import verify_gate_owner, host_capacity
    verify_gate_owner()
    caps = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in caps):
        raise ValueError('Two disjoint CPUs and native caps of two are required')
    host_capacity(2, minimum_ram_gib=96)
    plan = api.read_plan(out)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    device_MiB = float(subprocess.check_output(['nvidia-smi', '--id=' + os.environ['CUDA_VISIBLE_DEVICES'],
        '--query-gpu=memory.total', '--format=csv,noheader,nounits'], text=True).strip())
    guard = json.loads((out / 'host_guard/status.json').read_text())
    if (device_MiB < plan['minimum_GPU_MiB'] or guard['pid'] != os.getpid()
            or guard['maximum_RSS_GiB'] != plan['maximum_RSS_GiB'] or guard['minimum_free_GiB'] != 32):
        raise ValueError('The declared device and same-PID host guard are required')
    raw = Path(plan['raw_root']); raw.mkdir(parents=True, exist_ok=False)
    record = dict(status='isolated_DF_derivative_probe_running', started=time.time(),
        host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), GPU=os.environ['CUDA_VISIBLE_DEVICES'],
        SCF_performed=False, full_Hessian_computed=False, native_G_qualified=False)
    atomic_write_json(out / 'status.json', record)
    memory = None
    try:
        import numpy as np
        import cupy as cp
        from pyscf.scf import chkfile
        from gpu4pyscf.df.hessian import rhf
        from redox.qm import dft
        mol, state = chkfile.load_scf(plan['checkpoint'])
        source_settings = json.loads(Path(plan['settings']).read_text())
        if (mol.charge, mol.spin, mol.natm) != (0, 0, 139):
            raise ValueError('Unexpected saved molecule/state')
        if not np.allclose(mol.atom_coords(), source_settings['coordinates_Bohr'], rtol=0, atol=1e-10):
            raise ValueError('SCF molecule and durable source coordinates differ')
        coeff, occ = np.asarray(state['mo_coeff']), np.asarray(state['mo_occ'])
        if (coeff.shape != (mol.nao_nr(), mol.nao_nr()) or occ.shape != (mol.nao_nr(),)
                or not np.isfinite(coeff).all() or not np.isfinite(occ).all()
                or not np.isin(occ, [0., 2.]).all() or int(occ.sum()) != mol.nelectron):
            raise ValueError('Invalid restored molecular orbitals/occupation')
        overlap = mol.intor_symmetric('int1e_ovlp')
        orthogonality = float(np.max(abs(coeff.T @ overlap @ coeff - np.eye(mol.nao_nr()))))
        if orthogonality > 1e-8:
            raise ValueError('Restored orbitals are not orthonormal in the saved AO basis')
        engine = api.load_module('isolated_DF_probe_memory', ROOT / 'scripts/pipeline/run_pi_stacking.py')
        adapter = api.load_module('isolated_DF_probe_checked_operator', ADAPTER)
        with engine.MemoryMonitor() as memory, adapter.compact_smd_execution(raw / 'execution', 28,
                occupied_cap=plan['occupied_cap']), (raw / 'calculation.log').open('x') as log, redirect_stdout(log), redirect_stderr(log):
            mf = dft._build_smd_mf(mol, 'r2scan', 'acetonitrile', 'gpu', disp='d4')
            mf.with_df.build(build_cderi=False)
            result = rhf._get_veff(mf.with_df.intopt, cp.asarray(coeff), cp.asarray(occ), 1., 0.)
            if result.shape != (mol.natm, 3, mol.nao_nr(), int((occ > 0).sum())):
                raise ValueError('Unexpected DF derivative dimensions')
            value = cp.asnumpy(result)
            if not np.isfinite(value).all():
                raise ValueError('Nonfinite fixed-orbital derivative')
            np.save(raw / 'DF_derivative_MO_occ.npy', value)
            atomic_write_json(raw / 'fixed_state.json', dict(source_checkpoint=plan['checkpoint'],
                source_E_Ha=float(state['e_tot']), coordinates_Bohr=mol.atom_coords().tolist(),
                basis=repr(mol.basis), charge=0, mult=1, occupied_cap=plan['occupied_cap'],
                orbital_orthogonality_maximum_error=orthogonality, SCF_performed=False,
                nocc=int((occ > 0).sum()), shape=list(value.shape)))
        if not 0 < memory.peak <= plan['maximum_GPU_MiB']:
            raise ValueError('Isolated subroutine resource ceiling exceeded')
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
        record.update(status='isolated_DF_derivative_resource_probe_complete_comparison_pending',
            output_sha256={str(p): api.digest(p) for p in raw.rglob('*') if p.is_file()},
            orbital_orthogonality_maximum_error=orthogonality)
    except BaseException as error:
        record.update(status='isolated_DF_derivative_probe_failed_preserved', error=repr(error)); raise
    finally:
        record.update(finished=time.time(), observed_peak_GPU_MiB=memory.peak if memory else None)
        atomic_write_json(out / 'status.json', record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--source', type=Path)
    parser.add_argument('--controls', type=Path)
    parser.add_argument('--occupied-cap', type=int, choices=(4, 32), default=32)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        if not args.source or not args.controls:
            parser.error('prepare requires --source and --controls')
        prepare(args.source.resolve(), args.controls.resolve(), args.occupied_cap, args.out.resolve())
    else:
        run(args.out.resolve())


if __name__ == '__main__':
    main()
