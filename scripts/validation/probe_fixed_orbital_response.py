#!/usr/bin/env python
"""Test one atom's response using exact saved SCF orbitals and H1.

This is an isolated diagnostic. It runs no SCF, optimization, full molecular
Hessian or free-energy calculation. All saved intermediates remain unqualified.
"""
import argparse
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import inspect
import json
import os
from pathlib import Path
import signal
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / 'scripts/pipeline/residual_response_solver.py'


def controls(root):
    import numpy as np
    proof = {}
    profile_path=ROOT/'environments/cutensor22/runtime.json'
    profile=json.loads(profile_path.read_text())
    for mult in (1, 2):
        p = root / 'controls' / str(mult) / 'status.json'
        row = json.loads(p.read_text())
        if row['status'] != 'full_original_compact_equivalence_passed' or not row.get('response_solver_checks_passed'):
            raise ValueError('Both complete residual-solver spin controls required')
        api.verify_sources(row['source_sha256'])
        api.verify_sources(row['evidence_sha256'])
        runtime = root / 'runtime' / str(mult) / 'status.json'
        rt = json.loads(runtime.read_text())
        if (rt['status'] != 'task_returned_under_checked_cutensor_runtime'
                or rt['engine'] != 'cutensor' or rt['pid'] != row['pid']
                or rt['versions'] != profile['versions'] or rt['library_sha256'] != profile['library_sha256']
                or rt['profile_sha256'] != api.digest(profile_path)
                or rt['source_sha256'] != api.digest(ROOT/'scripts/fleet/tensor_runtime_worker.py')):
            raise ValueError('Exact completed runtime must match the full control')
        directory=p.parent
        for name in ('gradient.npy','hessian.npy'):
            difference=float(np.max(abs(np.load(directory/'original'/name)-np.load(directory/'compact'/name))))
            if not np.isfinite(difference) or difference>1e-6:
                raise ValueError('Independent full control-array comparison failed')
        for name in ('E_smd_eV','E_gas_eV','geometry_level_E_Ha'):
            if abs(row['arms']['original'][name]-row['arms']['compact'][name])>1e-6:
                raise ValueError('Independent full control-energy comparison failed')
        guard=root/'guards'/str(mult)
        owner=json.loads((guard/'owner.json').read_text());residency=json.loads((guard/'state.json').read_text())
        if (owner['pid'] != row['pid'] or not owner['guard_stopped']
                or not residency['own_residency_observed'] or residency['foreign_resident_PIDs']
                or (guard/'interruption.json').exists()):
            raise ValueError('Full control did not retain exclusive residency')
        proof.update(row['source_sha256']); proof.update(row['evidence_sha256'])
        proof.update({str(f): api.digest(f) for f in (p,runtime,guard/'owner.json',guard/'state.json')})
    return proof


def first_atom_source(source):
    changes = (
        ('mo1s = np.zeros(h1mo.shape)', 'mo1s = np.zeros((1,) + h1mo.shape[1:])'),
        ('e1s = np.zeros((natm, 3, nocc, nocc))', 'e1s = np.zeros((1, 3, nocc, nocc))'),
        ('for i0, i1 in lib.prange(0, natm, blksize):', 'for i0, i1 in [(0, 1)]:'))
    for old, new in changes:
        if source.count(old) != 1:
            raise ValueError('Installed first-atom response template changed')
        source = source.replace(old, new, 1)
    return source


def prepare(source, control_root, out):
    plan = api.read_plan(source)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    status = source / 'streamed_native/status.json'
    row = json.loads(status.read_text())
    owner = api.load_module('fixed_response_source_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    if ('Krylov subspace iterations diverge' not in row.get('error', '')
            or not row.get('finished') or owner.owner_alive(row['owner'])):
        raise ValueError('A preserved, stopped Krylov failure is required')
    raw = Path(plan['raw_root']) / 'streamed_execution'
    paths = [source / 'plan.json', status, source / 'SCF_checkpoints/SCF_0000/SCF.h5',
        source / 'SCF_checkpoints/SCF_0000/settings.json', raw / 'rks_make_h1_state_identity.json',
        raw / 'rks_make_h1_0.npy', ROOT / 'environments/cutensor22/runtime.json']
    inputs = {str(p): api.digest(p) for p in paths}
    out.mkdir(parents=True, exist_ok=False)
    proof = controls(control_root)
    producers = {**proof, **{str(p): api.digest(p) for p in (Path(__file__).resolve(), ADAPTER)}}
    atomic_write_json(out / 'plan.json', dict(input_sha256=inputs, producer_sha256=producers,
        source=str(source), checkpoint=str(paths[2]), settings=str(paths[3]), H1_identity=str(paths[4]),
        H1=str(paths[5]), runtime_profile=str(paths[6]), execution_driver=str(Path(__file__).resolve()),
        raw_root=str(ROOT / 'calcs/validation' / out.parent.name / out.name),
        maximum_host_RSS_GiB=80, minimum_GPU_MiB=32768, maximum_GPU_MiB=31744,
        atom_index=0, physical_residual_maximum_Ha=1e-8, SCF_recomputed=False,
        full_Hessian_qualified=False, native_G_qualified=False,
        purpose='First-atom fixed-state CPHF diagnostic, original physical equations and exact H1/MOs'))


def run(campaign):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    raw = Path(plan['raw_root']); raw.mkdir(parents=True, exist_ok=False)
    out = campaign / 'response'; out.mkdir(exist_ok=False)
    row = dict(status='isolated_fixed_response_probe_running', started=time.time(),
        owner=dict(host=os.uname().nodename, pid=os.getpid(), process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        SCF_recomputed=False, full_Hessian_qualified=False, native_G_qualified=False,
        atom_index=0, CPU_affinity=sorted(os.sched_getaffinity(0)))
    atomic_write_json(out / 'status.json', row)
    def stop(signum, frame):
        raise KeyboardInterrupt('Owned response probe stop; intermediate records preserved')
    previous = signal.signal(signal.SIGTERM, stop)
    memory = None
    try:
        import cupy as cp
        import numpy as np
        from pyscf.scf import chkfile
        from gpu4pyscf.hessian import rhf
        from redox.qm import dft
        mol, state = chkfile.load_scf(plan['checkpoint'])
        identity = json.loads(Path(plan['H1_identity']).read_text())
        if (mol.natm, mol.charge, mol.spin) != (139, 0, 0):
            raise ValueError('Unexpected saved state')
        for name in ('mo_coeff','mo_occ','mo_energy'):
            array = np.asarray(state[name])
            if (list(array.shape) != identity[name]['shape'] or str(array.dtype) != identity[name]['dtype']
                    or hashlib.sha256(array.tobytes()).hexdigest() != identity[name]['sha256']):
                raise ValueError('Saved H1 and restored molecular orbitals differ: '+name)
        if (not identity['SCF_converged'] or abs(float(state['e_tot'])-identity['SCF_E_Ha']) > 1e-12
                or not np.allclose(mol.atom_coords(),identity['coordinates_Bohr'],rtol=0,atol=1e-12)):
            raise ValueError('Saved derivative and exact baseline state differ')
        coeff, occ, energies = (np.asarray(state[n]) for n in ('mo_coeff','mo_occ','mo_energy'))
        overlap = mol.intor_symmetric('int1e_ovlp')
        if np.max(abs(coeff.T@overlap@coeff-np.eye(coeff.shape[1]))) > 1e-8:
            raise ValueError('Saved MOs are not orthonormal')
        H1 = np.load(plan['H1'],mmap_mode='r')
        if H1.shape != (139,3,coeff.shape[1],int((occ>0).sum())):
            raise ValueError('Saved H1 dimensions differ')
        adapter = api.load_module('fixed_state_response_checked_adapter', ADAPTER)
        observer = api.load_module('fixed_state_response_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        with observer.MemoryMonitor() as memory, (raw/'calculation.log').open('x') as log, redirect_stdout(log), redirect_stderr(log), adapter.compact_smd_execution(raw/'execution',28):
            mf = dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
            mf.mo_coeff, mf.mo_occ, mf.mo_energy = (cp.asarray(v) for v in (coeff,occ,energies))
            mf.e_tot=float(state['e_tot']); mf.converged=True
            mf.conv_tol=1e-12; mf.conv_tol_grad=1e-8; mf.grids.level=3
            mf.grids.build(); mf.cphf_grids=mf.grids.copy()
            mf.with_df.build()
            mf.with_solvent._get_vind(mf.make_rdm1())
            body = first_atom_source(inspect.getsource(rhf.solve_mo1))
            snapshot=raw/'first_atom_solve_mo1.py';snapshot.write_text(body)
            namespace=dict(rhf.solve_mo1.__globals__);local={}
            exec(compile(body,str(snapshot),'exec'),namespace,local)
            result, e1 = local['solve_mo1'](mf,mf.mo_energy,mf.mo_coeff,mf.mo_occ,H1,max_cycle=100)
            # Separate AO overlap derivative and original unpreconditioned
            # CPHF equation: (ea-ei)U + H1 - S1*ei + V[U] = 0.
            oi, vi = occ>0, occ==0
            ip=-mol.intor('int1e_ipovlp',comp=3)
            s=np.zeros_like(ip);a,b=mol.aoslice_by_atom()[0,2:]
            s[:,a:b]+=ip[:,a:b];s[:,:,a:b]+=ip[:,a:b].transpose(0,2,1)
            s1=np.einsum('xpq,pi,qj->xij',s,coeff,coeff[:,oi],optimize=True)
            induced=cp.asnumpy(rhf.gen_vind(mf,mf.mo_coeff,mf.mo_occ)(cp.asarray(result))).reshape(result.shape)
            hs=np.array(H1[:1])-s1[None]*energies[oi]
            residual=result[:,:,vi]*(energies[vi,None]-energies[oi])+hs[:,:,vi]+induced[:,:,vi]
            error=float(np.max(abs(residual)))
            np.save(raw/'first_atom_MO_response.npy',result)
            np.save(raw/'first_atom_orbital_energy_response.npy',e1)
            atomic_write_json(raw/'physical_residual.json',dict(maximum_Ha=error,
                acceptance_Ha=plan['physical_residual_maximum_Ha'],passed=bool(np.isfinite(error) and error<=plan['physical_residual_maximum_Ha']),
                all_atoms_solved=False,first_atom_only=True))
            if not np.isfinite(error) or error>plan['physical_residual_maximum_Ha']:
                raise ValueError('Independent unpreconditioned response residual failed')
        if not 0<memory.peak<=plan['maximum_GPU_MiB']:
            raise ValueError('Fixed-state response device ceiling exceeded')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        row.update(status='first_atom_response_equations_checked_full_native_QC_pending',maximum_physical_residual_Ha=error,
            output_sha256={str(p):api.digest(p) for p in raw.rglob('*') if p.is_file()})
    except BaseException as error:
        row.update(status='isolated_fixed_response_probe_failed_preserved',error=repr(error));raise
    finally:
        signal.signal(signal.SIGTERM,previous)
        row.update(finished=time.time(),observed_peak_GPU_MiB=memory.peak if memory else None)
        atomic_write_json(out/'status.json',row)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run'))
    parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--controls',type=Path)
    args=parser.parse_args()
    if args.action=='prepare':
        if not args.source or not args.controls:parser.error('prepare requires source and controls')
        prepare(args.source.resolve(),args.controls.resolve(),args.campaign.resolve())
    else:run(args.campaign.resolve())


if __name__=='__main__':
    main()
