#!/usr/bin/env python
"""Evaluate a specified same-basin conformer pair with native-SMD response.

Explicit geometry paths allow paired-conformer controls and new registered
motifs without relying on a dated discovery table. The output is an electronic
solution relaxation diagnostic, never a polymer constraint or experimental
accuracy assertion. Wrap execution in the fleet reservation/preflight workflow.
"""
from __future__ import annotations
import argparse
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT,HARTREE_EV
from redox.core.protocol import ACTIVE_SP
from redox.properties.solution_reorg import directional_relaxation,summarize_relaxations
from redox.qm.dft import _ecp_for,OPT_XC,OPT_DISP,OPT_BASIS,OPT_BASIS_ANION,_basis_for
from redox.qm.nonequilibrium_pcm import host_array
from redox.qm.nonequilibrium_smd import smd_object,smd_dielectrics,attach_frozen_slow_smd
from redox.qm.native_thermal import load_native_thermal,cache_path

spec=importlib.util.spec_from_file_location('reference_solver',ROOT/'scripts/validation/physical_reorg.py')
solver=importlib.util.module_from_spec(spec);spec.loader.exec_module(solver)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mol_at(state,geometry,basis):
    atoms=read(str(geometry));symbols=atoms.get_chemical_symbols()
    return gto.M(atom=list(zip(symbols,atoms.positions.tolist())),basis=basis,ecp=_ecp_for(symbols,basis),
                 charge=state['charge'],spin=state['mult']-1,verbose=0,max_memory=32000)


def atomic(path,rec):
    tmp=path.with_suffix(f'.{os.getpid()}.tmp');tmp.write_text(json.dumps(rec,indent=2,allow_nan=False));tmp.replace(path)


def release_solver_memory(backend):
    # Frozen-response closures can retain solvent arrays in reference cycles.
    # Releasing each completed SCF avoids accumulating them across vertical legs.
    gc.collect()
    if backend=='gpu':
        from redox.qm.dft import _free_gpu_pool
        _free_gpu_pool()


def evaluate(pair,out,backend,cache_pool):
    out.mkdir(parents=True,exist_ok=True)
    states=pair['states'];constants=smd_dielectrics()
    if len(states)!=2 or states[0]['charge']-states[1]['charge']!=1:
        raise ValueError('pair must specify oxidized then reduced adjacent electronic states')
    if any(s['charge']!=int(s['charge']) or s['mult']<1 for s in states):raise ValueError('invalid charge/multiplicity')
    if (read(states[0]['geometry']).get_chemical_symbols()!=read(states[1]['geometry']).get_chemical_symbols()):
        raise ValueError('Paired geometries have different atom ordering/composition')
    for state in states:
        if state['name'] in pair.get('require_native_minimum_states',[]):
            sd=Path(state.get('state_dir',Path(state['geometry']).parent))
            if (not (sd/'opt.xyz').exists() or sha(sd/'opt.xyz')!=sha(state['geometry'])
                or load_native_thermal(sd,state['charge'],state['mult']) is None):
                raise RuntimeError(f'Required native minimum is not verified: {state["name"]}')
    identity=dict(pair=pair,input_sha256={s['geometry']:sha(s['geometry']) for s in states},
        method=dict(ACTIVE_SP,solvent_model='native_SMD_frozen_slow',**constants,sasa_ng=590),
        source_sha256={str(p):sha(p) for p in [Path(__file__),ROOT/'src/redox/qm/nonequilibrium_smd.py',
            ROOT/'scripts/validation/physical_reorg.py']})
    atomic(out/'identity.json',identity)
    # Existing complete records can be used only on the same geometry, species,
    # uniform production protocol and native SMD static/optical constants.
    existing=[]
    if cache_pool:
        for p in Path(cache_pool).glob('*.json'):
            r=json.loads(p.read_text());k=r.get('key',{});method=k.get('method',{})
            if (r.get('converged') is True and np.isfinite([r.get('e_Ha',np.nan),r.get('e_eV',np.nan)]).all()
                and method==identity['method']):
                existing.append((p,r))
    def scf(state,at,reference=None,dm0=None):
        geometry=Path(at['geometry']);content=geometry.read_bytes()
        if sha(geometry)!=identity['input_sha256'][str(geometry)]:raise RuntimeError('Pair input changed before SCF')
        for path,r in existing:
            k=r['key'];species=k['identity'];g=k['geometry'];ref=k['reference']
            if (species['charge']!=state['charge'] or species['mult']!=state['mult']
                or g['geom_sha1']!=hashlib.sha1(content).hexdigest()[:16]):continue
            if (reference is None)!=(ref is None):continue
            if dm0 is not None and k.get('refdm_sha256')!=hashlib.sha256(host_array(dm0).tobytes()).hexdigest():continue
            return dict(r,reused_from=str(path)),np.load(path.with_suffix('.npz'))['dm']
        key=dict(state=state,geometry_sha256=sha(geometry),reference=reference,
                 reference_density_sha256=sha_bytes(dm0) if dm0 is not None else None,identity=identity)
        digest=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest()[:20]
        path=out/f'scf_{digest}.json'
        if path.exists():
            r=json.loads(path.read_text())
            if r.get('key')==key and r.get('converged') is True:return r,np.load(path.with_suffix('.npz'))['dm']
            raise RuntimeError(f'Existing unverified cache: {path}')
        mol=mol_at(state,geometry,ACTIVE_SP['basis'])
        mf=solver.make_mf(mol,backend,ACTIVE_SP['xc'],ACTIVE_SP['nlc'])
        mf=mf.SMD(smd_object(mol,backend=backend)) if dm0 is None else attach_frozen_slow_smd(mf,dm0,backend=backend)
        print(f'PAIR_SCF_START {state["name"]} at{at["name"]} reference={reference}',flush=True)
        best,guesses=solver.solve(mf);e,dm,orbitals,spin=best
        if sha(geometry)!=identity['input_sha256'][str(geometry)]:raise RuntimeError('Input changed during SCF')
        total=host_array(dm).sum(axis=0) if dm.ndim==3 else host_array(dm)
        trace=float(np.einsum('ij,ji',total,mol.intor('int1e_ovlp')))
        if abs(trace-mol.nelectron)>1e-6:raise RuntimeError('Electron count failed')
        r=dict(key=key,e_Ha=e,e_eV=e*HARTREE_EV,converged=True,scf_guesses=guesses,
            cds_Ha=float(mf.with_solvent.get_cds()),electron_trace=trace,s2=spin.get('s2'),mo_energy_Ha=orbitals)
        density=host_array(dm)
        np.savez_compressed(path.with_suffix('.npz'),dm=density);atomic(path,r)
        del mf,best,dm
        release_solver_memory(backend)
        return r,density
    endpoint=[];electronic=[]
    for state in states:
        mol=mol_at(state,state['geometry'],_basis_for(state['charge'],OPT_BASIS,OPT_BASIS_ANION))
        mf=solver.make_mf(mol,backend,OPT_XC);mf.disp=OPT_DISP;mf.grids.level=3
        mf=mf.SMD(smd_object(mol,backend=backend));best,_=solver.solve(mf)
        e=float(mf.kernel(dm0=best[1]));grad=mf.nuc_grad_method();grad.grid_response=True
        g=host_array(grad.kernel())
        if not mf.converged or abs(e-best[0])>1e-7:raise RuntimeError('Endpoint SCF restoration failed')
        endpoint.append(dict(state=state['name'],minimum_verified=False,conformer_policy=pair['conformer_policy'],
            gradient_max_atomic_norm_Ha_Bohr=float(np.linalg.norm(g,axis=1).max()),
            stationary_gradient_pass=bool(np.linalg.norm(g,axis=1).max()<4.5e-4)))
        if not endpoint[-1]['stationary_gradient_pass']:
            atomic(out/'endpoint_gradient_failure.json',dict(identity=identity,endpoint_evidence=endpoint))
            raise RuntimeError('Paired geometry is nonstationary with conservative full-grid response; reoptimize before lambda')
        sd=Path(state.get('state_dir',Path(state['geometry']).parent))
        if (sd/'opt.xyz').exists() and sha(sd/'opt.xyz')==identity['input_sha256'][state['geometry']]:
            thermal=load_native_thermal(sd,state['charge'],state['mult'])
            if thermal:
                native,_=cache_path(sd,state['charge'],state['mult'])
                endpoint[-1].update(minimum_verified=True,minimum_source=str(native),
                    minimum_source_sha256=sha(native),minimum_frequency_cm=min(thermal['frequencies_cm']),
                    minimum_energy_surface=thermal['identity']['protocol'],active_SP_stationarity_verified=False)
        del mf,grad,best
        release_solver_memory(backend)
    def electronic_evidence(state,at,response,record):
        # Audit all equilibrium and vertical legs, including inherited records.
        mol=mol_at(state,at['geometry'],ACTIVE_SP['basis']);energies=np.asarray(record['mo_energy_Ha'])
        if energies.ndim==2:
            na,nb=mol.nelec;homo=float(max(energies[0,na-1],energies[1,nb-1]))
        else:homo=float(energies[mol.nelectron//2-1])
        flags=[]
        if state['charge']<0 and homo>0:flags.append('positive_solution_homo')
        spin=(state['mult']-1)/2;s2=record.get('s2')
        contamination=None if s2 is None else s2-spin*(spin+1)
        if contamination is not None and contamination>.1:flags.append('spin_contamination')
        electronic.append(dict(state=state['name'],at=at['name'],solvent_response=response,
            solution_homo_eV=homo*HARTREE_EV,s2=s2,spin_contamination=contamination,flags=flags))
    finals={s['name']:scf(s,s) for s in states};directions=[]
    for s in states:electronic_evidence(s,s,'equilibrium',finals[s['name']][0])
    for initial,final in (states,states[::-1]):
        ieq,dm0=finals[initial['name']];feq,_=finals[final['name']]
        cross,_=scf(final,initial);neq,_=scf(final,initial,initial,dm0)
        electronic_evidence(final,initial,'equilibrium',cross)
        electronic_evidence(final,initial,'frozen_slow',neq)
        row=directional_relaxation(ieq['e_eV'],cross['e_eV'],neq['e_eV'],feq['e_eV'])
        row.update(initial=initial['name'],final=final['name']);directions.append(row)
    for s in states:
        if sha(s['geometry'])!=identity['input_sha256'][s['geometry']]:raise RuntimeError('Pair input changed')
    summary=summarize_relaxations(*directions,endpoint_evidence=endpoint,method_identity=identity['method'])
    summary['qualification_convention']='composite active-SP energies at geometry-level solution minima'
    summary['active_SP_stationary_endpoints_proven']=False
    flags=sorted({flag for e in electronic for flag in e['flags']})
    if flags:summary['physical_flags'].extend(flags);summary['endpoint_model_verified']=False
    atomic(out/'result.json',dict(identity=identity,directional=directions,endpoint_evidence=endpoint,
        electronic_state_evidence=electronic,summary=summary,
        validation_status='paired_conformer_sensitivity',polymer_constraint_proven=False))


def sha_bytes(array):return hashlib.sha256(host_array(array).tobytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--pair',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--cache-pool',type=Path)
    p.add_argument('--backend',choices=['cpu','gpu'],default='gpu');p.add_argument('--dry-run',action='store_true')
    a=p.parse_args();pair=json.loads(a.pair.read_text())
    if a.dry_run:
        print(json.dumps(dict(pair=pair,input_sha256={s['geometry']:sha(s['geometry']) for s in pair['states']})));return
    if a.backend=='gpu':
        from redox.workflows.resources import verify_gate_owner
        verify_gate_owner()
    def terminate(sig,frame):raise KeyboardInterrupt(f'termination signal{sig}')
    signal.signal(signal.SIGTERM,terminate)
    a.out.mkdir(parents=True,exist_ok=True);claim=a.out/'claim';claim.mkdir()
    atomic(claim/'owner.json',dict(host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().split()[21]))
    try:evaluate(pair,a.out,a.backend,a.cache_pool)
    finally:(claim/'owner.json').unlink();claim.rmdir()


if __name__=='__main__':main()
