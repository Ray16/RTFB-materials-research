#!/usr/bin/env python
"""Prepare/run a bounded native-SMD frozen-slow cohort in validation only.

Preparing a queue uses current geometry, charge/spin and active energy identities.
Running MUST be externally wrapped in the shared GPU reservation and resource
preflight. This module never claims a GPU, changes production or clears flags.
Three controls precede the flagged-discovery extension. Endpoint forces and
minimum/conformer evidence remain separate quality fields.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import time

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT,DFT,HARTREE_EV,read_result
from redox.core.protocol import ACTIVE_SP,geom_sha1,record_path
from redox.properties.solution_reorg import directional_relaxation,summarize_relaxations
from redox.qm.dft import _ecp_for,OPT_XC,OPT_DISP,OPT_BASIS,OPT_BASIS_ANION,_basis_for
from redox.qm.nonequilibrium_pcm import host_array
from redox.qm.nonequilibrium_smd import (smd_object,smd_dielectrics,attach_frozen_slow_smd,
                                        independent_hcore_smd_reference)

OUT=ROOT/'results/validation/physical_reorg_2026-10-11'
CONTROL=[('anthraquinone_parent','neu','red1'),('anthraquinone_parent','red1','red2'),
         ('anthraquinone','neu','red1')]
spec=importlib.util.spec_from_file_location('validated_noneq_solver',ROOT/'scripts/validation/physical_reorg.py')
solver=importlib.util.module_from_spec(spec);spec.loader.exec_module(solver)


def atomic(path,record):
    p=path.with_suffix(f'.{os.getpid()}.tmp');p.write_text(json.dumps(record,indent=2,allow_nan=False));p.replace(path)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def state_identity(gid,state):
    p=DFT/gid/state/'result.json';geom=p.with_name('opt.xyz')
    info=json.loads(p.read_text());q,m=int(info['charge']),int(info['mult'])
    energy=record_path(p.parent,q,m)
    active=read_result(gid,state)
    if active is None or active.get('e_smd_eV') is None:
        raise RuntimeError(f'No validated active solution energy for {gid}/{state}')
    return dict(id=gid,state=state,charge=q,mult=m,geom_sha1=geom_sha1(geom),
        raw_result_sha256=sha(p),active_record=str(energy),active_record_sha256=sha(energy),
        e_smd_eV=active['e_smd_eV'],thermal_qc=('missing' if active.get('g_thermal_eV') is None
             else 'imag_modes' if active.get('n_imag',0)>0 else 'ok'),
        thermal_qc_is_not_DFT_minimum_proof=True)


def prepare(table,with_matched_benchmark=False):
    rows=list(csv.DictReader(Path(table).open()))
    flagged=[(r['id'],*r['couple'].split('->')) for r in rows if r['id'].startswith('d3_') and r['flag']]
    cases=[]
    benchmark=[('bm_2_6_dimethyl_1_4_benzoquinone','neu','red1')] if with_matched_benchmark else []
    for gid,o,r in CONTROL+benchmark+flagged:
        cases.append(dict(id=gid,oxidized=o,reduced=r,kind='control' if (gid,o,r) in CONTROL
                         else 'matched_experimental_benchmark' if (gid,o,r) in benchmark
                         else 'flagged_discovery',states=[state_identity(gid,o),state_identity(gid,r)],
                         original_flag=next((x['flag'] for x in rows if x['id']==gid
                                            and x['couple']==f'{o}->{r}'),'')))
    rec=dict(created=time.time(),status='prepared_not_submitted',cases=cases,
        method=dict(ACTIVE_SP,solvent_model='native_SMD_frozen_slow',**smd_dielectrics(),
                    sasa_ng=590,cds_assumption='geometry_only_native_functional'),
        source_table=str(table),source_table_sha256=sha(table),
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),
            ROOT/'src/redox/qm/nonequilibrium_smd.py',ROOT/'src/redox/qm/nonequilibrium_pcm.py',
            ROOT/'src/redox/properties/solution_reorg.py',ROOT/'scripts/validation/physical_reorg.py']},
        cpu_slots_per_worker=2,max_gpu_workers=1,
        dependency='Run only after existing thermal/spin 2GPU campaign releases a slot and root authorizes launch',
        endpoint_contract=dict(minimum_evidence='same-method conservative-gradient Hessian or equivalent',
            conformer_policy='original minima retained; no flags cleared until same-basin/ensemble validation',
            production_promotion=False),
        original_gas_inner_values_preserved=True)
    OUT.mkdir(parents=True,exist_ok=True)
    atomic(OUT/'native_smd_queue.json',rec)
    print(json.dumps(dict(cases=len(cases),controls=3,matched_benchmarks=len(benchmark),
                          flagged=len(flagged),status=rec['status'])))


def molecule(identity,geom_identity):
    atoms=read(str(DFT/geom_identity['id']/geom_identity['state']/'opt.xyz'))
    sym=atoms.get_chemical_symbols()
    return gto.M(atom=list(zip(sym,atoms.positions.tolist())),basis=ACTIVE_SP['basis'],
        ecp=_ecp_for(sym,ACTIVE_SP['basis']),charge=identity['charge'],spin=identity['mult']-1,
        verbose=3,max_memory=32000)


def independent_reference(backend):
    """Gate new backend against an independent native-SMD SCF assembly."""
    atom='C 0 0 0; O 0 0 1.21; H .94 0 -.58; H -.94 0 -.58'
    neutral=gto.M(atom=atom,basis='6-31+g*',verbose=0)
    eq=solver.make_mf(neutral,backend,'pbe').SMD(smd_object(neutral,backend=backend,order=11))
    best,_=solver.solve(eq);dm0=host_array(best[1])
    ion=gto.M(atom=atom,basis='6-31+g*',charge=1,spin=1,verbose=0)
    neq=attach_frozen_slow_smd(solver.make_mf(ion,backend,'pbe'),dm0,backend=backend,order=11)
    a,_=solver.solve(neq)
    independent,constant=independent_hcore_smd_reference(solver.make_mf(ion,'cpu','pbe'),dm0,order=11)
    b,_=solver.solve(independent)
    residual=(a[0]-b[0]-constant)*HARTREE_EV
    rec=dict(status='passed' if abs(residual)<1e-5 else 'failed',backend=backend,
             residual_eV=residual,method='PBE/6-31+G* native SMD',
             wrapper_Ha=a[0],independent_hcore_Ha=b[0]+constant)
    atomic(OUT/f'native_smd_{backend}_reference.json',rec)
    if rec['status']!='passed':raise RuntimeError(f'Native SMD backend gate failed: {rec}')


def endpoint_gradient(identity,backend):
    """Conservative optimization-level force evidence; a force is not a Hessian."""
    path=OUT/'endpoint_forces'/f'{identity["id"]}_{identity["state"]}_{identity["geom_sha1"]}.json'
    path.parent.mkdir(exist_ok=True)
    key=dict(state=identity,method=dict(xc=OPT_XC,disp=OPT_DISP,
        basis=_basis_for(identity['charge'],OPT_BASIS,OPT_BASIS_ANION),
        solvent='native_SMD_MeCN',grid_response=True,grid_level=3,scf_conv_tol=1e-9),
        source_sha256=sha(Path(__file__)))
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('key')==key:return old
        raise RuntimeError(f'Existing force record differs: {path}')
    atoms=read(str(DFT/identity['id']/identity['state']/'opt.xyz'));symbols=atoms.get_chemical_symbols()
    basis=key['method']['basis']
    mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),charge=identity['charge'],
        spin=identity['mult']-1,basis=basis,ecp=_ecp_for(symbols,basis),verbose=0)
    mf=solver.make_mf(mol,backend,OPT_XC);mf.disp=OPT_DISP;mf.grids.level=3
    mf=mf.SMD(smd_object(mol,backend=backend))
    best,_=solver.solve(mf)
    # Restore the exact lowest converged density by a deterministic SCF restart;
    # gradients cannot be evaluated on the last higher independent guess.
    energy=float(mf.kernel(dm0=best[1]))
    if not mf.converged or abs(energy-best[0])>1e-7:
        raise RuntimeError('Could not restore chosen endpoint SCF for force')
    gradient=mf.nuc_grad_method();gradient.grid_response=True
    g=host_array(gradient.kernel())
    if not np.isfinite(g).all():raise RuntimeError('Nonfinite endpoint gradient')
    rec=dict(key=key,e_Ha=energy,gradient_Ha_Bohr=g.tolist(),
        gradient_max_Ha_Bohr=float(abs(g).max()),gradient_rms_Ha_Bohr=float(np.sqrt(np.mean(g*g))),
        gradient_max_atomic_norm_Ha_Bohr=float(np.linalg.norm(g,axis=1).max()),
        gradient_rms_atomic_norm_Ha_Bohr=float(np.sqrt(np.mean(np.sum(g*g,axis=1)))),
        gradient_max_threshold_Ha_Bohr=4.5e-4,
        stationary_gradient_pass=bool(np.linalg.norm(g,axis=1).max()<4.5e-4),
        minimum_verified=False,reason='Hessian minimum evidence still required')
    atomic(path,rec);return rec


def scf(identity,geom_identity,backend,reference=None,refdm=None):
    key=dict(identity=identity,geometry=geom_identity,reference=reference,
        method=dict(ACTIVE_SP,solvent_model='native_SMD_frozen_slow',**smd_dielectrics(),sasa_ng=590),
        source={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),
            ROOT/'src/redox/qm/nonequilibrium_smd.py',ROOT/'scripts/validation/physical_reorg.py']})
    if refdm is not None:key['refdm_sha256']=hashlib.sha256(host_array(refdm).tobytes()).hexdigest()
    digest=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest()[:20]
    p=OUT/'native_smd_scf'/f'{identity["id"]}_{identity["state"]}_at_{geom_identity["state"]}_{digest}.json'
    p.parent.mkdir(exist_ok=True)
    if p.exists():
        old=json.loads(p.read_text())
        if old.get('key')==key and old.get('converged') is True:
            return old,np.load(p.with_suffix('.npz'))['dm']
        raise RuntimeError(f'Existing unverified cache: {p}')
    mol=molecule(identity,geom_identity)
    mf=solver.make_mf(mol,backend,ACTIVE_SP['xc'],ACTIVE_SP['nlc'])
    mf=(mf.SMD(smd_object(mol,backend=backend)) if refdm is None
        else attach_frozen_slow_smd(mf,refdm,backend=backend))
    started=time.time();print(f'SMD_SCF_START {p.name}',flush=True)
    best,guesses=solver.solve(mf);energy,dm,orbitals,spin=best
    # Explicit CDS bookkeeping; native attachment includes it exactly once.
    cds=float(mf.with_solvent.get_cds())
    rec=dict(key=key,e_Ha=energy,e_eV=energy*HARTREE_EV,converged=True,
        cds_Ha=cds,scf_guesses=guesses,s2=spin.get('s2'),mo_energy_Ha=orbitals,
        elapsed_seconds=time.time()-started)
    density=host_array(dm);total=density.sum(axis=0) if density.ndim==3 else density
    electron_trace=float(np.einsum('ij,ji',total,mol.intor('int1e_ovlp')))
    rec.update(electron_trace=electron_trace,expected_electrons=mol.nelectron,
               electronic_state_flags=[])
    if abs(electron_trace-mol.nelectron)>1e-6 or not np.isfinite(density).all():
        raise RuntimeError('Density trace/finiteness audit failed')
    energies=np.asarray(orbitals)
    if energies.ndim==2:
        na,nb=mol.nelec
        homo=float(max(energies[0,na-1],energies[1,nb-1]))
    else:homo=float(energies[mol.nelectron//2-1])
    rec['solution_homo_eV']=homo*HARTREE_EV
    if identity['charge']<0 and homo>0:
        rec['electronic_state_flags'].append('positive_solution_homo')
    if rec['s2'] is not None:
        spin_quantum=(identity['mult']-1)/2
        rec['spin_contamination']=rec['s2']-spin_quantum*(spin_quantum+1)
        if rec['spin_contamination']>.1:
            rec['electronic_state_flags'].append('spin_contamination')
    if refdm is None and identity['state']==geom_identity['state']:
        delta=rec['e_eV']-identity['e_smd_eV'];rec['native_production_replay_delta_eV']=delta
        if abs(delta)>1e-4:
            raise RuntimeError(f'Native SMD active-energy replay disagrees by{delta}eV')
    np.savez_compressed(p.with_suffix('.npz'),dm=host_array(dm));atomic(p,rec)
    print(f'SMD_SCF_DONE {p.name}',flush=True)
    return rec,host_array(dm)


def run(plan,backend,limit):
    expected=plan['code_sha256']
    for file,h in expected.items():
        if sha(ROOT/file)!=h:raise RuntimeError(f'Prepared queue producer changed: {file}; reprepare queue')
    owner=dict(host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().split()[21],
        cpus=sorted(os.sched_getaffinity(0)),created=time.time())
    atomic(OUT/'native_smd_worker.json',owner)
    def termination(sig,frame):raise KeyboardInterrupt(f'termination signal{sig}')
    signal.signal(signal.SIGTERM,termination)
    independent_reference(backend)
    for case in plan['cases'][:limit]:
        gid=case['id'];o,r=case['oxidized'],case['reduced'];identities=case['states']
        for item in identities:
            if state_identity(gid,item['state'])!=item:
                raise RuntimeError(f'State changed after queue preparation: {gid}/{item["state"]}')
        p=OUT/f'native_smd_{gid}_{o}_{r}.json'
        if p.exists():
            old=json.loads(p.read_text())
            if old.get('case')==case and old.get('code_sha256')==expected:continue
            raise RuntimeError(f'Changed identity must not replace old record: {p}')
        claim=p.with_suffix('.claim');claim.mkdir();atomic(claim/'owner.json',owner)
        try:
            force_evidence=[endpoint_gradient(i,backend) for i in identities]
            finals={}
            for item in identities:finals[item['state']]=scf(item,item,backend)
            directions=[]
            electronic_evidence=[dict(species=s,geometry=s,solvent_response='equilibrium',
                flags=e[0]['electronic_state_flags']) for s,e in finals.items()]
            for initial,final in (identities,identities[::-1]):
                ieq,dm0=finals[initial['state']]
                feq,_=finals[final['state']]
                cross,_=scf(final,initial,backend)
                vertical,_=scf(final,initial,backend,initial,dm0)
                for response,record in [('equilibrium',cross),('frozen_slow',vertical)]:
                    electronic_evidence.append(dict(species=final['state'],geometry=initial['state'],
                        solvent_response=response,flags=record['electronic_state_flags']))
                row=directional_relaxation(ieq['e_eV'],cross['e_eV'],vertical['e_eV'],feq['e_eV'])
                row.update(initial=initial['state'],final=final['state'],
                    cds_relaxation_eV=(vertical['cds_Ha']-feq['cds_Ha'])*HARTREE_EV)
                directions.append(row)
            evidence=[dict(state=i['state'],minimum_verified=False,conformer_policy='different_minima',
                thermal_qc=i['thermal_qc'],stationary_gradient_pass=f['stationary_gradient_pass'],
                gradient_max_Ha_Bohr=f['gradient_max_atomic_norm_Ha_Bohr']) for i,f in zip(identities,force_evidence)]
            summary=summarize_relaxations(*directions,endpoint_evidence=evidence,
                                         method_identity=plan['method'])
            state_flags=sorted({flag for e in electronic_evidence for flag in e['flags']})
            if state_flags:
                summary['physical_flags'].extend(state_flags);summary['endpoint_model_verified']=False
            atomic(p,dict(case=case,code_sha256=expected,method=plan['method'],directional=directions,
                          summary=summary,endpoint_evidence=evidence,electronic_state_evidence=electronic_evidence,
                          validation_status='numerically_checked_solution_model_estimate'))
        except Exception as exc:
            atomic(p.with_name(p.stem+'_failure.json'),dict(case=case,error=repr(exc),owner=owner))
            raise
        finally:(claim/'owner.json').unlink();claim.rmdir()


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--prepare',action='store_true')
    p.add_argument('--with-matched-benchmark',action='store_true',
                   help='Place existing2,6-dimethylBQ experimental anchor after three controls')
    p.add_argument('--table',default=str(ROOT/'results/review_2026-10-11_physical_checkpoint/reorganization.csv'))
    p.add_argument('--execute',action='store_true')
    p.add_argument('--backend',choices=['cpu','gpu'],default='gpu')
    p.add_argument('--limit',type=int,default=3,help='Default only three controls; extension requires explicit limit')
    a=p.parse_args()
    if a.prepare:prepare(a.table,a.with_matched_benchmark)
    elif a.execute:
        if a.limit<1:raise SystemExit('positive limit required')
        run(json.loads((OUT/'native_smd_queue.json').read_text()),a.backend,a.limit)
    else:p.error('--prepare or --execute required')


if __name__=='__main__':main()
