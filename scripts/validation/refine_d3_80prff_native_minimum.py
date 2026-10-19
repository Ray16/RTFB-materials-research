#!/usr/bin/env python
"""Bounded conservative refinement away from a diagnosed sharp CDS contact."""
import fcntl,hashlib,json,os,shutil,socket,time
from pathlib import Path
import numpy as np
from ase.io import read,write
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from redox.core.common import ROOT,DFT,read_manifest,smiles_to_xyz_map
from redox.core.protocol import load_record,record_path,ACTIVE_SP
from redox.qm import dft as D,native_thermal as N
from redox.qm.sp import ensure,_ENERGY_KEYS
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10/d3_80prff_neutral_hessian_repair/conservative_refinement'
BYTES=Path(__file__).read_bytes()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,r):
 p=Path(p);t=p.with_suffix('.tmp');t.write_text(json.dumps(r,indent=2));t.replace(p)
def archive(src,dst):
 shutil.copytree(src,dst);s={str(p.relative_to(dst)):sha(p) for p in dst.rglob('*') if p.is_file()};save(dst.parent/(dst.name+'_SHA256.json'),s);assert all(sha(dst/k)==v for k,v in s.items())
def main():
 assert sorted(os.sched_getaffinity(0))==[22,23] and os.getenv('CUDA_VISIBLE_DEVICES') and os.getloadavg()[0]+2<os.cpu_count()
 BASE.mkdir(parents=True,exist_ok=True);(BASE/'producer_source.py').write_bytes(BYTES)
 sd=DFT/'d3_80prff/neu';raw=json.loads((sd/'result.json').read_text());q,m=raw['charge'],raw['mult'];assert (q,m)==(0,1)
 source=dict(geometry_sha256=sha(sd/'opt.xyz'),raw_sha256=sha(sd/'result.json'),SP_sha256=sha(record_path(sd,q,m)))
 status=dict(status='running',source=source,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],gpu=os.getenv('CUDA_VISIBLE_DEVICES'),cpus=[22,23],started=time.time(),branches=[],unchanged_native_QC=True)
 def checkpoint():save(BASE/'status.json',status)
 checkpoint();atoms=read(sd/'opt.xyz');symbols=atoms.get_chemical_symbols();smiles=next(r['smiles'] for r in read_manifest() if r['id']=='d3_80prff' and r['state']=='neu');factory=D._build_smd_mf
 def tight(*args,**kwargs):
  mf=factory(*args,**kwargs);mf.conv_tol_cpscf=1e-10;return mf
 try:
  with (ROOT/'logs/pipeline_locks/d3_80prff_native_thermal.lock').open('a+') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);good=[]
   for sign in [-1,1]:
    work=BASE/('minus' if sign<0 else 'plus');work.mkdir();item=dict(sign=sign,status='optimizing');status['branches'].append(item);checkpoint()
    try:
     seed=atoms.copy();seed.positions[28]+=sign*.02*np.array([1.,0.,1.])/np.sqrt(2);write(work/'seed.xyz',seed)
     mol=gto.M(atom=list(zip(symbols,seed.positions.tolist())),basis='def2-svp',charge=0,spin=0,verbose=0,max_memory=6000)
     mf=factory(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.kernel()
     if not mf.converged:raise RuntimeError('Conservative refinement SCF failed')
     gradient=mf.nuc_grad_method();gradient.grid_response=True
     final=optimize(gradient.as_scanner(),maxsteps=120,assert_convergence=True,convergence_energy=1e-8,convergence_grms=3e-5,convergence_gmax=5e-5,convergence_drms=1e-4,convergence_dmax=2e-4)
     out=atoms.copy();out.positions=final.atom_coords(unit='Angstrom');write(work/'opt.xyz',out);del mf,gradient;D._free_gpu_pool()
     if smiles_to_xyz_map(smiles,work/'opt.xyz') is None:raise RuntimeError('Refinement changes chemical connectivity')
     r=dict(raw,optimized=True,g_thermal_eV=None,n_imag=None,conservative_refinement_source=source)
     for key in _ENERGY_KEYS:r[key]=None
     save(work/'result.json',r);D._build_smd_mf=tight
     try:_,thermal=N.ensure_native_thermal('d3_80prff','neu',state_dir=work)
     finally:D._build_smd_mf=factory
     _,energy=ensure('d3_80prff','neu',work,q,m,backend='gpu');energy=load_record(work,q,m)
     if energy is None:raise RuntimeError('Refined active SP failed')
     r.update({k:energy.get(k) for k in _ENERGY_KEYS});r.update(sp_xc=ACTIVE_SP['xc'],sp_basis=ACTIVE_SP['basis'],sp_disp=ACTIVE_SP['disp'],sp_nlc=ACTIVE_SP['nlc'],xc=ACTIVE_SP['xc'],basis=ACTIVE_SP['basis']);save(work/'result.json',r)
     p,_=N.cache_path(work,q,m);record=json.loads(p.read_text());record['numerical_response_controls']=dict(conv_tol_cpscf=1e-10,physical_protocol_unchanged=True);save(p,record);assert N.load_native_thermal(work,q,m)
     item.update(status='verified',geometry_sha256=sha(work/'opt.xyz'),e_smd_eV=energy['e_smd_eV'],g_thermal_eV=thermal['g_thermal_total_eV'],min_frequency_cm=min(thermal['frequencies_cm']),maximum_asymmetry=thermal['maximum_asymmetry'],max_derivative_error=max(c['relative_error'] for c in thermal['derivative_checks']));good.append((energy['e_smd_eV'],work))
    except Exception as exc:item.update(status='failed_preserved',error=repr(exc))
    checkpoint();D._free_gpu_pool()
   if not good:raise RuntimeError('Both conservative refinement branches failed unchanged native physical gate')
   _,selected=min(good,key=lambda t:t[0])
   with (ROOT/'logs/pipeline_locks/d3_80prff.lock').open('a+') as motif:
    fcntl.flock(motif,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if sha(sd/'opt.xyz')!=source['geometry_sha256'] or sha(sd/'result.json')!=source['raw_sha256'] or sha(record_path(sd,q,m))!=source['SP_sha256']:raise RuntimeError('Canonical source changed; do not promote')
    archive(sd,BASE/'superseded_canonical');reorg=DFT/'d3_80prff/reorg'
    if reorg.exists():archive(reorg,BASE/'superseded_reorg')
    prepared=BASE/'prepared_canonical';shutil.copytree(selected,prepared);N.select_policy(prepared);assert load_record(prepared,q,m) and N.load_native_thermal(prepared,q,m)
    retired=BASE/'retired_live_canonical';sd.rename(retired)
    try:prepared.rename(sd)
    except Exception:retired.rename(sd);raise
    if reorg.exists():
     for p in reorg.glob('*.json'):p.unlink()
   status.update(status='verified_native_minimum_promoted',selected=str(selected),current_geometry_sha256=sha(sd/'opt.xyz'),adopted_native=str(N.cache_path(sd,q,m)[0]),cross_refresh_required=True,finished=time.time());checkpoint()
 except Exception as exc:status.update(status='failed_preserved',error=repr(exc),finished=time.time());checkpoint();raise
 finally:D._build_smd_mf=factory
if __name__=='__main__':main()
