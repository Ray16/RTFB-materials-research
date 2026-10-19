"""Bounded conservative minimum repair, then SHA-archived canonical promotion.

Only entered for an actual native gradient/stationarity failure. Every proposed
branch must pass the same native Hessian and connectivity contract before adoption.
"""
import fcntl,hashlib,json,shutil,time,os,socket
from contextlib import contextmanager
from pathlib import Path
import numpy as np
from redox.core.common import ROOT,DFT,read_manifest,smiles_to_xyz_map
from redox.qm.native_thermal import ensure_native_thermal,cache_path,load_native_thermal,select_policy,NATIVE_THERMAL
from redox.core.protocol import load_record
from redox.core.io import atomic_write_json

REPAIR_SOURCE_BYTES=Path(__file__).read_bytes()
REPAIR_SOURCE_SHA256=hashlib.sha256(REPAIR_SOURCE_BYTES).hexdigest()


def promotion_metadata(record,canonical,archive,expected_source_sha256,producer_snapshot):
 """Bind original geometry to its verified archive without inventing SP provenance."""
 original=Path(archive)/'opt.xyz'
 if hashlib.sha256(original.read_bytes()).hexdigest()!=expected_source_sha256:
  raise RuntimeError('Archived repair source geometry changed')
 snapshot=Path(producer_snapshot)
 if hashlib.sha256(snapshot.read_bytes()).hexdigest()!=REPAIR_SOURCE_SHA256:
  raise RuntimeError('Repair producer snapshot changed')
 updated=dict(record)
 # A repair copies a separately addressed SP; it is not a new direct producer
 # of the geometry-bound energy proof in result.json.
 updated.pop('energy_protocol',None);updated.pop('energy_geom_sha1',None)
 updated.update(opt_out=str(Path(canonical)/'opt.xyz'),
                minimum_repair_source=str(original),
                minimum_repair_source_sha256=expected_source_sha256,
                legacy_thermal_archive=str(Path(archive)/'result.json'),
                repair_producer_snapshot=str(snapshot),repair_producer_sha256=REPAIR_SOURCE_SHA256)
 return updated


@contextmanager
def promotion_guard(gid):
 """Acquire the motif lock or verify the same local runner's live inherited ownership."""
 path=ROOT/'logs/pipeline_locks'/f'{gid}.lock'
 with path.open('a+') as lock:
  try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:
   pid=int(os.environ.get('REDOX_PIPELINE_OWNER_PID','0'))
   start=os.environ.get('REDOX_PIPELINE_OWNER_START')
   proc=Path(f'/proc/{pid}')
   if pid<=0 or not proc.exists() or (proc/'stat').read_text().split()[21]!=start:
    raise RuntimeError('Another writer owns the motif; promotion refused')
   command=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
   if 'scripts/pipeline/run_motifs.py' not in command:raise RuntimeError('Unverified parent motif lock owner')
   inode=path.stat().st_ino
   owned=False
   for fd in (proc/'fd').iterdir():
    try:
     info=(proc/'fdinfo'/fd.name).read_text()
     if fd.stat().st_ino==inode and 'FLOCK' in info and 'WRITE' in info:owned=True
    except OSError:pass
   if not owned:raise RuntimeError('Parent does not own the declared exclusive motif lock')
  yield


def repair_minimum(gid,state,reason):
 import inspect,importlib
 import redox.qm.native_thermal as native
 # An already running cohort may have imported the pre-adapter function. Reload
 # only the directory adapter; the scientific protocol and equations are identical.
 if 'state_dir' not in inspect.signature(native.ensure_native_thermal).parameters:
  importlib.reload(native)
 ensure_native_thermal=native.ensure_native_thermal
 from ase.io import read,write
 from pyscf import gto
 from pyscf.hessian import thermo
 from pyscf.geomopt.geometric_solver import optimize
 from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
 from redox.qm.sp import ensure as ensure_sp,_ENERGY_KEYS
 sd=DFT/gid/state;raw=json.loads((sd/'result.json').read_text());q,m=raw['charge'],raw['mult']
 source_sha=hashlib.sha256((sd/'opt.xyz').read_bytes()).hexdigest()
 base=ROOT/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs'/gid/state/str(time.time_ns());base.mkdir(parents=True)
 (base/'repair_producer_source.py').write_bytes(REPAIR_SOURCE_BYTES)
 status=dict(status='running',id=gid,state=state,reason=reason,source_sha256=source_sha,branches=[],max_branches=2)
 def save():
  p=base/'status.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
 save();atoms=read(sd/'opt.xyz');symbols=atoms.get_chemical_symbols();nat=len(atoms)
 basis=NATIVE_THERMAL['anion_basis'] if q<0 else NATIVE_THERMAL['basis']
 def build(coords):
  mol=gto.M(atom=list(zip(symbols,coords.tolist())),charge=q,spin=m-1,basis=basis,ecp=_ecp_for(symbols,basis),verbose=0,max_memory=6000)
  mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
  return mol,mf
 p,_=cache_path(sd,q,m);h=p.with_suffix('')/'hessian_cartesian.npy'
 if h.exists():
  mol,_=build(atoms.positions);mat=np.load(h);H=mat.reshape(nat,3,nat,3).transpose(0,2,1,3)
  data=thermo.harmonic_analysis(mol,H,imaginary_freq=False);k=int(np.argmin(data['freq_wavenumber']))
  mode=np.asarray(data['norm_mode'][k]);mode/=max(np.linalg.norm(mode,axis=1));seeds=[(-1,-.1*mode),(1,.1*mode)]
  status['source_min_frequency_cm']=float(data['freq_wavenumber'][k])
 else:seeds=[(0,np.zeros_like(atoms.positions))]
 row=next(r for r in read_manifest() if r['id']==gid and r['state']==state)
 previous=[]
 for saved in base.parent.glob('*/status.json'):
  try:
   evidence=json.loads(saved.read_text())
   if saved.parent!=base and evidence.get('source_sha256')==source_sha:
    previous.append(saved.parent)
  except (OSError,ValueError):pass
 good=[]
 for sign,delta in seeds:
  work=base/f'branch_{sign}';work.mkdir();record=dict(sign=sign,status='optimizing');status['branches'].append(record);save()
  try:
   prior=next((folder/f'branch_{sign}'/'opt.xyz' for folder in reversed(previous)
               if (folder/f'branch_{sign}'/'opt.xyz').exists()),None)
   if prior is not None:
    # Only a geometry checkpoint is reused. SCF/Hessian/SP/graph are revalidated below.
    shutil.copy2(prior,work/'opt.xyz');record['resumed_geometry_checkpoint']=str(prior)
   else:
    seed=atoms.copy();seed.positions+=delta;write(work/'seed.xyz',seed)
    _,mf=build(seed.positions);mf.kernel()
    if not mf.converged:raise RuntimeError('Repair SCF unconverged')
    grad=mf.nuc_grad_method();grad.grid_response=True
    finalmol=optimize(grad.as_scanner(),maxsteps=120,assert_convergence=True,
        convergence_energy=1e-8,convergence_grms=3e-5,convergence_gmax=5e-5,convergence_drms=1e-4,convergence_dmax=2e-4)
    final=atoms.copy();final.positions=finalmol.atom_coords(unit='Angstrom');write(work/'opt.xyz',final)
    del mf,grad;_free_gpu_pool()
   if smiles_to_xyz_map(row['smiles'],work/'opt.xyz') is None:raise RuntimeError('Repair changes declared molecular connectivity')
   trial=dict(raw,optimized=True,g_thermal_eV=None,n_imag=None,minimum_repair_source=str(sd/'opt.xyz'),minimum_repair_source_sha256=source_sha)
   # Old energies refer to the original geometry and must never be adoptable.
   for key in _ENERGY_KEYS:trial[key]=None
   (work/'result.json').write_text(json.dumps(trial,indent=2))
   kind,th=ensure_native_thermal(gid,state,state_dir=work)
   energy_kind,energy=ensure_sp(gid,state,work,q,m,backend='gpu')
   if not load_record(work,q,m):raise RuntimeError('Repair active SP incomplete')
   trial.update({key:energy.get(key) for key in _ENERGY_KEYS})
   proto=energy['protocol']
   trial.update(sp_xc=proto['xc'],sp_basis=proto['basis'],sp_disp=proto['disp'],sp_nlc=proto['nlc'],
                xc=proto['xc'],basis=proto['basis'],software=energy.get('software'))
   trial.update(opt_out=str(work/'opt.xyz'),optimized=True,legacy_thermal_archive=str(base/'superseded_canonical/result.json'))
   (work/'result.json').write_text(json.dumps(trial,indent=2))
   record.update(status='verified',geometry=str(work/'opt.xyz'),g_thermal_eV=th['g_thermal_total_eV'],e_smd_eV=energy['e_smd_eV'])
   good.append((energy['e_smd_eV'],work,th))
  except Exception as exc:record.update(status='failed',error=repr(exc))
  save()
 if not good:status['status']='repair_failed';save();raise RuntimeError(f'Bounded minimum repair failed; evidence {base}')
 _,selected,thermal=min(good,key=lambda x:x[0])
 with promotion_guard(gid):
  if hashlib.sha256((sd/'opt.xyz').read_bytes()).hexdigest()!=source_sha:raise RuntimeError('Canonical source changed during bounded repair')
  archive=base/'superseded_canonical';shutil.copytree(sd,archive)
  reorg=DFT/gid/'reorg'
  if reorg.exists():shutil.copytree(reorg,base/'superseded_reorg')
  sums={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in [archive,base/'superseded_reorg'] if folder.exists() for p in folder.rglob('*') if p.is_file()}
  (base/'SHA256.json').write_text(json.dumps(sums,indent=2))
  assert all(hashlib.sha256((base/k).read_bytes()).hexdigest()==v for k,v in sums.items())
  # Archive first; retire only geometry-specific canonical files after the verified replacement is ready.
  replacement=promotion_metadata(json.loads((selected/'result.json').read_text()),sd,archive,source_sha,base/'repair_producer_source.py')
  shutil.copy2(selected/'opt.xyz',sd/'opt.xyz')
  atomic_write_json(sd/'result.json',replacement)
  for name in ['native_thermal','sp']:
   for src in (selected/name).iterdir():
    dst=sd/name/src.name
    if dst.exists():raise RuntimeError('Replacement addressed cache unexpectedly already exists')
    if src.is_dir():shutil.copytree(src,dst)
    else:shutil.copy2(src,dst)
  if reorg.exists():
   for old in reorg.glob('*.json'):old.unlink()
  select_policy(sd)
  assert load_native_thermal(sd,q,m) and load_record(sd,q,m)
 status.update(status='verified_minimum_promoted',selected=str(selected),new_geometry_sha256=hashlib.sha256((sd/'opt.xyz').read_bytes()).hexdigest(),archive=str(archive),finished=time.time());save()
 return 'repaired_minimum',thermal
