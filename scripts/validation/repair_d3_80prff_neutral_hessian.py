#!/usr/bin/env python
"""Reproduce and tighten CPHF for one failed native neutral without relaxing QC."""
import fcntl,hashlib,json,os,shutil,socket,time
from pathlib import Path
import numpy as np
from redox.core.common import ROOT,DFT
from redox.qm import native_thermal as N,dft as D
from redox.core.protocol import load_record,record_path
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10/d3_80prff_neutral_hessian_repair'
SOURCE_BYTES=Path(__file__).read_bytes()

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(path,data):
 path=Path(path);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2));tmp.replace(path)

def main():
 assert sorted(os.sched_getaffinity(0))==[22,23] and os.getenv('CUDA_VISIBLE_DEVICES')
 assert os.getloadavg()[0]+2<os.cpu_count()
 BASE.mkdir(parents=True,exist_ok=True);(BASE/'producer_source.py').write_bytes(SOURCE_BYTES)
 sd=DFT/'d3_80prff/neu';raw=json.loads((sd/'result.json').read_text());q,m=raw['charge'],raw['mult'];assert (q,m)==(0,1)
 p,ident=N.cache_path(sd,q,m);assert load_record(sd,q,m) and N.load_native_thermal(sd,q,m) is None
 source=dict(identity=ident,raw_sha256=sha(sd/'result.json'),SP_sha256=sha(record_path(sd,q,m)))
 status=dict(status='running',source=source,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],gpu=os.getenv('CUDA_VISIBLE_DEVICES'),cpus=[22,23],started=time.time(),stages=[],gate_relaxed=False)
 def checkpoint():save(BASE/'status.json',status)
 checkpoint();factory=D._build_smd_mf
 try:
  with (ROOT/'logs/pipeline_locks/d3_80prff_native_thermal.lock').open('a+') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
   for label,tolerance,cycles in [('default_reproduction',None,100),('tight_response_retry',1e-10,200)]:
    work=BASE/label
    if work.exists():raise RuntimeError('Preserve previous diagnostic outputs before repeating a stage')
    work.mkdir();shutil.copy2(sd/'opt.xyz',work/'opt.xyz');shutil.copy2(sd/'result.json',work/'result.json')
    stage=dict(label=label,requested_response_tolerance=tolerance,max_cycle=cycles,started=time.time(),status='computing');status['stages'].append(stage);checkpoint()
    def build(*args,**kwargs):
     mf=factory(*args,**kwargs)
     if tolerance is not None:mf.conv_tol_cpscf=tolerance
     stage['actual_conv_tol_cpscf']=float(mf.conv_tol_cpscf)
     original_hessian=mf.Hessian
     def make_hessian(*a,**kw):
      # SMD undo_solvent copies instance attributes; remove this adapter
      # temporarily so its vacuum Hessian dispatch cannot recurse.
      del mf.__dict__['Hessian']
      try:h=original_hessian(*a,**kw)
      finally:mf.Hessian=make_hessian
      kernel=h.kernel
      def captured_kernel(*ka,**kk):
       h.max_cycle=cycles;value=kernel(*ka,**kk)
       host=np.asarray(value.get() if hasattr(value,'get') else value)
       matrix=host.transpose(0,2,1,3).reshape(3*mf.mol.natm,3*mf.mol.natm)
       np.save(work/'actual_full_response_hessian.npy',host);np.save(work/'actual_full_response_matrix.npy',matrix)
       finite=np.isfinite(matrix);diagnostic=dict(finite=bool(finite.all()),nonfinite_count=int((~finite).sum()),maximum_asymmetry=float(abs(matrix-matrix.T).max()) if finite.all() else None,actual_conv_tol_cpscf=float(mf.conv_tol_cpscf),response_max_cycle=h.max_cycle,grid_response=h.grid_response,cphf_grid_level=mf.cphf_grids.level)
       if finite.all():
        from pyscf.hessian import thermo
        spectrum=thermo.harmonic_analysis(mf.mol,host,imaginary_freq=False)['freq_wavenumber'];diagnostic.update(frequencies_cm=np.asarray(spectrum).tolist(),minimum_frequency_cm=float(min(spectrum)))
       save(work/'hessian_diagnostics.json',diagnostic);stage['hessian_diagnostics']=diagnostic;checkpoint()
       return value
      h.kernel=captured_kernel;return h
     mf.Hessian=make_hessian;return mf
    D._build_smd_mf=build
    try:
     _,record=N.ensure_native_thermal('d3_80prff','neu',state_dir=work)
     stage.update(status='native_record_verified',g_thermal_total_eV=record['g_thermal_total_eV'],native_record=str(N.cache_path(work,q,m)[0]),max_independent_derivative_error=max(c['relative_error'] for c in record['derivative_checks']))
    except Exception as exc:stage.update(status='failed_preserved',error=repr(exc))
    finally:D._build_smd_mf=factory
    stage['finished']=time.time();checkpoint();D._free_gpu_pool()
   retry=status['stages'][-1]
   if retry['status']!='native_record_verified':raise RuntimeError('Tighter response retry did not pass the unchanged physical/numerical gate')
   work=BASE/'tight_response_retry';native_path,new_ident=N.cache_path(work,q,m);record=json.loads(native_path.read_text())
   record.update(numerical_response_controls=dict(conv_tol_cpscf=1e-10,max_cycle=200,physical_protocol_unchanged=True),response_repair_diagnostics=str(work/'hessian_diagnostics.json'),response_repair_source=str(BASE/'status.json'))
   save(native_path,record);assert N.load_native_thermal(work,q,m)
   if N.identity(sd,q,m)!=source['identity'] or sha(sd/'result.json')!=source['raw_sha256'] or sha(record_path(sd,q,m))!=source['SP_sha256']:raise RuntimeError('Canonical source changed; do not adopt this result')
   with (ROOT/'logs/pipeline_locks/d3_80prff.lock').open('a+') as motif:
    fcntl.flock(motif,fcntl.LOCK_EX|fcntl.LOCK_NB)
    old_work=p.with_suffix('');archive=BASE/'superseded_partial';assert old_work.exists() and not p.exists()
    shutil.copytree(old_work,archive);sums={str(f.relative_to(archive)):sha(f) for f in archive.rglob('*') if f.is_file()};save(BASE/'superseded_partial_SHA256.json',sums)
    assert all(sha(archive/key)==value for key,value in sums.items())
    # Only the old incomplete addressed work directory is retired, after SHA archive.
    shutil.rmtree(old_work);shutil.copytree(native_path.with_suffix(''),old_work);shutil.copy2(native_path,p)
    assert N.load_native_thermal(sd,q,m)
   status.update(status='native_neutral_repaired_and_adopted',adopted_native_record=str(p),partial_archive=str(archive),source_geometry_unchanged=True,source_SP_unchanged=True,finished=time.time());checkpoint()
 except Exception as exc:status.update(status='failed_preserved',error=repr(exc),finished=time.time());checkpoint();raise
 finally:D._build_smd_mf=factory
if __name__=='__main__':main()
