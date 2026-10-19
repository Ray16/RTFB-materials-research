#!/usr/bin/env python
"""Archive and promote the conservative Fc+ minimum, retaining exact provenance."""
import fcntl,hashlib,json,shutil,time
from pathlib import Path
import numpy as np
from redox.core.common import ROOT,DFT,read_manifest,smiles_to_xyz_map
from redox.core.protocol import ACTIVE_SP,geom_sha1,record_path,load_record
from redox.qm.sp import _ENERGY_KEYS
from redox.qm.native_thermal import cache_path,implementation_controls,load_native_thermal,select_policy
from redox.qm.thermal import spin_degeneracy_eV
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 with (ROOT/'logs/pipeline_locks/ferrocene.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  status=json.loads((BASE/'fc_conservative_minimum_refinement/status.json').read_text())
  assert status['status']=='refined_reference_verified'
  branch=min(status['branches'],key=lambda b:b['electronic_opt_energy_Ha'])
  assert branch['minimum_and_derivatives_verified'] and branch['minimum_frequency_cm']>0
  source=Path(branch['geometry']);source_dir=source.parent
  occupation=json.loads((BASE/f'refined_fc_smd_occupation_checks/{source_dir.name}.json').read_text())
  assert occupation['geometry_sha256']==sha(source)
  best=occupation['best'];assert best['converged'] and abs(best['smd_diagnostics']['spin_contam'])<.1
  row=next(r for r in read_manifest() if r['id']=='ferrocene' and r['state']=='ox')
  assert smiles_to_xyz_map(row['smiles'],source) is not None
  sd=DFT/'ferrocene/ox';old=json.loads((sd/'result.json').read_text())
  if old.get('minimum_promotion_source')==str(source):raise RuntimeError('Already promoted; preserve current state')
  assert sha(sd/'opt.xyz')==status['source_sha256']
  archive=BASE/f'fc_minimum_adoption_{time.time_ns()}';archive.mkdir()
  shutil.copytree(sd,archive/'canonical_ox');shutil.copytree(DFT/'ferrocene/reorg',archive/'dependent_reorg')
  checksums={str(p.relative_to(archive)):sha(p) for p in archive.rglob('*') if p.is_file()}
  assert all(sha(archive/k)==v for k,v in checksums.items())
  (archive/'SHA256.json').write_text(json.dumps(checksums,indent=2))
  energy=dict(branch['active_protocol_energy']);energy.update(e_smd_Ha=best['e_Ha'],e_smd_eV=best['e_eV'],
      smd_occupation_guesses=['rm0','rm1'],converged_smd=True)
  energy['scf_guess_energies_smd'] += [dict(guess=g['guess'],e_Ha=g['e_Ha'],converged=g['converged']) for g in occupation['guesses']]
  for k,v in best['smd_diagnostics'].items():energy['smd_'+k]=v
  energy['dG_solv_eV']=energy['e_smd_eV']-energy['e_gas_eV']
  result=dict(energy,optimized=True,id='ferrocene',state='ox',n_e=old.get('n_e'),opt_out=str(sd/'opt.xyz'),
    g_thermal_eV=None,n_imag=None,thermal_method=None,legacy_thermal_archive=str(archive/'canonical_ox/result.json'),
    minimum_promotion_source=str(source),minimum_promotion_source_sha256=sha(source),
    conservative_gradient_max=branch['gradient_max'],e_opt_smd_Ha=branch['electronic_opt_energy_Ha'])
  shutil.copy2(source,sd/'opt.xyz');(sd/'result.json').write_text(json.dumps(result,indent=2))
  rp=record_path(sd,1,2);assert not rp.exists()
  rec=dict(protocol=ACTIVE_SP,protocol_hash=rp.stem,geom=str(sd/'opt.xyz'),geom_sha1=geom_sha1(sd/'opt.xyz'),
      charge=1,mult=2,id='ferrocene',state='ox',source='verified_refined_minimum:phase_matched_occupation_checks',
      software=branch['software'],occupation_proof=str(BASE/f'refined_fc_smd_occupation_checks/{source_dir.name}.json'),
      **{k:energy.get(k) for k in _ENERGY_KEYS})
  rp.parent.mkdir(exist_ok=True);rp.write_text(json.dumps(rec,indent=2));assert load_record(sd,1,2)
  # Retire all old Fc cross files after verifying their complete archive. Current hashes differ.
  for p in (DFT/'ferrocene/reorg').glob('*.json'):p.unlink()
  native,ident=cache_path(sd,1,2);work=native.with_suffix('');work.mkdir(parents=True)
  for name in ['hessian','thermo_25.log','thermo_50.log','thermo_100.log']:shutil.copy2(source_dir/name,work/name)
  shutil.copy2(source,work/'mol.xyz');matrix=np.load(source_dir/'hessian_matrix.npy');np.save(work/'hessian_cartesian.npy',matrix)
  g=branch['nuclear_by_cutoff_eV']['50'];spin=spin_degeneracy_eV(2)
  checks=[dict(atom_index=c['atom'],axis=c['axis'],step_Bohr=c['step_Bohr'],relative_error=c['relative_error']) for c in branch['derivative_checks']]
  th=dict(id='ferrocene',state='ox',status='derivative_checked_solution_thermal',identity=ident,
      stationary_minimum=True,converged=True,grid_response=True,gradient_max=branch['gradient_max'],
      maximum_asymmetry=float(abs(matrix-matrix.T).max()),frequencies_cm=branch['frequencies_cm'],
      derivative_checks=checks,smd_spin_diagnostics=branch['smd_spin_diagnostics'],software=branch['software'],
      implementation_controls=implementation_controls(),g_nuclear_qrrho_eV=g,g_spin_partition_eV=spin,
      g_thermal_total_eV=g+spin,nuclear_by_cutoff_eV=branch['nuclear_by_cutoff_eV'],
      adopted_from_verified_minimum=str(source_dir),experimental_accuracy_validated=False)
  native.write_text(json.dumps(th,indent=2));assert load_native_thermal(sd,1,2)
  select_policy(sd)
  report=dict(status='verified_minimum_promoted',archive=str(archive),old_geometry_sha256=status['source_sha256'],
      new_geometry_sha256=sha(sd/'opt.xyz'),selected_branch=source_dir.name,source_branch=branch,
      corrected_active_smd_eV=best['e_eV'],native_thermal_record=str(native),native_thermal_total_eV=g+spin,
      archive_file_count=len(checksums),dependent_cross_records_archived=True)
  (BASE/'fc_minimum_promotion.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='source_branch'},indent=2))
if __name__=='__main__':main()
