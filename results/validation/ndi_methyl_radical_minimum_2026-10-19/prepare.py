"""Prepare a bounded two-sided unstable-mode refinement without canonical writes."""
from pathlib import Path
import hashlib,json,shutil,time,os
from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
import numpy as np
from ase.io import read,write
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT,read_manifest,smiles_to_xyz_map
from redox.core.io import atomic_write_json
b=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
gid='search_ndi_573f69a7895b';state='red1';canonical=ROOT/'calcs/dft'/gid/state
case=ROOT/'results/validation/candidate_property_continuation_2026-10-18/NDI_Me_radical'
assert json.loads((case/'status.json').read_text())['status']=='failed_preserved'
source=case/'work';H=source/'native_thermal/4932cfe17c4f7e66/hessian_cartesian.npy'
for name in ['result.json','opt.xyz']:assert sha(source/name)==sha(canonical/name),name
row=next(r for r in read_manifest() if r['id']==gid and r['state']==state);atoms=read(source/'opt.xyz');n=len(atoms);q,m=-1,2
matrix=np.load(H);assert matrix.shape==(3*n,3*n) and abs(matrix-matrix.T).max()<=1e-6
mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis='def2-svpd',charge=q,spin=m-1,verbose=0)
a=thermo.harmonic_analysis(mol,matrix.reshape(n,3,n,3).transpose(0,2,1,3),imaginary_freq=False)
f=np.asarray(a['freq_wavenumber']);i=int(f.argmin());assert f[i]<-5
mode=np.asarray(a['norm_mode'][i]);mode*=.1/np.linalg.norm(mode,axis=1).max()
archive=b/'source_archive';archive.mkdir(exist_ok=False)
for name in ['result.json','opt.xyz']:shutil.copy2(canonical/name,archive/name)
template=ROOT/'results/validation/remaining_radical_signed_modes_2026-10-18/ndi_methoxyethyl'
for name in ['branch.py','qualify.py']:shutil.copy2(template/name,b/name)
seeds={}
for label,sign in [('minus',-1),('plus',1)]:
 work=b/label;work.mkdir();trial=atoms.copy();trial.positions+=sign*mode;write(work/'seed.xyz',trial)
 assert smiles_to_xyz_map(row['smiles'],work/'seed.xyz') is not None
 seeds[label]=sha(work/'seed.xyz')
inputs={str(p):sha(p) for p in [canonical/'opt.xyz',canonical/'result.json',source/'opt.xyz',source/'result.json',H]}
producers={str(b/n):sha(b/n) for n in ['branch.py','qualify.py']}
for relative in ['qm/native_thermal.py','qm/dft.py','qm/sp.py','qm/hessian_memory.py','core/common.py','core/protocol.py','qm/qrrho.py','qm/xtb_kernel.py','qm/projected_thermal.py']:
 p=ROOT/'src/redox'/relative;producers[str(p)]=sha(p)
plan=dict(id=gid,state=state,charge=q,mult=m,smiles=row['smiles'],source_state_dir=str(canonical),source_hessian=str(H),source_minimum_frequency_cm=float(f[i]),seed_max_displacement_A=.1,mode_index=i,seed_sha256=seeds,inputs_sha256=inputs,producer_sha256=producers,maximum_signed_branches=2,maximum_optimizer_steps_per_branch=120,production_admission=False,scientific_thresholds_unchanged=True,created=time.time())
atomic_write_json(b/'source_plus.json',plan);print('Prepared two exact-graph seeds; minimum frequency',float(f[i]))
