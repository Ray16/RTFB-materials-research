#!/usr/bin/env python
"""Bounded private repair of the proven cap-methyl rotor saddle.

Exact +/-60 degree rotation takes a threefold methyl maximum towards its
neighbouring minima without disturbing initial heavy-atom coordinates. Both
branches are optimized with conservative full-grid-response gradients; the
lowest optimized branch must pass the addressed native Hessian contract.
Original calculations and production directories are not modified.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import time

import numpy as np
from ase.io import read,write
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from pyscf import gto
from pyscf.hessian import thermo
from pyscf.geomopt.geometric_solver import optimize

from redox.core.common import ROOT,read_manifest,smiles_to_xyz_map
from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
from redox.qm.native_thermal import ensure_native_thermal,load_native_thermal,cache_path,NATIVE_THERMAL

SD=ROOT/'results/validation/physical_reorg_2026-10-10/conformer_controls/local_unconstrained'
OUT=ROOT/'results/validation/physical_reorg_2026-10-11/cap_rotor_repair'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(r):
    tmp=OUT/'status.tmp';tmp.write_text(json.dumps(r,indent=2,allow_nan=False));tmp.replace(OUT/'status.json')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'source_archive').exists():raise RuntimeError('Existing repair must be inspected and resumed explicitly')
    shutil.copytree(SD,OUT/'source_archive')
    failure=ROOT/'results/validation/physical_reorg_2026-10-11/graft_local_hessian_status.json'
    shutil.copy2(failure,OUT/'source_hessian_failure.json')
    sums={str(p.relative_to(OUT/'source_archive')):sha(p) for p in (OUT/'source_archive').rglob('*') if p.is_file()}
    (OUT/'source_archive_SHA256.json').write_text(json.dumps(sums,indent=2))
    if any(sha(OUT/'source_archive'/p)!=v for p,v in sums.items()):raise RuntimeError('Archive hash verification failed')
    raw=json.loads((SD/'result.json').read_text());q,m=int(raw['charge']),int(raw['mult'])
    atoms=read(SD/'opt.xyz');symbols=atoms.get_chemical_symbols();n=len(atoms)
    hfile=next((SD/'native_thermal').glob('*/hessian_cartesian.npy'));matrix=np.load(hfile)
    dummy=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),basis='sto-3g',charge=q,spin=m-1,verbose=0)
    modes=thermo.harmonic_analysis(dummy,matrix.reshape(n,3,n,3).transpose(0,2,1,3),imaginary_freq=False)
    k=int(np.argmin(modes['freq_wavenumber']));mode=np.asarray(modes['norm_mode'][k])
    molecule=Chem.MolFromXYZFile(str(SD/'opt.xyz'));rdDetermineBonds.DetermineConnectivity(molecule)
    weights=np.sum(mode*mode,axis=1);dominant=int(np.argmax(weights))
    if symbols[dominant]!='H':raise RuntimeError('Proven instability is not the targeted methyl rotor')
    neighbours=[a.GetIdx() for a in molecule.GetAtomWithIdx(dominant).GetNeighbors()]
    if len(neighbours)!=1 or symbols[neighbours[0]]!='C':raise RuntimeError('Target hydrogen lacks unique methyl carbon')
    carbon=neighbours[0]
    hydrogens=[a.GetIdx() for a in molecule.GetAtomWithIdx(carbon).GetNeighbors() if a.GetSymbol()=='H']
    heavy=[a.GetIdx() for a in molecule.GetAtomWithIdx(carbon).GetNeighbors() if a.GetSymbol()!='H']
    if len(hydrogens)!=3 or len(heavy)!=1:raise RuntimeError('Target carbon is not terminal CH3')
    if weights[hydrogens].sum()/weights.sum()<.8:raise RuntimeError('Instability not localized sufficiently on cap CH3')
    centre=atoms.positions[carbon];axis=centre-atoms.positions[heavy[0]];axis/=np.linalg.norm(axis)
    r=dict(status='optimizing_bounded_rotor_branches',host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().split()[21],cpus=sorted(os.sched_getaffinity(0)),
        gpu=os.getenv('CUDA_VISIBLE_DEVICES'),started=time.time(),source_geometry_sha256=sha(SD/'opt.xyz'),
        source_frequency_cm=float(modes['freq_wavenumber'][k]),methyl_carbon=carbon,
        methyl_hydrogens=hydrogens,max_branches=2,branches=[],production_modified=False)
    save(r)
    def terminate(sig,frame):raise KeyboardInterrupt(f'termination signal{sig}')
    signal.signal(signal.SIGTERM,terminate)
    smiles=next(x['smiles'] for x in read_manifest() if x['id']=='anthraquinone')
    basis=NATIVE_THERMAL['anion_basis'] if q<0 else NATIVE_THERMAL['basis']
    candidates=[]
    for sign in (-1,1):
        work=OUT/f'branch_{sign}';work.mkdir();angle=sign*np.pi/3
        seed=atoms.copy()
        for i in hydrogens:
            v=atoms.positions[i]-centre
            seed.positions[i]=centre+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+axis*np.dot(axis,v)*(1-np.cos(angle))
        write(work/'seed.xyz',seed)
        mol=gto.M(atom=list(zip(symbols,seed.positions.tolist())),basis=basis,ecp=_ecp_for(symbols,basis),
                  charge=q,spin=m-1,verbose=0,max_memory=6000)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
        mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
        mf.kernel()
        if not mf.converged:raise RuntimeError('Rotor seed SCF failed')
        grad=mf.nuc_grad_method();grad.grid_response=True
        final=optimize(grad.as_scanner(),maxsteps=120,assert_convergence=True,
            convergence_energy=1e-9,convergence_grms=1e-6,convergence_gmax=3e-6,
            convergence_drms=3e-5,convergence_dmax=5e-5)
        geom=atoms.copy();geom.positions=final.atom_coords(unit='Angstrom');write(work/'opt.xyz',geom)
        if smiles_to_xyz_map(smiles,work/'opt.xyz') is None:raise RuntimeError('Rotor repair changes declared connectivity')
        trial=dict(raw,optimized=True,opt_out=str(work/'opt.xyz'),g_thermal_eV=None,n_imag=None,
                   rotor_repair_source=str(SD/'opt.xyz'),rotor_repair_source_sha256=r['source_geometry_sha256'])
        for key in list(trial):
            if key.startswith(('e_smd','e_gas','scf_guess_energies','gas_','dG_solv')):trial[key]=None
        (work/'result.json').write_text(json.dumps(trial,indent=2))
        finalmol=gto.M(atom=list(zip(symbols,geom.positions.tolist())),basis=basis,ecp=_ecp_for(symbols,basis),
                       charge=q,spin=m-1,verbose=0,max_memory=6000)
        fm=_build_smd_mf(finalmol,'r2scan','acetonitrile','gpu',disp='d4');fm.conv_tol=1e-12;fm.conv_tol_grad=1e-8;fm.kernel()
        if not fm.converged:raise RuntimeError('Optimized rotor branch SCF failed')
        item=dict(sign=sign,status='optimized',geometry=str(work/'opt.xyz'),opt_level_e_Ha=float(fm.e_tot))
        r['branches'].append(item);save(r);candidates.append((float(fm.e_tot),work,item))
        del mf,grad,fm;_free_gpu_pool()
    for energy,work,item in sorted(candidates):
        try:
            r['status']='verifying_selected_rotor_minimum';r['selected']=str(work);save(r)
            _,thermal=ensure_native_thermal('anthraquinone','red1',state_dir=work)
            if load_native_thermal(work,q,m) is None:raise RuntimeError('Native minimum record did not validate')
            p,_=cache_path(work,q,m);item.update(status='verified_minimum',native_record=str(p),minimum_frequency_cm=min(thermal['frequencies_cm']))
            r.update(status='verified_private_rotor_minimum',selected=str(work),selected_native_record=str(p),finished=time.time())
            save(r);return
        except Exception as exc:item.update(status='minimum_failed',error=repr(exc));save(r)
    r.update(status='bounded_repair_failed',finished=time.time());save(r)
    raise RuntimeError('Both rotor branches failed the physical minimum gate')


if __name__=='__main__':main()
