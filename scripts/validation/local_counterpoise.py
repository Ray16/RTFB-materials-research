"""Gas-only counterpoise with a fixed full-complex quadrature and RI space.

Ghost centers add AO functions without nuclei or electrons. All five terms use
one real-complex integration grid and one full-complex auxiliary basis, so the
own/ghost difference isolates AO-basis expansion within that numerical model.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

SETTINGS=dict(xc='wb97m-v',basis='def2-tzvpd',nlc='vv10',density_fit=True,only_dfj=False,
    grid_level=3,nlc_grid_level=3,small_rho_cutoff=0.,SCF_conv_tol=1e-10,SCF_conv_tol_grad=1e-7,
    SCF_max_cycle=200,medium='gas',quadrature='One CPU-generated full-real-complex grid for every term',
    auxiliary_space='One default PySCF full-real-complex RI-JK auxiliary molecule for every term')
ROLES=('AB','A_own','B_own','A_ghost','B_ghost')
EV_KJ=1.602176634e-19*6.02214076e23/1000
HA_EV=27.211386245988


def molecules(symbols,coordinates,A):
    from pyscf import gto
    coordinates=np.asarray(coordinates,dtype=float);n=len(symbols);A=set(A)
    if coordinates.shape!=(n,3) or not np.isfinite(coordinates).all() or not A or not A<set(range(n)):
        raise ValueError('A finite complex and two nonempty complementary atom sets are required')
    output={}
    for role in ROLES:
        real=set(range(n)) if role=='AB' else A if role.startswith('A_') else set(range(n))-A
        selected=range(n) if role=='AB' or role.endswith('ghost') else sorted(real)
        atoms=[(symbols[i] if i in real else 'ghost-'+symbols[i],coordinates[i].tolist()) for i in selected]
        output[role]=gto.M(atom=atoms,unit='Angstrom',basis=SETTINGS['basis'],charge=0,spin=0,verbose=0,max_memory=6000)
    if (output['AB'].nelectron!=output['A_ghost'].nelectron+output['B_ghost'].nelectron
            or any(output[r].nao_nr()!=output['AB'].nao_nr() for r in ('A_ghost','B_ghost'))
            or output['A_own'].nao_nr()+output['B_own'].nao_nr()!=output['AB'].nao_nr()):
        raise ValueError('Ghost nuclei/electrons or AO coverage differ from counterpoise definition')
    return output


def sums(energies):
    if set(energies)!=set(ROLES) or not all(type(v) in (int,float) and math.isfinite(v) for v in energies.values()):
        raise ValueError('Five finite common-grid common-auxiliary energies required')
    own=math.fsum([energies['AB'],-energies['A_own'],-energies['B_own']])
    cp=math.fsum([energies['AB'],-energies['A_ghost'],-energies['B_ghost']])
    shift=math.fsum([energies['A_own'],-energies['A_ghost'],energies['B_own'],-energies['B_ghost']])
    if abs(cp-own-shift)>1e-10:raise ValueError('Counterpoise identity does not close')
    return dict(uncorrected_common_grid_kJ_mol=own*HA_EV*EV_KJ,
        counterpoise_gas_kJ_mol=cp*HA_EV*EV_KJ,AO_basis_expansion_shift_kJ_mol=shift*HA_EV*EV_KJ,
        solvent_counterpoise_computed=False,thermal_free_energy_computed=False)


def compute(symbols,coordinates,A,backend,out):
    from pyscf import dft,df
    from pyscf.dft import gen_grid
    from redox.core.io import atomic_write_json
    mols=molecules(symbols,coordinates,A);full=mols['AB']
    grid=gen_grid.Grids(full);grid.level=3;grid.build(with_non0tab=False)
    auxiliary=df.addons.make_auxmol(full,None)
    # CPU arrays form a single persisted numerical quadrature for both engines.
    out.mkdir(parents=True,exist_ok=False)
    np.save(out/'grid_coordinates_AU.npy',grid.coords);np.save(out/'grid_weights.npy',grid.weights)
    grid_hashes={k:hashlib.sha256((out/k).read_bytes()).hexdigest() for k in ['grid_coordinates_AU.npy','grid_weights.npy']}
    if backend=='gpu':
        import cupy as cp
        from gpu4pyscf import dft as engine
        array=cp.asarray
    elif backend=='cpu':engine=dft;array=np.asarray
    else:raise ValueError('Explicit CPU or GPU backend required')
    results={}
    for role in ROLES:
        folder=out/role;folder.mkdir();mol=mols[role]
        mf=engine.RKS(mol).density_fit(only_dfj=False);mf.with_df.auxmol=auxiliary.copy()
        mf.xc=SETTINGS['xc'];mf.nlc=SETTINGS['nlc'];mf.small_rho_cutoff=0.
        mf.conv_tol=SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=SETTINGS['SCF_conv_tol_grad'];mf.max_cycle=200
        mf.grids.level=3;mf.nlcgrids.level=3
        for grids in (mf.grids,mf.nlcgrids):
            grids.coords=array(grid.coords.copy());grids.weights=array(grid.weights.copy())
        mf.chkfile=str(folder/'SCF.h5')
        mf.callback=lambda env:atomic_write_json(folder/'progress.json',dict(cycle=int(env['cycle'])+1,E_Ha=float(env['e_tot']),checked=time.time(),qualified=False))
        started=time.time();mf.kernel();value=float(mf.e_tot)
        if not mf.converged or not math.isfinite(value):raise ValueError('Counterpoise SCF failed: '+role)
        record=dict(role=role,E_Ha=value,converged=True,real_nuclear_charges=mol.atom_charges().tolist(),
            electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()),auxiliary_functions=int(mf.with_df.auxmol.nao_nr()),
            auxiliary_basis=repr(mf.with_df.auxmol.basis),common_grid_hashes=grid_hashes,settings=SETTINGS,seconds=time.time()-started,
            coordinates_Bohr=mol.atom_coords().tolist(),atom_symbols=[mol.atom_symbol(i) for i in range(mol.natm)])
        atomic_write_json(folder/'energy.json',record);results[role]=record
        del mf
        import gc
        gc.collect()
        if backend=='gpu':cp.get_default_memory_pool().free_all_blocks()
    summary=dict(status='gas_counterpoise_five_terms_computed_independent_review_pending',backend=backend,settings=SETTINGS,
        grid_points=len(grid.weights),terms=results,sums=sums({k:v['E_Ha'] for k,v in results.items()}),
        common_grid_hashes=grid_hashes,source_symbols=symbols,source_coordinates_A=np.asarray(coordinates).tolist(),A_indices=sorted(A),
        solvent_counterpoise_computed=False,primary_candidate_admission=False)
    atomic_write_json(out/'summary.json',summary);return summary
