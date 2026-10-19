#!/usr/bin/env python
"""Independent trace/units/provenance audit of the new PCM diagnostic records."""
import csv
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT, DFT, HARTREE_EV
from redox.core.protocol import geom_sha1
from redox.qm.dft import _ecp_for

OUT=ROOT/'results/validation/physical_reorg_2026-10-10'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=OUT)
    parser.add_argument('--allow-preserved-solver',action='store_true')
    args=parser.parse_args();destination=args.out;destination.mkdir(parents=True,exist_ok=True)
    launch=json.loads((OUT/'launch_manifest.json').read_text())
    worker_hash=launch['producer_sha256']['scripts/validation/physical_reorg.py']
    allowed={worker_hash}
    if args.allow_preserved_solver:
        before=ast.parse((OUT/'physical_reorg_v2_preserved.py').read_text())
        after=ast.parse((ROOT/'scripts/validation/physical_reorg.py').read_text())
        for name in ('make_mf','solve','state_info'):
            a=next(n for n in before.body if isinstance(n,ast.FunctionDef) and n.name==name)
            b=next(n for n in after.body if isinstance(n,ast.FunctionDef) and n.name==name)
            if ast.dump(a,include_attributes=False)!=ast.dump(b,include_attributes=False):
                raise RuntimeError(f'Preserved producer differs: {name}')
        allowed.add(hashlib.sha256((OUT/'physical_reorg_v2_preserved.py').read_bytes()).hexdigest())
        allowed.add(hashlib.sha256((ROOT/'scripts/validation/physical_reorg.py').read_bytes()).hexdigest())
    rows=[]
    overlap_cache={}
    for path in sorted((OUT/'scf').glob('*.json')):
        rec=json.loads(path.read_text()); k=rec['key']
        row=dict(record=path.name,id=k['id'],state=k['state'],geometry=k['at_state'],
                 reference=k['reference_state'],status='passed')
        if k.get('worker_code_sha256') not in allowed:
            row['status']='excluded_superseded_worker_revision'
            rows.append(row);continue
        geom=DFT/k['id']/k['at_state']/'opt.xyz'
        errors=[]
        if geom_sha1(geom)!=k['geom_sha1']:
            errors.append('geometry_hash_changed')
        if abs(rec['e_eV']-rec['e_Ha']*HARTREE_EV)>1e-8:
            errors.append('energy_unit_mismatch')
        converged=[g for g in rec['guesses'] if g['converged']]
        if not converged or abs(rec['e_Ha']-min(g['e_Ha'] for g in converged))>1e-10:
            errors.append('energy_not_lowest_converged_guess')
        symbols=read(str(geom)).get_chemical_symbols()
        geo_key=(str(geom),k['geom_sha1'],k['protocol']['basis'])
        if geo_key not in overlap_cache:
            atom=read(str(geom))
            # Overlap is independent of charge/spin but mol parity must be valid.
            mol=gto.M(atom=list(zip(symbols,atom.positions.tolist())),
                      charge=k['charge'],spin=k['mult']-1,basis=k['protocol']['basis'],
                      ecp=_ecp_for(symbols,k['protocol']['basis']),verbose=0)
            overlap_cache[geo_key]=(mol.intor('int1e_ovlp'),sum(mol.atom_charges()))
        overlap,nuclear_electrons=overlap_cache[geo_key]
        with np.load(path.with_suffix('.npz')) as npz:
            dm=npz['dm']
        dm_total=dm.sum(axis=0) if dm.ndim==3 else dm
        count=float(np.einsum('ij,ji',dm_total,overlap))
        expected=int(nuclear_electrons-k['charge'])
        row['electron_trace']=count;row['electron_expected']=expected
        row['electron_trace_error']=count-expected
        if abs(count-expected)>1e-6 or not np.isfinite(dm).all():
            errors.append('density_electron_count_or_finiteness')
        if not np.allclose(dm_total,dm_total.T,atol=1e-10):
            errors.append('density_not_hermitian')
        energies=np.asarray(rec['mo_energy_Ha'])
        if energies.ndim==2:
            nalpha=(expected+k['mult']-1)//2;nbeta=expected-nalpha
            row['homo_eV']=float(max(energies[0,nalpha-1],energies[1,nbeta-1]))*HARTREE_EV
        else:
            row['homo_eV']=float(energies[expected//2-1])*HARTREE_EV
        if rec.get('s2') is not None:
            spin=(k['mult']-1)/2
            row['spin_contamination']=rec['s2']-spin*(spin+1)
        if errors:
            row['status']='failed';row['errors']=';'.join(errors)
        rows.append(row)
    columns=['record','id','state','geometry','reference','status','electron_trace',
             'electron_expected','electron_trace_error','homo_eV','spin_contamination','errors']
    with (destination/'independent_record_audit.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    summary=dict(passed=sum(r['status']=='passed' for r in rows),
                 failed=sum(r['status']=='failed' for r in rows),
                 superseded=sum(r['status'].startswith('excluded') for r in rows),
                 expected_cases=3,
                 cases_present=sum((OUT/f'{g}_{o}_{r}.json').exists() for g,o,r in [
                     ('anthraquinone_parent','neu','red1'),
                     ('anthraquinone_parent','red1','red2'),
                     ('anthraquinone','neu','red1')]),
                 production_promotion=False,allowed_worker_sha256=sorted(allowed),
                 solver_AST_reuse_verified=args.allow_preserved_solver)
    (destination/'independent_record_audit.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary))
    if summary['failed']:
        raise SystemExit(2)


if __name__=='__main__':main()
