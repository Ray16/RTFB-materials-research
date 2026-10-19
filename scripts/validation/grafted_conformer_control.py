#!/usr/bin/env python
"""Same-basin and linker-fixed SMD controls for grafted anthraquinone.

The constraints are an explicit sensitivity experiment. They do not establish
that every polymer fixes the graft in this conformation and do not clear flags.
Results remain in the validation directory, never in production calcs/dft.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import time

import numpy as np
from ase.io import read
from rdkit import Chem

from redox.core.common import DFT, ROOT, read_manifest, read_result, smiles_to_xyz_map
from redox.core.protocol import ACTIVE_SP, geom_sha1
from redox.qm.dft import dft_smd


OUT=ROOT/'results/validation/physical_reorg_2026-10-10/conformer_controls'


def write_linker_constraints(smiles,geometry,path):
    mol=Chem.MolFromSmiles(smiles)
    tether=set(mol.GetSubstructMatch(Chem.MolFromSmarts('[CH3]c1ccc([CH2])cc1')))
    if not tether:
        raise ValueError('Merrifield surrogate not identified')
    # Find graft-side heavy bonds that connect tether to an external atom.
    bridges=[]
    for bond in mol.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if (a in tether)!=(b in tether):
            bridges.append((a,b) if a in tether else (b,a))
    if len(bridges)!=1:
        raise ValueError(f'Expected one attachment bond; found {bridges}')
    start,_=bridges[0]
    # Anthraquinone carbonyl-bearing ring excludes isolated aryl substituents.
    from reorg_conformation_audit import localized_core_indices
    core_xyz,_,_=localized_core_indices(smiles,geometry)
    mapping=smiles_to_xyz_map(smiles,geometry)
    core={i for i,j in enumerate(mapping) if j in set(core_xyz)}
    route=min((Chem.GetShortestPath(mol,start,c) for c in core),key=len)
    actual_heavy=[i for i,s in enumerate(read(str(geometry)).get_chemical_symbols()) if s!='H']
    lines=['$freeze']
    for a,b in zip(route,route[1:]):
        na=[n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx()!=b]
        nb=[n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetIdx()!=a]
        if not na or not nb or mol.GetBondBetweenAtoms(a,b).IsInRing():
            continue
        indices=[actual_heavy[mapping[i]]+1 for i in (na[0],a,b,nb[0])]
        lines.append('dihedral '+' '.join(map(str,indices)))
    if len(lines)<2:
        raise ValueError('No linker torsion selected')
    path.write_text('\n'.join(lines)+'\n')
    return len(lines)-1


def rmsd(a,b):
    a=a-a.mean(0);b=b-b.mean(0)
    u,_,vt=np.linalg.svd(b.T@a)
    rotation=u@np.diag([1,1,np.sign(np.linalg.det(u@vt))])@vt
    return float(np.sqrt(np.mean(np.sum((b@rotation-a)**2,axis=1))))


def main():
    p=argparse.ArgumentParser();p.add_argument('--backend',choices=['cpu','gpu'],default='cpu')
    p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    gid='anthraquinone';initial='neu';final='red1'
    seed=DFT/gid/initial/'opt.xyz';global_final=DFT/gid/final/'opt.xyz'
    smiles=next(r['smiles'] for r in read_manifest() if r['id']==gid)
    final_info=json.loads((DFT/gid/final/'result.json').read_text())
    OUT.mkdir(parents=True,exist_ok=True)
    constraint=OUT/'anthraquinone_linker.freeze'
    n=write_linker_constraints(smiles,seed,constraint)
    identity=dict(id=gid,initial=initial,final=final,charge=final_info['charge'],mult=final_info['mult'],
        initial_geometry_sha1=geom_sha1(seed),global_final_geometry_sha1=geom_sha1(global_final),
        constraints_sha256=hashlib.sha256(constraint.read_bytes()).hexdigest(),
        worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        protocol=ACTIVE_SP,n_frozen_linker_torsions=n,max_opt_steps=80,
        assert_convergence=True,backend=a.backend)
    (OUT/'plan.json').write_text(json.dumps(identity,indent=2))
    if a.dry_run:
        print(json.dumps(identity));return
    def terminate(signum,frame):raise KeyboardInterrupt(f'termination signal {signum}')
    signal.signal(signal.SIGTERM,terminate)
    for mode in ['local_unconstrained','linker_fixed']:
        target=OUT/mode;target.mkdir(exist_ok=True)
        recfile=target/'result.json'
        if recfile.exists():
            old=json.loads(recfile.read_text())
            if old.get('identity')==identity and old.get('converged_smd') is True:
                continue
            raise RuntimeError(f'Existing unverified diagnostic must be preserved: {recfile}')
        if (target/'opt.xyz').exists():
            raise RuntimeError(f'Partial geometry exists; inspect and archive before restart: {target}')
        claim=target/'claim';claim.mkdir()
        (claim/'owner.json').write_text(json.dumps(dict(host=os.uname().nodename,pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().split()[21]),indent=2))
        try:
            rec=dft_smd(seed,int(final_info['charge']),int(final_info['mult']),do_opt=True,
                do_gas=False,do_smd=True,do_freq=False,opt_out=target/'opt.xyz',backend=a.backend,
                constraints=str(constraint) if mode=='linker_fixed' else None,
                assert_convergence=True,max_opt_steps=80)
            if rec.get('converged_smd') is not True:
                raise RuntimeError('Production-level SMD SCF did not converge')
            original=read_result(gid,final)
            optr=read(str(target/'opt.xyz'));seedr=read(str(seed));globalr=read(str(global_final))
            mask=np.asarray(optr.get_chemical_symbols())!='H'
            rec.update(identity=identity,validation_status='diagnostic_only',
                optimized_geometry_sha1=geom_sha1(target/'opt.xyz'),
                relative_to_original_final_active_eV=rec['e_smd_eV']-original['e_smd_eV'],
                rmsd_to_neutral_A=rmsd(optr.positions[mask],seedr.positions[mask]),
                rmsd_to_global_red1_A=rmsd(optr.positions[mask],globalr.positions[mask]),
                completed=time.time(),interpretation='local/fixed-linker sensitivity, no polymer assertion')
            tmp=recfile.with_suffix('.tmp');tmp.write_text(json.dumps(rec,indent=2));tmp.replace(recfile)
            print(json.dumps({k:rec[k] for k in ['validation_status','relative_to_original_final_active_eV',
                  'rmsd_to_neutral_A','rmsd_to_global_red1_A']}),flush=True)
        except Exception as exc:
            (target/'failure.json').write_text(json.dumps(dict(identity=identity,error=repr(exc)),indent=2))
            raise
        finally:
            (claim/'owner.json').unlink();claim.rmdir()


if __name__=='__main__':main()
