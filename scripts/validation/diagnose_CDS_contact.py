#!/usr/bin/env python
"""CPU-only scan of legacy CDS energy at a suspected carbon-sphere contact.

This does not repair a geometry or qualify a Hessian/free energy. It records
one-sided energy slopes and finite-step curvature without assuming smoothness.
"""
import argparse
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',type=Path,required=True)
    parser.add_argument('--atoms',type=int,nargs=2,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    import numpy as np
    from ase.io import read
    from pyscf import gto
    from gpu4pyscf.solvent import smd
    xyz=args.geometry.resolve();atoms=read(xyz);i,j=args.atoms
    if atoms[i].symbol!='C' or atoms[j].symbol!='C':
        raise ValueError('This declared contact test is scoped to a carbon pair')
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),unit='Angstrom',
        basis='def2-svp',charge=0,spin=0,verbose=0)
    original=mol.atom_coords();bohr=0.529177210903
    direction=original[j]-original[i];initial=float(np.linalg.norm(direction))*bohr
    direction/=np.linalg.norm(direction)
    obj=SimpleNamespace(mol=mol,solvent_descriptors=None,solvent='acetonitrile')
    points=[]
    offsets=[-.01,-.005,-.002,-.001,-.0005,-.0002,-.0001,-.00004,-.00002,-.00001,
        0.,.00001,.00002,.00004,.0001,.0002,.0005,.001,.002,.005,.01]
    for delta in offsets:
        coordinates=original.copy();coordinates[j]=coordinates[i]+direction*((4.2+delta)/bohr)
        mol.set_geom_(coordinates,unit='Bohr');energy,gradient=smd.get_cds_legacy(obj)
        points.append(dict(offset_A=delta,distance_A=4.2+delta,E_CDS_Ha=float(energy),
            legacy_radial_gradient_Ha_Bohr=float(np.dot(gradient[j],direction))))
    left=[r for r in points if -.000041<r['offset_A']<0]
    right=[r for r in points if 0<r['offset_A']<.000041]
    def slope(rows):
        x=np.array([r['offset_A']/bohr for r in rows]);y=np.array([r['E_CDS_Ha'] for r in rows])
        return float(np.dot(x-x.mean(),y-y.mean())/np.dot(x-x.mean(),x-x.mean()))
    slopes=dict(left_Ha_Bohr=slope(left),right_Ha_Bohr=slope(right))
    energies={r['offset_A']:r['E_CDS_Ha'] for r in points}
    curvatures={str(h):float((energies[h]-2*energies[0.]+energies[-h])/(h/bohr)**2)
        for h in (.001,.0005,.0002,.0001,.00004,.00002,.00001)}
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),xyz,Path(smd.__file__))}
    result=dict(status='legacy_CDS_carbon_contact_energy_scan_complete',checked=time.time(),
        source_geometry=str(xyz),atom_indices=[i,j],initial_distance_A=initial,
        tested_contact_distance_A=4.2,points=points,one_sided_slopes=slopes,
        slope_jump_Ha_Bohr=slopes['right_Ha_Bohr']-slopes['left_Ha_Bohr'],
        centered_curvature_Ha_Bohr2_by_displacement_A=curvatures,source_sha256=sources,
        full_DFT_energy_evaluated=False,geometry_or_native_G_qualification_changed=False,
        interpretation='Measured legacy-CDS component behavior only; differentiability and a replacement protocol require separate review')
    args.out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(args.out,result)
    print(json.dumps({k:result[k] for k in ('status','initial_distance_A','one_sided_slopes','slope_jump_Ha_Bohr')}))


if __name__=='__main__':
    main()
