#!/usr/bin/env python
"""CPU-only local energy/gradient and scale diagnostics for first-site states."""
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
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--cached-scan',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    import numpy as np
    from ase.io import read
    from pyscf import gto
    from gpu4pyscf.solvent import smd
    from redox.qm.cds_finite import energy_gradient
    source=args.source.resolve();report=json.loads(source.read_text());api.verify_sources(report['source_sha256'])
    cached={r['role']:r for r in json.loads(args.cached_scan.read_text())['states']} if args.cached_scan else {}
    states={'empty':report['empty_reference']['state_dir'],'capped':report['local_capped_reference']['state_dir'],
        **{f'site{r["site"]}':r['state_dir'] for r in report['singles']}}
    sources={str(p):api.digest(p) for p in (source,Path(__file__).resolve(),Path(smd.__file__),ROOT/'src/redox/qm/cds_finite.py')}
    if args.cached_scan:sources[str(args.cached_scan.resolve())]=api.digest(args.cached_scan)
    rows=[]
    for role,directory in states.items():
        xyz=Path(directory)/'opt.xyz';atoms=read(xyz);sources[str(xyz)]=api.digest(xyz)
        mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),unit='Angstrom',basis='def2-svp',charge=0,spin=0,verbose=0)
        coords=mol.atom_coords();obj=SimpleNamespace(mol=mol,solvent_descriptors=None,solvent='acetonitrile')
        _,legacy=smd.get_cds_legacy(obj)
        def energy(x):obj.mol.set_geom_(x,unit='Bohr');return smd.get_cds_legacy(obj)[0]
        if role in cached:
            if cached[role]['geometry_sha256']!=api.digest(xyz):raise ValueError('Cached component scan has another geometry')
            values={k:np.array(v) for k,v in cached[role]['gradients_by_step_Ha_Bohr'].items()}
        else:values={str(h):energy_gradient(coords,energy,h) for h in (1e-4,5e-4,1e-3)}
        error=float(np.max(abs(np.asarray(legacy)-values['0.0001'])))
        rows.append(dict(role=role,state_dir=directory,geometry_sha256=api.digest(xyz),
            local_legacy_vs_energy_gradient_maximum_error_Ha_Bohr=error,
            scale_differences_Ha_Bohr={f'{a}:{b}':float(np.max(abs(values[str(a)]-values[str(b)])))
                for a,b in ((1e-4,5e-4),(5e-4,1e-3))},
            gradients_by_step_Ha_Bohr={k:v.tolist() for k,v in values.items()},
            local_match=bool(error<=1e-6)))
    api.verify_sources(sources)
    result=dict(status='first_site_local_CDS_gradients_checked_larger_scale_sensitivity_retained',checked=time.time(),
        states=rows,all_local_gradient_matches=all(r['local_match'] for r in rows),
        maximum_local_error_Ha_Bohr=max(r['local_legacy_vs_energy_gradient_maximum_error_Ha_Bohr'] for r in rows),
        cached_component_scan_reused=bool(args.cached_scan),source_sha256=sources,GPU_calculations_started=0,
        primary_native_admission_changed=False,global_thermal_or_conformer_convergence_inferred=False,
        scope='All Cartesian components of the legacy CDS term; not a replacement for full-energy or full-Hessian checks')
    args.out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(args.out,result)
    print(json.dumps({k:result[k] for k in ('status','all_local_gradient_matches','maximum_local_error_Ha_Bohr')}))


if __name__=='__main__':main()
