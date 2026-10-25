#!/usr/bin/env python
"""Independent fixed-geometry review of the experimental smooth CDS study.

This binds the actual CPU legacy library and recomputes representative values.
Reaction-level numbers are component energy sensitivities, never grafting G.
"""
import argparse
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    from ase.io import read
    from pyscf import gto
    from pyscf.solvent import smd,smd_experiment
    root=args.root.resolve();source=root/'summary.json';record=json.loads(source.read_text())
    api.verify_sources(record['source_sha256'])
    if record['status']!='experimental_smooth_CDS_component_study_complete':
        raise ValueError('Completed component study required')
    xyz=next(Path(p) for p in record['source_sha256'] if p.endswith('/rest/opt.xyz') and 'CDS_geometry_repair_interruption' in p)
    atoms=read(xyz);mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),
        unit='Angstrom',basis='def2-svp',charge=0,spin=0,verbose=0)
    obj=SimpleNamespace(mol=mol,solvent='acetonitrile',_solvent='acetonitrile',
        solvent_descriptors=None,sasa_ng=2030)
    legacy=float(smd.get_cds_legacy(obj)[0]);smooth=float(smd_experiment.get_cds(obj))
    stored=record['states']['contact'];expected=stored['grid_points'][-1]['smooth_E_Ha']
    if abs(smooth-expected)>1e-12 or abs(legacy-stored['legacy_E_Ha'])>1e-12:
        raise ValueError('Independent component energy replay differs')
    states=record['states'];sensitivities=[]
    for site in ('site2','site3'):
        for grid in (590,974,2030):
            def change(key):
                row=states[key];new=next(v['smooth_E_Ha'] for v in row['grid_points'] if v['grid']==grid)
                return new-row['legacy_E_Ha']
            free=sum(change(k)*v for k,v in {site:1,'empty':-1,'precursor':-1,'base':-1,'protonated':1,'chloride':1}.items())
            pair=sum(change(k)*v for k,v in {site:1,'empty':-1,'precursor':-1,'base':-1,'pair':1}.items())
            penalty=sum(change(k)*v for k,v in {site:1,'local':1,'empty':-1,'capped':-1}.items())
            sensitivities.append(dict(site=site,CDS_grid=grid,
                free_ion_reaction_component_shift_kJ_mol=free*2625.4996394799,
                paired_reaction_component_shift_kJ_mol=pair*2625.4996394799,
                polymer_minus_local_component_shift_kJ_mol=penalty*2625.4996394799))
    library=Path(smd.libsolvent._name)
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),source,library,
        Path(smd.__file__),Path(smd_experiment.__file__),xyz)}
    result=dict(status='smooth_CDS_fixed_geometry_component_review_passed',checked=time.time(),
        actual_legacy_library=str(library),legacy_energy_replay_error_Ha=abs(legacy-stored['legacy_E_Ha']),
        smooth_energy_replay_error_Ha=abs(smooth-expected),fixed_geometry_component_sensitivities=sensitivities,
        source_sha256=sources,full_DFT_electronic_energy_recomputed=False,nuclear_G_recomputed=False,
        original_geometry_changed=False,primary_records_changed=False,
        scope='CDS component sensitivity at fixed geometries only; geometry and nuclear terms must be recomputed for any alternate-model G')
    atomic_write_json(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],finest_grid=[r for r in sensitivities if r['CDS_grid']==2030])))


if __name__=='__main__':
    main()
