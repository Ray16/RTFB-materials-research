#!/usr/bin/env python
"""Independently check a prepared local-motif gas counterpoise geometry."""
import argparse
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def review(case,out):
    plan=api.read_plan(case);api.verify_sources(plan['input_sha256']);api.verify_sources(plan['producer_sha256'])
    if (plan['task']['atoms']!=38 or plan['task']['charge']!=0 or plan['task']['mult']!=1
            or plan['minimum_GPU_MiB']!=32768 or plan['primary_candidate_admission']
            or plan['solvent_counterpoise_computed']):raise ValueError('Wrong prospective point scope')
    control=json.loads(Path(plan['control_review']).read_text());api.verify_sources(control['source_sha256'])
    if control['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':raise ValueError('Unqualified numerical control')
    geom=json.loads(Path(plan['input_geometry']).read_text())
    manifests=[Path(p) for p in plan['input_sha256'] if Path(p).name=='manifest.json']
    if len(manifests)!=1:raise ValueError('Unique parent construction manifest required')
    manifest=json.loads(manifests[0].read_text());matches=[r for r in manifest['entries'] if r['key']==plan['source_construction_key']]
    if len(matches)!=1:raise ValueError('Parent construction key is not unique')
    entry=matches[0];molecules={r['id']:r for r in entry['molecules']}
    reader=api.load_module('review_motif_counterpoise_xyz',ROOT/'scripts/pipeline/prepare_local_counterpoise_point.py')
    a,b,pair=(reader.read_xyz(Path(molecules[k]['geometry'])) for k in ('A','B','AB'))
    if (a[0]+b[0]!=pair[0] or a[1]+b[1]!=pair[1] or geom['symbols']!=pair[0]
            or geom['coordinates_A']!=pair[1] or geom['A_indices']!=list(range(19))):
        raise ValueError('Frozen coordinates or fragment partition differ')
    core=api.load_module('review_prepared_counterpoise_pyscf',ROOT/'scripts/validation/local_counterpoise.py')
    mols=core.molecules(geom['symbols'],geom['coordinates_A'],geom['A_indices'])
    electrons={k:int(v.nelectron) for k,v in mols.items()};basis={k:int(v.nao_nr()) for k,v in mols.items()}
    if (electrons!={'AB':168,'A_own':84,'B_own':84,'A_ghost':84,'B_ghost':84}
            or basis['AB']!=1026 or basis['A_ghost']!=1026 or basis['B_ghost']!=1026):
        raise ValueError('Unexpected fragment electron or AO count')
    if (case/'launch.json').exists() or (case/'counterpoise/status.json').exists():
        raise ValueError('Prepared-only review requires an unstarted case')
    proof={**plan['input_sha256'],**plan['producer_sha256'],str(case/'plan.json'):api.digest(case/'plan.json'),
           str(Path(__file__).resolve()):api.digest(__file__)}
    result=dict(status='prepared_motif_counterpoise_geometry_independently_checked_no_scientific_start',
        checked=time.time(),campaign=str(case),parent_key=plan['source_construction_key'],
        atoms=38,electrons=electrons,basis_functions=basis,source_sha256=proof,
        solvent_counterpoise_computed=False,primary_candidate_admission=False)
    if out.exists():raise FileExistsError('Preserve existing review')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review(a.case.resolve(),a.out.resolve());print(json.dumps({k:result[k] for k in ('status','atoms','electrons','basis_functions')}))


if __name__=='__main__':main()
