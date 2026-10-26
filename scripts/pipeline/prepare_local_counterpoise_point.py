#!/usr/bin/env python
"""Prepare a source-bound frozen motif counterpoise point after control review.

The resulting plan is gas-only, uses the existing five-term calculator, and
does not admit a candidate or modify a solvent interaction result.
"""
import argparse
import json
import math
import os
import re
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

CORE=ROOT/'scripts/validation/local_counterpoise.py'
WORKER=ROOT/'scripts/pipeline/run_local_counterpoise.py'
GUARDED_WORKER=ROOT/'scripts/pipeline/run_local_counterpoise_point.py'
PROFILE=ROOT/'environments/cutensor22/runtime.json'


def read_xyz(path):
    lines=path.read_text().splitlines()
    if len(lines)<3:raise ValueError('Missing XYZ geometry')
    n=int(lines[0]);atoms=[];coordinates=[]
    if len(lines)!=n+2:raise ValueError('XYZ atom count differs')
    for line in lines[2:]:
        cols=line.split()
        if len(cols)!=4 or not re.fullmatch(r'[A-Z][a-z]?',cols[0]):raise ValueError('Invalid XYZ atom')
        xyz=[float(x) for x in cols[1:]]
        if not all(math.isfinite(v) for v in xyz):raise ValueError('Nonfinite XYZ atom')
        atoms.append(cols[0]);coordinates.append(xyz)
    return atoms,coordinates


def prepare(construction,construction_review,control_review,key,out):
    if out.exists():raise FileExistsError('Preserve prior attempts')
    manifest_path=construction/'manifest.json';manifest=json.loads(manifest_path.read_text())
    api.verify_sources(manifest['source_sha256'])
    verified=json.loads(construction_review.read_text());api.verify_sources(verified['source_sha256'])
    if (verified['status']!='fixed_parent_local_motif_construction_independently_checked'
            or verified['primary_candidate_admission']
            or sum(r['key']==key and r['parent_coordinates_preserved'] for r in verified['models'])!=1):
        raise ValueError('Independently checked construction for this exact case is required')
    control=json.loads(control_review.read_text());api.verify_sources(control['source_sha256'])
    if control['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':
        raise ValueError('A checked CPU/GPU gas counterpoise control is required')
    if control['solvent_counterpoise_qualified'] or control['primary_candidate_admission']:
        raise ValueError('Gas control has an incompatible scope')
    rows=[r for r in manifest['entries'] if r['key']==key]
    if len(rows)!=1 or not re.fullmatch(r'[a-zA-Z0-9_]+',key):raise ValueError('Unique safe declared case key required')
    entry=rows[0];molecules={r['id']:r for r in entry['molecules']}
    if set(molecules)!={'A','B','AB'} or any(r['charge']!=0 for r in molecules.values()):
        raise ValueError('Exact neutral two-fragment construction required')
    geometries={}
    for label,row in molecules.items():
        path=Path(row['geometry'])
        if api.digest(path)!=row['geometry_sha256']:raise ValueError('Frozen construction coordinates changed')
        geometries[label]=read_xyz(path)
    a,b,ab=(geometries[k] for k in ('A','B','AB'))
    if a[0]+b[0]!=ab[0] or a[1]+b[1]!=ab[1] or len(ab[0])!=38 or len(a[0])!=19:
        raise ValueError('Ordered AB must equal the two exact 19-atom frozen fragments')
    core=api.load_module('motif_counterpoise_protocol',CORE)
    if control['settings']!=core.SETTINGS:raise ValueError('Counterpoise protocol differs from qualified control')
    # Prospective unknown 38-center ghost jobs require a >=32-GB device. The
    # worker guard retains a 30-GiB observed peak ceiling on that device.
    out.mkdir(parents=True);case=out/key;case.mkdir()
    geometry=dict(symbols=ab[0],coordinates_A=ab[1],A_indices=list(range(len(a[0]))))
    atomic_write_json(case/'geometry.json',geometry)
    files=[WORKER,GUARDED_WORKER,Path(__file__).resolve(),CORE,ROOT/'scripts/pipeline/local_motif_energy.py',ROOT/'scripts/pipeline/run_pi_stacking.py',
        ROOT/'scripts/fleet/exclusive_python_worker.py',ROOT/'scripts/fleet/tensor_runtime_worker.py',
        ROOT/'scripts/fleet/bounded_host_python_worker.py']
    inputs=[case/'geometry.json',PROFILE,control_review,construction_review,manifest_path]+[Path(v['geometry']) for v in molecules.values()]
    plan=dict(task=dict(key=key,atoms=38,charge=0,mult=1,role='frozen_local_motif_gas_counterpoise'),
        settings=core.SETTINGS,input_geometry=str(case/'geometry.json'),
        input_sha256={str(p):api.digest(p) for p in inputs},producer_sha256={str(p):api.digest(p) for p in files},
        runtime_profile=str(PROFILE),execution_driver=str(GUARDED_WORKER),maximum_host_RSS_GiB=32,
        minimum_GPU_MiB=32768,maximum_GPU_MiB=30720,
        raw_root=str(ROOT/'calcs/local_motif_counterpoise'/out.name/key),
        control_review=str(control_review),source_construction_key=key,
        primary_candidate_admission=False,solvent_counterpoise_computed=False)
    atomic_write_json(case/'plan.json',plan)
    atomic_write_json(out/'campaign.json',dict(status='prepared_gas_motif_counterpoise_no_scientific_start',
        entries=[dict(key=key,campaign=str(case),plan_sha256=api.digest(case/'plan.json'))],
        gas_control_review_sha256=api.digest(control_review),minimum_GPU_MiB=32768,
        preparer_sha256=api.digest(__file__),source_sha256={str(p):api.digest(p) for p in files+inputs},
        primary_candidate_admission=False,solvent_counterpoise_computed=False))
    return case


def launch_32gb(case,gpu,preflight):
    from redox.workflows.resources import free_slots
    import socket
    launcher=api.load_module('motif_counterpoise_existing_reservation_worker',ROOT/'scripts/pipeline/local_motif_energy.py')
    launcher.__file__=str(GUARDED_WORKER)
    device_module=api.load_module('motif_counterpoise_device_probe',ROOT/'scripts/pipeline/local_motif_energy.py')
    actual=device_module.device(gpu);host=socket.gethostname().split('.')[0]
    if actual['memory_MiB']<32768:raise ValueError('Prospective ghost calculation requires a 32-GB device')
    if (host,gpu) not in free_slots([host]):return dict(status='atomic_placement_unavailable')
    return launcher.launch(case,gpu,preflight)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','launch'])
    p.add_argument('--construction',type=Path);p.add_argument('--construction-review',type=Path)
    p.add_argument('--control-review',type=Path);p.add_argument('--key');p.add_argument('--out',type=Path)
    p.add_argument('--case',type=Path);p.add_argument('--gpu',type=int);p.add_argument('--preflight',type=Path)
    a=p.parse_args()
    if a.action=='prepare':
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
        print(prepare(a.construction.resolve(),a.construction_review.resolve(),a.control_review.resolve(),a.key,a.out.resolve()))
    else:
        print(json.dumps(launch_32gb(a.case.resolve(),a.gpu,a.preflight.resolve())))


if __name__=='__main__':main()
