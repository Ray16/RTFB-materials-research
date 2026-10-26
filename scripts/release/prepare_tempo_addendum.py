#!/usr/bin/env python
"""Stage the four source-bound TEMPO/Fc states for isolated cached replay."""
import argparse
import ast
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT,DFT,read_manifest
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path
from redox.qm.native_thermal import cache_path,load_native_thermal
from redox.workflows import polymer_campaign as api


def prepare(source,out):
    report=json.loads(source.read_text());api.verify_sources(report['source_sha256'])
    if (report['status']!='TEMPO_local_oxidation_native_partition_independently_reconstructed'
            or report['experimental_exact_graft_potential_claimed'] or report['free_energy_of_grafting_computed']
            or report['all_axis_eligibility']):
        raise ValueError('Expected checked local-model result without promoted claims')
    original={(r['id'],r['state']):r for r in read_manifest()}
    definitions=[('tempo_rad','tempo','rad','tempo__rad'),('tempo_ox','tempo','ox','tempo__ox'),
                 ('fc_neutral','ferrocene','neu','ferrocene__neu'),('fc_oxidized','ferrocene','ox','ferrocene__ox')]
    if out.exists():raise FileExistsError('Preserve existing package')
    out.mkdir(parents=True);sources={};exported=[]

    def copy(path,target):
        path=Path(path);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
        if api.digest(path)!=api.digest(target):raise ValueError('Copied source bytes differ')
        sources[str(path)]=api.digest(path)

    copy(ROOT/'examples/cached_thermochemistry/reproduce.py',out/'cached_reader.py')
    copy(ROOT/'examples/tempo_replay/reproduce.py',out/'reproduce.py')
    copy(ROOT/'examples/control_replay/requirements-cpu.lock.txt',out/'requirements-cpu.lock.txt')
    copy(Path(__file__).resolve(),out/'prepare_source.py')
    copy(source,out/'tempo_source.json')
    derived={}
    for script,function,filename,imports in [
        (ROOT/'scripts/validation/review_smooth_quantum.py','projected_modes','projected_modes.py','import numpy as np'),
        (ROOT/'scripts/validation/audit/native_partition.py','nuclear_equations','nuclear_equations.py','import math\nimport numpy as np'),
        (ROOT/'examples/neighbor_cycle_replay/reproduce.py','native_replay','native_engine.py','import importlib.util\nimport json\nimport math')]:
        text=script.read_text();tree=ast.parse(text)
        node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==function)
        (out/filename).write_text('"""Frozen CPU replay function; original source is identified in the manifest."""\n'+imports+'\n\n'+ast.get_source_segment(text,node)+'\n')
        sources[str(script)]=api.digest(script)
        derived[filename]=dict(source=str(script),source_sha256=api.digest(script),function=function)
    for name,gid,state,key in definitions:
        state_dir=DFT/gid/state;raw=json.loads((state_dir/'result.json').read_text())
        q,m=raw['charge'],raw['mult'];active=load_record(state_dir,q,m);native=load_native_thermal(state_dir,q,m)
        if active is None or native is None or (q,m)!=(report['states'][key]['charge'],report['states'][key]['multiplicity']):
            raise ValueError('State lacks a matching qualified electronic/native parent: '+name)
        total=active['e_smd_eV']+native['g_thermal_total_eV']
        if abs(total-report['states'][key]['G_eV'])>1e-8:raise ValueError('State differs from checked TEMPO result')
        native_path,_=cache_path(state_dir,q,m)
        paths=dict(geometry=state_dir/'opt.xyz',result=state_dir/'result.json',energy=record_path(state_dir,q,m),
                   native=native_path,partition=Path(native['thermal_postprocessing_path']),
                   hessian=native_path.with_suffix('')/'hessian_cartesian.npy')
        files={}
        for role,path in paths.items():
            target=out/name/(role+path.suffix);copy(path,target);files[role]=str(target.relative_to(out))
        exported.append(dict(id=name,smiles=original[gid,state]['smiles'],charge=q,mult=m,partition_kind='projected_DFT',
                             files=files,checked_G_eV=total,original_state_directory=str(state_dir)))
    (out/'README.md').write_text('''# TEMPO local-graft raw-state addendum

Run `python3 -I -S /path/to/addendum/reproduce.py --data /path/to/addendum`
to reconstruct the cached potential. With the locked NumPy environment, add
`--full-native` to reconstruct four Hessian mode sets and nuclear sums.
Original cluster paths are descriptive; replay reads only package-relative files.
This is not a re-run of SCF, stationarity, spin, reaction kinetics or experiment.
It does not establish grafting free energy, polymer interactions, conformer
convergence, material eligibility or a public archival DOI.
''')
    api.verify_sources(sources)
    data=dict(schema=1,temperature_K=298.15,states=exported,reactions=[],
        potential=dict(oxidized='tempo_ox',reduced='tempo_rad',reference_oxidized='fc_oxidized',reference_reduced='fc_neutral'),
        source_energy_conversion_kJ_mol_per_eV=96.48533212331002,SI_energy_conversion_kJ_mol_per_eV=96.48533212331002,
        source_review='tempo_source.json',source_sha256=sources,derived_functions=derived,
        experimental_accuracy_inferred=False,material_eligibility_inferred=False,grafting_deltaG_computed=False,
        conformer_convergence_established=False,created=time.time(),
        files_sha256={str(p.relative_to(out)):api.digest(p) for p in out.rglob('*') if p.is_file()})
    atomic_write_json(out/'manifest.json',data)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    r=prepare(a.source.resolve(),a.out.resolve());print(json.dumps({'status':'TEMPO_raw_addendum_prepared_replay_pending','states':len(r['states'])}))
