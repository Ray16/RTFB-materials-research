#!/usr/bin/env python3
"""Export the corrected primary cohort using the existing frozen CPU replay."""
import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import shutil

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def export(original,cohort,reviews,output):
    if output.exists():raise FileExistsError('Preserve earlier primary molecular exports')
    inputs=json.loads((cohort/'scalar_inputs.json').read_text());cycles=json.loads((cohort/'cycles.json').read_text())
    manifest=deepcopy(inputs['package']);outstates=[];sources={};output.mkdir(parents=True)
    def copy(path,destination):
        path=Path(path);target=output/destination;target.parent.mkdir(parents=True,exist_ok=True)
        digest=sha256_file(path);shutil.copyfile(path,target)
        if sha256_file(target)!=digest or sha256_file(path)!=digest:raise ValueError('Exported input bytes differ')
        sources[str(path)]=digest;return str(destination)
    for state in manifest['states']:
        key=state['id'];destination=Path('states')/key;files={}
        component=inputs['components'][key]
        if len(component['symbols'])==1:
            for role in ('geometry','result','partition'):files[role]=copy(original/state['files'][role],destination/Path(state['files'][role]).name)
            bound=component['sources'];energy_path=next(Path(p) for p in bound if p.endswith('/corrected/result.json'))
        else:
            report_path=reviews/(key+'.json');report=json.loads(report_path.read_text())
            if (report['key']!=key or not report['native_component_requalified']
                    or report['actual_SCF_wavefunction_spins_checked']!=9
                    or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
                raise ValueError('Exact checked primary molecular report required')
            spec=json.loads(Path(report['namespace']).read_text());case=Path(spec['case']);produced=json.loads((case/'result.json').read_text())
            native_path=Path(produced['native_record_path']);native=json.loads(native_path.read_text())
            partitions=[p for p in (case/'state/projected_thermal').glob('*.json') if json.loads(p.read_text()).get('identity',{}).get('raw_native_sha256')==sha256_file(native_path)]
            if len(partitions)!=1:raise ValueError('One bound projected partition required')
            paths=dict(geometry=case/'state/opt.xyz',result=case/'state/result.json',native=native_path,
                partition=partitions[0],hessian=native_path.with_suffix('')/'hessian_cartesian.npy')
            for role,path in paths.items():files[role]=copy(path,destination/(role+path.suffix))
            points=sorted((case/'SCF_evaluations').glob('*/*/settings.json'))
            if len(points)!=9:raise ValueError('All nine saved SCF/gradient points required')
            settings=[json.loads(p.read_text()) for p in points]
            coords=np.asarray([p['coordinates_Bohr'] for p in settings]);gradients=np.asarray([np.load(p.with_name('gradient.npy'),allow_pickle=False) for p in points])
            checks=[]
            for check in native['derivative_checks']:
                item=dict(check)
                for label,sign in (('minus',-1),('plus',1)):
                    target=coords[0].copy();target[check['atom_index'],check['axis']]+=sign*check['step_Bohr']
                    matches=np.flatnonzero(np.max(abs(coords-target),axis=(1,2))<3e-9)
                    if len(matches)!=1:raise ValueError('Exact unique displaced gradient required')
                    item[label+'_index']=int(matches[0])
                checks.append(item)
            metadata=[dict(settings=s,SCF_result=json.loads(p.with_name('SCF_result.json').read_text())) for p,s in zip(points,settings)]
            target=output/destination/'fd_arrays.npz';np.savez_compressed(target,coordinates_Bohr=coords,gradient_Ha_Bohr=gradients);files['fd_arrays']=str(target.relative_to(output))
            target=target.with_name('fd_records.json');atomic_write_json(target,dict(checks=checks,points=metadata));files['fd_records']=str(target.relative_to(output))
            energy_path=Path(spec['sources'])/'corrected_electronic.json';sources[str(report_path)]=sha256_file(report_path)
        energy=json.loads(energy_path.read_text());energy['protocol']=energy['energy_protocol'];energy['geom_sha1']=energy['energy_geom_sha1']
        target=output/destination/'energy.json';atomic_write_json(target,energy);files['energy']=str(target.relative_to(output));sources[str(energy_path)]=sha256_file(energy_path)
        outstates.append(dict(state,files=files,checked_G_by_cutoff_eV=component['G_by_cutoff_eV'],finite_difference_replay=len(component['symbols'])>1))
    for name in ('projected_modes.py','nuclear_equations.py','requirements.lock.txt'):copy(original/name,Path(name))
    source=original/'replay.py';body=source.read_text();old="if len(values)!=10 or fd_states!=6:raise ValueError('Ten exact states, including six saved-gradient panels, required')"
    if body.count(old)!=1:raise ValueError('Existing frozen replay scope changed; review the new exporter')
    body=body.replace(old,"if len(values)!=10 or fd_states!=9:raise ValueError('Ten exact states, including nine saved-gradient panels, required')")
    (output/'replay.py').write_text(body);copy(Path(__file__),Path('exporter_source.py'));sources[str(source)]=sha256_file(source)
    expected={r['id']:r['corrected_by_cutoff'] for r in cycles['cycles']}
    for reaction in manifest['reactions']:reaction['expected_by_cutoff']=expected[reaction['id']]
    manifest.update(states=outstates,source_sha256=sources,scope_id='corrected_primary_native_viologen_cohort_2026-11-02',
                    replay_scope='Nine full saved-gradient panels; cached modes/partitions/cycles only; no SCF or orbital-spin reconstruction')
    (output/'README.md').write_text('''# Corrected primary exact-viologen CPU replay

Run `python replay.py replay --package .` with the hash-locked NumPy dependency.
The ten-state cohort contains nine newly repeated primary molecular Hessians,
all nine saved-gradient panels and analytic monatomic chloride. It reconstructs
two Fc-referenced potentials, three precursor-specific quaternizations and
disproportionation at all three rotor cutoffs. The runner reuses the original
frozen example implementation, changing only its required FD-panel count to nine.

Only cached numerical data are read. No SCF/GPU/xTB calculation, original project
path or installed redox package is required. Actual SCF orbitals, resource review
and independent spin reconstructions have separate source-bound receipts/raw
storage; they are not rerun by this NumPy cache. Code, modes, partition equations,
geometry and energy/protocol records are included. No experimental/speciation,
conformer-ensemble or material admission or full historical executable closure
is inferred. The original six-panel example remains preserved separately.
''')
    manifest['files_sha256']={str(p.relative_to(output)):sha256_file(p) for p in output.rglob('*') if p.is_file()}
    atomic_write_json(output/'manifest.json',manifest);return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('original','cohort','reviews','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=export(a.original.resolve(),a.cohort.resolve(),a.reviews.resolve(),a.out.resolve())
    print(json.dumps(dict(states=len(r['states']),files=len(r['files_sha256']),scope=r['scope_id'])))
