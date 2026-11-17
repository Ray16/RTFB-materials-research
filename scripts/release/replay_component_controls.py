#!/usr/bin/env python3
"""Replay a relocated electronic-control archive on CPU without original paths."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import types
import h5py


def sha(path):
    result=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):result.update(block)
    return result.hexdigest()


def local(root,name):
    p=Path(name)
    if p.is_absolute() or '..' in p.parts:raise ValueError('Package-relative paths required')
    result=(root/p).resolve()
    if not result.is_relative_to(root):raise ValueError('Package path escape')
    return result


def load(path,name):
    module=types.ModuleType(name);module.__file__=str(path)
    exec(compile(path.read_bytes(),str(path),'exec'),module.__dict__);return module


def replay(root):
    root=root.resolve();manifest=json.loads((root/'manifest.json').read_text())
    if any(p.is_symlink() for p in root.rglob('*')):raise ValueError('Package symlinks forbidden')
    files={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if files!=set(manifest['files_sha256'])|{'manifest.json'}:raise ValueError('Complete package inventory required')
    for name,digest in manifest['files_sha256'].items():
        if sha(local(root,name))!=digest:raise ValueError('Package input changed')
    frame_checker=load(local(root,'code/scf_frame.py'),'portable_actual_scf_frame')
    # Supply only the recorded conversion constants to the unchanged analysis
    # source. No installed redox package or original checkout is imported.
    common=types.ModuleType('redox.core.common')
    common.HARTREE_EV=manifest['hartree_eV'];common.EV_KJ=manifest['eV_kJ_mol']
    sys.modules['redox']=types.ModuleType('redox');sys.modules['redox.core']=types.ModuleType('redox.core');sys.modules['redox.core.common']=common
    analysis=load(local(root,'code/component_sensitivity.py'),'portable_component_sensitivity')
    values={};nuclei={};frames=0
    for key,row in manifest['cases'].items():
        folder=local(root,row['folder']);plan=json.loads((folder/'case/plan.json').read_text())
        xyz=(folder/'sources/input.xyz').read_text().splitlines();atoms=[line.split() for line in xyz[2:]]
        if len(atoms)!=int(xyz[0]) or any(len(a)!=4 for a in atoms):raise ValueError('Complete input XYZ required')
        symbols=[a[0] for a in atoms];coords=[[float(v)/manifest['bohr_A'] for v in a[1:]] for a in atoms]
        nuclei[key]=Counter(symbols);state=plan['source_state'];energies={}
        for variant in ('installed','corrected'):
            result=json.loads((folder/'case'/variant/'result.json').read_text());solved=[]
            if (result['energy_protocol']!=manifest['electronic_protocol'] or result['converged_smd'] is not True
                    or (result['charge'],result['mult'])!=(state['charge'],state['mult'])
                    or not math.isfinite(result['e_smd_Ha']) or abs(result['smd_spin_contam'])>.1):
                raise ValueError('Matched converged method/state and producer spin gate required')
            for checkpoint in sorted((folder/'case'/variant/'SCF_evaluations').glob('*/*/SCF.h5')):
                settings=json.loads(checkpoint.with_name('settings.json').read_text());record=json.loads(checkpoint.with_name('SCF_result.json').read_text())
                if (settings['charge'],settings['mult'])!=(state['charge'],state['mult']):continue
                if (settings['basis']!=manifest['electronic_protocol']['basis'] or settings['xc']!=manifest['electronic_protocol']['xc']
                        or settings['solvent_class']!='SMD' or settings['conv_tol']!=1e-9):raise ValueError('SCF method declaration differs')
                frame_checker.checkpoint_frame(checkpoint,charge=state['charge'],mult=state['mult'],basis=manifest['electronic_protocol']['basis'],coordinates_Bohr=coords,symbols=symbols)
                with h5py.File(checkpoint,'r') as h:
                    energy=float(h['scf/e_tot'][()])
                if not math.isfinite(energy) or abs(energy-record['energy_Ha'])>1e-9:raise ValueError('Checkpoint energy differs')
                frames+=1
                if record['converged']:solved.append(energy)
            if not solved or abs(min(solved)-result['e_smd_Ha'])>1e-9:raise ValueError('Selected solution lacks checkpoint support')
            energies[variant]=result['e_smd_Ha']
        values[key]=dict(stored_Ha=plan['reference_energy_Ha'],installed_Ha=energies['installed'],corrected_Ha=energies['corrected'],charge=state['charge'],mult=state['mult'])
    expected=json.loads(local(root,'expected_replay.json').read_text());panels={}
    for name,row in manifest['panels'].items():
        package=json.loads(local(root,row['package']).read_text());mapping=row['case_keys']
        for reaction in package['reactions']:
            elements=set(e for state in reaction['coefficients'] for e in nuclei[mapping[state]])
            if any(sum(c*nuclei[mapping[k]][e] for k,c in reaction['coefficients'].items())!=0 for e in elements):
                raise ValueError('Actual reaction geometry inventories do not balance')
        panels[name]=analysis.reconstruct(package,{state:values[case] for state,case in mapping.items()})
    if panels!=expected:raise ValueError('Electronic contributions differ from reviewed source')
    return dict(status='relocated_component_checkpoints_and_electronic_shifts_replayed',unique_execution_namespaces=len(values),actual_checkpoint_frames=frames,
        panel_states={name:r['states'] for name,r in panels.items()},panel_cycles={name:len(r['cycles']) for name,r in panels.items()},
        nuclear_terms_recomputed=False,SCF_calculations_run=False,producer_spin_diagnostics_independently_recomputed=False,
        new_free_energies_qualified=False,public_deposition_complete=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('package',type=Path);a=p.parse_args();print(json.dumps(replay(a.package)))
