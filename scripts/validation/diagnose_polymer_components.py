#!/usr/bin/env python3
"""Separate a saved polymer FD error into surface/dispersion and remainder.

Only the density-independent CDS or D4 energy/gradient is reevaluated on
CPU. No SCF, GPU, native free energy or scientific admission is produced.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows.archived_native_namespace import verify_spec


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def split_residual(plus,minus,full_column,cds_plus,cds_minus,cds_column,step):
    """Compare gradient differences before dividing by the FD displacement."""
    arrays=[np.asarray(a,dtype=float).ravel() for a in (plus,minus,full_column,cds_plus,cds_minus,cds_column)]
    if len({a.shape for a in arrays})!=1 or not all(np.isfinite(a).all() for a in arrays) or not np.isfinite(step) or step<=0:
        raise ValueError('Finite matching derivative arrays and a positive step required')
    plus,minus,full_column,cds_plus,cds_minus,cds_column=arrays
    full=plus-minus-2*step*full_column
    cds=cds_plus-cds_minus-2*step*cds_column
    remainder=full-cds
    denominator=np.linalg.norm((plus-minus)/(2*step))
    return dict(full_gradient_difference_residual_norm=float(np.linalg.norm(full)),
        CDS_gradient_difference_residual_norm=float(np.linalg.norm(cds)),
        remaining_gradient_difference_residual_norm=float(np.linalg.norm(remainder)),
        full_relative_error=float(np.linalg.norm(full)/(2*step)/max(denominator,1e-14)),
        remaining_error_over_full_residual=float(np.linalg.norm(remainder)/max(np.linalg.norm(full),1e-30))),full,cds,remainder


def diagnose(namespace,output,component='CDS'):
    if output.exists():raise FileExistsError('Preserve prior component diagnostics')
    spec=verify_spec(namespace);reader=NamespaceSources(spec);case=Path(spec['attempt'])
    plan=json.loads((case/'plan.json').read_text());reader.verify_all(plan['producer_sha256']);reader.verify_all(plan['input_sha256'])
    if len(plan['tasks'])!=1 or (plan['tasks'][0]['charge'],plan['tasks'][0]['mult'])!=(0,1):
        raise ValueError('This diagnostic is scoped to one neutral singlet')
    if (plan['experimental_model']['CDS_points_per_sphere']!=2030
            or plan['experimental_model']['solvent']!='acetonitrile'
            or plan['experimental_native_settings']['FD_steps_Bohr']!=[.0005,.001]):
        raise ValueError('Expected unchanged experimental MeCN model and FD steps')
    status_path=case/'experimental_native/status.json';status=json.loads(status_path.read_text())
    if status['status']!='experimental_native_failed_preserved' or 'finished' not in status:
        raise ValueError('Only a closed failed native attempt may be diagnosed')
    raw=reader.path(spec['original_raw_root']);checks_path=raw/'derivative_checks.json'
    checks=json.loads(checks_path.read_text())['checks'];hessian_path=raw/'hessian.npy';hessian=np.load(hessian_path,allow_pickle=False)
    points=[]
    for path in sorted((raw/'SCF_evaluations').glob('*/evaluation.json')):
        row=json.loads(path.read_text())
        if not row['SCF_converged']:raise ValueError('A displaced SCF was not converged')
        points.append((path,np.asarray(row['coordinates_Bohr']),np.asarray(row['gradient_Ha_Bohr'])))
    if not points or not checks:raise ValueError('Saved base/displaced gradients and derivative checks required')
    coordinates=points[0][1]
    def point(target):
        matches=[p for p in points if np.max(np.abs(p[1]-target))<1e-9]
        if len(matches)!=1:raise ValueError('One exact saved SCF at each displacement required')
        return matches[0]
    geometry=raw/'geometry.xyz';symbols=[s.split()[0] for s in geometry.read_text().splitlines()[2:]]
    if coordinates.shape!=(len(symbols),3) or hessian.shape!=(coordinates.size,coordinates.size):raise ValueError('Geometry/derivative dimensions differ')
    output.mkdir(parents=True);sources=output/'sources';sources.mkdir();shutil.copyfile(__file__,sources/'diagnostic.py')
    if component not in ('CDS','D4'):raise ValueError('Choose an explicit density-independent component')
    qualifier=reader.path(Path(spec['prior_root'])/'scripts/validation/qualify_smooth_cds.py')
    api=load(qualifier,'diagnostic_original_CDS_derivative_loader')
    from pyscf.solvent import smd_experiment as smooth
    from pyscf import gto
    gradient,gradient_source=api.gradient_functions(smooth,sources/'experimental_gradient_source.py')
    profile=json.loads((case/'model/model.json').read_text())['profile']
    reader.verify_all(profile['implementation_sha256'])
    for path in (Path(smooth.__file__),gradient_source):
        if profile['implementation_sha256'].get(str(path))!=sha256_file(path):raise ValueError('Installed CDS implementation differs from the original attempt')
    molecule=gto.M(atom=list(zip(symbols,coordinates)),unit='Bohr',basis='sto-3g',charge=0,spin=0,verbose=0)
    obj=SimpleNamespace(mol=molecule,solvent='acetonitrile',_solvent='acetonitrile',solvent_descriptors=None,sasa_ng=2030)
    if component=='D4':
        from pyscf.dispersion import dftd4
        from pyscf.scf.dispersion import parse_disp
        method,version,atm=parse_disp('r2scan','d4')
        if (method,version,atm)!=('r2scan','d4',True):raise ValueError('D4 parameters differ from the native r2SCAN-D4 definition')
    cache={}
    def evaluate(coords):
        key=coords.tobytes()
        if key not in cache:
            index=len(cache);molecule.set_geom_(coords,unit='Bohr')
            if component=='CDS':energy=float(smooth.get_cds(obj));derivative=np.asarray(gradient(obj))
            else:
                value=dftd4.DFTD4Dispersion(molecule,xc=method,atm=atm).get_dispersion(grad=True)
                energy=float(value['energy']);derivative=np.asarray(value['gradient'])
            if not np.isfinite(energy) or derivative.shape!=coordinates.shape or not np.isfinite(derivative).all():raise ValueError('Invalid CPU CDS energy/gradient')
            path=output/f'{component}_{index:03d}.npz';np.savez(path,coordinates_Bohr=coords,gradient_Ha_Bohr=derivative,energy_Ha=energy)
            cache[key]=(energy,derivative,path)
            atomic_write_json(output/'progress.json',dict(component=component,component_evaluations=len(cache),complete=False,SCF_computations=0,GPU_jobs=0))
        return cache[key]
    evaluate(coordinates);rows=[];fine={};residuals={}
    for row in checks:
        atom,axis,step=row['atom'],row['axis'],row['step_Bohr'];selector=(atom,axis)
        if selector not in fine:
            grads=[]
            for sign in (-1,1):
                target=coordinates.copy();target[atom,axis]+=sign*1e-5;grads.append(evaluate(target)[1])
            fine[selector]=((grads[1]-grads[0])/(2e-5)).ravel()
        plus=coordinates.copy();minus=coordinates.copy();plus[atom,axis]+=step;minus[atom,axis]-=step
        pp,pm=point(plus),point(minus)
        values,full,cds,remainder=split_residual(pp[2],pm[2],hessian[:,3*atom+axis],evaluate(pp[1])[1],evaluate(pm[1])[1],fine[selector],step)
        if abs(values['full_relative_error']-row['relative_error'])>1e-10:raise ValueError('Saved full derivative error does not reconstruct')
        rows.append(dict(atom=atom,axis=axis,step_Bohr=step,**values))
        label=f'{atom}_{axis}_{step}';residuals[label+'_full']=full;residuals[label+'_'+component]=cds;residuals[label+'_remaining']=remainder
        atomic_write_json(output/'component_checks.json',dict(checks=rows,complete=False,full_native_G_qualified=False))
    np.savez(output/'residuals.npz',**residuals)
    paths=[namespace,case/'plan.json',status_path,checks_path,hessian_path,geometry,qualifier,Path(smooth.__file__),gradient_source,
        case/'model/model.json',*sorted(sources.glob('*')),*sorted(output.glob('*.npz')),*[p[0] for p in points]]
    hashes={str(p):sha256_file(p) for p in paths};reader.recheck()
    if component=='D4':
        for path in (Path(dftd4.__file__),Path(dftd4.libdftd4._name)):hashes[str(path)]=sha256_file(path)
        for row in rows:row['D4_gradient_difference_residual_norm']=row.pop('CDS_gradient_difference_residual_norm')
    result=dict(status='saved_polymer_FD_error_split_into_component_and_remainder',key=spec['key'],component=component,
        method_id=spec['method_id'],checks=rows,component_evaluations=len(cache),component_Hessian_column_step_Bohr=1e-5,
        CPU_affinity=sorted(os.sched_getaffinity(0)),quantum_SCF_computations=0,GPU_jobs=0,
        full_native_G_qualified=False,primary_candidate_admission=False,source_sha256=hashes)
    atomic_write_json(output/'summary.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--component',choices=('CDS','D4'),default='CDS');a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(1))
    r=diagnose(a.spec.resolve(),a.out.resolve(),a.component);print(json.dumps(dict(status=r['status'],component=a.component,checks=r['checks'])))
