#!/usr/bin/env python3
"""Reconstruct cached molecular redox/grafting quantities using NumPy only."""
import argparse
from collections import Counter
import hashlib
import json
import math
import numbers
from pathlib import Path
import types

import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local(root,relative):
    p=Path(relative)
    if p.is_absolute() or '..' in p.parts:raise ValueError('Package path must be relative without traversal')
    result=(root/p).resolve()
    if not result.is_relative_to(root.resolve()):raise ValueError('Package path escapes its root')
    return result


def module(path,name):
    # A replay must not add bytecode files to its hash-locked input inventory.
    result=types.ModuleType(name);result.__file__=str(path)
    exec(compile(path.read_bytes(),str(path),'exec'),result.__dict__);return result


def close(a,b,tolerance=1e-8):
    if any(isinstance(v,(bool,np.bool_)) or not isinstance(v,numbers.Real) or not math.isfinite(v) for v in (a,b)) or abs(a-b)>tolerance:
        raise ValueError('Cached quantity does not reconstruct')


def replay(root):
    root=Path(root).resolve();manifest=json.loads((root/'manifest.json').read_text())
    if any(p.is_symlink() for p in root.rglob('*')):raise ValueError('Package symlinks are not allowed')
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if actual!=set(manifest['files_sha256'])|{'manifest.json'}:raise ValueError('Package inventory differs')
    for relative,digest in manifest['files_sha256'].items():
        path=local(root,relative)
        if path.is_symlink() or sha(path)!=digest:raise ValueError('Package content differs: '+relative)
    if (manifest['material_admission'] or manifest['global_conformer_convergence_proven']
            or manifest['experimental_exact_graft_accuracy_validated'] or manifest['full_historical_source_closure_revalidated']):
        raise ValueError('Numerical replay cannot promote scientific scope')
    if manifest['product_inchikey']!='PBNXFHGCSKYNHC-UHFFFAOYSA-N':raise ValueError('This example requires its exact source product identity')
    if (manifest['active_protocol']['basis']!='def2-tzvpd' or manifest['active_protocol']['xc']!='wb97m-v'
            or manifest['active_protocol']['solvent']!='acetonitrile' or manifest['native_protocol']['xc']!='r2scan'
            or manifest['native_protocol']['basis']!='def2-svp' or manifest['native_protocol']['temperature_K']!=298.15):
        raise ValueError('Reported molecular protocols differ from this exact example')
    modes=module(local(root,'projected_modes.py'),'portable_projected_modes')
    equations=module(local(root,'nuclear_equations.py'),'portable_nuclear_equations')
    values={};inventories={};charges={};records=[];fd_states=0
    for row in manifest['states']:
        name=row['id']
        if name in values:raise ValueError('Duplicate cached state')
        load=lambda role:json.loads(local(root,row['files'][role]).read_text())
        geometry=local(root,row['files']['geometry']);lines=geometry.read_text().splitlines();atoms=[s.split() for s in lines[2:]]
        if len(atoms)!=int(lines[0]) or any(len(a)!=4 for a in atoms):raise ValueError('Malformed cached geometry')
        xyz=np.array([[float(v) for v in a[1:]] for a in atoms])
        if not np.isfinite(xyz).all():raise ValueError('Nonfinite cached coordinates')
        raw=load('result');energy=load('energy');partition=load('partition')
        if any((r['charge'],r['mult'])!=(row['charge'],row['mult']) for r in (raw,energy)):
            raise ValueError('Cached charge/spin differs')
        if (energy['geom_sha1']!=hashlib.sha1(geometry.read_bytes()).hexdigest()[:16]
                or energy['converged_smd'] is not True or energy['protocol']!=manifest['active_protocol']):
            raise ValueError('Electronic protocol, geometry or convergence differs')
        if row['mult']<1 or type(row['mult']) is not int or type(row['charge']) is not int:raise ValueError('Integral electronic state required')
        frequency_error=0.;fd_checks=[]
        if row['partition_kind']=='monatomic':
            if len(atoms)!=1 or partition['internal_modes']!=0 or row['mult']!=1:raise ValueError('Invalid monatomic state')
            mass=np.array([partition['mass_amu']]);frequencies=np.array([]);sigma=1
            nuclear={c:equations.nuclear_equations(frequencies,xyz,mass,sigma,int(c)) for c in ('25','50','100')}
            for value in nuclear.values():close(value,partition['G_thermal_eV'])
            spin=0.
        else:
            native=load('native');hessian=np.load(local(root,row['files']['hessian']),allow_pickle=False)
            if (hessian.shape!=(3*len(atoms),3*len(atoms)) or not np.isfinite(hessian).all()
                    or np.max(abs(hessian-hessian.T))>1e-6 or native['identity']['protocol']!=manifest['native_protocol']
                    or native['converged'] is not True or native['stationary_minimum'] is not True
                    or native['grid_response'] is not True or not 0<=native['gradient_max']<=.00045
                    or not math.isfinite(native['smd_spin_diagnostics']['spin_contam'])
                    or abs(native['smd_spin_diagnostics']['spin_contam'])>.1):raise ValueError('Cached native surface gate differs')
            if (len(native['derivative_checks'])!=4 or any(not math.isfinite(r['relative_error']) or not 0<=r['relative_error']<=.01 for r in native['derivative_checks'])):
                raise ValueError('Reported native derivative gate differs')
            for source in (native,partition):
                if (source['identity']['geometry_sha256']!=sha(geometry)
                        or (source['identity']['charge'],source['identity']['mult'])!=(row['charge'],row['mult'])):
                    raise ValueError('Cached nuclear identity differs')
            if partition['identity']['raw_native_sha256']!=sha(local(root,row['files']['native'])):raise ValueError('Partition/native binding differs')
            mass=np.asarray(partition['atomic_masses_amu']);sigma=partition['symmetry_number']
            if mass.shape!=(len(atoms),) or np.any(mass<=0) or not np.isfinite(mass).all() or type(sigma) is not int or sigma<1:raise ValueError('Invalid masses or symmetry number')
            frequencies,rank=modes.projected_modes(hessian,xyz,mass)
            expected=np.asarray(partition['frequencies_cm'])
            if expected.shape!=frequencies.shape or not np.isfinite(expected).all() or rank!=partition['rigid_rank'] or min(frequencies)<-5:raise ValueError('Internal mode coverage or minimum differs')
            frequency_error=float(np.max(abs(frequencies-expected)))
            if frequency_error>1e-5:raise ValueError('Hessian modes do not reconstruct')
            nuclear={c:equations.nuclear_equations(frequencies,xyz,mass,sigma,int(c)) for c in ('25','50','100')}
            for c,value in nuclear.items():close(value,partition['nuclear_by_cutoff_eV'][c])
            spin=-8.617333262145e-5*298.15*math.log(row['mult']);close(spin,partition['g_spin_partition_eV'])
            # Scope explicitly distinguishes new saved-gradient reconstructions
            # from inherited reference records whose FD metadata alone survive.
            if row.get('finite_difference_replay'):
                saved=np.load(local(root,row['files']['fd_arrays']),allow_pickle=False);fd=load('fd_records')
                coords=saved['coordinates_Bohr'];gradients=saved['gradient_Ha_Bohr']
                if coords.shape!=gradients.shape or coords.shape[1:]!=(len(atoms),3) or not np.isfinite(coords).all() or not np.isfinite(gradients).all():raise ValueError('Malformed saved FD arrays')
                close(float(abs(gradients[0]).max()),native['gradient_max'],1e-10)
                if np.max(abs(coords[0]*.529177210903-xyz))>1e-7:raise ValueError('Base gradient and cached geometry differ')
                expected_checks=[{k:v for k,v in r.items() if k not in ('minus_index','plus_index')} for r in fd['checks']]
                if expected_checks!=native['derivative_checks'] or len(fd['points'])!=len(coords):raise ValueError('FD/native test binding differs')
                for index,point in enumerate(fd['points']):
                    s=point['settings']
                    if (point['SCF_result']['converged'] is not True
                            or (s['charge'],s['mult'])!=(row['charge'],row['mult'])
                            or s['xc']!='r2scan' or s['disp']!='d4' or s['conv_tol']!=1e-12 or s['conv_tol_grad']!=1e-8 or s['grid_level']!=3
                            or s['basis']!=manifest['native_protocol']['basis']
                            or not math.isfinite(point['SCF_result']['energy_Ha'])
                            or np.max(abs(coords[index]-np.asarray(s['coordinates_Bohr'])))>1e-12):
                        raise ValueError('Saved FD calculation state/settings differ')
                for entry in fd['checks']:
                    atom,axis,h=entry['atom_index'],entry['axis'],entry['step_Bohr']
                    if h not in (.0005,.001):raise ValueError('Unexpected FD step')
                    minus,plus=entry['minus_index'],entry['plus_index']
                    for index,sign in ((minus,-1),(plus,1)):
                        target=coords[0].copy();target[atom,axis]+=sign*h
                        if np.max(abs(coords[index]-target))>3e-9:raise ValueError('Saved gradient has the wrong displacement')
                    numerical=(gradients[plus]-gradients[minus]).ravel()/(2*h)
                    error=float(np.linalg.norm(numerical-hessian[3*atom+axis])/np.linalg.norm(numerical))
                    close(error,entry['relative_error'],1e-10)
                    if not math.isfinite(error) or error>.01:raise ValueError('Saved FD column fails')
                    fd_checks.append(dict(atom=atom,axis=axis,step_Bohr=h,relative_error=error))
                if len(fd_checks)!=4:raise ValueError('Two columns at two steps required')
                symbols=[a[0] for a in atoms]
                other=next((i for i,s in enumerate(symbols) if s=='O'),next((i for i,s in enumerate(symbols) if s=='Fe'),len(atoms)-1))
                required={(atom,axis,h) for atom,axis in [(0,1),(other,2)] for h in [.0005,.001]}
                if {(r['atom'],r['axis'],r['step_Bohr']) for r in fd_checks}!=required:raise ValueError('FD column coverage differs')
                fd_states+=1
        totals={c:math.fsum([energy['e_smd_eV'],nuclear[c],spin]) for c in nuclear}
        for c,value in totals.items():close(value,row['checked_G_by_cutoff_eV'][c])
        values[name]=totals;inventories[name]=Counter(a[0] for a in atoms);charges[name]=row['charge']
        records.append(dict(id=name,frequency_error_cm=frequency_error,
            nuclear_errors_eV={c:abs(nuclear[c]-(partition['G_thermal_eV'] if row['partition_kind']=='monatomic' else partition['nuclear_by_cutoff_eV'][c])) for c in nuclear},
            finite_difference_checks=fd_checks,FD_recomputed=bool(fd_checks)))
    if len(values)!=10 or fd_states!=6:raise ValueError('Ten exact states, including six saved-gradient panels, required')
    factor=1.602176634e-19*6.02214076e23/1000;outputs={}
    for c in ('25','50','100'):
        state={name:value[c] for name,value in values.items()};reactions={}
        for reaction in manifest['reactions']:
            if reaction['kind'] not in ('one_electron_potential_vs_Fc','grafting_free_energy','reaction_free_energy'):raise ValueError('Unknown balanced observable')
            coeff=reaction['coefficients'];balance=Counter()
            for name,nu in coeff.items():
                if type(nu) is not int:raise ValueError('Integral stoichiometry required')
                for symbol,n in inventories[name].items():balance[symbol]+=nu*n
            if any(balance.values()) or sum(charges[name]*nu for name,nu in coeff.items()) or sum(coeff.values())!=0:
                raise ValueError('Atom, charge or solute-count balance fails')
            dg=math.fsum(nu*state[name] for name,nu in coeff.items())
            if reaction['kind']=='one_electron_potential_vs_Fc' and (reaction.get('electrons')!=1 or coeff.get('fc_ox')!=1 or coeff.get('fc_neu')!=-1):raise ValueError('One-electron Fc reference required')
            reactions[reaction['id']]=(-dg if reaction['kind']=='one_electron_potential_vs_Fc' else dg*factor)
            close(reactions[reaction['id']],reaction['expected_by_cutoff'][c],1e-6)
        outputs[c]=reactions
    selected=min(manifest['precursor_ids'],key=lambda name:values[name]['50'])
    if selected!=manifest['selected_precursor_id']:raise ValueError('Primary lowest-precursor-G selection differs')
    return dict(status='exact_viologen_cached_redox_and_quaternization_reconstructed',states=len(values),records=records,
        by_cutoff=outputs,selected_precursor_id=selected,Hessians_reconstructed=9,FD_states_reconstructed=fd_states,
        monatomic_partitions_reconstructed=1,quantum_calculations_run=False,xTB_library_loaded=False,
        global_conformer_convergence_proven=False,material_admission=False,full_historical_source_closure_revalidated=False,
        source_manifest_sha256=sha(root/'manifest.json'),numpy_version=np.__version__)


if __name__=='__main__':
    here=Path(__file__).resolve().parent;default=here/'data' if (here/'data/manifest.json').is_file() else here
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('replay',));p.add_argument('--package',type=Path,default=default);p.add_argument('--out',type=Path);a=p.parse_args()
    result=replay(a.package.resolve())
    if a.out:
        with a.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(status=result['status'],states=result['states'],by_cutoff=result['by_cutoff'])))
