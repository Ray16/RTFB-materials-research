#!/usr/bin/env python
"""Extend the checked first-site raw package with the completed PMDI pair."""
import argparse
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.protocol import load_record,record_path
from redox.core.io import atomic_write_json
from redox.qm.native_thermal import cache_path,load_native_thermal
from redox.workflows import polymer_campaign as api


def prepare(first_site,source,out):
    original=json.loads((first_site/'manifest.json').read_text())
    for relative,expected in original['files_sha256'].items():
        path=(first_site/relative).resolve()
        if not path.is_relative_to(first_site) or api.digest(path)!=expected:
            raise ValueError('Original first-site package bytes differ or escape containment')
    report=json.loads(source.read_text());api.verify_sources(report['source_sha256'])
    if report['status']!='exact_candidate_neighbor_grafting_cycle_independently_reconstructed':
        raise ValueError('An independently reconstructed neighbor cycle is required')
    state=Path(report['states']['pair']['state_dir']);raw=json.loads((state/'result.json').read_text())
    charge,mult=raw['charge'],raw['mult'];energy=load_record(state,charge,mult);native=load_native_thermal(state,charge,mult)
    if energy is None or native is None:raise ValueError('Paired state has no accepted exact electronic/native parents')
    if abs(energy['e_smd_eV']+native['g_thermal_total_eV']-report['states']['pair']['G_eV'])>1e-8:
        raise ValueError('Paired raw G differs from the completed cycle')
    native_path,_=cache_path(state,charge,mult)
    paths=dict(geometry=state/'opt.xyz',result=state/'result.json',energy=record_path(state,charge,mult),
        native=native_path,partition=Path(native['thermal_postprocessing_path']),
        hessian=native_path.with_suffix('')/'hessian_cartesian.npy')
    out.mkdir(parents=True,exist_ok=False);sources={str(first_site/'manifest.json'):api.digest(first_site/'manifest.json')}
    def copy(path,target):
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        if api.digest(path)!=api.digest(target):raise ValueError('Copied archive bytes differ')
        sources[str(path)]=api.digest(path)
    for p in first_site.rglob('*'):
        if p.is_symlink():raise ValueError('First-site source contains a symlink')
        if p.is_file():copy(p,out/p.relative_to(first_site))
    shutil.copyfile(out/'manifest.json',out/'first_site_manifest.json')
    shutil.move(out/'reproduce.py',out/'cached_reader.py')
    shutil.move(out/'README.md',out/'first_site_README.md')
    copy(ROOT/'examples/neighbor_cycle_replay/reproduce.py',out/'reproduce.py')
    copy(Path(__file__).resolve(),out/'prepare_neighbor_source.py');copy(source,out/'neighbor_source.json')
    files={}
    for role,p in paths.items():
        target=out/'polymer_pair'/(role+p.suffix);copy(p,target);files[role]=str(target.relative_to(out))
    original['states'].append(dict(id='polymer_pair',smiles=report['states']['pair']['smiles'],charge=charge,mult=mult,
        partition_kind='projected_DFT',files=files,checked_G_eV=report['states']['pair']['G_eV'],original_state_directory=str(state)))
    for key,row in report['states'].items():
        name='polymer_pair' if key=='pair' else key
        expected=next(r for r in original['states'] if r['id']==name)
        if abs(expected['checked_G_eV']-row['G_eV'])>1e-8:raise ValueError('Inherited single/empty state differs from neighbor cycle')
    original['reactions'].append(dict(name='neighbor_nonadditivity',
        coefficients={'polymer_pair':1,'empty':1,'site2':-1,'site3':-1},delta_n_solutes=0,source_kJ_mol=report['deltaDeltaG_kJ_mol']))
    for mode in ('free_ions','contact_pair'):
        reagents={'protonated':1,'chloride':1} if mode=='free_ions' else {'pair':1}
        for site,other in ((2,3),(3,2)):
            coefficients={'polymer_pair':1,f'site{other}':-1,'precursor':-1,'base':-1,**reagents}
            original['reactions'].append(dict(name=f'graft{site}_given{other}_{mode}',coefficients=coefficients,
                delta_n_solutes=sum(coefficients.values()),source_kJ_mol=report['by_cutoff']['50']['models'][mode][f'graft_{site}_given_{other}_kJ_mol']))
        coefficients={'polymer_pair':1,'empty':-1,'precursor':-2,'base':-2,**{k:2*v for k,v in reagents.items()}}
        original['reactions'].append(dict(name='direct_double_'+mode,coefficients=coefficients,
            delta_n_solutes=sum(coefficients.values()),source_kJ_mol=report['by_cutoff']['50']['models'][mode]['direct_double_graft_kJ_mol']))
    original.update(neighbor_source='neighbor_source.json',source_sha256={**original['source_sha256'],**sources},
        conformer_convergence_established=False,synthesis_admission=False,isolated_motif_interaction_computed=False,
        electronic_coupling_computed=False,scope='Finite exact-PMDI native neighbor cycle and conditional DBU grafting; cached arithmetic, not independent nuclear integration or predictive validation')
    (out/'README.md').write_text('''# Exact-PMDI neighbor-cycle raw addendum

Run from any directory with standard Python and no site packages:

`python3 -I -S /path/to/package/reproduce.py --data /path/to/package`

Fourteen cached states include the empty, two singly grafted and doubly grafted
chains, local/reagent references and a molecular redox/Fc example. Each non-atomic
state includes its original geometry, electronic record, native Hessian, modes
and nuclear partition. Fourteen balanced reactions replay at the primary cutoff;
the complete conditional neighbor cycle additionally replays at 25/50/100 cm^-1.
Both attachment orders close. The replay verifies contained paths and hashes and
never opens the original cluster paths. It does not run QM, reintegrate nuclear
terms, establish converged ensembles, isolate motif binding or calculate coupling.
This extends the immutable first-site package without changing it. The archive
is locally prepared; public data licensing, author metadata and DOI remain open.
''')
    original['files_sha256']={str(p.relative_to(out)):api.digest(p) for p in out.rglob('*') if p.is_file() and p!=out/'manifest.json'}
    api.verify_sources(sources);atomic_write_json(out/'manifest.json',original)
    return dict(status='raw_neighbor_addendum_prepared_replay_pending',states=len(original['states']),
        reactions=len(original['reactions']),prepared=time.time(),publicly_deposited=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('first-site','source','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(a.first_site.resolve(),a.source.resolve(),a.out.resolve())))


if __name__=='__main__':main()
