#!/usr/bin/env python
"""Export raw first-site polymer states for relocated standard-library replay.

Extends the earlier reviewer package with the newer first-site result. Original
state bytes are copied without editing certificates or invoking quantum jobs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT,DFT,read_manifest
from redox.core.protocol import load_record,record_path
from redox.qm.native_thermal import cache_path,load_native_thermal


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    from redox.workflows import polymer_campaign as api
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    source=args.source.resolve();report=json.loads(source.read_text())
    api.verify_sources(report['source_sha256'])
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    reader=ROOT/'examples/cached_thermochemistry/reproduce.py'
    shutil.copyfile(reader,out/'reproduce.py')
    shutil.copyfile(Path(__file__).resolve(),out/'prepare_source.py')
    shutil.copyfile(source,out/'first_site_source.json')
    sources={str(p):digest(p) for p in (source,reader,Path(__file__).resolve())}
    selected={key:dict(row) for key,row in report['reagent_states'].items()}
    selected['empty']=report['empty_reference'];selected['capped']=report['local_capped_reference']
    for row in report['singles']:
        selected['site'+str(row['site'])]=dict(row,smiles=row['task']['smiles'])
    manifest={(r['id'],r['state']):r for r in read_manifest()}
    for key,gid,state in [('fc_neutral','ferrocene','neu'),('fc_oxidized','ferrocene','ox'),
            ('capped_reduced','design_pmdi_ethyl','red1')]:
        selected[key]=dict(state_dir=str(DFT/gid/state),smiles=manifest[gid,state]['smiles'])
    exported=[]
    for name,row in selected.items():
        state=Path(row['state_dir']);raw=json.loads((state/'result.json').read_text())
        charge,mult=raw['charge'],raw['mult'];active=load_record(state,charge,mult)
        if active is None:raise ValueError('Unqualified electronic parent: '+name)
        paths=dict(geometry=state/'opt.xyz',result=state/'result.json',energy=record_path(state,charge,mult))
        if (state/'monatomic_thermal.json').exists():
            paths['partition']=state/'monatomic_thermal.json'
            thermal=json.loads(paths['partition'].read_text())['G_thermal_eV'];kind='monatomic'
        else:
            native=load_native_thermal(state,charge,mult)
            if native is None:raise ValueError('Unqualified nuclear partition: '+name)
            path,_=cache_path(state,charge,mult)
            paths.update(native=path,partition=Path(native['thermal_postprocessing_path']),
                hessian=path.with_suffix('')/'hessian_cartesian.npy')
            thermal=native['g_thermal_total_eV'];kind='projected_DFT'
        total=active['e_smd_eV']+thermal
        if 'G_eV' in row and abs(total-row['G_eV'])>1e-8:
            raise ValueError('Raw state differs from checked first-site result: '+name)
        folder=out/name;folder.mkdir();files={}
        for role,path in paths.items():
            target=folder/(role+path.suffix);shutil.copyfile(path,target)
            if digest(path)!=digest(target):raise ValueError('Raw-state copy differs')
            sources[str(path)]=digest(path);files[role]=str(target.relative_to(out))
        exported.append(dict(id=name,smiles=row['smiles'],charge=charge,mult=mult,
            partition_kind=kind,files=files,checked_G_eV=total,original_state_directory=str(state)))
    reactions=[]
    for row in report['singles']:
        site='site'+str(row['site']);values=row['by_cutoff']['50']
        for mode in ('free_ions','contact_pair'):
            coefficients={site:1,'empty':-1,'precursor':-1,'base':-1}
            coefficients.update({'protonated':1,'chloride':1} if mode=='free_ions' else {'pair':1})
            reactions.append(dict(name=site+'_'+mode,coefficients=coefficients,
                delta_n_solutes=sum(coefficients.values()),source_kJ_mol=values['first_graft_kJ_mol'][mode]))
        reactions.append(dict(name=site+'_polymer_minus_local',coefficients={site:1,'local':1,'empty':-1,'capped':-1},
            delta_n_solutes=0,source_kJ_mol=values['polymer_minus_capped_kJ_mol']))
    reactions.append(dict(name='site2_minus_site3',coefficients={'site2':1,'site3':-1},delta_n_solutes=0,
        source_kJ_mol=report['placement']['site2_minus_site3_deltaG_kJ_mol']))
    data=dict(schema=1,temperature_K=report['temperature_K'],states=exported,reactions=reactions,
        potential=dict(oxidized='capped',reduced='capped_reduced',reference_oxidized='fc_oxidized',reference_reduced='fc_neutral'),
        source_energy_conversion_kJ_mol_per_eV=report['conversion_kJ_mol_per_eV'],
        SI_energy_conversion_kJ_mol_per_eV=report['conversion_kJ_mol_per_eV'],
        source_sha256=sources,scope='Primary-cutoff raw first-site arithmetic; no new QM, ensemble convergence or external validation',
        files_sha256={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()})
    api.verify_sources(sources)
    (out/'manifest.json').write_text(json.dumps(data,indent=2)+'\n')
    (out/'README.md').write_text('''# First-site PMDI raw-state addendum

Run with standard Python, from any directory:

`python3 -I -S /path/to/addendum/reproduce.py --data /path/to/addendum`

The replay verifies original cached state/partition/geometry identities, one
PMDI potential and seven atom/charge-balanced first-site or relative reactions.
All calculation objects retain original bytes. Absolute cluster paths are
descriptive provenance; the replay opens only relative packaged files.

The two polymer single-site states postdate the v5 reviewer archive. This
addendum complements that earlier 243-state package. It does not supply the
missing paired-state free energy, new conformer ensembles, reaction barriers,
material eligibility or a public archival DOI. Original proton-transfer and
speciation qualifications remain. Code retains its MIT license; public data
license and author metadata remain unconfirmed.
''')
    print(json.dumps(dict(status='first_site_raw_addendum_prepared_replay_pending',states=len(exported),reactions=len(reactions),
        package=str(out),prepared=time.time(),publicly_deposited=False)))


if __name__=='__main__':
    main()
