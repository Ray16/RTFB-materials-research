"""Exact local-graft chloride-contact hypotheses, separate from redox states.

The oxidized viologen core remains formally +2; adding Cl- gives a +1 singlet
microstate. It is not the isolated +1 doublet obtained by one-electron reduction.
These geometries test ion-convention sensitivity, not synthesis eligibility.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

import numpy as np
from ase import Atoms
from ase.io import read,write
from rdkit import Chem
from pyscf import gto

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def launch(namespace,gpu,output,pilot_review=None):
    from redox.workflows.molecular_native_namespace import verify,launch as launch_process
    from redox.workflows.screening_control import probe
    from redox.workflows.response_diagnostic import observe_budget
    import time
    spec=verify(namespace);plan=json.loads((Path(spec['case'])/'plan.json').read_text())
    if plan['role']!='primary_counterion_geometry_control' or plan['maximum_host_RSS_GiB']!=32:
        raise ValueError('Exact primary chloride-contact geometry resource pilot required')
    if plan['source_id']!='N0_above':
        if pilot_review is None:raise ValueError('First independently checked chloride-geometry pilot required')
        reference=json.loads(pilot_review.read_text())
        if (reference['status']!='primary_counterion_geometry_independently_checked' or reference['key']!='N0_above'
                or reference['nao']<plan['nao'] or reference['native_G_qualified']
                or any(sha256_file(p)!=h for p,h in reference['source_sha256'].items())):
            raise ValueError('Exact physically checked chloride pilot required before wider geometry release')
    probe(namespace);child,result=launch_process(namespace,gpu,output)
    observer_source=Path(observe_budget.__code__.co_filename);shutil.copyfile(observer_source,output/'budget_observer_source.py')
    result.update(maximum_wall_seconds=plan['maximum_wall_seconds'],maximum_durable_output_bytes=plan['maximum_durable_output_bytes'],
        budget_observer_owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        budget_observer_source_sha256=sha256_file(output/'budget_observer_source.py'))
    atomic_write_json(output/'launch.json',result)
    result=observe_budget(child,spec,result,output);return result


def seeds(atoms):
    symbols=atoms.get_chemical_symbols();nitrogen=[i for i,s in enumerate(symbols) if s=='N']
    if len(nitrogen)!=2:raise ValueError('Exact two-nitrogen viologen parent required')
    a,b=atoms.positions[nitrogen];axis=b-a;axis/=np.linalg.norm(axis)
    carbons=sorted([i for i,s in enumerate(symbols) if s=='C'],key=lambda i:np.linalg.norm(atoms.positions[i]-a))
    normal=next(np.cross(axis,atoms.positions[i]-a) for i in carbons if np.linalg.norm(np.cross(axis,atoms.positions[i]-a))>1e-5)
    normal/=np.linalg.norm(normal)
    for site,origin in enumerate((a,b)):
        for label,direction in [('above',normal),('below',-normal),('outer',axis*(1 if site else -1))]:
            accepted=None
            for distance in (3.8,4.4,5.0,5.6):
                point=origin+distance*direction
                if min(np.linalg.norm(atoms.positions-point,axis=1))>=2.0:accepted=(point,distance);break
            if accepted is None:raise ValueError('Declared chloride seed has unresolved atomic overlap')
            yield f'N{site}_{label}',accepted[0],dict(nitrogen_atom=nitrogen[site],direction=label,distance_A=accepted[1])
    yield 'separated_start',(a+b)/2+9*normal,dict(direction='core_normal',distance_from_N_midpoint_A=9.)


def prepare(archive,package,output):
    from redox.validation.protonation import hydrogen_topology_status
    from redox.core.runtime_paths import gpu4pyscf_source_root
    manifest=json.loads((package/'manifest.json').read_text());parent=next(s for s in manifest['states'] if s['id']=='graft_ox2')
    geometry=package/parent['files']['geometry']
    if sha256_file(geometry)!=manifest['files_sha256'][parent['files']['geometry']] or (parent['charge'],parent['mult'])!=(2,1):
        raise ValueError('Source-bound exact oxidized graft geometry required')
    smiles=parent['smiles']+'.[Cl-]';mol=Chem.MolFromSmiles(smiles)
    if Chem.GetFormalCharge(mol)!=1:raise ValueError('The contact pair must be a net +1 species')
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive,package)):
        raise FileExistsError('Fresh durable counterion-control namespace outside source trees required')
    output.mkdir(parents=True);atoms=read(geometry);rows=[];installed=gpu4pyscf_source_root()
    for key,position,definition in seeds(atoms):
        root=output/'states'/key;case=root/'case';sources=root/'sources';case.mkdir(parents=True);sources.mkdir()
        paired=atoms.copy();paired+=Atoms('Cl',positions=[position]);write(case/'input.xyz',paired,format='xyz')
        if hydrogen_topology_status(smiles,case/'input.xyz')!='ok':raise ValueError('Actual seed differs from the exact ionic graph')
        molecule=gto.M(atom=list(zip(paired.get_chemical_symbols(),paired.positions.tolist())),basis='def2-svp',charge=1,spin=0,verbose=0)
        if molecule.nao_nr()>512 or molecule.nelectron%2:raise ValueError('Explicit small closed-shell geometry pilot dimension required')
        for original,name in [(ROOT/'scripts/validation/run_counterion_geometry.py','worker.py'),
                (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
                (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),
                (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
                (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
                (package/'manifest.json','parent_manifest.json'),
                (Path(__file__),'preparer.py'),(geometry,'parent_geometry.xyz')]:shutil.copyfile(original,sources/name)
        files=[*sources.glob('*'),case/'input.xyz',
            *(archive/'src/redox/qm').glob('*.py'),*(archive/'src/redox/core').glob('*.py'),
            *(archive/'config').glob('*.py'),*(archive/'config').glob('*.json'),
            archive/'src/redox/workflows/resources.py',archive/'scripts/pipeline/run_pi_stacking.py',
            archive/'scripts/fleet/exclusive_python_worker.py',archive/'scripts/fleet/tensor_runtime_worker.py',
            archive/'environments/cutensor22/runtime.json',installed/'dft/gen_grid.py',installed/'dft/numint.py',installed/'lib/libgdft.so']
        plan=dict(role='primary_counterion_geometry_control',source_id=key,state='oxidized_core_chloride_contact',
            parent_smiles=parent['smiles'],smiles=smiles,parent_formal_core_charge=2,counterion_charge=-1,charge=1,mult=1,
            molecular_state_distinct_from_reduced_radical=True,seed=definition,atoms=len(paired),nao=molecule.nao_nr(),
            method=dict(xc='r2scan',basis='def2-svp',disp='d4',solvent='acetonitrile',solvent_model='primary SMD with legacy CDS',
                SCF_conv_tol=1e-12,SCF_conv_tol_grad=1e-8,SCF_max_cycle=200,grid_level=3,gradient_grid_response=True),
            maximum_host_RSS_GiB=32,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,maximum_wall_seconds=14400,
            maximum_durable_output_bytes=2*2**30,minimum_initial_storage_GiB=64,minimum_live_storage_GiB=32,
            source_sha256={str(p):sha256_file(p) for p in files},native_G_qualified=False,geometry_qualified=False,
            property_values_transferred=False,material_admission=False,experimental_speciation_weights=False,
            purpose='Incoming-chloride association sensitivity in the exact local oxidized graft; iodide remains an omitted separated spectator')
        atomic_write_json(case/'plan.json',plan)
        spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),source_id=key,
            state=plan['state'],scientific_sources_read_only=True,maximum_host_RSS_GiB=32,minimum_free_GiB=32,plan_sha256=sha256_file(case/'plan.json'))
        atomic_write_json(root/'namespace.json',spec)
        rows.append(dict(key=key,namespace=str(root/'namespace.json'),namespace_sha256=sha256_file(root/'namespace.json')))
    result=dict(status='seven_explicit_chloride_microstate_geometry_inputs_prepared',states=rows,pilot_key='N0_above',
        same_covalent_graft_core=True,formal_core_charge=2,total_charge=1,multiplicity=1,
        geometry_only=True,native_G_qualified=False,experimental_populations_inferred=False,quantum_jobs_submitted=0)
    atomic_write_json(output/'manifest.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',nargs='?',choices=('prepare','launch'),default='prepare')
    for name in ('archive','package','out','spec','pilot-review'):p.add_argument('--'+name,type=Path)
    p.add_argument('--gpu',type=int)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    if a.action=='prepare':
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
        if any(v is None for v in (a.archive,a.package,a.out)):p.error('Preparation requires archive, exact source package and fresh output')
        result=prepare(a.archive.resolve(),a.package.resolve(),a.out.resolve());print(json.dumps({k:result[k] for k in ('status','pilot_key','total_charge','multiplicity')}))
    else:
        if a.spec is None or a.gpu is None or a.out is None:p.error('Launch requires exact namespace, device and fresh placement')
        result=launch(a.spec.resolve(),a.gpu,a.out.resolve(),None if a.pilot_review is None else a.pilot_review.resolve())
        print(json.dumps(dict(status=result['status'],returncode=result['returncode'],budget_stop=result['budget_stop'])))
