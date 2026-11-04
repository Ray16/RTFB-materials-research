"""Implementation-impact controls for exact states of a portable reaction package.

Reuses the established electronic worker, reservation dispatcher, resource-pilot
review and independent SCF-checkpoint reviewer. No experimental endpoints are
invented and no new electronic energy is promoted to a total free energy.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT, HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.workflows.screening_energy_panel import output_budget


def inventory(package):
    from ase.io import read
    from pyscf import gto
    data=json.loads((package/'manifest.json').read_text());rows=[];seen=set()
    prefix=data.get('source_alias_prefix','exact_viologen')
    if not isinstance(prefix,str) or not prefix or any(c in prefix for c in (':','/','\\')):
        raise ValueError('A single explicit source alias namespace is required')
    for state in data['states']:
        key=state['id']
        if key in seen or not key or Path(key).name!=key:raise ValueError('Distinct safe package state IDs required')
        seen.add(key)
        paths={name:package/state['files'][name] for name in ('geometry','energy')}
        for name,path in paths.items():
            if not path.resolve().is_relative_to(package.resolve()) or sha256_file(path)!=data['files_sha256'][state['files'][name]]:
                raise ValueError('Bound package input changed')
        energy=json.loads(paths['energy'].read_text())
        if (energy['protocol']!=ACTIVE_SP or energy['charge']!=state['charge'] or energy['mult']!=state['mult']
                or type(state['charge']) is not int or type(state['mult']) is not int
                or energy['converged_smd'] is not True
                or energy['geom_sha1']!=hashlib.sha1(paths['geometry'].read_bytes()).hexdigest()[:16]):
            raise ValueError('Exact geometry, electronic protocol and charge/spin required')
        value=energy.get('e_smd_Ha',energy['e_smd_eV']/HARTREE_EV)
        if (any(type(v) not in (int,float) or not math.isfinite(v) for v in (value,energy['e_smd_eV']))
                or abs(value*HARTREE_EV-energy['e_smd_eV'])>1e-8):
            raise ValueError('Finite consistent energy units required')
        atoms=read(paths['geometry']);symbols=atoms.get_chemical_symbols()
        molecule=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),basis=ACTIVE_SP['basis'],
                       charge=state['charge'],spin=state['mult']-1,verbose=0)
        rows.append(dict(key=key,source_state=dict(id=key,aliases=[prefix+':'+key],
            smiles=state['smiles'],charge=state['charge'],mult=state['mult'],original_state_directory=str(paths['geometry'].parent)),
            geometry=str(paths['geometry']),energy=str(paths['energy']),reference_energy_Ha=value,
            atoms=len(atoms),nao=molecule.nao,has_iron='Fe' in symbols,
            output_budget_bytes=output_budget(molecule.nao,state['mult'],'Fe' in symbols)))
    if not rows:raise ValueError('No package states')
    return data,rows


def prepare(archive,package,output):
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Fresh durable namespace outside both scientific checkouts required')
    data,rows=inventory(package)
    output.mkdir(parents=True);inputs=output/'inputs';inputs.mkdir()
    shutil.copyfile(package/'manifest.json',inputs/'package_manifest.json')
    module=Path(importlib.util.find_spec('gpu4pyscf').origin).parent;prepared=[]
    for row in rows:
        root=output/'states'/row['key'];case=root/'case';sources=root/'sources'
        case.mkdir(parents=True);sources.mkdir()
        for name in ('dft','locks','repairs'):(case/name).mkdir()
        copies=[(ROOT/'scripts/validation/run_screening_energy_control.py','worker.py'),
            (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),(ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
            (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
            (Path(__file__),'preparer.py'),(Path(row['geometry']),'input.xyz'),(Path(row['energy']),'energy_reference.json')]
        for source,name in copies:shutil.copyfile(source,sources/name)
        files=[*sources.glob('*'),*inputs.glob('*'),*(archive/'src/redox/core').glob('*.py'),
            archive/'src/redox/workflows/resources.py',archive/'src/redox/qm/dft.py',
            archive/'scripts/pipeline/run_pi_stacking.py',archive/'scripts/fleet/exclusive_python_worker.py',
            archive/'scripts/fleet/tensor_runtime_worker.py',archive/'environments/cutensor22/runtime.json',
            module/'dft/gen_grid.py',module/'dft/numint.py',module/'lib/libgdft.so']
        plan=dict(role='screening_energy_control',energy_control=True,source_id=row['key'],state='rest',
            source_state=row['source_state'],electronic_protocol=ACTIVE_SP,reference_energy_Ha=row['reference_energy_Ha'],
            atoms=row['atoms'],nao=row['nao'],has_iron=row['has_iron'],source_sha256={str(p):sha256_file(p) for p in files},
            maximum_host_RSS_GiB=40,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
            maximum_wall_seconds=7200,maximum_durable_output_bytes=row['output_budget_bytes'],
            minimum_initial_storage_GiB=32+row['output_budget_bytes']/2**30,minimum_live_storage_GiB=32,
            comparison_absolute_tolerance_Ha=1e-7,native_G_qualified=False,chemical_free_energy_assembled=False,
            purpose='Exact candidate/reaction electronic implementation sensitivity; no new qualified potential or grafting G')
        atomic_write_json(case/'plan.json',plan)
        spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),source_id=row['key'],state='rest',
            maximum_host_RSS_GiB=40,minimum_free_GiB=32,scientific_sources_read_only=True,plan_sha256=sha256_file(case/'plan.json'))
        atomic_write_json(root/'namespace.json',spec)
        prepared.append(dict(row,namespace=str(root/'namespace.json'),namespace_sha256=sha256_file(root/'namespace.json')))
    result=dict(status='exact_component_screening_panel_prepared_no_jobs_submitted',target_states=len(prepared),states=prepared,
        workload=data.get('scope_id','exact_viologen_package'),
        source_package=str(package),source_package_manifest_sha256=sha256_file(package/'manifest.json'),
        maximum_total_durable_output_bytes=sum(r['output_budget_bytes'] for r in rows),
        source_sha256={str(inputs/'package_manifest.json'):sha256_file(inputs/'package_manifest.json')},
        preparer_sha256=sha256_file(__file__),native_G_qualified=False,new_redox_potentials_qualified=False,quantum_jobs_submitted=0)
    atomic_write_json(output/'manifest.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','package','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=prepare(a.archive.resolve(),a.package.resolve(),a.out.resolve())
    print(json.dumps({k:r[k] for k in ['status','target_states','maximum_total_durable_output_bytes']}))
