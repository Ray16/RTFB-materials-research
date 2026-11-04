"""Source-exact benchmark implementation controls and guarded case preparation."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT,HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def output_budget(nao,mult,has_iron):
    """Bound both variants' checkpoint orbitals, including three SCF fallbacks."""
    if type(nao) is not int or nao<1 or type(mult) is not int or mult<1 or type(has_iron) is not bool:
        raise ValueError('Explicit positive AO/state dimensions and element flag required')
    guesses=6 if has_iron and mult==2 else 3 if mult>1 else 1
    spins=2 if mult>1 else 1
    return 16*2**20+2*guesses*3*spins*nao*nao*8


def benchmark_inventory(benchmark,data,package):
    """Use every state of every frozen endpoint and its exact Fc reference."""
    from ase.io import read
    from pyscf import gto
    from redox.core.protocol import ACTIVE_SP
    aliases={}
    for row in data['states']:
        for alias in row['aliases']:
            if alias in aliases:raise ValueError('Ambiguous source state alias')
            aliases[alias]=row
    required={'benchmark:ferrocene:ox','benchmark:ferrocene:neu'}
    for row in benchmark['outliers']:
        names=row['event'].split('->')
        if len(names)!=2:raise ValueError('One explicit redox-state pair per endpoint required')
        required.update('benchmark:'+row['id']+':'+name for name in names)
    if benchmark['admitted_events']!=len(benchmark['outliers']) or not required<=set(aliases):
        raise ValueError('Frozen endpoint coverage is incomplete')
    selected={aliases[a]['id']:aliases[a] for a in required}
    inventory=[]
    for key,row in sorted(selected.items()):
        paths={name:package/row['files'][name] for name in ('geometry','energy')}
        for name,path in paths.items():
            if sha256_file(path)!=data['files_sha256'][row['files'][name]]:
                raise ValueError('Frozen state object changed: '+str(path))
        energy=json.loads(paths['energy'].read_text())
        if (energy['protocol']!=ACTIVE_SP or energy['charge']!=row['charge'] or energy['mult']!=row['mult']
                or energy['converged_smd'] is not True or energy['geom_sha1']!=hashlib.sha1(paths['geometry'].read_bytes()).hexdigest()[:16]
                or type(row['charge']) is not int or type(row['mult']) is not int):
            raise ValueError('Exact geometry/state/active-protocol parent required')
        reference=energy.get('e_smd_Ha',energy['e_smd_eV']/HARTREE_EV)
        if not math.isfinite(reference) or abs(reference*HARTREE_EV-energy['e_smd_eV'])>1e-8:
            raise ValueError('Finite consistent stored electronic energy required')
        atoms=read(paths['geometry']);symbols=atoms.get_chemical_symbols()
        molecule=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),basis=ACTIVE_SP['basis'],
                       charge=row['charge'],spin=row['mult']-1,verbose=0)
        inventory.append(dict(key=key,source_state={k:row[k] for k in ('id','aliases','smiles','charge','mult','original_state_directory')},
            geometry=str(paths['geometry']),energy=str(paths['energy']),reference_energy_Ha=reference,
            atoms=len(atoms),nao=molecule.nao,has_iron='Fe' in symbols,
            output_budget_bytes=output_budget(molecule.nao,row['mult'],'Fe' in symbols)))
    return inventory


def prepare(archive,benchmark_path,data_path,package,output):
    from redox.core.protocol import ACTIVE_SP
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Use a fresh durable panel outside both scientific checkouts')
    benchmark=json.loads(benchmark_path.read_text());data=json.loads(data_path.read_text())
    rows=benchmark_inventory(benchmark,data,package)
    output.mkdir(parents=True);inputs=output/'inputs';inputs.mkdir()
    shutil.copyfile(benchmark_path,inputs/'benchmark.json');shutil.copyfile(data_path,inputs/'source_data.json')
    package_dir=Path(importlib.util.find_spec('gpu4pyscf').origin).parent
    prepared=[]
    for row in rows:
        root=output/'states'/row['key'];case=root/'case';sources=root/'sources';case.mkdir(parents=True);sources.mkdir()
        for name in ('dft','locks','repairs'):(case/name).mkdir()
        copies=[(ROOT/'scripts/validation/run_screening_energy_control.py','worker.py'),
            (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),
            (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
            (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
            (Path(__file__),'preparer.py'),(Path(row['geometry']),'input.xyz'),(Path(row['energy']),'energy_reference.json')]
        for source,name in copies:shutil.copyfile(source,sources/name)
        files=[*sources.glob('*'),*inputs.glob('*'),*[p for p in (archive/'src/redox/core').glob('*.py')],
            archive/'src/redox/workflows/resources.py',archive/'src/redox/qm/dft.py',
            archive/'scripts/pipeline/run_pi_stacking.py',archive/'scripts/fleet/exclusive_python_worker.py',
            archive/'scripts/fleet/tensor_runtime_worker.py',archive/'environments/cutensor22/runtime.json',
            package_dir/'dft/gen_grid.py',package_dir/'dft/numint.py',package_dir/'lib/libgdft.so']
        plan=dict(role='screening_energy_control',energy_control=True,source_id=row['key'],state='rest',
            source_state=row['source_state'],electronic_protocol=ACTIVE_SP,reference_energy_Ha=row['reference_energy_Ha'],
            atoms=row['atoms'],nao=row['nao'],has_iron=row['has_iron'],source_sha256={str(p):sha256_file(p) for p in files},
            maximum_host_RSS_GiB=40,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
            maximum_wall_seconds=7200,maximum_durable_output_bytes=row['output_budget_bytes'],
            minimum_initial_storage_GiB=32+row['output_budget_bytes']/2**30,minimum_live_storage_GiB=32,
            comparison_absolute_tolerance_Ha=1e-7,native_G_qualified=False,chemical_free_energy_assembled=False,
            purpose='Matched source-frozen benchmark electronic implementation sensitivity; no new qualified E0 or grafting G')
        atomic_write_json(case/'plan.json',plan)
        spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),source_id=row['key'],state='rest',
            maximum_host_RSS_GiB=40,minimum_free_GiB=32,scientific_sources_read_only=True,plan_sha256=sha256_file(case/'plan.json'))
        atomic_write_json(root/'namespace.json',spec)
        prepared.append(dict(row,namespace=str(root/'namespace.json'),namespace_sha256=sha256_file(root/'namespace.json')))
    manifest=dict(status='matched_benchmark_screening_energy_panel_prepared_no_jobs_submitted',
        target_endpoints=benchmark['admitted_events'],target_states=len(rows),states=prepared,
        maximum_total_durable_output_bytes=sum(r['output_budget_bytes'] for r in rows),
        source_sha256={str(p):sha256_file(p) for p in [benchmark_path,data_path,Path(__file__),*inputs.glob('*')]},
        native_G_qualified=False,new_redox_potentials_qualified=False,quantum_jobs_submitted=0)
    atomic_write_json(output/'manifest.json',manifest)
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare',))
    for name in ('archive','benchmark','data','package','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare(args.archive.resolve(),args.benchmark.resolve(),args.data.resolve(),args.package.resolve(),args.out.resolve())
    print(json.dumps({k:result[k] for k in ('status','target_endpoints','target_states','maximum_total_durable_output_bytes')}))
