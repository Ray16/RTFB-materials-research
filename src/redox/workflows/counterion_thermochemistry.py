"""Prospective primary-model thermochemistry of checked chloride microstates.

The first completed contact geometry is a resource/physical pilot. This stage
reuses the molecular native equations and the uniform diffuse electronic basis.
It produces components, not ion-pair populations or a grafting-yield prediction.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil

from ase.io import read
from pyscf import gto

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.workflows.molecular_native_namespace import verify, load


def checked_component(path):
    pilot = json.loads(path.read_text())
    if (pilot.get('status') != 'primary_counterion_thermochemistry_independently_checked'
            or pilot.get('native_G_qualified') is not True
            or (pilot.get('charge'), pilot.get('mult'), pilot.get('formal_core_charge')) != (1, 1, 2)
            or pilot.get('dimensions') != {'def2-svp': 412, 'def2-tzvpd': 1003}
            or pilot.get('association_free_energy_computed') is not False
            or pilot.get('material_admission') is not False or pilot.get('conformer_ensemble_qualified') is not False
            or pilot.get('actual_native_SCF_frames_checked') != 9
            or not pilot.get('source_sha256')
            or any(sha256_file(p) != h for p, h in pilot['source_sha256'].items())):
        raise ValueError('Complete independent chloride thermochemistry pilot required')
    for key, low, high in [('peak_GPU_MiB', 1, 31744), ('peak_host_RSS_bytes', 1, 40*2**30),
                           ('minimum_frequency_cm', -5, float('inf')), ('maximum_column_error', 0, .01)]:
        value = pilot.get(key)
        if type(value) not in (int, float) or not math.isfinite(value) or not low <= value <= high:
            raise ValueError('Full pilot physical/resource quantities required')
    spec = verify(Path(pilot['namespace']))
    case = Path(spec['case']); plan = json.loads((case/'plan.json').read_text())
    raw = json.loads((case/'result.json').read_text())
    required = {str(p) for p in (Path(pilot['namespace']), case/'plan.json', case/'result.json',
                case/'electronic.json', Path(raw['native_record_path']))}
    if (plan['role'] != 'primary_counterion_thermochemistry'
            or plan['source_id'] != pilot['key'] or not required <= pilot['source_sha256'].keys()
            or any(pilot[key] != plan[key] for key in ('dimensions','electronic_protocol','native_protocol','thermal_protocol'))
            or pilot.get('temperature_K') != plan['thermal_protocol']['temperature_K']
            or pilot.get('standard_state') != plan['thermal_protocol']['pressure']
            or set(pilot.get('G_by_cutoff_eV', {})) != {'25','50','100'}
            or any(not math.isfinite(v) or abs(v-raw['G_by_cutoff_eV'][c]) > 1e-8 for c, v in pilot['G_by_cutoff_eV'].items())):
        raise ValueError('The pilot component and declared method do not match their actual raw values')
    return pilot, spec


def checked_pilot(path):
    report, spec = checked_component(path)
    if report['key'] != 'N0_above':
        raise ValueError('The first declared contact geometry is the thermochemistry resource pilot')
    return report, spec


def followup_host_bound(pilot):
    peak=pilot['peak_host_RSS_bytes']
    if type(peak) not in (int,float) or not math.isfinite(peak) or not 0<peak<=40*2**30:
        raise ValueError('Measured finite full-pilot host peak required')
    return 16 if 1.25*peak<=16*2**30 else 40


def prepare(archive, geometry_review, output, pilot_review=None):
    from redox.core.runtime_paths import gpu4pyscf_source_root
    from redox.validation.protonation import hydrogen_topology_status
    review = json.loads(geometry_review.read_text())
    if (review['status'] != 'primary_counterion_geometry_independently_checked'
            or (review['total_charge'], review['multiplicity'], review['formal_core_charge']) != (1, 1, 2)
            or review['native_G_qualified'] or review['association_free_energy_computed']
            or any(sha256_file(p) != h for p, h in review['source_sha256'].items())):
        raise ValueError('Exact independently checked primary chloride geometry required')
    original = verify(Path(review['namespace']))
    prior = Path(original['case'])
    geometry = prior / 'optimizer_return.xyz'
    if sha256_file(geometry) != review['geometry_sha256'] or hydrogen_topology_status(review['smiles'], geometry) != 'ok':
        raise ValueError('Actual chloride/graft geometry differs from its checked graph')
    if review['key'] != 'N0_above' and pilot_review is None:
        raise ValueError('The first full chloride thermochemistry pilot must pass before wider preparation')
    if pilot_review is not None:
        pilot, pilot_spec = checked_pilot(pilot_review)
    atoms = read(geometry)
    dimensions = {basis: gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
        basis=basis, charge=1, spin=0, verbose=0).nao_nr() for basis in ('def2-svp', ACTIVE_SP['basis'])}
    if len(atoms) != 42 or dimensions != {'def2-svp': 412, 'def2-tzvpd': 1003}:
        raise ValueError('Only the exact prepared closed-shell chloride-contact dimension class is authorized')
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT, archive)):
        raise FileExistsError('Fresh durable namespace outside both source trees required')
    output.mkdir(parents=True)
    case, sources = output / 'case', output / 'sources'
    case.mkdir(); sources.mkdir(); (case / 'state').mkdir(); (case / 'library').mkdir()
    for name in ('dft', 'locks', 'repairs'):
        (case / name).mkdir()
    for source, name in [(ROOT/'scripts/validation/run_counterion_thermochemistry.py', 'worker.py'),
            (ROOT/'src/redox/core/hashing.py', 'hashing.py'),
            (ROOT/'src/redox/qm/scf_checkpoints.py', 'checkpoints.py'),
            (ROOT/'src/redox/qm/screening_barrier.py', 'helper.py'),
            (ROOT/'tests/fixtures/screen_index_v1_8_1.cu', 'kernel.cu'),
            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py', 'host_guard.py'),
            (geometry_review, 'geometry_review.json'), (geometry, 'reference_geometry.xyz'),
            (prior/'optimizer_result.json', 'optimizer_result.json'), (Path(__file__), 'preparer.py')]:
        shutil.copyfile(source, sources / name)
    if pilot_review is not None:
        shutil.copyfile(pilot_review, sources / 'thermochemistry_pilot_review.json')
    shutil.copyfile(geometry, case / 'state/opt.xyz')
    shutil.copyfile(prior/'optimizer_result.json', case / 'state/result.json')
    identifier = 'exact_viologen_chloride_' + review['key']
    with (case/'library/manifest.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['id', 'state', 'smiles', 'charge', 'mult_hint'])
        writer.writeheader()
        writer.writerow(dict(id=identifier, state='oxidized_contact', smiles=review['smiles'], charge=1, mult_hint=1))
    native = load(archive/'src/redox/qm/native_thermal.py', 'counterion_native_definition')
    thermal = load(archive/'src/redox/qm/qrrho.py', 'counterion_thermal_definition')
    if pilot_review is not None:
        if (pilot['electronic_protocol'] != ACTIVE_SP or pilot['native_protocol'] != native.NATIVE_THERMAL
                or pilot['thermal_protocol'] != thermal.PROTOCOL or pilot['smiles'] != review['smiles']):
            raise ValueError('Followup state or thermodynamic context differs from the pilot')
        for name in ('worker.py','hashing.py','checkpoints.py','helper.py','kernel.cu','host_guard.py'):
            if sha256_file(sources/name) != sha256_file(Path(pilot_spec['sources'])/name):
                raise ValueError('Followup calculation code differs from the qualified pilot')
    installed = gpu4pyscf_source_root()
    files = [*sources.glob('*'), case/'state/opt.xyz', case/'state/result.json', case/'library/manifest.csv',
        *(archive/'src/redox/qm').glob('*.py'), *(archive/'src/redox/core').glob('*.py'),
        *(archive/'config').glob('*.py'), *(archive/'config').glob('*.json'),
        archive/'src/redox/workflows/resources.py', archive/'scripts/pipeline/run_pi_stacking.py',
        archive/'scripts/validation/audit/native_partition.py', archive/'scripts/fleet/exclusive_python_worker.py',
        archive/'scripts/fleet/tensor_runtime_worker.py', archive/'environments/cutensor22/runtime.json',
        installed/'dft/gen_grid.py', installed/'dft/numint.py', installed/'lib/libgdft.so']
    plan = dict(role='primary_counterion_thermochemistry', source_id=review['key'], state='oxidized_contact',
        native_id=identifier, native_state='oxidized_contact', smiles=review['smiles'], charge=1, mult=1,
        parent_formal_core_charge=2, counterion_charge=-1, atoms=42, dimensions=dimensions,
        geometry_sha256=review['geometry_sha256'], electronic_protocol=ACTIVE_SP,
        native_protocol=native.NATIVE_THERMAL, thermal_protocol=thermal.PROTOCOL,
        native_memory_cap_GiB=6, maximum_host_RSS_GiB=40 if pilot_review is None else followup_host_bound(pilot), minimum_free_GiB=32,
        maximum_peak_GPU_MiB=31744, maximum_wall_seconds=14400, maximum_durable_output_bytes=2*2**30,
        minimum_initial_storage_GiB=34, minimum_live_storage_GiB=32,
        source_sha256={str(p): sha256_file(p) for p in files},
        pilot_review=None if pilot_review is None else str(sources/'thermochemistry_pilot_review.json'),
        geometry_repair_permitted=False, native_G_qualified=False, association_free_energy_computed=False,
        material_admission=False, experimental_populations_inferred=False,
        purpose='Exact primary molecular microstate component for a matched incoming-chloride association cycle; iodide remains an omitted separated spectator')
    atomic_write_json(case/'plan.json', plan)
    spec = dict(prior_root=str(ROOT), archive=str(archive), case=str(case), sources=str(sources),
        source_id=review['key'], state='oxidized_contact', maximum_host_RSS_GiB=plan['maximum_host_RSS_GiB'], minimum_free_GiB=32,
        scientific_sources_read_only=True, library_overlay=str(case/'library'), plan_sha256=sha256_file(case/'plan.json'))
    atomic_write_json(output/'namespace.json', spec)
    return spec


def launch(namespace, gpu, output):
    from redox.workflows.bounded_molecular import launch as execute
    spec = verify(namespace)
    plan = json.loads((Path(spec['case'])/'plan.json').read_text())
    if plan['role'] != 'primary_counterion_thermochemistry' or (plan['source_id'] != 'N0_above' and plan['pilot_review'] is None):
        raise ValueError('Exact prospective pilot or independently released thermochemistry input required')
    if plan['pilot_review'] is not None:
        pilot,_=checked_pilot(Path(plan['pilot_review']))
        if plan['maximum_host_RSS_GiB']!=followup_host_bound(pilot):raise ValueError('Declared host envelope differs from the measured pilot release')
    elif plan['maximum_host_RSS_GiB']!=40:
        raise ValueError('Unqualified pilot retains the original host envelope')
    return execute(namespace, gpu, output, caller_source=Path(__file__))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'launch', 'analyze'))
    for key in ('archive', 'geometry-review', 'pilot-review', 'out', 'spec', 'primary-inputs', 'geometry-manifest'):
        p.add_argument('--'+key, type=Path)
    p.add_argument('--gpu', type=int)
    p.add_argument('--component-review', type=Path, action='append', default=[])
    a = p.parse_args()
    if a.action == 'prepare':
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
        if any(v is None for v in (a.archive, a.geometry_review, a.out)):
            p.error('Preparation needs archive, checked geometry review and fresh output')
        r = prepare(a.archive.resolve(), a.geometry_review.resolve(), a.out.resolve(),
                    None if a.pilot_review is None else a.pilot_review.resolve())
        print(json.dumps(dict(status='primary_counterion_thermochemistry_prepared', namespace=str(a.out/'namespace.json'))))
    elif a.action == 'launch':
        if a.spec is None or a.out is None or a.gpu is None:
            p.error('Launch needs exact namespace, GPU and fresh placement')
        r = launch(a.spec.resolve(), a.gpu, a.out.resolve())
        print(json.dumps({k:r[k] for k in ('status', 'returncode', 'budget_stop')}))
    else:
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        from redox.workflows.counterion_association import analyze
        with allocation_guard():
            host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
        if a.primary_inputs is None or a.geometry_manifest is None or a.out is None:
            p.error('Analysis needs checked primary inputs, declared geometry manifest and fresh output')
        r = analyze(a.primary_inputs.resolve(), a.geometry_manifest.resolve(),
                    [p.resolve() for p in a.component_review], a.out.resolve())
        print(json.dumps({k:r[k] for k in ('status','declared_microstate_starts','computed_microstate_cycles')}))
