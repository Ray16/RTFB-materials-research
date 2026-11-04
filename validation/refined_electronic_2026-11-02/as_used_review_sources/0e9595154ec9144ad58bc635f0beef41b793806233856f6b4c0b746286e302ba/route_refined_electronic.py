"""First high-level scoring of independently checked phenothiazine refinements.

Reuses the installed/corrected molecular solver and actual-checkpoint reviewer.
No stored electronic reference exists at these new geometries. Both neutral
branches remain distinct inputs until their own electronic/native checks pass.
Top-level imports stay standard-library-only for the copied admission profile.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil

PROFILE = 'phenothiazine_refined_first_electronic_v1'
SMILES = 'Cc1ccc(CN2c3ccccc3Sc3ccccc32)cc1'
STATES = {'phenothiazine_neu_mode_minus': (0, 1),
          'phenothiazine_neu_mode_plus': (0, 1), 'phenothiazine_ox_stationarity': (1, 2)}
# Frozen resource-class declarations; checked against the maintained calculator
# during preparation and tests. No new workflow module is needed in the archive.
OUTPUT_BYTES = {1: 62500864, 2: 291119104}


def validate_plan(plan):
    from redox.core.protocol import ACTIVE_SP
    state = plan.get('source_state', {}); key = plan.get('source_id')
    expected = dict(role='screening_energy_control', energy_control=True,
        first_route_refined_electronic=True, route_electronic_profile=PROFILE,
        candidate_alias='phenothiazine', electronic_protocol=ACTIVE_SP,
        reference_energy_Ha=None, stored_energy_comparison_performed=False,
        atoms=39, nao=976, has_iron=False, maximum_host_RSS_GiB=40,
        minimum_free_GiB=32, minimum_GPU_memory_MiB=31744,
        maximum_peak_GPU_MiB=31744, maximum_wall_seconds=7200,
        minimum_live_storage_GiB=32, comparison_absolute_tolerance_Ha=1e-7,
        native_G_qualified=False, chemical_free_energy_assembled=False,
        redox_potential_qualified=False, material_admission=False)
    if (key not in STATES or any(plan.get(k) != v for k, v in expected.items())
            or (state.get('charge'), state.get('mult')) != STATES.get(key)
            or state.get('id') != key or state.get('smiles') != SMILES):
        raise ValueError('Exact new-geometry phenothiazine scoring with no historical energy reference required')
    budget = OUTPUT_BYTES[state['mult']]
    if (plan.get('maximum_durable_output_bytes') != budget
            or plan.get('minimum_initial_storage_GiB') != 32+budget/2**30):
        raise ValueError('The existing electronic retained-checkpoint resource envelope is required')
    return expected


def validate_inputs(plan, sources):
    validate_plan(plan)
    report = json.loads((sources/'geometry_review.json').read_text())
    gate = json.loads((sources/'route_gate.json').read_text())
    digest = hashlib.sha256((sources/'input.xyz').read_bytes()).hexdigest()
    if (report['status'] != 'route_primary_geometry_independently_checked'
            or report['key'] != plan['source_id'] or report['smiles'] != SMILES
            or (report['total_charge'], report['multiplicity']) != STATES[plan['source_id']]
            or report['geometry_sha256'] != digest or plan['geometry_sha256'] != digest
            or report['native_G_qualified'] or report['Hessian_class_qualified']
            or not gate['proceed_to_full_redox'] or not gate['exact_route_binding_verified']):
        raise ValueError('Exact independently checked refined geometry and prospective route gate required')
    return report


def prepare(archive, reviews, output):
    from ase.io import read
    from pyscf import gto
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    from redox.core.io import atomic_write_json
    from redox.core.protocol import ACTIVE_SP
    from redox.core.runtime_paths import gpu4pyscf_source_root
    from redox.validation.protonation import hydrogen_topology_status
    from redox.workflows.molecular_native_namespace import verify
    from redox.workflows.route_geometry import validate_plan as geometry_plan, visible_bindings
    from redox.workflows.screening_energy_panel import output_budget
    reports = [json.loads(p.read_text()) for p in reviews]
    if (len(reports) != 3 or {r['key'] for r in reports} != set(STATES)
            or output.exists() or any(output.is_relative_to(p) for p in (ROOT, archive))):
        raise ValueError('Retain both neutral branches and the oxidized state in a fresh durable campaign')
    for report in reports:
        if (report['status'] != 'route_primary_geometry_independently_checked'
                or report['native_G_qualified'] or any(sha256_file(p) != h for p, h in report['source_sha256'].items())):
            raise ValueError('Complete independently checked refined geometry backing required')
    installed = gpu4pyscf_source_root(); rows = []
    for review, report in zip(reviews, reports):
        prior = verify(Path(report['namespace'])); parent = Path(prior['case'])
        old = json.loads((parent/'plan.json').read_text()); geometry_plan(old)
        geometry = parent/'optimizer_return.xyz'; atoms = read(geometry)
        if (sha256_file(geometry) != report['geometry_sha256'] or atoms[15].symbol != 'S'
                or hydrogen_topology_status(SMILES, geometry) != 'ok'):
            raise ValueError('The accepted exact geometry and sulfur atom ordering must remain unchanged')
        q, m = STATES[report['key']]
        molecule = gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
                         basis=ACTIVE_SP['basis'], charge=q, spin=m-1, verbose=0)
        if (len(atoms), molecule.nao, report['total_charge'], report['multiplicity']) != (39, 976, q, m):
            raise ValueError('Exact refined molecular dimensions and electronic state required')
        root = output/'states'/report['key']; case = root/'case'; sources = root/'sources'
        case.mkdir(parents=True); sources.mkdir()
        for name in ('dft', 'locks', 'repairs'): (case/name).mkdir()
        for source, name in [(ROOT/'scripts/validation/run_screening_energy_control.py', 'worker.py'),
                (ROOT/'src/redox/qm/screening_barrier.py', 'helper.py'),
                (ROOT/'src/redox/qm/scf_checkpoints.py', 'checkpoints.py'),
                (ROOT/'tests/fixtures/screen_index_v1_8_1.cu', 'kernel.cu'),
                (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py', 'host_guard.py'),
                (Path(__file__), 'admission.py'), (geometry, 'input.xyz'), (review, 'geometry_review.json'),
                (Path(prior['sources'])/'route_gate.json', 'route_gate.json')]:
            shutil.copyfile(source, sources/name)
        evidence = visible_bindings({**report['source_sha256'], str(review): sha256_file(review)}, sources, ROOT)
        files = [*[p for p in sources.iterdir() if p.is_file()], *(archive/'src/redox/core').glob('*.py'),
            archive/'src/redox/workflows/resources.py', archive/'src/redox/qm/dft.py',
            archive/'scripts/pipeline/run_pi_stacking.py', archive/'scripts/fleet/exclusive_python_worker.py',
            archive/'scripts/fleet/tensor_runtime_worker.py', archive/'environments/cutensor22/runtime.json',
            installed/'dft/gen_grid.py', installed/'dft/numint.py', installed/'lib/libgdft.so']
        budget = output_budget(976, m, False)
        if budget != OUTPUT_BYTES[m]: raise ValueError('The frozen electronic resource class differs from its maintained budget definition')
        plan = dict(role='screening_energy_control', energy_control=True, first_route_refined_electronic=True,
            route_electronic_profile=PROFILE, candidate_alias='phenothiazine', source_id=report['key'], state='rest',
            source_state=dict(id=report['key'], aliases=['phenothiazine:'+report['key']], smiles=SMILES,
                              charge=q, mult=m, original_state_directory=str(parent)),
            geometry_sha256=sha256_file(geometry), electronic_protocol=ACTIVE_SP, reference_energy_Ha=None,
            stored_energy_comparison_performed=False, atoms=39, nao=976, has_iron=False,
            maximum_host_RSS_GiB=40, minimum_free_GiB=32, minimum_GPU_memory_MiB=31744,
            maximum_peak_GPU_MiB=31744, maximum_wall_seconds=7200, maximum_durable_output_bytes=budget,
            minimum_initial_storage_GiB=32+budget/2**30, minimum_live_storage_GiB=32,
            comparison_absolute_tolerance_Ha=1e-7, native_G_qualified=False, chemical_free_energy_assembled=False,
            redox_potential_qualified=False, material_admission=False,
            purpose='First scoring at the accepted refined geometry; installed/corrected comparison only; no stored-geometry energy or nuclear transfer',
            source_sha256={**evidence, **{str(p):sha256_file(p) for p in files}})
        validate_inputs(plan, sources); atomic_write_json(case/'plan.json', plan)
        spec = dict(prior_root=str(ROOT), archive=str(archive), case=str(case), sources=str(sources),
            source_id=report['key'], state='rest', maximum_host_RSS_GiB=40, minimum_free_GiB=32,
            scientific_sources_read_only=True, plan_sha256=sha256_file(case/'plan.json'))
        atomic_write_json(root/'namespace.json', spec); verify(root/'namespace.json')
        rows.append(dict(key=report['key'], namespace=str(root/'namespace.json'),
                         namespace_sha256=sha256_file(root/'namespace.json'),
                         launcher_module='redox.workflows.route_refined_electronic'))
    result = dict(status='three_refined_first_electronic_states_prepared', states=rows,
        stored_energy_comparison_performed=False, scientific_submissions=0, native_G_qualified=False)
    atomic_write_json(output/'manifest.json', result); return result


def qualified_rows(manifest):
    from redox.core.hashing import sha256_file
    from redox.workflows.molecular_native_namespace import verify
    rows = json.loads(manifest.read_text())['states']
    if len(rows) != 3 or {r['key'] for r in rows} != set(STATES):
        raise ValueError('Both neutral branches and the oxidized input must remain explicit')
    for row in rows:
        if (row.get('launcher_module') != 'redox.workflows.route_refined_electronic'
                or sha256_file(row['namespace']) != row['namespace_sha256']):
            raise ValueError('Exact prepared namespace and maintained scoring launcher required')
        spec = verify(Path(row['namespace']))
        validate_inputs(json.loads((Path(spec['case'])/'plan.json').read_text()), Path(spec['sources']))
    return rows


def launch(namespace, gpu, output):
    """Keep the existing budget observer alive for the owned scientific child."""
    from redox.workflows.molecular_native_namespace import verify
    from redox.workflows.screening_control import launch as start, observe_budget
    spec=verify(namespace)
    validate_inputs(json.loads((Path(spec['case'])/'plan.json').read_text()),Path(spec['sources']))
    child,record=start(namespace,gpu,output)
    return observe_budget(child,spec,record,output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('action', choices=('prepare','launch','queue'))
    for name in ('archive','out','spec','manifest','storage-peers'): p.add_argument('--'+name, type=Path)
    p.add_argument('--geometry-review', type=Path, action='append', default=[])
    p.add_argument('--gpu', type=int); p.add_argument('--python'); p.add_argument('--hours', type=float, default=12)
    p.add_argument('--resume', action='store_true'); a = p.parse_args()
    if a.action != 'launch':
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard(): host_capacity(1,minimum_ram_gib=4); os.sched_setaffinity(0,cpu_allocation(1))
    if a.action == 'prepare':
        if a.archive is None or a.out is None or not a.geometry_review: p.error('Archive, all geometry reviews and fresh output required')
        r = prepare(a.archive.resolve(), [f.resolve() for f in a.geometry_review], a.out.resolve())
    elif a.action == 'launch':
        if a.spec is None or a.gpu is None or a.out is None: p.error('Exact namespace, GPU and placement required')
        r=launch(a.spec.resolve(),a.gpu,a.out.resolve())
    else:
        from redox.workflows.prepared_namespace_queue import run
        if any(v is None for v in (a.manifest,a.out,a.storage_peers,a.python)): p.error('Manifest, output, storage roster and interpreter required')
        manifest=a.manifest.resolve();rows=qualified_rows(manifest)
        r=run(rows,a.out.resolve(),a.storage_peers.resolve(),a.python,a.hours,
            finished_status='screening_energy_control_finished_review_pending',
            validate=lambda row:qualified_rows(manifest),sources=[Path(__file__),manifest],
            resume=a.resume,label='route_refined_electronic')
    print(json.dumps(dict(status=r['status'],native_G_qualified=False)))
