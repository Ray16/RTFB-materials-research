"""Prepare the complete corrected sampling cohort with compact SCF evidence.

The existing gas/SMD hybrid single-point equations and model are unchanged.
Old geometries are explicit starting inputs; historical E/G never enter this
new cohort. Empty, single and double grafts require separate resource pilots.
"""
import argparse
import hashlib
import json
import os
import math
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP


def profile(task):
    occupied = len(task['occupied_sites'])
    if occupied not in (0, 1, 2) or (task['charge'], task['mult']) != (0, 1):
        raise ValueError('Only the declared neutral empty/single/double controls are supported')
    return dict(maximum_host_RSS_GiB=240 if occupied == 2 else 160,
        host_memory_profile='corrected_hybrid_energy',
        minimum_free_GiB=32, maximum_peak_GPU_MiB=31744,
        minimum_GPU_memory_MiB=31 * 1024, maximum_wall_seconds=14400,
        maximum_durable_output_bytes=256 * 2**20,
        minimum_initial_storage_GiB=32.25, minimum_live_storage_GiB=32)


def validate_inventory(data):
    rows = data['states']
    if (data['status'] != 'complete_corrected_sampling_start_inventory'
            or data['declared_inputs'] != 240 or data['selected_inputs'] != 240
            or data['missing'] or data['historical_values_transferred']
            or data['experimental_occupancy_weights'] or len(rows) != 240
            or len({r['key'] for r in rows}) != 240):
        raise ValueError('The full source-bound twelve-by-twenty roster is required')
    cells = {}
    for row in rows:
        task = row['task']; profile(task)
        if (task['key'] != row['key'] or task['candidate_id'] != 'exp_phthalimide_control'
                or task['role'] != 'interaction_method_control; not synthesis admission'
                or task['stereoisomer'] not in ('all_same', 'alternating')):
            raise ValueError('Exact declared graph/configuration and method-control identity required')
        cells.setdefault((task['stereoisomer'], tuple(task['occupied_sites'])), []).append(row['key'])
    expected = {(cfg, sites) for cfg in ('all_same', 'alternating')
                for sites in ((), (2,), (3,), (4,), (2, 3), (2, 4))}
    if set(cells) != expected or any(len(keys) != 20 for keys in cells.values()):
        raise ValueError('Every declared graph/configuration needs all twenty starting IDs')
    return rows


def prepare(inventory_path, archive, output):
    from redox.build.geometry import geometry_stereo_status
    from redox.validation.protonation import hydrogen_topology_status
    from redox.core.runtime_paths import gpu4pyscf_source_root
    data = json.loads(inventory_path.read_text()); rows = validate_inventory(data)
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT, archive)):
        raise FileExistsError('Fresh durable output outside both scientific checkouts required')
    for path, digest in data['source_sha256'].items():
        if sha256_file(path) != digest:
            raise ValueError('Sampling source inventory changed: ' + path)
    template_path = archive / 'results/polymer_models/campaigns/2026-10-25_smooth_energy_point80/all_same__occ_none__conf001/plan.json'
    template = json.loads(template_path.read_text())
    model = template['experimental_model']; method_id = template['experimental_method_id']
    if (hashlib.sha256(json.dumps(model, sort_keys=True).encode()).hexdigest() != method_id
            or model['solvent'] != 'acetonitrile'):
        raise ValueError('The declared original smooth-surface MeCN model is required')
    package = gpu4pyscf_source_root()
    # Verify every actual geometry before creating namespaces. The prior native
    # review is not used as a substitute for exact graph/stereochemistry here.
    for row in rows:
        xyz = Path(row['geometry']); task = row['task']
        if (sha256_file(xyz) != row['geometry_sha256']
                or hydrogen_topology_status(task['smiles'], xyz) != 'ok'
                or geometry_stereo_status(task['smiles'], xyz) not in ('ok', 'not_applicable')):
            raise ValueError('Actual sampling XYZ differs from the declared graph/stereo: ' + row['key'])
    output.mkdir(parents=True); prepared = []
    # The worker's source view overlays ROOT with the archived checkout. Bind
    # new inputs from this durable namespace, where they remain visible.
    inputs = output / 'inputs'; inputs.mkdir()
    shutil.copyfile(inventory_path, inputs / 'inventory.json')
    for row in rows:
        task = row['task']; root = output / 'states' / row['key']
        case = root / 'case'; sources = root / 'sources'
        case.mkdir(parents=True); sources.mkdir()
        for original, name in [
            (ROOT / 'scripts/validation/run_corrected_sampling_energy.py', 'worker.py'),
            (ROOT / 'src/redox/qm/occupied_checkpoints.py', 'occupied.py'),
            (ROOT / 'src/redox/qm/screening_barrier.py', 'helper.py'),
            (ROOT / 'tests/fixtures/screen_index_v1_8_1.cu', 'kernel.cu'),
            (ROOT / 'scripts/fleet/bounded_peak_host_python_worker.py', 'host_guard.py'),
            (Path(__file__), 'preparer.py')]:
            shutil.copyfile(original, sources / name)
        shutil.copyfile(row['geometry'], case / 'input.xyz')
        atomic_write_json(sources / 'starting_geometry.json', row)
        bound = [*sources.glob('*'), case / 'input.xyz', inputs / 'inventory.json', template_path,
            *(archive / 'src/redox/qm').glob('*.py'), *(archive / 'src/redox/core').glob('*.py'),
            *(archive / 'config').glob('*.py'), *(archive / 'config').glob('*.json'),
            archive / 'src/redox/workflows/polymer_campaign.py', archive / 'src/redox/workflows/resources.py',
            archive / 'scripts/pipeline/run_pi_stacking.py', archive / 'scripts/pipeline/smooth_cds_control.py',
            archive / 'scripts/pipeline/smd_gc_workspace.py', archive / 'scripts/pipeline/smooth_polymer_energy.py',
            archive / 'scripts/pipeline/smooth_polymer_energy_bounded.py',
            archive / 'scripts/fleet/exclusive_python_worker.py', archive / 'scripts/fleet/tensor_runtime_worker.py',
            archive / 'environments/cutensor22/runtime.json', package / 'dft/gen_grid.py',
            package / 'dft/numint.py', package / 'lib/libgdft.so']
        hashes = {str(p): sha256_file(p) for p in bound}
        plan = dict(schema='corrected_sampling_energy.v1', role='corrected_sampling_energy',
            corrected_sampling_energy=True, source_id=row['key'], state='rest', tasks=[task],
            raw_root=str(case / 'raw'), scored_geometry_sha256=row['geometry_sha256'],
            solvent='acetonitrile', experimental_model=model, experimental_method_id=method_id,
            electronic_protocol=ACTIVE_SP,
            quantum_settings=template['quantum_settings'], runtime_profile=template['runtime_profile'],
            producer_sha256=hashes, input_sha256={}, source_sha256=hashes,
            energy_only_execution=True, geometry_only_execution=False, native_pipeline_enabled=False,
            native_G_computed=False, primary_candidate_admission=False, prospective_resource_pilot=True,
            historical_energy_or_thermal_values_transferred=False, experimental_occupancy_weights=False,
            checkpoint_format='redox.restricted_occupied_SCF.v1', checkpoint_virtual_coefficients=False,
            checkpoint_density_restart=True, checkpoint_response_restart=False,
            purpose='Complete corrected same-model gas/SMD electronic cohort; native/ensemble qualification remains separate',
            **profile(task))
        atomic_write_json(case / 'plan.json', plan)
        spec = dict(prior_root=str(ROOT), archive=str(archive), case=str(case), sources=str(sources),
            source_id=row['key'], state='rest', scientific_sources_read_only=True,
            maximum_host_RSS_GiB=plan['maximum_host_RSS_GiB'], minimum_free_GiB=32,
            plan_sha256=sha256_file(case / 'plan.json'))
        atomic_write_json(root / 'namespace.json', spec)
        prepared.append(dict(key=row['key'], namespace=str(root / 'namespace.json'),
            namespace_sha256=sha256_file(root / 'namespace.json'),
            occupied_count=len(task['occupied_sites']), nao=task['active_basis_functions']))
    pilots = [min([r for r in prepared if r['occupied_count'] == n], key=lambda r: (-r['nao'], r['key']))['key'] for n in (0, 1, 2)]
    manifest = dict(status='complete_corrected_sampling_energy_cohort_prepared_no_jobs_submitted',
        target_states=240, states=prepared, pilot_keys=pilots, experimental_method_id=method_id,
        source_inventory=str(inventory_path), source_inventory_sha256=sha256_file(inventory_path),
        quantum_jobs_submitted=0, native_G_qualified=False, occupied_ensemble_cycles_qualified=False,
        historical_values_transferred=False, primary_candidate_admission=False,
        resource_release='Independent corrected gas/SMD and compact-checkpoint reviews for all three occupancy-size pilots required')
    atomic_write_json(output / 'manifest.json', manifest)
    return manifest


def dispatch(manifest_path, interpreter, reviews=(), pilots_only=False, storage_peers=()):
    from redox.workflows.molecular_native_namespace import verify
    from redox.workflows.charged_native_followup import dispatch_rows
    data = json.loads(manifest_path.read_text()); rows = data['states']
    inventory_path=Path(data['source_inventory'])
    if sha256_file(inventory_path)!=data['source_inventory_sha256']:
        raise ValueError('The complete declared sampling source inventory changed')
    inputs={r['key']:r for r in validate_inventory(json.loads(inventory_path.read_text()))}
    if data['target_states'] != 240 or len(rows) != 240 or len({r['key'] for r in rows}) != 240:
        raise ValueError('The full declared sampling energy cohort must remain intact')
    expected = {min([r for r in rows if r['occupied_count'] == n], key=lambda r: (-r['nao'], r['key']))['key'] for n in (0, 1, 2)}
    if expected != set(data['pilot_keys']):
        raise ValueError('Separate exact empty/single/double size pilots required')
    selected = [r for r in rows if r['key'] in expected] if pilots_only else rows
    checked = [json.loads(path.read_text()) for path in reviews]
    pilot_plans = {}
    if not pilots_only:
        if len(checked) != 3 or {r['key'] for r in checked} != expected:
            raise ValueError('All three exact independently reviewed pilot receipts required')
        for report in checked:
            row = next(r for r in rows if r['key'] == report['key'])
            bound=profile(inputs[row['key']]['task'])['maximum_host_RSS_GiB']
            if (report['status'] != 'corrected_sampling_energy_independently_checked'
                    or report['namespace'] != row['namespace'] or not report['resources_checked']
                    or not report['actual_checkpoint_frames_verified'] or not report['occupied_wavefunctions_checked']
                    or report['native_G_qualified'] or report['primary_candidate_admission']
                    or report['portable_method_id'] != data['experimental_method_id']
                    or not isinstance(report['peak_host_RSS_bytes'],(int,float))
                    or not math.isfinite(report['peak_host_RSS_bytes'])
                    or not 0<1.25*report['peak_host_RSS_bytes']<=bound*2**30
                    or not isinstance(report['peak_GPU_MiB'],(int,float))
                    or not math.isfinite(report['peak_GPU_MiB']) or not 0<report['peak_GPU_MiB']<=31744
                    or not report['source_sha256']
                    or any(sha256_file(p) != h for p, h in report['source_sha256'].items())):
                raise ValueError('Complete same-model physical/resource pilot review required')
            spec = verify(Path(row['namespace']))
            pilot_plans[row['occupied_count']] = (spec, json.loads((Path(spec['case']) / 'plan.json').read_text()))
    for row in selected:
        if sha256_file(row['namespace']) != row['namespace_sha256']:
            raise ValueError('Sampling namespace changed')
        spec = verify(Path(row['namespace'])); plan = json.loads((Path(spec['case']) / 'plan.json').read_text())
        if (plan['role'] != 'corrected_sampling_energy' or len(plan['tasks']) != 1
                or plan['tasks'][0]!=inputs[row['key']]['task']
                or plan['scored_geometry_sha256']!=inputs[row['key']]['geometry_sha256']
                or plan['electronic_protocol'] != ACTIVE_SP
                or plan['experimental_method_id'] != data['experimental_method_id']
                or plan['tasks'][0]['active_basis_functions'] != row['nao']
                or any(plan.get(k) != v for k, v in profile(plan['tasks'][0]).items())):
            raise ValueError('Prepared energy method/dimension/resource profile changed')
        if not pilots_only:
            ps, pp = pilot_plans[row['occupied_count']]
            if (row['nao'] > pp['tasks'][0]['active_basis_functions']
                    or plan['experimental_model'] != pp['experimental_model']):
                raise ValueError('Followup exceeds its declared pilot dimension/model class')
            for name in ('worker.py', 'occupied.py', 'helper.py', 'kernel.cu', 'host_guard.py'):
                if sha256_file(Path(spec['sources']) / name) != sha256_file(Path(ps['sources']) / name):
                    raise ValueError('Followup execution differs from the checked pilots')
    return dispatch_rows(selected, manifest_path.parent, interpreter,
        [manifest_path, *reviews, Path(__file__)], storage_peers=storage_peers)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', nargs='?', choices=('prepare', 'dispatch'), default='prepare')
    for name in ('archive', 'inventory', 'out', 'manifest', 'storage-peers'):
        p.add_argument('--' + name, type=Path)
    p.add_argument('--python'); p.add_argument('--pilots-only', action='store_true')
    p.add_argument('--pilot-review', type=Path, action='append', default=[])
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    if a.action == 'prepare':
        if any(v is None for v in (a.inventory, a.archive, a.out)):
            p.error('Preparation requires the complete inventory, archive and fresh output')
        result = prepare(a.inventory.resolve(), a.archive.resolve(), a.out.resolve())
        print(json.dumps({k: result[k] for k in ('status', 'target_states', 'pilot_keys')}))
    else:
        if a.manifest is None or not a.python:
            p.error('Dispatch requires the prepared manifest and interpreter')
        peers = [] if a.storage_peers is None else json.loads(a.storage_peers.read_text())
        result = dispatch(a.manifest.resolve(), a.python, [p.resolve() for p in a.pilot_review], a.pilots_only, peers)
        print(json.dumps(dict(status=result['status'], launches=len(result['launches']), refusals=result['refusals'], unsubmitted=result['unsubmitted'])))
