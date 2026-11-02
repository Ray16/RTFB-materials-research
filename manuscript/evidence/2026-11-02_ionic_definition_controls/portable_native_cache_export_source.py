#!/usr/bin/env python3
"""Export or replay the raw numerical inputs of checked polymer native states.

Replay uses only the package and NumPy. It reconstructs saved Hessians,
displaced-gradient checks and thermochemistry; it does not rerun SCF, inspect
live GPU ownership, or verify every historical executable and scratch artifact.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            result.update(block)
    return result.hexdigest()


def local(root, relative):
    if Path(relative).is_absolute() or '..' in Path(relative).parts:
        raise ValueError('Portable artifact path must remain within the package')
    path = (root/relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Portable artifact escapes the package')
    return path


def function_source(path, names):
    text = path.read_text()
    functions = {node.name: ast.get_source_segment(text, node)
                 for node in ast.parse(text).body if isinstance(node, ast.FunctionDef)}
    return 'import math\nimport numpy as np\n\n\n'+'\n\n\n'.join(functions[name] for name in names)+'\n'


def prepare(inputs, prior, archive, output):
    from redox.validation.archived_sources import ArchivedSources
    from redox.validation.namespace_sources import NamespaceSources
    root = Path(__file__).resolve().parents[2]
    if output.exists():
        raise FileExistsError('Preserve earlier raw-state packages')
    output.mkdir(parents=True)
    objects = output/'objects'
    objects.mkdir()
    records, sources, keys = [], {}, set()

    def store(path, expected):
        if sha(path) != expected:
            raise ValueError('Source content changed before raw export: '+str(path))
        relative = 'objects/'+expected
        target = output/relative
        if not target.exists():
            shutil.copyfile(path, target)
        if sha(target) != expected or sha(path) != expected:
            raise ValueError('Raw artifact changed during export')
        return relative

    config = json.loads(inputs.read_text())
    sources[str(inputs)] = sha(inputs)
    for replay_path in map(Path, config['replays']):
        data = json.loads(replay_path.read_text())
        sources[str(replay_path)] = sha(replay_path)
        namespaces = {json.loads(Path(b['namespace']).read_text())['key']: b for b in data.get('bindings', [])}
        for row in data['records']:
            if row['key'] in keys or row['status'] != 'relocated_native_physical_replay_passed':
                raise ValueError('Distinct physically reconstructed source states required')
            keys.add(row['key'])
            if row['key'] in namespaces:
                binding = namespaces[row['key']]
                if sha(binding['namespace']) != binding['namespace_sha256']:
                    raise ValueError('Source namespace changed')
                reader = NamespaceSources(json.loads(Path(binding['namespace']).read_text()))
            else:
                reader = ArchivedSources(prior, archive, relative_to_prior=True)
            report_path = reader.verify(row['review'], row['review_sha256'])
            report = json.loads(report_path.read_text())
            if (report['task']['key'] != row['key'] or report['status'] != 'alternate_native_state_independently_reconstructed'
                    or report['primary_candidate_admission'] or report['experimental_model']['solvent'] != 'acetonitrile'):
                raise ValueError('Source review identity or scope differs')
            campaign = Path(report['campaign'])
            files, original = dict(review=store(report_path, row['review_sha256'])), {}

            def artifact(role, uri):
                uri = str(uri)
                digest = report['source_sha256'][uri]
                path = reader.verify(uri, digest)
                files[role] = store(path, digest)
                original[role] = dict(path=uri, sha256=digest)
                return path

            plan = json.loads(artifact('plan', campaign/'plan.json').read_text())
            state = json.loads(artifact('state', campaign/'experimental_native/status.json').read_text())
            raw = Path(state['raw_root'])
            for role, uri in [('partition', raw/'experimental_partition.json'), ('geometry', raw/'geometry.xyz'),
                              ('hessian', raw/'hessian.npy'), ('gradient', raw/'gradient.npy'),
                              ('symmetry', raw/'symmetry_thermo.log'), ('parent', Path(plan['electronic_review'])),
                              ('residency_owner', campaign/'residency_guard/owner.json'),
                              ('residency_state', campaign/'residency_guard/state.json'),
                              ('host_guard', campaign/'host_guard/status.json'),
                              ('tensor_runtime', campaign/'tensor_runtime/status.json')]:
                artifact(role, uri)
            if reader.path(campaign/'residency_guard/interruption.json').exists():
                raise ValueError('Interrupted source state cannot enter the checked raw cache')
            evaluations = sorted(reader.path(raw/'SCF_evaluations').glob('*/evaluation.json'))
            if not evaluations:
                raise ValueError('Saved displaced-gradient records required')
            for path in evaluations:
                name = path.parent.name
                if not name.isdigit():
                    raise ValueError('Unexpected saved evaluation name')
                artifact('evaluation_'+name, raw/'SCF_evaluations'/name/'evaluation.json')
            records.append(dict(key=row['key'], expected=row, original_raw_root=str(raw), files=files,
                                original_selected_sources=original, original_review_sha256=row['review_sha256']))
            reader.recheck()
    generated = {}
    for name, source, functions in [
        ('mode_checks.py', root/'scripts/validation/replay_archived_polymer_native.py',
         ['spin_valid', 'projected_modes', 'finite_difference_checks']),
        ('nuclear_equations.py', root/'scripts/validation/audit/native_partition.py', ['nuclear_equations'])]:
        body = function_source(source, functions)
        (output/name).write_text(body)
        generated[name] = dict(source=str(source), source_sha256=sha(source), functions=functions,
                               function_bodies_changed=False, generated_sha256=sha(output/name))
    shutil.copyfile(__file__, output/'replay.py')
    lock = root/'environments/reviewer_primary_cpu.lock.txt'
    numpy_line = [line for line in lock.read_text().splitlines() if line.startswith('numpy @ ')]
    if len(numpy_line) != 1 or '--hash=sha256:' not in numpy_line[0]:
        raise ValueError('An exact hash-locked NumPy wheel is required')
    (output/'requirements.lock.txt').write_text('# CPython 3.11 Linux x86_64; numerical replay only.\n'+numpy_line[0]+'\n')
    (output/'README.md').write_text(
        '# Portable native-state numerical cache\n\n'
        'Create a CPython 3.11 environment and install `pip install --require-hashes -r requirements.lock.txt`. '
        'Run `python -I replay.py replay --package . --out replay.json`. '
        'Only NumPy and files in this directory are used.\n\n'
        'The cache includes geometries, Cartesian Hessians, gradients, displaced-gradient SCF records, '
        'electronic parents, three-cutoff partitions, and saved execution reviews. Original URIs are descriptive '
        'metadata and are never opened during replay. Files are deduplicated by content hash.\n\n'
        'This repeats numerical checks on saved records, not quantum calculations or live execution. '
        'The complete historical source closure includes executable binaries and large scratch tensors that '
        'are not bundled here. Their original hashes remain in each unchanged original review. '
        'This cache does not claim their availability, a fresh GPU deployment, experimental accuracy, '
        'or global conformer convergence. Public licensing, metadata and deposition remain pending.\n')
    if any(sha(path) != digest for path, digest in sources.items()):
        raise ValueError('Source replay changed while exporting')
    manifest = dict(schema=1, status='portable_native_numerical_inputs_exported', states=records,
        source_replays_sha256=sources, exported_functions=generated,
        files_sha256={str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()},
        full_historical_source_closure_supplied=False, external_software_or_scratch_revalidated=False,
        public_archive_deposited=False, global_ensemble_convergence_proven=False)
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return dict(status=manifest['status'], states=len(records), files=len(manifest['files_sha256']))


def load_code(path, name):
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module


def replay(root):
    import numpy as np
    root = root.resolve()
    manifest = json.loads((root/'manifest.json').read_text())
    if (manifest['schema'] != 1 or manifest['full_historical_source_closure_supplied'] is not False
            or manifest['external_software_or_scratch_revalidated'] is not False):
        raise ValueError('Native numerical replay scope differs')
    for relative, expected in manifest['files_sha256'].items():
        if sha(local(root, relative)) != expected or (relative.startswith('objects/') and Path(relative).name != expected):
            raise ValueError('Portable artifact differs: '+relative)
    modes = load_code(local(root, 'mode_checks.py'), 'portable_native_mode_checks')
    equations = load_code(local(root, 'nuclear_equations.py'), 'portable_native_equations')
    rows, keys = [], set()
    for row in manifest['states']:
        files = row['files']
        def bound(role):
            path = local(root, files[role])
            original = row['original_selected_sources'].get(role)
            if original and sha(path) != original['sha256']:
                raise ValueError('Raw object differs from its original review binding')
            return path
        def read(role):
            return json.loads(bound(role).read_text())
        report, plan, state, native, parent = [read(role) for role in ('review', 'plan', 'state', 'partition', 'parent')]
        expected = row['expected']
        if (row['key'] in keys or row['key'] != report['task']['key'] or sha(bound('review')) != row['original_review_sha256']
                or report['status'] != 'alternate_native_state_independently_reconstructed'
                or report['primary_candidate_admission'] is not False
                or plan['tasks'] != [report['task']] or (report['task']['charge'], report['task']['mult']) != (0, 1)
                or report['experimental_model']['solvent'] != 'acetonitrile'
                or state['status'] != 'experimental_native_finished_independent_review_pending'
                or plan['experimental_model'] != report['experimental_model']
                or plan['experimental_method_id'] != report['experimental_method_id']
                or hashlib.sha256(json.dumps(report['experimental_model'], sort_keys=True).encode()).hexdigest() != report['experimental_method_id']
                or native['experimental_method_id'] != report['experimental_method_id']
                or native['native_settings'] != plan['experimental_native_settings']
                or native['symmetry_source'] != 'xTB autodetection on exact geometry'
                or not native['derivative_checks_passed'] or native['primary_candidate_admission'] is not False
                or not math.isfinite(native['spin_partition_eV']) or abs(native['spin_partition_eV']) > 1e-12
                or state['raw_root'] != row['original_raw_root']
                or expected['method_id'] != report['experimental_method_id']
                or expected['graph_sha256'] != report['task']['graph_sha256']
                or parent['status'] != 'alternate_electronic_geometry_and_resource_checked_native_G_pending'
                or parent['task'] != report['task'] or parent['experimental_model'] != report['experimental_model']
                or parent['portable_method_id'] != report['experimental_method_id']
                or not math.isfinite(expected['E_eV']) or abs(expected['E_eV']-parent['E_smd_eV']) > 1e-9
                or sha(bound('geometry')) != plan['scored_geometry_sha256']):
            raise ValueError('Exact raw native identity, model, geometry or scope differs')
        for role, binding in row['original_selected_sources'].items():
            if report['source_sha256'][binding['path']] != binding['sha256']:
                raise ValueError('Selected raw source does not belong to the original review')
        keys.add(row['key'])
        lines = bound('geometry').read_text().splitlines()
        count = int(lines[0]); atom_lines = [line.split() for line in lines[2:] if line.strip()]
        if len(atom_lines) != count or any(len(line) != 4 for line in atom_lines):
            raise ValueError('Malformed native coordinates')
        symbols = [line[0] for line in atom_lines]
        xyz = np.asarray([[float(x) for x in line[1:]] for line in atom_lines])
        matrix = np.load(bound('hessian'), allow_pickle=False)
        gradient = np.load(bound('gradient'), allow_pickle=False)
        if (matrix.shape != (3*count, 3*count) or gradient.shape != (count, 3)
                or not np.isfinite(matrix).all() or not np.isfinite(gradient).all() or not np.isfinite(xyz).all()
                or abs(matrix-matrix.T).max() > 1e-6 or abs(gradient).max() > .00045):
            raise ValueError('Saved Hessian or stationarity gate fails')
        evaluation_roles = {role.removeprefix('evaluation_'): role for role in files if role.startswith('evaluation_')}
        raw = Path(row['original_raw_root'])
        def evaluation(uri):
            relative = uri.relative_to(raw/'SCF_evaluations')
            if len(relative.parts) != 2 or relative.name != 'evaluation.json':
                raise ValueError('Unexpected displaced-gradient reference')
            return read(evaluation_roles[relative.parts[0]])
        checks, columns, required_count = modes.finite_difference_checks(raw, native, matrix, gradient,
            symbols, xyz/.529177210903, report['experimental_method_id'], evaluation)
        if len(evaluation_roles) != required_count:
            raise ValueError('Saved displaced-gradient inventory differs')
        nuclear = native['nuclear']; mass = np.asarray(nuclear['atomic_masses_amu'], dtype=float)
        if mass.shape != (count,) or not np.isfinite(mass).all() or min(mass) <= 0:
            raise ValueError('Invalid atom-matched masses')
        frequencies, rank = modes.projected_modes(matrix, xyz, mass)
        saved_frequencies = np.asarray(nuclear['frequencies_cm'])
        if saved_frequencies.shape != frequencies.shape or not np.isfinite(saved_frequencies).all():
            raise ValueError('Saved mode list has a different dimension or nonfinite value')
        error = float(np.max(abs(frequencies-saved_frequencies)))
        if rank != 6 or len(frequencies) != 3*count-6 or not np.isfinite(frequencies).all() or not math.isfinite(error) or error > 1e-5 or min(frequencies) < -5:
            raise ValueError('Saved internal-mode projection or minimum gate fails')
        sigma = int(re.search(r'rotational number\s+(\d+)', bound('symmetry').read_text())[1])
        if sigma < 1 or sigma != native['symmetry_number']:
            raise ValueError('Saved symmetry assignment differs')
        values, thermal_errors = {}, {}
        for rotor in ('25', '50', '100'):
            value = float(equations.nuclear_equations(frequencies, xyz, mass, sigma, int(rotor)))
            thermal_errors[rotor] = abs(value-nuclear['nuclear_by_cutoff_eV'][rotor])
            values[rotor] = float(parent['E_smd_eV']+value+native['spin_partition_eV'])
            stored_values = [nuclear['nuclear_by_cutoff_eV'][rotor], native['G_by_cutoff_eV'][rotor],
                             report['G_by_cutoff_eV'][rotor], expected['G_by_cutoff_eV'][rotor]]
            if (not all(math.isfinite(v) for v in stored_values) or not math.isfinite(parent['E_smd_eV'])
                    or not math.isfinite(values[rotor]) or not math.isfinite(value) or thermal_errors[rotor] > 1e-8
                    or abs(values[rotor]-native['G_by_cutoff_eV'][rotor]) > 1e-8
                    or abs(values[rotor]-report['G_by_cutoff_eV'][rotor]) > 1e-8
                    or abs(values[rotor]-expected['G_by_cutoff_eV'][rotor]) > 1e-8):
                raise ValueError('Reconstructed nuclear or free-energy term differs')
        rows.append(dict(key=row['key'], configuration=expected['configuration'], occupied_sites=expected['occupied_sites'],
                         G_by_cutoff_eV=values, E_eV=parent['E_smd_eV'], frequency_error_cm=error,
                         nuclear_errors_eV=thermal_errors, sampled_columns=len(columns), maximum_column_error=max(checks)))
    return dict(status='packaged_native_numerical_inputs_reconstructed', states=len(rows), records=rows,
                quantum_calculations_run=False, original_cluster_paths_opened=False, xTB_library_loaded=False,
                full_historical_source_closure_revalidated=False, live_execution_revalidated=False,
                experimental_accuracy_inferred=False, global_ensemble_convergence_proven=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'replay'))
    for name in ('inputs', 'prior', 'archive', 'package', 'out'):
        parser.add_argument('--'+name, type=Path)
    args = parser.parse_args()
    if args.action == 'prepare':
        from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
        import os
        with allocation_guard():
            host_capacity(1, minimum_ram_gib=4)
            os.sched_setaffinity(0, cpu_allocation(1))
        result = prepare(args.inputs.resolve(), args.prior.resolve(), args.archive.resolve(), args.out.resolve())
    else:
        result = replay(args.package.resolve())
        if args.out:
            with args.out.open('x') as stream:
                json.dump(result, stream, indent=2)
                stream.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'records'}))
