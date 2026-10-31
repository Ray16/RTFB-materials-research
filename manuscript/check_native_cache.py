"""Bind the portable numerical cache to the paper's existing 75 native states."""
import hashlib
import json
import math


def audit(folder):
    load = lambda name: json.loads((folder/name).read_text())
    expected = [r for name in ('polymer_native_replay.json', 'polymer_native_continuations.json',
                'native_addition_five.json', 'native_addition_three.json') for r in load(name)['records']]
    replay = load('portable_native_cache_replay.json')
    isolation = load('portable_native_cache_isolation.json')
    compression = load('portable_native_cache_compression.json')
    manifest = load('portable_native_cache_manifest.json')
    if (replay['status'] != 'packaged_native_numerical_inputs_reconstructed' or replay['states'] != 75
            or len(expected) != 75 or len(replay['records']) != 75
            or len({r['key'] for r in replay['records']}) != 75
            or replay['quantum_calculations_run'] or replay['xTB_library_loaded']
            or replay['full_historical_source_closure_revalidated'] or replay['live_execution_revalidated']
            or replay['experimental_accuracy_inferred'] or replay['global_ensemble_convergence_proven']):
        raise ValueError('Portable native cache scope or coverage differs')
    originals = {r['key']: r for r in expected}
    for row in replay['records']:
        old = originals[row['key']]
        if row['configuration'] != old['configuration'] or row['occupied_sites'] != old['occupied_sites']:
            raise ValueError('Portable native cache changed chemical identity')
        values = [row['E_eV'], row['frequency_error_cm'], row['maximum_column_error'],
                  *row['G_by_cutoff_eV'].values(), *row['nuclear_errors_eV'].values()]
        if (not all(math.isfinite(v) for v in values) or abs(row['E_eV']-old['E_eV']) > 1e-9
                or set(row['G_by_cutoff_eV']) != {'25', '50', '100'}
                or set(row['nuclear_errors_eV']) != {'25', '50', '100'}
                or row['sampled_columns'] != old['sampled_columns']
                or any(abs(v-old['G_by_cutoff_eV'][k]) > 1e-8 for k,v in row['G_by_cutoff_eV'].items())
                or not 0 <= row['frequency_error_cm'] <= 1e-5
                or not 0 <= row['maximum_column_error'] <= .01
                or not 0 <= max(row['nuclear_errors_eV'].values()) <= 1e-8):
            raise ValueError('Portable raw quantities differ from the existing native evidence')
    if (isolation['status'] != 'extracted_native_cache_fresh_numpy_replay_checked' or isolation['states'] != 75
            or isolation['archive_sha256'] != compression['archive_sha256']
            or isolation['package_manifest_sha256'] != compression['manifest_sha256']
            or hashlib.sha256((folder/'portable_native_cache_manifest.json').read_bytes()).hexdigest() != compression['manifest_sha256']
            or isolation['original_project_installed'] or isolation['quantum_calculations_run']
            or isolation['full_historical_source_closure_revalidated'] or isolation['public_archive_deposited']
            or len(isolation['blocked_original_paths']) != 4
            or isolation['scientific_dependency'] != 'numpy==2.2.6'
            or isolation['scientific_requirements_sha256'] != manifest['files_sha256']['requirements.lock.txt']
            or isolation['archive_files'] != len(manifest['files_sha256'])+1
            or len(manifest['states']) != 75 or manifest['full_historical_source_closure_supplied']
            or manifest['external_software_or_scratch_revalidated']
            or {r['key'] for r in manifest['states']} != set(originals)):
        raise ValueError('Portable archive identity or fresh-environment scope differs')
    for field, value in [('maximum_frequency_error_cm', max(r['frequency_error_cm'] for r in replay['records'])),
                         ('maximum_nuclear_error_eV', max(max(r['nuclear_errors_eV'].values()) for r in replay['records']))]:
        if not math.isfinite(isolation[field]) or abs(isolation[field]-value) > 1e-16:
            raise ValueError('Portable reconstruction maxima do not reconstruct')
    return dict(states=75, extracted_archive_fresh_numpy_replay=True,
                full_historical_source_closure_supplied=False, quantum_calculations_run=False,
                experimental_accuracy_inferred=False, public_archive_deposited=False)
