"""Reconstruct the current sampling figure and the tested backend limitation."""
from collections import Counter
import hashlib
import json
import math


def audit(folder):
    load = lambda name: json.loads((folder/name).read_text())
    previous = load('polymer_native_replay.json')
    new = load('polymer_native_continuations.json')
    figure = load('polymer_native_sampling_current.json')
    if (new['status'] != 'namespace_native_physical_records_replayed' or new['states'] != 7
            or new['source_bindings'] != 10937 or new['quantum_calculations_run']
            or new['primary_candidate_admission'] or new['global_ensemble_convergence_proven']):
        raise ValueError('Continuation reconstruction scope differs')
    rows = previous['records']+new['records']
    expected = {('all_same', ()): 20, ('alternating', ()): 16,
                ('all_same', (2,)): 18, ('alternating', (2,)): 13}
    if (len(rows) != 67 or len({r['key'] for r in rows}) != 67
            or len({r['method_id'] for r in rows}) != 1
            or Counter((r['configuration'], tuple(r['occupied_sites'])) for r in rows) != expected
            or len(figure['panels']) != 4 or figure['global_convergence_proven']):
        raise ValueError('Combined coverage, identities or solvent model differ')
    for row in new['records']:
        if (row['status'] != 'relocated_native_physical_replay_passed'
                or not math.isfinite(row['minimum_frequency_cm']) or row['minimum_frequency_cm'] < -5
                or not 0 <= row['maximum_column_error'] <= .01
                or not 0 <= row['frequency_error_cm'] <= 1e-5
                or set(row['G_by_cutoff_eV']) != {'25', '50', '100'}
                or any(not math.isfinite(v) for v in row['G_by_cutoff_eV'].values())
                or any(not 0 <= v <= 1e-8 for v in row['nuclear_errors_eV'].values())):
            raise ValueError('Continuation physical reconstruction fails')
    conv = 1.602176634e-19*6.02214076e23/1000
    for panel in figure['panels']:
        key = panel['configuration'], tuple(panel['occupied_sites'])
        subset = [r for r in rows if (r['configuration'], tuple(r['occupied_sites'])) == key]
        if not subset or len({r['graph_sha256'] for r in subset}) != 1:
            raise ValueError('One exact chemical graph per panel required')
        indices = [int(r['key'].split('conf')[-1]) for r in subset]
        values = [r['G_by_cutoff_eV']['50'] for r in subset]
        early = [v for i, v in zip(indices, values) if i <= 10]
        change = (min(values)-min(early))*conv
        if (panel['checked'] != len(subset) or panel['missing_indices'] != sorted(set(range(1, 21))-set(indices))
                or abs(panel['observed_range_kJ_mol']-(max(values)-min(values))*conv) > 1e-8
                or abs(panel['lowest_checked_minus_lowest_checked_in_first_ten_kJ_mol']-change) > 1e-8):
            raise ValueError('Current sampling figure does not reconstruct')
    digest = hashlib.sha256((folder/'polymer_native_sampling_current.pdf').read_bytes()).hexdigest()
    if digest != next(v for p, v in figure['outputs_sha256'].items() if p.endswith('.pdf')):
        raise ValueError('Current sampling figure changed')
    if figure['continuation_sha256'] != hashlib.sha256((folder/'polymer_native_continuations.json').read_bytes()).hexdigest():
        raise ValueError('Current figure has a different continuation source')
    nitro = load('nitro_repeatability.json')
    standard, tested = [nitro['modes'][mode] for mode in ('reference', 'nitro')]
    if (nitro['production_calculator_qualified'] or nitro['full_panel_qualified']
            or nitro['atoms'] != 42 or not standard['repeat_gate_passed'] or tested['repeat_gate_passed']
            or standard['checkpoint_sha256'] != tested['checkpoint_sha256']
            or not standard['initial_to_final_force_error_eV_A'] <= 1e-4
            or not tested['initial_to_final_force_error_eV_A'] > tested['force_tolerance_eV_A'] == 1e-4):
        raise ValueError('Tested Nitro limitation differs')
    return dict(checked_states=67, additional_states=7, additional_source_bindings=10937,
                complete_occupied_batch_convergence_proven=False,
                nitro_production_qualified=False,
                nitro_force_repeat_error_eV_A=tested['initial_to_final_force_error_eV_A'])
