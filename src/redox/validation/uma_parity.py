"""Compare frozen UMA states without treating execution as numerical admission."""
import math

import numpy as np


def compare_points(plan, reference, candidate):
    """Require every planned electronic state and geometry, energy and force.

    Loader versions may differ; weights and both atomic reference files may not.
    Tolerances come from the pre-execution panel, never from measured residuals.
    A diagnostic behavior override is always ineligible for production adoption.
    """
    expected = plan['model_identity']
    for record in (reference, candidate):
        identity = record['model_identity']
        if any(identity[k] != expected[k] for k in ('model', 'checkpoint_sha256', 'reference_sha256', 'task_name')):
            raise ValueError('Checkpoint, atomic references or task differ')
        if record['status'] != 'UMA_point_evaluations_complete_cross_mode_review_pending':
            raise ValueError('Point evaluations are not complete')

    ids = [row['id'] for row in plan['cases']]
    if len(set(ids)) != len(ids) or not ids:
        raise ValueError('Frozen panel must contain unique cases')

    def by_id(record):
        rows = record['cases']
        if len(rows) != len(ids) or {r['id'] for r in rows} != set(ids):
            raise ValueError('Missing, duplicate or unexpected case')
        return {r['id']: r for r in rows}

    left, right = by_id(reference), by_id(candidate)
    comparisons = []
    for case in plan['cases']:
        a, b = left[case['id']], right[case['id']]
        for record in (a, b):
            if any(record[k] != case[k] for k in ('charge', 'mult', 'atoms')):
                raise ValueError('Electronic state or atomic inventory differs')
            if len(record['frames']) != len(case['frames']):
                raise ValueError('Missing or unexpected geometry frame')
        tolerance = max(plan['minimum_energy_tolerance_eV'],
                        case['atoms'] * plan['energy_tolerance_eV_per_atom'])
        frames = []
        for index, (x, y) in enumerate(zip(a['frames'], b['frames'])):
            energies = [x['energy_eV'], y['energy_eV']]
            f, g = np.asarray(x['forces_eV_A']), np.asarray(y['forces_eV_A'])
            if (not all(math.isfinite(e) for e in energies)
                    or f.shape != (case['atoms'], 3) or g.shape != f.shape
                    or not np.isfinite(f).all() or not np.isfinite(g).all()):
                raise ValueError('Invalid energy or force values')
            energy_error = abs(energies[1] - energies[0])
            rmse, maximum = float(np.sqrt(np.mean((g - f)**2))), float(np.max(abs(g - f)))
            checks = dict(energy=energy_error <= tolerance,
                          force_RMSE=rmse <= plan['force_RMSE_tolerance_eV_A'],
                          force_max=maximum <= plan['force_max_tolerance_eV_A'])
            frames.append(dict(frame=index, signed_energy_difference_eV=energies[1]-energies[0],
                energy_error_eV=energy_error, force_RMSE_eV_A=rmse,
                force_max_eV_A=maximum, checks=checks, passed=all(checks.values())))
        comparisons.append(dict(id=case['id'], charge=case['charge'], mult=case['mult'], atoms=case['atoms'],
            energy_tolerance_eV=tolerance, maximum_energy_error_eV=max(f['energy_error_eV'] for f in frames),
            maximum_force_RMSE_eV_A=max(f['force_RMSE_eV_A'] for f in frames),
            maximum_force_error_eV_A=max(f['force_max_eV_A'] for f in frames),
            passed=all(f['passed'] for f in frames), frames=frames))
    passed = all(row['passed'] for row in comparisons)
    return dict(status='numerical_comparison_passed' if passed else 'numerical_comparison_failed',
        all_energy_and_force_checks_pass=passed, cases=len(comparisons),
        frames=sum(len(r['frames']) for r in comparisons), passed_cases=sum(r['passed'] for r in comparisons),
        diagnostic_only=candidate.get('diagnostic_only', False), production_migration_qualified=False,
        scope='Frozen fixed-geometry energy/force comparison; optimizer, thermochemistry and speedup require separate checks',
        maximum_energy_error_eV=max(r['maximum_energy_error_eV'] for r in comparisons),
        maximum_force_error_eV_A=max(r['maximum_force_error_eV_A'] for r in comparisons),
        comparisons=comparisons)
