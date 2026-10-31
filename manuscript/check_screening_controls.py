"""Reconstruct the frozen fixed-density grid-order control without CUDA."""
import json
import math
from pathlib import Path


def audit(directory):
    directory = Path(directory)
    results = {}
    for arch in ('V100', 'H100'):
        path = directory / ('screening_' + arch + '_integration.json')
        data = json.loads(path.read_text())
        review = json.loads((directory / ('screening_' + arch + '_review.json')).read_text())
        if (data['native_G_qualified'] or review['native_G_qualified']
                or review['SCF_gradient_Hessian_qualification_complete']
                or review['status'] != 'screening_masks_and_fixed_density_integration_independently_checked'
                or review['upstream_absolute_limit'] != 2e-9):
            raise ValueError('Fixed-density numerical scope or upstream tolerance differs')
        errors = {}
        for row in data['rows']:
            values = row['values']
            if [v['sorted'] for v in values] != [True, False] or values[0]['grid_points'] != values[1]['grid_points']:
                raise ValueError('Matched sorted/unsorted integration records required')
            numbers = dict(density_integral_difference=values[0]['numeric_density_integral'] - values[1]['numeric_density_integral'],
                XC_energy_difference_Ha=values[0]['XC_energy_Ha'] - values[1]['XC_energy_Ha'])
            if any(not math.isfinite(v) or abs(v - review['integration_differences'][row['variant']][k]) > 1e-15
                   for k, v in numbers.items()):
                raise ValueError('Frozen control arithmetic does not reconstruct')
            errors[row['variant']] = max(abs(v) for v in numbers.values())
        if set(errors) != {'installed', 'corrected'} or errors['installed'] <= 2e-9 or errors['corrected'] > 2e-9:
            raise ValueError('Frozen grid-order control no longer supports its stated conclusion')
        results[arch] = errors
    return dict(status='frozen_fixed_density_screening_controls_reconstructed', variants=results,
                SCF_gradient_Hessian_qualification_complete=False, native_G_qualified=False)
