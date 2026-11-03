"""Check arithmetic/scope of the frozen independently reviewed ionic component."""
import json
import math


def audit(folder):
    row=json.loads((folder/'ionic_native_component_review.json').read_text())
    if (row['status']!='charged_polymer_native_component_independently_reconstructed' or not row['native_G_qualified']
            or row['graph_key']!='all_same__occ_2' or (row['charge'],row['mult'],row['atoms'])!=(2,1,114)
            or row['primary_model_admission'] or row['conformer_ensemble_qualified'] or row['material_admission']
            or not 0<=row['maximum_column_error']<=.01 or row['minimum_frequency_cm'] < -5
            or not 0<row['peak_GPU_MiB']<=31744):
        raise ValueError('Exact scoped native pilot qualification differs')
    for cutoff in ('25','50','100'):
        actual=row['E_geometry_level_Ha']*27.211386245988+row['nuclear_by_cutoff_eV'][cutoff]
        if not math.isfinite(actual) or abs(actual-row['G_by_cutoff_eV'][cutoff])>1e-8:
            raise ValueError('Pilot component electronic/nuclear sum differs')
    return dict(checked_components=1,atoms=114,charge=2,mult=1,
        minimum_internal_frequency_cm=row['minimum_frequency_cm'],maximum_sampled_column_error=row['maximum_column_error'],
        cached_component_arithmetic_replayed=True,original_full_physical_review_repeated=False,
        primary_model_admission=False,occupied_free_energy_cycle_qualified=False)
