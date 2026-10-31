"""Reconstruct the frozen initial charged-state execution comparison."""
import json
import math
from pathlib import Path


def audit(directory):
    directory=Path(directory)
    a=json.loads((directory/'ionic_compact_review.json').read_text())
    b=json.loads((directory/'ionic_direct_review.json').read_text())
    expected=json.loads((directory/'ionic_execution_agreement.json').read_text())
    fields=('candidate_alias','graph_key','geometry_sha256','charge','mult','solvent_model','calculation_ion_convention')
    if any(a[k]!=b[k] for k in fields) or (a['atoms'],a['charge'],a['mult'])!=(137,4,1):
        raise ValueError('The exact declared resting ionic graph/state differs')
    if a['native_G_qualified'] or b['native_G_qualified'] or a['solvent_optimization_complete'] or b['solvent_optimization_complete']:
        raise ValueError('An initial execution pilot cannot imply thermochemistry or optimization')
    if any(len(row['gradient_Ha_Bohr'])!=137 or any(len(v)!=3 for v in row['gradient_Ha_Bohr']) for row in (a,b)):
        raise ValueError('Full atom-matched Cartesian gradients required')
    energy=b['E_Ha']-a['E_Ha'];force=max(abs(x-y) for aa,bb in zip(a['gradient_Ha_Bohr'],b['gradient_Ha_Bohr']) for x,y in zip(aa,bb))
    if (not all(math.isfinite(v) for v in [energy,force]) or abs(energy)>1e-7 or force>1e-6
            or energy!=expected['energy_difference_Ha'] or force!=expected['maximum_gradient_difference_Ha_Bohr']
            or not expected['initial_state_execution_agreement_passed']
            or (a['peak_GPU_MiB'],b['peak_GPU_MiB'])!=(29964,87026)):
        raise ValueError('Frozen ionic energy/force/resource comparison does not reconstruct')
    return dict(status='initial_ionic_execution_comparison_reconstructed',atoms=137,charge=4,mult=1,
        energy_difference_Ha=energy,maximum_gradient_difference_Ha_Bohr=force,
        compact_peak_GPU_MiB=a['peak_GPU_MiB'],direct_peak_GPU_MiB=b['peak_GPU_MiB'],
        native_G_qualified=False,polymer_interaction_strength_established=False)
