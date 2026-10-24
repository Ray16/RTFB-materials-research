"""Strict, fixed-state UMA Nitro calculator for gas-phase molecular warm starts.

The inference engine is checked before and after evaluation. Unsupported devices
or state changes fail explicitly; a non-Nitro fallback cannot enter these records.
"""
from __future__ import annotations

MODEL = 'uma-s-1p2p1'
TASK = 'omol'


def state_key(numbers, charge, multiplicity):
    if type(charge) is not int or type(multiplicity) is not int:
        raise ValueError('Explicit integer charge and spin multiplicity required')
    if not numbers or any(type(z) is not int or not 1 <= z <= 118 for z in numbers):
        raise ValueError('An ordered atomic-number inventory is required')
    electrons = sum(numbers) - charge
    if (not -100 <= charge <= 100 or not 1 <= multiplicity <= 100
            or electrons < multiplicity-1 or (electrons-(multiplicity-1)) % 2):
        raise ValueError('Charge/spin is incompatible with the declared electron count')
    return tuple(numbers), charge, multiplicity, TASK


def check_actual_nitro(predictor):
    if getattr(predictor.inference_settings, 'nitro', None) is not True:
        raise RuntimeError('Requested Nitro fell back to another inference mode')
    model = getattr(predictor.model, 'module', predictor.model)
    if 'uma_nitro' not in type(model).__module__:
        raise RuntimeError('The loaded model is not the standalone Nitro executor')
    return dict(actual_mode='nitro', executor=type(model).__module__+'.'+type(model).__name__)


def make_calculator(numbers, charge, multiplicity, model=MODEL, seed=41):
    if model != MODEL:
        raise ValueError('This execution profile is scoped to the exact UMA-S-1.2.1 checkpoint')
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    import numpy as np
    import torch
    from ase.calculators.calculator import Calculator, all_changes
    from fairchem.core import FAIRChemCalculator, pretrained_mlip

    key = state_key([int(z) for z in numbers], charge, multiplicity)
    capability = torch.cuda.get_device_capability()
    if capability[0] not in (8, 9, 10, 11, 12):
        raise RuntimeError('Nitro requires CUDA compute capability 8.x through 12.x')
    predictor = pretrained_mlip.get_predict_unit(model, device='cuda',
        inference_settings='nitro', seed=seed)
    check_actual_nitro(predictor)
    underlying = FAIRChemCalculator(predictor, task_name=TASK)

    class FixedStateNitro(Calculator):
        implemented_properties = ['energy', 'forces']

        def __init__(self):
            super().__init__()
            self.predictor = predictor
            self.execution = dict(requested_mode='nitro', task=TASK, model=model,
                charge=charge, multiplicity=multiplicity, capability=list(capability), seed=seed)

        def calculate(self, atoms=None, properties=('energy','forces'), system_changes=all_changes):
            if ('charge' not in atoms.info or 'spin' not in atoms.info
                    or state_key([int(z) for z in atoms.numbers], atoms.info['charge'], atoms.info['spin']) != key):
                raise ValueError('Nitro state changed; construct a new calculator for that state')
            check_actual_nitro(predictor)
            super().calculate(atoms, properties, system_changes)
            underlying.calculate(atoms, ['energy','forces'], system_changes)
            self.execution.update(check_actual_nitro(predictor))
            energy = float(underlying.results['energy'])
            forces = np.asarray(underlying.results['forces'], dtype=float)
            if not np.isfinite(energy) or forces.shape != (len(atoms),3) or not np.isfinite(forces).all():
                raise ValueError('Nitro returned invalid energy or forces')
            self.results = dict(energy=energy, forces=forces.copy())

    return FixedStateNitro()
