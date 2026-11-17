"""Balanced neutral/singly oxidized molecular potentials against Fc/Fc+.

Physical QC and source verification belong to the input reader. Different
qualified numerical differentiation steps are retained, not converted into
different physical solvent models or silently discarded.
"""
from collections import Counter
from copy import deepcopy
import math

from redox.properties.context import combine_properties

COEFFICIENTS = dict(graft_ox=1, graft_neu=-1, fc_ox=-1, fc_neu=1)
CUTOFFS = ('25', '50', '100')


def neutral_oxidation_vs_fc(components):
    if set(components) != set(COEFFICIENTS):
        raise ValueError('Exactly four graft/Fc redox components required')
    elements = Counter(); charge = 0
    for key, row in components.items():
        expected = (1, 2) if key.endswith('_ox') else (0, 1)
        if ((row['state']['charge'], row['state']['multiplicity']) != expected
                or row.get('native_G_qualified') is not True
                or row.get('thermodynamic_entities') != 1
                or set(row['G_by_cutoff_eV']) != set(CUTOFFS)
                or not row.get('symbols')):
            raise ValueError('Qualified neutral-singlet/oxidized-doublet components required')
        policy = row.get('native_numerical_policy', {})
        step = policy.get('Hessian_gradient_step_Bohr')
        if (policy.get('component') != 'primary legacy SMD CDS'
                or type(step) not in (int, float) or not math.isfinite(step) or step <= 0
                or any(policy.get(k) is not False for k in
                       ('energy_function_changed', 'gradient_function_changed', 'global_defaults_changed'))):
            raise ValueError('Explicit unchanged-model CDS numerical policy required')
        for symbol in row['symbols']:
            elements[symbol] += COEFFICIENTS[key]
        charge += COEFFICIENTS[key] * row['state']['charge']
    if charge or any(elements.values()):
        raise ValueError('Atom and charge balance required')
    for key in ('fc_neu', 'fc_ox'):
        if Counter(components[key]['symbols']) != Counter(C=10, H=10, Fe=1):
            raise ValueError('Exact unsubstituted ferrocene reference required')
    if components['graft_neu']['state']['graph_identity'] != components['graft_ox']['state']['graph_identity']:
        raise ValueError('Exact same graft connectivity required')
    values, assemblies = {}, {}
    for cutoff in CUTOFFS:
        inputs = {}
        for key, row in components.items():
            item = deepcopy(row)
            item['value'] = item['G_by_cutoff_eV'][cutoff]
            item['context']['thermal_definition']['selected_rotor_cutoff_cm'] = int(cutoff)
            inputs[key] = item
        combined = combine_properties(inputs, COEFFICIENTS)
        if (combined['context']['quantity'] != 'composite_solution_G'
                or combined['context']['unit'] != 'eV'):
            raise ValueError('Absolute composite free energies in eV required')
        # G(graft+) + G(Fc) - G(graft) - G(Fc+), divided by one electron.
        values[cutoff] = combined['value']
        assemblies[cutoff] = combined
    return dict(status='qualified_component_finite_oxidation_cycle_assembled',
        definition='[G(graft+) + G(Fc) - G(graft) - G(Fc+)] / e',
        unit='V versus Fc/Fc+', electron_count=1, coefficients=dict(COEFFICIENTS),
        potential_by_cutoff_V=values, delta_n_solutes=0,
        one_atm_to_one_M_correction_V=0.0,
        components=deepcopy(components), property_assemblies=assemblies,
        numerical_policy_variation_retained=True,
        conformer_ensemble_qualified=False, experimental_accuracy_validated=False,
        material_admission=False)
