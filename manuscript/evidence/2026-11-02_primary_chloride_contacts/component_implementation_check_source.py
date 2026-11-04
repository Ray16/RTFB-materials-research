"""Reconstruct candidate/reaction electronic implementation shifts for the SI."""
import json
from redox.validation.component_sensitivity import reconstruct


def audit(folder):
    inputs=json.loads((folder/'component_implementation_inputs.json').read_text())
    expected=json.loads((folder/'component_implementation_replay.json').read_text())
    actual={name:reconstruct(data['package'],data['energies']) for name,data in inputs.items()}
    if actual!=expected:raise ValueError('Electronic component sensitivity does not reconstruct')
    if actual['viologen']['states']!=10 or len(actual['viologen']['cycles'])!=6 or actual['base']['states']!=21 or len(actual['base']['cycles'])!=25:
        raise ValueError('Declared exact-state/reaction coverage differs')
    for report in actual.values():
        if report['nuclear_terms_recomputed'] or report['new_potentials_or_reaction_G_qualified']:
            raise ValueError('Electronic comparison promoted to full thermochemistry')
    pots=[r for r in actual['viologen']['cycles'] if r['unit']=='V']
    grafts=[r for r in actual['viologen']['cycles'] if r['kind']=='grafting_free_energy']
    largest=lambda rows,key:max(abs(r['electronic_contribution_to_shift'][key]) for r in rows)
    return dict(viologen_states=10,viologen_cycles=6,base_states=21,base_cycles=25,
        overlapping_reference_states_not_counted_as_new_unique_states=True,
        maximum_electronic_potential_shift_V=largest(pots,'corrected_minus_installed'),
        maximum_viologen_electronic_grafting_shift_kJ_mol=largest(grafts,'corrected_minus_installed'),
        maximum_base_electronic_shift_kJ_mol=largest(actual['base']['cycles'],'corrected_minus_installed'),
        maximum_base_repeat_vs_stored_shift_kJ_mol=largest(actual['base']['cycles'],'corrected_minus_stored'),
        nuclear_terms_recomputed=False,total_free_energies_requalified=False)
