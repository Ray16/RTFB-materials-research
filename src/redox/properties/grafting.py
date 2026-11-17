"""Defined separated-ion quaternization cycles, distinct from HCl/base routes."""
import math
from collections import Counter

from rdkit import Chem

from redox.core.common import EV_KJ


def quaternization_delta_g(reaction, components, *, property_records=None):
    """Sum matched component G values; no base or concentration offset is fitted.

    This is a local-site thermodynamic model. It does not predict reaction rates,
    graft conversion, ion pairing, counterion exchange or polymer accessibility.
    Physical record admission is the caller's responsibility.
    """
    roles=('precursor','site','product','byproduct');cutoffs={'25','50','100'}
    if (reaction.get('reaction_class')!='N_quaternization' or reaction.get('base_correction_required') is not False
            or not reaction.get('atom_balanced') or not reaction.get('charge_balanced') or set(components)!=set(roles)):
        raise ValueError('A complete balanced quaternization cycle is required')
    canonical=lambda s:Chem.MolToSmiles(Chem.MolFromSmiles(s))
    for role in roles:
        row=components[role];expected=reaction[role+'_smiles'];mol=Chem.MolFromSmiles(expected)
        if (canonical(row['smiles'])!=canonical(expected) or type(row['charge']) is not int or row['charge']!=Chem.GetFormalCharge(mol)
                or set(row['G_by_cutoff_eV'])!=cutoffs
                or any(type(v) not in (int,float) or not math.isfinite(v) for v in row['G_by_cutoff_eV'].values())):
            raise ValueError('Reaction component identity, charge or complete free energy differs')
    def inventory(role):return Counter(a.GetAtomicNum() for a in Chem.AddHs(Chem.MolFromSmiles(components[role]['smiles'])).GetAtoms())
    if (inventory('product')+inventory('byproduct')!=inventory('precursor')+inventory('site')
            or components['product']['charge']+components['byproduct']['charge']!=components['precursor']['charge']+components['site']['charge']):
        raise ValueError('The actual component graphs must conserve every atom and total charge')
    if canonical(reaction['byproduct_smiles'])!='[Cl-]':raise ValueError('Explicit separated chloride required')
    matched = {}
    if property_records is not None:
        from redox.properties.context import combine_properties
        if set(property_records) != cutoffs:
            raise ValueError('Context-bound records required at every cutoff')
        for cutoff in cutoffs:
            records = property_records[cutoff]
            required = dict(quantity='G', unit='eV', solvent='acetonitrile', solvent_model='SMD',
                temperature_K=298.15, standard_state='1atm_translation',
                reference='absolute_component_free_energy', rotor_cutoff_cm=int(cutoff))
            if any(any(records[r]['context'].get(k) != v for k,v in required.items()) for r in records):
                raise ValueError('Reaction context must match the declared MeCN G convention')
            if set(records) != set(roles) or any(records[r]['value'] != components[r]['G_by_cutoff_eV'][cutoff] for r in roles):
                raise ValueError('Context-bound values differ from reaction components')
            if any(records[r]['state']['charge'] != components[r]['charge'] or records[r]['state']['multiplicity'] != components[r]['mult'] or canonical(records[r]['state']['identity']) != canonical(components[r]['smiles']) for r in roles):
                raise ValueError('Context-bound chemical states differ from reaction components')
            matched[cutoff] = combine_properties(records, dict(product=1, byproduct=1, precursor=-1, site=-1))
    values={c:(matched[c]['value'] if property_records is not None else math.fsum([
        components['product']['G_by_cutoff_eV'][c],components['byproduct']['G_by_cutoff_eV'][c],
        -components['precursor']['G_by_cutoff_eV'][c],-components['site']['G_by_cutoff_eV'][c]]))*EV_KJ
        for c in sorted(cutoffs,key=int)}
    return dict(deltaG_kJ_mol_by_cutoff=values,conversion_kJ_mol_per_eV=EV_KJ,
        property_context_verified=property_records is not None, property_context=matched,
        delta_n_solutes=0,standard_state_correction_kJ_mol=0.,base_term_kJ_mol=0.,
        solvent='acetonitrile',temperature_K=298.15,standard_concentration_M=1.,
        ion_model='separated ions; unmodeled iodide cancels as a separated spectator only',
        reaction_class='N_quaternization',material_admission=False,
        ion_pairing_computed=False,counterion_exchange_computed=False,experimental_reflux_modelled=False)
