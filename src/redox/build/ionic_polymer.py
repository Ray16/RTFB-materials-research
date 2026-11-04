"""Explicit ionic occupancy graphs and dry-salt composition ceilings.

Graph construction, state hypotheses and a mass convention do not qualify QM,
experimental graft loading, ion speciation or a material selection.
"""
import math

from rdkit import Chem
from rdkit.Chem import Descriptors

from redox.build.grafting import precursor
from redox.build.polymer_cases import matched_graphs
from redox.properties.polymer_capacity import FARADAY_C_MOL,integer


def ionic_definition(fragment,model):
    molecule,reaction=precursor(fragment);attached=Chem.MolFromSmiles(fragment)
    charge=Chem.GetFormalCharge(attached)
    if reaction!='N_quaternization' or charge<=0:raise ValueError('This ionic model requires a positively charged quaternizing motif')
    if model.get('calculation_convention')!='bare_charged_polymer_with_separated_omitted_spectators':
        raise ValueError('Declare the calculation ion convention explicitly')
    if model.get('explicit_counterions_in_coordinates') is not False:
        raise ValueError('This preparation handles omitted separated ions only')
    reference=model.get('composition_reference',{})
    if reference.get('state')!='fully_oxidized_dry_salt' or not reference.get('label'):
        raise ValueError('Declare an oxidized dry-salt mass reference')
    ions=reference.get('counterions');counter_charge=0;counter_mass=0.;records=[]
    if not isinstance(ions,list) or not ions:raise ValueError('Explicit counterion identities/stoichiometries required')
    for ion in ions:
        number=integer(ion['per_graft'],'counterions per graft',1);mol=Chem.MolFromSmiles(ion['smiles'])
        if mol is None or len(Chem.GetMolFrags(mol))!=1 or Chem.GetFormalCharge(mol)>=0:
            raise ValueError('A connected anionic counterion is required')
        counter_charge+=number*Chem.GetFormalCharge(mol);counter_mass+=number*Descriptors.MolWt(mol)
        records.append(dict(smiles=Chem.MolToSmiles(mol),per_graft=number,charge=Chem.GetFormalCharge(mol)))
    if charge+counter_charge!=0:raise ValueError('The declared dry-salt composition is not electroneutral')
    return dict(charge_per_graft=charge,precursor_smiles=Chem.MolToSmiles(molecule),reaction_class=reaction,
        counterions=records,counterion_mass_per_graft_g_mol=counter_mass,composition_reference_label=reference['label'])


def capacity_bound(fragment,model,n,handles,grafts,electrons,additional_inactive_mass_g_mol=0.):
    n=integer(n,'repeat count',1);handles=integer(handles,'reactive handle count',0)
    grafts=integer(grafts,'graft count',0);electrons=integer(electrons,'electron budget',1)
    if not grafts<=handles<=n:raise ValueError('Grafts must occupy declared handles')
    if type(additional_inactive_mass_g_mol) not in (int,float) or not math.isfinite(additional_inactive_mass_g_mol) or additional_inactive_mass_g_mol<0:
        raise ValueError('Inactive mass must be finite and nonnegative')
    definition=ionic_definition(fragment,model);mol=Chem.MolFromSmiles(definition['precursor_smiles']);table=Chem.GetPeriodicTable()
    carbon,hydrogen,chloride=[table.GetAtomicWeight(i) for i in (6,1,17)]
    # Quaternization loses Cl from the site, without removing a precursor H.
    polymer=n*(8*carbon+8*hydrogen)+handles*(carbon+hydrogen+chloride)+2*hydrogen+grafts*(Descriptors.MolWt(mol)-chloride)
    counterion=grafts*definition['counterion_mass_per_graft_g_mol']
    total=polymer+counterion+additional_inactive_mass_g_mol
    return dict(polymer_cation_mass_g_mol=polymer,counterion_mass_g_mol=counterion,mass_g_mol=total,
        bare_polymer_charge=grafts*definition['charge_per_graft'],reference_salt_charge=0,
        capacity_upper_bound_mAh_g=grafts*electrons*FARADAY_C_MOL/(3.6*total),capacity_lower_bound_mAh_g=0.,
        electrons_per_graft=electrons,composition_reference_label=definition['composition_reference_label'],
        additional_inactive_mass_g_mol=additional_inactive_mass_g_mol,experimental_loading_inferred=False,
        accessible_capacity_established=False,mobile_ion_mass_during_redox_established=False)


def charged_panel(spec,fragment):
    if len(spec['reactive_sites'])!=len(set(spec['reactive_sites'])):raise ValueError('Reactive sites must be unique')
    model=spec['ionic_model'];definition=ionic_definition(fragment,model);cases,cycles=matched_graphs(spec,fragment)
    manifolds=model.get('proposed_global_states_by_graft_count',{})
    for case in cases:
        grafts=len(case['occupied_sites'])
        if case['charge']!=grafts*definition['charge_per_graft']:raise ValueError('Graph charge does not scale with its graft count')
        states=manifolds.get(str(grafts))
        if not isinstance(states,list) or not states or len({s['id'] for s in states})!=len(states):raise ValueError('Declare distinct global charge/spin hypotheses for each occupancy')
        numbers=sum(a.GetAtomicNum() for a in Chem.AddHs(Chem.MolFromSmiles(case['smiles'])).GetAtoms())
        for state in states:
            q,m=state['charge'],state['mult']
            if type(q) is not int or type(m) is not int or m<1 or m>numbers-q+1 or (numbers-q-(m-1))%2:
                raise ValueError('Global charge/spin is incompatible with the exact electron count')
        if not any(s['charge']==case['charge'] for s in states):raise ValueError('The resting graph charge is missing from its state hypotheses')
        case.update(proposed_global_states=states,charge_or_spin_localization_established=False,
            explicit_counterions_in_coordinates=False,calculation_ion_convention=model['calculation_convention'],
            composition_bounds={str(z):capacity_bound(fragment,model,spec['repeat_units'],len(spec['reactive_sites']),grafts,z) for z in spec['electrons_per_graft_scenarios']},
            quantum_calculation_admitted=False,native_G_qualified=False,material_admission=False)
    return cases,cycles,definition
