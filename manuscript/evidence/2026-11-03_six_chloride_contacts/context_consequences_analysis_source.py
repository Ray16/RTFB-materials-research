"""Reconstruct consequences of explicitly stated model transfers from saved data.

These are controlled interpretations of existing quantities, not new quantum
results, experimental decision accuracy or a comparison against a prior algorithm.
"""
import math
from redox.core.common import EV_KJ


def finite(value):
    if type(value) not in (int,float) or not math.isfinite(value):raise ValueError('Finite stated property required')
    return value


def assess(base,first,capacity,spacing,tempo):
    complete=[];missing=[]
    for row in base['candidates']:
        ref=row['reference_deltaG_kJ_mol']
        if ref is None:missing.append(row['id']);continue
        ref=finite(ref);dbu=finite(row['base_models']['dbu']['free_ions_kJ_mol'])
        shift=finite(base['corrections']['dbu']['free_ions_kJ_mol'])
        if abs((dbu-ref)-shift)>1e-6:raise ValueError('Within-base shared term does not reconstruct')
        complete.append(dict(id=row['id'],reference_kJ_mol=ref,DBU_kJ_mol=dbu,
            nonpositive_screen_changes=(ref<=0)!=(dbu<=0)))
    if not complete:raise ValueError('No complete reaction comparisons')
    # A shared term must preserve strict pairwise order at this fixed convention.
    pairs=0
    for i,a in enumerate(complete):
        for b in complete[i+1:]:
            if (a['reference_kJ_mol']<b['reference_kJ_mol'])!=(a['DBU_kJ_mol']<b['DBU_kJ_mol']):
                raise ValueError('A common base term changed pairwise rank')
            pairs+=1
    singles={r['site']:r for r in first['singles']}
    if set(singles)!={2,3}:raise ValueError('Exact site-2/site-3 comparison required')
    penalties=[];preferences={}
    for cutoff in ('25','50','100'):
        values={}
        for site,row in singles.items():
            entry=row['by_cutoff'][cutoff]
            actual=finite(entry['first_graft_kJ_mol']['free_ions'])-finite(entry['capped_graft_kJ_mol']['free_ions'])
            if abs(actual-finite(entry['polymer_minus_capped_kJ_mol']))>1e-7:raise ValueError('Polymer/local difference changed')
            pair=finite(entry['first_graft_kJ_mol']['contact_pair'])-finite(entry['capped_graft_kJ_mol']['contact_pair'])
            if abs(pair-actual)>1e-7:raise ValueError('Shared salt terms fail to cancel')
            penalties.append(dict(site=site,rotor_cutoff_cm=int(cutoff),penalty_kJ_mol=actual))
            values[site]=entry['first_graft_kJ_mol']['free_ions']
        preferences[cutoff]=values[2]-values[3]
    nominal=finite(capacity['reference_capped_two_e_nominal_mAh_g']);bounds={}
    if abs(nominal-2*EV_KJ*1000/(3.6*finite(capacity['reference_capped_mass_g_mol'])))>1e-6:
        raise ValueError('Molecular capacity does not reconstruct from mass and electrons')
    for row in capacity['exact_graphs_at_handle_density_0_6']:
        number=len(row['occupied_sites'])
        if not number:continue
        value=finite(row['two_e_upper_bound_mAh_g'])
        if abs(value-2*number*EV_KJ*1000/(3.6*finite(row['mass_g_mol'])))>1e-6:
            raise ValueError('Fragment capacity does not reconstruct from mass and graft count')
        if number in bounds and abs(bounds[number]-value)>1e-7:raise ValueError('Equal-composition capacity differs')
        bounds[number]=value
    if set(bounds)!={1,2,3}:raise ValueError('Three predeclared occupancy composition bounds required')
    contacts={}
    for row in spacing['states']:
        graph=row['graph'];distance=finite(row['contact']['centroid_distance_A'])
        contacts.setdefault(graph,[]).append(distance)
    if set(contacts)!={'occ_2_3','occ_2_4'} or any(len(v)!=2 for v in contacts.values()):
        raise ValueError('Matched two-pose spacing controls required')
    response=tempo['summary'];left=finite(response['lambda_reduction_eV']);right=finite(response['lambda_oxidation_eV'])
    if min(left,right)<=0:raise ValueError('Positive directional relaxation costs required')
    asym=abs(left-right)/((left+right)/2)
    if abs(asym-finite(response['directional_relative_asymmetry']))>1e-12:raise ValueError('Directional asymmetry does not reconstruct')
    return dict(status='saved_context_consequences_reconstructed',
        route_definition=dict(complete=len(complete),missing=missing,
            changes_at_declared_nonpositive_screen=sum(r['nonpositive_screen_changes'] for r in complete),
            comparisons=complete,unchanged_pairwise_orders=pairs,criterion_scope='Model deltaG <= 0; not an experimental feasibility or synthesis gate'),
        polymer_context=dict(penalties=penalties,positive_penalties=sum(r['penalty_kJ_mol']>0 for r in penalties),
            site2_minus_site3_by_cutoff_kJ_mol=preferences,site_preference_changes_sign=min(preferences.values())<0<max(preferences.values()),
            conformer_convergence_established=False),
        capacity=dict(molecular_two_e_mAh_g=nominal,fragment_two_e_bounds_mAh_g={str(k):v for k,v in sorted(bounds.items())},
            nominal_to_one_graft_bound_ratio=nominal/bounds[1],experimental_loading_inferred=False),
        spatial_contacts=dict(centroid_distances_A={k:sorted(v) for k,v in contacts.items()},
            contact_population_or_binding_G_established=False),
        directional_response=dict(reduction_eV=left,oxidation_eV=right,relative_asymmetry=asym,
            diagnostic_tolerance=response['equal_curvature_diagnostic_tolerance'],
            symmetric_parameter_qualified=asym<=response['equal_curvature_diagnostic_tolerance']),
        new_quantum_calculations=False,algorithmic_superiority_established=False,experimental_material_decisions_validated=False)
