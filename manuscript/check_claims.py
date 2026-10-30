"""Independently reconstruct the numerical statements in the frozen paper snapshot."""

import csv
import hashlib
import json
import math
import statistics as descriptive
from collections import Counter
from pathlib import Path


def audit_projected_portability(load, folder):
    if not (folder / "projected_portability.json").exists():
        return None
    data = load("projected_portability.json")
    if (data["status"] != "relocated_projected_partitions_and_physical_guards_independently_checked"
            or {r["id"] for r in data["controls"]} != {"formaldehyde", "methyl_radical"}
            or len(data["controls"]) != 2 or data["negative_controls"] != 10
            or data["original_checkout_file_accesses"] != 0 or data["xtb_CLI_executions"] != 6
            or data["quantum_calculations_run"] or data["native_Hessian_execution_requalified"]
            or data["general_platform_portability_inferred"]):
        raise ValueError("Projected relocation scope changed")
    for row in data["controls"]:
        if (row["maximum_frequency_difference_cm"] != 0
                or row["nuclear_difference_eV_by_cutoff"] != {"25": 0., "50": 0., "100": 0.}):
            raise ValueError("Exact relocated partition agreement changed")
    negative = load("projected_portability_adverse.json")
    if len(negative["checks"]) != 10 or not all(c["rejected"] for c in negative["checks"]):
        raise ValueError("Relocated physical guards failed")
    import xml.etree.ElementTree as ET
    suites = list(ET.parse(folder / "projected_portability_tests.xml").getroot().iter("testsuite"))
    if [sum(int(s.get(k, "0")) for s in suites) for k in ("tests", "failures", "errors", "skipped")] != [22, 0, 0, 0]:
        raise ValueError("Relocation tests incomplete")
    return dict(controls=2, negative_controls=10, original_checkout_file_accesses=0,
                native_GPU_Hessian_execution_requalified=False)


def audit_spacing_qc_and_bibliography(load, folder):
    if not (folder / "spacing_native_QC.json").exists():
        return None
    qc = load("spacing_native_QC.json")
    rows = qc["rows"]
    expected = {c["key"] for graph in load("finite_ensemble_states.json")["cases"] for c in graph["conformers"]}
    if len(rows) != 12 or {r["key"] for r in rows} != expected:
        raise ValueError("Physical-QC table does not cover all twelve exact poses")
    for row in rows:
        values = [row["maximum_gradient_Ha_Bohr"], row["maximum_Hessian_asymmetry"], row["lowest_projected_internal_frequency_cm"], row["largest_sampled_column_relative_error"]]
        if (not all(math.isfinite(v) for v in values) or not 0 <= values[0] <= .00045
                or not 0 <= values[1] <= 1e-6 or values[2] < -5 or not 0 <= values[3] <= .01):
            raise ValueError("Reported spacing physical checks failed")
    gradients = [r["maximum_gradient_Ha_Bohr"] for r in rows]
    if qc["maximum_gradient_range_Ha_Bohr"] != [min(gradients), max(gradients)]:
        raise ValueError("Gradient range does not reconstruct")
    current = (folder / "bibliography_current_audit.json").exists()
    bibliography = load("bibliography_current_audit.json" if current else "bibliography_audit.json")
    checks = bibliography["checks"]
    non_doi = load("bibliography_non_DOI_review.json")["records"]
    total, doi_total = (48, 43) if current else (42, 37)
    if (len(checks) != total or sum(r["status"] == "metadata_matches" for r in checks) != doi_total
            or {r["key"] for r in checks if r["status"] == "non_DOI_source_requires_direct_review"} != {r["key"] for r in non_doi}
            or len(non_doi) != 5):
        raise ValueError("Bibliography identity coverage changed")
    source_studies = None
    if current:
        import unicodedata
        frozen_bib = (folder / "bibliography_current_source.bib").read_text()
        current_bib = Path(__file__).with_name("references.bib").read_text()
        if unicodedata.normalize("NFC", frozen_bib) != unicodedata.normalize("NFC", current_bib):
            raise ValueError("Bibliography changed beyond canonical Unicode normalization")
        inventory = load("benchmark_source_inventory.json")
        endpoint_rows = load("benchmark.json")["outliers"]
        expected_counts = {doi: sum(row["doi"] == doi for row in endpoint_rows) for doi in {row["doi"] for row in endpoint_rows}}
        if ({r["doi"]: r["endpoints"] for r in inventory["studies"]} != expected_counts
                or inventory["total_endpoints"] != 39 or len(inventory["studies"]) != 7):
            raise ValueError("Seven-study benchmark inventory does not reconstruct")
        by_key = {r["key"]: r for r in checks}
        for study in inventory["studies"]:
            row = by_key[study["citation_key"]]
            if row["status"] != "metadata_matches" or row["citation"]["doi"].lower() != study["doi"].lower():
                raise ValueError("Benchmark primary source is not in the checked bibliography")
        source_studies = 7
    return dict(spacing_QC_rows=12, maximum_gradient_range_Ha_Bohr=[min(gradients), max(gradients)],
                bibliography_entries=total, DOI_identity_matches=doi_total, non_DOI_sources_reviewed=5,
                experimental_studies_cited=source_studies)


def audit_neighbor_archive(load):
    try:current=load('neighbor_archive.json')
    except FileNotFoundError:return None
    proof=load('neighbor_archive_relocation.json');manifest=load('neighbor_archive_manifest.json')
    if (proof['status']!='neighbor_archive_inventory_and_two_relocated_replays_checked'
            or current['archive_sha256']!=proof['archive_sha256']
            or current['package_manifest_sha256']!=proof['package_manifest_sha256']
            or proof['archive_regular_files']!=len(manifest['files_sha256'])+1
            or len(manifest['states'])!=14 or len(manifest['reactions'])!=14
            or proof['raw_states']!=14 or proof['balanced_reactions']!=14
            or not proof['nuclear_partition_reintegrated']
            or proof['public_archive_deposited'] or proof['conformer_convergence_established']):
        raise ValueError('Raw neighbor archive coverage, identity or scope differs')
    for key in ('conformer_convergence_established','synthesis_admission',
                'isolated_motif_interaction_computed','electronic_coupling_computed'):
        if manifest[key] or proof['arithmetic'][key]:raise ValueError('Neighbor replay scope was promoted')
    for mode,key in [('stdlib','isolation'),('native','native_isolation')]:
        row=proof[key]
        if (row['attempted_original_reads'] or row['forbidden_modules_loaded'] or not row['isolated_python']
                or (mode=='stdlib' and (not row['site_imports_disabled'] or row['numpy_version'] is not None))
                or (mode=='native' and row['numpy_version']!='2.2.6')):
            raise ValueError('Neighbor replay lacks isolated declared dependencies')
    native=proof['native_reconstruction'];states={r['id']:r for r in manifest['states']}
    expected={r['id'] for r in manifest['states'] if r['partition_kind']!='monatomic'}
    if (len(states)!=14 or native['Hessians_reconstructed']!=13 or native['monatomic_partitions_recomputed']!=1
            or len(native['states'])!=13 or {r['id'] for r in native['states']}!=expected
            or native['xTB_library_loaded'] or native['original_stationarity_spin_and_displaced_SCF_QC_repeated']):
        raise ValueError('Native archive reconstruction coverage differs')
    for row in native['states']:
        errors=[row['frequency_max_error_cm'],*row['nuclear_errors_eV'].values()]
        if (set(row['nuclear_errors_eV'])!={'25','50','100'}
                or not all(math.isfinite(v) and v>=0 for v in errors)
                or errors[0]>1e-5 or max(errors[1:])>1e-8):
            raise ValueError('Neighbor Hessian/nuclear replay exceeds its tolerances')
    if (len(native['monatomic_states'])!=1 or native['monatomic_states'][0]['id']!='chloride'
            or not 0<=native['monatomic_states'][0]['nuclear_error_eV']<=1e-8):
        raise ValueError('Monatomic replay differs')
    source=load('neighbor_cycle.json')
    for alias,name in [('polymer_pair','pair'),('empty','empty'),('site2','site2'),('site3','site3')]:
        value=states[alias]['checked_G_eV']
        if not math.isfinite(value) or abs(value-source['states'][name]['G_eV'])>1e-8:
            raise ValueError('Archived exact chain state differs from paper cycle')
    for cutoff in ('25','50','100'):
        values=proof['arithmetic']['by_cutoff'][cutoff];expected=source['by_cutoff'][cutoff]
        if not math.isfinite(values['deltaDeltaG_kJ_mol']) or abs(values['deltaDeltaG_kJ_mol']-expected['deltaDeltaG_kJ_mol'])>1e-6:
            raise ValueError('Archived neighbor cycle differs')
        for mode,rows in values['models'].items():
            for name,value in rows.items():
                if not math.isfinite(value) or abs(value-expected['models'][mode][name])>1e-6:
                    raise ValueError('Archived conditional reaction differs')
    return dict(raw_states=14,balanced_reactions=14,Hessians_reconstructed=13,
        monatomic_partitions_recomputed=1,original_project_reads_blocked=True,
        independent_nuclear_reconstruction=True,original_quantum_QC_repeated=False)


def audit_neighbor_cycle(load):
    try:
        row=load('neighbor_cycle.json')
    except FileNotFoundError:
        return None
    if (row['status']!='exact_candidate_neighbor_grafting_cycle_independently_reconstructed'
            or row['conformers_per_graph']!=1
            or any(row[k] for k in ('conformer_convergence_established','synthesis_admission',
                                   'isolated_motif_interaction_computed','electronic_coupling_computed',
                                   'graft_loading_distribution_inferred'))
            or set(row['states'])!={'empty','site2','site3','pair'}):
        raise ValueError('Finite matched cycle scope or four-state coverage changed')
    states=row['states'];first=load('pmdi_first_site_grafting.json')
    if row['candidate_id']!=first['candidate_id'] or states['empty']['smiles']!=first['empty_reference']['smiles']:
        raise ValueError('Neighbor cycle and first-site exact identities differ')
    if states['pair']['task']['occupied_sites']!=[2,3]:
        raise ValueError('Neighbor cycle must use the declared adjacent sequence sites')
    for site in (2,3):
        task=states['site'+str(site)]['task']
        if (task['occupied_sites']!=[site] or task['candidate_id']!=states['pair']['task']['candidate_id']
                or task['backbone_stereo_fingerprint']!=states['pair']['task']['backbone_stereo_fingerprint']):
            raise ValueError('Matched occupancy/backbone/candidate identity differs')
    factor=96.48533212331002
    if abs(row['conversion_kJ_mol_per_eV']-factor)>1e-12:
        raise ValueError('Neighbor cycle conversion changed')
    for state in states.values():
        if (not math.isfinite(state['G_eV']) or not math.isfinite(state['E_eV']) or not math.isfinite(state['minimum_frequency_cm'])
                or state['minimum_frequency_cm']<-5
                or abs(state['G_eV']-state['G_by_cutoff_eV']['50'])>1e-8):
            raise ValueError('Neighbor state lacks finite matched primary thermochemistry')
    for cutoff in ('25','50','100'):
        values={name:s['G_by_cutoff_eV'][cutoff] for name,s in states.items()}
        if not all(math.isfinite(v) for v in values.values()):raise ValueError('Nonfinite neighbor state G')
        recorded=row['by_cutoff'][cutoff]
        if any(abs(values[k]-recorded['state_G_eV'][k])>1e-8 for k in values):
            raise ValueError('Cutoff-specific state join differs')
        g0,g2,g3,g23=(values[k] for k in ('empty','site2','site3','pair'))
        delta=math.fsum([g23,g0,-g2,-g3])*factor
        if abs(delta-recorded['deltaDeltaG_kJ_mol'])>1e-7:
            raise ValueError('Neighbor nonadditivity does not reconstruct')
        for name,constant in row['shared_reagent_constants_eV'][cutoff].items():
            expected=dict(graft_2_on_empty_kJ_mol=math.fsum([g2,-g0,constant])*factor,
                graft_3_on_empty_kJ_mol=math.fsum([g3,-g0,constant])*factor,
                graft_2_given_3_kJ_mol=math.fsum([g23,-g3,constant])*factor,
                graft_3_given_2_kJ_mol=math.fsum([g23,-g2,constant])*factor,
                direct_double_graft_kJ_mol=math.fsum([g23,-g0,2*constant])*factor,
                deltaDeltaG_kJ_mol=delta)
            if any(not math.isfinite(recorded['models'][name][k])
                   or abs(v-recorded['models'][name][k])>1e-7 for k,v in expected.items()):
                raise ValueError('Conditional/direct neighboring-graft reaction does not reconstruct')
    if abs(row['deltaDeltaG_kJ_mol']-row['by_cutoff']['50']['deltaDeltaG_kJ_mol'])>1e-8:
        raise ValueError('Primary neighboring-graft value differs')
    electronic=math.fsum([states['pair']['E_eV'],states['empty']['E_eV'],-states['site2']['E_eV'],-states['site3']['E_eV']])*factor
    if not math.isfinite(row['electronic_deltaDeltaE_kJ_mol']) or abs(electronic-row['electronic_deltaDeltaE_kJ_mol'])>1e-7:
        raise ValueError('Electronic neighbor cycle differs')
    return dict(matched_states=4,cutoffs=3,ionic_representations=2,
        deltaDeltaG_kJ_mol=row['deltaDeltaG_kJ_mol'],conformer_convergence_inferred=False,
        isolated_binding_or_coupling_inferred=False)


def audit_control_addendum(load):
    try:
        current=load('control_addendum.json')
    except FileNotFoundError:
        return None
    proof=load('control_archive_relocation.json');full=load('control_full_replay.json')
    if (proof['status']!='archive_inventory_and_two_extracted_replays_checked'
            or (proof['raw_quantum_control_states'],proof['Fc_electronic_records'],proof['raw_SCF_evaluations'],
                proof['baseline_evaluations'],proof['displaced_evaluations'])!=(8,6,392,8,384)
            or full['quantum_states_reconstructed']!=8 or len(full['states'])!=8
            or proof['archive_sha256']!=current['archive_sha256']
            or proof['polymer_G_qualification_inferred'] or full['polymer_G_qualified']
            or full['quantum_calculations_run'] or full['CUDA_or_xTB_loaded']):
        raise ValueError('Control addendum coverage, archive identity or cached-only scope differs')
    for mode,row in proof['isolated_replays'].items():
        if (row['original_project_open_attempts'] or row['forbidden_quantum_modules_loaded']
                or row['validator_sha256']!=proof['validator_sha256']
                or (mode=='full' and row['numpy_version']!='2.2.6')
                or (mode=='electronic' and row['numpy_version'] is not None)):
            raise ValueError('Control replay did not retain isolated declared dependencies')
    if (full['files_verified']!=len(load('control_example_manifest.json')['files_sha256'])
            or proof['isolated_replays']['full']['example_manifest_sha256']!=current['CPU_example_manifest_sha256']):
        raise ValueError('Control replay belongs to another raw example inventory')
    for state in full['states']:
        values=[state['frequency_max_error_cm'],state['maximum_force_error_Ha_Bohr'],
                state['maximum_Hessian_column_relative_error'],*state['nuclear_errors_eV'].values()]
        if (not all(math.isfinite(v) and v>=0 for v in values)
                or state['frequency_max_error_cm']>1e-5 or state['maximum_force_error_Ha_Bohr']>1e-6
                or state['maximum_Hessian_column_relative_error']>.01
                or any(v>1e-8 for v in state['nuclear_errors_eV'].values())):
            raise ValueError('Relocated cached numerical control exceeds its declared tolerance')
    return dict(raw_quantum_records=8,Fc_electronic_records=6,total_SCF_evaluations=392,
                baseline_evaluations=8,displaced_evaluations=384,original_project_reads_blocked=True,
                fresh_quantum_or_polymer_G_inferred=False)


def audit_fc_functional_controls(load):
    """Rebuild electronic-only reference diagnostics from six stored states."""
    try:
        report = load('Fc_functional_review.json')
    except FileNotFoundError:
        return None
    if (report['status'] != 'frozen_geometry_Fc_functional_sensitivity_independently_checked'
            or any(report[k] for k in ('primary_reference_modified', 'primary_benchmark_modified',
                                       'nuclear_thermochemistry_borrowed', 'alternative_potentials_calculated'))):
        raise ValueError('Fc functional control scope was changed')
    methods = {}
    functionals = {'wb97mv': 'wb97m-v', 'b3lyp': 'b3lyp', 'pbe0': 'pbe0'}
    for name, xc in functionals.items():
        states = {}
        for state, charge, mult in [('neu', 0, 1), ('ox', 1, 2)]:
            row = load(f'Fc_functional_{name}_{state}.json')
            if (row['charge'] != charge or row['mult'] != mult or row['sp_xc'] != xc
                    or row['sp_basis'] != 'def2-tzvpd' or row['solvent'] != 'acetonitrile'
                    or not row['converged_smd'] or not row['converged_gas']
                    or any(row[k] for k in ('primary_reference_modified', 'geometry_reoptimized',
                                            'thermochemistry_computed', 'potential_computed'))):
                raise ValueError('Fc electronic control protocol or state identity differs')
            values = [row[k] for k in ('e_smd_eV', 'e_gas_eV', 'smd_spin_contam', 'gas_spin_contam')]
            if not all(math.isfinite(v) for v in values) or max(abs(v) for v in values[2:]) > .1:
                raise ValueError('Fc electronic control contains invalid energy or spin data')
            states[state] = row
        methods[name] = {medium: states['ox'][f'e_{medium}_eV'] - states['neu'][f'e_{medium}_eV']
                         for medium in ('smd', 'gas')}
    errors = [r['error_V'] for r in load('benchmark.json')['outliers']]
    for name, gaps in methods.items():
        row = report['methods'][name]
        shift = gaps['smd'] - methods['wb97mv']['smd']
        residuals = [e - shift for e in errors]
        metrics = dict(bias_V=math.fsum(residuals) / len(residuals),
                       MAE_V=math.fsum(abs(e) for e in residuals) / len(residuals),
                       RMSE_V=math.sqrt(math.fsum(e * e for e in residuals) / len(residuals)))
        expected = {'gas_electronic_gap_eV': gaps['gas'], 'SMD_electronic_gap_eV': gaps['smd'],
                    'differential_solvation_eV': gaps['smd'] - gaps['gas'],
                    'SMD_reference_electronic_shift_eV': shift,
                    'gas_gap_shift_eV': gaps['gas'] - methods['wb97mv']['gas'],
                    'solvation_gap_shift_eV': (gaps['smd'] - gaps['gas'])
                    - (methods['wb97mv']['smd'] - methods['wb97mv']['gas'])}
        if row['thermochemistry_recomputed'] or row['full_alternate_potential_qualified']:
            raise ValueError('An electronic reference control cannot establish an alternative potential')
        if any(not math.isfinite(row[k]) or abs(row[k] - v) > 1e-8 for k, v in expected.items()):
            raise ValueError('Fc functional electronic gaps do not reconstruct')
        if any(not math.isfinite(row['reference_only_counterfactual'][k])
               or abs(row['reference_only_counterfactual'][k] - v) > 1e-9 for k, v in metrics.items()):
            raise ValueError('Fc reference-only counterfactual does not reconstruct')
    return dict(electronic_states=6, functionals=3, gaps_and_counterfactuals_reconstructed=True,
                alternative_potential_qualified=False, primary_reference_modified=False)


def audit_surface_controls(load):
    """Reconstruct component curvatures and verify the stated control scope."""
    try:
        legacy=load('CDS_contact.json')
    except FileNotFoundError:
        return None
    smooth=load('CDS_smooth_component.json')
    qualification=load('CDS_smooth_qualification.json')
    if (qualification['status']!='experimental_neutral_CHNOCl_model_controls_qualified'
            or qualification['primary_candidate_admission'] or qualification['full_molecular_polymer_G_qualified']
            or len(qualification['quantum_controls'])!=4):
        raise ValueError('Surface controls cannot establish polymer G or primary candidate admission')
    method=qualification['method']
    if method['CDS_points_per_sphere']!=2030 or method['electrostatic_surface_points']!=590:
        raise ValueError('Declared alternate surface grids changed')
    series={}
    for label,energies,curvatures in [
        ('legacy',{row['offset_A']:row['E_CDS_Ha'] for row in legacy['points']},legacy['centered_curvature_Ha_Bohr2_by_displacement_A']),
        ('smooth',{float(k):v for k,v in smooth['states']['contact']['contact_scan']['energy_Ha_by_offset_A'].items()},
            smooth['states']['contact']['contact_scan']['centered_curvature_Ha_Bohr2'])]:
        values={}
        for text,reported in curvatures.items():
            h=float(text);value=(energies[h]-2*energies[0.]+energies[-h])/(h/.529177210903)**2
            if not math.isfinite(value) or abs(value-reported)>1e-8:
                raise ValueError('Surface component curvature does not reconstruct')
            values[h]=value
        series[label]=values
    if series['legacy'][1e-5]/series['legacy'][1e-3]<50:
        raise ValueError('Legacy contact-step sensitivity no longer supports the stated diagnostic')
    if (max(series['smooth'].values())-min(series['smooth'].values()))/max(series['smooth'].values())>.01:
        raise ValueError('Smooth component curvature is not stable at the tested steps')
    gradient=load('CDS_smooth_full_gradients.json')
    if gradient['status']!='smooth_amine_imide_full_gradients_checked' or set(gradient['states'])!={'precursor','base','protonated','pair'}:
        raise ValueError('Full imide/amine gradient panel is incomplete')
    for row in gradient['states'].values():
        if row['tested_components']!=row['total_components'] or len(row['checks'])!=3*row['total_components']:
            raise ValueError('Claimed full Cartesian component checks are incomplete')
        if not all(math.isfinite(c[k]) for c in row['checks'] for k in ('analytic_Ha_Bohr','numerical_Ha_Bohr')):
            raise ValueError('Nonfinite surface component-gradient data')
        actual=max(abs(c['analytic_Ha_Bohr']-c['numerical_Ha_Bohr']) for c in row['checks'])
        if abs(actual-row['maximum_error_Ha_Bohr'])>1e-14 or actual>1e-6:
            raise ValueError('Full smooth component-gradient evidence does not reconstruct')
    controls=[]
    for name in ('CDS_smooth_quantum_CHO.json','CDS_smooth_quantum_NCl.json'):
        panel=load(name)
        if panel['status']!='experimental_smooth_quantum_controls_independently_checked' or panel['primary_candidate_admission']:
            raise ValueError('Full quantum surface controls have an invalid status or scope')
        if panel['portable_method_id']!=qualification['portable_method_id']:
            raise ValueError('Different surface models were mixed')
        for row in panel['controls']:
            arm=row['arms']['smooth']
            if not all(math.isfinite(arm[k]) for k in ('maximum_force_error_Ha_Bohr',
                    'maximum_Hessian_column_relative_error','minimum_internal_frequency_cm')):
                raise ValueError('Nonfinite quantum surface-control data')
            if (arm['tested_force_components']!=12 or arm['tested_Hessian_columns']!=2
                    or arm['steps_Bohr']!=[.0005,.001] or arm['maximum_force_error_Ha_Bohr']>1e-6
                    or arm['maximum_Hessian_column_relative_error']>.01 or arm['minimum_internal_frequency_cm']<-5):
                raise ValueError('Quantum surface-control gates failed or were overstated')
            controls.append(arm)
    return dict(component_curvatures_reconstructed=True,neutral_small_quantum_controls=len(controls),
        all_axis_component_states=len(gradient['states']),
        maximum_small_control_force_error_Ha_Bohr=max(v['maximum_force_error_Ha_Bohr'] for v in controls),
        maximum_small_control_Hessian_relative_error=max(v['maximum_Hessian_column_relative_error'] for v in controls),
        polymer_G_qualification_inferred=False,primary_values_modified=False)


def audit_first_site_grafting(load):
    """Reconstruct exact-candidate site and capped/polymer differences offline."""
    try:
        report = load('pmdi_first_site_grafting.json')
    except FileNotFoundError:
        return None
    pair = load('pmdi_pair_electronic.json')
    if (report['status'] != 'available_first_site_grafting_independently_reconstructed'
            or report['solvent'] != 'acetonitrile' or report['temperature_K'] != 298.15
            or report['base'] != 'DBU' or report['synthesis_admission']
            or report['conformer_convergence_established'] or report['paired_G_eV'] is not None
            or report['neighbor_deltaDeltaG_kJ_mol'] is not None):
        raise ValueError('First-site evidence scope or missing neighboring G changed')
    factor = 1.602176634e-19 * 6.02214076e23 / 1000
    if report['candidate_id'] != 'design_pmdi_ethyl' or abs(report['conversion_kJ_mol_per_eV'] - factor) > 1e-12:
        raise ValueError('Exact-candidate identity or energy conversion changed')
    correction = -.00831446261815324 * 298.15 * math.log(.08205736608096 * 298.15)
    if abs(correction - report['contact_pair_standard_state_correction_kJ_mol']) > 1e-10:
        raise ValueError('The contact-pair standard-state correction does not reconstruct')
    sites = {row['site']: row for row in report['singles']}
    if len(report['singles']) != 2 or set(sites) != {2, 3}:
        raise ValueError('Exactly the two declared single-site graphs are required')
    if any(row['task']['candidate_id'] != report['candidate_id'] for row in sites.values()):
        raise ValueError('First-site components belong to another exact candidate')
    results = {}; placement = {}
    for cutoff in ('25', '50', '100'):
        g = {name: row['G_by_cutoff_eV'][cutoff] for name, row in report['reagent_states'].items()}
        constants = dict(free_ions=g['protonated'] + g['chloride'] - g['base'] - g['precursor'],
            contact_pair=g['pair'] - g['base'] - g['precursor'] + correction / factor)
        empty = report['empty_reference']['G_by_cutoff_eV'][cutoff]
        cap = report['local_capped_reference']['G_by_cutoff_eV'][cutoff]
        for name, value in constants.items():
            if abs(value - report['shared_constant_eV_by_cutoff'][cutoff][name]) > 1e-8:
                raise ValueError('The shared first-graft reagent term does not reconstruct')
        for site, row in sites.items():
            product = row['state_G_by_cutoff_eV'][cutoff]
            target = row['by_cutoff'][cutoff]
            relative = (product + g['local'] - empty - cap) * factor
            if abs(relative - target['polymer_minus_capped_kJ_mol']) > 1e-7:
                raise ValueError('Base-canceling local/polymer grafting difference does not reconstruct')
            for name, value in constants.items():
                first = (product - empty + value) * factor
                local = (cap - g['local'] + value) * factor
                if (not math.isfinite(first) or abs(first - target['first_graft_kJ_mol'][name]) > 1e-7
                        or abs(local - target['capped_graft_kJ_mol'][name]) > 1e-7
                        or abs(first - local - relative) > 1e-7):
                    raise ValueError('The first-site/capped reaction sums or closure differ')
            if not all(v['atom_balanced'] and v['charge_balanced'] for v in row['balances'].values()):
                raise ValueError('A first-site reaction lacks checked atom/charge balance')
            if cutoff == '50':
                results[str(site)] = dict(first_graft_kJ_mol=target['first_graft_kJ_mol'],
                    polymer_minus_capped_kJ_mol=relative)
        placement[cutoff] = (sites[2]['state_G_by_cutoff_eV'][cutoff] - sites[3]['state_G_by_cutoff_eV'][cutoff]) * factor
    deltaE = (sites[2]['E_eV'] - sites[3]['E_eV']) * factor
    if (abs(deltaE - report['placement']['site2_minus_site3_deltaE_kJ_mol']) > 1e-7
            or abs(placement['50'] - report['placement']['site2_minus_site3_deltaG_kJ_mol']) > 1e-7):
        raise ValueError('Site-placement energy/free-energy differences do not reconstruct')
    if (pair['status'] != 'state_records_verified' or pair['physical_QC'] != 'passed'
            or pair['task']['occupied_sites'] != [2, 3]
            or pair['task']['candidate_id'] != report['candidate_id']):
        raise ValueError('The paired electronic reference is not qualified within scope')
    cycle = (pair['E_smd_eV'] + report['empty_reference']['E_eV'] - sites[2]['E_eV'] - sites[3]['E_eV']) * factor
    if abs(cycle - report['electronic_neighbor_cycle_kJ_mol']) > 1e-7:
        raise ValueError('The separate electronic neighboring-site cycle does not reconstruct')
    return dict(checked_first_site_reactions=4, sites=results, placement_deltaE_kJ_mol=deltaE,
        placement_deltaG_by_cutoff_kJ_mol=placement, paired_native_G_available=False,
        base_term_cancels_from_local_polymer_difference=True)


def audit_spatial_conformers(load):
    """Recompute frozen centroid distances and within-graph free-energy gaps."""
    try:
        figure = load('spacing_conformers.json')
    except FileNotFoundError:
        return None
    baseline = {row['key']:row for row in load('finite_ensemble.json')['states']}
    expected = {f'occ_{sites}__pose{pose}' for sites in ('2_3','2_4') for pose in (0,1)}
    rows = figure['states']
    if (len(rows) != 4 or {r['key'] for r in rows} != expected
            or figure['conformer_convergence_established'] or figure['electronic_coupling_computed']
            or not figure['proper_rotation_only']):
        raise ValueError('Spatial figure coverage, mapping or physical scope changed')
    distances, gaps = {}, {}
    factor = 1.602176634e-19 * 6.02214076e23 / 1000
    for row in rows:
        if row['graph'] != row['key'].split('__')[0]:
            raise ValueError('Spatial graph identity changed')
        rotation = row['coordinate_rotation']
        if (len(rotation) != 3 or any(len(r) != 3 for r in rotation)
                or any(not math.isfinite(v) for r in rotation for v in r)):
            raise ValueError('Invalid spatial-view rotation')
        determinant = (rotation[0][0]*(rotation[1][1]*rotation[2][2]-rotation[1][2]*rotation[2][1])
            - rotation[0][1]*(rotation[1][0]*rotation[2][2]-rotation[1][2]*rotation[2][0])
            + rotation[0][2]*(rotation[1][0]*rotation[2][1]-rotation[1][1]*rotation[2][0]))
        if (abs(determinant-1) > 1e-10 or any(abs(sum(rotation[i][k]*rotation[j][k] for k in range(3))
                -(1 if i == j else 0)) > 1e-10 for i in range(3) for j in range(3))):
            raise ValueError('The spatial view must be a proper rotation, without reflection or rescaling')
        if abs(row['G_eV']-baseline[row['key']]['G_eV']) > 1e-9:
            raise ValueError('Spatial figure uses another native free-energy parent')
        cores = row['core_coordinates_A']
        if len(cores) != 2 or any(len(c) < 5 for c in cores):
            raise ValueError('Two complete conjugated core coordinate groups required')
        if any(len(p) != 3 or not all(math.isfinite(x) for x in p) for c in cores for p in c):
            raise ValueError('Invalid core coordinates')
        centroids = [[math.fsum(p[j] for p in c)/len(c) for j in range(3)] for c in cores]
        distance = math.dist(*centroids)
        if abs(distance-row['contact']['centroid_distance_A']) > 1e-9:
            raise ValueError('Spatial centroid distance does not reconstruct')
        minimum = min(r['G_eV'] for r in rows if r['graph'] == row['graph'])
        delta = (row['G_eV']-minimum)*factor
        if abs(delta-row['deltaG_conformer_kJ_mol']) > 1e-9:
            raise ValueError('Within-graph conformer free-energy gap does not reconstruct')
        distances[row['key']] = distance
        gaps[row['key']] = delta
    return dict(checked_conformers=4, centroid_distances_A=distances,
        within_graph_deltaG_kJ_mol=gaps, isolated_motif_binding_inferred=False,
        converged_ensemble_inferred=False, electron_transfer_inferred=False)


def audit_quantum_deployment(load):
    """Reconstruct deployment comparisons without using a positive status alone."""
    try:
        report = load('quantum_environment.json')
    except FileNotFoundError:
        return None
    states = load('quantum_environment_states.json')['states']
    changes = load('quantum_environment_compatibility.json')['third_party_changes']
    if (report['status'] != 'fresh_compatible_GPU_deployment_independently_checked'
            or report['original_environment_modified'] or report['method_packages_changed']
            or report['full_portability_claimed'] or len(changes) != 14
            or report['compatible_CUDA_runtime_changes'] != len(changes)
            or report['libxtb_sha256'] != 'c405797f34ed96470f23158406e6d8a8e102681d3c3cbbaea16692f514f22853'
            or any(not r['name'].startswith('nvidia-') or not r['name'].endswith('-cu12')
                or r['recorded'] == r['replacement'] for r in changes)):
        raise ValueError('Deployment profile identity or qualification scope changed')
    expected = {'formaldehyde': 1, 'methyl_radical': 2}
    if len(states) != 2 or {r['case'] for r in states} != set(expected):
        raise ValueError('Both deployment spin controls are required')
    comparisons = {r['case']: r for r in report['controls']}
    for row in states:
        case = row['case']
        old, new = row['original'], row['fresh']
        if (old['mult'] != expected[case] or new['mult'] != expected[case]
                or old['symbols'] != new['symbols']
                or sorted(old['symbols']) != sorted(['C','H','H','O'] if case == 'formaldehyde' else ['C','H','H','H'])
                or any(r['status'] != 'normal_primary_electronic_and_native_admission_passed' for r in (old, new))):
            raise ValueError('Deployment control identity/admission changed')
        for state in (old, new):
            if (len(state['coordinates_A']) != len(state['symbols'])
                    or any(len(x) != 3 or not all(math.isfinite(v) for v in x) for x in state['coordinates_A'])
                    or len(state['thermal']['frequencies_cm']) != 3*len(state['symbols'])-6
                    or not all(math.isfinite(v) for v in state['thermal']['frequencies_cm'])):
                raise ValueError('Deployment geometry or mode inventory is incomplete')
            spin = -8.617333262145e-5 * 298.15 * math.log(state['mult'])
            G = state['energy']['e_smd_eV'] + state['thermal']['nuclear_by_cutoff_eV']['50'] + spin
            if not math.isfinite(G) or abs(G - state['G_eV']) > 1e-9:
                raise ValueError('Deployment free-energy reconstruction failed')
        errors = dict(energy_eV=max(abs(old['energy'][k]-new['energy'][k]) for k in ('e_smd_eV','e_gas_eV')),
            free_energy_eV=abs(old['G_eV']-new['G_eV']),
            frequency_cm=max(abs(a-b) for a,b in zip(sorted(old['thermal']['frequencies_cm']), sorted(new['thermal']['frequencies_cm']))))
        direct_rms = math.sqrt(sum(sum((a-b)**2 for a,b in zip(x,y))
            for x,y in zip(old['coordinates_A'],new['coordinates_A'])) / len(old['symbols']))
        if direct_rms > report['tolerances']['aligned_RMSD_Angstrom']:
            raise ValueError('Deployment direct-coordinate bound failed')
        if any(not math.isfinite(v) or v > report['tolerances'][k] for k,v in errors.items()):
            raise ValueError('Deployment comparison exceeds the declared bounds')
        if any(abs(v-comparisons[case]['differences'][k]) > 1e-9 for k,v in errors.items()):
            raise ValueError('Deployment reported comparison does not reconstruct')
    return dict(controls=2, external_distributions=report['external_distributions'],
        compatible_CUDA_runtime_changes=len(changes), universal_portability_claimed=False)


def audit_reviewer_package(load):
    """Check release coverage against its frozen inventory and replay receipts."""
    try:
        current = load("review_package.json")
    except FileNotFoundError:
        return None  # Earlier immutable epochs predate this packaging report.
    package = load("review_package_manifest.json")
    data = load("review_package_data.json")
    replay = load("review_package_replay.json")
    arithmetic = load("review_package_arithmetic.json")
    hessians = load("review_package_hessians.json")
    states = {r["id"]: r for r in data["states"]}
    expected_hessians = {k for k, v in states.items() if v["partition_kind"] != "monatomic"}
    actual_hessians = {r["id"] for r in hessians["states"]}
    electronic = data["auxiliary"]["electronic"]["records"]
    if (len(states) != len(data["states"]) or len(states) != package["raw_states"]
            or len(states) != current["raw_states"] or len(states) != arithmetic["raw_states"]
            or len(electronic) != package["electronic_only_records"]
            or "package_manifest.json" in package["files_sha256"]
            or len(package["files_sha256"]) + 1 != current["archive_files"]):
        raise ValueError("Reviewer package inventory counts differ")
    if (actual_hessians != expected_hessians or len(actual_hessians) != len(hessians["states"])
            or len(actual_hessians) != hessians["checked"]
            or len(actual_hessians) != current["CPU_Hessian_states"]):
        raise ValueError("Reviewer package Hessian coverage differs")
    if (replay["package_manifest_sha256"] != current["package_manifest_sha256"]
            or replay["arithmetic_interpreter"] != "/usr/bin/python3 -I -S"
            or not replay["original_checkout_and_calculation_opens_blocked_by_audit_hook"]
            or not replay["optional_Hessian_reconstruction"]
            or current["public_archive_deposited"] != package["public_archive_deposited"]
            or current["public_archive_deposited"] != replay["code_and_dataset_publicly_deposited"]
            or current["public_DOI_created"] or replay["full_GPU_environment_restored"]):
        raise ValueError("Reviewer package replay/public-deposition scope differs")
    if (arithmetic["original_cluster_paths_opened"] or arithmetic["quantum_calculations_run"]
            or arithmetic["independent_physical_QC_reestablished"] or arithmetic["experimental_accuracy_inferred"]
            or hessians["original_cluster_paths_opened"] or hessians["quantum_calculations_run"]
            or hessians["xTB_library_loaded"] or hessians["original_physical_QC_repeated"]
            or hessians["experimental_accuracy_inferred"]):
        raise ValueError("Cached replay was promoted to fresh physical or experimental validation")
    coverage = current["primary_coverage"]
    for key in ("benchmark_endpoints", "reference_reactions", "base_assisted_values",
                "relative_basicity_controls", "finite_graph_partitions"):
        if arithmetic[key] != coverage[key]:
            raise ValueError("Reviewer package primary replay coverage differs")
    for key in ("benchmark_endpoints", "reference_reactions", "base_assisted_values"):
        if package[key] != coverage[key]:
            raise ValueError("Reviewer package manuscript inventory differs")
    if arithmetic["auxiliary_scientific_records"] != current["auxiliary_coverage"]:
        raise ValueError("Reviewer package auxiliary replay coverage differs")
    negative = replay["negative_checks"]
    if len(negative) != current["negative_checks"] or len({r["case"] for r in negative}) != len(negative) or not all(r["rejected"] for r in negative):
        raise ValueError("Reviewer package negative checks differ")
    frequency_error = max(r["frequency_max_error_cm"] for r in hessians["states"])
    nuclear_error = max(abs(x) for r in hessians["states"] for x in r["nuclear_errors_eV"].values())
    if (frequency_error != current["maximum_frequency_reconstruction_error_cm"]
            or nuclear_error != current["maximum_nuclear_reconstruction_error_eV"]):
        raise ValueError("Reviewer package reconstruction errors differ")
    try:
        draft = load("review_package_draft.json")
    except FileNotFoundError:
        draft = None
    if draft is not None:
        assets = {r["name"]: r for r in draft["assets"]}
        archive = Path(current["compressed_archive"]).name
        if (draft["draft"] is not True or draft["published_at"] is not None
                or draft["public_raw_archive_deposited"] or draft["DOI_created"]
                or assets[archive]["sha256"] != current["archive_sha256"]
                or assets[archive]["bytes"] != current["compressed_bytes"]
                or assets["package_manifest.json"]["sha256"] != current["package_manifest_sha256"]):
            raise ValueError("Reviewer archive draft status or asset identity differs")
    return dict(raw_states=len(states), electronic_only_records=len(electronic),
                archive_files=len(package["files_sha256"]) + 1, Hessian_states=len(actual_hessians),
                negative_checks=len(negative), public_deposition=False, full_GPU_restoration=False)


def audit_thermal_settings(load):
    """Recompute published cutoff aggregates; retain the primary 50 cm^-1 values."""
    try:
        result = load("thermal_cutoff_sensitivity.json")
    except FileNotFoundError:
        return None
    if (result["quantum_calculations_run"] or result["primary_results_modified"]
            or result["experimental_accuracy_inferred"] or result["candidate_prediction_intervals_established"]
            or result["settings"]["cutoffs_cm"] != [25, 50, 100]
            or result["settings"]["primary_cutoff_cm"] != 50):
        raise ValueError("Cutoff scenarios were promoted beyond thermal-setting sensitivity")
    rows = result["benchmark"]["endpoints"]
    original = {(r["id"], r["event"]): r for r in load("benchmark.json")["outliers"]}
    if len(rows) != len(original) or {(r["id"], r["event"]) for r in rows} != set(original):
        raise ValueError("Cutoff benchmark coverage or identity differs")
    for row in rows:
        source = original[(row["id"], row["event"])]
        if (row["exp_V"] != source["exp_V"] or abs(row["E_by_cutoff_V"]["50"]-source["calc_V"]) > 1e-9
                or row["wave"] != (1 if "neu" in row["event"].split("->") else 2)):
            raise ValueError("Primary benchmark or wave labels changed in cutoff analysis")
        span = 1000*(max(row["E_by_cutoff_V"].values())-min(row["E_by_cutoff_V"].values()))
        if abs(span-row["range_mV"]) > 1e-9:
            raise ValueError("Endpoint cutoff span differs")
    for key in ("25", "50", "100"):
        errors = [r["E_by_cutoff_V"][key]-r["exp_V"] for r in rows]
        computed = dict(MAE_V=math.fsum(abs(v) for v in errors)/len(errors),
            RMSE_V=math.sqrt(math.fsum(v*v for v in errors)/len(errors)), bias_V=math.fsum(errors)/len(errors))
        if any(abs(value-result["benchmark"]["metrics_by_cutoff"][key][name]) > 1e-12 for name,value in computed.items()):
            raise ValueError("Cutoff benchmark statistics do not reconstruct")
    largest = max(r["range_mV"] for r in rows)
    if abs(largest-result["benchmark"]["max_endpoint_range_mV"]) > 1e-9:
        raise ValueError("Maximum cutoff benchmark span differs")
    reactions = result["grafting"]["reactions"]
    if len({r["name"] for r in reactions}) != 87 or len(reactions) != 87:
        raise ValueError("Expected the 67 reference and 20 base-assisted cutoff reactions")
    for row in reactions:
        values = row["SI_deltaG_by_cutoff_kJ_mol"].values()
        if abs(max(values)-min(values)-row["range_kJ_mol"]) > 1e-9:
            raise ValueError("Reaction cutoff span differs")
    crossing = [r["name"] for r in reactions if min(r["SI_deltaG_by_cutoff_kJ_mol"].values()) < 0 < max(r["SI_deltaG_by_cutoff_kJ_mol"].values())]
    max_reaction = max(r["range_kJ_mol"] for r in reactions)
    if crossing != result["grafting"]["zero_crossings"] or abs(max_reaction-result["grafting"]["max_range_kJ_mol"]) > 1e-9:
        raise ValueError("Cutoff reaction sign/range inventory differs")
    review = load("thermal_cutoff_review.json")
    if (review["status"] != "thermal_cutoff_full_state_reconstruction_checked"
            or review["potential_sums"] != 3*len(rows) or review["reaction_sums"] != 3*len(reactions)
            or review["maximum_potential_reconstruction_error_V"] > 1e-8
            or review["maximum_reaction_reconstruction_error_kJ_mol"] > 1e-6):
        raise ValueError("Independent thermal-cutoff source review differs")
    return dict(benchmark_endpoints=len(rows), reaction_models=len(reactions),
                maximum_endpoint_span_mV=largest, maximum_reaction_span_kJ_mol=max_reaction,
                reaction_sign_crossings=len(crossing), prediction_intervals_established=False)


def audit_decision_followup(load):
    """Recompute finite partitions, capacity bounds and decision counts independently."""
    baseline = load("finite_ensemble.json")
    states = {r["key"]: r for r in baseline["states"]}
    if len(states) != 12 or baseline["checked_native_G"] != 12 or baseline["pending_keys"]:
        raise ValueError("Incomplete twelve-state finite baseline")
    conversion = 1.602176634e-19 * 6.02214076e23 / 1000
    rt = 8.31446261815324e-3 * 298.15 / conversion
    partitions = {}
    for pair in baseline["finite_pairs"]:
        result = {}
        for cutoff in ("25", "50", "100"):
            energies = [states[k]["G_by_rotor_cutoff_eV"][cutoff] for k in pair["conformer_keys"]]
            minimum = min(energies)
            g = minimum - rt * math.log(math.fsum(math.exp(-(e-minimum)/rt) for e in energies))
            if abs(g-pair["independent_G_eV_by_rotor_cutoff_cm"][cutoff]) > 1e-8:
                raise ValueError("Finite conformer partition does not reconstruct")
            result[cutoff] = g
        partitions[pair["graph_key"]] = result
    for row in baseline["cycles"]:
        i, j = row["sites"]
        for cutoff in ("25", "50", "100"):
            dg = math.fsum((partitions[f"occ_{i}_{j}"][cutoff], partitions["occ_none"][cutoff],
                           -partitions[f"occ_{i}"][cutoff], -partitions[f"occ_{j}"][cutoff])) * conversion
            if abs(dg-row["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"][cutoff]) > 1e-6:
                raise ValueError("Finite ensemble occupancy cycle does not reconstruct")
        for terms in row["conditional_DBU_grafting"].values():
            for second, first in (("graft_j_given_i_kJ_mol", "graft_j_on_empty_kJ_mol"),
                                  ("graft_i_given_j_kJ_mol", "graft_i_on_empty_kJ_mol")):
                if abs(terms[second]-terms[first]-row["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"]["50"]) > 1e-6:
                    raise ValueError("Finite conditional grafting does not close")
    capacity = load("pmdi_capacity.json")
    if capacity["accessible_capacity_established"] or capacity["experimental_loading_inferred"]:
        raise ValueError("Composition bounds were promoted to material performance")
    weights = dict(C=12.011, H=1.008, N=14.007, O=15.999, Cl=35.453)
    faraday = conversion * 1000
    graph_bounds = {}
    for graph in capacity["exact_graphs_at_handle_density_0_6"]:
        grafts = len(graph["occupied_sites"])
        atoms = dict(C=43+12*grafts, H=45+7*grafts, Cl=3-grafts, N=2*grafts, O=4*grafts)
        formula = "".join(k+(str(v) if v != 1 else "") for k, v in atoms.items() if v)
        mass = math.fsum(weights[k]*v for k, v in atoms.items())
        if graph["formula"] != formula or abs(mass-graph["mass_g_mol"]) > 1e-8:
            raise ValueError("Exact occupancy mass does not reconstruct")
        for n, label in ((1, "one_e"), (2, "two_e")):
            expected = grafts*n*faraday/(3.6*mass)
            if abs(expected-graph[label+"_upper_bound_mAh_g"]) > 1e-8:
                raise ValueError("Capacity upper bound does not reconstruct")
        graph_bounds[tuple(graph["occupied_sites"])] = graph["two_e_upper_bound_mAh_g"]
    combinations = set()
    for row in capacity["composition_grid"]:
        n, h, g = (row[k] for k in ("repeat_units", "reactive_handles", "grafts"))
        if not 0 <= g <= h <= n or (n, h, g) in combinations:
            raise ValueError("Invalid composition grid")
        combinations.add((n, h, g))
        # Independently sum atomic compositions; do not import the producer's mass model.
        atoms = dict(C=8*n+h+12*g, H=8*n+h+2+7*g, Cl=h-g, N=2*g, O=4*g)
        mass = math.fsum(weights[k]*v for k, v in atoms.items())
        for electrons in (1, 2):
            result = row[f"{electrons}e"]
            expected = g*electrons*faraday/(3.6*mass)
            if (abs(result["mass_g_mol"]-mass) > 1e-8
                    or abs(result["capacity_upper_bound_mAh_g"]-expected) > 1e-8):
                raise ValueError("Composition capacity grid does not reconstruct")
    summary = load("all_axes_summary.json")
    decisions = load("all_axes_decisions.json")["candidates"]
    required = [r["id"] for r in load("all_axes_policy.json")["required_axes"]]
    if len({r["product_key"] for r in decisions}) != len(decisions) or len(decisions) != summary["exact_products"]:
        raise ValueError("Decision identities or counts differ")
    if dict(Counter(r["status"] for r in decisions)) != summary["decisions"]:
        raise ValueError("Decision status counts differ")
    eligible = []
    for row in decisions:
        axes = {r["axis"]: r["status"] for r in row["axes"]}
        if set(axes) != set(required) or len(row["axes"]) != len(required):
            raise ValueError("Required decision axes differ")
        passed = all(axes[k] == "pass" for k in required)
        if row["all_required_axes_pass"] != passed or row["eligible_for_Pareto"] != passed:
            raise ValueError("Candidate bypassed all-axis admission")
        if passed:
            eligible.append(row["id"])
    if sorted(eligible) != sorted(summary["eligible_ids"]):
        raise ValueError("Eligible candidate counts differ")
    precursor = load("pmdi_precursor_review.json")
    route = precursor["prospective_route_decision"]
    if (not route["evidence"]["exact_precursor_access_reviewed"] or route["route_status"] != "hold"
            or route["proceed_to_full_redox"] or precursor["exact_polymer_grafting_demonstrated"]):
        raise ValueError("Precursor access was confused with a completed attachment route")
    deformation = load("segment_relaxation.json")
    for pair in deformation["pairs"]:
        if abs(math.fsum(f["deformation_kJ_mol"] for f in pair["fragments"])
               - pair["deformation_kJ_mol_CODATA"]) > 1e-6:
            raise ValueError("Signed fragment deformation does not close")
    origin = None
    try:
        origin = load("basicity_reference_origin.json")
    except FileNotFoundError:
        pass  # Earlier immutable epochs predate this source-origin review.
    if origin is not None:
        if (origin["independent_experimental_pKaH_MeCN"] is not None
                or origin["source_role"]["independent_experimental_origin_eligible"]
                or origin["source"]["origin_class"] != "computed_PCM_estimate"
                or origin["primary_grafting_results_modified"] or origin["candidate_decisions_modified"]):
            raise ValueError("Computed DIPEA source was promoted to experimental validation")
        capture = -.00831446261815324*298.15*math.log(10)*(origin["source"]["reported_value"]-load("acid_base.json")["acid_anchor"]["pKa_MeCN"])
        if abs(capture-origin["computed_literature_capture_kJ_mol"]) > 1e-10:
            raise ValueError("DIPEA model sensitivity capture does not reconstruct")
        for row in origin["candidate_model_sensitivities"]:
            if abs(row["reference_deltaG_kJ_mol"]+capture-row["fixed_reference_plus_computed_literature_capture_kJ_mol"]) > 1e-9 or row["calibrated_full_route_prediction"]:
                raise ValueError("DIPEA sensitivity changed or acquired a calibration claim")
    return dict(finite_native_G=len(states), finite_graph_partitions=len(partitions),
                finite_cycles=len(baseline["cycles"]), capacity_graphs=len(graph_bounds),
                capacity_compositions=len(combinations),
                all_axis_graphs=len(decisions), all_axis_eligible=len(eligible),
                all_axis_status_counts=summary["decisions"], signed_deformation_pairs=len(deformation["pairs"]))


def audit_tempo_exact_support(load):
    try:
        route=load('tempo_route.json')
    except FileNotFoundError:
        return None
    screen=load('tempo_route_screen.json');native=load('tempo_native_potential.json')
    decision=load('tempo_decision.json');summary=load('tempo_summary.json')
    if (route['candidate']['id']!='tempo' or screen['candidates'][0]['route_status']!='proceed'
            or screen['candidates'][0]['evidence']!=route['evidence']
            or not all(route['NaH_overall'][key] for key in ('atom_balanced','charge_balanced','radical_electrons_balanced'))
            or route['NaH_overall']['radical_electrons']!=1
            or native['product_smiles']!=route['candidate']['product_smiles']
            or native['experimental_exact_graft_potential_claimed']
            or native['free_energy_of_grafting_computed'] or native['all_axis_eligibility']):
        raise ValueError('TEMPO route, identity or scientific scope differs')
    states=native['states'];expected={'tempo__rad':(0,2),'tempo__ox':(1,1),'ferrocene__neu':(0,1),'ferrocene__ox':(1,2)}
    if set(states)!=set(expected):raise ValueError('TEMPO/Fc state coverage differs')
    for key,(charge,mult) in expected.items():
        s=states[key]
        if (s['charge'],s['multiplicity'])!=(charge,mult) or abs(s['G_eV']-s['active_E_smd_eV']-s['native_thermal_eV'])>1e-8:
            raise ValueError('TEMPO state free-energy reconstruction differs')
    values={}
    for cutoff in ('25','50','100'):
        values[cutoff]=math.fsum((states['tempo__ox']['G_by_cutoff_eV'][cutoff],
            -states['tempo__rad']['G_by_cutoff_eV'][cutoff],-states['ferrocene__ox']['G_by_cutoff_eV'][cutoff],
            states['ferrocene__neu']['G_by_cutoff_eV'][cutoff]))
        if abs(values[cutoff]-native['potential_V_vs_Fc_by_rotor_cutoff'][cutoff])>1e-9:
            raise ValueError('TEMPO/Fc potential arithmetic differs')
    result=decision['decision']
    if (result['status']!='pending' or result['all_required_axes_pass'] or len(result['pending_axes'])!=9
            or next(a for a in result['axes'] if a['axis']=='attachment_route')['status']!='pass'
            or summary['actual_route_deltaG_computed'] or summary['exact_graft_experimental_potential_known']):
        raise ValueError('TEMPO all-axis decision or uncertainty scope differs')
    return dict(route='proceed',potential_V_vs_Fc_by_cutoff=values,material_decision='pending',
                pending_axes=len(result['pending_axes']),new_quantum_jobs=0)


def audit_tempo_archive(load,digest):
    try:
        archive=load('tempo_archive.json')
    except FileNotFoundError:
        return None
    relocation=load('tempo_archive_relocation.json');package=load('tempo_archive_manifest.json')
    if (archive['status']!='TEMPO_raw_archive_and_two_relocated_replays_checked'
            or relocation['status']!='TEMPO_archive_and_two_relocated_replays_checked'
            or archive['raw_states']!=4 or archive['Hessians_reconstructed']!=4
            or relocation['archive_sha256']!=archive['archive_sha256']
            or digest('tempo_archive_manifest.json')!=archive['package_manifest_sha256']
            or relocation['package_manifest_sha256']!=archive['package_manifest_sha256']
            or len(package['states'])!=4 or archive['public_archive_deposited'] or archive['public_DOI_created']):
        raise ValueError('TEMPO raw archive coverage or scope differs')
    for key in ('stdlib_isolation','native_isolation'):
        isolation=relocation[key]
        if (not isolation['isolated_python'] or isolation['original_project_open_attempts']
                or isolation['forbidden_modules_loaded']):
            raise ValueError('TEMPO relocated replay isolation differs')
    native=relocation['native_result']
    if (not native['nuclear_partition_reintegrated'] or native['native_reconstruction']['Hessians_reconstructed']!=4
            or native['experimental_accuracy_inferred'] or native['grafting_deltaG_computed']
            or native['material_eligibility_inferred']):
        raise ValueError('TEMPO Hessian replay or interpretation differs')
    expected=load('tempo_native_potential.json')['potential_V_vs_Fc_by_rotor_cutoff']
    if any(abs(native['oxidation_V_vs_Fc_by_cutoff'][c]-expected[c])>1e-9 for c in ('25','50','100')):
        raise ValueError('TEMPO archived potential differs from manuscript evidence')
    return dict(raw_states=4,Hessians_reconstructed=4,original_project_reads_blocked=True,
                quantum_calculations_repeated=False,public_archive_deposited=False)


def audit_local_motif(load):
    try:
        panel = load('local_motif_panel.json')
    except FileNotFoundError:
        return None
    if (panel['status'] != 'all_twelve_frozen_local_interaction_triplets_independently_reconstructed'
            or panel['completed_triplets'] != 12 or panel['declared_triplets'] != 12
            or panel['checked_electronic_geometries'] != 36 or len(panel['states']) != 12
            or any(panel[key] for key in ('free_energies_computed', 'electronic_coupling_computed',
                                          'conformer_convergence_established', 'basis_superposition_correction_computed'))):
        raise ValueError('Local-motif panel scope differs')
    rows = {}
    for row in panel['states']:
        key = (row['parent_key'], round(row['cap_length_A'], 2))
        if key in rows or row['missing_states']:
            raise ValueError('Duplicate or incomplete local-motif triplet')
        if row['free_energy_computed'] or row['basis_superposition_correction_computed']:
            raise ValueError('Local-motif result was promoted beyond electronic energy')
        for medium in ('gas', 'acetonitrile'):
            e = row['exact_state_energies_eV']
            actual = (e['AB'][medium]-e['A'][medium]-e['B'][medium])*96.48533212331002
            if not math.isfinite(actual) or abs(actual-row['interaction_energy_kJ_mol_by_medium'][medium]) > 1e-7:
                raise ValueError('Local-motif interaction does not reconstruct')
        rows[key] = row
    parents = ('occ_2_3__pose0', 'occ_2_3__pose1', 'occ_2_4__pose0', 'occ_2_4__pose1')
    if set(rows) != {(p, c) for p in parents for c in (1.07, 1.09, 1.11)}:
        raise ValueError('Local-motif coverage differs')
    central = {p: rows[p, 1.09]['interaction_energy_kJ_mol_by_medium']['acetonitrile'] for p in parents}
    for parent, expected in zip(parents, (-24.7615, 0.0477, 0.0613, -26.2196)):
        if abs(central[parent]-expected) > 0.0001:
            raise ValueError('Local-motif manuscript result differs')
    expected_ranges = ((-24.81, -24.71), (0.047, 0.049), (0.060, 0.062), (-26.24, -26.20))
    for parent, (expected_min, expected_max) in zip(parents, expected_ranges):
        values = [rows[parent, c]['interaction_energy_kJ_mol_by_medium']['acetonitrile']
                  for c in (1.07, 1.09, 1.11)]
        tolerance = 0.0051 if expected_min < 0 else 0.00051
        if abs(min(values)-expected_min) > tolerance or abs(max(values)-expected_max) > tolerance:
            raise ValueError('Local-motif cap sensitivity differs from manuscript')
    return dict(triplets=12, electronic_geometries=36, central_mECN_kJ_mol=central,
                free_energy_computed=False, electronic_coupling_computed=False,
                basis_superposition_correction_computed=False)


def audit():
    root = Path(__file__).resolve().parents[1]
    pointer = json.loads((root / "manuscript/evidence/current.json").read_text())
    folder = (root / pointer["manifest"]).parent
    load = lambda name: json.loads((folder / name).read_text())
    replay = None
    if (folder/'polymer_native_replay.json').exists():
        from check_polymer_replay import audit as audit_polymer_replay
        replay = audit_polymer_replay(folder)
    benchmark = load("benchmark.json")
    panel = load("benchmark_panel.json")
    rows = list(csv.DictReader((folder / "benchmark_points.csv").open()))
    errors = [float(r["calc_V"]) - float(r["exp_V"]) for r in rows]
    keys = [(r["id"], r["event"], r["doi"]) for r in rows]
    targets = {(r["id"], r["event"], r["doi"]): r for r in panel["targets"]}
    if len(set(keys)) != len(keys) or not set(keys) <= set(targets):
        raise ValueError("Duplicated or non-panel manuscript benchmark point")
    if any(abs(float(r["exp_V"]) - targets[key]["exp_V"]) > 1e-12 for r, key in zip(rows, keys)):
        raise ValueError("Manuscript experimental label differs from frozen panel")
    statistics = dict(MAE_V=math.fsum(abs(e) for e in errors) / len(errors),
                      RMSE_V=math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
                      bias_V=math.fsum(errors) / len(errors), maximum_absolute_error_V=max(abs(e) for e in errors))
    if any(abs(statistics[key] - benchmark["metrics"][key]) > 1e-12 for key in statistics):
        raise ValueError("Manuscript benchmark arithmetic mismatch")
    if (len(rows) != benchmark["admitted_events"] or len(targets) != benchmark["target_events"]
            or len({r["id"] for r in rows}) != benchmark["admitted_molecules"]
            or len({r["doi"] for r in rows}) != benchmark["original_studies"]):
        raise ValueError("Manuscript coverage mismatch")
    bias = load("reference_bias.json")
    if bias["primary_potentials_modified"] or bias["external_validation_claimed"]:
        raise ValueError("Diagnostic offsets replaced the primary benchmark")
    cv_checks = 0
    for model in bias["cross_validation"]:
        residuals = []
        for index, row in enumerate(rows):
            label = model["holdout"]
            members = [j for j, other in enumerate(rows)
                       if (j != index if label == "endpoint" else other[label] != row[label])
                       and (model["model"] == "constant_offset"
                            or targets[keys[j]]["wave"] == targets[keys[index]]["wave"])]
            if not members:
                continue
            offset = math.fsum(errors[j] for j in members) / len(members)
            residuals.append(errors[index] - offset)
        checked = dict(MAE_V=math.fsum(abs(e) for e in residuals)/len(residuals),
                       RMSE_V=math.sqrt(math.fsum(e*e for e in residuals)/len(residuals)),
                       bias_V=math.fsum(residuals)/len(residuals))
        if any(abs(checked[k]-model["metrics"][k]) > 1e-12 for k in checked):
            raise ValueError("Grouped offset diagnostic does not reconstruct")
        cv_checks += 1
    for interval in bias["uncertainty"]:
        primary = next(r for r in benchmark["uncertainty"] if r["cluster"] == interval["cluster"])
        if interval["intervals"] != primary["intervals"]:
            raise ValueError("Printed bootstrap interval differs from primary evidence")
    base = load("base_assisted.json")
    base_values = 0
    for candidate in base["candidates"]:
        for name, modes in candidate["base_models"].items():
            for mode, value in modes.items():
                if value is None:
                    continue
                if abs(value - candidate["reference_deltaG_kJ_mol"] - base["corrections"][name][mode]) > 1e-7:
                    raise ValueError("Manuscript base-assisted cycle mismatch")
                base_values += 1
    dipea_experiment = None
    if (folder / "dipea_experimental_anchor.json").exists():
        anchor = load("dipea_experimental_anchor.json")
        c = anchor["comparison"]
        if (anchor["status"] != "primary_DIPEA_salt_NMR_estimate_source_checked_conditional_comparison"
                or anchor["source"]["doi"] != "10.1039/d3sc02611f"
                or anchor["source"]["counterion"] != "hexafluorophosphate"
                or anchor["source"]["reported_pKa"] != 18.7
                or anchor["calibrated_full_route_prediction"]
                or anchor["primary_grafting_results_modified"] or anchor["candidate_decisions_modified"]):
            raise ValueError("DIPEA experimental-source identity or scope changed")
        acid = load("acid_base.json")
        old = load("basicity_reference_origin.json")
        factor = 8.31446261815324*acid["temperature_K"]*math.log(10)/1000
        nmr = anchor["source"]["reported_pKa"]
        hcl = acid["acid_anchor"]["pKa_MeCN"]
        dbu = next(r["reference"]["pKaH_MeCN"] for r in acid["rows"] if r["base"] == "dbu")
        dft = base["corrections"]["dipea"]["free_ions_kJ_mol"]
        reference = -factor*(nmr-hcl)
        relative = factor*(dbu-nmr)
        expected = dict(RTln10_kJ_mol=factor,pKa_implied_capture_kJ_mol=reference,
            DFT_separated_ion_capture_kJ_mol=dft,DFT_minus_pKa_implied_capture_kJ_mol=dft-reference,
            pKa_implied_DIPEA_DBU_transfer_kJ_mol=relative,
            DFT_minus_pKa_implied_relative_kJ_mol=old["DFT_relative_DIPEA_DBU_transfer_kJ_mol"]-relative)
        if any(abs(c[k]-value)>1e-9 for k,value in expected.items()):
            raise ValueError("DIPEA experimental comparison does not reconstruct")
        frozen = {r["id"]:r for r in base["candidates"]}
        for row in anchor["candidate_hybrid_sensitivities"]:
            original=frozen[row["id"]]
            if row["hybrid_is_calibrated_reaction_prediction"]:raise ValueError("Hybrid sensitivity was promoted to a prediction")
            if original["reference_deltaG_kJ_mol"] is None:
                if row["hybrid_capture_sensitivity_kJ_mol"] is not None:raise ValueError("Missing grafting value was filled")
            elif abs(row["hybrid_capture_sensitivity_kJ_mol"]-original["reference_deltaG_kJ_mol"]-reference)>1e-9:
                raise ValueError("DIPEA hybrid sensitivity does not reconstruct")
        dipea_experiment=dict(source_pKa=nmr,capture_difference_kJ_mol=dft-reference,
            relative_difference_kJ_mol=old["DFT_relative_DIPEA_DBU_transfer_kJ_mol"]-relative,
            exact_chloride_experimental_match=False,hybrid_is_calibrated_prediction=False)
    dipea_uv = None
    if (folder / "dipea_uv_anchor.json").exists():
        uv = load("dipea_uv_anchor.json"); reviewed = load("dipea_uv_review.json")
        if (uv["status"] != "primary_DIPEA_UV_protonation_constant_source_checked_conditional_comparison"
                or uv["source"]["doi"] != "10.1039/d1ra04738h"
                or uv["source"]["reported_value"] != 18.64
                or uv["source"]["reported_notation"] != "18.64(1)"
                or uv["source"]["temperature_C"] != 25.0
                or reviewed["status"] != "independent_UV_source_and_comparison_review_checked"
                or reviewed["UV_reported_logKH"] != "18.64(1)"
                or uv["calibrated_full_route_prediction"] or uv["candidate_decisions_modified"]):
            raise ValueError("DIPEA UV source identity or scope changed")
        sources=[h for path,h in reviewed["source_sha256"].items() if Path(path).name == "summary.json"]
        if len(sources)!=1 or hashlib.sha256((folder/"dipea_uv_anchor.json").read_bytes()).hexdigest()!=sources[0]:
            raise ValueError("Frozen UV assessment differs from its independent source review")
        acid = load("acid_base.json"); old = load("basicity_reference_origin.json");c=uv["comparison"]
        factor=8.31446261815324*acid["temperature_K"]*math.log(10)/1000
        hcl=acid["acid_anchor"]["pKa_MeCN"]
        dbu=next(r["reference"]["pKaH_MeCN"] for r in acid["rows"] if r["base"]=="dbu")
        value=uv["source"]["reported_value"]
        capture=-factor*(value-hcl);relative=factor*(dbu-value)
        dft=base["corrections"]["dipea"]["free_ions_kJ_mol"]
        expected=dict(pKa_implied_capture_kJ_mol=capture,DFT_capture_kJ_mol=dft,
            DFT_minus_pKa_implied_capture_kJ_mol=dft-capture,
            pKa_implied_relative_base_transfer_kJ_mol=relative,
            DFT_minus_pKa_implied_relative_kJ_mol=old["DFT_relative_DIPEA_DBU_transfer_kJ_mol"]-relative)
        if any(abs(c[k]-v)>1e-9 for k,v in expected.items()):raise ValueError("DIPEA UV comparison does not reconstruct")
        frozen={r["id"]:r for r in base["candidates"]}
        for row in uv["candidate_hybrid_sensitivities"]:
            original=frozen[row["id"]]
            if row["calibrated_reaction_prediction"]:raise ValueError("UV hybrid sensitivity was promoted to a prediction")
            if original["reference_deltaG_kJ_mol"] is None:
                if row["hybrid_capture_sensitivity_kJ_mol"] is not None:raise ValueError("Missing UV sensitivity was filled")
            elif abs(row["hybrid_capture_sensitivity_kJ_mol"]-original["reference_deltaG_kJ_mol"]-capture)>1e-9:
                raise ValueError("DIPEA UV hybrid sensitivity does not reconstruct")
        dipea_uv=dict(source_logKH=value,capture_difference_kJ_mol=dft-capture,
            relative_difference_kJ_mol=old["DFT_relative_DIPEA_DBU_transfer_kJ_mol"]-relative,
            matching_activity_conventions_established=False,hybrid_is_calibrated_prediction=False)
    polymer = load("polymer_validation.json")
    cases = {r["key"]: r for r in polymer["cases"]}
    conversion_differences = []
    for cycle in polymer["cycles"]:
        difference_eV = math.fsum(cases[key]["E_eV"] * weight for key, weight in cycle["coefficients"].items())
        # Existing source-frozen records use this documented five-digit factor.
        # Also record its difference from the SI conversion rather than hiding it.
        value = difference_eV * 96.485
        conversion_differences.append(abs(difference_eV * 96.48533212 - value))
        if abs(value - cycle["electronic_deltaDeltaE_kJ_mol"]) > 1e-6:
            raise ValueError("Manuscript electronic occupancy cycle mismatch")
        if any(cases[key]["G_eV"] is None for key in cycle["coefficients"]):
            if cycle["free_energy_deltaDeltaG_kJ_mol"] is not None:
                raise ValueError("Missing native thermochemistry was filled in")
        else:
            dg = math.fsum(cases[key]["G_eV"] * weight for key, weight in cycle["coefficients"].items()) * 96.485
            if abs(dg - cycle["free_energy_deltaDeltaG_kJ_mol"]) > 1e-6:
                raise ValueError("Manuscript native occupancy cycle mismatch")
            for cutoff, expected in cycle["rotor_cutoff_sensitivity_kJ_mol"].items():
                value = math.fsum(cases[key]["cutoff_G_eV"][cutoff] * weight for key, weight in cycle["coefficients"].items()) * 96.485
                if abs(value - expected) > 1e-6:
                    raise ValueError("Manuscript rotor-cutoff sensitivity mismatch")
            for terms in cycle["conditional_grafting_kJ_mol"].values():
                for second, first in (("second_j_given_i", "first_j"), ("second_i_given_j", "first_i")):
                    if abs(terms[second] - terms[first] - dg) > 1e-6:
                        raise ValueError("Conditional grafting difference does not close")
    routes = load("exact_routes.json")
    if routes["completed"] != routes["expected"] or sum(routes["counts"].values()) != routes["completed"]:
        raise ValueError("Exact route handoff coverage mismatch")
    acid = load("acid_base.json")
    for row in acid["rows"]:
        if row["reference"] is None:
            if row["external_capture_kJ_mol"] is not None:
                raise ValueError("Missing exact-base anchor filled in")
            continue
        expected = -.00831446261815324 * acid["temperature_K"] * math.log(10) * (row["reference"]["pKaH_MeCN"] - acid["acid_anchor"]["pKa_MeCN"])
        if abs(expected - row["external_capture_kJ_mol"]) > 1e-10:
            raise ValueError("External acid-capture conversion mismatch")
    concurrency = load("cpu_concurrency.json")
    if concurrency["run_count"] != len(concurrency["runs"]) or concurrency["serial_discrepancies"]:
        raise ValueError("Concurrency experiment count or independent comparison differs")
    for condition in concurrency["conditions"]:
        runs = [r for r in concurrency["runs"] if r["condition"] == condition["condition"]]
        if len(runs) != condition["repeats"]:
            raise ValueError("Concurrency repeat count differs")
        for key in ("launch_to_observed_completion_seconds", "workers_executing_chunks"):
            values = [r[key] for r in runs]
            expected = dict(median=descriptive.median(values), minimum=min(values), maximum=max(values))
            if expected != condition[key]:
                raise ValueError("Concurrency summary does not reconstruct from actual runs")
    fresh = load("fresh_checkout.json")
    if not fresh["figure_CSV_exact_match"] or not fresh["figure_RGB_pixels_exact_match"] or fresh["serial_discrepancies"]:
        raise ValueError("Fresh-checkout example did not reproduce the checked results")
    transfer = load("proton_transfer.json")
    components = {r["id"]: r["G_eV"] for r in transfer["states"]}
    references = load("proton_transfer_sources.json")
    for row, anchor in zip(transfer["comparisons"], references["comparisons"]):
        if row["base"] != anchor["base"]:
            raise ValueError("Relative-basicity source identity mismatch")
        key = row["base"]
        dg = (components[key+"_protonated"] + components["dbu_neutral"]
              - components[key+"_neutral"] - components["dbu_protonated"]) * 96.48533212331002
        external = .00831446261815324 * 298.15 * math.log(10) * (
            anchor["reference"]["pKaH_MeCN"]-anchor["source"]["pKaH_MeCN"])
        if abs(dg-row["DFT_deltaG_kJ_mol"]) > 1e-6 or abs(external-row["external_deltaG_kJ_mol"]) > 1e-9:
            raise ValueError("Relative-basicity comparison does not reconstruct")
    restart = load("quantum_restart.json")
    for control in restart["controls"]:
        if any(control["differences"][key] > tolerance for key, tolerance in restart["tolerances"].items()):
            raise ValueError("Quantum restart exceeds its predeclared tolerance")
    basis = load("neutral_basis.json")
    if basis["checked"] != 5 or basis["pending"] or basis["primary_values_modified"]:
        raise ValueError("Neutral basis panel is incomplete or modified primary values")
    for row in basis["rows"]:
        if abs(row["delta_thermal_eV"] + row["delta_geometry_SP_eV"]
               - row["delta_first_reduction_V_holding_anion_and_Fc_fixed"]) > 1e-8:
            raise ValueError("Basis sensitivity components do not close")
    exact = load("exact_case_join.json")
    for row in exact["cases"]:
        g = row["state_G_eV"]
        if (abs(g["neu"] - g["red1"] - exact["Fc_absolute_V"] - row["E1_V_vs_Fc"]) > 1e-9
                or abs(g["red1"] - g["red2"] - exact["Fc_absolute_V"] - row["E2_V_vs_Fc"]) > 1e-9
                or row["route"]["proceed_to_full_redox"] or row["end_to_end_complete"]):
            raise ValueError("Exact-case values or pending route status changed")
    cached = load("cached_states_manifest.json")
    states = {row["id"]: row["checked_G_eV"] for row in cached["states"]}
    reproduced = load("cached_reproduction.json")
    p = cached["potential"]
    potential = states[p["oxidized"]] - states[p["reduced"]] - states[p["reference_oxidized"]] + states[p["reference_reduced"]]
    if abs(potential - reproduced["PMDI_first_reduction_V_vs_Fc"]) > 1e-9:
        raise ValueError("Cached example potential changed")
    for recipe, value in zip(cached["reactions"], reproduced["reactions"]):
        dg = math.fsum(states[key] * nu for key, nu in recipe["coefficients"].items())
        if abs(dg - value["deltaG_eV"]) > 1e-9:
            raise ValueError("Cached example grafting sum changed")
    added = audit_decision_followup(load) if (folder / "finite_ensemble.json").exists() else None
    release = audit_reviewer_package(load)
    cutoff = audit_thermal_settings(load)
    deployment = audit_quantum_deployment(load)
    spatial = audit_spatial_conformers(load)
    first_site = audit_first_site_grafting(load)
    surface = audit_surface_controls(load)
    fc_functional = audit_fc_functional_controls(load)
    neighbor = audit_neighbor_cycle(load)
    controls_addendum = audit_control_addendum(load)
    neighbor_archive = audit_neighbor_archive(load)
    if neighbor_archive is not None:
        if hashlib.sha256((folder/'neighbor_archive_manifest.json').read_bytes()).hexdigest()!=load('neighbor_archive.json')['package_manifest_sha256']:
            raise ValueError('Frozen neighbor manifest differs from the extracted package')
        import xml.etree.ElementTree as ET
        suites=list(ET.parse(folder/'neighbor_archive_tests.xml').getroot().iter('testsuite'))
        totals=[sum(int(s.get(k,'0')) for s in suites)
                for k in ('tests','failures','errors','skipped')]
        if totals!=[14,0,0,0]:raise ValueError('Neighbor archive integration checks did not all pass')
    if cutoff is not None:
        if hashlib.sha256((folder / "thermal_cutoff_sensitivity.json").read_bytes()).hexdigest() != load("thermal_cutoff_review.json")["result_sha256"]:
            raise ValueError("Frozen thermal-cutoff result differs from its independent source review")
    if release is not None:
        package = load("review_package_manifest.json")
        current = load("review_package.json")
        digest = lambda name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if (digest("review_package_manifest.json") != current["package_manifest_sha256"]
                or digest("review_package_data.json") != package["files_sha256"]["data/manifest.json"]
                or package["code_commit"] != current["code_snapshot_commit"]):
            raise ValueError("Frozen reviewer inventory differs from the tested package")
    return dict(status="numerical_snapshot_reconstructed", polymer_native_replay=replay, benchmark=dict(admitted=len(rows), target=len(targets),
                molecules=len({r["id"] for r in rows}), studies=len({r["doi"] for r in rows}), **statistics),
                base_assisted_values=base_values, dipea_experimental_anchor=dipea_experiment, dipea_uv_anchor=dipea_uv,
                qualified_spacing_G=sum(r["G_eV"] is not None for r in cases.values()),
                spacing_graphs=len(cases), exact_route_entries=len(load("route_evidence.json")["candidates"]),
                exact_routes_checked=routes["completed"], route_decisions=routes["counts"],
                acid_base_scope=acid["status"], conditional_native_cycles=len(polymer["cycles"]),
                checked_concurrency_runs=len(concurrency["runs"]), fresh_checkout_designs=fresh["exact_designs"],
                source_energy_conversion_kJ_mol_per_eV=96.485,
                maximum_SI_conversion_difference_kJ_mol=max(conversion_differences),
                checked_offset_models=cv_checks, relative_basicity_controls=len(transfer["comparisons"]),
                quantum_restart_controls=len(restart["controls"]),
                neutral_basis_controls=basis["checked"], exact_cases_with_explicit_holds=len(exact["cases"]),
                decision_followup=added,
                reviewer_package=release,
                thermal_cutoff_sensitivity=cutoff,
                quantum_environment=deployment,
                spatial_conformers=spatial,
                first_site_grafting=first_site,
                surface_controls=surface,
                Fc_functional_controls=fc_functional,
                neighbor_cycle=neighbor,
                control_addendum=controls_addendum,
                neighbor_archive=neighbor_archive,
                tempo_exact_support=audit_tempo_exact_support(load),
                tempo_archive=audit_tempo_archive(load,lambda name:hashlib.sha256((folder/name).read_bytes()).hexdigest()),
                local_motif=audit_local_motif(load),
                projected_partition_relocation=audit_projected_portability(load, folder),
                spacing_QC_and_bibliography=audit_spacing_qc_and_bibliography(load, folder),
                experimental_accuracy_inferred=False, submission_ready=False)


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
