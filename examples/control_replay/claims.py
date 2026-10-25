"""Frozen independent replay functions; provenance is in the manifest."""
import math

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
