"""Reconstruct the exact TEMPO directional result from six frozen SCF cells."""
import hashlib
import json
import math
from pathlib import Path


def close(actual,expected,tolerance=1e-8):
    if (any(type(v) not in (float,int) or not math.isfinite(v) for v in (actual,expected))
            or abs(actual-expected)>tolerance):
        raise ValueError('Frozen TEMPO solution quantity does not reconstruct')


def audit(folder):
    load=lambda name:json.loads((folder/name).read_text())
    data=load('tempo_solution_numeric_inputs.json');review=load('tempo_solution_review.json')
    candidate=load('tempo_solution_candidate.json');figure=load('tempo_solution_figure.json')
    if (review['status']!='exact_TEMPO_directional_solution_descriptor_independently_checked'
            or not review['summary']['calculation_verified'] or not review['summary']['endpoint_model_verified']
            or review['material_admission'] or review['solution_rate_validated']
            or review['global_conformer_convergence_proven'] or len(data['cells'])!=6):
        raise ValueError('TEMPO calculation scope or coverage differs')
    if candidate['product_inchikey']!=review['product_inchikey'] or candidate['all_axis_decision']['status']!='pending':
        raise ValueError('Exact candidate identity or all-axis decision differs')
    cells={}
    geometries={Path(path).parent.name:digest for path,digest in data['input_sha256'].items()}
    if set(geometries)!={'ox','rad'} or data['method']['solvent_model']!='native_SMD_frozen_slow':
        raise ValueError('Exact geometry pair or solvent-response model differs')
    for cell in data['cells']:
        tag=(cell['state'],cell['at'],cell['reference'])
        if tag in cells or not cell['converged'] or cell['state'] not in ('ox','rad'):
            raise ValueError('Duplicated or unqualified energy cell')
        charge,mult=(1,1) if cell['state']=='ox' else (0,2)
        if (cell['charge'],cell['mult'])!=(charge,mult):raise ValueError('State charge/spin differs')
        if cell['at'] not in geometries or cell['geometry_sha256']!=geometries[cell['at']]:
            raise ValueError('SCF geometry does not match its declared pair role')
        close(cell['e_eV'],cell['e_Ha']*27.211386245988)
        close(cell['e_Ha'],min(g['e_Ha'] for g in cell['scf_guesses'] if g['converged']),1e-10)
        cells[tag]=cell
    expected={('ox','ox',None),('rad','rad',None),('rad','ox',None),('rad','ox','ox'),('ox','rad',None),('ox','rad','rad')}
    if set(cells)!=expected:raise ValueError('The six equilibrium/frozen-slow roles differ')
    values=[]
    for initial,final in [('ox','rad'),('rad','ox')]:
        ref=cells[initial,initial,None];cross=cells[final,initial,None];neq=cells[final,initial,initial];eq=cells[final,final,None]
        if neq['reference_density_sha256']!=ref['density_matrix_sha256']:
            raise ValueError('Frozen solvent reference density differs')
        row=dict(total_solution_relaxation_eV=neq['e_eV']-eq['e_eV'],
            equilibrium_geometry_component_eV=cross['e_eV']-eq['e_eV'],
            frozen_slow_increment_eV=neq['e_eV']-cross['e_eV'])
        for key,value in row.items():close(review['directional'][len(values)][key],value)
        values.append(row['total_solution_relaxation_eV'])
    summary=review['summary'];close(summary['lambda_reduction_eV'],values[0]);close(summary['lambda_oxidation_eV'],values[1])
    asymmetry=abs(values[0]-values[1])/(sum(values)/2)
    close(summary['directional_relative_asymmetry'],asymmetry)
    if (summary['equal_curvature_diagnostic_tolerance']!=.05 or asymmetry<=.05
            or 'equal_curvature_approximation_unverified' not in summary['physical_flags']
            or candidate['qualified_symmetric_single_body_eV'] is not None
            or candidate['grafting_free_energy_kJ_mol'] is not None
            or candidate['occupancy_nonadditivity_kJ_mol'] is not None):
        raise ValueError('Unsupported scalar or material property promotion')
    for key,value in zip(('reduction','oxidation'),values):close(candidate['solution_relaxation_eV'][key],value)
    if (figure['figure_sha256']!=hashlib.sha256((folder/'tempo_solution_relaxation.pdf').read_bytes()).hexdigest()
            or figure['source_sha256']!=hashlib.sha256((folder/'tempo_solution_candidate.json').read_bytes()).hexdigest()
            or figure['statistical_intervals_plotted'] or figure['symmetric_scalar_qualified']
            or len(figure['totals_eV'])!=2 or figure['labels']!=['reduction','oxidation']):
        raise ValueError('Figure source or interpretation differs')
    for actual,value in zip(figure['totals_eV'],values):close(actual,value)
    return dict(states=2,SCF_cells=6,reduction_eV=values[0],oxidation_eV=values[1],
                directional_asymmetry=asymmetry,symmetric_scalar_qualified=False,
                experimental_rate_inferred=False,all_axis_decision='pending')
