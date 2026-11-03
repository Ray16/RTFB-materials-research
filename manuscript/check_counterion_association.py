"""Replay finite primary chloride association and duplicate scope from frozen data."""
from copy import deepcopy
from decimal import Decimal
import json
import math

import numpy as np

from redox.core.hashing import sha256_file
from redox.properties.counterion import association_cycle
from redox.properties.conformer_geometry import alignment_molecule,symmetry_rmsd_matrix
from redox.properties.conformer_ensemble import cluster_minima

FACTOR=Decimal('1.602176634e-19')*Decimal('6.02214076e23')/Decimal(1000)


def close(actual,expected,tolerance=1e-7):
    if (type(actual) not in (int,float) or not math.isfinite(actual)
            or abs(actual-float(expected))>tolerance):
        raise ValueError('Frozen chloride association arithmetic differs')


def reconstruct(folder):
    prefix='counterion_v2_' if (folder/'counterion_v2_inputs.json').exists() else 'counterion_'
    def artifact(name):return folder/(prefix+name)
    catalog=json.loads(artifact('inputs.json').read_text())
    parents=json.loads((artifact('primary_inputs.json')).read_text())
    if (folder/'primary_native_cohort_inputs.json').exists():
        if parents!=json.loads((folder/'primary_native_cohort_inputs.json').read_text()):
            raise ValueError('Counterion and manuscript primary parents differ')
    thermal=json.loads((artifact('parent_thermal_protocol.json')).read_text())
    if (catalog['checked_outcomes'] not in (4,6) or len(catalog['outcomes'])!=catalog['checked_outcomes'] or catalog['declared_original_starts']!=7
            or catalog['conformer_ensemble_qualified'] or catalog['experimental_accuracy_validated']
            or catalog['material_admission'] or catalog['reaction_rate_or_yield_inferred']):
        raise ValueError('Exact finite four-outcome scope required')
    context=deepcopy(parents['components']['graft_ox2']['context'])
    if context!=parents['components']['chloride']['context'] or context['temperature_K']!=298.15:
        raise ValueError('Matched primary parent context and temperature required')
    reference=dict(energy=context['reference']['energy'],
        ions='explicit incoming chloride; omitted separated iodide spectator on both sides',
        numerical_implementation=context['reference']['numerical_implementation'])
    context['reference']=reference
    states={row['id']:row for row in parents['package']['states']}
    shared={}
    for source_id,role in [('graft_ox2','graft'),('chloride','chloride')]:
        component=parents['components'][source_id]
        if component['qc']!='independently_checked' or not component['native_component_requalified']:
            raise ValueError('Independently requalified primary parents required')
        shared[role]=dict(state=dict(identity=source_id,smiles=states[source_id]['smiles'],
            charge=component['charge'],multiplicity=component['mult']),context=deepcopy(context),
            G_by_cutoff_eV=component['G_by_cutoff_eV'],qc=component['qc'],
            sources=component['sources'],thermodynamic_entities=1)
    rows=[];molecules=[];electronic=[]
    for item in catalog['outcomes']:
        report=json.loads((folder/item['review']).read_text());geometry=folder/item['geometry']
        if (report['status']!='primary_counterion_thermochemistry_independently_checked'
                or (report['charge'],report['mult'],report['formal_core_charge'])!=(1,1,2)
                or report['native_G_qualified'] is not True or report['conformer_ensemble_qualified']
                or report['experimental_populations_inferred'] or report['material_admission']
                or report['actual_native_SCF_frames_checked']!=9
                or report['electronic_protocol']!=context['electronic_protocol']
                or report['native_protocol']!=context['thermal_definition']['native_protocol']
                or report['thermal_protocol']!=thermal or report['temperature_K']!=298.15
                or report['numerical_implementation']!=reference['numerical_implementation']
                or report['standard_state']!=thermal['pressure']
                or report['geometry_sha256']!=sha256_file(geometry)
                or not -5<=report['minimum_frequency_cm']<float('inf')
                or not 0<=report['maximum_column_error']<=.01):
            raise ValueError('Exact checked primary chloride component/state/model required')
        for cutoff in ('25','50','100'):
            close(report['G_by_cutoff_eV'][cutoff],report['electronic_eV']+
                  report['nuclear_by_cutoff_eV'][cutoff]+report['spin_partition_eV'],1e-8)
        pair=dict(state=dict(identity='oxidized_graft_chloride:'+report['key'],smiles=report['smiles'],
            charge=1,multiplicity=1),context=deepcopy(context),G_by_cutoff_eV=report['G_by_cutoff_eV'],
            qc='independently_checked',sources={item['review']:sha256_file(folder/item['review'])},thermodynamic_entities=1)
        computed=association_cycle(dict(paired=pair,**shared))
        ledger=json.loads((folder/item['ledger']).read_text())
        recorded=next(row for row in ledger['states'] if row['key']==report['key'])
        if (ledger['reaction_reference']!=reference or ledger['delta_n_solutes']!=-1
                or ledger['ion_pair_populations_inferred'] or ledger['covalent_quaternization_cycle_modified']
                or recorded['coefficients']!={'paired':1,'graft':-1,'chloride':-1}
                or recorded['delta_n_solutes']!=-1 or recorded['population_weights_inferred']):
            raise ValueError('Original association scope and balanced two-to-one convention required')
        temperature=Decimal(str(context['temperature_K']))
        correction=-Decimal('0.00831446261815324')*temperature*(Decimal('0.08205736608096')*temperature).ln()
        for cutoff in ('25','50','100'):
            source=(Decimal(str(report['G_by_cutoff_eV'][cutoff]))
                    -Decimal(str(shared['graft']['G_by_cutoff_eV'][cutoff]))
                    -Decimal(str(shared['chloride']['G_by_cutoff_eV'][cutoff])))*FACTOR
            close(computed['source_convention_deltaG_kJ_mol'][cutoff],source)
            close(recorded['source_convention_deltaG_kJ_mol'][cutoff],source)
            close(computed['one_M_deltaG_kJ_mol'][cutoff],source+correction)
            close(recorded['one_M_deltaG_kJ_mol'][cutoff],source+correction)
        close(recorded['one_atm_to_one_M_correction_kJ_mol'],correction)
        rows.append(dict(key=item['key'],label=item['label'],minimum_frequency_cm=report['minimum_frequency_cm'],
            maximum_column_error=report['maximum_column_error'],G=recorded['one_M_deltaG_kJ_mol'],
            source_G=recorded['source_convention_deltaG_kJ_mol'],correction=float(correction)))
        molecules.append(alignment_molecule(report['smiles'],geometry))
        electronic.append(dict(key=item['key'],E_eV=report['electronic_eV']))
    expected_keys=['N0_minus','N0_plus','N1_above','N1_below']+(['N0_below_minus','N0_outer_minus'] if catalog['checked_outcomes']==6 else [])
    if [row['key'] for row in rows]!=expected_keys:
        raise ValueError('Four exact source outcomes must retain their individual identities')
    duplicate=json.loads((artifact('duplicate_review.json')).read_text())
    matrix,matches=symmetry_rmsd_matrix(molecules)
    if (not np.allclose(matrix,duplicate['proper_rotation_symmetry_heavy_RMSD_A'],atol=1e-8,rtol=0)
            or matches!=duplicate['automorphisms'] or duplicate['populations_inferred']
            or duplicate['conformer_ensemble_qualified']):
        raise ValueError('Actual frozen geometry duplicate comparison differs')
    for cutoff in (.05,.1,.2):
        groups=cluster_minima(electronic,matrix,cutoff,.2)
        if groups!=duplicate['groups_by_RMSD_cutoff_A'][str(cutoff)] or len(groups)!=len(rows)-1:
            raise ValueError('Four outcomes must not be counted as four independent conformational wells')
    return rows,float(matrix[0,1])


def values(rows):
    numbers=dict(ChlorideAssocMin=min(r['G']['50'] for r in rows),
        ChlorideAssocMax=max(r['G']['50'] for r in rows),ChlorideCorrection=rows[0]['correction'])
    numbers.update({name:row['G']['50'] for name,row in zip(
        ['ChlorideNZeroMinus','ChlorideNZeroPlus','ChlorideNOneAbove','ChlorideNOneBelow','ChlorideNZeroBelow','ChlorideNZeroOuter'],rows)})
    return ''.join('\\newcommand{\\'+name+'}{'+f'{value:.2f}'+'}\n' for name,value in numbers.items())


def table(rows):
    lines=[r'\begin{tabular}{lrrrrr}\toprule',
        r'Outcome & $\nu_{\min}$ & FD error (\%) & $\Delta G_{25}$ & $\Delta G_{50}$ & $\Delta G_{100}$\\\midrule']
    for row in rows:
        lines.append(row['label']+' & '+f"{row['minimum_frequency_cm']:.2f} & {100*row['maximum_column_error']:.4f}"+
            ''.join(' & '+f"{row['G'][c]:.2f}" for c in ('25','50','100'))+r'\\')
    return '\n'.join([*lines,r'\bottomrule\end{tabular}'])+'\n'


def audit(folder):
    rows,rmsd=reconstruct(folder)
    prefix='counterion_v2_' if (folder/'counterion_v2_inputs.json').exists() else 'counterion_'
    if (folder/(prefix+'values.tex')).read_text()!=values(rows) or (folder/(prefix+'table.tex')).read_text()!=table(rows):
        raise ValueError('Counterion numerical TeX inputs differ from frozen evidence')
    return dict(status='frozen_primary_contact_outcomes_replayed',checked_outcomes=len(rows),duplicate_groups=len(rows)-1,
        same_group_N0_heavy_RMSD_A=rmsd,minimum_primary_association_1M_kJ_mol=min(r['G']['50'] for r in rows),
        maximum_primary_association_1M_kJ_mol=max(r['G']['50'] for r in rows),
        original_physical_QC_repeated=False,geometry_duplicate_comparison_repeated=True,
        conformer_ensemble_qualified=False,experimental_accuracy_validated=False,material_admission=False)
