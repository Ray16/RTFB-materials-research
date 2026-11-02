"""Replay frozen primary-cohort arithmetic and its scoped reported controls."""
import json
import math

from redox.validation.native_component_sensitivity import reconstruct


def scientific(value):
    if value==0:return '$0$'
    mantissa,exponent=f'{value:.2e}'.split('e')
    return '$'+mantissa+r'\times10^{'+str(int(exponent))+'}$'


def table(rows):
    labels={'graft_ox2':r'Viologen $2+$','graft_ox1':r'Viologen $1+$','graft_neu':'Viologen neutral',
        'precursor_0':'Precursor seed 0','precursor_1':'Precursor seed 1','precursor_2':'Precursor seed 2',
        'fc_neu':'Fc','fc_ox':r'Fc$^+$','local':'Local benzyl chloride site'}
    order=list(labels)
    lines=[r'\begin{tabular}{lrrrrr}\toprule',
        r'State & $q$ & $2S+1$ & $\nu_{\min}$ & FD error (\%) & $\max|\delta G_{\mathrm{nuc}}|$ (eV)\\\midrule']
    indexed={r['key']:r for r in rows}
    for key in order:
        r=indexed[key];value=max(abs(v) for v in r['reference_nuclear_difference_eV'].values())
        lines.append(f"{labels[key]} & {r['charge']} & {r['mult']} & {r['minimum_frequency_cm']:.2f} & "+scientific(100*r['maximum_column_error'])+' & '+scientific(value)+r'\\')
    return '\n'.join([*lines,r'\bottomrule\end{tabular}'])+'\n'


def audit(folder):
    inputs=json.loads((folder/'primary_native_cohort_inputs.json').read_text())
    expected=json.loads((folder/'primary_native_cohort_cycles.json').read_text())
    actual=reconstruct(inputs['package'],inputs['components'])
    if actual!=expected or actual['states']!=10 or len(actual['cycles'])!=6:
        raise ValueError('Frozen full primary-native cohort does not reconstruct')
    rows=json.loads((folder/'primary_native_component_reviews.json').read_text())
    if (len(rows)!=9 or len({r['key'] for r in rows})!=9
            or sum(r['actual_SCF_checkpoint_frames_checked'] for r in rows)!=81
            or sum(r['actual_SCF_wavefunction_spins_checked'] for r in rows)!=81
            or any(r['status']!='primary_molecular_native_control_independently_reconstructed'
                   or not -5<=r['minimum_frequency_cm']<math.inf or not 0<=r['maximum_column_error']<=.01
                   or not 0<=r['maximum_reconstructed_spin_contamination']<=.1 for r in rows)):
        raise ValueError('Declared native physical/spin control coverage differs')
    if (folder/'primary_native_cohort_table.tex').read_text()!=table(rows):
        raise ValueError('Primary-native control table differs from its source reviews')
    potentials=max(abs(v) for r in actual['cycles'] if r['unit']=='V' for v in r['corrected_minus_source_by_cutoff'].values())
    reactions=max(abs(v) for r in actual['cycles'] if r['unit']=='kJ/mol' for v in r['corrected_minus_source_by_cutoff'].values())
    if potentials>2.45e-8 or reactions>1.55e-6:raise ValueError('Full-cohort changes exceed the stated rounded bounds')
    return dict(checked_states=10,new_molecular_Hessians=9,analytic_monatomic_components=1,
        checked_SCF_frames=81,independently_reconstructed_wavefunction_spins=81,complete_cycles=6,
        maximum_potential_change_V=potentials,maximum_reaction_change_kJ_mol=reactions,
        cached_cycle_arithmetic_replayed=True,original_raw_physical_review_repeated=False,
        experimental_accuracy_validated=False,occupied_polymer_free_energy_cycle_qualified=False)
