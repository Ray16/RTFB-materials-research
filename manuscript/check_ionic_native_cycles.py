"""Replay frozen finite-minimum ionic cycles without promoting missing states."""
import json
from decimal import Decimal
from redox.workflows.ionic_native_cycles import native_cycle,validate_component


def values(data):
    cycle=next(c for c in data['cycles'] if c['native_cycle_qualified'])
    g=cycle['delta_delta_G_by_cutoff_kJ_mol'];large=data['records']['all_same__occ_2_3']
    fields=dict(IonicAdjacentG=g['50'],IonicAdjacentE=cycle['delta_delta_E_kJ_mol'],
        IonicAdjacentGLow=min(g.values()),IonicAdjacentGHigh=max(g.values()),
        IonicLargeMinimum=large['minimum_frequency_cm'])
    return ''.join('\\newcommand{\\'+key+'}{'+f'{value:.2f}'+'}\n' for key,value in fields.items())


def table(data):
    labels={'all_same__occ_none':'Empty backbone','all_same__occ_2':'Site 2',
        'all_same__occ_3':'Site 3','all_same__occ_2_3':'Sites 2 and 3'}
    lines=[r'\begin{tabular}{lrrrr}\toprule',r'Component & $q$ & $2S+1$ & $\nu_{\min}$ & FD error (\%)\\\midrule']
    for key,label in labels.items():
        row=data['records'][key];err=100*row['maximum_column_error'];mantissa,exponent=f'{err:.2e}'.split('e')
        lines.append(f"{label} & {row['charge']} & {row['mult']} & {row['minimum_frequency_cm']:.2f} & "+
            '$'+mantissa+r'\times10^{'+str(int(exponent))+'}$'+r'\\')
    return '\n'.join([*lines,r'\bottomrule\end{tabular}'])+'\n'


def audit(folder):
    data=json.loads((folder/'ionic_native_cycles.json').read_text())
    preparation=json.loads((folder/'ionic_native_cycle_preparation.json').read_text())
    cases={r['case']['key']:r['case'] for r in preparation['cases']}
    if (data['checked_components']!=9 or len(data['records'])!=9 or data['computed_cycles']!=1
            or data['all_declared_native_cycles_qualified'] or data['conformer_ensemble_qualified']
            or data['primary_model_admission'] or data['material_admission'] or len(data['cycles'])!=4):
        raise ValueError('Exactly one scoped checked ionic cycle and three pending cycles required')
    for key,row in data['records'].items():
        validate_component(row,cases[key]['charge'])
        if (row['actual_SCF_checkpoint_frames_checked']<1 or not 0<=row['maximum_column_error']<=.01
                or row['minimum_frequency_cm'] < -5):raise ValueError('Declared component physical/frame QC differs')
    large=data['records']['all_same__occ_2_3']
    if (large['atoms'],large['charge'],large['mult'],large['actual_SCF_checkpoint_frames_checked'])!=(137,4,1,25):
        raise ValueError('Exact largest ionic component and all 25 frames required')
    for row in data['cycles']:
        actual=native_cycle(row['definition'],data['records'],cases)
        if actual!={k:v for k,v in row.items() if k!='definition'}:
            raise ValueError('Frozen ionic cycle does not reconstruct from its exact components')
    cycle=next(c for c in data['cycles'] if c['native_cycle_qualified'])
    if cycle['definition']['configuration']!='all_same' or cycle['definition']['sites']!=[2,3]:
        raise ValueError('Only the first adjacent all-same ionic cycle is reported in this epoch')
    factor=Decimal('1.602176634e-19')*Decimal('6.02214076e23')/Decimal(1000)
    for cutoff,value in cycle['delta_delta_G_by_cutoff_kJ_mol'].items():
        independent=sum(Decimal(str(data['records'][key]['G_by_cutoff_eV'][cutoff]))*Decimal(weight)
                        for key,weight in cycle['definition']['coefficients'].items())*factor
        if abs(float(independent)-value)>1e-7:raise ValueError('Independent decimal four-state sum differs')
    if (folder/'ionic_native_cycle_values.tex').read_text()!=values(data) or (folder/'ionic_native_cycle_table.tex').read_text()!=table(data):
        raise ValueError('Ionic cycle macros/table differ from the cached checked inputs')
    return dict(checked_components=9,complete_finite_cycles=1,pending_cycles=3,
        primary_cutoff_cm=50,delta_delta_G_by_cutoff_kJ_mol=cycle['delta_delta_G_by_cutoff_kJ_mol'],
        delta_delta_E_kJ_mol=cycle['delta_delta_E_kJ_mol'],
        native_cycle_arithmetic_replayed=True,original_physical_QC_repeated=False,
        conformer_ensemble_qualified=False,isolated_motif_binding_energy=False,
        electronic_coupling_computed=False,experimental_accuracy_validated=False,material_admission=False)
