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
    spacer=next((c for c in data['cycles'] if c['native_cycle_qualified'] and c['definition']['configuration']=='all_same'
                 and c['definition']['sites']==[2,4]),None)
    if spacer is not None:
        sg=spacer['delta_delta_G_by_cutoff_kJ_mol']
        fields.update(IonicSpacerG=sg['50'],IonicSpacerE=spacer['delta_delta_E_kJ_mol'],
            IonicSpacerGLow=min(sg.values()),IonicSpacerGHigh=max(sg.values()))
    alternating=next((c for c in data['cycles'] if c['native_cycle_qualified'] and c['definition']['configuration']=='alternating'
                      and c['definition']['sites']==[2,4]),None)
    if alternating is not None:
        ag=alternating['delta_delta_G_by_cutoff_kJ_mol']
        fields.update(IonicAlternatingSpacerG=ag['50'],IonicAlternatingSpacerGLow=min(ag.values()),IonicAlternatingSpacerGHigh=max(ag.values()))
    return ''.join('\\newcommand{\\'+key+'}{'+f'{value:.2f}'+'}\n' for key,value in fields.items())


def table(data):
    labels={'all_same__occ_none':'Empty backbone','all_same__occ_2':'Site 2',
        'all_same__occ_3':'Site 3','all_same__occ_2_3':'Sites 2 and 3'}
    if 'all_same__occ_2_4' in data['records']:
        labels.update(all_same__occ_4='Site 4',all_same__occ_2_4='Sites 2 and 4')
    if 'alternating__occ_2_4' in data['records']:
        labels.update(alternating__occ_none='Alternating empty',alternating__occ_2='Alternating site 2',
            alternating__occ_4='Alternating site 4',alternating__occ_2_4='Alternating sites 2 and 4')
    lines=[r'\begin{tabular}{lrrrr}\toprule',r'Component & $q$ & $2S+1$ & $\nu_{\min}$ & FD error (\%)\\\midrule']
    for key,label in labels.items():
        row=data['records'][key];err=100*row['maximum_column_error'];mantissa,exponent=f'{err:.2e}'.split('e')
        lines.append(f"{label} & {row['charge']} & {row['mult']} & {row['minimum_frequency_cm']:.2f} & "+
            '$'+mantissa+r'\times10^{'+str(int(exponent))+'}$'+r'\\')
    return '\n'.join([*lines,r'\bottomrule\end{tabular}'])+'\n'


def audit(folder):
    suffix=next((s for s in ('_v3','_v2','') if (folder/('ionic_native_cycles'+s+'.json')).exists()))
    data=json.loads((folder/('ionic_native_cycles'+suffix+'.json')).read_text())
    preparation=json.loads((folder/'ionic_native_cycle_preparation.json').read_text())
    cases={r['case']['key']:r['case'] for r in preparation['cases']}
    complete=sum(c['native_cycle_qualified'] for c in data['cycles'])
    if (data['checked_components']!=len(data['records']) or data['computed_cycles']!=complete
            or data['all_declared_native_cycles_qualified']!=(complete==4) or data['conformer_ensemble_qualified']
            or data['primary_model_admission'] or data['material_admission'] or len(data['cycles'])!=4):
        raise ValueError('Exact checked component/cycle census and field-specific incomplete scopes required')
    for key,row in data['records'].items():
        validate_component(row,cases[key]['charge'])
        if (row['actual_SCF_checkpoint_frames_checked']<1 or not 0<=row['maximum_column_error']<=.01
                or row['minimum_frequency_cm'] < -5):raise ValueError('Declared component physical/frame QC differs')
        if cases[key]['charge']==4 and (row['atoms'],row['actual_SCF_checkpoint_frames_checked'])!=(137,25):
            raise ValueError('Every retained double-graft minimum requires all 25 actual SCF frames')
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
    for checked in (c for c in data['cycles'] if c['native_cycle_qualified']):
        for cutoff,value in checked['delta_delta_G_by_cutoff_kJ_mol'].items():
            independent=sum(Decimal(str(data['records'][key]['G_by_cutoff_eV'][cutoff]))*Decimal(weight)
                            for key,weight in checked['definition']['coefficients'].items())*factor
            if abs(float(independent)-value)>1e-7:raise ValueError('Independent decimal four-state sum differs')
    if (folder/('ionic_native_cycle_values'+suffix+'.tex')).read_text()!=values(data) or (folder/('ionic_native_cycle_table'+suffix+'.tex')).read_text()!=table(data):
        raise ValueError('Ionic cycle macros/table differ from the cached checked inputs')
    return dict(checked_components=len(data['records']),complete_finite_cycles=complete,pending_cycles=4-complete,
        primary_cutoff_cm=50,delta_delta_G_by_cutoff_kJ_mol=cycle['delta_delta_G_by_cutoff_kJ_mol'],
        delta_delta_E_kJ_mol=cycle['delta_delta_E_kJ_mol'],
        native_cycle_arithmetic_replayed=True,original_physical_QC_repeated=False,
        conformer_ensemble_qualified=False,isolated_motif_binding_energy=False,
        electronic_coupling_computed=False,experimental_accuracy_validated=False,material_admission=False)
