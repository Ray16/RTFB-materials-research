"""Independently reconstruct the exact viologen ladder in the frozen paper."""
import hashlib
import json
import math


def qc_table(qc):
    names={'ox2':'Dication','ox1':'Radical cation','neu':'Neutral'}
    rows={r['state']:r for r in qc['states']}
    text='\\begin{tabular}{lrrrr}\\toprule\nState & Charge & Multiplicity & Initial $\\nu_{\\min}$ & Accepted $\\nu_{\\min}$\\\\\\midrule\n'
    for state in ('ox2','ox1','neu'):
        r=rows[state]
        text+=f"{names[state]} & {r['charge']} & {r['mult']} & {r['initial_minimum_frequency_cm']:.2f} & {r['final_minimum_frequency_cm']:.2f}\\\\\n"
    return text+'\\bottomrule\\end{tabular}\n'


def close(actual,expected,tolerance=1e-8):
    if (any(type(v) not in (int,float) or not math.isfinite(v) for v in (actual,expected))
            or abs(actual-expected)>tolerance):
        raise ValueError('Frozen viologen quantity does not reconstruct')


def audit(folder,figures=None):
    load=lambda name:json.loads((folder/name).read_text())
    candidate=load('viologen_native_candidate.json');numeric=load('viologen_native_numeric_inputs.json')
    qc=load('viologen_native_QC.json');figure=load('viologen_native_figure.json')
    index=load('current_registry_identity_index.json')['graphs']
    if (candidate['status']!='exact_molecular_native_redox_published_with_field_QC'
            or candidate['product_inchikey']!='PBNXFHGCSKYNHC-UHFFFAOYSA-N'
            or candidate['all_axis_decision']['status']!='pending' or candidate['material_admission']
            or candidate['experimental_exact_graft_accuracy_validated'] or candidate['global_conformer_convergence_proven']
            or any(candidate[k] is not None for k in ('grafting_free_energy_kJ_mol','occupancy_nonadditivity_kJ_mol','electronic_coupling_eV'))):
        raise ValueError('Exact candidate identity or assessment scope differs')
    if (len(index)!=904 or len({r['graph_key'] for r in index})!=904
            or any(hashlib.sha256(r['product_smiles'].encode()).hexdigest()!=r['graph_key'] for r in index)
            or not any(r['graph_key']==candidate['registry_graph_key'] and r['inchikey']==candidate['product_inchikey'] for r in index)):
        raise ValueError('Maintained registry identity coverage differs')
    rows=numeric['states'];models={};states={};reference={}
    if len(rows)!=5:raise ValueError('Three candidate and two matched reference states required')
    if len(candidate['redox']['events'])!=2 or candidate['redox']['fitted_offset_V']!=0.:
        raise ValueError('Two unfitted reduction events required')
    for row in rows:
        if (row['active_protocol']['basis']!='def2-tzvpd' or row['active_protocol']['xc']!='wb97m-v'
                or row['active_protocol']['solvent']!='acetonitrile' or row['active_protocol']['solvent_model']!='SMD'
                or row['native_protocol']['xc']!='r2scan' or row['native_protocol']['basis']!='def2-svp'
                or row['native_protocol']['disp']!='d4' or row['native_protocol']['temperature_K']!=298.15):
            raise ValueError('Candidate/reference protocols differ from the reported model')
        for key in ('active_protocol','native_protocol'):
            encoded=json.dumps(row[key],sort_keys=True)
            if models.setdefault(key,encoded)!=encoded:raise ValueError('State protocols are not matched')
        if set(row['nuclear_by_cutoff_eV'])!={'25','50','100'} or set(row['G_by_cutoff_eV'])!={'25','50','100'}:
            raise ValueError('Three rotor comparisons required for every state')
        close(row['spin_partition_eV'],-8.617333262145e-5*298.15*math.log(row['mult']))
        for cutoff,value in row['G_by_cutoff_eV'].items():
            close(value,math.fsum([row['active_E_smd_eV'],row['nuclear_by_cutoff_eV'][cutoff],row['spin_partition_eV']]))
        target=states if row['role']=='candidate' else reference if row['role']=='reference' else None
        if target is None or row['state'] in target:raise ValueError('Duplicated or unknown native role')
        if row['role']=='candidate' and row['product_inchikey']!=candidate['product_inchikey']:
            raise ValueError('A native component belongs to another exact candidate')
        target[row['state']]=row
    if (set(states)!={'ox2','ox1','neu'} or set(reference)!={'neu','ox'}
            or [(states[s]['charge'],states[s]['mult']) for s in ('ox2','ox1','neu')]!=[(2,1),(1,2),(0,1)]
            or [(reference[s]['charge'],reference[s]['mult']) for s in ('neu','ox')]!=[(0,1),(1,2)]):
        raise ValueError('Charge/spin assignment differs')
    values=[]
    for i,(ox,red) in enumerate([('ox2','ox1'),('ox1','neu')]):
        event=candidate['redox']['events'][i]
        if event['event']!=ox+'->'+red or event['electrons']!=1:raise ValueError('Redox event differs')
        for cutoff in ('25','50','100'):
            potential=math.fsum([states[ox]['G_by_cutoff_eV'][cutoff],-states[red]['G_by_cutoff_eV'][cutoff],
                -reference['ox']['G_by_cutoff_eV'][cutoff],reference['neu']['G_by_cutoff_eV'][cutoff]])
            close(event['potential_V_vs_Fc_by_cutoff'][cutoff],potential)
        values.append(event['potential_V_vs_Fc_by_cutoff']['50'])
    factor=1.602176634e-19*6.02214076e23/1000
    for cutoff in ('25','50','100'):
        spacing=candidate['redox']['events'][0]['potential_V_vs_Fc_by_cutoff'][cutoff]-candidate['redox']['events'][1]['potential_V_vs_Fc_by_cutoff'][cutoff]
        close(candidate['redox']['spacing_V_by_cutoff'][cutoff],spacing)
        close(candidate['redox']['disproportionation_kJ_mol_by_cutoff'][cutoff],spacing*factor)
    if (len(qc['states'])!=3 or {r['state'] for r in qc['states']}!=set(states)
            or qc['conformer_ensemble_convergence_proven']):raise ValueError('QC coverage/scope differs')
    for row in qc['states']:
        for key,lower,upper in [('final_minimum_frequency_cm',-5,math.inf),('maximum_gradient_Ha_Bohr',0,.00045),
                ('maximum_Hessian_asymmetry',0,1e-6),('maximum_sampled_column_error',0,.01)]:
            value=row[key]
            if not math.isfinite(value) or not lower<=value<=upper:raise ValueError('Reported native QC failed')
    if (folder/'viologen_native_QC_table.tex').exists() and (folder/'viologen_native_QC_table.tex').read_text()!=qc_table(qc):
        raise ValueError('Generated viologen QC table differs')
    pdf=(folder if figures is None else figures)/'viologen_native_redox.pdf'
    if (hashlib.sha256(pdf.read_bytes()).hexdigest()!=figure['figure_sha256']
            or hashlib.sha256((folder/'viologen_native_candidate.json').read_bytes()).hexdigest()!=figure['source_sha256']
            or figure['statistical_intervals_plotted'] or figure['experimental_exact_graft_accuracy_claimed']):
        raise ValueError('Figure source or interpretation differs')
    for actual,expected in zip(figure['primary_potentials_V'],values):close(actual,expected)
    macros={'MolecularFirstReduction':f'{values[0]:+.3f}','MolecularSecondReduction':f'{values[1]:+.3f}',
        'MolecularSpacing':f"{candidate['redox']['spacing_V_by_cutoff']['50']:.3f}",
        'MolecularDisproportionation':f"{candidate['redox']['disproportionation_kJ_mol_by_cutoff']['50']:.2f}"}
    expected=''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in macros.items())
    if (folder/'viologen_native_values.tex').read_text()!=expected:raise ValueError('Generated viologen macros differ')
    return dict(candidate_states=3,matched_reference_states=2,registry_graphs=904,potentials_V_vs_Fc=values,
        disproportionation_kJ_mol=candidate['redox']['disproportionation_kJ_mol_by_cutoff']['50'],
        fitted_offset_V=0.,material_admission=False)
