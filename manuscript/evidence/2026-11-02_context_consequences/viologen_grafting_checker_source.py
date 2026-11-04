"""Independent arithmetic for the frozen exact-viologen quaternization cycle."""
import json
import math


def close(a,b):
    if any(type(v) not in (int,float) or not math.isfinite(v) for v in (a,b)) or abs(a-b)>1e-8:
        raise ValueError('Frozen quaternization value does not reconstruct')


def values_tex(result):
    values=result['deltaG_kJ_mol_by_cutoff']
    return ''.join('\\newcommand{\\'+name+'}{'+f'{value:.2f}'+'}\n' for name,value in [
        ('ViologenGraftDG',values['50']),('ViologenGraftRotorLow',min(values.values())),
        ('ViologenGraftRotorHigh',max(values.values())),('ViologenGraftSeedSpan',result['finite_precursor_deltaG_span_kJ_mol'])])


def qc_table(result):
    best=min(r['G_by_cutoff_eV']['50'] for r in result['precursor_conformers']);factor=1.602176634e-19*6.02214076e23/1000
    text='\\begin{tabular}{rrrr}\\toprule\nSeed & Relative precursor $G$ & Grafting $\\Delta G$ & $\\nu_{\\min}$\\\\\\midrule\n'
    for row,cycle in zip(result['precursor_conformers'],result['per_precursor_cycle']):
        text+=f"{row['conformer_index']} & {(row['G_by_cutoff_eV']['50']-best)*factor:.2f} & {cycle['deltaG_kJ_mol_by_cutoff']['50']:.2f} & {row['minimum_internal_frequency_cm']:.2f}\\\\\n"
    return text+'\\bottomrule\\end{tabular}\n'


def audit(folder):
    path=folder/'viologen_grafting_result.json'
    if not path.exists():return None
    result=json.loads(path.read_text());components=result['components'];precursors=result['precursor_conformers']
    if (result['status']!='exact_local_quaternization_native_cycle_checked_finite_sampling'
            or result['product_inchikey']!='PBNXFHGCSKYNHC-UHFFFAOYSA-N'
            or result['reaction_class']!='N_quaternization' or result['delta_n_solutes']!=0
            or result['base_term_kJ_mol']!=0 or result['standard_state_correction_kJ_mol']!=0
            or result['material_admission'] or result['global_conformer_convergence_proven']
            or result['ion_pairing_computed'] or result['counterion_exchange_computed'] or result['experimental_reflux_modelled']
            or result['experimental_grafting_accuracy_validated'] or result['temperature_K']!=298.15 or result['solvent']!='acetonitrile'):
        raise ValueError('Quaternization identity, reaction convention or scope differs')
    expected={'precursor':('C[n+]1ccc(-c2ccncc2)cc1',1),'site':('Cc1ccc(CCl)cc1',0),
        'product':('Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1',2),'byproduct':('[Cl-]',-1)}
    if set(components)!=set(expected) or len(precursors)!=3 or len(result['per_precursor_cycle'])!=3:
        raise ValueError('Four reaction components and three complete precursor seeds required')
    for role,row in components.items():
        if (row['smiles'],row['charge'])!=expected[role] or row['mult']!=1:raise ValueError('Quaternization species differs')
    if sorted(r['conformer_index'] for r in precursors)!=[0,1,2]:raise ValueError('Three distinct fixed seeds required')
    for row in [*components.values(),*precursors]:
        if set(row['G_by_cutoff_eV'])!={'25','50','100'} or row['mult']!=1:raise ValueError('Three matched singlet nuclear partitions required')
        for c,v in row['G_by_cutoff_eV'].items():
            close(v,math.fsum([row['active_E_smd_eV'],row['nuclear_by_cutoff_eV'][c],row['spin_partition_eV']]))
    for row in precursors:
        if (row['smiles'],row['charge'])!=expected['precursor'] or not 0<=row['maximum_sampled_column_error']<=.01:
            raise ValueError('Precursor identity or sampled derivative gate differs')
        if not -5<=row['minimum_internal_frequency_cm'] or not 0<=row['maximum_gradient_Ha_Bohr']<=.00045:
            raise ValueError('Precursor minimum or stationary gate failed')
    selected=min(precursors,key=lambda r:r['G_by_cutoff_eV']['50'])
    if result['selected_precursor_conformer_index']!=selected['conformer_index'] or components['precursor']!=selected:
        raise ValueError('Lowest precursor G50 selection differs')
    factor=1.602176634e-19*6.02214076e23/1000;close(result['conversion_kJ_mol_per_eV'],factor)
    for row,cycle in zip(precursors,result['per_precursor_cycle']):
        for c in ('25','50','100'):
            value=math.fsum([components['product']['G_by_cutoff_eV'][c],components['byproduct']['G_by_cutoff_eV'][c],
                -components['site']['G_by_cutoff_eV'][c],-row['G_by_cutoff_eV'][c]])*factor
            close(cycle['deltaG_kJ_mol_by_cutoff'][c],value)
            if row['conformer_index']==selected['conformer_index']:close(result['deltaG_kJ_mol_by_cutoff'][c],value)
    span=max(c['deltaG_kJ_mol_by_cutoff']['50'] for c in result['per_precursor_cycle'])-min(c['deltaG_kJ_mol_by_cutoff']['50'] for c in result['per_precursor_cycle'])
    close(result['finite_precursor_deltaG_span_kJ_mol'],span)
    if (folder/'viologen_grafting_values.tex').read_text()!=values_tex(result) or (folder/'viologen_grafting_QC_table.tex').read_text()!=qc_table(result):
        raise ValueError('Frozen quaternization LaTeX numbers differ')
    return dict(status='exact_viologen_quaternization_frozen_numbers_reconstructed',
        deltaG_kJ_mol_by_cutoff=result['deltaG_kJ_mol_by_cutoff'],precursor_seed_span_kJ_mol=span,
        material_admission=False,global_conformer_convergence_proven=False)
