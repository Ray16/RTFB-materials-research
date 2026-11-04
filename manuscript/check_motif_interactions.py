"""Replay frozen motif/whole-chain E definitions; do not repeat physical QC."""
from decimal import Decimal
import json
from pathlib import Path

from redox.core.hashing import sha256_file

FACTOR=Decimal('27.211386245988')*Decimal('1.602176634e-19')*Decimal('6.02214076e23')/Decimal(1000)
KEYS=('all_same__occ_2_3','all_same__occ_2_4','alternating__occ_2_3','alternating__occ_2_4')


def reconstruct(folder):
    collection=json.loads((folder/'motif_energy_collection.json').read_text())
    whole=json.loads((folder/'motif_parent_stationary.json').read_text())
    expected=json.loads((folder/'motif_definition_comparison.json').read_text())
    if collection['checked_triplets']!=12 or collection['independent_review_errors'] or len(collection['states'])!=12:
        raise ValueError('All twelve independently checked frozen triplets are required')
    local={};ranges={};frames=0
    for item in collection['states']:
        key=item['key'];report=json.loads((folder/('motif_review_'+key+'.json')).read_text())
        if (report['status']!='frozen_charged_motif_triplet_independently_checked' or report['actual_SCF_frames_checked']!=6
                or report['native_G_qualified'] or report['basis_superposition_corrected'] or report['electronic_coupling_computed']
                or report['solvent_model']!=whole['model'] or report['numerical_implementation']!=whole['numerical_implementation']):
            raise ValueError('Exact E-only scope, frame count and model identity required')
        positions={}
        for label,charge,count in [('A',2,28),('B',2,28),('AB',4,56)]:
            path=folder/('motif_'+key+'_'+label+'.xyz');lines=path.read_text().splitlines();entries=[line.split() for line in lines[2:]]
            if int(lines[0])!=count or len(entries)!=count:raise ValueError('Frozen motif XYZ coverage differs')
            positions[label]=entries
            for medium in ('gas','acetonitrile'):
                state=report['exact_energies_Ha'][label][medium]
                if (state['charge'],state['mult'],state['SCF_converged'])!=(charge,1,True) or sha256_file(path)!=state['geometry_sha256']:
                    raise ValueError('Frozen electronic component state or geometry differs')
        if positions['AB']!=positions['A']+positions['B']:raise ValueError('Pair and monomer frames are not exact concatenations')
        values={}
        for medium in ('gas','acetonitrile'):
            energy=report['exact_energies_Ha'];value=(Decimal(str(energy['AB'][medium]['E_Ha']))
                -Decimal(str(energy['A'][medium]['E_Ha']))-Decimal(str(energy['B'][medium]['E_Ha'])))*FACTOR
            if abs(float(value)-report['interaction_E_kJ_mol_by_medium'][medium])>1e-6:
                raise ValueError('Independent decimal local interaction sum differs')
            values[medium]=report['interaction_E_kJ_mol_by_medium'][medium]
        local[(report['parent_key'],report['cap_length_A'])]=values;frames+=6
    result=[]
    for key in KEYS:
        for length in (1.07,1.09,1.11):
            if (key,length) not in local:raise ValueError('Missing frame or cap-sensitivity point')
        cycle=next(row for row in whole['cycles'] if row['definition']['configuration']==key.split('__')[0]
            and row['definition']['sites']==list(map(int,key.split('__occ_')[1].split('_'))))
        value=sum(Decimal(str(whole['records'][state]['E_Ha']))*Decimal(weight)
            for state,weight in cycle['definition']['coefficients'].items())*FACTOR
        if abs(float(value)-cycle['delta_delta_E_kJ_mol'])>1e-6:raise ValueError('Independent whole-chain E sum differs')
        central=local[(key,1.09)];span=max(local[(key,l)]['acetonitrile'] for l in (1.07,1.09,1.11))-min(local[(key,l)]['acetonitrile'] for l in (1.07,1.09,1.11))
        row=next(row for row in expected['states'] if row['parent_key']==key)
        if (row['whole_chain_nonadditivity_E_kJ_mol']!=cycle['delta_delta_E_kJ_mol']
                or row['frozen_methyl_motif_interaction_E_kJ_mol']!=central['acetonitrile']):
            raise ValueError('Definition-comparison table differs from its raw-backed E parents')
        result.append(dict(key=key,whole_E=cycle['delta_delta_E_kJ_mol'],motif_E=central['acetonitrile'],gas_E=central['gas'],cap_span=span))
    return result,frames


def values(rows):
    numbers=dict(MotifAdjacentE=rows[0]['motif_E'],MotifMinimumE=min(r['motif_E'] for r in rows),
        MotifMaximumE=max(r['motif_E'] for r in rows),MotifCapSpan=max(r['cap_span'] for r in rows))
    return ''.join('\\newcommand{\\'+key+'}{'+f'{value:.3f}'+'}\n' for key,value in numbers.items())


def table(rows):
    labels=['All-same, adjacent','All-same, one spacer','Alternating, adjacent','Alternating, one spacer']
    lines=[r'\begin{tabular}{lrrrr}\toprule',r'Parent frame & Whole-chain $\Delta\Delta E$ & Motif $E_{\rm int}^{\rm MeCN}$ & Motif $E_{\rm int}^{\rm gas}$ & Cap span (MeCN)\\\midrule']
    for label,row in zip(labels,rows):lines.append(f"{label} & {row['whole_E']:.2f} & {row['motif_E']:.2f} & {row['gas_E']:.2f} & {row['cap_span']:.4f}"+r'\\')
    return '\n'.join([*lines,r'\bottomrule\end{tabular}'])+'\n'


def audit(folder):
    rows,frames=reconstruct(folder)
    if (folder/'motif_interaction_values.tex').read_text()!=values(rows) or (folder/'motif_interaction_table.tex').read_text()!=table(rows):
        raise ValueError('Frozen motif TeX values/table differ')
    return dict(checked_triplets=12,SCF_frames_represented=frames,parent_frames=4,
        opposite_sign_frames=sum((r['whole_E']>=0)!=(r['motif_E']>=0) for r in rows),
        minimum_motif_MeCN_E_kJ_mol=min(r['motif_E'] for r in rows),maximum_motif_MeCN_E_kJ_mol=max(r['motif_E'] for r in rows),
        maximum_cap_span_MeCN_kJ_mol=max(r['cap_span'] for r in rows),original_physical_QC_repeated=False,
        binding_free_energy_supplied=False,unique_backbone_decomposition=False,electronic_coupling_computed=False)
