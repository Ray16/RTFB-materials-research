#!/usr/bin/env python3
"""Replay the exact-PMDI finite neighbor cycle from contained cached records."""
import argparse
import importlib.util
import json
import math
from pathlib import Path


def replay(root):
    root=Path(root).resolve()
    manifest=json.loads((root/'manifest.json').read_text())
    if (manifest['synthesis_admission'] or manifest['conformer_convergence_established']
            or manifest['isolated_motif_interaction_computed'] or manifest['electronic_coupling_computed']):
        raise ValueError('Cached finite-cycle replay cannot promote its scientific scope')
    if any(type(row['checked_G_eV']) not in (int,float) or not math.isfinite(row['checked_G_eV'])
           for row in manifest['states']):raise ValueError('Finite checked state G required')
    names=[row['id'] for row in manifest['states']]
    if len(names)!=14 or len(set(names))!=14:raise ValueError('Fourteen distinct raw state identities required')
    def number(value):
        if type(value) not in (int,float) or not math.isfinite(value):raise ValueError('Finite source value required')
        return value
    spec=importlib.util.spec_from_file_location('contained_cached_state_reader',root/'cached_reader.py')
    reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    primary=reader.reconstruct(root)
    def load(relative):
        path=(root/relative).resolve()
        if not path.is_relative_to(root):raise ValueError('Cached artifact escapes package')
        return json.loads(path.read_text())
    source=load(manifest['neighbor_source'])
    if (source['candidate_id']!='design_pmdi_ethyl' or source['solvent']!='acetonitrile'
            or source['conformer_convergence_established'] or source['synthesis_admission']):
        raise ValueError('Expected finite exact-PMDI cycle source required')
    aliases={'pair':'polymer_pair','empty':'empty','site2':'site2','site3':'site3'}
    states={row['id']:row for row in manifest['states']}
    conversion=number(manifest['SI_energy_conversion_kJ_mol_per_eV'])
    if conversion<=0 or number(manifest['temperature_K'])!=298.15:raise ValueError('Expected positive conversion and 298.15 K')
    standard=-.00831446261815324*manifest['temperature_K']*math.log(.08205736608096*manifest['temperature_K'])
    cutoffs={}
    for cutoff in ('25','50','100'):
        values={}
        for key,row in states.items():
            energy=load(row['files']['energy']);nuclear=load(row['files']['partition'])
            correction=(nuclear['G_thermal_eV'] if row['partition_kind']=='monatomic'
                else nuclear['nuclear_by_cutoff_eV'][cutoff]+nuclear['g_spin_partition_eV'])
            value=energy['e_smd_eV']+correction
            if not math.isfinite(value):raise ValueError('Nonfinite cutoff-specific G')
            values[key]=value
        expected=source['by_cutoff'][cutoff]
        for name,key in aliases.items():
            if abs(values[key]-number(expected['state_G_eV'][name]))>1e-8:raise ValueError('Cycle state differs')
        delta=math.fsum([values['polymer_pair'],values['empty'],-values['site2'],-values['site3']])*conversion
        if abs(delta-number(expected['deltaDeltaG_kJ_mol']))>1e-6:raise ValueError('Neighbor nonadditivity differs')
        models={}
        for mode in ('free_ions','contact_pair'):
            reagent=math.fsum([values['protonated'],values['chloride'],-values['precursor'],-values['base']]) if mode=='free_ions' else math.fsum([values['pair'],-values['precursor'],-values['base']])
            correction=0. if mode=='free_ions' else standard
            computed={}
            for site,other in ((2,3),(3,2)):
                computed[f'graft_{site}_on_empty_kJ_mol']=math.fsum([values[f'site{site}'],-values['empty'],reagent])*conversion+correction
                computed[f'graft_{site}_given_{other}_kJ_mol']=math.fsum([values['polymer_pair'],-values[f'site{other}'],reagent])*conversion+correction
            computed['direct_double_graft_kJ_mol']=math.fsum([values['polymer_pair'],-values['empty'],2*reagent])*conversion+2*correction
            for name,value in computed.items():
                if abs(value-number(expected['models'][mode][name]))>1e-6:raise ValueError('Conditional reaction differs')
            for site,other in ((2,3),(3,2)):
                closure=computed[f'graft_{site}_on_empty_kJ_mol']+computed[f'graft_{other}_given_{site}_kJ_mol']
                if abs(closure-computed['direct_double_graft_kJ_mol'])>1e-6:raise ValueError('Attachment-order cycle does not close')
            models[mode]=computed
        cutoffs[cutoff]=dict(deltaDeltaG_kJ_mol=delta,models=models)
    return dict(status='portable_exact_PMDI_neighbor_cycle_reconstructed',states=primary['states'],
        balanced_reactions=len(primary['reactions']),PMDI_first_reduction_V_vs_Fc=primary['PMDI_first_reduction_V_vs_Fc'],
        by_cutoff=cutoffs,quantum_calculations_run=False,nuclear_partition_reintegrated=False,
        conformer_convergence_established=False,synthesis_admission=False,
        isolated_motif_interaction_computed=False,electronic_coupling_computed=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--out',type=Path);args=parser.parse_args();result=replay(args.data)
    if args.out:
        with args.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
