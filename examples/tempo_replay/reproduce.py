#!/usr/bin/env python3
"""Reconstruct the exact local TEMPO oxidation potential from cached states."""
import argparse
import importlib.util
import json
import math
from pathlib import Path


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def replay(root,full_native=False):
    root=Path(root).resolve();manifest=json.loads((root/'manifest.json').read_text())
    expected={'tempo_rad':(0,2),'tempo_ox':(1,1),'fc_neutral':(0,1),'fc_oxidized':(1,2)}
    rows=manifest['states']
    if (len(rows)!=4 or {r['id'] for r in rows}!=set(expected)
            or any((r['charge'],r['mult'])!=expected[r['id']] for r in rows)
            or manifest['experimental_accuracy_inferred'] or manifest['material_eligibility_inferred']
            or manifest['grafting_deltaG_computed'] or manifest['conformer_convergence_established']):
        raise ValueError('Unexpected state identity or promoted scientific scope')
    reader=load_module(root/'cached_reader.py','tempo_cached_reader')
    primary=reader.reconstruct(root)
    source_path=(root/manifest['source_review']).resolve()
    if not source_path.is_relative_to(root):raise ValueError('Source review escapes package')
    source=json.loads(source_path.read_text())
    if (source['status']!='TEMPO_local_oxidation_native_partition_independently_reconstructed'
            or source['experimental_exact_graft_potential_claimed'] or source['free_energy_of_grafting_computed']):
        raise ValueError('TEMPO source scope differs')
    states={r['id']:r for r in rows};cutoffs={}
    aliases={'tempo_rad':'tempo__rad','tempo_ox':'tempo__ox','fc_neutral':'ferrocene__neu','fc_oxidized':'ferrocene__ox'}
    for cutoff in ('25','50','100'):
        g={}
        for key,row in states.items():
            def load(relative):
                p=(root/relative).resolve()
                if not p.is_relative_to(root):raise ValueError('State file escapes package')
                return json.loads(p.read_text())
            energy=load(row['files']['energy']);partition=load(row['files']['partition'])
            value=energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][cutoff]+partition['g_spin_partition_eV']
            expected_value=source['states'][aliases[key]]['G_by_cutoff_eV'][cutoff]
            if not math.isfinite(value) or abs(value-expected_value)>1e-8:
                raise ValueError('TEMPO/Fc state differs from checked source')
            g[key]=value
        value=math.fsum((g['tempo_ox'],-g['tempo_rad'],-g['fc_oxidized'],g['fc_neutral']))
        if abs(value-source['potential_V_vs_Fc_by_rotor_cutoff'][cutoff])>1e-9:
            raise ValueError('TEMPO potential differs from checked source')
        cutoffs[cutoff]=value
    if abs(primary['PMDI_first_reduction_V_vs_Fc']-cutoffs['50'])>1e-9:
        raise ValueError('Primary cached-state reader disagrees')
    result=dict(status='portable_TEMPO_local_oxidation_reconstructed',states=4,
        oxidation_V_vs_Fc_by_cutoff=cutoffs,quantum_calculations_run=False,
        nuclear_partition_reintegrated=False,experimental_accuracy_inferred=False,
        grafting_deltaG_computed=False,material_eligibility_inferred=False)
    if full_native:
        native=load_module(root/'native_engine.py','tempo_native_engine').native_replay(root,manifest)
        if native['Hessians_reconstructed']!=4 or native['monatomic_partitions_recomputed']!=0:
            raise ValueError('Expected four molecular Hessian reconstructions')
        native['source_algorithm_status']=native['status']
        native['status']='TEMPO_and_Fc_Hessians_modes_and_nuclear_sums_reconstructed'
        result['native_reconstruction']=native
        result['nuclear_partition_reintegrated']=True
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--out',type=Path);p.add_argument('--full-native',action='store_true');a=p.parse_args()
    result=replay(a.data,a.full_native)
    if a.out:
        with a.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
