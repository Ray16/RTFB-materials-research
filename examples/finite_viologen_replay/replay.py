#!/usr/bin/env python3
"""Standard-library replay of the paper's finite viologen comparisons."""
import ast
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
EV_KJ=Decimal('1.602176634e-19')*Decimal('6.02214076e23')/Decimal(1000)
HA_EV=Decimal('27.211386245988')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def number(value):
    if type(value) not in (int,float) or not math.isfinite(value):raise ValueError('Finite component required')
    return Decimal(str(value))
def close(actual,expected,tolerance=1e-7):
    if abs(float(actual)-float(number(expected)))>tolerance:raise ValueError('Cached arithmetic differs')

def replay(root=ROOT):
    manifest=json.loads((root/'manifest.json').read_text())
    for name,digest in manifest['files_sha256'].items():
        path=(root/name).resolve()
        if not path.is_relative_to(root.resolve()) or sha(path)!=digest:raise ValueError('Packaged input changed')
    folder=root/'data';data=json.loads((folder/'ionic_native_cycles_v4.json').read_text())
    preparation=json.loads((folder/'ionic_native_cycle_preparation.json').read_text())
    cases={r['case']['key']:r['case'] for r in preparation['cases']};records=data['records']
    if (len(records)!=12 or len(data['cycles'])!=4 or not data['all_declared_native_cycles_qualified']
            or data['conformer_ensemble_qualified'] or data['primary_model_admission'] or data['material_admission']):
        raise ValueError('The exact finite twelve-component/four-cycle scope is required')
    common=None
    for key,r in records.items():
        context=(r['model'],r['quantum_settings'])
        if common is not None and context!=common:raise ValueError('Mixed electronic/solvent contexts')
        common=context
        if (not r['native_G_qualified'] or r['mult']!=1 or r['charge']!=cases[key]['charge']
                or r['atoms']!=cases[key]['atoms'] or number(r['minimum_frequency_cm']) < -5
                or not 0<=number(r['maximum_column_error'])<=Decimal('.01')):raise ValueError('Cached component state or QC differs')
        for cutoff in ('25','50','100'):
            close(number(r['E_geometry_level_Ha'])*HA_EV+number(r['nuclear_by_cutoff_eV'][cutoff]),r['G_by_cutoff_eV'][cutoff])
    cycles=[];common_thermal=None;common_reference=None
    for cycle in data['cycles']:
        definition=cycle['definition'];configuration=definition['configuration'];i,j=definition['sites']
        expected={f'{configuration}__occ_{i}_{j}':1,f'{configuration}__occ_none':1,
                  f'{configuration}__occ_{i}':-1,f'{configuration}__occ_{j}':-1}
        if definition['coefficients']!=expected or not cycle['native_cycle_qualified'] or cycle['missing']:
            raise ValueError('Exact complete balanced occupancy cycle required')
        for field in ('charge','atoms'):
            if sum(records[k][field]*c for k,c in expected.items())!=0:raise ValueError('Charge or atom-count imbalance')
        energy=sum((number(records[k]['E_geometry_level_Ha'])*c for k,c in expected.items()),Decimal(0))*HA_EV*EV_KJ
        close(energy,cycle['delta_delta_E_kJ_mol']);free={}
        for cutoff in ('25','50','100'):
            assembly=cycle['property_assemblies'][cutoff];ctx=assembly['context']
            if (not assembly['context_verified'] or ctx['unit']!='eV' or ctx['temperature_K']!=298.15
                    or ctx['electronic_protocol']!=common[1] or ctx['solvent_model']!=common[0]
                    or ctx['thermal_definition']['selected_rotor_cutoff_cm']!=int(cutoff)
                    or ctx['reference']['explicit_spectator_ions']):raise ValueError('Thermal/reference context differs')
            thermal=ctx['thermal_definition']['protocol']
            if common_thermal is not None and thermal!=common_thermal:raise ValueError('Mixed thermal protocols')
            common_thermal=thermal
            reference=(ctx['reference'],ctx['standard_state'])
            if common_reference is not None and reference!=common_reference:raise ValueError('Mixed reference or standard state')
            common_reference=reference
            value=sum((number(records[k]['G_by_cutoff_eV'][cutoff])*c for k,c in expected.items()),Decimal(0))*EV_KJ
            close(value,cycle['delta_delta_G_by_cutoff_kJ_mol'][cutoff]);free[cutoff]=float(value)
        cycles.append(dict(configuration=configuration,sites=[i,j],E_kJ_mol=float(energy),G_by_cutoff_kJ_mol=free))
    # Preserve the existing independent motif checker verbatim. Replace its sole
    # project dependency (file hashing) with the identical standard-library adapter.
    source=root/'motif_checker_source.py';tree=ast.parse(source.read_text());removed=[];body=[]
    for node in tree.body:
        if isinstance(node,ast.ImportFrom) and node.module and node.module.startswith('redox'):
            if node.module!='redox.core.hashing' or [a.name for a in node.names]!=['sha256_file']:
                raise ValueError('Unexpected project dependency in frozen checker')
            removed.append(node)
        else:body.append(node)
    if len(removed)!=1:raise ValueError('Exactly one declared hashing adapter required')
    scope={'sha256_file':sha,'__file__':str(source)}
    exec(compile(ast.Module(body=body,type_ignores=[]),str(source),'exec'),scope)
    motifs,frames=scope['reconstruct'](folder)
    for row in motifs:
        config,sites=row['key'].split('__occ_');matched=next(c for c in cycles if c['configuration']==config and c['sites']==list(map(int,sites.split('_'))))
        close(matched['E_kJ_mol'],row['whole_E'])
    return dict(status='packaged_finite_cycle_and_motif_arithmetic_replayed',components=12,cycles=cycles,
        motif_triplets=12,previously_reviewed_SCF_frames_represented=frames,motif_comparisons=motifs,
        opposite_sign_frames=sum((r['whole_E']<0)!=(r['motif_E']<0) for r in motifs),
        context=dict(electronic=common[1],solvent_model=common[0],thermal=common_thermal,
                     spectator_ions='omitted separated spectators',temperature_K=298.15),
        raw_quantum_QC_repeated=False,SCF_or_Hessian_calculations_run=False,
        conformer_ensemble_qualified=False,experimental_accuracy_validated=False,material_admission=False)

if __name__=='__main__':print(json.dumps(replay(),indent=2))
