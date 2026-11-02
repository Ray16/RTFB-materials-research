#!/usr/bin/env python3
"""Review six actual SCFs before assembling a charged capped-motif interaction."""
import argparse
from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path

from redox.core.common import ROOT, HARTREE_EV, EV_KJ
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.properties.context import combine_properties
from redox.validation.execution_review import review_closure, freeze_source, SOURCE_BYTES as CLOSURE_BYTES
from redox.validation import primary_checkpoints
from redox.validation.scf_frame import freeze_reviewer
from redox.workflows.molecular_native_namespace import verify, load
from redox.workflows.ionic_stationary_cycles import screening_signature, validate_cycle_balance

SOURCE_BYTES=Path(__file__).read_bytes()


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read(path):used.append(Path(path));return json.loads(Path(path).read_text())
    plan=read(case/'plan.json')
    if plan['role']!='frozen_charged_motif_energy' or plan['native_G_qualified'] or plan['geometry_optimized']:
        raise ValueError('Exact frozen charged-motif E-only scope required')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','motif_energy_exact_owner')
    stage=review_closure(case,read,owners.alive,
        expected_status='frozen_charged_motif_energy_finished_independent_review_pending',host_GiB=16,
        GPU_MiB=31744,worker=sources/'worker.py')
    used.extend([freeze_reviewer(case),freeze_source(case,Path(__file__),SOURCE_BYTES),
        freeze_source(case,ROOT/'src/redox/validation/execution_review.py',CLOSURE_BYTES),
        freeze_source(case,Path(primary_checkpoints.__file__),primary_checkpoints.SOURCE_BYTES),
        freeze_source(case,ROOT/'src/redox/properties/context.py')])
    construction=read(sources/'construction_review.json')
    if construction['status']!='fixed_parent_local_motif_construction_independently_checked' or any(sha256_file(p)!=h for p,h in construction['source_sha256'].items()):
        raise ValueError('The exact full-precision caps/pose construction audit changed')
    result=read(case/'result.json');model=read(case/'model/model.json');actual_model=model['profile']
    if (model['status']!='experimental_smooth_CDS_context_finished_QC_separate'
            or set(actual_model)!=set(plan['solvent_model'])|{'implementation_sha256'}
            or any(actual_model[k]!=v for k,v in plan['solvent_model'].items())
            or result['method']!=plan['method'] or result['solvent_model']!=plan['solvent_model']):
        raise ValueError('Actual gas/MeCN implementation or declared physics differs')
    from redox.validation.protonation import hydrogen_topology_status
    tasks={t['id']:t for t in plan['tasks']};records={medium:{} for medium in plan['media']};counts=[]
    validate_cycle_balance(dict(coefficients={'AB':1,'A':-1,'B':-1}),tasks)
    for label,task in tasks.items():
        geometry=Path(task['geometry']);used.append(geometry)
        if sha256_file(geometry)!=task['geometry_sha256'] or hydrogen_topology_status(task['smiles'],geometry)!='ok':
            raise ValueError('Exact frozen atom frame or full hydrogen graph differs')
        for medium in plan['media']:
            expected=dict(plan['method']);expected.pop('density_fit')
            expected.update(charge=task['charge'],mult=1,solvent_class='SMD' if medium=='acetonitrile' else None)
            series=primary_checkpoints.checkpoint_series(case/'SCF_evaluations'/label/medium,geometry,expected,
                read,used,expected_count=1,fixed_geometry=True)
            value=result['energies'][label][medium]
            if (not value['SCF_converged'] or not math.isfinite(value['E_Ha'])
                    or abs(value['E_Ha']-series['energies_Ha'][0])>1e-10
                    or (value['charge'],value['mult'])!=(task['charge'],1)
                    or value['geometry_sha256']!=task['geometry_sha256']):
                raise ValueError('Actual frozen motif energy/state does not reconstruct')
            context=dict(quantity='electronic_energy',unit='Ha',electronic_protocol=plan['method'],
                thermal_definition='none: electronic energies only',solvent=medium,
                solvent_model=plan['solvent_model'] if medium=='acetonitrile' else 'vacuum',
                temperature_K=298.15,standard_state='not applied to fixed-geometry electronic energies',
                reference='absolute electronic components at the same frozen relative pose and declared medium')
            records[medium][label]=dict(value=value['E_Ha'],context=context,
                state=dict(identity=task['smiles'],charge=task['charge'],multiplicity=1),
                qc='independently_checked',sources={str(geometry):sha256_file(geometry),str(case/'result.json'):sha256_file(case/'result.json')})
            counts.append(series['actual_SCF_frames_checked'])
    sums={medium:combine_properties(rows,{'AB':1,'A':-1,'B':-1}) for medium,rows in records.items()}
    values={medium:row['value']*HARTREE_EV*EV_KJ for medium,row in sums.items()}
    implementation=case/'screening_implementation/implementation.json';kernel=implementation.with_name('screening_kernel.cu')
    signature=screening_signature(read(implementation),sources,kernel);used.append(kernel)
    host=read(case/'host_guard/status.json');verify(namespace)
    return dict(status='frozen_charged_motif_triplet_independently_checked',namespace=str(namespace),key=plan['source_id'],
        parent_key=plan['parent_key'],cap_length_A=plan['cap_length_A'],actual_SCF_frames_checked=sum(counts),
        method=plan['method'],solvent_model=plan['solvent_model'],numerical_implementation=signature,
        exact_energies_Ha=result['energies'],interaction_E_kJ_mol_by_medium=values,
        solvent_response_to_interaction_kJ_mol=values['acetonitrile']-values['gas'],checked_sums=sums,
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'],peak_host_RSS_bytes=max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes']),
        native_G_qualified=False,geometry_optimized=False,basis_superposition_corrected=False,
        electronic_coupling_computed=False,deformation_energy_computed=False,material_admission=False,
        scope='Frozen methyl-capped local motif E interaction including cap and continuum response; not binding G, pure pi energy, coupling, occupancy population or complete backbone decomposition',
        source_sha256={str(p):sha256_file(p) for p in used},producer_sha256=hashlib.sha256(SOURCE_BYTES).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier independent triplet reviews')
    result=review(a.spec.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','key','actual_SCF_frames_checked','interaction_E_kJ_mol_by_medium')}))
