#!/usr/bin/env python
"""Review gas counterpoise control states, grids and CPU/GPU parity."""
import argparse
import ast
import json
import math
import os
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def compare(cpu,gpu):
    if (cpu['settings']!=gpu['settings'] or cpu['common_grid_hashes']!=gpu['common_grid_hashes']
            or cpu['source_symbols']!=gpu['source_symbols'] or cpu['source_coordinates_A']!=gpu['source_coordinates_A']
            or cpu['A_indices']!=gpu['A_indices'] or set(cpu['terms'])!=set(gpu['terms'])):
        raise ValueError('Counterpoise controls differ in protocol, geometry, partition or common grid')
    errors={role:abs(cpu['terms'][role]['E_Ha']-gpu['terms'][role]['E_Ha']) for role in cpu['terms']}
    if not all(math.isfinite(v) and v<=1e-6 for v in errors.values()):raise ValueError('CPU/GPU counterpoise energy parity failed')
    error=abs(cpu['sums']['counterpoise_gas_kJ_mol']-gpu['sums']['counterpoise_gas_kJ_mol'])
    if not math.isfinite(error) or error>.01:raise ValueError('CPU/GPU corrected interaction differs')
    return dict(maximum_energy_difference_Ha=max(errors.values()),per_role_energy_differences_Ha=errors,
        counterpoise_interaction_difference_kJ_mol=error)


def review(cpu_root,campaign,out):
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    state_path=campaign/'counterpoise/status.json';state=json.loads(state_path.read_text())
    owners=api.load_module('counterpoise_stopped_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if state['status']!='gas_counterpoise_GPU_complete_CPU_comparison_and_review_pending' or owners.owner_alive(state['owner']):
        raise ValueError('Completed stopped GPU counterpoise control required')
    api.verify_sources(state['output_sha256']);gpu_root=Path(state['raw_root'])/'GPU'
    core=api.load_module('counterpoise_independent_state_join',ROOT/'scripts/validation/local_counterpoise.py')
    documents=[json.loads((root/'summary.json').read_text()) for root in (cpu_root,gpu_root)]
    provenance_path=cpu_root.with_name(cpu_root.name+'_provenance.json')
    provenance=json.loads(provenance_path.read_text());api.verify_sources(provenance['source_sha256'])
    if provenance['status']!='fixed_RI_CPU_reference_computed_review_pending' or owners.owner_alive(provenance['owner']):
        raise ValueError('Completed stopped fixed-RI CPU reference required')
    parity=compare(*documents);proof={**plan['producer_sha256'],**plan['input_sha256'],**state['output_sha256']}
    from pyscf.lib import chkfile
    for root,doc in zip((cpu_root,gpu_root),documents):
        expected=core.molecules(doc['source_symbols'],doc['source_coordinates_A'],doc['A_indices'])
        if doc['settings']!=plan['settings'] or set(doc['terms'])!=set(core.ROLES):raise ValueError('Wrong five-term counterpoise protocol')
        if doc['sums']!=core.sums({k:v['E_Ha'] for k,v in doc['terms'].items()}):raise ValueError('Counterpoise arithmetic does not reconstruct')
        aux_sizes=set();aux_bases=[]
        for role,point in doc['terms'].items():
            mol=chkfile.load_mol(str(root/role/'SCF.h5'));scf=chkfile.load(str(root/role/'SCF.h5'),'scf');target=expected[role]
            if (not point['converged'] or point['settings']!=doc['settings']
                    or list(mol.atom_charges())!=list(target.atom_charges())
                    or not np.allclose(mol.atom_coords(),target.atom_coords(),rtol=0,atol=1e-11)
                    or mol.nelectron!=target.nelectron or mol.nao_nr()!=target.nao_nr()
                    or not math.isfinite(float(scf['e_tot'])) or abs(float(scf['e_tot'])-point['E_Ha'])>1e-10
                    or not math.isfinite(float(np.sum(scf['mo_occ']))) or abs(float(np.sum(scf['mo_occ']))-target.nelectron)>1e-10
                    or point['common_grid_hashes']!=doc['common_grid_hashes']):
                raise ValueError('Counterpoise checkpoint nucleus, electron, basis or energy binding failed')
            aux_sizes.add(point['auxiliary_functions']);aux_bases.append(ast.literal_eval(point['auxiliary_basis']))
        if len(aux_sizes)!=1 or min(aux_sizes)<=0 or any(v!=aux_bases[0] for v in aux_bases):raise ValueError('Common RI space changed between counterpoise terms')
        for own,ghost in [('A_own','A_ghost'),('B_own','B_ghost')]:
            if doc['terms'][ghost]['E_Ha']>doc['terms'][own]['E_Ha']+1e-7:raise ValueError('Expanded-basis monomer energy did not satisfy its variational check')
        for file in root.rglob('*'):
            if file.is_file() and file.suffix in ('.json','.npy','.h5'):proof[str(file)]=api.digest(file)
    for name,filename in [('coordinates','grid_coordinates_AU.npy'),('weights','grid_weights.npy')]:
        if not np.array_equal(np.load(cpu_root/filename),np.load(gpu_root/filename)):raise ValueError('Actual common grid arrays differ')
        for root,doc in zip((cpu_root,gpu_root),documents):
            if api.digest(root/filename)!=doc['common_grid_hashes'][filename]:raise ValueError('Stored common grid bytes differ')
    guard=campaign/'residency_guard';owner=json.loads((guard/'owner.json').read_text());resident=json.loads((guard/'state.json').read_text())
    host_path=campaign/'host_guard/status.json';host=json.loads(host_path.read_text());runtime_path=campaign/'tensor_runtime/status.json';runtime=json.loads(runtime_path.read_text());profile=json.loads(Path(plan['runtime_profile']).read_text())
    if (not owner['guard_stopped'] or resident['foreign_resident_PIDs'] or not resident['own_residency_observed']
            or (guard/'interruption.json').exists() or host['status']!='program_returned_review_required'
            or not 0<host['peak_RSS_bytes']<=32*2**30 or not 0<state['observed_peak_GPU_MiB']<=14336
            or runtime['status']!='task_returned_under_checked_cutensor_runtime' or runtime['program']!=plan['execution_driver']
            or runtime['profile_sha256']!=api.digest(Path(plan['runtime_profile'])) or runtime['versions']!=profile['versions']):
        raise ValueError('Counterpoise runtime/resource gate failed')
    proof.update(provenance['source_sha256'])
    for p in (Path(__file__).resolve(),provenance_path,campaign/'plan.json',state_path,guard/'owner.json',guard/'state.json',host_path,runtime_path):proof[str(p)]=api.digest(p)
    result=dict(status='gas_counterpoise_CPU_GPU_control_independently_checked',checked=time.time(),settings=core.SETTINGS,
        roles=list(core.ROLES),CPU=documents[0]['sums'],GPU=documents[1]['sums'],parity=parity,
        exact_ghost_nuclei_electrons_and_checkpoint_binding=True,common_grid_and_auxiliary_space=True,
        observed_peak_GPU_MiB=state['observed_peak_GPU_MiB'],peak_host_RSS_bytes=host['peak_RSS_bytes'],source_sha256=proof,
        solvent_counterpoise_qualified=False,thermal_free_energy_computed=False,primary_candidate_admission=False,
        scope='Gas five-term numerical counterpoise validation; no solvent correction or experimental interaction accuracy claim')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cpu',type=Path,required=True);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    row=review(a.cpu.resolve(),a.campaign.resolve(),a.out.resolve());print(json.dumps(row['parity']))


if __name__=='__main__':main()
