"""Independently collect a completed source-bound molecular solution pair."""
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.workflows.host_local_review import command_for_host
from redox.workflows.solution_pair_namespace import verify,command,geometry_pair


class ReviewPending(RuntimeError):
    """An existing audit handoff cannot yet be treated as terminal."""


def density_spin(dm,overlap,charge_spin):
    """Recover electron/spin traces directly from the archived AO density."""
    electrons,multiplicity=charge_spin
    if (overlap.ndim!=2 or overlap.shape[0]!=overlap.shape[1]
            or not np.isfinite(dm).all() or not np.isfinite(overlap).all()):raise ValueError('Nonfinite density or overlap')
    if dm.ndim==2:
        if multiplicity!=1 or dm.shape!=overlap.shape:raise ValueError('Restricted density has the wrong shape/multiplicity')
        alpha=beta=dm/2
    elif dm.shape==(2,*overlap.shape):alpha,beta=dm
    else:raise ValueError('Density shape differs from the declared spin state')
    if abs(dm-np.swapaxes(dm,-1,-2)).max()>1e-8:raise ValueError('Density is not Hermitian')
    errors=[float(np.linalg.norm(value@overlap@value-value)/max(1.,np.linalg.norm(value))) for value in (alpha,beta)]
    if not np.isfinite(errors).all() or max(errors)>1e-6:raise ValueError('Density occupations are not idempotent')
    traces=[float(np.einsum('ij,ji',value,overlap)) for value in (alpha,beta)]
    if abs(sum(traces)-electrons)>1e-6 or abs(traces[0]-traces[1]-(multiplicity-1))>1e-6:
        raise ValueError('Density electron or spin populations differ')
    s2=(traces[0]-traces[1])**2/4+sum(traces)/2-float(np.trace(alpha@overlap@beta@overlap))
    spin=(multiplicity-1)/2
    if not math.isfinite(s2) or abs(s2-spin*(spin+1))>.1:raise ValueError('Reconstructed spin contamination fails')
    return dict(alpha_electrons=traces[0],beta_electrons=traces[1],s_squared=s2,
                idempotency_errors=errors)


def review(spec_path,output):
    if output.exists():raise FileExistsError('Preserve earlier solution-pair reviews')
    spec=verify(spec_path);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
    stage=json.loads((case/'status.json').read_text())
    if stage['status']!='solution_pair_finished_independent_review_pending':raise ValueError('Completed numerical result required')
    reader=ArchivedSources(spec['prior_root'],spec['archive'],relative_to_prior=True)
    guard=json.loads((case/'residency_guard/owner.json').read_text());residency=json.loads((case/'residency_guard/state.json').read_text())
    host=json.loads((case/'host_guard/status.json').read_text());runtime=json.loads((case/'tensor_runtime/status.json').read_text())
    if (not guard.get('guard_stopped') or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
            or (case/'residency_guard/interruption.json').exists()
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or host['maximum_RSS_GiB']!=plan['maximum_host_RSS_GiB'] or host['minimum_free_GiB']!=32
            or not host['pid']==guard['pid']==runtime['pid']==stage['owner']['pid']
            or not host['process_start']==guard['process_start']==stage['owner']['process_start']
            or not host['CPU_affinity']==guard['CPU_affinity']==runtime['CPU_affinity'] or len(host['CPU_affinity'])!=2
            or not 0<max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])<=plan['maximum_host_RSS_GiB']*2**30
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or runtime['program']!=str(Path(spec['sources'])/'worker.py')
            or not 0<stage['observed_peak_GPU_MiB']<=31744):
        raise ValueError('Closed execution/resource evidence did not pass')
    control=json.loads((case/'backend_control/native_smd_gpu_reference.json').read_text())
    if (control['status']!='passed' or not math.isfinite(control['residual_eV']) or abs(control['residual_eV'])>1e-5
            or not all(math.isfinite(control[k]) for k in ('wrapper_Ha','independent_hcore_Ha'))
            or abs((control['wrapper_Ha']-control['independent_hcore_Ha'])*27.211386245988-control['residual_eV'])>1e-8):
        raise ValueError('Fresh backend reference failed')
    result=json.loads((case/'pair_results/result.json').read_text())
    if result['identity']['pair']!=plan['pair'] or sha256_file(case/'pair_results/result.json')!=stage['result_sha256']:
        raise ValueError('The result belongs to a different exact pair')
    geometries=[reader.path(s['geometry']) for s in plan['pair']['states']]
    diagnostic=geometry_pair(geometries,plan['product_smiles'])
    runtime_profile=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    audits=case/'independent_checks'
    if not (audits/'solution_pair_audit.json').exists():
        if output.with_suffix('.log').exists():
            raise ReviewPending('Existing audit handoff preserved; no duplicate reviewer launched')
        argv=command(spec)+['env','PYTHONPATH='+spec['prior_root']+'/src',runtime_profile['python'],
            spec['prior_root']+'/scripts/validation/audit_solution_reorg_pair.py','--out',str(case/'pair_results'),
            '--report-prefix','solution_pair','--report-dir',str(audits),'--pair-only']
        argv=command_for_host(argv,stage['owner']['host'])
        with output.with_suffix('.log').open('x') as log:
            try:subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
            except subprocess.TimeoutExpired as error:
                raise ReviewPending('Audit observation expired; inspect its existing handoff before restart') from error
    audit=json.loads((audits/'solution_pair_audit.json').read_text())
    if audit['failed_checks'] or audit['fresh_SCF_records']!=6:raise ValueError('Six fresh independent SCF/density checks required')
    native=json.loads((Path(spec['sources'])/'native_review.json').read_text());energy_differences={};spin_checks=[]
    states={s['name']:s for s in plan['pair']['states']}
    for path in (case/'pair_results').glob('scf_*.json'):
        item=json.loads(path.read_text());key=item['key'];state=states[key['state']['name']]
        xyz=reader.path(next(s['geometry'] for s in states.values() if sha256_file(reader.path(s['geometry']))==key['geometry_sha256']))
        atoms=read(xyz);mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),
            charge=state['charge'],spin=state['mult']-1,basis=result['identity']['method']['basis'],verbose=0)
        with np.load(path.with_suffix('.npz'),allow_pickle=False) as saved:dm=saved['dm']
        spin=density_spin(dm,mol.intor('int1e_ovlp'),(mol.nelectron,state['mult']))
        if state['mult']>1 and (not math.isfinite(item['s2']) or abs(spin['s_squared']-item['s2'])>1e-6):
            raise ValueError('Saved and reconstructed spin differ')
        spin_checks.append(dict(record=str(path),**spin))
        if key['reference'] is None and key['geometry_sha256']==sha256_file(reader.path(state['geometry'])):
            reference=native['states']['tempo__'+state['name']]['active_E_smd_eV']
            delta=item['e_eV']-reference
            if not math.isfinite(delta) or abs(delta)>1e-4:raise ValueError('Own-geometry SP does not reproduce the active source')
            energy_differences[state['name']]=delta
    if set(energy_differences)!=set(states):raise ValueError('Both own-geometry source energies must be checked')
    source_files=[case/'plan.json',case/'status.json',case/'pair_results/result.json',case/'host_guard/status.json',
        case/'residency_guard/owner.json',case/'residency_guard/state.json',case/'tensor_runtime/status.json',
        case/'backend_control/native_smd_gpu_reference.json',audits/'solution_pair_audit.json',audits/'solution_pair_checks.json',
        Path(__file__).resolve(),ROOT/'src/redox/workflows/solution_pair_namespace.py',
        ROOT/'src/redox/validation/protonation.py',ROOT/'src/redox/build/geometry.py',
        ROOT/'src/redox/properties/solution_reorg.py',
        *list((case/'pair_results').glob('scf_*.json')),*list((case/'pair_results').glob('scf_*.npz'))]
    from redox.properties.solution_reorg import summarize_relaxations
    summary=summarize_relaxations(*result['directional'],endpoint_evidence=result['endpoint_evidence'],
        method_identity=result['identity']['method'],independent_validation=True,implementation_validation=True)
    summary.update(qualification_convention=result['summary']['qualification_convention'],
        active_SP_stationary_endpoints_proven=False)
    extra_flags=sorted({flag for row in result['electronic_state_evidence'] for flag in row['flags']})
    if extra_flags:summary['physical_flags'].extend(extra_flags);summary['endpoint_model_verified']=False
    reader.verify_all(plan['source_sha256']);reader.recheck()
    result_review=dict(status='exact_TEMPO_directional_solution_descriptor_independently_checked',
        candidate_id='tempo',product_smiles=plan['product_smiles'],product_inchikey=plan['product_inchikey'],
        source_sha256={**{str(reader.path(p)):h for p,h in plan['source_sha256'].items()},
                       **{str(p):sha256_file(p) for p in source_files}},namespace=str(spec_path),namespace_sha256=sha256_file(spec_path),
        method=result['identity']['method'],directional=result['directional'],summary=summary,source_result_summary=result['summary'],
        geometry_diagnostic=diagnostic,spin_checks=spin_checks,own_geometry_energy_differences_eV=energy_differences,
        endpoint_evidence=result['endpoint_evidence'],independent_checks=audit['checks'],fresh_SCF_records=6,
        resource_pilot_passed_for_this_case=True,solution_rate_validated=False,global_conformer_convergence_proven=False,
        grafting_free_energy_computed=False,material_admission=False,gas_inner_or_Born_outer_terms_added=False)
    atomic_write_json(output,result_review);return result_review
