"""Source-bound molecular native review with independent density/resource checks."""
import json
import math
from pathlib import Path,PurePosixPath
import shutil
import subprocess
import time

import numpy as np
from pyscf import gto

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.validation.solution_pair_review import density_spin
from redox.workflows.host_local_review import command_for_host
from redox.workflows.molecular_native_namespace import command,verify


class ReviewPending(RuntimeError):pass


def completed_review_record(path,timeout_seconds=30):
    """Allow bounded NFS visibility lag after the existing producer exits."""
    deadline=time.monotonic()+timeout_seconds
    while True:
        try:return json.loads(Path(path).read_text())
        except (FileNotFoundError,json.JSONDecodeError) as error:
            if time.monotonic()>=deadline:
                raise ReviewPending('Completed review handoff not yet visible; preserve it without launching a duplicate') from error
            time.sleep(.2)


def cross_device_uma_review(check,reference):
    """Recompute the reported errors from fixed-geometry energy/force arrays."""
    force=np.asarray(reference['forces_eV_A'],dtype=float)
    if (check['energy_tolerance_eV']!=1e-4 or check['force_tolerance_eV_A']!=1e-4
            or len(check['repeats'])!=3 or force.shape!=(len(reference['numbers']),3)
            or not np.isfinite(force).all()):raise ValueError('Incomplete fixed-geometry UMA control')
    errors=[]
    for row in check['repeats']:
        got=np.asarray(row['forces_eV_A'],dtype=float)
        if got.shape!=force.shape or not np.isfinite(got).all():raise ValueError('Malformed repeated UMA force array')
        energy=abs(row['energy_eV']-reference['energy_eV']);maximum=float(abs(got-force).max())
        if (not math.isfinite(energy) or energy>1e-4 or maximum>1e-4
                or abs(energy-row['energy_error_eV'])>1e-12 or abs(maximum-row['force_max_error_eV_A'])>1e-12):
            raise ValueError('Cross-device standard UMA energy/forces do not reconstruct or pass')
        errors.append(dict(energy_error_eV=energy,force_max_error_eV_A=maximum))
    return errors


def spin_record_agreement(spin,saved,mult):
    """Honor the producer's four decimals without weakening the physical gate."""
    S=(mult-1)/2;computed=spin['s_squared'];contamination=computed-S*(S+1)
    if not math.isfinite(computed) or abs(contamination)>.1:
        raise ValueError('Full-precision spin-contamination gate failed')
    half_unit=0.5e-4;roundoff=1e-10
    for actual,key in [(computed,'s_squared'),(contamination,'spin_contam')]:
        value=saved[key]
        if (type(value) not in (int,float) or not math.isfinite(value) or round(value,4)!=value
                or abs(actual-value)>half_unit+roundoff):
            raise ValueError('Independent density spin differs beyond the stored four-decimal precision')
    return dict(stored_decimal_places=4,full_precision_contamination=contamination,
                full_precision_limit=.1,stored_s_squared=saved['s_squared'])


class MolecularSources(ArchivedSources):
    def __init__(self,spec):
        super().__init__(spec['prior_root'],spec['archive'],relative_to_prior=True)
        prior=PurePosixPath(spec['prior_root']);case=Path(spec['case'])
        self.aliases=[(prior/'calcs/dft',case/'dft'),
            (prior/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs',case/'repairs')]

    def path(self,original):
        path=PurePosixPath(str(original))
        if '..' in path.parts:raise ValueError('Evidence paths cannot traverse parents')
        if not path.is_absolute():path=self.prior/path
        for prefix,target in self.aliases:
            try:return target.joinpath(*path.relative_to(prefix).parts)
            except ValueError:pass
        return super().path(path)


def review(spec_path,output,*,reaction_component=False):
    if output.exists():raise FileExistsError('Preserve existing molecular reviews')
    spec=verify(spec_path);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
    if reaction_component!=(plan.get('role')=='reaction_component'):
        raise ValueError('Candidate and reagent reviews must be requested separately')
    identity_key='inchikey' if reaction_component else 'product_inchikey'
    state=json.loads((case/'status.json').read_text())
    completed=('reaction_component_native_finished_independent_review_pending' if reaction_component
               else 'molecular_native_state_finished_independent_review_pending')
    if state['status']!=completed:raise ValueError('Completed native result required')
    guard=json.loads((case/'residency_guard/owner.json').read_text());residency=json.loads((case/'residency_guard/state.json').read_text())
    host=json.loads((case/'host_guard/status.json').read_text());runtime=json.loads((case/'tensor_runtime/status.json').read_text())
    if (not guard.get('guard_stopped') or not residency['own_residency_observed'] or residency['foreign_resident_PIDs']
            or (case/'residency_guard/interruption.json').exists()
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or not host['pid']==guard['pid']==runtime['pid']==state['owner']['pid']
            or not host['process_start']==guard['process_start']==state['owner']['process_start']
            or host['maximum_RSS_GiB']!=plan['maximum_host_RSS_GiB'] or host['minimum_free_GiB']!=32
            or not 0<max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])<=plan['maximum_host_RSS_GiB']*2**30
            or not host['CPU_affinity']==guard['CPU_affinity']==runtime['CPU_affinity'] or len(host['CPU_affinity'])!=2
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or runtime['program']!=str(Path(spec['sources'])/'worker.py')
            or not 0<state['observed_peak_GPU_MiB']<=plan['maximum_peak_GPU_MiB']):
        raise ValueError('Closed exact resource/runtime/residency chain failed')
    numerical_path=case/'numerical_review.json'
    if not numerical_path.exists():
        log=output.with_suffix('.log')
        if log.exists():raise ReviewPending('Existing numerical-review handoff retained; no duplicate launched')
        runtime_profile=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
        argv=command(spec)+['env','PYTHONPATH='+spec['prior_root']+'/src',runtime_profile['python'],
            str(Path(spec['sources'])/'reviewer.py'),'--plan',str(case/'plan.json'),'--out',str(numerical_path)]
        with log.open('x') as stream:
            try:subprocess.run(command_for_host(argv,state['owner']['host']),stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=300)
            except subprocess.TimeoutExpired as error:raise ReviewPending('Review observation expired; retain its existing handoff') from error
    numerical=completed_review_record(numerical_path);reader=MolecularSources(spec)
    reader.verify_all(numerical['source_sha256'])
    produced=json.loads((case/'result.json').read_text())
    if ((not reaction_component and sha256_file(case/'result.json')!=state['result_sha256'])
            or state.get('plan_sha256')!=sha256_file(case/'plan.json')
            or numerical[identity_key]!=plan[identity_key]
            or numerical['state']!=plan['state'] or numerical['charge']!=plan['charge'] or numerical['mult']!=plan['mult']
            or any(abs(v-produced['G_by_cutoff_eV'][c])>1e-8 for c,v in numerical['G_by_cutoff_eV'].items())):
        raise ValueError('Independent molecular result has different state/identity/energies')
    if reaction_component:
        post=json.loads((case/'post_UMA_runtime.json').read_text())
        if (produced.get('role')!='reaction_component' or numerical.get('role')!='reaction_component'
                or post['versions']!=runtime['versions'] or post['engine']!='cutensor'
                or post['library_sha256']!=runtime['library_sha256']
                or len(post['FP64_smoke_checks'])!=3
                or any(not math.isfinite(c['maximum_absolute_error']) or c['maximum_absolute_error']>1e-11
                       for c in post['FP64_smoke_checks'])):
            raise ValueError('Reagent role or post-UMA contraction-runtime gate failed')
    surface=json.loads(reader.path(produced['native_record_path']).read_text())
    matrix=np.load(reader.path(produced['native_record_path']).with_suffix('')/'hessian_cartesian.npy',allow_pickle=False)
    candidates=[]
    for path in (case/'SCF_evaluations').glob('*/*/hessian_tensor.npy'):
        h=np.load(path,allow_pickle=False)
        if h.ndim==4 and h.shape[1]==h.shape[0] and h.shape[2:]==(3,3):
            mat=h.transpose(0,2,1,3).reshape(matrix.shape)
            if np.max(np.abs(mat-matrix))<1e-12:candidates.append(path)
    if len(candidates)!=1:raise ValueError('One checkpointed base Hessian must bind the accepted native matrix')
    path=candidates[0];settings=json.loads(path.with_name('settings.json').read_text())
    mol=gto.M(atom=list(zip(settings['elements'],settings['coordinates_Bohr'])),unit='Bohr',
        charge=plan['charge'],spin=plan['mult']-1,basis=settings['basis'],verbose=0)
    if any(z>36 for z in mol.atom_charges()):raise ValueError('ECP density review needs an explicit additional binding')
    spin=density_spin(np.load(path.with_name('density.npy'),allow_pickle=False),mol.intor('int1e_ovlp'),(mol.nelectron,plan['mult']))
    spin_agreement=spin_record_agreement(spin,surface['smd_spin_diagnostics'],plan['mult'])
    snapshot=case/'review_sources'/str(time.time_ns());snapshot.mkdir(parents=True)
    copied=snapshot/'molecular_native_reviewer.py';shutil.copyfile(Path(__file__).resolve(),copied)
    density_source=snapshot/'density_spin_reviewer.py'
    shutil.copyfile(ROOT/'src/redox/validation/solution_pair_review.py',density_source)
    sources=[case/'plan.json',case/'status.json',case/'result.json',numerical_path,
        case/'residency_guard/owner.json',case/'residency_guard/state.json',case/'host_guard/status.json',
        case/'tensor_runtime/status.json',path,path.with_name('density.npy'),path.with_name('settings.json'),
        copied,density_source]
    reader.recheck()
    result=dict(status=('reaction_component_native_G_independently_checked' if reaction_component
                       else 'exact_molecular_native_G_independently_checked'),
        source_id=plan['source_id'],state=plan['state'],charge=plan['charge'],mult=plan['mult'],
        G_by_cutoff_eV=numerical['G_by_cutoff_eV'],
        density_spin=spin,spin_record_agreement=spin_agreement,resource_pilot_passed_for_this_state=True,
        repaired_geometry_used=(produced['final_geometry_sha256']!=produced['initial_optimized_geometry_sha256']
            if reaction_component else numerical['repaired_geometry_used']),namespace=str(spec_path),
        namespace_sha256=sha256_file(spec_path),global_conformer_convergence_proven=False,
        material_admission=False,experimental_accuracy_inferred=False,declared_spin_model=True,
        source_sha256={**{str(reader.path(p)):h for p,h in numerical['source_sha256'].items()},
                       **{str(p):sha256_file(p) for p in sources}})
    result[identity_key]=plan[identity_key]
    if reaction_component:
        result.update(role='reaction_component',route_candidate=plan['route_candidate'])
        for path in (case/'post_UMA_runtime.json',case/'uma/result.json'):
            result['source_sha256'][str(path)]=sha256_file(path)
        if plan.get('UMA_cross_device_reference'):
            path=case/'UMA_cross_device_check.json';cross=json.loads(path.read_text())
            reference=json.loads(Path(plan['UMA_cross_device_reference']).read_text())
            if cross['reference_sha256']!=sha256_file(plan['UMA_cross_device_reference']):
                raise ValueError('Cross-device standard UMA reference checks failed')
            result['UMA_cross_device_errors']=cross_device_uma_review(cross,reference)
            result['source_sha256'][str(path)]=sha256_file(path)
            for origin,digest in reference['source_sha256'].items():
                if sha256_file(origin)!=digest:raise ValueError('Original V100 fixed-geometry reference changed')
                result['source_sha256'][origin]=digest
    else:result['candidate_alias']=plan['candidate_alias']
    atomic_write_json(output,result);return result
