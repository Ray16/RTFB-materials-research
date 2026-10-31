"""Source-bound molecular native review with independent density/resource checks."""
import json
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


def review(spec_path,output):
    if output.exists():raise FileExistsError('Preserve existing molecular reviews')
    spec=verify(spec_path);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
    state=json.loads((case/'status.json').read_text())
    if state['status']!='molecular_native_state_finished_independent_review_pending':raise ValueError('Completed native result required')
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
    numerical=json.loads(numerical_path.read_text());reader=MolecularSources(spec)
    reader.verify_all(numerical['source_sha256'])
    produced=json.loads((case/'result.json').read_text())
    if (sha256_file(case/'result.json')!=state['result_sha256'] or numerical['product_inchikey']!=plan['product_inchikey']
            or numerical['state']!=plan['state'] or numerical['charge']!=plan['charge'] or numerical['mult']!=plan['mult']
            or any(abs(v-produced['G_by_cutoff_eV'][c])>1e-8 for c,v in numerical['G_by_cutoff_eV'].items())):
        raise ValueError('Independent molecular result has different state/identity/energies')
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
    if abs(spin['s_squared']-surface['smd_spin_diagnostics']['s_squared'])>1e-6:
        raise ValueError('Independent baseline density spin does not reconstruct')
    snapshot=case/'review_sources'/str(time.time_ns());snapshot.mkdir(parents=True)
    copied=snapshot/'molecular_native_reviewer.py';shutil.copyfile(Path(__file__).resolve(),copied)
    sources=[case/'plan.json',case/'status.json',case/'result.json',numerical_path,
        case/'residency_guard/owner.json',case/'residency_guard/state.json',case/'host_guard/status.json',
        case/'tensor_runtime/status.json',path,path.with_name('density.npy'),path.with_name('settings.json'),
        copied,ROOT/'src/redox/validation/solution_pair_review.py']
    reader.recheck()
    result=dict(status='exact_molecular_native_G_independently_checked',candidate_alias=plan['candidate_alias'],
        source_id=plan['source_id'],state=plan['state'],charge=plan['charge'],mult=plan['mult'],
        product_inchikey=plan['product_inchikey'],G_by_cutoff_eV=numerical['G_by_cutoff_eV'],
        density_spin=spin,resource_pilot_passed_for_this_state=True,
        repaired_geometry_used=numerical['repaired_geometry_used'],namespace=str(spec_path),
        namespace_sha256=sha256_file(spec_path),global_conformer_convergence_proven=False,
        material_admission=False,experimental_accuracy_inferred=False,declared_spin_model=True,
        source_sha256={**{str(reader.path(p)):h for p,h in numerical['source_sha256'].items()},
                       **{str(p):sha256_file(p) for p in sources}})
    atomic_write_json(output,result);return result
