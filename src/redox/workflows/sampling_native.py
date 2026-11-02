"""Prepare corrected native thermochemistry from an exact checked energy input.

Uses the original smooth-model native driver and its independent reconstruction.
The checked largest neutral Hessian supplies a resource witness, not a G value.
Every new input must pass its own derivative, mode and partition checks.
"""
import argparse
import ast
import json
import os
from pathlib import Path
import shutil

from ase.io import read
from pyscf import gto
from pyscf.lib import chkfile

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify, load

PILOT_KEYS={'all_same__occ_none__conf001','all_same__occ_2__conf001','all_same__occ_2_3__conf001'}


def profile():
    return dict(maximum_host_RSS_GiB=80,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
        minimum_GPU_memory_MiB=31*1024,maximum_wall_seconds=28800,
        maximum_durable_output_bytes=2*2**30,minimum_initial_storage_GiB=34,
        minimum_live_storage_GiB=32,response_available_memory_GiB=4)


def native_settings(path):
    tree=ast.parse(path.read_text())
    node=next(n.value for n in tree.body if isinstance(n,ast.Assign)
              and any(isinstance(t,ast.Name) and t.id=='NATIVE_SETTINGS' for t in n.targets))
    return {item.arg:ast.literal_eval(item.value) for item in node.keywords}


def resource_witness(review_path, model):
    from redox.validation.execution_review import review_closure
    report=json.loads(review_path.read_text())
    if (report['status']!='corrected_full_polymer_Hessian_independently_reconstructed'
            or report['atoms']!=119 or not report['equilibrium_solvent_H1_retained']
            or not report['complete_saved_Hessian_validated'] or report['maximum_column_error']>.01
            or report['minimum_frequency_cm'] < -5
            or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
        raise ValueError('Complete independently reconstructed corrected neutral Hessian required')
    checked=verify(Path(report['namespace']));checked_plan=json.loads((Path(checked['case'])/'plan.json').read_text())
    original=verify(Path(checked_plan['saved_full_Hessian']['namespace']));case=Path(original['case'])
    plan=json.loads((case/'plan.json').read_text());used=[review_path,Path(report['namespace']),Path(checked['case'])/'plan.json',
        Path(checked_plan['saved_full_Hessian']['namespace']),case/'plan.json']
    def read_json(p):used.append(p);return json.loads(p.read_text())
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','sampling_native_resource_owner')
    stage=review_closure(Path(checked['case']),read_json,owners.alive,
        expected_status='screening_polymer_Hessian_control_finished_review_pending',host_GiB=80,
        GPU_MiB=31744,worker=Path(checked['sources'])/'worker.py')
    construction=read_json(case/'status.json');guard=read_json(case/'residency_guard/owner.json')
    host=read_json(case/'host_guard/status.json');resident=read_json(case/'residency_guard/state.json')
    placements=list(case.parent.glob('placement*/completion.json'))
    if len(placements)!=1:raise ValueError('One exact Hessian construction placement required')
    completion=read_json(placements[0]);launch=read_json(placements[0].with_name('launch.json'))
    # This matrix completed before a recorded wall stop during followup work;
    # its separate fresh derivative/minimum review closed normally. Do not turn
    # the incomplete producer into a normal resource completion or G admission.
    if (construction['status']!='screening_polymer_Hessian_control_failed_preserved'
            or completion['budget_stop']!='wall_time_budget' or completion['returncode'] in (None,0)
            or owners.alive(construction['owner']) is not False or not guard.get('guard_stopped')
            or resident.get('foreign_resident_PIDs') or not resident.get('own_residency_observed')
            or any(guard.get(k)!=construction['owner'].get(k) for k in ('host','pid','process_start'))
            or not 0<max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])<=80*2**30
            or host['free_RAM_bytes']<32*2**30 or not 31*1024<=launch['device']['memory_MiB']<=33*1024):
        raise ValueError('Proven closed matrix-producing wall stop within the 32-GB/80-GiB class required')
    actual=read_json(case/'model/model.json')
    actual_model=actual['profile'];implementation=actual_model.get('implementation_sha256',{})
    if (set(actual_model)!=set(model)|{'implementation_sha256'} or any(actual_model[k]!=v for k,v in model.items())
            or not implementation or plan['response_available_memory_GiB']!=4
            or (case/'points/0000/evaluation.json').read_text()==''):
        raise ValueError('Neutral native witness has a different solvent/response model')
    from redox.validation.archived_sources import ArchivedSources
    reader=ArchivedSources(ROOT,Path(checked['archive']),relative_to_prior=True)
    for path,digest in implementation.items():used.append(reader.verify(path,digest))
    checkpoint=case/'SCF_evaluations/0000/0000/SCF.h5';mol=chkfile.load_mol(str(checkpoint));used.append(checkpoint)
    if (mol.natm,mol.nao_nr(),mol.nelec,mol.charge,mol.spin)!=(119,1193,(235,235),0,0):
        raise ValueError('Exact largest neutral RKS dimension witness required')
    return dict(status='corrected_neutral_Hessian_dimension_witness_native_pilot_required',atoms=119,nao=1193,nocc=235,
        response_available_memory_GiB=4,maximum_peak_GPU_MiB=31744,
        derivative_validation_peak_GPU_MiB=stage['observed_peak_GPU_MiB'],
        construction_device_memory_MiB=launch['device']['memory_MiB'],
        construction_host_peak_RSS_bytes=max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes']),
        construction_GPU_peak_retained=False,new_complete_native_resource_pilot_required=True,source_G_transferred=False,
        solvent_implementation_sha256=implementation,
        source_sha256={**report['source_sha256'],**{str(p):sha256_file(p) for p in used}})


def prepare(energy_review, resource_review, archive, output, native_pilot=None):
    from redox.core.runtime_paths import gpu4pyscf_source_root
    from redox.validation.sampling_receipts import checked_energy
    preliminary=json.loads(energy_review.read_text())
    report=checked_energy(energy_review,Path(preliminary['namespace']),preliminary['key'],preliminary['portable_method_id'])
    source=verify(Path(report['namespace']));parent=Path(source['case']);old=json.loads((parent/'plan.json').read_text());task=old['tasks'][0]
    if task!=report['task'] or old['experimental_method_id']!=report['portable_method_id']:
        raise ValueError('Exact electronic task and solvent model differ')
    prerequisite=None
    if native_pilot is not None:
        prerequisite=json.loads(native_pilot.read_text())
        if (prerequisite['status']!='corrected_sampling_native_independently_reconstructed'
                or prerequisite['key']!='all_same__occ_2_3__conf001' or not prerequisite['native_G_qualified']
                or prerequisite['experimental_model']!=old['experimental_model']
                or not 0<prerequisite['observed_peak_GPU_MiB']<=31744
                or not 0<prerequisite['peak_host_RSS_bytes']<=80*2**30
                or any(sha256_file(p)!=h for p,h in prerequisite['source_sha256'].items())):
            raise ValueError('Exact completed largest corrected native pilot required')
    elif task['key'] not in PILOT_KEYS:
        raise ValueError('Only the three declared native pilots may precede the largest complete native review')
    witness=resource_witness(resource_review,old['experimental_model'])
    atoms=read(parent/'input.xyz');molecule=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),
        unit='Angstrom',basis='def2-svp',charge=0,spin=0,verbose=0)
    if (task['charge'],task['mult'],task['candidate_id'])!=(0,1,'exp_phthalimide_control') or len(atoms)>119 or molecule.nao_nr()>1193 or molecule.nelec[0]>235:
        raise ValueError('Only the checked neutral phthalimide native dimension class is supported')
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Fresh durable native namespace outside both scientific checkouts required')
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir();sources.mkdir()
    for original,name in [(ROOT/'scripts/validation/run_corrected_sampling_native.py','worker.py'),
            (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
            (ROOT/'src/redox/qm/sampling_native_checkpoints.py','native_checkpoints.py'),
            (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),
            (ROOT/'src/redox/qm/response_diagnostic.py','response_helper.py'),
            (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
            (ROOT/'scripts/validation/review_corrected_sampling_native.py','reviewer.py'),
            (energy_review,'electronic_review.json'),(Path(__file__),'preparer.py')]:
        shutil.copyfile(original,sources/name)
    shutil.copyfile(parent/'input.xyz',case/'input.xyz');atomic_write_json(sources/'native_resource_witness.json',witness)
    if prerequisite is not None:shutil.copyfile(native_pilot,sources/'native_pilot_review.json')
    package=gpu4pyscf_source_root()
    paths=[*sources.glob('*'),case/'input.xyz',*(archive/'src/redox/qm').glob('*.py'),
        *(archive/'src/redox/core').glob('*.py'),*(archive/'config').glob('*.py'),*(archive/'config').glob('*.json'),
        archive/'src/redox/workflows/polymer_campaign.py',archive/'src/redox/workflows/resources.py',
        archive/'scripts/pipeline/smooth_polymer_native.py',archive/'scripts/pipeline/smooth_cds_control.py',
        archive/'scripts/pipeline/residual_response_solver.py',archive/'scripts/pipeline/run_pi_stacking.py',
        archive/'scripts/validation/audit/native_partition.py',archive/'scripts/fleet/exclusive_python_worker.py',
        archive/'scripts/fleet/tensor_runtime_worker.py',archive/'environments/cutensor22/runtime.json',
        package/'dft/gen_grid.py',package/'dft/numint.py',package/'lib/libgdft.so']
    hashes={str(p):sha256_file(p) for p in paths}
    plan=dict(role='corrected_sampling_native',corrected_sampling_native=True,source_id=task['key'],state='rest',
        tasks=[task],atoms=len(atoms),native_nao=molecule.nao_nr(),native_nocc=molecule.nelec[0],
        raw_root=str(case/'raw'),electronic_review=str(sources/'electronic_review.json'),
        experimental_method_id=old['experimental_method_id'],experimental_model=old['experimental_model'],
        electronic_protocol=old['electronic_protocol'],experimental_native_settings=native_settings(archive/'scripts/pipeline/smooth_polymer_native.py'),
        thermal_protocol=load(archive/'src/redox/qm/qrrho.py','sampling_native_thermal_definition').PROTOCOL,
        scored_geometry_sha256=sha256_file(case/'input.xyz'),runtime_profile=str(archive/'environments/cutensor22/runtime.json'),
        execution_driver=str(sources/'worker.py'),producer_sha256=hashes,input_sha256={},source_sha256=hashes,
        native_G_qualified=False,primary_candidate_admission=False,conformer_ensemble_qualified=False,
        prospective_resource_pilot=prerequisite is None,native_class_prerequisite=(None if prerequisite is None else str(sources/'native_pilot_review.json')),
        source_nuclear_or_G_transferred=False,**profile())
    if prerequisite is not None:
        pilot=verify(Path(prerequisite['namespace']))
        for name in ('worker.py','checkpoints.py','native_checkpoints.py','helper.py','response_helper.py','kernel.cu','host_guard.py'):
            if sha256_file(sources/name)!=sha256_file(Path(pilot['sources'])/name):
                raise ValueError('Native followup execution differs from the completed largest pilot')
    atomic_write_json(case/'plan.json',plan)
    spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),
        source_id=task['key'],state='rest',scientific_sources_read_only=True,
        maximum_host_RSS_GiB=80,minimum_free_GiB=32,plan_sha256=sha256_file(case/'plan.json'))
    atomic_write_json(output/'namespace.json',spec)
    return dict(status='corrected_sampling_native_prepared_no_jobs_submitted',key=task['key'],
        namespace=str(output/'namespace.json'),namespace_sha256=sha256_file(output/'namespace.json'),
        source_electronic_review=str(energy_review),native_G_qualified=False)


def dispatch(manifest_path,interpreter,storage_peers=(),avoid_hosts=()):
    from redox.workflows.charged_native_followup import dispatch_rows
    data=json.loads(manifest_path.read_text());rows=data['states']
    if data['declared_sampling_inputs']!=240 or len({r['key'] for r in rows})!=len(rows):
        raise ValueError('Distinct native rows within the complete 240-input sampling objective required')
    for row in rows:
        if sha256_file(row['namespace'])!=row['namespace_sha256']:raise ValueError('Native namespace changed')
        spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
        if (plan['role']!='corrected_sampling_native' or plan['tasks'][0]['key']!=row['key']
                or any(plan.get(k)!=v for k,v in profile().items())
                or (plan.get('native_class_prerequisite') is None and row['key'] not in PILOT_KEYS)):
            raise ValueError('Exact native pilot or independently released followup class required')
    return dispatch_rows(rows,manifest_path.parent,interpreter,[manifest_path,Path(__file__)],storage_peers=storage_peers,avoid_hosts=avoid_hosts)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',nargs='?',choices=('prepare','dispatch'),default='prepare')
    for name in ('energy-review','resource-review','archive','out','manifest','storage-peers'):p.add_argument('--'+name,type=Path)
    p.add_argument('--native-pilot',type=Path)
    p.add_argument('--python')
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':
        if any(v is None for v in (a.energy_review,a.resource_review,a.archive,a.out)):p.error('Preparation requires energy/resource reviews, archive and fresh output')
        print(json.dumps(prepare(a.energy_review.resolve(),a.resource_review.resolve(),a.archive.resolve(),a.out.resolve(),
            None if a.native_pilot is None else a.native_pilot.resolve())))
    else:
        if a.manifest is None or not a.python:p.error('Dispatch requires the exact prepared native manifest and interpreter')
        peers=[] if a.storage_peers is None else json.loads(a.storage_peers.read_text())
        result=dispatch(a.manifest.resolve(),a.python,peers)
        print(json.dumps(dict(status=result['status'],launches=len(result['launches']),refusals=result['refusals'],unsubmitted_count=len(result['unsubmitted']))))
