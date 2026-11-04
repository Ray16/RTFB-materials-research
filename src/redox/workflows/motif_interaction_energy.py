"""Source-bound electronic A/B/AB controls for frozen charged motif pairs.

This is an operational capped-fragment interaction, including continuum response.
It is separate from whole-chain nonadditivity, deformation, binding G and coupling.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil

from ase.io import read
from pyscf import gto

from redox.core.common import ROOT, HARTREE_EV, EV_KJ
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify, load
from redox.workflows.ionic_stationary_cycles import validate_cycle_balance

PILOT = 'all_same__occ_2_3__cap109'


def validate_construction_binding(spec,plan,report):
    data=json.loads((Path(spec['sources'])/'construction.json').read_text())
    matches=[r for r in data['entries'] if r['key']==plan['source_id']]
    if len(matches)!=1:raise ValueError('The exact declared parent/cap triplet is required')
    entry=matches[0];expected={r['id']:r for r in entry['molecules']}
    if (plan['parent_key']!=entry['parent_key'] or plan['cap_length_A']!=entry['cap_length_A']
            or report['parent_key']!=entry['parent_key'] or report['cap_length_A']!=entry['cap_length_A']
            or {r['id'] for r in plan['tasks']}!={'A','B','AB'}):
        raise ValueError('Parent-frame or cap-length labels differ from the audited construction')
    for task in plan['tasks']:
        declared=expected[task['id']]
        if (any(task[k]!=declared[k] for k in ('geometry_sha256','graph_sha256','charge','smiles','atoms'))
                or sha256_file(task['geometry'])!=declared['geometry_sha256']
                or sha256_file(task['molecular_graph'])!=declared['graph_sha256']):
            raise ValueError('The calculated A/B/AB input differs from the audited exact parent/cap frame')


def checked_triplet(path,expected_key):
    report=json.loads(path.read_text())
    if (report['status']!='frozen_charged_motif_triplet_independently_checked' or report['key']!=expected_key
            or report['actual_SCF_frames_checked']!=6 or report['native_G_qualified'] or report['geometry_optimized']
            or report['material_admission'] or not report['source_sha256']
            or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
        raise ValueError('Complete independent frozen-motif electronic pilot required')
    for key,maximum in [('peak_GPU_MiB',31744),('peak_host_RSS_bytes',16*2**30)]:
        value=report[key]
        if type(value) not in (int,float) or not math.isfinite(value) or not 0<value<=maximum:
            raise ValueError('The frozen-motif pilot must meet its declared resource class')
    spec=verify(Path(report['namespace']));case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text());raw=json.loads((case/'result.json').read_text())
    if (plan['role']!='frozen_charged_motif_energy' or plan['source_id']!=expected_key
            or report['method']!=plan['method'] or report['solvent_model']!=plan['solvent_model']
            or report['exact_energies_Ha']!=raw['energies']):
        raise ValueError('Pilot protocol and actual electronic components differ')
    validate_construction_binding(spec,plan,report)
    for medium in ('gas','acetonitrile'):
        values={key:raw['energies'][key][medium]['E_Ha'] for key in ('A','B','AB')}
        expected=math.fsum([values['AB'],-values['A'],-values['B']])*HARTREE_EV*EV_KJ
        if not math.isfinite(expected) or abs(expected-report['interaction_E_kJ_mol_by_medium'][medium])>1e-7:
            raise ValueError('The checked interaction does not reproduce its exact raw electronic parents')
    return report,spec


def checked_pilot(path):return checked_triplet(path,PILOT)


def collect(manifest_path,output):
    if output.exists():raise FileExistsError('Preserve earlier frozen-motif collections')
    from redox.validation.execution_review import freeze_source
    manifest=json.loads(manifest_path.read_text());rows=manifest['states']
    if manifest['target_triplets']!=12 or len(rows)!=12 or len({r['key'] for r in rows})!=12:
        raise ValueError('All four frames and three cap lengths must remain in the ledger')
    reviewer=load(ROOT/'scripts/validation/review_motif_interaction_energy.py','local_motif_triplet_actual_review')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','local_motif_collection_owner')
    reviews=output.parent/'triplet_reviews';reviews.mkdir(parents=True,exist_ok=True)
    records=[];errors=[];producer=freeze_source(output.parent,Path(__file__))
    sources={str(manifest_path):sha256_file(manifest_path),str(producer):sha256_file(producer)}
    for row in rows:
        namespace=Path(row['namespace']);case=Path(json.loads(namespace.read_text())['case'])
        if sha256_file(namespace)!=row['namespace_sha256']:raise ValueError('Frozen motif namespace changed')
        record=dict(key=row['key'],namespace=str(namespace),interaction_E_kJ_mol_by_medium=None,
            status='not_scientifically_started',native_G_qualified=False)
        path=reviews/(row['key']+'.json')
        try:
            if not path.exists():
                status=case/'status.json'
                if not status.exists():records.append(record);continue
                stage=json.loads(status.read_text());record['status']=stage['status']
                if stage['status']!='frozen_charged_motif_energy_finished_independent_review_pending':records.append(record);continue
                if owners.alive(stage['owner']) is not False or not list(case.parent.glob('placement*/completion.json')):
                    record['status']='execution_closure_pending';records.append(record);continue
                atomic_write_json(path,reviewer.review(namespace))
            report,_=checked_triplet(path,row['key'])
            record.update(status='frozen_motif_E_independently_checked',parent_key=report['parent_key'],
                cap_length_A=report['cap_length_A'],interaction_E_kJ_mol_by_medium=report['interaction_E_kJ_mol_by_medium'],
                method=report['method'],solvent_model=report['solvent_model'],
                checked_sums=report['checked_sums'],review=str(path),review_sha256=sha256_file(path))
            sources[str(path)]=sha256_file(path)
        except Exception as error:
            record.update(status='independent_E_review_failed_preserved',error=repr(error));errors.append(dict(key=row['key'],error=repr(error)))
        records.append(record)
    sensitivity={}
    for parent in sorted({r.get('parent_key') for r in records if r.get('parent_key')}):
        group=[r for r in records if r.get('parent_key')==parent]
        complete={r['cap_length_A'] for r in group}=={1.07,1.09,1.11}
        sensitivity[parent]=dict(checked_cap_lengths=[r['cap_length_A'] for r in group],
            complete=complete,range_kJ_mol_by_medium=({m:max(r['interaction_E_kJ_mol_by_medium'][m] for r in group)
                -min(r['interaction_E_kJ_mol_by_medium'][m] for r in group) for m in ('gas','acetonitrile')} if complete else None))
    result=dict(status='twelve_triplet_frozen_charged_motif_E_collection',states=records,
        checked_triplets=sum(r['interaction_E_kJ_mol_by_medium'] is not None for r in records),
        declared_triplets=12,independent_review_errors=errors,cap_sensitivity=sensitivity,
        basis_superposition_corrected=False,free_energy_computed=False,electronic_coupling_computed=False,
        full_polymer_nonadditivity_replaced=False,deformation_or_pure_backbone_energy_supplied=False,
        material_admission=False,source_sha256=sources)
    atomic_write_json(output,result);return result


def prepare(archive, construction, output, key=PILOT, pilot_review=None):
    from redox.core.runtime_paths import gpu4pyscf_source_root
    summary=json.loads((construction/'summary.json').read_text())
    if summary['status']!='four_ionic_parent_frames_methyl_fragment_construction_checked' or summary['quantum_jobs_submitted']:
        raise ValueError('The four checked ionic parent-frame construction is required')
    path=Path(summary['construction']);proof=Path(summary['independent_review'])
    if sha256_file(path)!=summary['construction_sha256'] or sha256_file(proof)!=summary['independent_review_sha256']:
        raise ValueError('Construction or independent review changed')
    checked=json.loads(proof.read_text())
    if checked['status']!='fixed_parent_local_motif_construction_independently_checked' or any(sha256_file(p)!=h for p,h in checked['source_sha256'].items()):
        raise ValueError('Complete source-bound fragment/pose audit required')
    matches=[r for r in json.loads(path.read_text())['entries'] if r['key']==key]
    if len(matches)!=1:raise ValueError('One declared frozen motif triplet required')
    entry=matches[0]
    if key!=PILOT and pilot_review is None:raise ValueError('The full charged A/B/AB pilot must pass before followups')
    if pilot_review is not None:
        reference,pilot_spec=checked_pilot(pilot_review)
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Fresh durable namespace outside scientific source trees required')
    tasks=[]
    for molecule in entry['molecules']:
        geometry=Path(molecule['geometry']);graph=Path(molecule['molecular_graph'])
        if sha256_file(geometry)!=molecule['geometry_sha256'] or sha256_file(graph)!=molecule['graph_sha256']:
            raise ValueError('Exact frozen motif geometry or graph changed')
        atoms=read(geometry);charge=4 if molecule['id']=='AB' else 2
        if molecule['charge']!=charge or len(atoms)!=(56 if charge==4 else 28):raise ValueError('Only two methyl viologen dications and their pair are admitted')
        mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis='def2-svp',charge=charge,spin=0,verbose=0)
        tasks.append(dict(molecule,charge=charge,mult=1,nao=mol.nao_nr()))
    validate_cycle_balance(dict(coefficients={'AB':1,'A':-1,'B':-1}),{r['id']:r for r in tasks})
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir();sources.mkdir()
    for source,name in [(ROOT/'scripts/validation/run_motif_interaction_energy.py','worker.py'),
        (ROOT/'src/redox/core/hashing.py','hashing.py'),(ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
        (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),(ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
        (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),(path,'construction.json'),
        (proof,'construction_review.json'),(Path(__file__),'preparer.py')]:shutil.copyfile(source,sources/name)
    for task in tasks:
        for field,extension in [('geometry','xyz'),('molecular_graph','mol')]:
            destination=sources/(task['id']+'.'+extension);shutil.copyfile(task[field],destination);task[field]=str(destination)
    if pilot_review is not None:shutil.copyfile(pilot_review,sources/'pilot_review.json')
    model=load(archive/'scripts/pipeline/smooth_cds_control.py','motif_energy_solvent_definition')
    if model.MODEL!=summary['physical_context']['solvent_model']:raise ValueError('Parent and frozen motif solvent models differ')
    installed=gpu4pyscf_source_root()
    files=[*sources.glob('*'),*(archive/'src/redox/qm').glob('*.py'),*(archive/'src/redox/core').glob('*.py'),
        *(archive/'config').glob('*.py'),*(archive/'config').glob('*.json'),archive/'src/redox/workflows/resources.py',
        archive/'scripts/pipeline/smooth_cds_control.py',archive/'scripts/validation/qualify_smooth_cds.py',
        archive/'scripts/pipeline/run_pi_stacking.py',archive/'scripts/fleet/exclusive_python_worker.py',
        archive/'scripts/fleet/tensor_runtime_worker.py',archive/'environments/cutensor22/runtime.json',
        installed/'dft/gen_grid.py',installed/'dft/numint.py',installed/'lib/libgdft.so']
    plan=dict(role='frozen_charged_motif_energy',source_id=key,state='frozen',tasks=tasks,parent_key=entry['parent_key'],
        cap_length_A=entry['cap_length_A'],method=dict(xc='r2scan',basis='def2-svp',disp='d4',density_fit=True,
            conv_tol=1e-12,conv_tol_grad=1e-8,max_cycle=200,grid_level=3),solvent_model=model.MODEL,
        media=['gas','acetonitrile'],geometry_optimized=False,maximum_host_RSS_GiB=16,minimum_free_GiB=32,
        maximum_peak_GPU_MiB=31744,maximum_wall_seconds=7200,maximum_durable_output_bytes=512*2**20,
        minimum_initial_storage_GiB=32.5,minimum_live_storage_GiB=32,source_sha256={str(p):sha256_file(p) for p in files},
        prospective_resource_pilot=pilot_review is None,pilot_review=None if pilot_review is None else str(sources/'pilot_review.json'),
        native_G_qualified=False,electronic_coupling_computed=False,material_admission=False,
        purpose='Frozen methyl-capped dication pair E_AB-E_A-E_B; gas and alternate smooth-surface MeCN separately; no binding G, populations or complete backbone/cap decomposition')
    if pilot_review is not None:
        if reference['method']!=plan['method'] or reference['solvent_model']!=plan['solvent_model']:
            raise ValueError('Followup and pilot physical contexts differ')
        for name in ('worker.py','hashing.py','checkpoints.py','helper.py','kernel.cu','host_guard.py'):
            if sha256_file(sources/name)!=sha256_file(Path(pilot_spec['sources'])/name):
                raise ValueError('Followup calculation or guard differs from its qualified pilot')
    atomic_write_json(case/'plan.json',plan)
    spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),source_id=key,state='frozen',
        maximum_host_RSS_GiB=16,minimum_free_GiB=32,scientific_sources_read_only=True,plan_sha256=sha256_file(case/'plan.json'))
    atomic_write_json(output/'namespace.json',spec);return spec


def launch(namespace,gpu,output):
    spec=verify(namespace);plan=json.loads((Path(spec['case'])/'plan.json').read_text())
    if plan['role']!='frozen_charged_motif_energy' or (plan['source_id']!=PILOT and plan['pilot_review'] is None):
        raise ValueError('Exact charged triplet pilot or independently released followup required')
    if plan['pilot_review'] is not None:checked_pilot(Path(plan['pilot_review']))
    from redox.workflows.bounded_molecular import launch as execute
    return execute(namespace,gpu,output,caller_source=Path(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','launch','collect'))
    for name in ('archive','construction','out','spec','pilot-review','manifest'):p.add_argument('--'+name,type=Path)
    p.add_argument('--key',default=PILOT);p.add_argument('--gpu',type=int);a=p.parse_args()
    if a.action=='prepare':
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
        if a.archive is None or a.construction is None or a.out is None:p.error('Prepare requires archive, checked construction and fresh output')
        spec=prepare(a.archive.resolve(),a.construction.resolve(),a.out.resolve(),a.key,None if a.pilot_review is None else a.pilot_review.resolve())
        print(json.dumps(dict(status='frozen_charged_motif_triplet_prepared',namespace=str(a.out/'namespace.json'))))
    elif a.action=='launch':
        if a.spec is None or a.out is None or a.gpu is None:p.error('Launch requires exact namespace, GPU and placement')
        result=launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps({k:result[k] for k in ('status','returncode','budget_stop')}))
    else:
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
        if a.manifest is None or a.out is None:p.error('Collection requires the complete manifest and fresh output')
        result=collect(a.manifest.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],checked_triplets=result['checked_triplets'],errors=result['independent_review_errors'])))
