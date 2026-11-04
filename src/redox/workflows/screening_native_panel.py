"""Prepare primary-model molecular Hessian repeats for the exact lead package.

The existing native thermal and projected partition implementation is retained.
These are numerical method controls, including references/precursors, with no
new synthesis or material admission. The largest RKS/UKS repeats are pilots.
"""
import argparse
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil

from ase.io import read
from pyscf import gto

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.screening_component_panel import inventory
from redox.workflows.molecular_native_namespace import verify


def prepare(archive, package, energy_reviews, output):
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Fresh durable namespaces outside both scientific checkouts required')
    manifest, entries=inventory(package);states={s['id']:s for s in manifest['states']}
    candidates=[];atomic=[]
    for row in entries:
        state=states[row['key']];files=state['files'];atoms=read(row['geometry'])
        if len(atoms)==1:
            atomic.append(dict(key=row['key'],status='monatomic_partition_no_Hessian_repeat',state=state));continue
        for name in ('result','native','partition','hessian'):
            p=package/files[name]
            if sha256_file(p)!=manifest['files_sha256'][files[name]]:
                raise ValueError('Native reference input changed')
        raw=json.loads((package/files['result']).read_text())
        if ((raw['charge'],raw['mult'])!=(state['charge'],state['mult']) or not raw.get('optimized')
                or any(raw.get(k)!=v for k,v in dict(opt_xc='r2scan',opt_basis='def2-svp',opt_disp='d4',solvent='acetonitrile').items())):
            raise ValueError('Exact previously admitted primary geometry required')
        review_path=energy_reviews/(row['key']+'.json');review=json.loads(review_path.read_text())
        if (review['status']!='benchmark_screening_energy_control_independently_checked'
                or not review['actual_checkpoint_frames_verified'] or review['key']!=row['key']
                or (review['source_state']['charge'],review['source_state']['mult'])!=(state['charge'],state['mult'])
                or any(sha256_file(p)!=h for p,h in review['source_sha256'].items())):
            raise ValueError('Exact independent corrected-electronic review required')
        energy_ns=Path(review['namespace']);energy_spec=json.loads(energy_ns.read_text())
        energy_case=Path(energy_spec['case']);energy_plan=json.loads((energy_case/'plan.json').read_text())
        if sha256_file(Path(energy_spec['sources'])/'input.xyz')!=sha256_file(row['geometry']):
            raise ValueError('Corrected electronic energy and Hessian reference geometry differ')
        molecule=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),
                       charge=state['charge'],spin=state['mult']-1,basis='def2-svp',verbose=0)
        if molecule.nao>1193:raise ValueError('Molecular native control exceeds the established dimension class')
        candidates.append((row,state,raw,review,review_path,energy_case,molecule.nao))
    output.mkdir(parents=True);rows=[]
    for row,state,raw,review,review_path,energy_case,nao in candidates:
        root=output/'states'/row['key'];case=root/'case';sources=root/'sources'
        case.mkdir(parents=True);sources.mkdir();(case/'state').mkdir();(case/'library').mkdir()
        for name in ('dft','locks','repairs'):(case/name).mkdir()
        copies=[(ROOT/'scripts/validation/run_screening_native_control.py','worker.py'),
            (ROOT/'src/redox/core/hashing.py','hashing.py'),
            (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),
            (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
            (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
            (Path(__file__),'preparer.py'),(review_path,'corrected_electronic_review.json'),
            (energy_case/'corrected/result.json','corrected_electronic.json')]
        copies += [(package/state['files'][name],'reference_'+Path(state['files'][name]).name)
                   for name in ('geometry','result','native','partition','hessian')]
        for original,name in copies:shutil.copyfile(original,sources/name)
        shutil.copyfile(row['geometry'],case/'state/opt.xyz')
        shutil.copyfile(package/state['files']['result'],case/'state/result.json')
        # A private one-state manifest declares the same exact graph/state for
        # the existing native function; the production library is read-only.
        with (case/'library/manifest.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=['id','state','smiles','charge','mult_hint'])
            writer.writeheader();writer.writerow(dict(id=raw['id'],state=raw['state'],smiles=state['smiles'],charge=state['charge'],mult_hint=state['mult']))
        bound=[*sources.glob('*'),case/'state/opt.xyz',case/'state/result.json',case/'library/manifest.csv',
               *(archive/'src/redox/qm').glob('*.py'),*(archive/'src/redox/core').glob('*.py'),
               *(archive/'config').glob('*.py'),*(archive/'config').glob('*.json'),
               archive/'scripts/pipeline/run_pi_stacking.py',archive/'scripts/validation/audit/native_partition.py',
               archive/'scripts/fleet/exclusive_python_worker.py',archive/'scripts/fleet/tensor_runtime_worker.py',
               archive/'environments/cutensor22/runtime.json',archive/'src/redox/workflows/resources.py']
        plan=dict(role='primary_molecular_Hessian_control',primary_molecular_Hessian_control=True,
            source_id=row['key'],state='rest',native_id=raw['id'],native_state=raw['state'],smiles=state['smiles'],
            charge=state['charge'],mult=state['mult'],atoms=len(read(row['geometry'])),native_nao=nao,
            native_memory_cap_GiB=6,maximum_host_RSS_GiB=40,minimum_free_GiB=32,
            maximum_peak_GPU_MiB=31744,maximum_wall_seconds=14400,maximum_durable_output_bytes=2**30,
            minimum_initial_storage_GiB=33,minimum_live_storage_GiB=32,
            source_sha256={**review['source_sha256'],**{str(p):sha256_file(p) for p in bound}},
            reference_checked_G_by_cutoff_eV=state['checked_G_by_cutoff_eV'],
            purpose='Fixed-source primary-model native repeat under corrected screening; no synthesis eligibility or ensemble inference',
            native_G_qualified=False,material_admission=False,geometry_repair_permitted=False)
        atomic_write_json(case/'plan.json',plan)
        spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),
            source_id=row['key'],state='rest',plan_sha256=sha256_file(case/'plan.json'),
            maximum_host_RSS_GiB=40,minimum_free_GiB=32,scientific_sources_read_only=True,
            library_overlay=str(case/'library'))
        atomic_write_json(root/'namespace.json',spec)
        rows.append(dict(key=row['key'],namespace=str(root/'namespace.json'),namespace_sha256=sha256_file(root/'namespace.json'),nao=nao,mult=state['mult']))
    pilots={max([r for r in rows if (r['mult']>1)==open_shell],key=lambda r:r['nao'])['key'] for open_shell in (False,True)}
    result=dict(status='exact_primary_native_controls_prepared_no_jobs_submitted',states=rows,
        target_native_states=len(rows),monatomic_components=atomic,pilot_keys=sorted(pilots),
        source_package_manifest_sha256=sha256_file(package/'manifest.json'),source_package=str(package),
        quantum_jobs_submitted=0,new_complete_free_energy_cycles_qualified=False)
    atomic_write_json(output/'manifest.json',result);return result


def qualified_rows(manifest_path,reviews):
    data=json.loads(manifest_path.read_text());rows=data['states']
    if len(rows)!=data['target_native_states'] or len({r['key'] for r in rows})!=len(rows):
        raise ValueError('Distinct complete primary-native roster required')
    expected={max([r for r in rows if (r['mult']>1)==open_shell],key=lambda r:r['nao'])['key'] for open_shell in (False,True)}
    if expected!=set(data['pilot_keys']):raise ValueError('Exact largest RKS/UKS pilots required')
    checked=[json.loads(p.read_text()) for p in reviews]
    if len(checked)!=2 or {r['key'] for r in checked}!=expected:
        raise ValueError('Both largest primary-native pilot reviews required')
    pilots={}
    for report in checked:
        if (report['status']!='primary_molecular_native_control_independently_reconstructed'
                or not report['native_component_requalified'] or report['material_admission']
                or report['conformer_ensemble_qualified'] or report['actual_SCF_checkpoint_frames_checked']!=9
                or not 0<report['peak_GPU_MiB']<=31744 or not 0<=report['maximum_column_error']<=.01
                or not -5<=report['minimum_frequency_cm']<float('inf')
                or not report.get('source_sha256') or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
            raise ValueError('Complete physical/resource pilot receipts required')
        row=next(r for r in rows if r['key']==report['key'])
        if Path(row['namespace']).resolve()!=Path(report['namespace']).resolve():raise ValueError('Pilot namespace differs')
        spec=verify(Path(report['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
        raw=json.loads((Path(spec['case'])/'result.json').read_text())
        if ((report['charge'],report['mult'])!=(plan['charge'],plan['mult']) or report['mult']!=row['mult']
                or set(report['G_by_cutoff_eV'])!={'25','50','100'}
                or any(not math.isfinite(v) or abs(v-raw['G_by_cutoff_eV'][c])>1e-8 for c,v in report['G_by_cutoff_eV'].items())):
            raise ValueError('Pilot state and complete raw G values must match the independent review')
        pilots[report['mult']>1]=spec
    for row in rows:
        if sha256_file(row['namespace'])!=row['namespace_sha256']:raise ValueError('Prepared primary-native namespace changed')
        spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
        pilot=pilots[row['mult']>1];pp=json.loads((Path(pilot['case'])/'plan.json').read_text())
        if (plan['role']!='primary_molecular_Hessian_control' or plan['native_nao']>pp['native_nao']
                or plan['native_memory_cap_GiB']!=pp['native_memory_cap_GiB']
                or plan['mult']!=row['mult'] or plan['native_nao']!=row['nao'] or plan['charge'] not in (0,1,2)):
            raise ValueError('Primary native followup exceeds the declared class')
        for name in ('native_thermal.py','projected_thermal.py','qrrho.py','thermal.py','dft.py','cds_finite.py','hessian_memory.py','xtb_kernel.py','vibrations.py'):
            original=str(Path(spec['archive'])/'src/redox/qm'/name);reference=str(Path(pilot['archive'])/'src/redox/qm'/name)
            if plan['source_sha256'][original]!=pp['source_sha256'][reference]:
                raise ValueError('Primary native quantum/thermal implementation differs')
        for name in ('worker.py','helper.py','checkpoints.py','hashing.py','kernel.cu'):
            if sha256_file(Path(spec['sources'])/name)!=sha256_file(Path(pilot['sources'])/name):
                raise ValueError('Primary native worker/kernel implementation differs')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',nargs='?',choices=('prepare','dispatch'),default='prepare')
    for key in ('archive','package','energy-reviews','out','manifest','storage-peers'):p.add_argument('--'+key,type=Path)
    p.add_argument('--pilot-review',type=Path,action='append',default=[]);p.add_argument('--python')
    p.add_argument('--avoid-host',action='append',default=[])
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':
        if any(v is None for v in (a.archive,a.package,a.energy_reviews,a.out)):p.error('Preparation requires archive, package, energy reviews and fresh output')
        r=prepare(a.archive.resolve(),a.package.resolve(),a.energy_reviews.resolve(),a.out.resolve())
        print(json.dumps({k:r[k] for k in ('status','target_native_states','pilot_keys')}))
    else:
        if a.manifest is None or not a.python:p.error('Dispatch requires the prepared manifest, exact pilot reviews and interpreter')
        rows=qualified_rows(a.manifest.resolve(),[r.resolve() for r in a.pilot_review])
        peers=[] if a.storage_peers is None else json.loads(a.storage_peers.read_text())
        from redox.workflows.charged_native_followup import dispatch_rows
        r=dispatch_rows(rows,a.manifest.resolve().parent,a.python,
            [a.manifest.resolve(),*a.pilot_review,Path(__file__),*([a.storage_peers] if a.storage_peers else [])],a.avoid_host,peers)
        print(json.dumps(dict(status=r['status'],launches=len(r['launches']),refusals=r['refusals'],unsubmitted=r['unsubmitted'])))
