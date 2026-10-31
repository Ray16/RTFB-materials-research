#!/usr/bin/env python3
"""Export the checked exact-viologen redox/grafting case as a portable CPU example."""
import argparse
import json
import os
from pathlib import Path
import shutil

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.qm.native_thermal import NATIVE_THERMAL
from redox.validation.archived_sources import ArchivedSources
from redox.validation.molecular_native import MolecularSources


def export(candidate_path,grafting_path,reference_path,shared_package,archive,output):
    if output.exists():raise FileExistsError('Preserve earlier molecular examples')
    candidate=json.loads(candidate_path.read_text());graft=json.loads(grafting_path.read_text());reference=json.loads(reference_path.read_text())
    shared=json.loads((shared_package/'manifest.json').read_text())
    if (candidate['product_inchikey']!='PBNXFHGCSKYNHC-UHFFFAOYSA-N' or graft['product_inchikey']!=candidate['product_inchikey']
            or candidate['material_admission'] or graft['material_admission']):raise ValueError('Exact checked viologen with partial material scope required')
    for path,digest in graft['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Checked grafting source changed')
    output.mkdir(parents=True);rows=[];sources={};readers=[]
    def copy(path,target):
        path=Path(path);target=output/target;target.parent.mkdir(parents=True,exist_ok=True)
        digest=sha256_file(path);shutil.copyfile(path,target)
        if sha256_file(path)!=digest or sha256_file(target)!=digest:raise ValueError('Exported bytes differ')
        sources[str(path)]=digest;return str(target.relative_to(output))
    def state(name,smiles,q,m,paths,expected,fd=None):
        files={role:copy(path,Path('states')/name/(role+Path(path).suffix)) for role,path in paths.items()}
        row=dict(id=name,smiles=smiles,charge=q,mult=m,files=files,checked_G_by_cutoff_eV=expected,
            partition_kind='monatomic' if q==-1 and name=='chloride' else 'projected_DFT',finite_difference_replay=fd is not None)
        if fd is not None:
            case,native=fd;matrix=np.load(paths['hessian'],allow_pickle=False);base=[]
            for p in case.glob('SCF_evaluations/*/*/hessian_tensor.npy'):
                h=np.load(p,allow_pickle=False)
                if h.ndim==4 and h.transpose(0,2,1,3).shape[:2]==(h.shape[0],3):
                    if np.max(abs(h.transpose(0,2,1,3).reshape(matrix.shape)-matrix))<1e-12:base.append(p.parent)
            if len(base)!=1:raise ValueError('One checkpoint must bind the accepted Hessian')
            base_coords=np.asarray(json.loads((base[0]/'settings.json').read_text())['coordinates_Bohr']);points=[base[0]];checks=[]
            possible=[]
            for p in case.glob('SCF_evaluations/*/*/settings.json'):
                settings=json.loads(p.read_text())
                if (p.with_name('gradient.npy').exists() and p.with_name('SCF_result.json').exists()
                        and settings['charge']==q and settings['mult']==m and settings['conv_tol']==1e-12
                        and settings['xc']=='r2scan' and settings['basis']=='def2-svp'):
                    possible.append((p.parent,np.asarray(settings['coordinates_Bohr'])))
            for check in native['derivative_checks']:
                rowcheck=dict(check)
                for label,sign in [('minus',-1),('plus',1)]:
                    target=base_coords.copy();target[check['atom_index'],check['axis']]+=sign*check['step_Bohr']
                    matches=[p for p,x in possible if np.max(abs(x-target))<3e-9]
                    if len(matches)!=1:raise ValueError('One saved native gradient per selected displacement required')
                    points.append(matches[0]);rowcheck[label+'_index']=len(points)-1
                checks.append(rowcheck)
            coords=[];gradients=[];metadata=[]
            for point in points:
                settings=json.loads((point/'settings.json').read_text());solved=json.loads((point/'SCF_result.json').read_text())
                if solved['converged'] is not True:raise ValueError('Converged saved native SCF required')
                coords.append(settings['coordinates_Bohr']);gradients.append(np.load(point/'gradient.npy',allow_pickle=False))
                metadata.append(dict(settings=settings,SCF_result=solved,original_gradient=str(point/'gradient.npy'),original_gradient_sha256=sha256_file(point/'gradient.npy')))
                for f in ['settings.json','SCF_result.json','gradient.npy']:sources[str(point/f)]=sha256_file(point/f)
            target=output/'states'/name/'fd_arrays.npz';np.savez_compressed(target,coordinates_Bohr=np.asarray(coords),gradient_Ha_Bohr=np.asarray(gradients))
            row['files']['fd_arrays']=str(target.relative_to(output));target=target.with_name('fd_records.json')
            atomic_write_json(target,dict(checks=checks,points=metadata));row['files']['fd_records']=str(target.relative_to(output))
        rows.append(row)
    def molecular(path,name):
        report=json.loads(path.read_text());spec_path=Path(report['namespace']);spec=json.loads(spec_path.read_text());reader=MolecularSources(spec)
        if report['status'] not in ('exact_molecular_native_G_independently_checked','reaction_component_native_G_independently_checked') or not report['resource_pilot_passed_for_this_state']:
            raise ValueError('Checked molecular native/resource review required')
        reader.verify_all(report['source_sha256']);readers.append(reader);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text());result=json.loads((case/'result.json').read_text())
        native=reader.path(result['native_record_path']);partition=[reader.path(p) for p in report['source_sha256'] if Path(p).parent.name=='projected_thermal' and Path(p).suffix=='.json']
        if len(partition)!=1:raise ValueError('One exact projected nuclear partition required')
        paths=dict(geometry=reader.path(Path(plan['target_directory'])/'opt.xyz'),result=reader.path(Path(plan['target_directory'])/'result.json'),
            energy=reader.path(result['active_record_path']),native=native,partition=partition[0],hessian=native.with_suffix('')/'hessian_cartesian.npy')
        state(name,plan.get('product_smiles',plan.get('smiles')),report['charge'],report['mult'],paths,report['G_by_cutoff_eV'],(case,json.loads(native.read_text())))
        sources[str(path)]=sha256_file(path);sources[str(spec_path)]=sha256_file(spec_path)
    for row in candidate['native_state_reviews']:
        path=Path(row['path'])
        if sha256_file(path)!=row['sha256']:raise ValueError('Candidate native review changed')
        report=json.loads(path.read_text());molecular(path,'graft_'+report['state'])
    precursor_reviews=[Path(p) for p in graft['source_sha256'] if p.endswith('independent_review.json') and 'viologen_grafting_2026-10-31' in p and '/existing_product_' not in p]
    if len(precursor_reviews)!=3:raise ValueError('Three source-bound precursor reviews required')
    for path in precursor_reviews:
        report=json.loads(path.read_text());spec=json.loads(Path(report['namespace']).read_text());plan=json.loads((Path(spec['case'])/'plan.json').read_text());molecular(path,'precursor_'+str(plan['conformer_index']))
    reader=ArchivedSources(ROOT,archive,relative_to_prior=True);reader.verify_all(reference['source_sha256']);readers.append(reader)
    for label,q,m in [('neu',0,1),('ox',1,2)]:
        prefix=str(ROOT/'calcs/dft/ferrocene'/label);bound=reference['source_sha256']
        def unique(folder):
            found=[reader.path(p) for p in bound if p.startswith(prefix+'/') and Path(p).parent.name==folder and Path(p).suffix=='.json']
            if len(found)!=1:raise ValueError('One matched Fc surface/energy/partition required')
            return found[0]
        native=unique('native_thermal');paths=dict(geometry=reader.path(prefix+'/opt.xyz'),result=reader.path(prefix+'/result.json'),energy=unique('sp'),native=native,partition=unique('projected_thermal'),hessian=native.with_suffix('')/'hessian_cartesian.npy')
        fc=next(r for r in shared['states'] if r['id']==('fc_neutral' if label=='neu' else 'fc_oxidized'))
        state('fc_'+label,fc['smiles'],q,m,paths,reference['states']['ferrocene__'+label]['G_by_cutoff_eV'])
    for relative,digest in shared['files_sha256'].items():
        if sha256_file(shared_package/relative)!=digest:raise ValueError('Shared reference cache changed')
    for name,role in [('local','site'),('chloride','byproduct')]:
        row=next(r for r in shared['states'] if r['id']==name)
        state(name,row['smiles'],row['charge'],row['mult'],{k:shared_package/v for k,v in row['files'].items()},graft['components'][role]['G_by_cutoff_eV'])
    for name in ['projected_modes.py','nuclear_equations.py']:copy(shared_package/name,name)
    copy(shared_package/'requirements-cpu.lock.txt','requirements.lock.txt')
    copy(ROOT/'examples/exact_viologen_replay/replay.py','replay.py');copy(Path(__file__),'exporter_source.py')
    for p in (candidate_path,grafting_path,reference_path):sources[str(p)]=sha256_file(p)
    reactions=[]
    for i,(ox,red) in enumerate([('ox2','ox1'),('ox1','neu')]):
        reactions.append(dict(id='E'+str(i+1)+'_V_vs_Fc',kind='one_electron_potential_vs_Fc',electrons=1,
            coefficients={'graft_'+red:1,'fc_ox':1,'graft_'+ox:-1,'fc_neu':-1},expected_by_cutoff=candidate['redox']['events'][i]['potential_V_vs_Fc_by_cutoff']))
    for i,cycle in enumerate(graft['per_precursor_cycle']):
        index=graft['precursor_conformers'][i]['conformer_index'];reactions.append(dict(id='graft_precursor_'+str(index)+'_kJ_mol',kind='grafting_free_energy',
            coefficients={'graft_ox2':1,'chloride':1,'precursor_'+str(index):-1,'local':-1},expected_by_cutoff=cycle['deltaG_kJ_mol_by_cutoff']))
    reactions.append(dict(id='disproportionation_kJ_mol',kind='reaction_free_energy',coefficients={'graft_ox2':1,'graft_neu':1,'graft_ox1':-2},expected_by_cutoff=candidate['redox']['disproportionation_kJ_mol_by_cutoff']))
    (output/'README.md').write_text('# Exact viologen cached CPU example\n\nRun `python replay.py replay --package .` after installing the hash-locked NumPy requirement with CPython 3.11/Linux x86_64. This reconstructs two Fc-referenced potentials and a local quaternization from ten cached states. Nine Hessians/mode sets and nuclear partitions reconstruct; six recent states also retain displaced-gradient arrays. Historical Fc/site FD entries remain metadata-only. No SCF, xTB library or original project path is needed.\n\nThe example uses a capped local graft, MeCN, 298.15 K and separated chloride. Individual nuclear terms retain the source 1-atm convention; all displayed balanced reactions have zero net solute count, so the 1-M correction cancels. No ion-pairing, counterion-exchange, reflux-temperature or polymer correction is supplied. Spin-only Fc and declared candidate spin states are retained. This is finite sampling and does not admit a complete material or revalidate full historical executable/scratch closure.\n')
    for reader in readers:reader.recheck()
    if any(sha256_file(p)!=h for p,h in sources.items()):raise ValueError('A selected source changed during export')
    manifest=dict(schema=1,product_inchikey=candidate['product_inchikey'],product_smiles=candidate['product_smiles'],states=rows,reactions=reactions,
        active_protocol=ACTIVE_SP,native_protocol=NATIVE_THERMAL,precursor_ids=['precursor_'+str(i) for i in range(3)],selected_precursor_id='precursor_'+str(graft['selected_precursor_conformer_index']),
        source_sha256=sources,material_admission=False,global_conformer_convergence_proven=False,experimental_exact_graft_accuracy_validated=False,full_historical_source_closure_revalidated=False,
        files_sha256={str(p.relative_to(output)):sha256_file(p) for p in output.rglob('*') if p.is_file()})
    atomic_write_json(output/'manifest.json',manifest);return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('candidate','grafting','reference','shared-package','archive','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=export(a.candidate.resolve(),a.grafting.resolve(),a.reference.resolve(),a.shared_package.resolve(),a.archive.resolve(),a.out.resolve())
    print(json.dumps(dict(status='exact_molecular_CPU_example_exported',states=len(result['states']),files=len(result['files_sha256']))))
