#!/usr/bin/env python
"""Build a config-driven occupancy panel, sample MeCN conformers, freeze QM inputs."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.polymer_cases import panel


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _pin():
    permitted=sorted(os.sched_getaffinity(0));identity=mp.current_process()._identity
    os.sched_setaffinity(0,{permitted[((identity[0]-1) if identity else 0)%len(permitted)]})


def sample(payload):
    _pin()
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from ase.io import read
    from pyscf import gto
    from redox.build.polymer_contacts import starting_ensemble
    from redox.build.structures import to_xyz
    from redox.build.geometry import geometry_stereo_status
    from redox.core.protocol import ACTIVE_SP
    base,case,fragment,sampling=payload;base=Path(base)
    spec=importlib.util.spec_from_file_location('validated_seed',ROOT/'scripts/pipeline/prepare_pi_stacking.py')
    seed=importlib.util.module_from_spec(spec);spec.loader.exec_module(seed)
    audit_spec=importlib.util.spec_from_file_location('geometry_audit',ROOT/'scripts/validation/audit_polymer_geometries.py')
    audit=importlib.util.module_from_spec(audit_spec);audit_spec.loader.exec_module(audit)
    directory=base/'ensembles'/case['key'];receipt=directory/'ensemble.json'
    if receipt.exists():
        previous=json.loads(receipt.read_text());assert previous['case']==case
        for row in previous['selected']:assert sha(row['xyz'])==row['xyz_sha256']
        return previous
    directory.mkdir(parents=True,exist_ok=True);records=[];failures=[]
    modes=['unbiased','open','slipped'] if len(case['occupied_sites'])>=2 else ['unbiased']
    for mode in modes:
        try:
            mol,groups,ranked=starting_ensemble(case['smiles'],fragment,mode,
                count=sampling['embed_count_per_mode'],seed=sampling['random_seed'])
            assert len(groups)==len(case['occupied_sites'])
        except Exception as exc:failures.append(dict(mode=mode,stage='embedding',error=repr(exc)));continue
        for i,(_,cid) in enumerate(ranked[:sampling['xtb_seeds_per_mode']]):
            run=directory/f'{mode}_{i}'
            if run.exists():
                row=json.loads((run/'record.json').read_text()) if (run/'record.json').exists() else None
                if row is None:raise RuntimeError('Unreceipted seed attempt requires review: '+str(run))
                if row.get('qualified'):records.append(row)
                else:failures.append(row)
                continue
            run.mkdir();one=Chem.Mol(mol);one.RemoveAllConformers();one.AddConformer(mol.GetConformer(cid),assignId=True)
            (run/'input.xyz').write_text(to_xyz(one,'Exploratory constrained FF pose; constraints removed before xTB'))
            row=dict(mode=mode,seed=i,qualified=False,pi_groups=groups)
            try:
                row.update(seed.xtb(run,case['smiles'],groups))
                stereo=geometry_stereo_status(case['smiles'],row['xyz'])
                if stereo not in ('ok','not_applicable'):raise ValueError('Stereo gate: '+stereo)
                task=dict(case,key=case['key'],pi_groups=groups,model_type='covalent_oligomer')
                review=audit.audit(task,row['xyz'])
                if review['issues']:raise ValueError('Physical geometry gate: '+str(review['issues']))
                atomic_write_json(run/'geometry_review.json',review)
                row.update(qualified=True,xyz_sha256=sha(row['xyz']),stereochemistry=stereo)
                records.append(row)
            except Exception as exc:row['error']=repr(exc);failures.append(row)
            atomic_write_json(run/'record.json',row)
    if not records:
        atomic_write_json(directory/'hold.json',dict(case=case,failures=failures,status='no_qualified_conformer'))
        return dict(case=case,selected=[],records=[],failures=failures,status='preparation_held')
    # Pick the xTB energy leader and then a structurally distinct converged seed.
    # These are warm starts, not thermodynamic population weights.
    def rdkit_coordinates(path):
        observed=Chem.MolFromXYZFile(str(path));from rdkit.Chem import rdDetermineBonds
        rdDetermineBonds.DetermineConnectivity(observed,charge=0)
        from redox.build.geometry import _skeleton
        declared=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        mapping=_skeleton(observed).GetSubstructMatch(_skeleton(declared));assert len(mapping)==declared.GetNumAtoms()
        conf=Chem.Conformer(declared.GetNumAtoms())
        for a,b in enumerate(mapping):conf.SetAtomPosition(a,observed.GetConformer().GetAtomPosition(b))
        declared.AddConformer(conf);return Chem.RemoveHs(declared)
    ranked=sorted(records,key=lambda r:r['E_eV']);selected=[];molecules=[]
    for row in ranked:
        molecule=rdkit_coordinates(row['xyz'])
        if molecules and min(rdMolAlign.GetBestRMS(Chem.Mol(molecule),other,maxMatches=10000) for other in molecules)<.5:continue
        selected.append(row);molecules.append(molecule)
        if len(selected)==sampling['DFT_poses_per_graph']:break
    atoms=read(selected[0]['xyz']);qmol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),basis=ACTIVE_SP['basis'],charge=0,spin=0,unit='Angstrom',verbose=0)
    result=dict(case=case,records=records,selected=selected,failures=failures,status='conformers_prepared_not_DFT',
                active_basis_functions=qmol.nao_nr(),conformer_population_qualified=False)
    atomic_write_json(receipt,result);return result


def prepare(config,out,workers):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    if workers<1:raise ValueError('Positive CPU worker count required')
    spec=json.loads(config.read_text());out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'plan.json').exists():raise RuntimeError('Existing frozen QM plan must not be regenerated')
    feed=json.loads((ROOT/'results/candidates/current.json').read_text());motif=next(r for r in feed['original39']+feed['additional_candidates'] if r['id']==spec['motif_id'])
    cases,cycles=panel(spec,motif['frag'])
    if spec['quantum']['solvent']!='acetonitrile':raise ValueError('This validated backend currently admits MeCN only')
    snapshot=dict(config=spec,motif_id=motif['id'],fragment=motif['frag'],source_snapshot=feed['snapshot_path'],cases=cases,cycles=cycles,
                  scope=spec['scope'],status='graphs_prepared',experimental_loading_known=False)
    case_path=out/'cases.json'
    if case_path.exists():assert json.loads(case_path.read_text())==snapshot,'Existing exact graph specification changed'
    else:atomic_write_json(case_path,snapshot)
    for name in ('claims','done','failures','tasks','workers','launches','pilot','live_inputs'):(out/name).mkdir(exist_ok=True)
    with allocation_guard():
        resources=host_capacity(workers,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(workers))
    atomic_write_json(out/'preparation_owner.json',dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],cpus=sorted(os.sched_getaffinity(0)),resources=resources))
    ensembles=[]
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
        futures={pool.submit(sample,(str(out),case,motif['frag'],spec['sampling'])):case for case in cases}
        for future in as_completed(futures):
            ensemble=future.result();ensembles.append(ensemble)
            atomic_write_json(out/'preparation_progress.json',dict(graphs=len(cases),finished_graphs=len(ensembles),held=sum(not e['selected'] for e in ensembles)))
            print('Prepared',ensemble['case']['key'],len(ensemble['selected']),'qualified warm starts',flush=True)
    tasks=[]
    for ensemble in sorted(ensembles,key=lambda e:e['case']['key']):
        for i,row in enumerate(ensemble['selected']):
            case=ensemble['case'];key=case['key']+f'__pose{i}';directory=out/'live_inputs'/key;directory.mkdir(exist_ok=False)
            shutil.copyfile(row['xyz'],directory/'input.xyz')
            task=dict(case,key=key,candidate_id=spec['motif_id'],model_key=case['key'],model_type='covalent_oligomer',
                      repeats=case['repeat_units'],stereoisomer='common_mapped_backbone',target=f'pose{i}',fragment=motif['frag'],
                      pi_groups=row['pi_groups'],mult=1,role=spec['role'],active_basis_functions=ensemble['active_basis_functions'])
            tasks.append(task)
            atomic_write_json(directory/'task.json',dict(task,input_xyz=str(directory/'input.xyz'),
                input_sha256={str(directory/'input.xyz'):sha(directory/'input.xyz'),str(out/'ensembles'/case['key']/'ensemble.json'):sha(out/'ensembles'/case['key']/'ensemble.json')}))
    if not tasks:raise RuntimeError('No physically qualified inputs; inspect preparation holds')
    inputs=[config,case_path]+list((out/'live_inputs').glob('*/task.json'))
    sources=[config,Path(__file__),ROOT/'src/redox/build/polymer_cases.py',ROOT/'src/redox/build/polymer_contacts.py',ROOT/'scripts/validation/audit_polymer_geometries.py',ROOT/'scripts/pipeline/prepare_pi_stacking.py']
    atomic_write_json(out/'preparation.json',dict(schema=1,created=time.time(),tasks=tasks,cycles=cycles,graph_cases=cases,
        raw_root=str(ROOT/'calcs/polymer_models'/out.name),solvent='acetonitrile',
        pilot_key=max(tasks,key=lambda t:t['active_basis_functions'])['key'],maximum_memory_pilot_MiB=80*1024,
        producer_sha256={str(p):sha(p) for p in set(sources)},input_sha256={str(p):sha(p) for p in inputs},
        quantum_settings=spec['quantum'],electronic_protocol='Production r2SCAN-D4/SVP geometry and wB97M-V/TZVPD/SMD energy',
        thermochemistry=spec['thermochemistry'],grafting=spec['grafting'],motif_id=spec['motif_id'],
        property_scope=spec['scope'],prepared_graphs=sum(bool(e['selected']) for e in ensembles),held_graphs=[e['case']['key'] for e in ensembles if not e['selected']],
        maximum_attempts_per_item=1,consecutive_failure_limit=3,Hessians=False,
        minimum_large_GPU_MiB=80*1024,experimental_loading_known=False))
    print('Prepared',len(tasks),'DFT inputs from',len(cases),'matched graph cases; freeze and memory pilot required',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/'config/polymer_spacing_pilot.json')
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();prepare(args.config.resolve(),args.out,args.workers)
