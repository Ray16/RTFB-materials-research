#!/usr/bin/env python
"""Prepare all missing five-repeat occupancy/tacticity searches reproducibly."""
import argparse,hashlib,json,os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.polymer_cases import build_case
from redox.build.structures import to_xyz
from redox.build.geometry import geometry_stereo_status
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False);source=ROOT/'results/reviews/2026-10-22/spacing_native_complete/validation.json'
    original=json.loads(source.read_text());cases={tuple(r['occupied_sites']):r for r in original['cases']};seed=20261023;tasks=[]
    for configuration,tags in [('all_same',['CW']*4),('alternating',['CW','CCW','CW','CCW'])]:
        for occupied in ((2,),(3,),(4,),(2,3),(2,4)):
            case=build_case('[*:1]N1C(=O)c2ccccc2C1=O',5,[2,3,4],list(occupied),tags)
            if configuration=='all_same' and case['smiles']!=cases[occupied]['smiles']:raise ValueError('Reconstructed graph differs from the checked original')
            key=configuration+'__occ_'+'_'.join(map(str,occupied));directory=args.out/key;directory.mkdir()
            case.update(key=key,mult=1,role='Neutral phthalimide method control; not synthesis admission',initial_embedding_seed=seed)
            mol=Chem.AddHs(Chem.MolFromSmiles(case['smiles']));parameters=AllChem.ETKDGv3();parameters.randomSeed=seed;parameters.numThreads=1
            confs=list(AllChem.EmbedMultipleConfs(mol,numConfs=32,params=parameters))
            if not confs:raise ValueError('No embedded structure for '+key)
            optimized=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=2000)
            good=[(energy,cid) for cid,(status,energy) in zip(confs,optimized) if status==0]
            if not good:raise ValueError('No converged MMFF input for '+key)
            energy,cid=min(good);one=Chem.Mol(mol);one.RemoveAllConformers();one.AddConformer(mol.GetConformer(cid),assignId=True)
            xyz=directory/'input.xyz';xyz.write_text(to_xyz(one,'MMFF initial geometry only; broad CREST search follows'))
            if hydrogen_topology_status(case['smiles'],xyz)!='ok' or geometry_stereo_status(case['smiles'],xyz) not in ('ok','not_applicable'):
                raise ValueError('Prepared stereochemical graph mismatch')
            case.update(embedded=len(confs),converged_MMFF=len(good),initial_MMFF_energy_kcal_mol=energy)
            atomic_write_json(directory/'case.json',case)
            tasks.append(dict(key=key,case=str(directory/'case.json'),input=str(xyz),raw_out=str(ROOT/'calcs/polymer_models/2026-10-23_conformer_revision'/key),
                input_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),case_sha256=hashlib.sha256((directory/'case.json').read_bytes()).hexdigest()))
    atomic_write_json(args.out/'plan.json',dict(tasks=tasks,old_empty_backbone_searches=['calcs/polymer_models/2026-10-22_crest_all_same','calcs/polymer_models/2026-10-22_crest_alternating'],
        scope='Ten remaining neutral phthalimide graph/configuration searches; no atactic population model',
        low_level=dict(method='GFN2',solvent_model='ALPB',solvent='acetonitrile',energy_window_kcal_mol=8,threads=14,BLAS_threads=1,MD_seed='CREST default; trajectories retained'),
        DFT_followup=dict(status='not_yet_submitted',target_unique_starts_per_graph=20,initial_rerank_batch=10,expanded_batch=20,
            minimum_change_tolerance_kJ_mol=2,ensemble_change_tolerance_kJ_mol=2,
            limit='These within-search/batch criteria do not prove a global minimum or a complete conformer ensemble'),
        producer_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),source,ROOT/'src/redox/build/polymer_cases.py',ROOT/'scripts/pipeline/sample_polymer_conformers.py')}))
    print(json.dumps(dict(prepared=len(tasks),existing_reference_searches=2)))


if __name__=='__main__':main()
