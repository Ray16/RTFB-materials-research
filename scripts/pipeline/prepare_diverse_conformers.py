#!/usr/bin/env python
"""Prepare energy-prioritized, geometry-diverse warm starts from saved searches.

Completed and stopped searches remain distinct. Selected input geometries are
validated; they are not DFT minima, population weights or converged ensembles.
Earlier input plans and calculations remain immutable.
"""
import argparse
import inspect
import json
import math
import os
from pathlib import Path
import tempfile

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def records(entry,out,stager,cutoff,maximum_pool):
    from rdkit.Chem import rdMolAlign
    helper=api.load_module('diverse_input_chemical_mapping',ROOT/'scripts/pipeline/polymer_conformer_thermals.py')
    selector=api.load_module('bounded_saved_conformer_pool',ROOT/'scripts/pipeline/stage_interrupted_conformers.py')
    reader=api.load_module('diverse_strict_xyz_frames',ROOT/'scripts/pipeline/sample_polymer_crest.py')
    owner=api.load_module('diverse_search_terminal_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    raw=Path(entry['source_search']);status_path=raw/'status.json';status=json.loads(status_path.read_text())
    allowed={'low_level_ensemble_complete_DFT_pending','failed_or_interrupted_preserving_outputs'}
    if status['status'] not in allowed or not status.get('finished') or owner.owner_alive(status):
        raise ValueError('A known completed or walltime-stopped source must have ended')
    if status['status'].startswith('failed') and 'TimeoutExpired' not in status.get('error',''):
        raise ValueError('Other failed searches require separate source review')
    api.verify_sources(status['source_sha256'])
    case_path=Path(entry['case_path']);case=json.loads(case_path.read_text())
    if case!=status['case'] or (case['charge'],case['mult'])!=(0,1):raise ValueError('Exact source graph/state differs')
    argv=status['argv']
    for option,value in (('--alpb','acetonitrile'),('--chrg','0'),('--uhf','0')):
        if option not in argv or argv[argv.index(option)+1]!=value:raise ValueError('Source medium/state differs')
    ensemble=raw/'crest_conformers.xyz'
    ranked,count=selector.lowest_frames(ensemble,reader.xyz_frames,stager.plain_xyz_energy,maximum_pool)
    sources=dict(status['source_sha256'])
    sources.update({str(p):api.digest(p) for p in (status_path,ensemble,case_path,Path(__file__).resolve())})
    molecules=[];maps=None;references=[]
    for ref in entry.get('mandatory_references',[]):
        if api.digest(ref['path'])!=ref['sha256']:raise ValueError('Reference input geometry changed')
        molecules.append(helper.alignment_molecule(case['smiles'],Path(ref['path'])))
        references.append(ref['key'])
    chosen=[];rejected=[];lowest=None;tested=0
    target=20-len(references)
    saved=out/'selected_source_frames'/entry['key'];saved.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='redox_diverse_seed_') as scratch:
        path=Path(scratch)/'candidate.xyz'
        for energy,index,frame in ranked:
            if lowest is not None and (energy-lowest)*627.509451>8.:break
            path.write_text(frame);tested+=1
            try:mol=helper.alignment_molecule(case['smiles'],path)
            except (ValueError,RuntimeError):
                rejected.append(dict(index=index,reason='exact graph/stereo mapping failed'));continue
            if maps is None:
                matches=mol.GetSubstructMatches(mol,uniquify=False,useChirality=True,maxMatches=10001)
                if not matches or len(matches)>=10001:raise ValueError('Symmetry enumeration needs review')
                maps=[list(enumerate(m)) for m in matches]
            minimum=min((rdMolAlign.GetBestAlignmentTransform(mol,old,map=maps,reflect=False,numThreads=1)[0]
                         for old in molecules),default=float('inf'))
            # Retain the actual lowest eligible low-level structure, even when
            # close to an inherited reference. Additional starts must be diverse.
            first=lowest is None
            if first:lowest=energy
            if not first and minimum<cutoff:continue
            dest=saved/f'conf_{index:05d}.xyz';dest.write_text(frame)
            chosen.append(dict(index=index,path=str(dest),sha256=api.digest(dest),GFN2_energy_Ha=energy,
                source_graph='ok',source_stereo='ok',minimum_RMSD_to_prior_inputs_A=None if math.isinf(minimum) else float(minimum)))
            molecules.append(mol)
            if len(chosen)==target:break
    if len(chosen)!=target:
        atomic_write_json(saved/'insufficient_diversity.json',dict(selected=len(chosen),required=target,tested=tested,source_frames=count))
        raise ValueError('Source pool does not supply the declared diverse twenty-start set')
    snapshot=dict(status='diverse_input_geometries_checked_not_DFT_minima',source_frames=count,
        ranked_pool_size=len(ranked),tested_candidates=tested,energy_window_kcal=8.,
        additional_start_RMSD_cutoff_A=cutoff,known_references=references,selected=chosen,rejected=rejected,
        complete_source_search=status['status']=='low_level_ensemble_complete_DFT_pending',
        source_search_status=status['status'],source_sha256=sources,thermodynamic_weights_from_hit_counts=False)
    atomic_write_json(saved/'selection.json',snapshot);sources[str(saved/'selection.json')]=api.digest(saved/'selection.json')
    return case,chosen,rejected,sources


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--rmsd-cutoff',type=float,default=.5);parser.add_argument('--maximum-pool',type=int,default=2000)
    args=parser.parse_args()
    if args.rmsd_cutoff!=.5 or args.maximum_pool!=2000:parser.error('Use the declared 0.5 Å / 2000-frame input pilot')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(1))
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source_root.resolve();request=json.loads((source/'request.json').read_text())
    stager_path=ROOT/'scripts/pipeline/prepare_crest_refinement.py'
    stager=api.load_module('retained_twenty_start_input_builder',stager_path)
    request['raw_root']=str(ROOT/'calcs/polymer_models'/out.name)
    request['source_sha256'].update({str(p):api.digest(p) for p in (Path(__file__).resolve(),source/'request.json',stager_path)})
    request['input_selection']='Lowest eligible source plus energy-prioritized separated inputs and inherited reference anchors'
    atomic_write_json(out/'request.json',request)
    stager.source_records=lambda entry:records(entry,out,stager,args.rmsd_cutoff,args.maximum_pool)
    body=inspect.getsource(stager.stage_case).replace('final_CREST_conformer','diverse_saved_search_warm_start').replace('final_CREST_index_','diverse_saved_search_index_')
    namespace=dict(vars(stager));local={}
    exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local)
    rows=[]
    for entry in request['cases']:
        try:
            row=local['stage_case'](out,request,entry);case=out/entry['key'];plan=api.read_plan(case)
            selection=json.loads((out/'selected_source_frames'/entry['key']/'selection.json').read_text())
            plan.update(source_search_complete=selection['complete_source_search'],
                input_selection=selection,input_geometry_diversity_is_not_ensemble_convergence=True)
            plan['producer_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
            atomic_write_json(case/'plan.json',plan);rows.append(row)
        except (ValueError,RuntimeError,OSError) as error:
            rows.append(dict(key=entry['key'],status='input_preparation_requires_review_preserved',error=repr(error)))
        atomic_write_json(out/'preparation_status.json',dict(cases=rows,quantum_jobs_submitted=0,
            converged_ensemble_claimed=False,input_strategy=request['input_selection']))
    print(json.dumps(dict(prepared=sum(r.get('starts',0)==20 for r in rows),cases=rows)))


if __name__=='__main__':main()
