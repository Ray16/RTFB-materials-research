#!/usr/bin/env python
"""Freeze exact-graph DFT warm starts from stopped CREST checkpoints.

These are incomplete low-level searches, not final ensembles or populations.
Every selected frame receives graph/stereochemistry QC and later full DFT QC.
Original searches and their terminal failures remain unchanged.
"""
import argparse
import heapq
import inspect
import json
import math
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def lowest_frames(path, reader, energy, maximum=500):
    selected=[]; count=0
    for index,frame in enumerate(reader(path)):
        value=energy(frame)
        if not math.isfinite(value):
            raise ValueError('Checkpoint conformer energy is nonfinite')
        count+=1
        item=(-value,-index,frame)
        if len(selected)<maximum:heapq.heappush(selected,item)
        elif item>selected[0]:heapq.heapreplace(selected,item)
    if count==0:raise ValueError('Empty conformer checkpoint')
    return sorted([(-e,-i,frame) for e,i,frame in selected]),count


def checkpoint_records(entry, root, stager):
    from redox.build.geometry import geometry_stereo_status
    from redox.validation.protonation import hydrogen_topology_status
    from redox.build.polymer_cases import build_case
    raw=Path(entry['source_search']);status_path=raw/'status.json'
    status=json.loads(status_path.read_text())
    owner=api.load_module('interrupted_CREST_terminal_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if (status['status']!='failed_or_interrupted_preserving_outputs' or not status.get('finished')
            or 'TimeoutExpired' not in status.get('error','') or owner.owner_alive(status)):
        raise ValueError('A confirmed stopped walltime-limited search is required')
    api.verify_sources(status['source_sha256'])
    case_path=Path(entry['case_path']);case=json.loads(case_path.read_text())
    if case!=status['case'] or (case['charge'],case['mult'])!=(0,1):
        raise ValueError('Checkpoint search state/graph differs from requested case')
    rebuilt=build_case(entry['fragment'],case['repeat_units'],case['reactive_sites'],case['occupied_sites'],case['backbone_stereo_tags'])
    if rebuilt['smiles']!=case['smiles'] or rebuilt['backbone_stereo_fingerprint']!=case['backbone_stereo_fingerprint']:
        raise ValueError('Checkpoint is not the exact open-chain occupancy/configuration graph')
    argv=status['argv']
    for option,value in (('--alpb','acetonitrile'),('--chrg','0'),('--uhf','0')):
        if option not in argv or argv[argv.index(option)+1]!=value:
            raise ValueError('Checkpoint low-level medium/charge/spin differs')
    if '--gfn2' not in argv:raise ValueError('Checkpoint is not GFN2')
    ensemble=raw/'crest_conformers.xyz';log=raw/'crest.log'
    if 'MTD Simulations done' not in log.read_text(errors='replace'):
        raise ValueError('Required completed metadynamics phase is not recorded')
    reader=api.load_module('checkpoint_strict_XYZ_reader',ROOT/'scripts/pipeline/sample_polymer_crest.py')
    ranked,count=lowest_frames(ensemble,reader.xyz_frames,stager.plain_xyz_energy)
    saved=root/'checkpoint_inputs'/entry['key'];saved.mkdir(parents=True,exist_ok=False)
    good=[];bad=[]
    for value,index,frame in ranked:
        path=saved/f'conf_{index:05d}.xyz';path.write_text(frame)
        graph=hydrogen_topology_status(case['smiles'],path)
        stereo=geometry_stereo_status(case['smiles'],path)
        row=dict(index=index,path=str(path),sha256=api.digest(path),GFN2_energy_Ha=value,graph=graph,stereo=stereo)
        if graph=='ok' and stereo in ('ok','not_applicable'):good.append(row)
        else:bad.append(row)
        if len(good)>=20:break
    if len(good)<20:raise ValueError('Fewer than twenty exact-graph checkpoint starts survive QC')
    sources=dict(status['source_sha256'])
    sources.update({str(p):api.digest(p) for p in (status_path,ensemble,log,case_path,Path(__file__).resolve())})
    atomic_write_json(saved/'snapshot.json',dict(status='interrupted_search_warm_starts_checked',
        source_search_complete=False,metadynamics_phase_reported_complete=True,total_complete_frames=count,
        ranked_candidates_checked=len(good)+len(bad),selected=good,rejected=bad,source_sha256=sources,
        DFT_minima_or_ensemble_G_qualified=False,search_hit_counts_are_not_population_weights=True))
    sources[str(saved/'snapshot.json')]=api.digest(saved/'snapshot.json')
    return case,good,bad,sources


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source_root.resolve();request=json.loads((source/'request.json').read_text())
    stager_path=ROOT/'scripts/pipeline/prepare_crest_refinement.py'
    stager=api.load_module('unchanged_conformer_DFT_input_builder',stager_path)
    request['cases']=[e for e in request['cases'] if json.loads((Path(e['source_search'])/'status.json').read_text())['status']=='failed_or_interrupted_preserving_outputs']
    request['raw_root']=str(ROOT/'calcs/polymer_models'/out.name)
    request['source_sha256'].update({str(p):api.digest(p) for p in (Path(__file__).resolve(),stager_path,source/'request.json')})
    request['source_searches_complete']=False
    atomic_write_json(out/'request.json',request)
    stager.source_records=lambda entry:checkpoint_records(entry,out,stager)
    body=inspect.getsource(stager.stage_case).replace('final_CREST_conformer','interrupted_CREST_checkpoint_warm_start').replace('final_CREST_index_','checkpoint_CREST_index_')
    namespace=dict(vars(stager));local={}
    exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local)
    rows=[]
    for entry in request['cases']:
        try:
            row=local['stage_case'](out,request,entry)
            campaign=out/entry['key'];plan=api.read_plan(campaign)
            plan['producer_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
            plan.update(source_search_complete=False,initial_ensemble_role='stopped-search DFT warm starts, not an equilibrium ensemble')
            atomic_write_json(campaign/'plan.json',plan)
            rows.append(row)
        except (ValueError,RuntimeError,OSError) as error:
            rows.append(dict(key=entry['key'],status='preparation_requires_review_preserved',error=repr(error)))
        atomic_write_json(out/'preparation_status.json',dict(cases=rows,quantum_jobs_submitted=0,
            source_searches_complete=False,broader_conformer_convergence_established=False))
    print(json.dumps(dict(prepared=sum(r.get('starts',0)>0 for r in rows),cases=rows)))


if __name__=='__main__':
    main()
