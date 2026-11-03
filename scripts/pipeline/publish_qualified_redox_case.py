#!/usr/bin/env python3
"""Publish a checked exact redox case without filling unresolved material axes."""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.screening.qualified_cases import assemble


def publish(redox_path, registry_path, policy_path, output):
    if output.exists():raise FileExistsError('Preserve previous candidate publications')
    paths=dict(redox=redox_path,registry=registry_path,policy=policy_path)
    data={key:json.loads(path.read_text()) for key,path in paths.items()}
    bindings={key:dict(original_path=str(path),sha256=sha256_file(path),snapshot='inputs/'+key+'.json')
              for key,path in paths.items()}
    # Recheck the original independent-review backing; no fresh QM is implied.
    if any(sha256_file(path)!=digest for path,digest in data['redox']['source_sha256'].items()):
        raise ValueError('Qualified redox evidence backing changed')
    candidate=assemble(data['redox'],data['registry'],data['policy'],bindings)
    output.mkdir(parents=True);(output/'inputs').mkdir()
    for key,path in paths.items():
        target=output/bindings[key]['snapshot'];shutil.copyfile(path,target)
        if sha256_file(target)!=bindings[key]['sha256']:raise ValueError('Copied input differs')
    candidate['publisher_sha256']=sha256_file(__file__)
    (output/'candidate.json').write_text(json.dumps(candidate,indent=2)+'\n')
    return candidate


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('redox','registry','policy','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=publish(a.redox.resolve(),a.registry.resolve(),a.policy.resolve(),a.out.resolve())
    print(json.dumps(dict(candidate=result['candidate_alias'],first_oxidation=result['properties']['first_oxidation']['value'],
                          all_axis_status=result['all_axis_decision']['status'])))
