#!/usr/bin/env python
"""Assemble five checked terms with explicit per-role attempt locations.

Only the record locations differ from the existing collector. Method, common
grid, checkpoint identity and independent stationarity gates remain unchanged.
"""
import argparse
import inspect
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

ORIGINAL=ROOT/'scripts/validation/collect_counterpoise_cpu_terms.py'
ROLES={'AB','B_own','B_ghost'}


def locations(data):
    terms=data.get('remaining_terms',{})
    if data.get('schema')!=1 or set(terms)!=ROLES:
        raise ValueError('Three exact remaining-term locations required')
    raw=[]
    for role,item in terms.items():
        if set(item)!={'raw','review','stationarity'}:raise ValueError('Raw, checkpoint review and stationarity locations required')
        for field,path in item.items():
            p=Path(path).resolve()
            if not p.is_relative_to(ROOT/'calcs'):raise ValueError('Term records must remain on durable project storage')
        raw.append(str(Path(item['raw']).resolve()))
    if len(set(raw))!=3:raise ValueError('Distinct per-role attempt directories required')
    return terms


def adapted_collect(module,terms,proof):
    body=inspect.getsource(module.collect)
    old="        raw=remaining_root/role;state_path=raw/role/'energy.json';launch_path=raw/'launch.json';status_path=raw/'status.json'"
    new="        raw=Path(_term_locations[role]['raw']);state_path=raw/role/'energy.json';launch_path=raw/'launch.json';status_path=raw/'status.json'"
    if body.count(old)!=1:raise ValueError('Term-location template changed')
    body=body.replace(old,new,1)
    old="        review_path=reviews_root/role/'review.json';gradient_path=reviews_root/role/'stationarity/summary.json'"
    new="        review_path=Path(_term_locations[role]['review']);gradient_path=Path(_term_locations[role]['stationarity'])"
    if body.count(old)!=1:raise ValueError('Independent-review location template changed')
    body=body.replace(old,new,1)
    old="    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result"
    new="""    result['source_sha256'].update(_manifest_proof)
    result['term_record_locations']=_term_locations
    api.verify_sources(_manifest_proof)
"""+old
    if body.count(old)!=1:raise ValueError('Counterpoise output contract changed')
    namespace=dict(vars(module),_term_locations=terms,_manifest_proof=proof);local={}
    exec(compile(body.replace(old,new,1),__file__,'exec'),namespace,local)
    return local['collect']


def collect(manifest,out):
    data=json.loads(manifest.read_text());terms=locations(data)
    module=api.load_module('unchanged_counterpoise_gates_explicit_locations',ORIGINAL)
    proof={str(p):api.digest(p) for p in (manifest,Path(__file__).resolve(),ORIGINAL)}
    function=adapted_collect(module,terms,proof)
    paths=[]
    for field in ('case','geometry_review','own','own_stationarity','ghost','ghost_stationarity'):
        p=Path(data[field]).resolve()
        if not p.is_relative_to(ROOT):raise ValueError('Shared source inputs required')
        paths.append(p)
    return function(*paths,manifest.parent,manifest.parent,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    result=collect(a.manifest.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],sums=result['sums'])))
