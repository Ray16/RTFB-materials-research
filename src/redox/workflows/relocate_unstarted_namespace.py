"""Relocate byte-identical dependency bindings in a never-started namespace.

No code/geometry/state/method bytes are changed. Preserve the original namespace,
prelaunch refusals and owner records; require every recorded controller closed.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def relocate(namespace,output,mapping):
    spec=verify(namespace);old=namespace.parent;case=Path(spec['case'])
    if output.exists() or any(output.is_relative_to(p.resolve()) for p in (ROOT,Path(spec['archive']),old)):
        raise ValueError('A fresh destination outside the scientific/source namespaces is required')
    if (case/'status.json').exists() or (case/'residency_guard').exists() or any(old.glob('placement_*/launch.json')):
        raise ValueError('Only a namespace with no scientific start can be relocated')
    definition=importlib.util.spec_from_file_location('relocation_closed_owners',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    for receipt in (old/'dispatch_claim').glob('try_*.json'):
        if owners.alive(json.loads(receipt.read_text())['owner']) is not False:
            raise ValueError('Every original controller must be proven closed')
    plan=json.loads((case/'plan.json').read_text());bindings={}
    for original,replacement in mapping.items():
        if original not in plan['source_sha256'] or sha256_file(replacement)!=plan['source_sha256'][original]:
            raise ValueError('Relocated dependency must match the original bound bytes exactly')
        bindings[original]=dict(replacement=str(Path(replacement).resolve()),sha256=plan['source_sha256'][original])
    def replace(value):
        if isinstance(value,str):
            if value in bindings:return bindings[value]['replacement']
            if value==str(old) or value.startswith(str(old)+'/'):return str(output)+value[len(str(old)):]
            return value
        if isinstance(value,list):return [replace(x) for x in value]
        if isinstance(value,dict):return {replace(k):replace(v) for k,v in value.items()}
        return value
    output.mkdir(parents=True);shutil.copytree(case,output/'case');shutil.copytree(Path(spec['sources']),output/'sources')
    preserved=dict(source_namespace=str(namespace),source_namespace_sha256=sha256_file(namespace),
        source_plan_sha256=spec['plan_sha256'],dependency_relocations=bindings,scientific_parameters_changed=False,scientific_jobs_submitted=0)
    atomic_write_json(output/'sources/dependency_relocation.json',preserved)
    shutil.copyfile(Path(__file__),output/'sources/relocation_preparer.py')
    plan=replace(plan);plan['source_sha256'].update({str(output/'sources'/name):sha256_file(output/'sources'/name) for name in ('dependency_relocation.json','relocation_preparer.py')})
    atomic_write_json(output/'case/plan.json',plan)
    spec=replace(spec);spec['plan_sha256']=sha256_file(output/'case/plan.json');atomic_write_json(output/'namespace.json',spec)
    verify(output/'namespace.json')
    return dict(status='unstarted_namespace_dependencies_explicitly_relocated',namespace=str(output/'namespace.json'),**preserved)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mapping',type=Path,required=True);a=p.parse_args()
    print(json.dumps(relocate(a.spec.resolve(),a.out.resolve(),json.loads(a.mapping.read_text())),indent=2))
