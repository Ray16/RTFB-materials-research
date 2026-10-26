#!/usr/bin/env python
"""Use plan-declared local-energy RAM guards without changing the calculator."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def adapted_launch(worker):
    body=inspect.getsource(worker._launch)
    replacements={
        'minimum_ram_gib=64':"minimum_ram_gib=plan['maximum_host_RSS_GiB']+32",
        'bound=32*2**30':"bound=plan['maximum_host_RSS_GiB']*2**30",
        "'--maximum-rss-gib','32'":"'--maximum-rss-gib',str(plan['maximum_host_RSS_GiB'])"}
    for old,new in replacements.items():
        if body.count(old)!=1:raise ValueError('Declared local-resource placement template changed')
        body=body.replace(old,new,1)
    capacity="        capacity=host_capacity(2,minimum_ram_gib=plan['maximum_host_RSS_GiB']+32);promise=placement.remaining_promises()"
    if body.count(capacity)!=1:raise ValueError('Local capacity template changed')
    body=body.replace(capacity,"""        try:
            capacity=host_capacity(2,minimum_ram_gib=plan['maximum_host_RSS_GiB']+32)
        except RuntimeError as error:
            return dict(status='resource_capacity_unavailable',reason=str(error))
        promise=placement.remaining_promises()""",1)
    namespace={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(worker)),namespace)
    return namespace['_launch']


def main():
    worker=api.load_module('unchanged_local_electronic_calculator',ROOT/'scripts/pipeline/local_motif_energy.py')
    worker._launch=adapted_launch(worker);worker.main()


if __name__=='__main__':main()
