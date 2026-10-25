#!/usr/bin/env python
"""Build the PCM surface before DF storage selection in the fixed-state probe.

The original SCF driver builds in this order. No operator, orbital, grid or
threshold changes; the DF implementation chooses storage from actual free RAM.
The first diagnostic's allocation failure remains in its original directory.
"""
import inspect
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def reordered_source(body):
    old = '            mf.with_df.build()\n            mf.with_solvent._get_vind(mf.make_rdm1())'
    new = '            mf.with_solvent._get_vind(mf.make_rdm1())\n            dft._free_gpu_pool()\n            mf.with_df.build()'
    if body.count(old) != 1:
        raise ValueError('Saved-orbital probe construction template changed')
    return body.replace(old,new,1)


def main():
    source=ROOT/'scripts/validation/probe_fixed_orbital_response.py'
    module=api.load_module('preserved_fixed_state_response_probe',source)
    if sys.argv[1]=='run':
        namespace=dict(vars(module));local={}
        exec(compile(reordered_source(inspect.getsource(module.run)),str(Path(__file__).resolve()),'exec'),namespace,local)
        module.run=local['run']
    module.main()
    if sys.argv[1]=='prepare':
        campaign=Path(sys.argv[sys.argv.index('--campaign')+1]).resolve()
        plan=api.read_plan(campaign)
        plan['execution_driver']=str(Path(__file__).resolve())
        plan['producer_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
        plan['construction_order']='PCM surface, release unreachable GPU pool, DF storage selection'
        atomic_write_json(campaign/'plan.json',plan)


if __name__=='__main__':
    main()
