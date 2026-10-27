#!/usr/bin/env python
"""Place the first-case repeat diagnostic using the existing reservation contract."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

UPSTREAM=ROOT/'scripts/validation/wait_uma_pair.py'
DIAGNOSTIC=ROOT/'scripts/validation/qualify_uma_current_diagnostic.py'


def main():
    module=api.load_module('original_first_case_wait_contract',UPSTREAM)
    source=inspect.getsource(module.main)
    old="script=ROOT/'scripts/validation/qualify_uma_current.py'"
    new="script=ROOT/'scripts/validation/qualify_uma_current_diagnostic.py'"
    if source.count(old)!=1:raise ValueError('UMA reserved-worker launch template changed')
    source=source.replace(old,new,1)
    old="helper=ROOT.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'"
    new="helper=ROOT/'scripts/fleet/preflight_cpu_dependencies.py'"
    if source.count(old)!=1:raise ValueError('CPU preflight template changed')
    source=source.replace(old,new,1)
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),UPSTREAM,DIAGNOSTIC,
        ROOT/'scripts/fleet/preflight_cpu_dependencies.py')}
    original_verify=api.verify_sources
    def verify(mapping):
        original_verify(proof)
        original_verify(mapping)
    api.verify_sources=verify
    namespace=dict(vars(module),__file__=__file__);local={}
    try:
        exec(compile(source,__file__,'exec'),namespace,local)
        local['main']()
    finally:api.verify_sources=original_verify


if __name__=='__main__':main()
