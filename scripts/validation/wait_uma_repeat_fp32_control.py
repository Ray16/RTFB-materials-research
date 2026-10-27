#!/usr/bin/env python
"""Reserve one fixed-device precision diagnostic using the existing waiter."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    upstream=ROOT/'scripts/validation/wait_uma_pair.py'
    diagnostic=ROOT/'scripts/validation/qualify_uma_repeat_fp32_control.py'
    module=api.load_module('UMA_precision_control_reserved_wait',upstream)
    body=inspect.getsource(module.main)
    changes={"script=ROOT/'scripts/validation/qualify_uma_current.py'":"script=ROOT/'scripts/validation/qualify_uma_repeat_fp32_control.py'",
             "helper=ROOT.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'":"helper=ROOT/'scripts/fleet/preflight_cpu_dependencies.py'"}
    for old,new in changes.items():
        if body.count(old)!=1:raise ValueError('Reserved diagnostic template changed')
        body=body.replace(old,new,1)
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),upstream,diagnostic,
        ROOT/'scripts/validation/qualify_uma_current_diagnostic.py',ROOT/'scripts/fleet/preflight_cpu_dependencies.py')}
    verify=api.verify_sources
    def bound(mapping):verify(proof);verify(mapping)
    api.verify_sources=bound
    local={}
    try:
        exec(compile(body,__file__,'exec'),dict(vars(module),__file__=__file__),local)
        local['main']()
    finally:api.verify_sources=verify


if __name__=='__main__':main()
