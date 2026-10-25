#!/usr/bin/env python
"""Reuse reservation/preflight rules for two experimental-surface controls."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    source=ROOT/'scripts/validation/launch_response_controls.py'
    module=api.load_module('qualified_small_control_placement',source)
    text=inspect.getsource(module.main)
    replacements={
        'scripts/validation/qualify_residual_response.py':'scripts/validation/qualify_smooth_quantum.py',
        'scripts/pipeline/residual_response_solver.py':'scripts/pipeline/smooth_cds_control.py',
        'redox-residual-response-full-control':'redox-experimental-smooth-CDS-control',
        'Exactly one full singlet and one full doublet replay, 14 GiB device ceiling each':
            'Two distinct-model full quantum controls; force/Hessian/minimum gates, 14 GiB ceiling each'}
    for old,new in replacements.items():
        if text.count(old)!=1:raise ValueError('Small-control placement template changed')
        text=text.replace(old,new,1)
    namespace=dict(vars(module),__file__=str(Path(__file__).resolve()));local={}
    exec(compile(text,str(Path(__file__).resolve()),'exec'),namespace,local)
    local['main']()


if __name__=='__main__':
    main()
