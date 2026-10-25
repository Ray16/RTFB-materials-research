#!/usr/bin/env python
"""Use qualified resource placement for the two additional N/Cl controls."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    source=ROOT/'scripts/validation/launch_response_controls.py'
    module=api.load_module('NCl_control_qualified_placement',source)
    body=inspect.getsource(module.main)
    replacements={
        'scripts/validation/qualify_residual_response.py':'scripts/validation/qualify_smooth_NCl.py',
        'scripts/pipeline/residual_response_solver.py':'scripts/pipeline/smooth_cds_control.py',
        'redox-residual-response-full-control':'redox-smooth-CDS-NCl-controls',
        'Exactly one full singlet and one full doublet replay, 14 GiB device ceiling each':
            'Two full N/Cl model controls, one singlet and one doublet; 14 GiB device ceiling'}
    for old,new in replacements.items():
        if body.count(old)!=1:raise ValueError('Small control placement template changed')
        body=body.replace(old,new,1)
    namespace=dict(vars(module),__file__=str(Path(__file__).resolve()));local={}
    exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local)
    local['main']()


if __name__=='__main__':main()
