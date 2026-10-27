#!/usr/bin/env python
"""Test whether Nitro FP16 intermediates explain the retained repeat failure.

This is a separate implementation diagnostic using the installed executor's
explicit precision switches. No production profile or tolerance is changed.
"""
import json
import os
from pathlib import Path
import runpy
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

FLAGS=('FP16_RADIAL','FP16_MERGED_BWD','FP16_LN1','FP16_XEGRAD','FP16_GRADBLK',
       'FP16_XEDGE','FP16_GRAD_RAD','FP16_DH','FP16_DX','FP16_M0RAD')


def main():
    if any(name.startswith('fairchem') for name in sys.modules):
        raise RuntimeError('Set the precision control before importing the executor')
    out=Path(sys.argv[sys.argv.index('--out')+1]).resolve()
    source=ROOT/'scripts/validation/qualify_uma_current_diagnostic.py'
    flags={'UMA_NITRO_'+name:'0' for name in FLAGS}
    previous={name:os.environ.get(name) for name in flags}
    for name,value in flags.items():os.environ[name]=value
    record=dict(status='Nitro_FP32_intermediate_diagnostic_requested',started=time.time(),
                environment=flags,production_profile_changed=False,full_panel_qualified=False,
                source_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(),source)})
    try:
        runpy.run_path(str(source),run_name='__main__')
    finally:
        executors=[m for name,m in sys.modules.items() if name.endswith('uma_nitro._inference')]
        if len(executors)==1:
            executor=executors[0]
            record['effective_precision_flags']={name:getattr(executor,name,None) for name in FLAGS}
            record['all_requested_controls_observed']=all(v is False for v in record['effective_precision_flags'].values())
            record['executor_source_sha256']={str(Path(executor.__file__).resolve()):api.digest(executor.__file__)}
        else:record['all_requested_controls_observed']=False
        if out.exists():
            record['finished']=time.time();atomic_write_json(out/'precision_control.json',record)
        for name,value in previous.items():
            if value is None:os.environ.pop(name,None)
            else:os.environ[name]=value


if __name__=='__main__':main()
