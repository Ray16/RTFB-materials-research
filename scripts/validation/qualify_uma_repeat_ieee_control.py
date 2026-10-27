#!/usr/bin/env python
"""Process-local IEEE-matmul diagnostic for Nitro on the failed first case.

Adds a strict PyTorch matmul control to the FP16-off diagnostic. This changes an
installed internal precision policy only in this process and is not a qualified
production profile or an alteration of model weights or acceptance thresholds.
"""
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
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    if any(name.startswith('fairchem') for name in sys.modules):raise RuntimeError('Use a fresh diagnostic process')
    out=Path(sys.argv[sys.argv.index('--out')+1]).resolve();source=ROOT/'scripts/validation/qualify_uma_current_diagnostic.py'
    previous={name:os.environ.get('UMA_NITRO_'+name) for name in FLAGS}
    for name in FLAGS:os.environ['UMA_NITRO_'+name]='0'
    from fairchem.core.models.uma.uma_nitro import _inference as runtime
    if type(getattr(runtime,'_USE_TF32',None)) is not bool:raise ValueError('Installed precision policy changed')
    old=runtime._USE_TF32;runtime._USE_TF32=False
    record=dict(status='Nitro_internal_IEEE_precision_diagnostic_requested',started=time.time(),original_USE_TF32=old,
                effective_USE_TF32=runtime._USE_TF32,effective_precision_flags={name:getattr(runtime,name,None) for name in FLAGS},
                official_default_profile=False,production_profile_changed=False,full_panel_qualified=False,
                source_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(),source,Path(runtime.__file__).resolve())})
    if not all(v is False for v in record['effective_precision_flags'].values()):raise ValueError('FP16 control flags did not take effect')
    try:runpy.run_path(str(source),run_name='__main__')
    finally:
        record['finished']=time.time()
        if out.exists():atomic_write_json(out/'precision_control.json',record)
        runtime._USE_TF32=old
        for name,value in previous.items():
            if value is None:os.environ.pop('UMA_NITRO_'+name,None)
            else:os.environ['UMA_NITRO_'+name]=value


if __name__=='__main__':main()
