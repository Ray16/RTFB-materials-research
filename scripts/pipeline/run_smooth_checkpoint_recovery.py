#!/usr/bin/env python
"""Run the unchanged smooth geometry worker with bounded checkpoint-lock retries."""
import argparse
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoint_io import checkpoint_lock_retries

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('run',));p.add_argument('--campaign',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    worker=api.load_module('same_model_geometry_checkpoint_recovery',ROOT/'scripts/pipeline/run_smooth_geometry.py')
    with checkpoint_lock_retries(a.campaign.resolve()/'checkpoint_lock_retries.jsonl'):
        worker.run(a.campaign.resolve())
