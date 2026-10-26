#!/usr/bin/env python
"""Enforce prospective-device admission inside a reserved motif worker."""
import argparse
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def validate_device(plan,device):
    if (plan.get('minimum_GPU_MiB')!=32768 or device['memory_MiB']<32768
            or plan.get('primary_candidate_admission') or plan.get('solvent_counterpoise_computed')):
        raise ValueError('Motif counterpoise requires a reserved 32-GB-class device and gas-only scope')


def run(case):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    plan=api.read_plan(case);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if plan['execution_driver']!=str(Path(__file__).resolve()):raise ValueError('Unexpected point execution driver')
    device=api.load_module('counterpoise_inner_device_probe',ROOT/'scripts/pipeline/local_motif_energy.py')
    validate_device(plan,device.device(int(os.environ['CUDA_VISIBLE_DEVICES'])))
    worker=api.load_module('counterpoise_inner_existing_worker',ROOT/'scripts/pipeline/run_local_counterpoise.py')
    worker.run(case)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['run']);p.add_argument('--campaign',type=Path,required=True);a=p.parse_args()
    run(a.campaign.resolve())
