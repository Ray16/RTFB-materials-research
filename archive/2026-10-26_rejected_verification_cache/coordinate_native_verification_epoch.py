#!/usr/bin/env python
"""CPU-native coordinator with per-pass large-checkpoint verification reuse.

Remote launchers and scientific workers still use the original strict verifier.
The protocol, input plans, numerical code and resource bounds are unchanged.
"""
import inspect
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    upstream=ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py'
    cache_path=ROOT/'scripts/fleet/verification_epoch.py'
    cache_module=api.load_module('native_CPU_verification_epoch',cache_path)
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),upstream,cache_path)}
    original=api.verify_sources
    if '--batch' not in sys.argv:raise ValueError('Declare an exact existing native batch')
    batch=Path(sys.argv[sys.argv.index('--batch')+1]).resolve()
    epoch=cache_module.VerificationEpoch(original,report=batch/'CPU_verification_epochs.jsonl')
    def verify(mapping):original(proof);epoch.verify(mapping)
    api.verify_sources=verify
    module=api.load_module('unchanged_native_dispatch_with_epoch_verification',upstream)
    body=inspect.getsource(module.coordinate)
    anchor="            while time.time() < record['deadline']:"
    if body.count(anchor)!=1:raise ValueError('Native dispatch-pass template changed')
    body=body.replace(anchor,anchor+'\n                _verification_epoch.begin()',1)
    namespace={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(module),_verification_epoch=epoch),namespace)
    module.coordinate=namespace['coordinate']
    # SCRIPT deliberately remains the original coordinator: remote launchers
    # and fresh scientific interpreters receive no verification cache.
    try:module.main()
    finally:epoch.begin()


if __name__=='__main__':main()
