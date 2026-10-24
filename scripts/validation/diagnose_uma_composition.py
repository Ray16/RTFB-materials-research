#!/usr/bin/env python
"""Ablate the historical composition-mean rule without changing installed code.

This deliberately changes in-memory model behavior for diagnosis only. Its
outputs cannot qualify the official current model or a production migration.
"""
import argparse
import json
from pathlib import Path
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    verify_gate_owner()
    from fairchem import core

    original = core.pretrained_mlip.get_predict_unit
    changes = []

    def diagnostic_predictor(*a, **kw):
        predictor = original(*a, **kw)
        backbone = predictor.model.module.backbone
        old_id = backbone.model_id
        if old_id != 'UMA-S-1.2.1' or predictor.inference_settings.nitro:
            raise RuntimeError('The diagnostic requires the official 1.2.1 general model')
        # escn_moe.py uses this field only to select include_self in the
        # composition mean. UMA-S-1.2 includes the zero-initialized accumulator;
        # the official 1.2.1 path does not. No weights or references are changed.
        backbone.model_id = 'UMA-S-1.2'
        changes.append(dict(original_model_id=old_id, diagnostic_model_id=backbone.model_id,
            model_version=backbone.model_version, include_self=True))
        return predictor

    core.pretrained_mlip.get_predict_unit = diagnostic_predictor
    worker = api.load_module('V100_general_FP32_composition_diagnostic',
        ROOT/'scripts/validation/qualify_uma_v100.py')
    previous = sys.argv
    try:
        sys.argv = [str(worker.__file__), '--campaign', str(args.campaign), '--out', str(args.out)]
        worker.main()
        path = args.out/'status.json'
        result = json.loads(path.read_text())
        result.update(diagnostic_only=True, production_qualified=False,
            hypothesis='composition averaging include_self explains legacy/current drift',
            in_memory_changes=changes, diagnostic_source_sha256=api.digest(Path(__file__)),
            completed=time.time())
        atomic_write_json(path, result)
    finally:
        core.pretrained_mlip.get_predict_unit = original
        sys.argv = previous


if __name__ == '__main__':
    main()
