#!/usr/bin/env python
"""Retain the native physical pipeline with replayed tiled solvent response."""
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source = ROOT / 'scripts/pipeline/run_streamed_native.py'
    adapter = ROOT / 'scripts/pipeline/streamed_pcm_weighted.py'
    checker = ROOT / 'scripts/validation/qualify_streamed_weighted.py'
    module = api.load_module('preserved_streamed_native_physical_gates', source)
    module.ADAPTER, module.CHECKER = adapter, checker
    module.main()
    if sys.argv[1] == 'prepare':
        campaign = Path(sys.argv[sys.argv.index('--campaign') + 1]).resolve()
        plan = api.read_plan(campaign)
        plan.update(execution_driver=str(Path(__file__).resolve()),
            execution_profile='streamed_FP64_PCM_curvature_and_tiled_weighted_response')
        for path in (Path(__file__).resolve(), source, adapter, checker):
            plan['producer_sha256'][str(path)] = api.digest(path)
        atomic_write_json(campaign / 'plan.json', plan)


if __name__ == '__main__':
    main()
