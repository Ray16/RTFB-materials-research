#!/usr/bin/env python
"""Native thermochemistry with replayed occupied-orbital DF blocking."""
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source = ROOT / 'scripts/pipeline/run_streamed_native.py'
    adapter = ROOT / 'scripts/pipeline/df_occupied_workspace.py'
    checker = ROOT / 'scripts/validation/qualify_df_occupied_workspace.py'
    module = api.load_module('native_gates_occupied_block_workspace', source)
    module.ADAPTER, module.CHECKER = adapter, checker
    module.main()
    if sys.argv[1] == 'prepare':
        campaign = Path(sys.argv[sys.argv.index('--campaign') + 1]).resolve()
        plan = api.read_plan(campaign)
        plan.update(execution_driver=str(Path(__file__).resolve()),
            execution_profile='streamed_FP64_PCM_and_32_occupied_DF_blocks')
        for path in (Path(__file__).resolve(), source, adapter, checker):
            plan['producer_sha256'][str(path)] = api.digest(path)
        atomic_write_json(campaign / 'plan.json', plan)


if __name__ == '__main__':
    main()
