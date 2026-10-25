#!/usr/bin/env python
"""Source-bound cuTENSOR execution with the unchanged native physical gates."""
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source = ROOT / 'scripts/pipeline/run_streamed_native.py'
    adapter = ROOT / 'scripts/pipeline/df_occupied_workspace.py'
    checker = ROOT / 'scripts/validation/review_cutensor_controls.py'
    module = api.load_module('native_physical_gates_cutensor_runtime', source)
    module.ADAPTER, module.CHECKER = adapter, checker
    module.main()
    if sys.argv[1] == 'prepare':
        campaign = Path(sys.argv[sys.argv.index('--campaign') + 1]).resolve()
        plan = api.read_plan(campaign)
        comparison = Path(plan['streamed_controls']) / 'comparison.json'
        runtime = Path(__import__('json').loads(comparison.read_text())['runtime_profile'])
        plan.update(execution_driver=str(Path(__file__).resolve()),
            execution_profile='cutensor_FP64_streamed_PCM_and_32_occupied_DF_blocks', runtime_profile=str(runtime))
        plan['input_sha256'][str(runtime)] = api.digest(runtime)
        for path in (Path(__file__).resolve(), source, adapter, checker,
                ROOT / 'scripts/fleet/tensor_runtime_worker.py'):
            plan['producer_sha256'][str(path)] = api.digest(path)
        plan['workflow_source_sha256'] = {str(p): api.digest(p) for p in (
            ROOT / 'scripts/pipeline/dispatch_native_case.py', ROOT / 'scripts/pipeline/publish_candidate_native_component.py')}
        atomic_write_json(campaign / 'plan.json', plan)


if __name__ == '__main__':
    main()
