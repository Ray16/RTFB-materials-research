#!/usr/bin/env python
"""Use checked residual solves within the existing full native admission path."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source=ROOT/'scripts/pipeline/run_streamed_native.py'
    adapter=ROOT/'scripts/pipeline/residual_response_solver.py'
    checker=ROOT/'scripts/validation/review_residual_response.py'
    module=api.load_module('unchanged_native_physical_QC_residual_solver',source)
    module.ADAPTER,module.CHECKER=adapter,checker
    module.main()
    if sys.argv[1]=='prepare':
        campaign=Path(sys.argv[sys.argv.index('--campaign')+1]).resolve()
        plan=api.read_plan(campaign)
        comparison=Path(plan['streamed_controls'])/'comparison.json'
        runtime=Path(json.loads(comparison.read_text())['runtime_profile'])
        plan.update(execution_driver=str(Path(__file__).resolve()),runtime_profile=str(runtime),
            execution_profile='residual_checked_FP64_CPHF_with_streamed_PCM_and_occupied_DF')
        plan['input_sha256'][str(runtime)]=api.digest(runtime)
        for p in (Path(__file__).resolve(),source,adapter,checker,ROOT/'scripts/fleet/tensor_runtime_worker.py'):
            plan['producer_sha256'][str(p)]=api.digest(p)
        plan['workflow_source_sha256']={str(p):api.digest(p) for p in (
            ROOT/'scripts/pipeline/dispatch_native_case.py',ROOT/'scripts/pipeline/publish_candidate_native_component.py')}
        atomic_write_json(campaign/'plan.json',plan)


if __name__=='__main__':
    main()
