#!/usr/bin/env python
"""Full singlet/doublet Hessian replay of a residual-checked response solver."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source = ROOT / 'scripts/validation/qualify_smd_workspace.py'
    adapter = ROOT / 'scripts/pipeline/residual_response_solver.py'
    module = api.load_module('full_response_solver_replay_template', source)
    previous = module.adapted_body
    def body(text):
        text = previous(text)
        key = 'scripts/pipeline/smd_workspace_execution.py'
        if text.count(key) != 1:
            raise ValueError('Full spin replay template changed')
        return text.replace(key, 'scripts/pipeline/residual_response_solver.py')
    module.adapted_body = body
    observer = api.load_module('response_control_resource_observer', ROOT / 'scripts/pipeline/run_pi_stacking.py')
    with observer.MemoryMonitor() as memory:
        module.main()
    if not 0 < memory.peak <= 14336:
        raise ValueError('Small full response control exceeded its 14 GiB device ceiling')
    out = Path(sys.argv[sys.argv.index('--out') + 1])
    proof_path = out / 'compact/execution/response_solver.json'
    proof = json.loads(proof_path.read_text())
    if (proof['status'] != 'response_equations_solved_full_physical_QC_separate'
            or not proof['calls'] or any(call['true_residual_norm'] > call['solver_absolute_target']
                or call['solver_info'] != 0 for call in proof['calls'])):
        raise ValueError('Original-equation residuals must pass for every response block')
    row = json.loads((out / 'status.json').read_text())
    row['source_sha256'].update({str(p): api.digest(p) for p in (Path(__file__).resolve(), source, adapter, proof_path)})
    row.update(response_solver_checks_passed=True, full_native_G_qualified=False,
        observed_peak_GPU_MiB=memory.peak)
    atomic_write_json(out / 'status.json', row)


if __name__ == '__main__':
    main()
