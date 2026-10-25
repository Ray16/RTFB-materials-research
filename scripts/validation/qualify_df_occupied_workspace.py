#!/usr/bin/env python
"""Full four-atom spin/Hessian replays at production and stress occupied caps."""
import json
from pathlib import Path
import subprocess
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / 'scripts/pipeline/df_occupied_workspace.py'


def occupied_checks(folder, mult, cap):
    path = Path(folder) / 'compact/execution/occupied_workspace.json'
    row = json.loads(path.read_text())
    if (row['status'] != 'occupied_block_execution_finished_requires_replay_and_physical_checks'
            or row['occupied_cap'] != cap or row['precision'] != 'float64' or row['equations_changed']
            or len(row['checkpoints']) < 2):
        raise ValueError('Occupied-block execution and state-bound checkpoints must pass')
    if mult == 1 and (not row['batches'] or max(b['size'] for b in row['batches']) > cap
            or (cap == 4 and len(row['batches']) < 2)):
        raise ValueError('The restricted stress control must exercise multiple complete occupied blocks')
    api.verify_sources(row['source_sha256'])
    proof = {str(path): api.digest(path), **row['source_sha256']}
    for record in row['checkpoints']:
        if api.digest(record['path']) != record['sha256']:
            raise ValueError('Derivative state identity changed')
        state = json.loads(Path(record['path']).read_text())
        if not state['SCF_converged'] or not state['arrays']:
            raise ValueError('A complete converged derivative checkpoint is required')
        api.verify_sources(state['arrays'])
        proof[record['path']] = record['sha256']
    return proof


def checked_controls(root):
    root = Path(root); proof = {}
    base = api.load_module('occupied_cap_weighted_spin_checks', ROOT / 'scripts/validation/qualify_streamed_weighted.py')
    base.ADAPTER = ADAPTER
    solute = api.load_module('occupied_cap_solute_checks', ROOT / 'scripts/validation/qualify_solute_workspace.py')
    df = api.load_module('occupied_cap_DF_checks', ROOT / 'scripts/validation/qualify_df_derivative_workspace.py')
    for cap in (4, 32):
        folder = root / 'caps' / str(cap)
        proof.update(base.checked_controls(folder))
        for mult in (1, 2):
            control = folder / 'controls' / str(mult)
            row = json.loads((control / 'status.json').read_text())
            if not row.get('occupied_block_checks_passed'):
                raise ValueError('Both spin classes must pass at both occupied caps')
            proof.update(solute.offload_checks(control / 'compact/execution/solute_workspace.json', mult))
            proof.update(df.allocation_checks(control / 'compact/execution/DF_derivative_workspace.json', mult))
            proof.update(occupied_checks(control, mult, cap))
    return proof


def main():
    index = sys.argv.index('--occupied-cap')
    cap = int(sys.argv[index + 1]); del sys.argv[index:index + 2]
    if cap not in (4, 32):
        raise ValueError('Use the production or predeclared stress cap')
    source = ROOT / 'scripts/validation/qualify_smd_workspace.py'
    module = api.load_module('occupied_cap_full_spin_template', source)
    old = module.adapted_body
    def body(text):
        text = old(text)
        key = 'scripts/pipeline/smd_workspace_execution.py'
        call = 'budget_GiB=args.workspace_GiB)'
        if text.count(key) != 1 or text.count(call) != 1:
            raise ValueError('The full spin replay template changed')
        return text.replace(key, 'scripts/pipeline/df_occupied_workspace.py').replace(
            call, f'budget_GiB=args.workspace_GiB, occupied_cap={cap})')
    module.adapted_body = body
    observer = api.load_module('occupied_cap_resource_monitor', ROOT / 'scripts/pipeline/run_pi_stacking.py')
    with observer.MemoryMonitor() as memory:
        module.main()
    out = Path(sys.argv[sys.argv.index('--out') + 1])
    mult = int(sys.argv[sys.argv.index('--mult') + 1])
    weighted = api.load_module('occupied_cap_weighted_operators', ROOT / 'scripts/validation/qualify_streamed_weighted.py')
    solute = api.load_module('occupied_cap_solute_operators', ROOT / 'scripts/validation/qualify_solute_workspace.py')
    df = api.load_module('occupied_cap_DF_operators', ROOT / 'scripts/validation/qualify_df_derivative_workspace.py')
    evidence = weighted.operator_checks(out)
    evidence.update(solute.offload_checks(out / 'compact/execution/solute_workspace.json', mult))
    evidence.update(df.allocation_checks(out / 'compact/execution/DF_derivative_workspace.json', mult))
    evidence.update(occupied_checks(out, mult, cap))
    if not 0 < memory.peak <= 14336:
        raise ValueError('Four-atom control resource qualification exceeded its 14-GiB ceiling')
    path = out / 'status.json'; row = json.loads(path.read_text())
    row['source_sha256'].update(evidence)
    row['source_sha256'].update({str(p): api.digest(p) for p in (Path(__file__).resolve(), source, ADAPTER)})
    row.update(streamed_weighted_checks_passed=True, solute_workspace_checks_passed=True,
        DF_derivative_checks_passed=True, occupied_block_checks_passed=True, occupied_cap=cap,
        observed_peak_GPU_MiB=memory.peak, qualification_scope='Four-atom full replay only; full molecular workflow retains the 32-GiB device requirement',
        large_resource_qualified=False, native_G_production_qualified=False)
    path.write_text(json.dumps(row, indent=2) + '\n')


if __name__ == '__main__':
    main()
