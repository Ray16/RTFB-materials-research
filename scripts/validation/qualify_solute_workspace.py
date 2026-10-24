#!/usr/bin/env python
"""Full spin/Hessian controls for temporary pure-solute matrix suspension."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / 'scripts/pipeline/pcm_solute_workspace.py'
PARENT_CHECKER = ROOT / 'scripts/validation/qualify_streamed_weighted.py'


def offload_checks(path, mult):
    row = json.loads(path.read_text())
    label = 'rks' if mult == 1 else 'uks'
    if (row['status'] != 'solute_workspace_finished_requires_full_controls_and_physical_review'
            or min(row['calls'][label + '.' + key] for key in ('partial_hess_elec', 'make_h1')) < 1
            or row['suspensions'] < 2 or row['restorations'] != row['suspensions']
            or row['exact_matrix_round_trips'] != 3 * row['suspensions']
            or row['maximum_suspended_bytes'] <= 0 or row['equations_changed'] or row['precision_changed']):
        raise ValueError('Pure solute calls must restore every FP64 solvent matrix exactly')
    api.verify_sources(row['source_sha256'])
    return {**row['source_sha256'], str(path): api.digest(path)}


def checked_controls(root):
    root = Path(root)
    parent = api.load_module('solute_workspace_weighted_control_checks', PARENT_CHECKER)
    parent.ADAPTER = ADAPTER
    proof = parent.checked_controls(root)
    for mult in (1, 2):
        control = root / 'controls' / str(mult)
        row = json.loads((control / 'status.json').read_text())
        if not row.get('solute_workspace_checks_passed'):
            raise ValueError('Both full solute-workspace spin controls must pass')
        proof.update(offload_checks(control / 'compact/execution/solute_workspace.json', mult))
    return proof


def main():
    source = ROOT / 'scripts/validation/qualify_smd_workspace.py'
    module = api.load_module('solute_workspace_full_spin_template', source)
    old = module.adapted_body
    def body(text):
        text = old(text)
        key = 'scripts/pipeline/smd_workspace_execution.py'
        if text.count(key) != 1:
            raise ValueError('The original full spin replay template changed')
        return text.replace(key, 'scripts/pipeline/pcm_solute_workspace.py')
    module.adapted_body = body
    module.main()
    out = Path(sys.argv[sys.argv.index('--out') + 1])
    mult = int(sys.argv[sys.argv.index('--mult') + 1])
    parent = api.load_module('solute_workspace_weighted_operators', PARENT_CHECKER)
    evidence = parent.operator_checks(out)
    evidence.update(offload_checks(out / 'compact/execution/solute_workspace.json', mult))
    path = out / 'status.json'
    row = json.loads(path.read_text())
    row['source_sha256'].update(evidence)
    row['source_sha256'].update({str(p): api.digest(p) for p in (Path(__file__).resolve(), source, ADAPTER, PARENT_CHECKER)})
    row.update(streamed_weighted_checks_passed=True, solute_workspace_checks_passed=True,
        large_resource_qualified=False, native_G_production_qualified=False)
    path.write_text(json.dumps(row, indent=2) + '\n')


if __name__ == '__main__':
    main()
