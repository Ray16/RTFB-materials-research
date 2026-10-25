#!/usr/bin/env python
"""Full Hessian replay of streamed AO derivatives and bounded DF batches."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / 'scripts/pipeline/df_derivative_workspace.py'
SOLUTE_CHECKER = ROOT / 'scripts/validation/qualify_solute_workspace.py'


def allocation_checks(path, mult):
    row = json.loads(path.read_text())
    spin = 'rhf' if mult == 1 else 'uhf'
    allocations = [r for r in row['allocations'] if r['spin'] == spin]
    if (row['status'] != 'DF_derivative_workspace_finished_requires_full_replay_and_physical_checks'
            or row['calls'][spin] < 1 or not allocations
            or any(r['full_atom_AO_cache'] or r['blksize'] <= 0 or r['auxiliary_batch'] <= 0 for r in allocations)
            or {c['stage'] for c in row['checkpoints']} !=
                {('rks' if mult == 1 else 'uks') + '_' + k for k in ('partial_hess_elec', 'make_h1')}
            or row['equations_changed'] or row['precision_changed']):
        raise ValueError('The streamed DF path and both subroutine checkpoints must be exercised')
    api.verify_sources(row['source_sha256'])
    evidence = {str(path): api.digest(path), **row['source_sha256']}
    for record in row['checkpoints']:
        if api.digest(record['path']) != record['sha256']:
            raise ValueError('A DF derivative checkpoint changed')
        evidence[record['path']] = record['sha256']
    return evidence


def checked_controls(root):
    root = Path(root)
    parent = api.load_module('DF_underlying_weighted_replay', ROOT / 'scripts/validation/qualify_streamed_weighted.py')
    parent.ADAPTER = ADAPTER
    proof = parent.checked_controls(root)
    solute = api.load_module('DF_solute_round_trip_checks', SOLUTE_CHECKER)
    for mult in (1, 2):
        folder = root / 'controls' / str(mult)
        if not json.loads((folder / 'status.json').read_text()).get('DF_derivative_checks_passed'):
            raise ValueError('Both full DF derivative spin controls must pass')
        proof.update(solute.offload_checks(folder / 'compact/execution/solute_workspace.json', mult))
        proof.update(allocation_checks(folder / 'compact/execution/DF_derivative_workspace.json', mult))
    return proof


def main():
    source = ROOT / 'scripts/validation/qualify_smd_workspace.py'
    module = api.load_module('DF_full_spin_replay_template', source)
    old = module.adapted_body
    def body(text):
        text = old(text)
        key = 'scripts/pipeline/smd_workspace_execution.py'
        if text.count(key) != 1:
            raise ValueError('The full spin replay template changed')
        return text.replace(key, 'scripts/pipeline/df_derivative_workspace.py')
    module.adapted_body = body
    module.main()
    out = Path(sys.argv[sys.argv.index('--out') + 1])
    mult = int(sys.argv[sys.argv.index('--mult') + 1])
    weighted = api.load_module('DF_weighted_operator_controls', ROOT / 'scripts/validation/qualify_streamed_weighted.py')
    solute = api.load_module('DF_solute_matrix_controls', SOLUTE_CHECKER)
    evidence = weighted.operator_checks(out)
    evidence.update(solute.offload_checks(out / 'compact/execution/solute_workspace.json', mult))
    evidence.update(allocation_checks(out / 'compact/execution/DF_derivative_workspace.json', mult))
    path = out / 'status.json'
    row = json.loads(path.read_text())
    row['source_sha256'].update(evidence)
    row['source_sha256'].update({str(p): api.digest(p) for p in (Path(__file__).resolve(), source, ADAPTER, SOLUTE_CHECKER)})
    row.update(streamed_weighted_checks_passed=True, solute_workspace_checks_passed=True,
        DF_derivative_checks_passed=True, large_resource_qualified=False, native_G_production_qualified=False)
    path.write_text(json.dumps(row, indent=2) + '\n')


if __name__ == '__main__':
    main()
