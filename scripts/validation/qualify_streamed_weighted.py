#!/usr/bin/env python
"""Full original/streamed-weighted energy, geometry, gradient and Hessian replay."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / 'scripts/pipeline/streamed_pcm_weighted.py'


def operator_checks(out):
    folder = Path(out) / 'compact/execution'
    weighted = json.loads((folder / 'weighted_response.json').read_text())
    streamed = json.loads((folder / 'streamed_PCM.json').read_text())
    contracted = json.loads((folder / 'contracted_R.json').read_text())
    if (weighted['status'] != 'weighted_response_finished_requires_controls_and_physical_review'
            or min(weighted['transformed_calls'].values()) < 1
            or min(weighted[k] for k in ('weighted_matrices', 'matrix_tiles', 'bounded_response_contracts')) < 1
            or weighted['maximum_matrix_tile_bytes'] > 512 * 1024**2
            or weighted['precision'] != 'float64' or weighted['equations_changed']
            or streamed['status'] != 'streamed_PCM_finished_requires_full_replay_and_physical_QC'
            or min(streamed[k] for k in ('switch_contractions', 'switch_primitive_calls',
                'diagonal_contractions', 'offdiagonal_primitive_calls')) < 1
            or streamed['maximum_offdiagonal_primitive_bytes'] > 512 * 1024**2
            or streamed['peak_GPU_block_cache_bytes'] > 128 * 1024**2
            or contracted['R_applications'] < 1):
        raise ValueError('Every weighted response and streamed-curvature operator must be exercised')
    evidence = {}
    for name, row in (('weighted_response', weighted), ('streamed_PCM', streamed), ('contracted_R', contracted)):
        api.verify_sources(row['source_sha256'])
        evidence.update(row['source_sha256'])
        path = folder / (name + '.json')
        evidence[str(path)] = api.digest(path)
    return evidence


def checked_controls(root):
    base = api.load_module('weighted_full_first_derivative_checks',
        ROOT / 'scripts/validation/qualify_pcm_first_response.py')
    base.ADAPTER = ADAPTER
    root = Path(root)
    proof = base.checked_controls(root)
    for mult in (1, 2):
        path = root / 'controls' / str(mult) / 'status.json'
        row = json.loads(path.read_text())
        if not row.get('streamed_weighted_checks_passed') or row['mult'] != mult:
            raise ValueError('Both full spin/Hessian controls must pass')
        proof.update(operator_checks(path.parent))
    return proof


def main():
    source = ROOT / 'scripts/validation/qualify_smd_workspace.py'
    module = api.load_module('weighted_full_spin_template', source)
    old = module.adapted_body
    def body(text):
        text = old(text)
        key = 'scripts/pipeline/smd_workspace_execution.py'
        if text.count(key) != 1:
            raise ValueError('The original full spin replay template changed')
        return text.replace(key, 'scripts/pipeline/streamed_pcm_weighted.py')
    module.adapted_body = body
    module.main()
    out = Path(sys.argv[sys.argv.index('--out') + 1])
    evidence = operator_checks(out)
    path = out / 'status.json'
    row = json.loads(path.read_text())
    row['source_sha256'].update(evidence)
    row['source_sha256'].update({str(p): api.digest(p) for p in (Path(__file__).resolve(), source, ADAPTER)})
    row.update(streamed_weighted_checks_passed=True, large_resource_qualified=False,
        native_G_production_qualified=False)
    path.write_text(json.dumps(row, indent=2) + '\n')


if __name__ == '__main__':
    main()
