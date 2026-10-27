#!/usr/bin/env python3
"""Apply the declared 10-to-20 diagnostic to source-checked empty backbones."""
import argparse
import json
from pathlib import Path

from redox.properties.conformer_ensemble import batch_change
from analyze_empty_conformer_set import analyze, module

ROOT = Path(__file__).resolve().parents[2]


def checked_batch(root, manifest, saved):
    recomputed = analyze(manifest, root)
    if json.loads(saved.read_text()) != recomputed:
        raise ValueError('Saved finite partition differs from reconstructed inputs')
    checker = module(root/'scripts/validation/summarize_checked_empty_conformers.py', 'batch_inputs')
    reviewed = checker.summarize(json.loads(manifest.read_text()), root)
    return recomputed, reviewed


def compare(first, expanded, first_review, expanded_review):
    signature = ('configuration', 'graph_sha256', 'method_id', 'solvent')
    if any(first_review[k] != expanded_review[k] for k in signature):
        raise ValueError('Graph, configuration, solvent or method differs')
    for expected, finite, reviewed in ((10, first, first_review), (20, expanded, expanded_review)):
        if ([r['search_rank'] for r in reviewed['records']] != list(range(1, expected+1))
                or finite['keys'] != [r['key'] for r in reviewed['records']]):
            raise ValueError('Requires exactly the first 10 and first 20 ranked inputs')
    comparisons = {}
    for rmsd in ('0.05', '0.1', '0.2'):
        comparisons[rmsd] = {
            rotor: batch_change(first['partitions_by_rmsd_cutoff_A'][rmsd][rotor],
                                expanded['partitions_by_rmsd_cutoff_A'][rmsd][rotor], 2.0)
            for rotor in ('25', '50', '100')}
    return dict(status='checked_empty_backbone_10_to_20_batch_comparison',
                **{k: first_review[k] for k in signature},
                by_rmsd_A_and_rotor_cm=comparisons,
                global_ensemble_convergence_proven=False, polymer_cycle_computed=False,
                primary_model_results_modified=False, quantum_jobs_submitted=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    for name in ('first-manifest', 'first-finite', 'expanded-manifest', 'expanded-finite', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    root = a.root.resolve()
    first, fr = checked_batch(root, a.first_manifest.resolve(), a.first_finite.resolve())
    expanded, er = checked_batch(root, a.expanded_manifest.resolve(), a.expanded_finite.resolve())
    result = compare(first, expanded, fr, er)
    checker = module(root/'scripts/validation/summarize_checked_empty_conformers.py', 'batch_hashes')
    sources = [a.first_manifest, a.first_finite, a.expanded_manifest, a.expanded_finite,
               Path(__file__), root/'docs/CONFORMER_ENSEMBLES.md']
    result['source_sha256'] = {str(path.resolve()): checker.sha(path) for path in sources}
    result['scope'] = ('Within-protocol stability of one empty-backbone configuration in the alternate '
                       'acetonitrile model; not a complete polymer cycle or global sampling proof.')
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
