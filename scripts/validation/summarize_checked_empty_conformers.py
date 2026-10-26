#!/usr/bin/env python
"""Summarize a declared set of already reviewed empty-chain native states.

This is a within-set range diagnostic, not a conformer ensemble or cycle.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


EV_KJ_MOL = 96.48533212331002  # elementary charge times exact Avogadro constant


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(manifest, root):
    rows = []
    signature = None
    seen = set()
    for item in manifest['reviews']:
        path = root / item['path']
        if sha(path) != item['sha256']:
            raise ValueError('Review file differs from declared hash')
        record = json.loads(path.read_text())
        task = record['task']
        model = record['experimental_model']
        if (record['status'] != 'alternate_native_state_independently_reconstructed'
                or record['primary_candidate_admission'] is not False
                or task['occupied_sites'] != [] or model['solvent'] != 'acetonitrile'):
            raise ValueError('Not a checked empty-chain alternative-solvent-model state')
        current = (task['graph_sha256'], tuple(tuple(x) for x in task['backbone_stereo_fingerprint']),
                   task['configuration_label'], record['experimental_method_id'])
        if signature is not None and current != signature:
            raise ValueError('Graph, tacticity or method differs across records')
        signature = current
        rank = task['search_rank']
        if type(rank) is not int or rank < 1 or rank in seen:
            raise ValueError('Each input must have a distinct positive search rank')
        seen.add(rank)
        energies = {str(c): float(record['G_by_cutoff_eV'][str(c)]) for c in (25, 50, 100)}
        if not all(math.isfinite(g) for g in energies.values()):
            raise ValueError('Non-finite free energy')
        if abs(float(record['G_eV']) - energies['50']) > 1e-9:
            raise ValueError('Primary free energy does not match the 50 cm-1 value')
        rows.append(dict(key=task['key'], search_rank=rank, G_by_cutoff_eV=energies,
                         minimum_frequency_cm=record['minimum_frequency_cm'], review=item))
    if not rows:
        raise ValueError('At least one checked record required')
    rows.sort(key=lambda r: r['search_rank'])
    results = {}
    for cutoff in ('25', '50', '100'):
        energies = [r['G_by_cutoff_eV'][cutoff] for r in rows]
        low = min(energies)
        results[cutoff] = dict(range_kJ_mol=(max(energies) - low) * EV_KJ_MOL,
                               relative_to_lowest_kJ_mol=[(e - low) * EV_KJ_MOL for e in energies],
                               running_minimum_eV=[min(energies[:n]) for n in range(1, len(rows)+1)])
    return dict(status='checked_empty_chain_free_energy_range_not_ensemble',
                configuration=signature[2], graph_sha256=signature[0], method_id=signature[3],
                solvent='acetonitrile', records=rows, by_rotor_cutoff_cm=results,
                conversion_kJ_mol_per_eV=EV_KJ_MOL, comparison='within same empty graph and alternate model',
                finite_ensemble_calculated=False, graph_duplicate_review_performed=False,
                wider_sampling_converged=False, polymer_cycle_or_binding_inferred=False,
                primary_model_cache_modified=False, quantum_jobs_submitted=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    result = summarize(json.loads(a.manifest.read_text()), a.root)
    result.update(manifest_sha256=sha(a.manifest), producer_sha256=sha(__file__))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
