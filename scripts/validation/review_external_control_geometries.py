#!/usr/bin/env python
"""Recheck saved external-control seed coordinates and source identities."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import rdMolAlign


ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(directory):
    directory = Path(directory)
    plan_path = directory / 'plan.json'
    plan = json.loads(plan_path.read_text())
    if (plan['status'] != 'neutral_external_control_seed_geometries_prepared_not_quantum_qualified'
            or plan['charged_state_minima_computed'] or plan['quantum_jobs_submitted']):
        raise ValueError('Unexpected seed status or scope')
    if sha(plan['challenge']) != plan['challenge_sha256']:
        raise ValueError('Challenge source changed')
    producer = ROOT / 'scripts/validation/prepare_external_control_geometries.py'
    if sha(producer) != plan['producer_sha256']:
        raise ValueError('Seed producer changed')
    graph = Chem.MolFromSmiles(plan['source_smiles'])
    if graph is None or Chem.MolToInchiKey(graph) != plan['inchikey']:
        raise ValueError('Graph identity changed')
    manifest = directory / plan['proposed_manifest']
    if sha(manifest) != plan['proposed_manifest_sha256']:
        raise ValueError('Proposed state manifest changed')
    with manifest.open(newline='') as stream:
        states = list(csv.DictReader(stream))
    if ([(s['state'], int(s['charge']), int(s['mult_hint'])) for s in states]
            != [('neu', 0, 1), ('red1', -1, 2), ('red2', -2, 1)]
            or any(s['smiles'] != plan['source_smiles'] or int(s['n_conf']) != plan['selected'] for s in states)):
        raise ValueError('Proposed manifest does not match control graph and states')
    path = ROOT / 'scripts/pipeline/polymer_conformer_thermals.py'
    spec = importlib.util.spec_from_file_location('external_seed_graph_alignment', path)
    geometry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
    mols, distances, files = [], [], []
    for entry in plan['entries']:
        xyz = directory / entry['path']
        if sha(xyz) != entry['sha256']:
            raise ValueError('Seed geometry hash changed')
        atoms = read(xyz)
        if len(atoms) != Chem.AddHs(graph).GetNumAtoms() or not np.isfinite(atoms.positions).all():
            raise ValueError('Composition or coordinates invalid')
        pair = np.linalg.norm(atoms.positions[:, None, :] - atoms.positions[None, :, :], axis=2)
        np.fill_diagonal(pair, np.inf)
        minimum = float(pair.min())
        if minimum < 0.6:
            raise ValueError('Unphysical close contact in seed geometry')
        distances.append(minimum)
        mols.append(geometry.alignment_molecule(plan['source_smiles'], xyz))
        files.append(dict(path=entry['path'], sha256=entry['sha256'], min_atom_distance_A=minimum))
    if len(mols) != plan['selected'] or not mols:
        raise ValueError('Selected geometry count differs')
    rmsd = [rdMolAlign.GetBestRMS(Chem.Mol(mols[i]), mols[j], maxMatches=10000)
            for i in range(len(mols)) for j in range(i)]
    if rmsd and (not all(math.isfinite(v) for v in rmsd) or min(rmsd) + 1e-6 < plan['post_optimization_heavy_rmsd_A']):
        raise ValueError('Post-optimization geometry separation does not match policy')
    return dict(status='external_neutral_seed_geometry_review_passed_not_quantum',
                plan_sha256=sha(plan_path), producer_sha256=sha(producer), reviewer_sha256=sha(__file__),
                alignment_sha256=sha(path), selected=len(files), files=files,
                minimum_atom_distance_A=min(distances), minimum_pairwise_heavy_rmsd_A=min(rmsd) if rmsd else None,
                maximum_pairwise_heavy_rmsd_A=max(rmsd) if rmsd else None,
                checked_graph_and_stereochemistry=True, quantum_optimization_run=False,
                charged_state_geometry_qualified=False, experimental_accuracy_inferred=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = review(args.input.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'status': result['status'], 'selected': result['selected']}))


if __name__ == '__main__':
    main()
