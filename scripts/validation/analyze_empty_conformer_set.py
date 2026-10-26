#!/usr/bin/env python
"""Analyze a source-bound finite set of reviewed alternate-model empty states."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from redox.properties.conformer_ensemble import cluster_minima, summarize


ROOT = Path(__file__).resolve().parents[2]


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def analyze(manifest_path, root):
    arithmetic = module(root / 'scripts/validation/summarize_checked_empty_conformers.py', 'empty_range')
    geometry = module(root / 'scripts/pipeline/polymer_conformer_thermals.py', 'empty_geometry')
    manifest = json.loads(manifest_path.read_text())
    checked = arithmetic.summarize(manifest, root)
    structures, rows, inputs = [], [], {}
    for item in checked['records']:
        review_path = root / item['review']['path']
        report = json.loads(review_path.read_text())
        campaign = Path(report['campaign'])
        plan_path = campaign / 'plan.json'
        status_path = campaign / 'experimental_native/status.json'
        for p in (plan_path, status_path):
            if arithmetic.sha(p) != report['source_sha256'][str(p)]:
                raise ValueError('Native plan or status changed since independent review')
        plan = json.loads(plan_path.read_text())
        status = json.loads(status_path.read_text())
        xyz = Path(status['raw_root']) / 'geometry.xyz'
        parent_path = Path(plan['electronic_review'])
        for p in (xyz, parent_path):
            if arithmetic.sha(p) != report['source_sha256'][str(p)]:
                raise ValueError('Geometry or electronic review changed since independent review')
        if arithmetic.sha(xyz) != plan['scored_geometry_sha256']:
            raise ValueError('Native and electronic geometry differ')
        electronic = json.loads(parent_path.read_text())
        energy = float(electronic['E_smd_eV'])
        if not math.isfinite(energy):
            raise ValueError('Non-finite electronic energy')
        structures.append(geometry.alignment_molecule(report['task']['smiles'], xyz))
        rows.append(dict(key=item['key'], E_eV=energy, G_by_cutoff_eV=item['G_by_cutoff_eV']))
        inputs[item['key']] = dict(geometry=str(xyz), geometry_sha256=arithmetic.sha(xyz),
                                   electronic_review=str(parent_path), electronic_review_sha256=arithmetic.sha(parent_path),
                                   native_review=str(review_path), native_review_sha256=arithmetic.sha(review_path))
    n = len(rows)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            matrix[i, j] = matrix[j, i] = rdMolAlign.GetBestRMS(
                Chem.Mol(structures[i]), structures[j], maxMatches=10000)
    by_cutoff = {}
    for rmsd_cutoff in (0.05, 0.1, 0.2):
        groups = cluster_minima(rows, matrix, rmsd_cutoff_A=rmsd_cutoff, energy_tolerance_kJ_mol=0.2)
        by_rotor = {}
        for rotor in ('25', '50', '100'):
            rotor_rows = [dict(key=r['key'], E_eV=r['E_eV'], G_eV=r['G_by_cutoff_eV'][rotor]) for r in rows]
            by_rotor[rotor] = summarize(rotor_rows, groups, [r['key'] for r in rows], temperature_K=298.15)
        by_cutoff[str(rmsd_cutoff)] = by_rotor
    return dict(status='finite_checked_empty_conformer_set_not_converged', manifest_sha256=arithmetic.sha(manifest_path),
                source_sha256={str(p): arithmetic.sha(p) for p in (Path(__file__).resolve(),
                    root / 'scripts/validation/summarize_checked_empty_conformers.py',
                    root / 'scripts/pipeline/polymer_conformer_thermals.py',
                    root / 'src/redox/properties/conformer_ensemble.py')},
                inputs=inputs, keys=[r['key'] for r in rows], rmsd_A=matrix.tolist(),
                geometry_check='declared graph and stereochemistry checked before symmetry-aware heavy-atom alignment',
                maximum_alignment_matches=10000, energy_duplicate_tolerance_kJ_mol=0.2,
                partitions_by_rmsd_cutoff_A=by_cutoff,
                scope='Finite empty-chain set in one configuration and alternate acetonitrile model; no experimental loading or grafting-cycle inference',
                converged_ensemble=False, physical_QC_repeated=False, quantum_jobs_submitted=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = analyze(args.manifest.resolve(), args.root.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
