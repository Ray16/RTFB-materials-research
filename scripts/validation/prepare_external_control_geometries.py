#!/usr/bin/env python
"""Generate reproducible neutral-graph starting geometries for an external control.

These force-field structures are inputs only. They are not UMA or DFT results,
charged-state minima, native thermochemistry, or benchmark admission.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign

from redox.build.geometry import geometry_stereo_status
from redox.validation.protonation import hydrogen_topology_status


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(challenge_path, out, seed=20261026, attempts=32, top_k=10, rmsd_cutoff=0.15):
    challenge = json.loads(challenge_path.read_text())
    if challenge['status'] != 'source_checked_NDI_external_challenge_not_scored' or challenge['computed_potentials'] is not None:
        raise ValueError('Expected an unscored source-checked external control')
    if not (0 < seed < 2**31 and 1 <= top_k <= attempts <= 200 and 0 < rmsd_cutoff <= 1):
        raise ValueError('Invalid deterministic sampling parameters')
    state_map = {'neutral': 'neu', 'reduced_1': 'red1', 'reduced_2': 'red2'}
    expected = {'neutral': (0, 1), 'reduced_1': (-1, 2), 'reduced_2': (-2, 1)}
    states = challenge['proposed_states']
    if (len(states) != len(expected) or {row['label'] for row in states} != set(expected)
            or any((row['charge'], row['multiplicity']) != expected[row['label']] for row in states)):
        raise ValueError('Challenge states differ from the supported three-state input schema')
    molecule = Chem.MolFromSmiles(challenge['molecule']['smiles'])
    if molecule is None or Chem.MolToInchiKey(molecule) != challenge['molecule']['inchikey']:
        raise ValueError('Challenge molecular identity differs')
    if Chem.GetFormalCharge(molecule) != 0:
        raise ValueError('This preparer accepts neutral source graphs only')
    molecule = Chem.AddHs(molecule)
    if AllChem.MMFFGetMoleculeProperties(molecule) is None:
        raise ValueError('MMFF parameters unavailable; do not silently change force fields')
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    params.pruneRmsThresh = 0.1
    ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=attempts, params=params))
    if not ids:
        raise ValueError('No conformers embedded')
    optimized = AllChem.MMFFOptimizeMoleculeConfs(molecule, numThreads=1, maxIters=2000)
    ranked = sorted((float(energy), int(cid)) for cid, (status, energy) in zip(ids, optimized)
                    if status == 0 and math.isfinite(energy))
    if not ranked:
        raise ValueError('No MMFF optimization converged')
    heavy = Chem.RemoveHs(Chem.Mol(molecule))
    retained = []
    for energy, cid in ranked:
        if any(rdMolAlign.GetBestRMS(Chem.Mol(heavy), heavy, prbId=cid, refId=prior,
                                     maxMatches=10000) < rmsd_cutoff for _, prior in retained):
            continue
        retained.append((energy, cid))
        if len(retained) == top_k:
            break
    if not retained:
        raise ValueError('No geometries retained')
    out.mkdir(parents=True, exist_ok=False)
    (out / 'conformers').mkdir()
    entries = []
    with tempfile.TemporaryDirectory(dir=out) as temporary:
        for rank, (energy, cid) in enumerate(retained, 1):
            xyz = Chem.MolToXYZBlock(molecule, confId=cid, precision=8)
            scratch = Path(temporary) / f'{rank:03d}.xyz'
            scratch.write_text(xyz)
            topology = hydrogen_topology_status(challenge['molecule']['smiles'], scratch)
            stereo = geometry_stereo_status(challenge['molecule']['smiles'], scratch)
            if topology != 'ok' or stereo not in ('ok', 'not_applicable'):
                raise ValueError('Generated geometry failed declared graph or stereochemistry')
            path = out / 'conformers' / f'conf_{rank:03d}.xyz'
            path.write_text(xyz)
            entries.append(dict(rank=rank, conformer_id=cid, path=str(path.relative_to(out)), sha256=sha(path),
                                mmff_energy_kcal_mol=energy, topology=topology, stereochemistry=stereo))
    manifest_path = out / 'manifest.csv'
    with manifest_path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['id', 'name', 'family', 'state', 'charge',
            'mult_hint', 'n_e', 'n_conf', 'solv_preopt', 'smiles'])
        writer.writeheader()
        for state in challenge['proposed_states']:
            writer.writerow(dict(id='bm_ndi_styrene_external', name=challenge['molecule']['name'],
                family='NDI experimental method control', state=state_map[state['label']],
                charge=state['charge'], mult_hint=state['multiplicity'], n_e=state['charge'],
                n_conf=len(entries), solv_preopt=0, smiles=challenge['molecule']['smiles']))
    report = dict(status='neutral_external_control_seed_geometries_prepared_not_quantum_qualified',
        challenge=str(challenge_path.resolve()), challenge_sha256=sha(challenge_path),
        producer_sha256=sha(__file__), source_smiles=challenge['molecule']['smiles'],
        inchikey=challenge['molecule']['inchikey'], rdkit_version=rdBase.rdkitVersion,
        method='ETKDGv3 followed by MMFF94 optimization; neutral molecular graph only',
        seed=seed, requested_embeddings=attempts, embedded=len(ids), mmff_converged=len(ranked),
        optimization_iterations=2000, embedding_prune_A=0.1, post_optimization_heavy_rmsd_A=rmsd_cutoff,
        proposed_manifest='manifest.csv', proposed_manifest_sha256=sha(manifest_path),
        manifest_scope='Separate method-control input only; not merged into the project library or synthesis candidate registry.',
        selected=len(entries), requested_selected=top_k, entries=entries,
        proposed_state_labels=challenge['proposed_states'],
        downstream='Nitro where supported or validated V100 fallback for warm starts, then state-specific solution DFT and native thermal QC',
        charged_state_minima_computed=False, quantum_jobs_submitted=0, experimental_accuracy_inferred=False,
        limitations='Finite force-field seeds only; no solvent, redox potential, population or conformer convergence established.')
    (out/'plan.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--challenge', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20261026)
    parser.add_argument('--attempts', type=int, default=32)
    parser.add_argument('--top-k', type=int, default=10)
    args = parser.parse_args()
    result = prepare(args.challenge.resolve(), args.out.resolve(), args.seed, args.attempts, args.top_k)
    print(json.dumps({k: result[k] for k in ('status', 'embedded', 'mmff_converged', 'selected')}))


if __name__ == '__main__':
    main()
