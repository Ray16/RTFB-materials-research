#!/usr/bin/env python3
"""Prepare and replay the source-bound frozen local-motif electronic panel.

Replay uses only stored records and Python's standard library. It does not rerun
SCF, physical qualification or experimental validation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
PANEL = Path('results/reviews/2026-10-26/local_motif_complete/summary.json')
CONSTRUCTION = Path('results/polymer_models/campaigns/2026-10-26_local_motif_fragments/manifest.json')
CONVERSION = 1.602176634e-19 * 6.02214076e23 / 1000


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def prepare(root, out):
    root = root.resolve(); out = out.resolve()
    if out.exists():
        raise FileExistsError('Preserve existing raw packages')
    panel = json.loads((root/PANEL).read_text())
    if (panel['status'] != 'all_twelve_frozen_local_interaction_triplets_independently_reconstructed'
            or panel['completed_triplets'] != 12 or panel['checked_electronic_geometries'] != 36):
        raise ValueError('Expected complete independently reviewed panel')
    sources = {}
    for absolute, digest in panel['source_sha256'].items():
        source = Path(absolute).resolve()
        if not source.is_relative_to(root):
            raise ValueError('Panel source lies outside project root')
        relative = str(source.relative_to(root))
        if relative in sources or sha(source) != digest:
            raise ValueError('Duplicate or changed panel source: '+relative)
        sources[relative] = digest
    for extra in (PANEL, Path(panel['source_replay']).resolve().relative_to(root), Path(__file__).resolve().relative_to(root)):
        relative = str(extra)
        value = sha(root/relative)
        if relative in sources and sources[relative] != value:
            raise ValueError('Conflicting source digest')
        sources[relative] = value
    if sha(panel['source_replay']) != panel['source_replay_sha256'] or str(CONSTRUCTION) not in sources:
        raise ValueError('Missing matching source replay or construction manifest')
    out.mkdir(parents=True)
    for relative, digest in sorted(sources.items()):
        source = root/relative; destination = out/'project'/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if sha(destination) != digest or sha(source) != digest:
            raise ValueError('Source changed while copying: '+relative)
    data = dict(status='local_motif_raw_package_prepared_not_public', original_root=str(root),
                panel=str(PANEL), construction=str(CONSTRUCTION), files_sha256=sources,
                scope='Frozen capped-fragment electronic records; no new physical calculation',
                public_archive_deposited=False, public_DOI_created=False)
    (out/'manifest.json').write_text(json.dumps(data,indent=2)+'\n')
    (out/'README.md').write_text('''# Local-motif raw electronic package\n\nThe `project/` tree preserves the relative locations and hashes of the original panel's geometry, plans, electronic outputs, checkpoint files and reviews. To replay after extraction, run `python -I project/scripts/release/local_motif_raw_archive.py replay --package PACKAGE --out NEW_RESULT.json`. The original absolute path stored in historical records is used only to map paths into this extracted tree.\n\nReplay checks source hashes, reviewed geometry and protocol identity, cached energy records and balanced electronic interaction sums. It does not rerun SCF, verify original physical QC independently, establish basis-set convergence, infer association free energies or coupling, or establish experimental accuracy. No public deposition or DOI is implied.\n''')
    return data


def replay(package):
    package = package.resolve(); manifest = json.loads((package/'manifest.json').read_text())
    original = Path(manifest['original_root'])
    files = manifest['files_sha256']; root = package/'project'
    def translated(path):
        path = Path(path)
        relative = str(path.relative_to(original)) if path.is_absolute() else str(path)
        if relative not in files or '..' in Path(relative).parts:
            raise ValueError('Unpackaged source reference: '+relative)
        result = (root/relative).resolve()
        if not result.is_relative_to(root.resolve()):
            raise ValueError('Packaged source escapes extracted project root')
        return result
    for relative, digest in files.items():
        if sha(translated(relative)) != digest:
            raise ValueError('Package source hash differs: '+relative)
    panel = json.loads(translated(manifest['panel']).read_text())
    construction = json.loads(translated(manifest['construction']).read_text())
    if (panel['completed_triplets'] != 12 or panel['checked_electronic_geometries'] != 36
            or any(panel[k] for k in ('free_energies_computed','electronic_coupling_computed',
                                      'conformer_convergence_established','basis_superposition_correction_computed'))):
        raise ValueError('Invalid panel coverage or scope')
    catalog = {}
    for relative in files:
        if not relative.endswith('.json'):
            continue
        try:
            review = json.loads(translated(relative).read_text())
        except (ValueError, UnicodeDecodeError):
            continue
        if not isinstance(review, dict) or review.get('status') != 'frozen_local_electronic_state_and_resource_pilot_independently_checked':
            continue
        if review['primary_candidate_admission'] or review['thermochemistry_computed']:
            raise ValueError('Review scope differs')
        plan = json.loads(translated(str(Path(review['campaign'])/'plan.json')).read_text())
        geometry = plan['geometry_sha256']
        if geometry in catalog:
            raise ValueError('Ambiguous reviewed geometry')
        for medium in ('gas','acetonitrile'):
            candidates = [path for path in review['source_sha256'] if path.endswith('/'+medium+'/energy.json')]
            matching = []
            for path in candidates:
                energy = json.loads(translated(path).read_text())
                if (energy.get('geometry_sha256') == geometry and energy.get('converged') is True
                        and energy.get('medium') == medium and energy.get('settings') == review['settings']
                        and abs(energy['E_eV']-review['energies'][medium]) < 1e-9):
                    matching.append(path)
            if not matching:
                raise ValueError('No matching cached electronic energy: '+review['task']['key'])
        catalog[geometry] = review
    if len(catalog) != 36:
        raise ValueError('Expected 36 distinct reviewed local electronic geometries')
    reports = {row['key']:row for row in panel['states']}
    if len(reports) != 12 or len(construction['entries']) != 12:
        raise ValueError('Expected twelve unique construction triplets')
    checks = []
    for entry in construction['entries']:
        row = reports[entry['key']]; states = {}
        for molecule in entry['molecules']:
            if sha(translated(molecule['geometry'])) != molecule['geometry_sha256']:
                raise ValueError('Geometry hash differs')
            review = catalog[molecule['geometry_sha256']]
            if review['task']['atoms'] != molecule['atoms'] or review['task']['charge'] != 0 or review['task']['mult'] != 1:
                raise ValueError('State identity differs')
            states[molecule['id']] = review
        if set(states) != {'A','B','AB'} or any(states[s]['settings'] != states['A']['settings'] for s in states):
            raise ValueError('Incomplete or mismatched protocol')
        differences = {}
        for medium in ('gas','acetonitrile'):
            actual = math.fsum((states['AB']['energies'][medium],-states['A']['energies'][medium],-states['B']['energies'][medium]))*CONVERSION
            expected = row['interaction_energy_kJ_mol_by_medium'][medium]
            if not math.isfinite(actual) or abs(actual-expected) > 1e-7:
                raise ValueError('Electronic interaction arithmetic differs')
            differences[medium] = actual
        checks.append(dict(key=entry['key'], interaction_kJ_mol=differences))
    return dict(status='local_motif_raw_records_and_interactions_replayed', files_verified=len(files),
                geometries=36, triplets=12, media=2, results=checks,
                quantum_calculations_repeated=False, original_physical_QC_repeated=False,
                association_free_energy_inferred=False, experimental_accuracy_inferred=False)


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='mode', required=True)
    a = sub.add_parser('prepare'); a.add_argument('--root', type=Path, default=ROOT); a.add_argument('--out',type=Path,required=True)
    b = sub.add_parser('replay'); b.add_argument('--package',type=Path,required=True); b.add_argument('--out',type=Path,required=True)
    args = p.parse_args()
    if args.mode == 'prepare':
        result = prepare(args.root,args.out)
        print(json.dumps(dict(status=result['status'], files=len(result['files_sha256']))))
    else:
        result = replay(args.package)
        with args.out.open('x') as stream:
            json.dump(result,stream,indent=2); stream.write('\n')
        print(json.dumps(dict(status=result['status'],files=result['files_verified'],triplets=result['triplets'])))


if __name__ == '__main__':
    main()
