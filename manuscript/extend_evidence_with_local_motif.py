#!/usr/bin/env python3
"""Freeze the independently reviewed local-motif electronic interaction panel."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extend(name):
    if not name or Path(name).name != name:
        raise ValueError('One snapshot name required')
    pointer = json.loads((ROOT/'manuscript/evidence/current.json').read_text())
    parent = ROOT/pointer['manifest']
    prior = json.loads(parent.read_text())
    old_figures = ROOT/pointer['figures']
    evidence = ROOT/'manuscript/evidence'/name
    figures = ROOT/'results/figures/manuscript'/name
    if evidence.exists() or figures.exists():
        raise FileExistsError('Never overwrite an evidence epoch')
    mapping = {}
    for path, expected in prior['snapshot_sha256'].items():
        source = ROOT/path
        if sha(source) != expected:
            raise ValueError('Parent snapshot changed: '+path)
        if source.parent == parent.parent:
            destination = evidence/source.name
        elif source.parent == old_figures:
            destination = figures/source.name
        else:
            raise ValueError('Unexpected parent snapshot location: '+path)
        mapping[destination] = source
    panel = ROOT/'results/reviews/2026-10-26/local_motif_complete/summary.json'
    d = json.loads(panel.read_text())
    if (d['status'] != 'all_twelve_frozen_local_interaction_triplets_independently_reconstructed'
            or d['completed_triplets'] != 12 or d['checked_electronic_geometries'] != 36
            or any(d[k] for k in ('free_energies_computed', 'electronic_coupling_computed',
                                  'conformer_convergence_established', 'basis_superposition_correction_computed'))):
        raise ValueError('Local-motif panel is incomplete or its scope differs')
    if sha(d['source_replay']) != d['source_replay_sha256']:
        raise ValueError('Local-motif source replay changed')
    for path, digest in d['source_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Local-motif source changed: '+path)
    mapping[evidence/'local_motif_panel.json'] = panel
    identities = {str(p.relative_to(ROOT)): sha(p) for p in set(mapping.values())}
    evidence.mkdir(); figures.mkdir()
    for destination, source in mapping.items():
        shutil.copyfile(source, destination)
    if any(sha(ROOT/path) != digest for path, digest in identities.items()):
        raise RuntimeError('Source changed during extension')
    manifest = dict(kind='immutable_draft_evidence_snapshot',
        captured_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        parent_manifest=pointer['manifest'], parent_manifest_sha256=sha(parent),
        source_sha256=identities, snapshot_sha256={str(p.relative_to(ROOT)):sha(p) for p in mapping},
        scope='Prior evidence plus frozen phthalimide local-motif electronic interactions; no new quantum results')
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(manifest=str((evidence/'manifest.json').relative_to(ROOT)),figures=str(figures.relative_to(ROOT)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name', required=True)
    a = p.parse_args()
    print(json.dumps(extend(a.name), indent=2))
