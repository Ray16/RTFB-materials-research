#!/usr/bin/env python3
"""Publish a source-bound interpretation comparison without new QM or selection."""
import argparse
import json
from pathlib import Path
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.context_consequences import assess

NAMES=('base_assisted.json','pmdi_first_site_grafting.json','pmdi_capacity.json','spacing_conformers.json','tempo_solution_raw_result.json')


def collect(folder):
    manifest=json.loads((folder/'manifest.json').read_text());sources={}
    for name in NAMES:
        path=folder/name;key=str(path.relative_to(folder.parents[2]))
        # Accept exactly the current snapshot's declared local artifact, never a
        # similarly named file from another calculation tree.
        matching=[p for p,h in manifest['snapshot_sha256'].items() if Path(p).name==name and h==sha256_file(path)]
        if len(matching)!=1:raise ValueError('One unchanged frozen source required: '+name)
        sources[name]=sha256_file(path)
    result=assess(*(json.loads((folder/name).read_text()) for name in NAMES))
    return dict(result,source_files_sha256=sources,source_manifest_sha256=sha256_file(folder/'manifest.json'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('Preserve earlier interpretation comparisons')
    r=collect(a.evidence.resolve());atomic_write_json(a.out,r)
    print(json.dumps(dict(status=r['status'],reaction_sign_screen_changes=r['route_definition']['changes_at_declared_nonpositive_screen'],complete_reactions=r['route_definition']['complete'],positive_polymer_penalties=r['polymer_context']['positive_penalties'],site_preference_changes_sign=r['polymer_context']['site_preference_changes_sign'])))
