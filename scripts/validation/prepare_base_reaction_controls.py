#!/usr/bin/env python3
"""Extract exact headline base-reaction electronic inputs from the checked package.

Only geometry and active-energy bytes are copied. The output is an input package
for implementation sensitivity, not a portable thermochemistry package or new G.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import shutil
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def prepare(source,output):
    if output.exists():raise FileExistsError('Preserve earlier extracted inputs')
    manifest=source/'manifest.json';data=json.loads(manifest.read_text())
    base=[r for r in data['reactions'] if r['name'].startswith(('dbu:','dipea:'))]
    if len(base)!=20:raise ValueError('Exactly the twenty declared base-assisted models required')
    ids={r['name'].split(':')[-1] for r in base}
    refs=[r for r in data['reactions'] if r['name'] in {'reference:'+name for name in ids}]
    if len(refs)!=5:raise ValueError('Five exact matching HCl references required')
    required={key for r in base+refs for key in r['coefficients']}
    aliases={alias:row for row in data['states'] for alias in row['aliases']}
    fc={aliases['benchmark:ferrocene:'+state]['id']:'fc_'+state for state in ('neu','ox')}
    required.update(fc);mapping={key:fc.get(key,key) for key in required};rows=[];hashes={}
    if any(not isinstance(key,str) or not key or Path(key).name!=key or key in ('.','..') for key in mapping.values()):
        raise ValueError('Safe exact state IDs required')
    output.mkdir(parents=True)
    for row in data['states']:
        if row['id'] not in required:continue
        files={}
        for name in ('geometry','energy'):
            relative=row['files'][name];original=source/relative
            if not original.resolve().is_relative_to(source.resolve()):raise ValueError('Source input escapes package')
            if sha256_file(original)!=data['files_sha256'][relative]:raise ValueError('Checked source input changed')
            target=output/'states'/mapping[row['id']]/(name+original.suffix);target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(original,target);relative=str(target.relative_to(output))
            hashes[relative]=sha256_file(target);files[name]=relative
        rows.append(dict(id=mapping[row['id']],aliases=row['aliases'],smiles=row['smiles'],charge=row['charge'],mult=row['mult'],files=files,
            original_source_state=row['id'],original_state_directory=row['original_state_directory']))
    if {row['id'] for row in rows}!=set(mapping.values()):raise ValueError('Every declared state required')
    reactions=[dict(id=r['name'],kind='grafting_free_energy',coefficients={mapping[k]:v for k,v in r['coefficients'].items()},
                    target_delta_n_solutes=r['delta_n_solutes']) for r in refs+base]
    result=dict(schema=1,scope_id='headline_base_reaction_components',states=rows,reactions=reactions,files_sha256=hashes,
        source_manifest=str(manifest),source_manifest_sha256=sha256_file(manifest),
        native_thermochemistry_copied=False,new_G_qualified=False,scope='Electronic implementation contributions to 25 declared reaction sums only; no nuclear or standard-state correction is recomputed')
    atomic_write_json(output/'manifest.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=prepare(a.source.resolve(),a.out.resolve());print(json.dumps(dict(states=len(r['states']),target_cycles=len(r['reactions']),new_G_qualified=False)))
