"""Preserve immutable replay references when extending manuscript evidence."""
import hashlib
import json

import pytest

from redox.validation.evidence_epoch import append_epoch


def make_parent(root):
    parent=root/'manuscript/evidence/old'
    figures=root/'results/figures/manuscript/old'
    parent.mkdir(parents=True);figures.mkdir(parents=True)
    (parent/'value.json').write_text('{"value": 3}')
    (figures/'figure.pdf').write_bytes(b'fixed figure bytes')
    external=root/'raw/graph/input.xyz';external.parent.mkdir(parents=True);external.write_text('1\nfixed raw structure\nH 0 0 0\n')
    paths=[parent/'value.json',figures/'figure.pdf',external]
    hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (parent/'manifest.json').write_text(json.dumps(dict(snapshot_sha256=hashes)))
    pointer=dict(manifest='manuscript/evidence/old/manifest.json',figures='results/figures/manuscript/old')
    (root/'manuscript/evidence/current.json').write_text(json.dumps(pointer))
    return external,hashes


def test_extension_carries_exact_bound_raw_inputs_without_rewriting_them(tmp_path):
    external,hashes=make_parent(tmp_path)
    before=external.read_bytes()
    result=append_epoch(tmp_path,'new',{},'same scientific scope','source-commit')
    manifest=json.loads((tmp_path/result['manifest']).read_text())
    key=str(external.relative_to(tmp_path))
    assert manifest['snapshot_sha256'][key]==hashes[key]
    assert manifest['inherited_bound_references_sha256']=={key:hashes[key]}
    assert external.read_bytes()==before
    assert (tmp_path/'manuscript/evidence/new/value.json').read_bytes()==(tmp_path/'manuscript/evidence/old/value.json').read_bytes()


def test_changed_bound_raw_input_rejects_extension_before_writing(tmp_path):
    external,_=make_parent(tmp_path)
    external.write_text('changed graph')
    with pytest.raises(ValueError,match='Parent evidence changed'):
        append_epoch(tmp_path,'new',{},'same scientific scope','source-commit')
    assert not (tmp_path/'manuscript/evidence/new').exists()
