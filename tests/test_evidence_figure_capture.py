import hashlib
import json
from pathlib import Path

import pytest

from redox.validation.evidence_epoch import append_epoch


@pytest.fixture
def root(tmp_path):
    old=tmp_path/'manuscript/evidence/old';fig=tmp_path/'results/figures/manuscript/old'
    old.mkdir(parents=True);fig.mkdir(parents=True)
    (old/'data.json').write_text('{}');(fig/'old.pdf').write_bytes(b'old figure bytes')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    snapshot={str(p.relative_to(tmp_path)):sha(p) for p in [old/'data.json',fig/'old.pdf']}
    (old/'manifest.json').write_text(json.dumps(dict(snapshot_sha256=snapshot)))
    (tmp_path/'manuscript/evidence/current.json').write_text(json.dumps(dict(manifest='manuscript/evidence/old/manifest.json',figures='results/figures/manuscript/old')))
    (tmp_path/'new.pdf').write_bytes(b'new scientific figure bytes');return tmp_path


def test_added_figure_is_captured_in_the_actual_compiled_figure_collection(root):
    before=(root/'manuscript/evidence/old/manifest.json').read_bytes()
    result=append_epoch(root,'new',{},'pre-reviewed numeric addition','revision',{'new.pdf':'new.pdf'})
    snapshot=json.loads((root/result['manifest']).read_text())['snapshot_sha256']
    key='results/figures/manuscript/new/new.pdf'
    assert snapshot[key]==hashlib.sha256((root/'new.pdf').read_bytes()).hexdigest()
    assert (root/key).read_bytes()==(root/'new.pdf').read_bytes()
    assert (root/'manuscript/evidence/old/manifest.json').read_bytes()==before


@pytest.mark.parametrize('name',['old.pdf','../escape.pdf'])
def test_capture_cannot_replace_an_inherited_figure_or_escape_the_collection(root,name):
    with pytest.raises(ValueError):append_epoch(root,'new',{},'scope','revision',{name:'new.pdf'})
