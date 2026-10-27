import hashlib
import json
from pathlib import Path
import pytest
from redox.validation.evidence_epoch import append_epoch


def fixture(tmp_path):
    root=tmp_path;e=root/'manuscript/evidence/old';e.mkdir(parents=True)
    f=root/'results/figures/manuscript/old';f.mkdir(parents=True)
    (e/'raw.json').write_text('{"checked":true}');(f/'figure.pdf').write_bytes(b'preserved figure')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (e/'manifest.json').write_text(json.dumps(dict(snapshot_sha256={str(p.relative_to(root)):sha(p) for p in (e/'raw.json',f/'figure.pdf')})))
    pointer=root/'manuscript/evidence/current.json';pointer.write_text(json.dumps(dict(manifest=str((e/'manifest.json').relative_to(root)),figures=str(f.relative_to(root)))))
    (root/'added.json').write_text('{"new":true}')
    return root,pointer


def test_epoch_inherits_bytes_and_does_not_publish_its_own_pointer(tmp_path):
    root,pointer=fixture(tmp_path);before=pointer.read_bytes()
    result=append_epoch(root,'new',{'added.json':'added.json'},'Declared new evidence','revision')
    assert pointer.read_bytes()==before
    assert (root/'manuscript/evidence/new/raw.json').read_bytes()==(root/'manuscript/evidence/old/raw.json').read_bytes()
    manifest=json.loads((root/result['manifest']).read_text())
    assert manifest['scientific_qualification_inferred_from_capture'] is False
    with pytest.raises(FileExistsError):append_epoch(root,'new',{},'scope','revision')


@pytest.mark.parametrize('extras',[{'raw.json':'added.json'},{'../escape':'added.json'},{'x':'../outside'},{'manifest.json':'added.json'}])
def test_overwriting_or_escaping_evidence_is_rejected(tmp_path,extras):
    root,pointer=fixture(tmp_path)
    with pytest.raises(ValueError):append_epoch(root,'new',extras,'scope','revision')
    assert not (root/'manuscript/evidence/new').exists()


def test_tampered_parent_is_rejected_before_capture(tmp_path):
    root,pointer=fixture(tmp_path);(root/'manuscript/evidence/old/raw.json').write_text('changed')
    with pytest.raises(ValueError,match='changed'):append_epoch(root,'new',{},'scope','revision')
    assert not (root/'manuscript/evidence/new').exists()
