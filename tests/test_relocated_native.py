import hashlib
from pathlib import Path
import pytest
from redox.qm.relocated_native import controls_equal


def fixture(tmp_path):
    root = tmp_path/'new'; p = root/'control.json';root.mkdir();p.write_text('checked')
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    return root,[dict(path='/old/control.json',sha256=sha,status='checked')],[dict(path=str(p),sha256=sha,status='checked')]


def test_only_declared_prefix_changes_and_original_is_not_opened(tmp_path):
    root, old, expected = fixture(tmp_path);receipt={}
    assert controls_equal(old,expected,root,'/old',receipt)
    assert receipt['control.json']['original_path']=='/old/control.json'


@pytest.mark.parametrize('field,value',[('path','/other/control.json'),('sha256','changed'),('status','failed')])
def test_different_source_or_metadata_is_rejected(tmp_path,field,value):
    root,old,expected=fixture(tmp_path);old[0][field]=value
    assert not controls_equal(old,expected,root,'/old',{})


def test_missing_or_changed_local_certificate_evidence_is_rejected(tmp_path):
    root,old,expected=fixture(tmp_path);(root/'control.json').write_text('changed')
    assert not controls_equal(old,expected,root,'/old',{})
    (root/'control.json').unlink()
    assert not controls_equal(old,expected,root,'/old',{})


def test_escaping_current_proof_and_ambiguous_prefix_are_rejected(tmp_path):
    root,old,expected=fixture(tmp_path)
    with pytest.raises(ValueError):controls_equal(old,expected,root,'/old/../other',{})
    expected[0]['path']=str(tmp_path/'outside')
    assert not controls_equal(old,expected,root,'/old',{})


def test_reader_binding_is_restored_after_failure():
    from redox.qm import native_thermal
    from redox.qm.relocated_native import relocated_native_reader
    original=native_thermal._load_native_surface_record
    with pytest.raises(RuntimeError,match='intentional'):
        with relocated_native_reader('/old'):
            assert native_thermal._load_native_surface_record is not original
            raise RuntimeError('intentional')
    assert native_thermal._load_native_surface_record is original
