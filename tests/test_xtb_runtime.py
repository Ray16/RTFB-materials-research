import hashlib
from pathlib import Path

import pytest

from redox.qm import xtb_runtime


def test_interpreter_prefix_discovery_and_explicit_paths(tmp_path,monkeypatch):
    prefix=tmp_path/'environment';(prefix/'lib').mkdir(parents=True);(prefix/'bin').mkdir()
    (prefix/'lib/libxtb.so.6.7.1').write_bytes(b'qualified separately');(prefix/'bin/xtb').write_bytes(b'qualified separately')
    monkeypatch.setattr(xtb_runtime.sys,'prefix',str(prefix));monkeypatch.setattr(xtb_runtime.sys,'base_prefix',str(tmp_path/'base'))
    assert xtb_runtime.configured_library({})==prefix/'lib/libxtb.so.6.7.1'
    assert xtb_runtime.configured_binary({})==prefix/'bin/xtb'
    chosen=tmp_path/'explicit'
    assert xtb_runtime.configured_library({'REDOX_XTB_PREFIX':str(chosen)})==chosen/'lib/libxtb.so.6.7.1'
    assert xtb_runtime.configured_binary({'REDOX_XTB_PREFIX':str(chosen)})==chosen/'bin/xtb'
    assert xtb_runtime.configured_library({'REDOX_XTB_LIBRARY':str(tmp_path/'custom.so')})==tmp_path/'custom.so'
    assert xtb_runtime.configured_binary({'REDOX_XTB_EXECUTABLE':str(tmp_path/'custom_xtb')})==tmp_path/'custom_xtb'


@pytest.mark.parametrize('name',['REDOX_XTB_PREFIX','REDOX_XTB_LIBRARY','REDOX_XTB_EXECUTABLE'])
@pytest.mark.parametrize('value',['','relative/path'])
def test_ambiguous_runtime_overrides_are_refused(name,value):
    call={'REDOX_XTB_PREFIX':xtb_runtime.configured_prefix,'REDOX_XTB_LIBRARY':xtb_runtime.configured_library,'REDOX_XTB_EXECUTABLE':xtb_runtime.configured_binary}[name]
    with pytest.raises(ValueError,match='absolute'):call({name:value})


def test_prefix_rejects_changed_bytes_and_files_outside_prefix(tmp_path):
    prefix=tmp_path/'prefix';(prefix/'bin').mkdir(parents=True);(prefix/'lib').mkdir()
    lib=prefix/'lib/libxtb.so.6.7.1';binary=prefix/'bin/xtb';lib.write_bytes(b'lib');binary.write_bytes(b'bin');binary.chmod(0o755)
    h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    record=xtb_runtime.validate_prefix(prefix,h(lib),h(binary));assert record['prefix']==str(prefix)
    lib.write_bytes(b'changed')
    with pytest.raises(ValueError,match='qualified bytes'):xtb_runtime.validate_prefix(prefix,hashlib.sha256(b'lib').hexdigest(),h(binary))
    lib.unlink();outside=tmp_path/'outside';outside.write_bytes(b'lib');lib.symlink_to(outside)
    with pytest.raises(ValueError,match='inside'):xtb_runtime.validate_prefix(prefix,h(outside),h(binary))


def test_binding_restores_paths_and_environment_without_modifying_other_source(tmp_path,monkeypatch):
    from redox.qm import qrrho,native_thermal,projected_thermal
    old=(qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB)
    monkeypatch.setattr(qrrho,'_WORKER',None);monkeypatch.setattr(qrrho,'_LAST_WORKER',None)
    monkeypatch.setattr(xtb_runtime,'validate_prefix',lambda *a,**k:{'library':str(tmp_path/'lib.so'),'binary':str(tmp_path/'xtb')})
    def checked():return str(tmp_path/'lib.so')
    checked.cache_clear=lambda:None
    monkeypatch.setattr(qrrho,'_checked_library',checked);monkeypatch.setenv('LD_LIBRARY_PATH','before')
    with xtb_runtime.qualified_xtb_prefix(tmp_path):
        assert qrrho.LIBRARY==tmp_path/'lib.so' and native_thermal.XTB==projected_thermal.XTB==tmp_path/'xtb'
    assert (qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB)==old
    import os
    assert os.environ['LD_LIBRARY_PATH']=='before'
    with pytest.raises(RuntimeError,match='intentional'):
        with xtb_runtime.qualified_xtb_prefix(tmp_path):raise RuntimeError('intentional')
    assert (qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB)==old
    assert os.environ['LD_LIBRARY_PATH']=='before'


def test_binding_refuses_a_preexisting_thermodynamic_worker(tmp_path,monkeypatch):
    from redox.qm import qrrho
    monkeypatch.setattr(qrrho,'_WORKER',object())
    with pytest.raises(RuntimeError,match='before the first'):
        with xtb_runtime.qualified_xtb_prefix(tmp_path):pass


def test_paths_restore_even_if_helper_shutdown_fails(tmp_path,monkeypatch):
    from redox.qm import qrrho,native_thermal,projected_thermal
    old=(qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB)
    monkeypatch.setattr(qrrho,'_WORKER',None);monkeypatch.setattr(qrrho,'_LAST_WORKER',None)
    monkeypatch.setattr(xtb_runtime,'validate_prefix',lambda *a,**k:{'library':str(tmp_path/'lib.so'),'binary':str(tmp_path/'xtb')})
    def checked():return str(tmp_path/'lib.so')
    checked.cache_clear=lambda:None
    monkeypatch.setattr(qrrho,'_checked_library',checked)
    monkeypatch.setattr(qrrho,'_close_worker',lambda:(_ for _ in ()).throw(RuntimeError('shutdown failed')))
    monkeypatch.setenv('LD_LIBRARY_PATH','before')
    with pytest.raises(RuntimeError,match='shutdown failed'):
        with xtb_runtime.qualified_xtb_prefix(tmp_path):pass
    import os
    assert (qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB)==old
    assert os.environ['LD_LIBRARY_PATH']=='before'
