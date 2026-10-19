"""Reference synchronization must preserve unrelated project facts."""
import importlib.util
import json
from pathlib import Path

import pytest

from redox.properties import potentials


def updater():
    path=Path(__file__).resolve().parents[1]/'scripts/pipeline/set_fc_reference.py'
    spec=importlib.util.spec_from_file_location('reference_sync_test',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_sync_changes_only_the_derived_reference(tmp_path,monkeypatch):
    module=updater();path=tmp_path/'project.json'
    facts=dict(reference=dict(fc_abs_computed_V=4.1,other_scale='preserve'),solvent={'name':'MeCN'},note='保存')
    text=json.dumps(facts,indent=2,ensure_ascii=False)+'\n';path.write_text(text)
    monkeypatch.setattr(module,'PROJECT',path)
    monkeypatch.setattr(module,'g_solv',lambda state:{'ox':14.25,'neu':10.}[state])
    monkeypatch.setattr(potentials,'fc_live',lambda:4.25)
    monkeypatch.setattr('sys.argv',['set_fc_reference.py'])
    module.main()
    assert path.read_text()==text.replace('4.1','4.25',1)
    assert not list(tmp_path.glob('*.tmp'))


def test_sync_refuses_missing_matched_reference(tmp_path,monkeypatch):
    module=updater();path=tmp_path/'project.json';path.write_text('{"reference":{"fc_abs_computed_V":4.1}}')
    before=path.read_bytes()
    monkeypatch.setattr(module,'PROJECT',path);monkeypatch.setattr(potentials,'fc_live',lambda:None)
    monkeypatch.setattr('sys.argv',['set_fc_reference.py'])
    with pytest.raises(SystemExit):module.main()
    assert path.read_bytes()==before
