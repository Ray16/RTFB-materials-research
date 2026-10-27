import builtins
import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('uma_ieee_control',ROOT/'scripts/validation/qualify_uma_repeat_ieee_control.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def test_internal_precision_and_environment_restore_after_diagnostic_failure(tmp_path,monkeypatch):
    import redox.workflows.resources as resources
    monkeypatch.setattr(resources,'verify_gate_owner',lambda:None)
    out=tmp_path/'result';monkeypatch.setattr(m.sys,'argv',['diagnostic','pair','--out',str(out)])
    runtime=SimpleNamespace(_USE_TF32=True,__file__=m.__file__,**{name:False for name in m.FLAGS})
    original=builtins.__import__
    def importing(name,*args,**kwargs):
        if name=='fairchem.core.models.uma.uma_nitro':return SimpleNamespace(_inference=runtime)
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',importing)
    monkeypatch.setenv('UMA_NITRO_FP16_RADIAL','1')
    def run(*args,**kwargs):
        assert runtime._USE_TF32 is False and os.environ['UMA_NITRO_FP16_RADIAL']=='0'
        out.mkdir();raise RuntimeError('diagnostic failed')
    monkeypatch.setattr(m.runpy,'run_path',run)
    with pytest.raises(RuntimeError,match='diagnostic failed'):m.main()
    assert runtime._USE_TF32 is True and os.environ['UMA_NITRO_FP16_RADIAL']=='1'
    import json
    record=json.loads((out/'precision_control.json').read_text())
    assert record['official_default_profile'] is False and record['full_panel_qualified'] is False
