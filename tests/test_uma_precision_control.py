import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('uma_precision_control',ROOT/'scripts/validation/qualify_uma_repeat_fp32_control.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def test_precision_control_is_process_local_and_records_effective_flags(tmp_path,monkeypatch):
    out=tmp_path/'result';monkeypatch.setattr(m.sys,'argv',['test','pair','--out',str(out)])
    monkeypatch.setenv('UMA_NITRO_FP16_RADIAL','1')
    fake=SimpleNamespace(__file__=str(Path(m.__file__).resolve()),**{name:False for name in m.FLAGS})
    def run(*args,**kwargs):
        assert all(os.environ['UMA_NITRO_'+name]=='0' for name in m.FLAGS)
        out.mkdir();monkeypatch.setitem(m.sys.modules,'test.uma_nitro._inference',fake)
    monkeypatch.setattr(m.runpy,'run_path',run)
    m.main()
    assert os.environ['UMA_NITRO_FP16_RADIAL']=='1'
    import json
    record=json.loads((out/'precision_control.json').read_text())
    assert record['all_requested_controls_observed'] is True and record['full_panel_qualified'] is False
