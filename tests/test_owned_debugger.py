import json
import os

import pytest

from redox.workflows import owned_debugger as module


def test_wrong_debugger_executable_cannot_receive_permission(tmp_path,monkeypatch):
    (tmp_path/'debugger_owner.json').write_text(json.dumps(dict(pid=os.getpid(),process_start=module.process_start(os.getpid()))))
    monkeypatch.setattr(module.ctypes,'CDLL',lambda *a,**k:pytest.fail('Permission API must not be called for an unverified debugger'))
    with pytest.raises(ValueError,match='identity differs'):module.permit(tmp_path,'0'*64,timeout=1)


@pytest.mark.parametrize('change',['start','command'])
def test_an_unmatched_target_cannot_be_attached(tmp_path,monkeypatch,change):
    (tmp_path/'debugger_wait.json').write_text(json.dumps(dict(pid=os.getpid(),process_start='0' if change=='start' else module.process_start(os.getpid()))))
    monkeypatch.setattr(module.subprocess,'Popen',lambda *a,**k:pytest.fail('No debugger may spawn for an unverified target'))
    with pytest.raises(ValueError,match='exact waiting diagnostic'):module.start(tmp_path,{},'unmatched-case-plan.json')
