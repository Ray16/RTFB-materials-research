"""Workflow schema adaptation preserves canonical receipts and quantum bytes."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('energy_marker',ROOT/'scripts/pipeline/normalize_polymer_energy_handoff.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def setup(tmp_path,monkeypatch):
    (tmp_path/'plan.json').write_text(json.dumps(dict(tasks=[dict(key='k')],producer_sha256={},input_sha256={})))
    status=tmp_path/'status.json';status.write_text(json.dumps(dict(status='state_records_verified',physical_QC='passed')))
    done=tmp_path/'energy/done';done.mkdir(parents=True);receipt=done/'k.json'
    receipt.write_text(json.dumps(dict(status_path=str(status))))
    marker=dict(status='stationary_electronic_case_checked_native_G_pending',source_receipt=str(receipt),source_receipt_sha256=MODULE.api.digest(receipt),native_G_computed=False)
    (tmp_path/'workspace_energy_status.json').write_text(json.dumps(marker))
    calls=[];monkeypatch.setattr(MODULE,'verify_state_checkpoint',lambda p:calls.append(p))
    return status,receipt,marker,calls


def test_adapter_is_explicit_idempotent_and_leaves_quantum_parent_unchanged(tmp_path,monkeypatch):
    status,receipt,marker,calls=setup(tmp_path,monkeypatch)
    previous={p:p.read_bytes() for p in (status,receipt)}
    MODULE.normalize(tmp_path)
    data=json.loads((tmp_path/'workspace_energy_status.json').read_text())
    assert data['energy_receipt_sha256']==MODULE.api.digest(status)
    assert json.loads((tmp_path/'handoff_schema/original_marker.json').read_text())==marker
    assert data['schema_adapter']['quantum_records_changed'] is False
    assert calls==[receipt] and all(p.read_bytes()==value for p,value in previous.items())
    before=(tmp_path/'workspace_energy_status.json').read_bytes();MODULE.normalize(tmp_path)
    assert (tmp_path/'workspace_energy_status.json').read_bytes()==before


def test_changed_canonical_receipt_is_rejected_before_any_marker_edit(tmp_path,monkeypatch):
    status,receipt,marker,calls=setup(tmp_path,monkeypatch)
    marker['source_receipt_sha256']='wrong'
    path=tmp_path/'workspace_energy_status.json';path.write_text(json.dumps(marker));before=path.read_bytes()
    with pytest.raises(ValueError):MODULE.normalize(tmp_path)
    assert path.read_bytes()==before and not (tmp_path/'handoff_schema').exists()
