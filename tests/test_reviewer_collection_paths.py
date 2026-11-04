"""Explicit review outputs may live inside the canonical project."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('review_collection',ROOT/'scripts/release/reviewer_collection.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_fresh_canonical_destination_is_allowed_without_creating_it(tmp_path):
    archive=tmp_path/'archive';archive.mkdir()
    output=ROOT/'calcs/validation/review_destination_not_created_by_test'
    assert not output.exists()
    assert module.output_location(archive,output)==(archive,output)
    assert not output.exists()


@pytest.mark.parametrize('destination',['archive','archive/child','existing'])
def test_existing_or_nested_source_destinations_are_refused(tmp_path,destination):
    archive=tmp_path/'archive';archive.mkdir();(tmp_path/'existing').mkdir()
    with pytest.raises(FileExistsError):module.output_location(archive,tmp_path/destination)


def test_changed_relocated_archive_is_refused_before_any_output(tmp_path,monkeypatch):
    import hashlib
    import json
    root=tmp_path/'project';root.mkdir();archive=tmp_path/'inputs';archive.mkdir()
    monkeypatch.setattr(module,'ROOT',root)
    monkeypatch.setattr(module.subprocess,'check_output',lambda *a,**k:'fixture-source-commit\n')
    monkeypatch.setattr(module.subprocess,'Popen',lambda *a,**k:pytest.fail('Changed input reached source export'))
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in ('primary_review','first_site','controls','neighbor_cycle','tempo_addendum','local_motif'):
        payload=archive/(name+'.tar.zst');payload.write_bytes(b'checked fixture bytes')
        record=archive/'results/releases'/name/'current.json';record.parent.mkdir(parents=True)
        record.write_text(json.dumps({'compressed_archive':payload.name,'archive_sha256':digest(payload)}))
    new=root/'calcs/native_inputs';new.mkdir(parents=True)
    payload=new/'native.tar.zst';payload.write_bytes(b'original native archive');expected=digest(payload)
    old=tmp_path/'legacy_native_inputs'
    record=root/'validation/native_cache_2026-10-30/summary.json';record.parent.mkdir(parents=True)
    record.write_text(json.dumps({'archive':str(old/'native.tar.zst'),'archive_sha256':expected}))
    original_record=record.read_bytes()
    cpu=root/'results/redox_release_candidates/2026-10-31_exact_viologen_CPU_final.tar.zst'
    cpu.parent.mkdir(parents=True);cpu.write_bytes(b'checked CPU fixture')
    info=root/'validation/exact_molecular_CPU_2026-10-31/summary.json';info.parent.mkdir(parents=True)
    info.write_text(json.dumps({'archive_sha256':digest(cpu)}))
    payload.write_bytes(b'changed native archive')
    output=root/'results/new_review'
    mapping={'moved':[{'old':str(old),'new':str(new)}]}
    with pytest.raises(ValueError,match='checked package changed'):
        module.assemble(archive,output,'HEAD',relocation=mapping)
    assert not output.exists() and not old.exists() and record.read_bytes()==original_record
