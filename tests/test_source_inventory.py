import hashlib
import json
from pathlib import Path

from redox.validation.source_inventory import inventory


def write(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_nested_bindings_and_exact_archive_relocation(tmp_path):
    root=tmp_path/'root';archive=tmp_path/'archive';source=archive/'raw.json'
    digest=write(source,{'value':1});seed=root/'review.json'
    seed_hash=write(seed,{'records':[{'source_sha256':{str(root/'raw.json'):digest}}]})
    r=inventory(root,archive,{seed:seed_hash})
    assert not r['unresolved'] and len(r['files'])==2
    assert r['explicit_archive_relocations'][0]['path']==str(source)


def test_missing_and_wrong_known_copy_remain_unresolved(tmp_path):
    root=tmp_path/'root';seed=root/'review.json';wrong=tmp_path/'wrong.json'
    write(wrong,{'value':2});expected='a'*64
    digest=write(seed,{'source_sha256':{str(root/'absent.json'):expected}})
    r=inventory(root,tmp_path/'archive',{seed:digest},[dict(path=str(wrong),sha256=expected)])
    assert len(r['unresolved'])==1


def test_multiple_versions_of_same_logical_path_are_not_collapsed(tmp_path):
    root=tmp_path/'root';first=tmp_path/'first.json';second=tmp_path/'second.json'
    a=write(first,{'version':1});b=write(second,{'version':2});logical=root/'old.json'
    r=inventory(root,tmp_path/'archive',[(logical,a),(logical,b)],
                [dict(path=str(first),sha256=a),dict(path=str(second),sha256=b)])
    assert len(r['files'])==2 and len(r['explicit_archive_relocations'])==2
    assert not r['unresolved']


def test_later_discovered_exact_content_resolves_earlier_missing_binding(tmp_path):
    root=tmp_path/'root';source=root/'actual.json';payload=write(source,{'value':1})
    seed=root/'review.json'
    digest=write(seed,{'source_sha256':{str(source):payload,str(root/'old.json'):payload}})
    r=inventory(root,tmp_path/'archive',{seed:digest})
    assert not r['unresolved']
    assert any(x['original']==str(root/'old.json') for x in r['explicit_archive_relocations'])


def test_unrecognized_mapping_is_not_claimed_as_covered(tmp_path):
    root=tmp_path/'root';seed=root/'review.json'
    digest=write(seed,{'other_sha256':{str(root/'absent'):'a'*64}})
    r=inventory(root,tmp_path/'archive',{seed:digest})
    assert len(r['files'])==1 and not r['implicit_dependencies_covered']


def test_incremental_inventory_checks_new_versions_and_declares_inheritance(tmp_path):
    root=tmp_path/'root';archive=tmp_path/'archive';source=root/'raw.json'
    old=write(source,{'value':1})
    base=inventory(root,archive,{source:old})
    new=write(source,{'value':2})
    seed=root/'review.json'
    digest=write(seed,{'source_sha256':{str(source):new}})
    result=inventory(root,archive,{seed:digest},baseline=base)
    assert not result['unresolved']
    assert {r['sha256'] for r in result['files'] if r['path']==str(source)}=={old,new}
    assert result['inherited_files']==1 and result['new_files_hash_checked']==2
    assert result['baseline_raw_files_rehashed'] is False
    # The inherited old version cannot satisfy a different new expected value.
    bad=inventory(root,archive,{source:'0'*64},baseline=base)
    assert len(bad['unresolved'])==1


def test_incremental_inventory_rejects_incomplete_baseline(tmp_path):
    import pytest
    base=inventory(tmp_path,tmp_path/'archive',{tmp_path/'missing':'0'*64})
    with pytest.raises(ValueError,match='complete baseline'):
        inventory(tmp_path,tmp_path/'archive',{},baseline=base)


def test_new_alias_to_inherited_content_is_checked_again(tmp_path):
    root=tmp_path/'root';source=root/'actual.json';old=write(source,{'value':1})
    base=inventory(root,tmp_path/'archive',{source:old})
    write(source,{'value':2})
    result=inventory(root,tmp_path/'archive',{root/'new_alias.json':old},baseline=base)
    assert len(result['unresolved'])==1
