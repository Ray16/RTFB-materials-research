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
